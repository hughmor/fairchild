use crate::mna::MnaMatrix;

pub const K_BOLTZMANN: f64 = 1.380649e-23;
pub const Q_ELECTRON: f64 = 1.602176634e-19;

/// Index of a terminal in the MNA solution vector; `None` → ground (excluded from matrix).
pub type NodeId = Option<usize>;

/// One noise generator that injects into several places at once, all driven by
/// the SAME underlying random process.
///
/// The taps of one `CorrelatedNoise` add **coherently** —
/// `|Σ wₖ·(λ[posₖ] − λ[negₖ])|² · psd` — which is what separates it from
/// returning several entries from [`Device::noise_sources`], where each entry
/// is independent and the transfer magnitudes add in quadrature.
///
/// Laser RIN is the case that needs it: one intensity fluctuation `δP` lands on
/// both the `re` and `im` field wires, split by the emission phase.  At φ₀ = 0
/// the two forms happen to agree; at 45° the quadrature sum is √2 low.
pub struct CorrelatedNoise {
    /// One-sided PSD of the driving process, in the squared units of whatever
    /// the taps inject — A²/Hz for a current into a node, or the square of the
    /// enforced potential's unit for an injection into a branch row.
    pub psd: f64,
    /// `(pos, neg, weight)`.  A `neg` of `None` is ground.
    pub taps: Vec<(NodeId, NodeId, f64)>,
}

/// Simulator context passed to device model callbacks at every eval.
pub struct SimContext {
    pub temperature: f64, // Kelvin; default 300.15 K (27 °C, SPICE TNOM)
    pub omega_0: f64,     // rad/s carrier frequency for optical ports; 0 for electrical
    /// When true, device models apply junction-step limiters (pnjlim, fetlim).
    /// Mapped from `SimOptions::pnjlim` / `.options nopnjlim`.
    pub jlim_enabled: bool,
    /// Minimum conductance placed **across each pn junction** — SPICE's `GMIN`.
    ///
    /// A junction's `dI/dV` collapses to nothing in reverse bias, so a node that
    /// reaches the circuit only through reverse-biased junctions gets a row of
    /// almost zeros and the Jacobian goes near-singular. `gmin` in parallel keeps
    /// every junction conducting something.
    ///
    /// It reaches devices through the context because it belongs to the *solve*,
    /// not to the model card: `diode.rs`, `bjt.rs` and `mosfet1.rs` each used to
    /// hardcode their own `const GMIN: f64 = 1e-12`, so `.options gmin=` moved a
    /// nodal floor and left every junction alone.
    ///
    /// A junction is what it crosses, so `mosfet1.rs` — Level 1, no body diodes
    /// — has nowhere to put it and keeps it as a Jacobian-only channel floor;
    /// see the note there.
    pub gmin: f64,
    /// Centre wavelength of the photonic band of interest, in metres.  Devices
    /// use this to bootstrap their λ wire during the initial NR iterate (when
    /// the laser hasn't driven the wire yet — value at x=0).  After iteration 0
    /// the actual wire value wins.  Set via `.options lambda_center_nm=1310`
    /// for O-band designs, `lambda_center_nm=1550` (default) for C-band.
    pub lambda_center_m: f64,
    /// When true, optical bundles carry 5 wires per channel
    /// (re_fw, im_fw, re_bw, im_bw, λ) and every photonic device stamps
    /// both forward and backward paths.  When false (default), bundles
    /// carry 3 wires (re, im, λ) and devices only stamp the forward
    /// direction.  Sourced from `SimOptions::bidirectional_propagation`.
    pub bidirectional_propagation: bool,
    /// When true, devices exposing a group delay (optical waveguides) model it
    /// as a true delay line rather than an instantaneous transmission.  Sourced
    /// from `SimOptions::waveguide_delay`.
    pub waveguide_delay: bool,
    /// Absolute transient time of the step currently being evaluated (s).  The
    /// transient loop updates this before each device `eval`/`load_*` so delay
    /// lines can look up historical port values at `time_s - τ`.  Zero during
    /// DC operating-point and AC analyses.
    pub time_s: f64,
    /// How the integrator is discretising the current timestep, for devices
    /// that stamp their own reactive companion rather than declaring
    /// [`Device::reactive_branches`].
    ///
    /// `None` outside transient. Set immediately before each `eval`, which such
    /// a device is already required to have called before `load_*_tran`.
    pub discretisation: Option<Discretisation>,
}

/// The integrator's discretisation of one timestep.
///
/// The `alpha` handed to [`Device::load_residual_tran`] cannot express anything
/// but Backward Euler: Trapezoidal and BDF-2 need history terms there is no room
/// for in a single scalar. A device that stamps its own reactance therefore
/// needs the method itself, and gets it here rather than through 18 changed
/// signatures.
///
/// Feed it to [`crate::reactive::charge_current`] — the one place a method is
/// interpreted for charge-based reactance.
#[derive(Clone, Copy, Debug)]
pub struct Discretisation {
    pub mode: crate::tran::IntegratorMode,
    /// The step being taken, in seconds.
    pub h: f64,
    /// `Some(h_prev)` only when BDF-2 is permitted this step — GEAR, no recent
    /// rejection, sane step ratio. Absent, GEAR demotes to BE.
    pub gear2_h_prev: Option<f64>,
}

impl Default for SimContext {
    fn default() -> Self {
        SimContext {
            temperature: 300.15,
            omega_0: 0.0,
            jlim_enabled: true,
            // Matches `SimOptions::default().gmin`, asserted by
            // `the_two_defaults_for_gmin_are_one_number` in `options.rs` — two
            // defaults for one quantity is two chances to drift.
            gmin: 1e-12,
            lambda_center_m: 1.55e-6,
            bidirectional_propagation: false,
            waveguide_delay: false,
            time_s: 0.0,
            discretisation: None,
        }
    }
}

impl SimContext {
    /// Thermal voltage kT/q in volts.
    pub fn vt(&self) -> f64 {
        K_BOLTZMANN * self.temperature / Q_ELECTRON
    }

    /// Number of underlying SVEA wires that a single bundle channel occupies.
    /// 3 for unidirectional (re, im, λ); 5 for bidirectional (re_fw, im_fw,
    /// re_bw, im_bw, λ).  Photonic devices that derive `n_channels` from
    /// `terminals.len()` should query this to know the per-channel stride.
    pub fn wires_per_channel(&self) -> usize {
        if self.bidirectional_propagation {
            5
        } else {
            3
        }
    }
}

/// Kind of a device-internal reactive branch the integrator should manage
/// with its companion-model machinery.  See [`ReactiveBranchSpec`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReactiveKind {
    Capacitor,
    Inductor,
}

/// Specification for a linear (or quasi-linear) reactive branch contributed
/// by a device.  The integrator owns the companion-model state for each
/// branch — the device only declares "I have a capacitor between (pos, neg)
/// with current value C(V_op)" and the transient solver does the BE / TR /
/// BDF-2 companion stamping AND advance from the post-converged solution.
///
/// `value` is queried EVERY NR iteration via [`Device::reactive_branches`],
/// so bias-dependent capacitance (e.g., depletion C_j(V) on a PN junction)
/// works naturally — return `C_j(V_pn)` evaluated at the device's current
/// cached operating-point voltage.  Newton converges through the value's
/// nonlinearity the same way it converges through I(V) nonlinearity in
/// `load_residual` / `load_jacobian`.
///
/// For state-variable physics that doesn't reduce to a linear C·dV/dt
/// (carrier dynamics, thermal RC with self-heating from absorbed light),
/// the device instead owns its own state through `commit_timestep` and
/// stamps the discretised state equation directly in `load_jacobian_tran`.
/// See the L3 PN-PS implementation for the canonical example (when it lands).
/// One complex admittance entry at one frequency, from
/// [`Device::ac_stamps`].
///
/// `row`/`col` are MNA indices and `(re, im)` is the complex coefficient added
/// at that cell — so the real part lands in the `G` block and the imaginary
/// part in the susceptance block, exactly as `jωC` would.
#[derive(Clone, Copy, Debug)]
pub struct AcStamp {
    pub row: usize,
    pub col: usize,
    pub re: f64,
    pub im: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct ReactiveBranchSpec {
    pub kind: ReactiveKind,
    pub pos: NodeId,
    pub neg: NodeId,
    /// Current capacitance (F) or inductance (H) at the device's cached
    /// operating point.  Re-queried per NR iteration.
    pub value: f64,
    /// `dC/dv` (F/V), or `dL/di` (H/A), at that same operating point.  Zero for
    /// a branch whose value does not depend on its own state, which is most of
    /// them — and the default, so a device that has no bias dependence says
    /// nothing.
    ///
    /// **Load-bearing for the Jacobian, not for the residual.**  The branch
    /// carries `q = value·v`, so its current is `α·(C(v)·v − C_prev·v_prev)`
    /// and the true derivative is `α·(C + v·dC/dv)`.  Stamping `α·C` alone
    /// converges — Newton reaches the same fixed point by successive
    /// substitution on the missing term — so a forward run looks correct and
    /// only pays in iterations.  What it is not is the *true* `∂f/∂x`, and the
    /// adjoint needs that: `dL/dp = −λᵀ·∂f/∂p` is the total derivative only if
    /// `Jᵀλ = ∂L/∂x` was solved with the real `J`.  Measured on an MZI with a
    /// `fc_pn_ps_cap` arm, the missing term is 22 % of the drive node's
    /// diagonal and put the parameter gradient 16 % out.
    pub dvalue_dstate: f64,
    /// Stored charge (C) or flux (Wb) at the current iterate, when `value·state`
    /// is not it.
    ///
    /// **This is the state variable the integrator advances**, and for a
    /// junction it is not `C(v)·v`. The depletion charge is `q(v) = ∫C dv`, and
    /// with `C ∝ (1 − v/V_bi)^−m` the two differ by a factor of 2.3 over a 2 V
    /// step. Integrating `C(v)·v` makes the device *faster* than it is — a
    /// `fc_pn_ps_cap` arm that should take 46 ps to cross 10-90 through 50 Ω
    /// took 20 ps, and every eye, bandwidth and edge measured through it was
    /// wrong in the optimistic direction, with nothing to notice.
    ///
    /// `None` means `q = value·state`, which is exact for a linear branch and
    /// is what every constant `C` or `L` reports. A branch reporting a charge
    /// also gets its Jacobian from `value` alone (`∂i/∂v = α·dq/dv = α·C`), so
    /// it leaves [`Self::dvalue_dstate`] at zero: the correction that field
    /// exists for is a correction to the wrong charge model.
    pub charge: Option<f64>,
}

/// Bit-flags controlling which contributions `eval` should compute.
pub struct EvalFlags {
    /// Compute resistive (DC / quasi-static) contributions — I(V) and dI/dV.
    pub resistive: bool,
    /// Compute reactive (capacitive/inductive) contributions — dq/dt and dq/dV.
    pub transient: bool,
}

impl EvalFlags {
    pub fn dc() -> Self {
        EvalFlags {
            resistive: true,
            transient: false,
        }
    }

    pub fn tran() -> Self {
        EvalFlags {
            resistive: true,
            transient: true,
        }
    }
}

/// Unified interface for circuit device models — both OSDI-loaded and Rust-native.
///
/// Call order per Newton-Raphson iteration:
///   1. `eval`           — device evaluates physics at current `x`, caches Norton equivalent
///   2. `load_residual`  — stamps cached currents into `b`
///   3. `load_jacobian`  — stamps cached conductances into `mat.a`
pub trait Device: Send + Sync {
    fn num_terminals(&self) -> usize;

    /// Finalise model-level derived quantities (temperature-scaled params, etc.).
    /// Called once after construction before any `setup_instance`.
    fn setup_model(&mut self, ctx: &SimContext);

    /// Bind this instance to specific MNA node indices.
    /// `terminals` length must equal `num_terminals()`; `None` = connected to ground.
    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext);

    /// Evaluate device physics at operating point `x`.
    /// Results are cached internally for subsequent `load_residual` / `load_jacobian` calls.
    fn eval(&mut self, x: &[f64], flags: EvalFlags, ctx: &SimContext);

    /// Add the Norton-equivalent current source contribution to the residual vector `b`.
    fn load_residual(&self, b: &mut [f64]);

    /// Add the Norton-equivalent conductance contribution to the MNA Jacobian `mat`.
    fn load_jacobian(&self, mat: &mut MnaMatrix);

    /// Transient residual: stamp reactive + resistive contributions for BE companion.
    /// `alpha` = h_new/h_old (1.0 for fixed-step BE).
    /// Default falls back to DC residual (correct for purely resistive devices).
    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }

    /// Transient Jacobian: stamp resistive + α·reactive conductances.
    /// Default falls back to DC Jacobian.
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }

    /// Record the current solution as the previous-timestep reference for reactive terms.
    ///
    /// Called once per accepted timestep in the variable-step integrator.
    /// Default is a no-op (correct for resistive-only and built-in devices).
    fn commit_timestep(&mut self, _x: &[f64]) {}

    /// Linear (or bias-dependent linear) reactive branches the device
    /// contributes to the MNA matrix.  The integrator stamps companion-
    /// model (G_eq, I_hist) between (pos, neg) of each branch every NR
    /// iteration in transient analysis, using the device-reported `value`
    /// (typically C(V_op) at the current iterate).  After each successful
    /// timestep the integrator reads `V_C = x[pos] − x[neg]` from the
    /// converged solution and updates the history.
    ///
    /// Default empty.  Override for devices with linear C(V) or L(I)
    /// contributions (depletion C_j on a PN, parasitic c_par on a PD).
    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> {
        Vec::new()
    }

    /// Small-signal reactive contributions at the current operating point, for
    /// frequency-domain analyses (`.ac`, `.noise`).
    ///
    /// Reported as capacitances (∂Q/∂V, farads) and inductances (∂Φ/∂I, henries)
    /// between node pairs — the *same* physical reactances the device stamps in
    /// transient, whether through [`Device::reactive_branches`] (the integrator-
    /// companion path) or through its own [`Device::load_jacobian_tran`]
    /// companion stamping. AC and noise build their jωC / 1/(jωL) blocks from
    /// these, so device-internal caps (diode C_j, MOSFET Meyer/junction caps,
    /// photonic parasitics) are included rather than silently dropped.
    ///
    /// **Query after an [`Device::eval`] with [`EvalFlags::tran`] at the
    /// operating point**, so devices whose cap evaluation is transient-gated
    /// have populated their cached values. The default returns
    /// [`Device::reactive_branches`] — correct for devices that already expose
    /// their reactance through the integrator-companion path (e.g. photonic
    /// phase-shifter junction caps); devices that stamp companions themselves
    /// (diode, MOSFET) override this.
    fn small_signal_reactances(&self) -> Vec<ReactiveBranchSpec> {
        self.reactive_branches()
    }

    /// Stamp this device's reactive Jacobian ∂q/∂x directly into the
    /// small-signal capacitance matrix, for devices whose reactance is a
    /// general matrix rather than a set of two-terminal branches.
    ///
    /// `.ac` and `.noise` form their susceptance block as `ω·C − L/ω`, so what
    /// lands here is the frequency-domain twin of the `α·J_react` that
    /// [`Device::load_jacobian_tran`] stamps: same entries, same positions,
    /// α → jω. Indices are MNA rows/columns, as for
    /// [`Device::load_jacobian`].
    ///
    /// Default is a no-op. Override this **instead of**
    /// [`Device::small_signal_reactances`], not as well — whatever both report
    /// gets stamped, so a device answering to both double-counts. Native
    /// devices report two-terminal branches and are fine as they are; this
    /// exists for OSDI/Verilog-A, where ∂q_i/∂v_j need not equal ∂q_j/∂v_i
    /// (transcapacitance) and so cannot be expressed as reciprocal branches at
    /// all.
    fn load_reactive_jacobian(&self, _c_mat: &mut [crate::mna::SparseRow]) {}

    /// Complex admittance entries at one frequency that are **not** of the form
    /// `G + jωC + Λ/(jω)`.
    ///
    /// The `G`/`C`/`Λ` matrices cover every device whose small-signal behaviour
    /// is a rational function of `jω` with the poles the assembly already knows
    /// about. A delay is not one: its frequency response is `exp(−jωτ)`, which
    /// is transcendental and cannot be spelled in those three matrices at all.
    ///
    /// Such a device stamps the delayed coupling here instead of writing it to
    /// the residual, which is where the transient keeps it and where no
    /// frequency-domain analysis reads it (#110). Returning a non-empty list is
    /// also the declaration that this device has **no linear matrix pencil**,
    /// which is why `.pz` refuses a circuit containing one rather than solving a
    /// system it has silently truncated.
    ///
    /// Indices are MNA rows/columns, as for [`Device::load_jacobian`], and the
    /// entries are *added* to whatever `load_jacobian` already stamped.
    fn ac_stamps(&self, _omega: f64) -> Vec<AcStamp> {
        Vec::new()
    }

    /// Resolve the matrix cells this device stamps into, once, against the
    /// structural pattern the hot loop's matrix was built with.
    ///
    /// Finding a cell is a `binary_search` per stamp, per device, per Newton
    /// iteration, per timestep. The pattern never moves, so the answer is
    /// settled here instead of being re-derived hundreds of millions of times
    /// (issue #24). A device that implements this stores
    /// [`Pattern::id`](crate::mna::Pattern::id) alongside the cells and
    /// compares it to [`MnaMatrix::pattern_id`](crate::mna::MnaMatrix::pattern_id)
    /// before using them.
    ///
    /// Default is a no-op, and **not implementing it is always safe** — the
    /// device keeps addressing cells by `a[r][c]`, which is slower and
    /// identical. So is never calling it: a device with no resolved cells
    /// takes the same searching path. That asymmetry is deliberate. A stamp
    /// through a stale slot is a silent wrong answer, so every way of getting
    /// this wrong has to land on the slow path rather than the wrong cell.
    fn resolve_cells(&mut self, _pattern: &crate::mna::Pattern) {}

    /// Set a named real-valued parameter on this device instance.
    ///
    /// Returns `true` if the parameter was found and set; `false` if not supported
    /// or if `name` is not a parameter of this device (built-in devices return false
    /// by default — their parameters are set at construction time).
    fn set_real_param(&mut self, _name: &str, _value: f64) -> bool {
        false
    }

    /// Last word on whether this device's configuration makes sense, called
    /// once at construction after the terminals are bound *and* every parameter
    /// has been applied.
    ///
    /// This is the only point in a device's life at which both halves of its
    /// configuration are present: `setup_instance` knows the terminals but not
    /// the parameters, and `set_real_param` knows one parameter but not the
    /// geometry. Anything that needs both — a power budget that must sum to one,
    /// a transfer matrix that must be passive, a port count that follows from a
    /// parameter — can only be judged here.
    ///
    /// Devices used to reach for one of two workarounds instead, and both are
    /// worse than an error:
    ///
    /// - **Defer to the first `eval`** behind a "have I been checked yet" flag.
    ///   `eval` cannot return an error, so the check had to `assert!`, which
    ///   means a mis-typed power budget arrived as a panic from inside the
    ///   solve — through `pyo3`, for a Python caller — rather than as a
    ///   diagnostic naming the element.
    /// - **Assert inside `setup_instance`**, which fires before the parameter
    ///   that would have made the configuration legal has been applied.
    ///
    /// The `Err` string says what is wrong with the device; the caller
    /// (`build_devices`) prefixes the element and model name, so an
    /// implementation should not repeat them.
    ///
    /// Returning `Ok` is always safe: a device with nothing to check does not
    /// implement this.
    /// The largest timestep this device will accept for the *next* step, if it
    /// has an opinion.
    ///
    /// Defaults to `None`. Native models have no way to ask: their reactive
    /// branches are declared, so the integrator already knows their time
    /// constants and the LTE controller handles the rest.
    ///
    /// It exists for compiled models, where Verilog-A's `$bound_step` is the
    /// model saying "do not step past this or you will miss something". LTE alone
    /// cannot cover that, because it measures the error of a step already taken.
    fn requested_max_timestep(&self) -> Option<f64> {
        None
    }

    /// Whether the *last* `eval` produced a result the model itself considers
    /// valid.
    ///
    /// Defaults to `Ok(())`, which is right for every native model: they evaluate
    /// a closed form and have no way to disown the answer.
    ///
    /// It exists for compiled models. OSDI's `eval` returns a flag word, and
    /// `EVAL_RET_FLAG_FATAL` means "this bias point is not one I can evaluate".
    /// Stamping the result anyway invents an answer out of numbers the model
    /// disowned. The Newton loop treats a device reporting `Err` here the way it
    /// treats a clamped step: this iterate cannot be the converged one, and if a
    /// device is still saying so when the iteration budget runs out, the message
    /// becomes the diagnosis rather than a bare non-convergence.
    fn eval_status(&self) -> Result<(), String> {
        Ok(())
    }

    fn validate(&mut self) -> Result<(), String> {
        Ok(())
    }

    /// Scale this device's INDEPENDENT source output for source-stepping
    /// homotopy, where `scale` ramps 0 → 1.
    ///
    /// Netlist `V`/`I` sources are ramped by the solver directly (see
    /// `mna::stamp_netlist_scaled_in_place`), but a device that is itself an
    /// independent source — an optical laser, most importantly — is invisible
    /// to that machinery, so without this hook it stays at full output for the
    /// whole ramp. In a photonic circuit the optical power IS the dominant
    /// excitation, which left source-stepping with no homotopy path at all:
    /// every trial point sat in the fully-illuminated nonlinearity.
    ///
    /// Default is a no-op, which is correct for every device that is not an
    /// independent source.
    fn set_source_scale(&mut self, _scale: f64) {}

    /// Internal noise current sources for `.noise` analysis.
    ///
    /// Called AFTER a DC-OP `eval()` so the device can read its cached
    /// bias-point quantities (e.g. Id, gm) when building the PSDs.
    ///
    /// Returns one or more `(pos, neg, S_i)` tuples where `S_i` is the
    /// one-sided current PSD in A²/Hz of an uncorrelated noise source
    /// connected between those two terminal node indices.  Default is empty
    /// (resistor thermal noise is iterated as `Element::Resistor` in
    /// `noise_analysis`, not through this hook).
    ///
    /// `freq` is the analysis frequency in Hz.  Every native device here is
    /// flat and ignores it — shot noise and RIN do not care — but an OSDI model
    /// may call `flicker_noise()`, whose density is a function of frequency, so
    /// the argument has to reach the hook.  `.noise` passes the sweep point.
    /// Transient noise realises a *white* sample sequence and cannot represent
    /// a sloped density at all, so it probes mid-band and warns if the source
    /// turns out to vary; see [`crate::noise::TransientNoise`].
    fn noise_sources(&self, _ctx: &SimContext, _freq: f64) -> Vec<(NodeId, NodeId, f64)> {
        Vec::new()
    }

    /// Noise generators whose one random process reaches the circuit at more
    /// than one place at once — see [`CorrelatedNoise`].  Default is empty.
    fn correlated_noise_sources(&self, _ctx: &SimContext) -> Vec<CorrelatedNoise> {
        Vec::new()
    }

    /// Device node indices whose potential is a **temperature in kelvin**
    /// rather than a voltage, so the solver can bound them with `temptol`
    /// instead of `vntol`.
    ///
    /// Numbered as the device's own nodes are: `0..terminals.len()` are the
    /// terminals it was handed, and the indices above that are the extra rows
    /// `num_extra_nodes` asked for, in the same order. `push_device` translates
    /// both into MNA rows, which is why one list covers a self-heating internal
    /// node and a shared thermal port alike.
    ///
    /// A Verilog-A model declares this by declaring its node `thermal`, and
    /// nothing else has to agree: OSDI carries the discipline's units through to
    /// the descriptor, so this is read off the model rather than restated in the
    /// deck. Default empty — a device with no temperature unknown says nothing.
    ///
    /// # What a thermal node is not
    ///
    /// It is not a different kind of row. The potential is kelvin and the flow
    /// is watts, so KCL over that node is conservation of power, and a SPICE
    /// `R h1 h2 4.2e4` across two of them stamps `ΔT/R` watts — a thermal
    /// resistance of 42 kK/W, which is exactly right and is how an electrical
    /// element earns its place in a thermal network. That is why there is no
    /// discipline check here refusing R/C/I on a thermal net, as there is for
    /// optical wires: the electrical primitives *are* the thermal primitives,
    /// and refusing them would ban thermal crosstalk between two devices.
    fn thermal_nodes(&self) -> Vec<usize> {
        Vec::new()
    }

    /// How many extra MNA rows this device needs beyond its terminal nodes.
    ///
    /// Verilog-A models with potential contributions (`V(port) <+ value;`)
    /// implicitly introduce branch-flow nodes that OSDI exposes as
    /// `num_nodes − num_terminals`.  The topology allocates that many rows
    /// per device and then calls `bind_extra_nodes` with the starting index.
    /// Returns 0 by default — native Rust devices that stamp directly into
    /// terminal nodes have no extras.
    fn num_extra_nodes(&self) -> usize {
        0
    }

    /// Tell the device which contiguous MNA rows (starting at `first_idx`)
    /// have been allocated for its internal nodes.  Called once after
    /// construction and before the first `eval()`.  Default is a no-op.
    fn bind_extra_nodes(&mut self, _first_idx: usize) {}

    /// MNA rows/columns this device stamps into *beyond* its own terminals and
    /// the extra rows bound by `bind_extra_nodes`.
    ///
    /// The sparsity pattern (`mna::Pattern`) is derived structurally from every
    /// device's terminals plus extras, on the invariant that a device only
    /// couples nodes it was handed.  Devices that reach further — a behavioural
    /// source whose expression references arbitrary `V(node)` / `I(vsrc)` —
    /// must report those rows here, or their entries get silently dropped from
    /// the sparse solve.  Debug builds catch a violation via
    /// `MnaMatrix::debug_assert_covers`.  Default is empty.
    fn extra_stamp_rows(&self) -> Vec<usize> {
        Vec::new()
    }

    /// The exact `(row, col)` cells this device stamps, when it knows them.
    ///
    /// By default `mna::Pattern` takes a **clique** over every row a device
    /// owns, because a stamper handed a set of rows may couple any of them to
    /// any other. That is the tightest footprint derivable without asking, and
    /// it is fine for a two-terminal device — but it costs `O(rows²)` and some
    /// photonic devices own a great many rows for a structurally sparse stamp.
    /// An N×N `fc_awgr` owns `9N²` rows: at N = 8 the clique is 332 k cells
    /// standing in for a true footprint of 6 k, and it grows as `N⁴` where the
    /// device is `N³`.
    ///
    /// Return `Some(pairs)` to declare the footprint exactly; anything stamped
    /// outside it is silently dropped by the sparse solve, so a device that
    /// implements this must keep it in sync with its stamping (debug builds
    /// catch a violation via `MnaMatrix::debug_assert_covers`). `None` — the
    /// default — keeps the conservative clique.
    fn stamp_pairs(&self) -> Option<Vec<(usize, usize)>> {
        None
    }

    /// MNA columns whose `∂f/∂x` this device deliberately does **not** stamp.
    ///
    /// A device may linearise about a coefficient it freezes at the previous
    /// Newton iterate rather than differentiating it.  Successive substitution
    /// reaches the same fixed point, so the converged answer is identical — and
    /// for a term like `∂φ/∂λ ≈ 2.7e9 rad/m` freezing is the only thing that
    /// converges at all.  Newton is content with any iteration matrix that
    /// contracts; it does not need the exact derivative.
    ///
    /// The adjoint method does.  `dL/dp = −λᵀ·∂f/∂p` is the total derivative
    /// only if `Jᵀλ = ∂L/∂x` used the true `J`, so a frozen block does not make
    /// the gradient approximate — it makes every path through that block
    /// silently contribute **zero**.  For an electro-optic device that is the
    /// path the user cares about most.
    ///
    /// Report those columns here.  `crate::adjoint` re-derives them numerically
    /// for the adjoint solve alone and leaves the Newton iteration matrix
    /// untouched, so declaring a column costs nothing at solve time.
    ///
    /// λ wires do not need declaring, and no longer *can* be: a wavelength is
    /// resolved before the solve rather than solved for, so there is no λ
    /// column to freeze.  Declare the *electrical* columns an optical device
    /// reads: drive voltages, control wires, thermal nodes.
    ///
    /// `crate::adjoint::jacobian_check` is the oracle for this: any mismatch in
    /// an undeclared column is a missing declaration or a wrong stamp.
    fn frozen_jacobian_columns(&self) -> Vec<usize> {
        Vec::new()
    }

    /// How a wavelength label moves through this device, as
    /// `(from_terminal, to_terminal)` pairs in the device's own terminal
    /// numbering.
    ///
    /// A wavelength is a label, not a state: measured across every photonic deck
    /// in the tree, each λ row read exactly a source's wavelength and never
    /// anything computed (`tests/lambda_is_a_label.rs`). Resolving it before the
    /// solve is what stopped it being an MNA unknown — which took 864 of 2840
    /// rows off the giona RNN and deleted `LambdaSelect`'s latch, the
    /// `lambdatol` tolerance class and λ's trust-region exemption with them.
    ///
    /// It has to be *declared* rather than read off the assembled matrix,
    /// because the matrix is exactly what is going away. This is the same shape
    /// as bundle arity: knowledge that used to be inferred from a structure the
    /// device happens to build, moved next to the device that knows it.
    ///
    /// Default is empty — correct for anything with no optical ports, and for a
    /// terminator like a photodetector that ends the path.
    fn lambda_routing(&self) -> Vec<(usize, usize)> {
        Vec::new()
    }

    /// The resolved wavelength at each of this device's terminals, in metres,
    /// in terminal order.
    ///
    /// Called once per build, after `setup_instance` and before the first
    /// `eval`. λ is resolved before the solve (`crate::lambda`), so a device
    /// that needs a channel's wavelength — to evaluate a propagation phase, a
    /// filter passband, a router's grid — reads it from here instead of from a
    /// matrix row. Non-λ terminals carry the band centre; a device indexes the
    /// λ positions of its own layout and ignores the rest.
    ///
    /// Default is a no-op: correct for every device whose physics does not
    /// depend on wavelength, which is all of the electrical ones.
    fn set_resolved_lambda(&mut self, _per_terminal: &[f64]) {}

    /// Every terminal of this device that carries a wavelength label.
    ///
    /// This is what makes the λ net set *total*: a label's value is decided
    /// before the solve, so every wire that carries one has to be enumerable
    /// without looking at a matrix. It used to be read off net names — anything
    /// ending `_wl_<k>` that a `.optical_port` had declared — which cannot see a
    /// PCell that hand-wires its bundle (`examples/photonic/pcells/source_bank.sp`
    /// names them `a1w`), and every real deck in this tree does exactly that.
    ///
    /// The default derives it from the two declarations that already exist, so a
    /// device that routes or emits a label says nothing extra. Override only
    /// when a device *reads* a label on a terminal it neither routes from nor
    /// emits at: `fc_awgr` evaluates each input port's passband at that port's
    /// own λ while mirroring only one port's tag onto its outputs, so the other
    /// ports' λ terminals appear nowhere in its routing.
    fn lambda_terminals(&self) -> Vec<usize> {
        let mut t: Vec<usize> = self
            .lambda_routing()
            .into_iter()
            .flat_map(|(a, b)| [a, b])
            .chain(self.lambda_emitted().into_iter().map(|(t, _)| t))
            .collect();
        t.sort_unstable();
        t.dedup();
        t
    }

    /// Wavelengths this device *originates*, as `(terminal, λ in metres)`.
    ///
    /// A source is where resolution starts. Native emitters already hold this as
    /// a parameter and merely deliver it through the matrix, so exposing it costs
    /// nothing; a model that writes its λ wire in `analog` and declares nothing
    /// cannot be resolved, and is worth a diagnostic rather than a silent zero.
    fn lambda_emitted(&self) -> Vec<(usize, f64)> {
        Vec::new()
    }
}
