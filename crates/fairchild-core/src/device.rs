use crate::mna::MnaMatrix;

pub const K_BOLTZMANN: f64 = 1.380649e-23;
pub const Q_ELECTRON: f64 = 1.602176634e-19;

/// Index of a terminal in the MNA solution vector; `None` → ground (excluded from matrix).
pub type NodeId = Option<usize>;

/// Simulator context passed to device model callbacks at every eval.
pub struct SimContext {
    pub temperature: f64, // Kelvin; default 300.15 K (27 °C, SPICE TNOM)
    pub omega_0: f64,     // rad/s carrier frequency for optical ports; 0 for electrical
    /// When true, device models apply junction-step limiters (pnjlim, fetlim).
    /// Mapped from `SimOptions::pnjlim` / `.options nopnjlim`.
    pub jlim_enabled: bool,
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
}

impl Default for SimContext {
    fn default() -> Self {
        SimContext {
            temperature: 300.15,
            omega_0: 0.0,
            jlim_enabled: true,
            lambda_center_m: 1.55e-6,
            bidirectional_propagation: false,
            waveguide_delay: false,
            time_s: 0.0,
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
#[derive(Clone, Copy, Debug)]
pub struct ReactiveBranchSpec {
    pub kind: ReactiveKind,
    pub pos: NodeId,
    pub neg: NodeId,
    /// Current capacitance (F) or inductance (H) at the device's cached
    /// operating point.  Re-queried per NR iteration.
    pub value: f64,
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
    /// timestep the integrator reads V_C = x[pos] − x[neg] from the
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
    fn load_reactive_jacobian(&self, _c_mat: &mut [Vec<f64>]) {}

    /// Set a named real-valued parameter on this device instance.
    ///
    /// Returns `true` if the parameter was found and set; `false` if not supported
    /// or if `name` is not a parameter of this device (built-in devices return false
    /// by default — their parameters are set at construction time).
    fn set_real_param(&mut self, _name: &str, _value: f64) -> bool {
        false
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
    fn noise_sources(&self, _ctx: &SimContext) -> Vec<(NodeId, NodeId, f64)> {
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
}
