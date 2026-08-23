//! Host-driven transient stepping.
//!
//! [`TranStepper`] is the fixed-step transient loop turned inside out: instead
//! of running to `stop` and handing back a [`TranResult`], it holds the
//! integrator state between timesteps so an external program can drive the
//! clock, read node voltages, and write source values at every point.  That is
//! what mixed-signal co-simulation needs — the digital side advances the analog
//! side one step, samples it, and drives it back.
//!
//! [`crate::tran::tran_nr_with_registry_opts`] is implemented on top of this,
//! so there is exactly one fixed-step integrator in the tree and the batch
//! goldens cover the stepping path too.
//!
//! ```no_run
//! # use fairchild_core::{DeviceRegistry, SimOptions, TranStepper};
//! # fn demo(netlist: fairchild_parser::Netlist) -> Result<(), fairchild_core::SimError> {
//! let mut registry = DeviceRegistry::new();
//! registry.register_builtin_models(&netlist.models);
//! let opts = SimOptions::from_netlist(&netlist);
//! let mut st = TranStepper::new(netlist, &registry, &opts, 1e-9)?;
//! while st.time() < 1e-6 {
//!     let comparator = st.node("out")? > 0.9;          // analog -> digital
//!     st.set_source("vdrive", if comparator { 0.0 } else { 1.8 })?; // digital -> analog
//!     st.step()?;
//! }
//! # Ok(())
//! # }
//! ```

use std::collections::HashMap;

use fairchild_parser::{Element, Netlist, Waveform};

use crate::device::{Device, EvalFlags, SimContext};
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::mna::{CircuitTopology, MnaMatrix, RowFloor, StampPlan};
use crate::newton::{build_devices_with_footprints, dc_op_nr_with_registry_opts};
use crate::options::SimOptions;
use crate::reactive::{stamp_device_branches, ReactiveState};
use crate::solver::{Factorisation, LinearSolver};
use crate::tran::IntegratorMode;

/// A transient analysis paused between timesteps.
///
/// The step size is fixed for the lifetime of the stepper.
///
// ponytail: fixed h only — but no longer for any deep reason.  Reactive history
// is physical now (`crate::reactive`), so companions rebuild correctly for any
// `h`; landing on an arbitrary externally-chosen event time needs only a
// `step_to(h)` that varies it, plus a decision about what the LTE controller
// should do with a host-imposed step.  Add it when a digital side actually
// needs off-grid event times; a clock-locked host does not.
pub struct TranStepper {
    /// Owned so `set_source` can rewrite source waveforms between steps.
    netlist: Netlist,
    opts: SimOptions,
    ctx: SimContext,
    mode: IntegratorMode,
    topo: CircuitTopology,
    x: Vec<f64>,
    devices: Vec<Box<dyn Device>>,
    plan: StampPlan,
    solver: Box<dyn LinearSolver>,
    fact: Option<Box<dyn Factorisation>>,
    mat: MnaMatrix,
    /// Every reactive branch's history plus the companions derived from it —
    /// netlist C/L and device-declared alike, one representation, shared with
    /// the variable-step integrator.
    reactive: ReactiveState,
    step: f64,
    t: f64,
    /// Trapezoidal deliberately takes its first step with Backward Euler, for
    /// stability across the t=0 discontinuity.  Step control, not integration
    /// method, so it lives here rather than in `ReactiveState`.
    first_step: bool,
    /// Lower-cased V/I source name → index into `netlist.elements`, so
    /// `set_source` is a lookup rather than a scan of every element.
    sources: HashMap<String, usize>,
    /// Per-row Newton step tolerances — not every unknown is a volt.  Built
    /// once here because the netlist and topology cannot change under a
    /// stepper; `set_source` only rewrites a waveform.  See `crate::tolerance`.
    tol: crate::tolerance::Tolerances,
    /// `None` unless `.options trannoise=1`; see `crate::noise::TransientNoise`.
    noise: Option<crate::noise::TransientNoise>,
}

impl TranStepper {
    /// Build the t = 0 state: DC operating point (or `.ic` under UIC), devices,
    /// sparsity pattern, linear solver, and the reactive companion models.
    ///
    /// `step` is clamped to `opts.max_step`, matching the batch path.
    pub fn new(
        netlist: Netlist,
        registry: &DeviceRegistry,
        opts: &SimOptions,
        step: f64,
    ) -> Result<Self, SimError> {
        Self::build(netlist, registry, opts, step, None)
    }

    /// [`TranStepper::new`] with `t = 0` supplied rather than solved for.
    ///
    /// The adjoint's parameter replay needs a stepper whose trajectory is
    /// *imposed*: it re-walks an already-solved run with one parameter nudged,
    /// holding every `x_k` frozen, so the DC operating point that `new` would
    /// compute is both wasted work and the wrong point to start from — the
    /// whole construction differentiates the residual at a fixed `x`.
    pub(crate) fn seeded(
        netlist: Netlist,
        registry: &DeviceRegistry,
        opts: &SimOptions,
        step: f64,
        x0: &[f64],
    ) -> Result<Self, SimError> {
        Self::build(netlist, registry, opts, step, Some(x0))
    }

    fn build(
        netlist: Netlist,
        registry: &DeviceRegistry,
        opts: &SimOptions,
        step: f64,
        seed: Option<&[f64]>,
    ) -> Result<Self, SimError> {
        if opts.sanity_check && opts.uic {
            crate::sanity::check_netlist_sanity(&netlist);
        }
        crate::connectivity::check_connectivity(&netlist)?;
        let ctx = opts.sim_context();
        let mode = opts.method;

        // With UIC: skip DC OP, seed x from `.ic` (or 0 where unspecified).
        // Without UIC (the default): use DC OP as t=0 condition.
        let (mut topo, mut x) = if let Some(x0) = seed {
            // Reproduce the row layout the operating-point path arrives at.
            // `push_device` appends fresh rows on every call, so building the
            // devices twice — once for the DC solve, once here — is what fixes
            // where each device's internal nodes live.  The seeded path skips
            // the solve but must not skip the allocation, or it would be
            // indexing different unknowns than the run it is replaying.
            //
            // ponytail: those first rows are then dead weight in every
            // non-UIC transient — 6 of 21 on a modest photonic netlist, and the
            // DC operating point's internal-node values are discarded with
            // them.  Worth removing, but that moves every matrix index in the
            // tree and belongs in its own change.
            let mut t = CircuitTopology::build_resolved(&netlist, &ctx, registry);
            let _ = build_devices_with_footprints(&netlist, &mut t, &ctx, registry)?;
            (t, x0.to_vec())
        } else if opts.uic {
            let topo = CircuitTopology::build_resolved(&netlist, &ctx, registry);
            let mut x = vec![0.0f64; topo.size];
            for (name, value) in &netlist.ic {
                if let Some(&i) = topo.node_index.get(name) {
                    x[i] = *value;
                }
            }
            (topo, x)
        } else {
            let dc = dc_op_nr_with_registry_opts(&netlist, registry, opts)?;
            // DC OP already allocated extras; reuse its topology so the row
            // layout (and matrix size) stays consistent through the transient.
            (dc.topo, dc.x)
        };

        let (mut devices, footprints) =
            build_devices_with_footprints(&netlist, &mut topo, &ctx, registry)?;
        // After build_devices, topo.size is final.  Pad the initial x vector
        // to match — when we came via UIC, x was sized to the pre-allocation
        // topology; OSDI internal nodes get a zero initial guess.
        x.resize(topo.size, 0.0);
        // Topology is fixed across the whole transient — build the linear
        // solver and the structural sparsity pattern once.
        let solver = opts.linear_solver(topo.size);
        let plan = StampPlan::new(&topo, &netlist, &footprints);
        plan.resolve_device_cells(&mut devices);

        // Seed x_tprev from DC OP (or UIC initial conditions) so reactive
        // history is defined before the first step.
        for dev in &mut devices {
            dev.commit_timestep(&x);
        }

        // Honour opts.max_step as an upper bound on the step size.
        let step = step.min(opts.max_step);

        let reactive = ReactiveState::new(&netlist, &topo, &mut devices, &ctx, &x);
        let mat = MnaMatrix::with_pattern(topo.size, plan.pattern.clone());

        // The parser already lower-cases element names, so this is keyed as-is.
        let sources: HashMap<String, usize> = netlist
            .elements
            .iter()
            .enumerate()
            .filter_map(|(i, el)| match el {
                Element::VoltageSource { name, .. } | Element::CurrentSource { name, .. } => {
                    Some((name.to_lowercase(), i))
                }
                _ => None,
            })
            .collect();

        let tol = crate::tolerance::Tolerances::build(&topo, opts);
        let noise = crate::noise::TransientNoise::new(&netlist, &topo, opts);

        Ok(TranStepper {
            netlist,
            opts: opts.clone(),
            ctx,
            mode,
            topo,
            x,
            devices,
            plan,
            solver,
            fact: None,
            mat,
            reactive,
            step,
            t: 0.0,
            first_step: true,
            sources,
            tol,
            noise,
        })
    }

    // ── stepping ──────────────────────────────────────────────────────────

    /// Advance one timestep.  Returns the new simulation time.
    ///
    /// On `Err(NoConvergence)` the stepper is left on the previous accepted
    /// timepoint — the failed iterate is not committed, so a host can lower a
    /// drive level and retry.
    pub fn step(&mut self) -> Result<f64, SimError> {
        let t_next = self.t + self.step;
        self.solve_at(t_next)?;
        self.commit(t_next);
        self.advance_history();
        Ok(self.t)
    }

    /// Step until the simulation time reaches `t_target`.
    ///
    /// The step size is fixed, so the final timepoint lands on the first grid
    /// point at or past `t_target`; the returned time says where that was.
    /// Already at or past `t_target` ⇒ no steps taken.
    pub fn advance_to(&mut self, t_target: f64) -> Result<f64, SimError> {
        // Tolerance against float drift in `t += h` accumulation, so a target
        // that is an exact multiple of h doesn't buy one extra step.
        let eps = self.step * 1e-9;
        while self.t < t_target - eps {
            self.step()?;
        }
        Ok(self.t)
    }

    /// One Newton-Raphson solve for the timepoint `t_next`, leaving the result
    /// in `self.x`.  Does not commit history — see [`Self::commit`].
    ///
    /// `alpha` comes from the configured step, not from `t_next - self.t`: the
    /// companion models were built for `self.step` and the batch driver's final
    /// clamped step relies on exactly this behaviour.
    pub(crate) fn solve_at(&mut self, t_next: f64) -> Result<(), SimError> {
        self.prepare(t_next, self.step_mode(), self.step);
        // One noise realisation for this step, drawn at the previous accepted
        // bias and then held: shot noise follows the current, and the current
        // is whatever the circuit last settled at.  Devices are re-evaluated at
        // `self.x` first so the very first step (where UIC leaves them fresh)
        // is not drawn from a zero bias.  Never redrawn inside the Newton loop
        // — see `TransientNoise::draw`.  Drawn here rather than in `prepare`,
        // which the adjoint calls a second time at a perturbed `h`: that must
        // re-stamp the same realisation, not a fresh one.
        if let Some(noise) = self.noise.as_mut() {
            for dev in &mut self.devices {
                dev.eval(&self.x, EvalFlags::tran(), &self.ctx);
            }
            noise.draw(&self.devices, &self.ctx, self.step);
        }
        // Solve into a scratch vector so a non-converging step leaves the last
        // accepted solution intact for the caller to inspect or retry from.
        let mut x_try = self.x.clone();

        for _iter in 0..self.opts.itl4 {
            self.stamp_at(&x_try);
            // Added here rather than inside `stamp_at` so the adjoint's
            // re-stamps stay a pure function of the frozen state: the held
            // realisation belongs to this timepoint's forward solve only.
            if let Some(noise) = self.noise.as_ref() {
                for (b, n) in self.mat.b.iter_mut().zip(noise.rhs()) {
                    *b += n;
                }
            }

            let x_new = if let Some(f) = self.fact.as_mut() {
                f.refactor_and_solve_mat(&self.mat)?
            } else {
                let mut f = self.solver.factorise_mat(&self.mat)?;
                let r = f.refactor_and_solve_mat(&self.mat)?;
                self.fact = Some(f);
                r
            };

            let max_dv = x_new
                .iter()
                .zip(x_try.iter())
                .take(self.topo.n_nodes())
                .map(|(n, o)| (n - o).abs())
                .fold(0.0f64, f64::max);

            let x_next: Vec<f64> = if max_dv > self.opts.vmax {
                let scale = self.opts.vmax / max_dv;
                x_try
                    .iter()
                    .zip(x_new.iter())
                    .map(|(o, n)| o + scale * (n - o))
                    .collect()
            } else {
                x_new
            };

            let converged = self.tol.converged(&x_next, &x_try);

            x_try = x_next;
            if converged {
                self.x = x_try;
                return Ok(());
            }
        }

        Err(SimError::NoConvergence {
            iters: self.opts.itl4,
        })
    }

    /// Everything one timestep needs that does *not* change per Newton
    /// iteration: the companion models for `h` under `mode`, and the absolute
    /// time devices resolve their sources and delays against.
    ///
    /// Split out from [`Self::stamp_at`] because the adjoint re-stamps the same
    /// timepoint at a second step size to read the charge Jacobian out of the
    /// difference, and because rebuilding companions per Newton iteration would
    /// be pure waste.  Reactive history is physical, so calling this again with
    /// the step's own `h` restores exactly what the last call left behind.
    pub(crate) fn prepare(&mut self, t: f64, mode: IntegratorMode, h: f64) {
        self.ctx.time_s = t;
        // Companions for this step, derived from physical history rather than
        // advanced in place.  One rebuild per step buys a single representation
        // that any step size can consume — see `crate::reactive`.
        self.reactive.build(&self.devices, mode, h, None);
        // Devices that stamp their own reactance (OSDI/Verilog-A `ddt`) need the
        // method, not just `alpha`, which can only express Backward Euler. Same
        // `mode` and the same absent BDF-2 history as the branch stamper in
        // `stamp_at`, so one decision still reaches everything.
        self.ctx.discretisation = Some(crate::device::Discretisation {
            mode,
            h,
            gear2_h_prev: None,
        });
    }

    /// Stamp the whole transient system at `probe`, leaving `A` and `b` in
    /// `self.mat` — one Newton iteration with the solve taken out.
    ///
    /// [`Self::prepare`] must have run for this timepoint.
    pub(crate) fn stamp_at(&mut self, probe: &[f64]) {
        let d = self
            .ctx
            .discretisation
            .expect("stamp_at without a preceding prepare");
        let alpha = 1.0 / d.h;
        let t = self.ctx.time_s;
        crate::mna::stamp_netlist_in_place(
            &mut self.mat,
            &self.topo,
            &self.netlist,
            t,
            &self.reactive.cap_state,
            &self.reactive.ind_state,
            Some(&self.plan),
            crate::mna::InductorDc::Short,
        );

        for dev in &mut self.devices {
            dev.set_source_scale(1.0);
            dev.eval(probe, EvalFlags::tran(), &self.ctx);
            dev.load_residual_tran(&mut self.mat.b, alpha);
            dev.load_jacobian_tran(&mut self.mat, alpha);
        }
        // Stamp integrator-managed reactive companions for every
        // device-declared linear reactive branch (uses the device's current
        // bias-dependent value AND the history from the previous accepted
        // timestep).
        // BDF-2 needs two-timepoint history, which only the variable-step
        // integrator carries; GEAR demotes to BE here, as it always has for the
        // built-in C/L on this path too.
        stamp_device_branches(
            &self.devices,
            &self.reactive.dev_state,
            &mut self.mat,
            probe,
            d.h,
            d.mode,
            None,
        );

        self.topo
            .stamp_gmin(&mut self.mat.a, self.opts.gmin, RowFloor::PinEmptyRows);
    }

    /// Stamp the *DC* system at `probe`, the way the `t = 0` operating point was
    /// solved.  The adjoint needs it because the initial condition is a
    /// constraint like any other timestep, and its Jacobian is the DC one.
    pub(crate) fn stamp_dc_at(&mut self, probe: &[f64]) {
        self.ctx.discretisation = None;
        let _ = crate::newton::residual_l2(
            &mut self.mat,
            &self.topo,
            &self.netlist,
            &mut self.devices,
            &self.ctx,
            &self.opts,
            1.0,
            0.0,
            Some(&self.plan),
            probe,
        );
    }

    /// The system as last stamped.
    pub(crate) fn matrix(&self) -> &MnaMatrix {
        &self.mat
    }

    pub(crate) fn devices(&self) -> &[Box<dyn Device>] {
        &self.devices
    }

    /// The conductance each device-declared reactive branch is currently
    /// stamping, labelled by its owning element.
    ///
    /// Only `crate::adjoint_tran::charge_lag` wants this: it rebuilds the
    /// companion after evaluating at the step's own solution and compares, which
    /// is how you see a bias-dependent capacitance being stamped one `eval`
    /// stale.
    pub(crate) fn device_branch_conductances(&self) -> Vec<(String, f64)> {
        let names = crate::newton::device_element_names(&self.netlist);
        let mut out = Vec::new();
        for (d, dev) in self.devices.iter().enumerate() {
            let n_br = dev.reactive_branches().len();
            for b in 0..n_br {
                let label = names
                    .get(d)
                    .map_or_else(|| format!("device[{d}]"), |n| n.clone());
                let g = self.reactive.dev_state[d][b].0;
                out.push((format!("{label}#{b}"), g));
            }
        }
        out
    }

    /// Impose a solution rather than solving for one — the replay's whole point.
    pub(crate) fn force_solution(&mut self, x: &[f64]) {
        self.x.copy_from_slice(x);
    }

    /// Retune a live device, then re-seed the reactive history from it.
    ///
    /// The re-seed is the part that is easy to forget: `set_real_param` can move
    /// a device's declared capacitance, and the `t = 0` history was seeded from
    /// the old one.  Leaving it stale would make the gradient with respect to
    /// that parameter miss its entire history path.
    ///
    /// This used to re-`eval` every device by hand first, because a
    /// bias-dependent capacitance reports whatever its last eval cached and
    /// seeding from a stale one is silently wrong.  `18b5744` hit the same
    /// hazard from the other direction — a transient starting from a DC point
    /// its devices had never been evaluated at — and fixed it inside
    /// `ReactiveState::new`, which now takes `&mut` and evaluates for itself.
    /// One place interprets it, so the loop here is gone.
    ///
    /// Returns whether the device recognised the parameter.
    pub(crate) fn set_device_param(&mut self, i: usize, name: &str, value: f64) -> bool {
        let ok = self.devices[i].set_real_param(name, value);
        self.ctx.discretisation = None;
        self.reactive = ReactiveState::new(
            &self.netlist,
            &self.topo,
            &mut self.devices,
            &self.ctx,
            &self.x,
        );
        ok
    }

    /// Accept the solved iterate as the state at `t_next`.
    pub(crate) fn commit(&mut self, t_next: f64) {
        for dev in &mut self.devices {
            dev.commit_timestep(&self.x);
        }
        self.t = t_next;
    }

    /// Roll the reactive companion history forward so the next step integrates
    /// from the just-accepted solution.
    pub(crate) fn advance_history(&mut self) {
        self.reactive.accept(&self.devices, &self.x);
        self.first_step = false;
    }

    /// The method this step actually integrates with.
    ///
    /// Trapezoidal takes its first step with Backward Euler, for stability
    /// across the t=0 discontinuity; `accept` then records the capacitor current
    /// that step implied, which is exactly what TR needs to continue.
    pub(crate) fn step_mode(&self) -> IntegratorMode {
        if self.first_step && matches!(self.mode, IntegratorMode::Trapezoidal) {
            IntegratorMode::BackwardEuler
        } else {
            self.mode
        }
    }

    // ── reading the analog state ──────────────────────────────────────────

    /// Current simulation time, in seconds.
    pub fn time(&self) -> f64 {
        self.t
    }

    /// The fixed step size, after clamping to `opts.max_step`.
    pub fn step_size(&self) -> f64 {
        self.step
    }

    /// Voltage at `node` for the current timepoint.
    pub fn node(&self, node: &str) -> Result<f64, SimError> {
        self.topo.node_voltage(node, &self.x)
    }

    /// Current through voltage source `name` for the current timepoint.
    pub fn vsrc_current(&self, name: &str) -> Result<f64, SimError> {
        self.topo.vsrc_current(name, &self.x)
    }

    /// Every solvable node name, in MNA row order.
    pub fn node_names(&self) -> impl Iterator<Item = &str> {
        self.topo.node_index.keys().map(|s| s.as_str())
    }

    /// The raw MNA solution vector for the current timepoint.
    pub fn solution(&self) -> &[f64] {
        &self.x
    }

    pub(crate) fn topology(&self) -> &CircuitTopology {
        &self.topo
    }

    // ── driving the analog state ──────────────────────────────────────────

    /// Hold voltage or current source `name` at `value` from the next step on.
    ///
    /// This is a zero-order hold: the source keeps `value` until the next
    /// `set_source`, which is the right semantics for a digital driver and
    /// matches what ngspice's `GetVSRCData` callback provides.
    ///
    /// A change in value is a slope discontinuity the integrator cannot see
    /// coming, exactly like a `PULSE` edge — with a fixed step there is no
    /// predictor to protect, but the step immediately after a change resolves
    /// the edge no better than `h` allows.  Pick `h` accordingly.
    pub fn set_source(&mut self, name: &str, value: f64) -> Result<(), SimError> {
        // ponytail: the lower-casing allocates on every call (~36 ns measured,
        // against ≥560 ns of solver work for the smallest useful circuit).  Only
        // worth removing if a host crosses a wide bus every step, and the fix
        // then is a resolve-once handle (`fc_source_handle`), not a cleverer
        // lookup — an exact-match fast path was measured and it was slower,
        // because parsed names are already lower case so it always missed.
        let idx = *self.sources.get(&name.to_lowercase()).ok_or_else(|| {
            SimError::ParameterError(format!("no voltage or current source named '{name}'"))
        })?;
        match &mut self.netlist.elements[idx] {
            Element::VoltageSource { waveform, .. } | Element::CurrentSource { waveform, .. } => {
                *waveform = Waveform::Dc(value);
            }
            _ => unreachable!("sources map only holds V/I source indices"),
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tran::tran_nr_with_registry_opts;
    use fairchild_parser::parse_spice;

    fn registry_for(netlist: &Netlist) -> DeviceRegistry {
        let mut registry = DeviceRegistry::new();
        registry.register_builtin_models(&netlist.models);
        registry
    }

    /// The stepper and the batch path are the same integrator, so driving the
    /// stepper over a plain netlist must reproduce the batch waveform exactly —
    /// not approximately.  This is the guard on the tran.rs refactor.
    #[test]
    fn stepper_matches_batch_bit_for_bit() {
        let src =
            "* rc\nV1 in 0 PULSE(0 1 0 1n 1n 50n 100n)\nR1 in out 1k\nC1 out 0 1p\n.tran 1n 100n\n";
        let netlist = parse_spice(src).unwrap();
        let registry = registry_for(&netlist);
        let opts = SimOptions::from_netlist(&netlist);

        let batch = tran_nr_with_registry_opts(&netlist, 1e-9, 100e-9, &registry, &opts).unwrap();

        // 100 steps of 1 ns, plus the t=0 operating point.  The batch path can
        // emit one extra trailing point because its final step is clamped to
        // land exactly on `stop` after float accumulation has already put it
        // there — pre-existing behaviour, and not what this test is about.
        const N: usize = 100;
        let mut st = TranStepper::new(netlist.clone(), &registry, &opts, 1e-9).unwrap();
        let mut stepped = vec![(st.time(), st.node("out").unwrap())];
        for _ in 0..N {
            st.step().unwrap();
            stepped.push((st.time(), st.node("out").unwrap()));
        }

        let batch_out = &batch.node_voltages["out"];
        assert!(
            batch_out.len() > N,
            "batch produced {} points",
            batch_out.len()
        );
        for (i, (t, v)) in stepped.iter().enumerate() {
            assert_eq!(*t, batch.time[i], "timepoint {i} time diverged");
            assert_eq!(*v, batch_out[i], "timepoint {i} value diverged");
        }
    }

    /// The whole point of the stepper: a value written between steps changes
    /// the circuit's future, and the analog state can be read back to decide
    /// what to write next.
    #[test]
    fn set_source_drives_the_circuit() {
        let src = "* rc\nV1 in 0 DC 0\nR1 in out 1k\nC1 out 0 1p\n.tran 1n 100n\n";
        let netlist = parse_spice(src).unwrap();
        let registry = registry_for(&netlist);
        let opts = SimOptions::from_netlist(&netlist);

        let mut st = TranStepper::new(netlist, &registry, &opts, 1e-10).unwrap();
        assert!(st.node("out").unwrap().abs() < 1e-12, "starts discharged");

        // Bang-bang: drive high below 0.5 V, low above it.  A comparator in the
        // host program, which is the mixed-signal loop in miniature.
        let mut flips = 0;
        let mut driving_high = true;
        st.set_source("V1", 1.0).unwrap();
        for _ in 0..2000 {
            st.step().unwrap();
            let out = st.node("out").unwrap();
            let want_high = out < 0.5;
            if want_high != driving_high {
                driving_high = want_high;
                flips += 1;
                st.set_source("V1", if want_high { 1.0 } else { 0.0 })
                    .unwrap();
            }
        }
        assert!(flips >= 2, "expected the loop to toggle, got {flips} flips");
        let out = st.node("out").unwrap();
        assert!(
            (out - 0.5).abs() < 0.05,
            "bang-bang should hold out near the 0.5 V threshold, got {out}"
        );
    }

    #[test]
    fn advance_to_lands_on_the_grid_and_never_walks_backwards() {
        let src = "* rc\nV1 in 0 DC 1\nR1 in out 1k\nC1 out 0 1p\n.tran 1n 100n\n";
        let netlist = parse_spice(src).unwrap();
        let registry = registry_for(&netlist);
        let opts = SimOptions::from_netlist(&netlist);
        let mut st = TranStepper::new(netlist, &registry, &opts, 1e-9).unwrap();

        assert!((st.advance_to(10e-9).unwrap() - 10e-9).abs() < 1e-18);
        // Already past: no-op, not a rewind and not an extra step.
        assert!((st.advance_to(5e-9).unwrap() - 10e-9).abs() < 1e-18);
        // Off-grid target rounds up to the next grid point.
        assert!((st.advance_to(10.5e-9).unwrap() - 11e-9).abs() < 1e-18);
    }

    #[test]
    fn set_source_rejects_unknown_names() {
        let src = "* rc\nV1 in 0 DC 1\nR1 in out 1k\nC1 out 0 1p\n.tran 1n 100n\n";
        let netlist = parse_spice(src).unwrap();
        let registry = registry_for(&netlist);
        let opts = SimOptions::from_netlist(&netlist);
        let mut st = TranStepper::new(netlist, &registry, &opts, 1e-9).unwrap();
        assert!(st.set_source("vnope", 1.0).is_err());
        assert!(st.set_source("R1", 1.0).is_err(), "not a source");
        assert!(
            st.set_source("v1", 1.0).is_ok(),
            "name match is case-insensitive"
        );
    }
}
