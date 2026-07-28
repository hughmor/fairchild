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

use indexmap::IndexMap;
use std::collections::HashMap;

use fairchild_parser::{Element, Netlist, Waveform};

use crate::device::{Device, EvalFlags, SimContext};
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::mna::{CircuitTopology, MnaMatrix, StampPlan};
use crate::newton::{build_devices_with_footprints, dc_op_nr_with_registry_opts};
use crate::options::SimOptions;
use crate::solver::{Factorisation, LinearSolver};
use crate::tran::{
    advance_companions, advance_device_reactive_state, init_companions, init_device_reactive_state,
    stamp_device_reactive_companions, IntegratorMode,
};

/// A transient analysis paused between timesteps.
///
/// The step size is fixed for the lifetime of the stepper — the reactive
/// companion models (and the device-internal ones) are built for one `h` and
/// advanced incrementally from it.
///
// ponytail: fixed h only. Landing on an arbitrary externally-chosen event time
// needs companions rebuilt from raw state for a new h — the machinery is in
// `tran::tran_nr_with_registry_var_opts` (see its `cap_v` / `ind_i` raw maps).
// Lift that in if the digital side ever needs off-grid event times; a
// clock-locked host does not.
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
    cap_state: IndexMap<String, (f64, f64)>,
    ind_state: IndexMap<String, (f64, f64)>,
    dev_reactive_state: Vec<Vec<(f64, f64)>>,
    step: f64,
    t: f64,
    first_tr: bool,
    /// Lower-cased V/I source name → index into `netlist.elements`, so
    /// `set_source` is a lookup rather than a scan of every element.
    sources: HashMap<String, usize>,
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
        if opts.sanity_check && opts.uic {
            crate::sanity::check_netlist_sanity(&netlist);
        }
        crate::connectivity::check_connectivity(&netlist)?;
        let ctx = opts.sim_context();
        let mode = opts.method;

        // With UIC: skip DC OP, seed x from `.ic` (or 0 where unspecified).
        // Without UIC (the default): use DC OP as t=0 condition.
        let (mut topo, mut x) = if opts.uic {
            let topo = CircuitTopology::build(&netlist);
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

        // Seed x_tprev from DC OP (or UIC initial conditions) so reactive
        // history is defined before the first step.
        for dev in &mut devices {
            dev.commit_timestep(&x);
        }

        // Honour opts.max_step as an upper bound on the step size.
        let step = step.min(opts.max_step);

        let (cap_state, ind_state) = init_companions(&netlist, &topo, step, &x, mode);
        let dev_reactive_state = init_device_reactive_state(&devices, &x, step, mode);
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
            cap_state,
            ind_state,
            dev_reactive_state,
            step,
            t: 0.0,
            first_tr: true,
            sources,
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
        let alpha = 1.0 / self.step;
        // Expose the absolute time of this step to devices (delay lines look up
        // historical port values at `time_s − τ`).
        self.ctx.time_s = t_next;
        // Solve into a scratch vector so a non-converging step leaves the last
        // accepted solution intact for the caller to inspect or retry from.
        let mut x_try = self.x.clone();

        for _iter in 0..self.opts.itl4 {
            crate::mna::stamp_netlist_in_place(
                &mut self.mat,
                &self.topo,
                &self.netlist,
                t_next,
                &self.cap_state,
                &self.ind_state,
                Some(&self.plan),
            );

            for dev in &mut self.devices {
                dev.set_source_scale(1.0);
                dev.eval(&x_try, EvalFlags::tran(), &self.ctx);
                dev.load_residual_tran(&mut self.mat.b, alpha);
                dev.load_jacobian_tran(&mut self.mat, alpha);
            }
            // Stamp integrator-managed reactive companions for every
            // device-declared linear reactive branch (uses the device's
            // current bias-dependent value AND the history from the
            // previous accepted timestep).
            stamp_device_reactive_companions(
                &self.devices,
                &self.dev_reactive_state,
                &mut self.mat,
                self.step,
            );

            self.topo.stamp_gmin(&mut self.mat.a, self.opts.gmin);

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

            let converged = x_next
                .iter()
                .zip(x_try.iter())
                .all(|(n, o)| (n - o).abs() < self.opts.vntol + self.opts.reltol * n.abs());

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
        advance_companions(
            &self.netlist,
            &self.topo,
            self.step,
            &self.x,
            &mut self.cap_state,
            &mut self.ind_state,
            self.mode,
            self.first_tr,
        );
        advance_device_reactive_state(
            &self.devices,
            &self.x,
            &mut self.dev_reactive_state,
            self.step,
        );
        self.first_tr = false;
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
        let src = "* rc\nV1 in 0 PULSE(0 1 0 1n 1n 50n 100n)\nR1 in out 1k\nC1 out 0 1p\n.tran 1n 100n\n.end\n";
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
        let src = "* rc\nV1 in 0 DC 0\nR1 in out 1k\nC1 out 0 1p\n.tran 1n 100n\n.end\n";
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
        let src = "* rc\nV1 in 0 DC 1\nR1 in out 1k\nC1 out 0 1p\n.tran 1n 100n\n.end\n";
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
        let src = "* rc\nV1 in 0 DC 1\nR1 in out 1k\nC1 out 0 1p\n.tran 1n 100n\n.end\n";
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
