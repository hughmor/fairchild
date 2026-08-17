/// Transient integrators: fixed-step BE/TR/GEAR and variable-step BE+LTE.
///
/// `tran_nr` / `tran_nr_tr` — fixed-step, Newton-Raphson.
/// `tran_nr_var`            — variable-step with LTE control.
///
/// All of them go through `TranStepper`, which handles a linear circuit as the
/// degenerate case of a nonlinear one: two Newton iterations, the second only
/// confirming convergence, over a matrix the solver's factorisation cache then
/// recognises as unchanged.
use indexmap::IndexMap;
use std::collections::HashSet;

use fairchild_parser::{Element, Netlist};

use crate::device::EvalFlags;
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::mna::CircuitTopology;
use crate::newton::{build_devices_with_footprints, dc_op_nr_with_registry_opts};
use crate::options::SimOptions;
use crate::tran_step::TranStepper;

/// Transient integration method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntegratorMode {
    /// Backward Euler (BDF-1): first-order, unconditionally stable, single-history.
    BackwardEuler,
    /// Trapezoidal Rule (TR / BDF-2-like): second-order, A-stable, minimal overhead.
    Trapezoidal,
    /// GEAR / BDF-2: second-order, L-stable, two-step history.  First step and
    /// the step after any rejection demote to BE (order control 1↔2).
    ///
    /// Applies uniformly to netlist L/C and device-declared branches — every
    /// method does, since `crate::reactive` interprets `mode` in one place.
    /// Only the variable-step integrator carries the two-timepoint history
    /// BDF-2 needs, so the fixed-step path demotes GEAR to BE throughout.
    Gear,
}

/// Output of a transient simulation.
pub struct TranResult {
    pub time: Vec<f64>,
    /// Node voltages over time: node name → time series.
    pub node_voltages: IndexMap<String, Vec<f64>>,
    /// Voltage-source branch currents over time: source name → time series.
    pub vsrc_currents: IndexMap<String, Vec<f64>>,
}

impl TranResult {
    /// Drop every timepoint before `tstart` — `.tran`'s third argument.
    ///
    /// SPICE's tstart selects what is *saved*, not where integration begins, so
    /// this runs on the finished result rather than gating the solver. A tstart
    /// past the end leaves the last point, because returning an empty waveform
    /// would be a worse answer than a short one.
    pub fn trim_before(&mut self, tstart: f64) {
        if tstart <= 0.0 || self.time.is_empty() {
            return;
        }
        let keep_from = self
            .time
            .partition_point(|&t| t < tstart)
            .min(self.time.len() - 1);
        if keep_from == 0 {
            return;
        }
        self.time.drain(..keep_from);
        for v in self.node_voltages.values_mut() {
            v.drain(..keep_from);
        }
        for i in self.vsrc_currents.values_mut() {
            i.drain(..keep_from);
        }
    }

    /// Node voltage at a specific time, with linear interpolation.
    /// Returns None if the node is unknown or t is out of range.
    pub fn voltage_at(&self, node: &str, t: f64) -> Option<f64> {
        if node == "0" || node == "gnd" {
            return Some(0.0);
        }
        let v_series = self.node_voltages.get(node)?;
        interp(&self.time, v_series, t)
    }

    /// Voltage-source current at a specific time, with linear interpolation.
    pub fn isrc_at(&self, vsrc_name: &str, t: f64) -> Option<f64> {
        let i_series = self.vsrc_currents.get(vsrc_name)?;
        interp(&self.time, i_series, t)
    }

    /// Write all waveforms as a Nutmeg ASCII rawfile (ngspice `rawread` format).
    ///
    /// One value per line; point index only on the first variable's line per point.
    pub fn write_nutmeg<W: std::io::Write>(&self, mut w: W, title: &str) -> std::io::Result<()> {
        let n_vars = 1 + self.node_voltages.len() + self.vsrc_currents.len();
        let n_pts = self.time.len();
        writeln!(w, "Title: {title}")?;
        writeln!(w, "Plotname: Transient Analysis")?;
        writeln!(w, "Flags: real")?;
        writeln!(w, "No. Variables: {n_vars}")?;
        writeln!(w, "No. Points: {n_pts}")?;
        writeln!(w, "Variables:")?;
        writeln!(w, "\t0\ttime\ttime")?;
        let mut idx = 1usize;
        for name in self.node_voltages.keys() {
            writeln!(w, "\t{idx}\tv({name})\tvoltage")?;
            idx += 1;
        }
        for name in self.vsrc_currents.keys() {
            writeln!(w, "\t{idx}\ti({name})\tcurrent")?;
            idx += 1;
        }
        writeln!(w, "Values:")?;
        for (ti, &t) in self.time.iter().enumerate() {
            // ngspice format: point index on the first variable's line only;
            // remaining variables each get their own tab-indented line.
            writeln!(w, " {ti}\t{t:.6e}")?;
            for v in self.node_voltages.values() {
                writeln!(w, "\t{:.6e}", v[ti])?;
            }
            for i in self.vsrc_currents.values() {
                writeln!(w, "\t{:.6e}", i[ti])?;
            }
        }
        Ok(())
    }

    /// Write all waveforms to CSV.
    ///
    /// Columns: `time`, then `V(<node>)` for every node, then `I(<vsrc>)` for every
    /// voltage source. Values are written in scientific notation.
    pub fn write_csv<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        // Header
        write!(w, "time")?;
        for name in self.node_voltages.keys() {
            write!(w, ",V({name})")?;
        }
        for name in self.vsrc_currents.keys() {
            write!(w, ",I({name})")?;
        }
        writeln!(w)?;
        // Rows
        for (ti, &t) in self.time.iter().enumerate() {
            write!(w, "{t:.6e}")?;
            for v in self.node_voltages.values() {
                write!(w, ",{:.6e}", v[ti])?;
            }
            for i in self.vsrc_currents.values() {
                write!(w, ",{:.6e}", i[ti])?;
            }
            writeln!(w)?;
        }
        Ok(())
    }
}

fn interp(xs: &[f64], ys: &[f64], x: f64) -> Option<f64> {
    if xs.is_empty() {
        return None;
    }
    if x <= xs[0] {
        return Some(ys[0]);
    }
    if x >= *xs.last().unwrap() {
        return Some(*ys.last().unwrap());
    }
    let i = xs.partition_point(|&xi| xi <= x).saturating_sub(1);
    let t = (x - xs[i]) / (xs[i + 1] - xs[i]);
    Some(ys[i] + t * (ys[i + 1] - ys[i]))
}

// ---------------------------------------------------------------------------
// Helpers shared by the transient solvers
// ---------------------------------------------------------------------------

/// Currents through a coupled inductor pair after a solved step.
///
/// The standalone update `i = G_eq·v + I_hist` misses the mutual term, so using
/// it on a coupled pair lets the stored current drift from what was actually
/// stamped.  The mutual form is
///
/// ```text
///   I_L1 = G11·V_L1 + G12·V_L2 + I_hist1
///   I_L2 = G12·V_L1 + G22·V_L2 + I_hist2
/// ```
///
/// `scale` is the conductance scale the standalone companions used — `G_eq · L`,
/// which is `h` under BE and `1/α` under BDF-2, so the same expression covers
/// both.  `i_hist` are the history terms that were *stamped* for this step, so
/// whatever method produced them is already folded in.
///
/// Returns `None` for a degenerate pair (k → 1), which the stamper also skips.
pub(crate) fn coupled_inductor_currents(
    l1_value: f64,
    l2_value: f64,
    coupling: f64,
    scale: f64,
    vl1: f64,
    vl2: f64,
    i_hist1: f64,
    i_hist2: f64,
) -> Option<(f64, f64)> {
    let m = coupling * (l1_value * l2_value).sqrt();
    let det = l1_value * l2_value - m * m; // = L1·L2·(1−k²)
    if det.abs() < 1e-40 {
        return None;
    }
    let g11 = scale * l2_value / det;
    let g22 = scale * l1_value / det;
    let g12 = -scale * m / det;
    Some((
        g11 * vl1 + g12 * vl2 + i_hist1,
        g12 * vl1 + g22 * vl2 + i_hist2,
    ))
}

/// Append one time-point to a TranResult.
pub(crate) fn push_timepoint(result: &mut TranResult, t: f64, topo: &CircuitTopology, x: &[f64]) {
    result.time.push(t);
    for (name, &idx) in &topo.node_index {
        result.node_voltages.get_mut(name).unwrap().push(x[idx]);
    }
    let n = topo.n_nodes();
    for (name, &idx) in &topo.vsrc_index {
        result.vsrc_currents.get_mut(name).unwrap().push(x[n + idx]);
    }
}

// ---------------------------------------------------------------------------
// Nonlinear transient solver
// ---------------------------------------------------------------------------

/// Fixed-step Backward Euler transient with Newton-Raphson and a pre-built device registry.
///
/// Honours `.options` directives from the netlist (but forces `method=be`).
pub fn tran_nr_with_registry(
    netlist: &Netlist,
    step: f64,
    stop: f64,
    registry: &DeviceRegistry,
) -> Result<TranResult, SimError> {
    let mut opts = SimOptions::from_netlist(netlist);
    opts.method = IntegratorMode::BackwardEuler;
    tran_nr_with_registry_opts(netlist, step, stop, registry, &opts)
}

/// Fixed-step transient with explicit `SimOptions`.
///
/// The integration method (`BackwardEuler` or `Trapezoidal`) comes from
/// `opts.method`.  Tolerances, max NR iterations, gmin, vmax, etc. all read
/// from `opts`.
///
/// This is a driver over [`TranStepper`] — the integrator itself lives there so
/// that batch runs and host-driven mixed-signal stepping cannot drift apart.
pub fn tran_nr_with_registry_opts(
    netlist: &Netlist,
    step: f64,
    stop: f64,
    registry: &DeviceRegistry,
    opts: &SimOptions,
) -> Result<TranResult, SimError> {
    // The stepper owns its netlist so it can rewrite source waveforms between
    // steps.  One clone per transient run, against thousands of timesteps.
    let mut st = TranStepper::new(netlist.clone(), registry, opts, step)?;
    let step = st.step_size();

    let n_steps = ((stop / step).ceil() as usize) + 2;
    let topo = st.topology();
    let mut result = TranResult {
        time: Vec::with_capacity(n_steps),
        node_voltages: topo
            .node_index
            .keys()
            .map(|k| (k.clone(), Vec::with_capacity(n_steps)))
            .collect(),
        vsrc_currents: topo
            .vsrc_index
            .keys()
            .map(|k| (k.clone(), Vec::with_capacity(n_steps)))
            .collect(),
    };

    // Store t = 0 from DC OP.
    push_timepoint(&mut result, 0.0, st.topology(), st.solution());

    // The first timepoint is `step` even when that overshoots `stop` (a
    // stop < step run still produces one solved point); every later one is
    // clamped so the run lands exactly on `stop`.
    let mut t_next = step;
    loop {
        st.solve_at(t_next)?;
        st.commit(t_next);
        push_timepoint(&mut result, st.time(), st.topology(), st.solution());
        if st.time() >= stop {
            break;
        }
        st.advance_history();
        t_next = (st.time() + step).min(stop);
    }

    result.trim_before(opts.tstart);
    Ok(result)
}

/// Fixed-step Backward Euler transient using only built-in models from `.model` cards.
pub fn tran_nr(netlist: &Netlist, step: f64, stop: f64) -> Result<TranResult, SimError> {
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&netlist.models);
    tran_nr_with_registry(netlist, step, stop, &registry)
}

/// Fixed-step Trapezoidal Rule transient with Newton-Raphson and a pre-built registry.
///
/// Honours `.options` directives from the netlist (but forces `method=tr`).
pub fn tran_nr_with_registry_tr(
    netlist: &Netlist,
    step: f64,
    stop: f64,
    registry: &DeviceRegistry,
) -> Result<TranResult, SimError> {
    let mut opts = SimOptions::from_netlist(netlist);
    opts.method = IntegratorMode::Trapezoidal;
    tran_nr_with_registry_opts(netlist, step, stop, registry, &opts)
}

/// Fixed-step Trapezoidal Rule transient using only built-in models from `.model` cards.
pub fn tran_nr_tr(netlist: &Netlist, step: f64, stop: f64) -> Result<TranResult, SimError> {
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&netlist.models);
    tran_nr_with_registry_tr(netlist, step, stop, &registry)
}

// ---------------------------------------------------------------------------
// Variable-step BE + LTE nonlinear transient solver
// ---------------------------------------------------------------------------

/// Variable-step Backward Euler transient with Newton-Raphson and LTE timestep control.
///
/// Honours `.options` directives from the netlist.
pub fn tran_nr_with_registry_var(
    netlist: &Netlist,
    step: f64,
    stop: f64,
    registry: &DeviceRegistry,
) -> Result<TranResult, SimError> {
    tran_nr_with_registry_var_opts(
        netlist,
        step,
        stop,
        registry,
        &SimOptions::from_netlist(netlist),
    )
}

/// Variable-step transient with explicit `SimOptions`.
///
/// `step` is the maximum allowed timestep (upper bound, further capped by
/// `opts.max_step`).  Every accepted internal step is stored;
/// `TranResult::voltage_at` interpolates between them.
///
/// Timestep control uses a first-order predictor-corrector LTE estimate:
///   LTE_norm = max |x_corr − x_pred| * 0.5 / (vntol + reltol·|x|)
/// Accept if LTE_norm ≤ 1; adjust h_new = h·(0.9/LTE_norm)^0.5, clamped to [0.1h, 4h]·min(step).
pub fn tran_nr_with_registry_var_opts(
    netlist: &Netlist,
    step: f64,
    stop: f64,
    registry: &DeviceRegistry,
    opts: &SimOptions,
) -> Result<TranResult, SimError> {
    // Adaptive step control and injected noise do not mix: the LTE estimator
    // reads a fresh random sample as a fast signal and shrinks the step to
    // chase it, and the step size then becomes correlated with the noise, which
    // biases the spectrum it was meant to reproduce.  Fixed steps are what SDE
    // solvers use, for this reason.  Refusing beats quietly returning a
    // plausible waveform with the wrong noise in it.
    if opts.trannoise {
        return Err(SimError::ParameterError(
            "transient noise needs a fixed timestep; set `.options variable_step=0` \
             (the LTE controller would chase the noise and bias its spectrum)"
                .into(),
        ));
    }
    if opts.sanity_check && opts.uic {
        crate::sanity::check_netlist_sanity(netlist);
    }
    crate::connectivity::check_connectivity(netlist)?;
    let mut ctx = opts.sim_context();
    let step = step.min(opts.max_step);
    let (topo, mut x) = if opts.uic {
        let topo = CircuitTopology::build(netlist);
        let mut x = vec![0.0f64; topo.size];
        for (name, value) in &netlist.ic {
            if let Some(&i) = topo.node_index.get(name) {
                x[i] = *value;
            }
        }
        (topo, x)
    } else {
        let dc = dc_op_nr_with_registry_opts(netlist, registry, opts)?;
        (dc.topo, dc.x)
    };
    let mut topo = topo;
    let (mut devices, footprints) =
        build_devices_with_footprints(netlist, &mut topo, &ctx, registry)?;
    // Pad x for any OSDI internal-node rows allocated by build_devices.
    x.resize(topo.size, 0.0);
    let solver = opts.linear_solver(topo.size);
    // Structural sparsity pattern: same for every timestep and every NR
    // iteration within them, so build it once here.
    let plan = crate::mna::StampPlan::new(&topo, netlist, &footprints);

    let n_nodes = topo.n_nodes();
    let h_min = step * 1e-6;
    // Not every unknown is a volt — see `crate::tolerance`.  Serves both the NR
    // step test and the LTE norm below, so the step controller weighs a λ error
    // on the same scale that decides convergence.
    let tol = crate::tolerance::Tolerances::build(netlist, &topo, opts);

    // Nodes constrained by voltage sources are excluded from LTE: their voltages
    // change due to the source waveform, not integration error.
    let forced_nodes: HashSet<usize> = netlist
        .elements
        .iter()
        .filter_map(|el| {
            if let Element::VoltageSource { pos, neg, .. } = el {
                Some([
                    topo.node_index.get(pos).copied(),
                    topo.node_index.get(neg).copied(),
                ])
            } else {
                None
            }
        })
        .flatten()
        .flatten()
        .collect();

    for dev in &mut devices {
        dev.commit_timestep(&x);
    }

    // All reactive history — netlist C/L and device-declared branches alike —
    // in the one representation `crate::reactive` owns.  Physical state, so the
    // companion can be rebuilt for whatever `h` this step turns out to use.
    let mut reactive = crate::reactive::ReactiveState::new(netlist, &topo, &mut devices, &ctx, &x);
    // Track previous accepted timestep for non-uniform BDF-2 (h_prev_accepted
    // is the step that took us from t_{n-1} to t_n; distinct from h_prev which
    // is the trial step for predictor extrapolation).
    let mut h_prev_accepted = 0.0_f64;

    let n_hint = ((stop / step).ceil() as usize) + 2;
    let mut result = TranResult {
        time: Vec::with_capacity(n_hint),
        node_voltages: topo
            .node_index
            .keys()
            .map(|k| (k.clone(), Vec::with_capacity(n_hint)))
            .collect(),
        vsrc_currents: topo
            .vsrc_index
            .keys()
            .map(|k| (k.clone(), Vec::with_capacity(n_hint)))
            .collect(),
    };
    push_timepoint(&mut result, 0.0, &topo, &x);

    let mut t = 0.0_f64;
    let mut h = step;
    let mut h_prev = 0.0_f64;
    let mut x_prev = x.clone();
    let mut consecutive_rejects = 0usize;

    // Set to `true` whenever an accepted step lands on a known waveform
    // breakpoint (slope discontinuity).  The linear predictor / LTE
    // estimator cannot see through a kink — they assume the system
    // evolves smoothly, so using x_prev from the steep-slope side
    // pushes x_pred far away from x_try and the LTE refuses to drop
    // below 1 no matter how small h becomes (the error isn't from h,
    // it's from the predictor extrapolating across the kink).
    //
    // When the flag is set, the *next* step uses zeroth-order
    // prediction (x_pred = x) and bypasses the LTE check entirely,
    // accepting whatever NR converges to.  This matches SPICE's
    // standard breakpoint-restart behaviour.
    let mut just_crossed_breakpoint = false;

    // Cached factorisation across the variable-step transient run —
    // see the fixed-step path above for the same pattern; the sparsity
    // is fixed across the entire integration so one symbolic factor-
    // isation is amortised over thousands of refactors.
    let mut fact: Option<Box<dyn crate::solver::Factorisation>> = None;

    // Reusable MnaMatrix — see fixed-step path for rationale.
    let mut mat = crate::mna::MnaMatrix::with_pattern(topo.size, plan.pattern.clone());

    'outer: loop {
        if t >= stop {
            break;
        }

        // Clamp h to the next waveform slope discontinuity so we always land on it.
        let next_bp: Option<f64> = netlist
            .elements
            .iter()
            .filter_map(|el| match el {
                Element::VoltageSource { waveform, .. }
                | Element::CurrentSource { waveform, .. } => waveform.next_breakpoint(t),
                _ => None,
            })
            .reduce(f64::min);
        let mut h_want = h.min(stop - t);
        if let Some(bp) = next_bp {
            let to_bp = bp - t;
            if to_bp > h_min {
                h_want = h_want.min(to_bp);
            }
        }
        let h_actual = h_want.max(h_min);
        let t_next = t + h_actual;
        // Expose the trial-step time to devices (delay lines look up historical
        // port values at `time_s − τ`).
        ctx.time_s = t_next;

        // Order control: GEAR-2 needs two accepted steps of history.  Demote
        // to BE on the first two steps, after any rejection (history stale),
        // and on extreme step ratios where BDF-2 would amplify noise.
        let step_ratio = if h_prev_accepted > 0.0 {
            h_actual / h_prev_accepted
        } else {
            1.0
        };
        let use_gear2 = matches!(opts.method, IntegratorMode::Gear)
            && h_prev_accepted > 0.0
            && reactive.gear2_ready()
            && consecutive_rejects == 0
            && step_ratio > 0.25
            && step_ratio < 4.0;
        let gear2 = if use_gear2 {
            Some(h_prev_accepted)
        } else {
            None
        };
        // Trapezoidal takes its first step with Backward Euler, matching the
        // fixed-step path — and available here at all only because history is
        // physical rather than companion-shaped.
        let step_mode =
            if h_prev_accepted == 0.0 && matches!(opts.method, IntegratorMode::Trapezoidal) {
                IntegratorMode::BackwardEuler
            } else {
                opts.method
            };

        // Companions for this trial `h`, for every reactive branch at once.
        reactive.build(&devices, step_mode, h_actual, gear2);

        // Predictor: linear extrapolation. Zero-order on the first step
        // (no history) and immediately after crossing a waveform
        // breakpoint (history is on the wrong side of a slope kink).
        let x_pred: Vec<f64> = if h_prev > 0.0 && !just_crossed_breakpoint {
            let scale = h_actual / h_prev;
            x.iter()
                .zip(x_prev.iter())
                .map(|(xi, xp)| xi + scale * (xi - xp))
                .collect()
        } else {
            x.clone()
        };

        // NR corrector starting from x_pred.
        let alpha = 1.0 / h_actual;
        // Devices that stamp their own reactance (OSDI/Verilog-A `ddt`) need the
        // method, not just `alpha` — which can only express Backward Euler.
        // Same `step_mode` and `gear2` the branch stamper below uses, so one
        // decision still reaches everything.
        ctx.discretisation = Some(crate::device::Discretisation {
            mode: step_mode,
            h: h_actual,
            gear2_h_prev: gear2,
        });
        let mut x_try = x_pred.clone();
        let mut nr_converged = false;

        for _iter in 0..opts.itl4 {
            crate::mna::stamp_netlist_in_place(
                &mut mat,
                &topo,
                netlist,
                t_next,
                &reactive.cap_state,
                &reactive.ind_state,
                Some(&plan),
                crate::mna::InductorDc::Short,
            );

            for dev in devices.iter_mut() {
                dev.set_source_scale(1.0);
                dev.eval(&x_try, EvalFlags::tran(), &ctx);
                dev.load_residual_tran(&mut mat.b, alpha);
                dev.load_jacobian_tran(&mut mat, alpha);
            }
            // Device-declared reactive branches, same method and order control
            // as the netlist C/L above — one decision, applied to everything.
            crate::reactive::stamp_device_branches(
                &devices,
                &reactive.dev_state,
                &mut mat,
                &x_try,
                h_actual,
                step_mode,
                gear2,
            );

            topo.stamp_gmin(&mut mat.a, opts.gmin);

            let x_new = if let Some(f) = fact.as_mut() {
                f.refactor_and_solve_mat(&mat)?
            } else {
                let mut f = solver.factorise_mat(&mat)?;
                let r = f.refactor_and_solve_mat(&mat)?;
                fact = Some(f);
                r
            };

            let max_dv = x_new
                .iter()
                .zip(x_try.iter())
                .take(n_nodes)
                .map(|(n, o)| (n - o).abs())
                .fold(0.0f64, f64::max);

            let x_next: Vec<f64> = if max_dv > opts.vmax {
                let scale = opts.vmax / max_dv;
                x_try
                    .iter()
                    .zip(x_new.iter())
                    .map(|(o, n)| o + scale * (n - o))
                    .collect()
            } else {
                x_new
            };

            let converged = tol.converged(&x_next, &x_try);

            x_try = x_next;
            if converged {
                nr_converged = true;
                break;
            }
        }

        if !nr_converged {
            consecutive_rejects += 1;
            if consecutive_rejects > opts.max_rejections {
                return Err(SimError::NoConvergence { iters: opts.itl4 });
            }
            h = (h_actual * 0.5).max(h_min);
            continue 'outer;
        }

        // LTE estimate: skipped on the first step (no predictor history)
        // and on the step right after crossing a waveform breakpoint
        // (predictor is intentionally zeroth-order there, so the
        // difference `x_try - x_pred` reflects the slope discontinuity
        // rather than integration error).
        let lte_norm: f64 = if h_prev > 0.0 && !just_crossed_breakpoint {
            x_try
                .iter()
                .zip(x_pred.iter())
                .enumerate()
                .take(n_nodes)
                .filter(|(idx, _)| !forced_nodes.contains(idx))
                .map(|(idx, (xc, xp))| (xc - xp).abs() * 0.5 / tol.bound(idx, *xc))
                .fold(0.0f64, f64::max)
        } else {
            0.0
        };

        if lte_norm <= 1.0 || h_actual <= h_min {
            // Accept step.
            consecutive_rejects = 0;
            x_prev = std::mem::replace(&mut x, x_try);
            h_prev = h_actual;
            t = t_next;

            // Did this step land on (or close to) a waveform breakpoint?
            // If so, the next step starts on the other side of a slope
            // discontinuity — disable the linear predictor + LTE for
            // that step.  `next_bp` was computed at the top of this
            // iteration relative to the *old* `t`; t_next equals it
            // (within float epsilon) exactly when the breakpoint
            // clamp determined h_actual.
            just_crossed_breakpoint = match next_bp {
                Some(bp) => (t - bp).abs() <= h_actual * 1e-9 + h_min,
                None => false,
            };

            h_prev_accepted = h_actual;

            for dev in &mut devices {
                dev.commit_timestep(&x);
            }

            // Roll every reactive branch's history forward — netlist C/L, the
            // coupled-pair mutual correction, and device branches, in one call.
            // Needs neither the method nor `h`: the companions it reads are the
            // ones that were stamped, so all of that is already folded in.
            reactive.accept(&devices, &x);

            push_timepoint(&mut result, t, &topo, &x);

            h = if lte_norm < 1e-10 {
                h_actual * 2.0
            } else {
                h_actual * (0.9 / lte_norm).sqrt()
            };
            h = h.clamp(h_actual * 0.1, h_actual * 4.0).min(step);
        } else {
            consecutive_rejects += 1;
            if consecutive_rejects > opts.max_rejections {
                return Err(SimError::NoConvergence { iters: opts.itl4 });
            }
            h = (h_actual * (0.9 / lte_norm).sqrt()).max(h_min);
        }
    }

    result.trim_before(opts.tstart);
    Ok(result)
}

/// Variable-step BE + LTE transient using only built-in models from `.model` cards.
pub fn tran_nr_var(netlist: &Netlist, step: f64, stop: f64) -> Result<TranResult, SimError> {
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&netlist.models);
    tran_nr_with_registry_var(netlist, step, stop, &registry)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

    /// Backward Euler through the Newton path. `tran_nr` takes its method from
    /// the netlist, which defaults to trapezoidal, and several tests below pin
    /// BE's error behaviour specifically.
    fn tran_be(netlist: &Netlist, step: f64, stop: f64) -> TranResult {
        let mut registry = DeviceRegistry::new();
        registry.register_builtin_models(&netlist.models);
        let opts = SimOptions {
            method: IntegratorMode::BackwardEuler,
            ..SimOptions::from_netlist(netlist)
        };
        tran_nr_with_registry_opts(netlist, step, stop, &registry, &opts).unwrap()
    }

    // ---------- tran_nr tests ----------

    #[test]
    fn tran_nr_diode_steady_state() {
        // R-D series, constant V=5V, no reactive elements.
        // tran_nr result at t=1µs must match dc_op_nr within 0.1%.
        let netlist_str = "* Diode DC via transient\nVdd a 0 DC 5\nR1 a b 10k\nD1 b 0 myd\n\
             .model myd D (Is=1e-14 N=1)\n.tran 1u 10u\n.end\n";
        let netlist = parse_spice(netlist_str).unwrap();

        let dc = crate::newton::dc_op_nr(&netlist).unwrap();
        let vb_dc = dc.node_voltage("b").unwrap();

        let tr = tran_nr(&netlist, 1e-6, 10e-6).unwrap();
        // At the last time step the solution should still be the DC OP.
        let vb_tran = tr.voltage_at("b", 10e-6).unwrap();
        assert!(
            (vb_tran - vb_dc).abs() < 1e-6,
            "V(b) tran_nr={vb_tran:.6e}  dc_op_nr={vb_dc:.6e}"
        );
    }

    // ---------- fixed-step regression tests ----------

    #[test]
    fn rc_step_response_shape() {
        // R=1k C=1µF τ=1ms, step to 1V at t=0. V(out) should be ~0.632 at t=τ.
        let netlist = parse_spice(
            "* RC step\nV1 in 0 PULSE(0 1 0 1n 1n 10m 20m)\nR1 in out 1k\nC1 out 0 1u\n.tran 1u 5m\n.end\n"
        ).unwrap();

        let result = tran_be(&netlist, 1e-6, 5e-3);

        let v_1tau = result.voltage_at("out", 1e-3).unwrap();
        let v_5tau = result.voltage_at("out", 5e-3).unwrap();

        // At t=τ=1ms: exact = 1-e^-1 ≈ 0.6321. BE with h=1µs ≈ 0.6314 (0.1% error).
        assert!((v_1tau - 0.6321).abs() < 0.01, "v(out) at 1τ = {v_1tau:.4}");
        // At t=5τ: essentially fully charged.
        assert!(v_5tau > 0.99, "v(out) at 5τ = {v_5tau:.4}");
        // Monotonically increasing (no oscillation).
        assert!(v_1tau < v_5tau);
    }

    #[test]
    fn rl_step_response_shape() {
        // R=1k L=1H τ=1ms, 1V step. I_L ramps up; V(R) = R*I_L(t).
        // V(out) = V_R = 1 * (1 - e^{-t/τ}). Same shape as RC.
        let netlist = parse_spice(
            "* RL step\nV1 in 0 PULSE(0 1 0 1n 1n 10m 20m)\nR1 in out 1k\nL1 out 0 1\n.tran 1u 5m\n.end\n"
        ).unwrap();

        let result = tran_be(&netlist, 1e-6, 5e-3);

        let v_1tau = result.voltage_at("out", 1e-3).unwrap();
        let v_5tau = result.voltage_at("out", 5e-3).unwrap();

        // V(out) = V across L = Vs - I*R = 1 - (1-e^-t/τ) = e^{-t/τ}
        // At t=τ: e^-1 ≈ 0.368
        // V(in)-V(out) = I*R = 1*(1-e^-1) ≈ 0.632 → V(out) ≈ 0.368
        assert!((v_1tau - 0.3679).abs() < 0.01, "v(out) at 1τ = {v_1tau:.4}");
        // At t=5τ: V(out) ≈ e^-5 ≈ 0.0067 (nearly 0, almost all voltage on R)
        assert!(v_5tau < 0.01, "v(out) at 5τ = {v_5tau:.4}");
    }

    // ---------- Trapezoidal Rule tests ----------

    #[test]
    fn tr_rc_more_accurate_than_be() {
        // Use a large step (h = τ/5) so BE error is visible.
        // RC: R=1kΩ C=1µF → τ=1ms; test at t=τ.
        // Exact V(t=τ) = 1 − e^−1 ≈ 0.6321.
        // The source must *start* at 0: SPICE computes an operating point
        // before the transient, and a plain `DC 1` would charge the cap through
        // its open circuit at t=0, leaving nothing to integrate.
        let netlist = parse_spice(
            "* RC step\nV1 in 0 PULSE(0 1 0 1n 1n 10m 20m)\nR1 in out 1k\nC1 out 0 1u\n.tran 200u 5m\n.end\n",
        )
        .unwrap();
        let h = 200e-6; // 5 steps per τ — large enough to show BE error
        let exact = 1.0 - (-1.0_f64).exp(); // ≈ 0.6321

        let r_be = tran_be(&netlist, h, 5e-3);
        let r_tr = tran_nr_tr(&netlist, h, 5e-3).unwrap();

        let v_be = r_be.voltage_at("out", 1e-3).unwrap();
        let v_tr = r_tr.voltage_at("out", 1e-3).unwrap();

        let err_be = (v_be - exact).abs();
        let err_tr = (v_tr - exact).abs();

        assert!(err_tr < err_be, "TR should be more accurate than BE at same step: be_err={err_be:.4e} tr_err={err_tr:.4e}");
        assert!(
            err_tr < 0.01,
            "TR error at t=τ should be < 1%: {err_tr:.4e}"
        );
    }

    // ---------- Variable-step tests ----------

    #[test]
    fn tran_nr_var_rc_step_response() {
        // RC: R=1k C=1µF τ=1ms, step to 1V at t=0 (PULSE so DC OP has V=0).
        // With a large hint step the variable-step solver should take fewer steps
        // but still hit V(τ) ≈ 0.6321 within 1%.
        let netlist = parse_spice(
            "* RC var-step\nV1 in 0 PULSE(0 1 0 1n 1n 10m 20m)\nR1 in out 1k\nC1 out 0 1u\n.tran 500u 5m\n.end\n"
        ).unwrap();
        let result = tran_nr_var(&netlist, 500e-6, 5e-3).unwrap();

        let exact = 1.0 - (-1.0_f64).exp();
        let v_1tau = result.voltage_at("out", 1e-3).unwrap();
        let v_5tau = result.voltage_at("out", 5e-3).unwrap();

        // BE with h = τ/2 has ~6% inherent discretisation error; allow 15%.
        assert!(
            (v_1tau - exact).abs() < 0.15,
            "V(1τ) = {v_1tau:.4}, exact = {exact:.4}"
        );
        assert!(v_5tau > 0.95, "V(5τ) = {v_5tau:.4}");
        assert!(result.time.len() > 2, "must have at least a few timepoints");
    }

    #[test]
    fn tran_nr_var_matches_fixed_step() {
        // Variable-step result should agree with fixed-step BE within 0.1% at t=1ms.
        let netlist = parse_spice(
            "* RC\nV1 in 0 PULSE(0 1 0 1n 1n 10m 20m)\nR1 in out 1k\nC1 out 0 1u\n.tran 1u 2m\n.end\n"
        ).unwrap();
        let r_fixed = tran_nr(&netlist, 1e-6, 2e-3).unwrap();
        let r_var = tran_nr_var(&netlist, 1e-6, 2e-3).unwrap();

        let v_fixed = r_fixed.voltage_at("out", 1e-3).unwrap();
        let v_var = r_var.voltage_at("out", 1e-3).unwrap();
        assert!(
            (v_fixed - v_var).abs() < 1e-3,
            "fixed={v_fixed:.6}  var={v_var:.6}"
        );
    }

    /// The variable-step path used to drop every device-internal reactive
    /// branch: it never called init/stamp/advance_device_reactive_state, so a
    /// PN phase shifter's depletion C_j simply wasn't in the matrix.  The
    /// junction then settled within one timestep instead of over its RC, and
    /// nothing warned.
    ///
    /// Driven through 10 k with a near-open junction (g_pn = 1e-9), the ONLY
    /// thing setting the settling time is the device's own C_j: with c_j0 = 10 pF
    /// and a 1 V reverse step, tau ~ 100 ns at zero bias.  Both integrators must
    /// show that lag; an integrator that drops C_j jumps straight to the final
    /// value.
    #[test]
    fn var_step_includes_device_internal_capacitance() {
        let src = "* device-internal C_j\n\
             .optical_port lam\n.optical_port psout\n\
             Xlaser lam fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
             Xps lam psout a 0 fc_pn_ps_cap L_um=500 V_pi_L=2e-3 g_pn=1e-9 \
             c_j0=10p v_bi=0.917\n\
             VS src 0 PULSE(0 -1 20n 100p 100p 2u 4u)\n\
             RS src a 10k\n\
             .tran 5n 400n\n.end\n";
        let netlist = parse_spice(src).unwrap();
        let mut registry = DeviceRegistry::new();
        registry.register_builtin_models(&netlist.models);
        let opts = SimOptions::from_netlist(&netlist);

        let fixed = tran_nr_with_registry_opts(&netlist, 5e-9, 400e-9, &registry, &opts).unwrap();
        let var = tran_nr_with_registry_var_opts(&netlist, 5e-9, 400e-9, &registry, &opts).unwrap();

        // 25 ns after the edge — well inside one tau — the junction must still be
        // far from its final value in BOTH runs.
        let vf = fixed.voltage_at("a", 45e-9).unwrap();
        let vv = var.voltage_at("a", 45e-9).unwrap();
        let settled = fixed.voltage_at("a", 390e-9).unwrap();
        // Not exactly -1 V: g_pn plus the diagonal gmin form a slight divider.
        assert!(
            (settled + 1.0).abs() < 0.01,
            "fixed-step should reach essentially the full -1 V eventually, got {settled}"
        );
        assert!(
            vf / settled < 0.5,
            "fixed-step should be mid-slew at 45 ns, got {vf} of {settled}"
        );
        assert!(
            vv / settled < 0.5,
            "variable-step dropped the device's C_j: {vv} of {settled} at 45 ns \
             (fixed-step has {vf}) — the junction settled instantly"
        );
        // And the two integrators should agree closely across the whole edge.
        for &t in &[30e-9, 45e-9, 70e-9, 120e-9, 250e-9] {
            let (a, b) = (
                fixed.voltage_at("a", t).unwrap(),
                var.voltage_at("a", t).unwrap(),
            );
            assert!(
                (a - b).abs() < 0.05,
                "t={t:e}: fixed={a:.6} var={b:.6} differ by more than 50 mV"
            );
        }
    }

    /// A device-internal reactive branch and an equivalent netlist element are
    /// the same circuit, so every integrator must produce the same numbers to
    /// round-off — under every method.
    ///
    /// `m_j = 0` makes the PN phase shifter's C_j bias-independent, so
    /// `fc_pn_ps_cap` with `c_j0 = C` is exactly `fc_pn_ps` (no cap) plus a
    /// discrete `C` across the same nodes.  This is the strongest statement of
    /// what device-declared branches are supposed to mean, and it caught two
    /// real bugs: the variable-step path dropping them entirely, and the
    /// fixed-step path integrating them with Backward Euler no matter what
    /// `opts.method` said (a 23 mV, O(h) error on the *default* TR setting).
    #[test]
    fn device_branch_equals_equivalent_netlist_element() {
        const C: &str = "10p";
        let common = "L_um=500 V_pi_L=2e-3 g_pn=1e-9";
        let dev = format!(
            "* device-internal C_j\n.optical_port lam\n.optical_port psout\n\
             Xlaser lam fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
             Xps lam psout a 0 fc_pn_ps_cap {common} c_j0={C} v_bi=0.917 m_j=0\n\
             VS src 0 PULSE(0 -1 20n 100p 100p 2u 4u)\nRS src a 10k\n.tran 5n 300n\n.end\n"
        );
        let refr = format!(
            "* explicit C\n.optical_port lam\n.optical_port psout\n\
             Xlaser lam fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
             Xps lam psout a 0 fc_pn_ps {common}\n\
             VS src 0 PULSE(0 -1 20n 100p 100p 2u 4u)\nRS src a 10k\nC1 a 0 {C}\n\
             .tran 5n 300n\n.end\n"
        );
        let (nd, nr) = (parse_spice(&dev).unwrap(), parse_spice(&refr).unwrap());

        for mode in [
            IntegratorMode::BackwardEuler,
            IntegratorMode::Trapezoidal,
            IntegratorMode::Gear,
        ] {
            let mut opts = SimOptions::from_netlist(&nd);
            opts.method = mode;
            let reg = DeviceRegistry::new();

            for variable_step in [false, true] {
                let run = |n: &Netlist| {
                    if variable_step {
                        tran_nr_with_registry_var_opts(n, 5e-9, 300e-9, &reg, &opts).unwrap()
                    } else {
                        tran_nr_with_registry_opts(n, 5e-9, 300e-9, &reg, &opts).unwrap()
                    }
                };
                let (rd, rr) = (run(&nd), run(&nr));
                for &t in &[25e-9, 45e-9, 70e-9, 120e-9, 200e-9, 290e-9] {
                    let (a, b) = (
                        rd.voltage_at("a", t).unwrap(),
                        rr.voltage_at("a", t).unwrap(),
                    );
                    assert!(
                        (a - b).abs() < 1e-9,
                        "{mode:?} variable_step={variable_step} t={t:e}: \
                         device-internal C gives {a:.9}, explicit C gives {b:.9} — \
                         the integrator is not treating them as the same circuit"
                    );
                }
                // Sanity: the capacitance is actually doing something, so the
                // assertion above is not passing on two identical flat lines.
                let mid = rd.voltage_at("a", 45e-9).unwrap();
                let end = rd.voltage_at("a", 290e-9).unwrap();
                assert!(
                    mid / end < 0.6,
                    "{mode:?} variable_step={variable_step}: expected an RC lag, \
                     got {mid} at 45 ns vs {end} settled"
                );
            }
        }
    }

    /// The same "two spellings of one circuit" probe, aimed at the devices that
    /// stamp their own reactance instead of declaring a branch.
    ///
    /// `load_residual_tran` receives only `alpha = 1/h`, which is Backward Euler
    /// and cannot express anything else — so the diode's `Cj` and the MOSFET's
    /// Meyer caps were integrated with BE under *every* method, while a netlist
    /// `C` in the same circuit honoured `.options method`. TR is the default, so
    /// this was the common case: 19.4 mV (diode) and 9.9 mV (MOSFET) on a 2 V
    /// swing, and exactly 0 under `be`/`gear`, which is what hid it.
    ///
    /// Both devices are set up so their internal cap is *exactly* a linear
    /// capacitor — `m=0` makes the diode's depletion charge `Q = CJO·V` in
    /// closed form, and with `cox=0` the MOSFET's Meyer caps are pure overlap —
    /// so any difference is the integrator, not the physics.
    #[test]
    fn self_stamped_device_caps_equal_an_equivalent_netlist_c() {
        const C: &str = "10p";
        // 10 kΩ · 10 pF = 100 ns, sampled across a 300 ns window.
        let cases: [(&str, &str, &str); 3] = [
            (
                "diode Cj",
                &format!(
                    "* device-internal Cj\n\
                     .model dm D (IS=1e-14 N=1 CJO={C} VJ=0.917 M=0)\n\
                     V1 src 0 PULSE(0 -2 20n 100p 100p 2u 4u)\nR1 src a 10k\nD1 a 0 dm\n\
                     .tran 5n 300n\n.end\n"
                ),
                &format!(
                    "* explicit C\n\
                     .model dm D (IS=1e-14 N=1)\n\
                     V1 src 0 PULSE(0 -2 20n 100p 100p 2u 4u)\nR1 src a 10k\nD1 a 0 dm\n\
                     C1 a 0 {C}\n.tran 5n 300n\n.end\n"
                ),
            ),
            (
                "MOSFET Cgs",
                // cgso · W = 1e-6 · 10e-6 = 10 pF, and cox=0 keeps the
                // region-dependent channel caps out of it.
                "* device-internal Cgs\n\
                 .model nm NMOS (VTO=0.7 KP=100u CGSO=1e-6)\n\
                 VDD dd 0 2\nRD dd d 1k\n\
                 V1 src 0 PULSE(0 2 20n 100p 100p 2u 4u)\nR1 src a 10k\n\
                 M1 d a 0 0 nm w=10u l=1u\n.tran 5n 300n\n.end\n",
                &format!(
                    "* explicit C\n\
                     .model nm NMOS (VTO=0.7 KP=100u)\n\
                     VDD dd 0 2\nRD dd d 1k\n\
                     V1 src 0 PULSE(0 2 20n 100p 100p 2u 4u)\nR1 src a 10k\n\
                     M1 d a 0 0 nm w=10u l=1u\nC1 a 0 {C}\n.tran 5n 300n\n.end\n"
                ),
            ),
            (
                // MJE=0 makes the B-E depletion charge exactly linear, and TF=0
                // keeps the transit-time charge out of it. Emitter grounded, so
                // CJE is precisely a capacitor from the base to 0.
                "BJT Cje",
                "* device-internal Cje\n\
                 .model qm NPN (IS=1e-16 BF=100 CJE=10p VJE=0.75 MJE=0)\n\
                 VCC cc 0 5\nRC cc c 1k\n\
                 V1 src 0 PULSE(0 -2 20n 100p 100p 2u 4u)\nR1 src a 10k\n\
                 Q1 c a 0 qm\n.tran 5n 300n\n.end\n",
                &format!(
                    "* explicit C\n\
                     .model qm NPN (IS=1e-16 BF=100)\n\
                     VCC cc 0 5\nRC cc c 1k\n\
                     V1 src 0 PULSE(0 -2 20n 100p 100p 2u 4u)\nR1 src a 10k\n\
                     Q1 c a 0 qm\nC1 a 0 {C}\n.tran 5n 300n\n.end\n"
                ),
            ),
        ];

        for (what, dev, refr) in cases {
            let (nd, nr) = (parse_spice(dev).unwrap(), parse_spice(refr).unwrap());

            for mode in [
                IntegratorMode::BackwardEuler,
                IntegratorMode::Trapezoidal,
                IntegratorMode::Gear,
            ] {
                let mut opts = SimOptions::from_netlist(&nd);
                opts.method = mode;

                for variable_step in [false, true] {
                    let run = |n: &Netlist| {
                        // Per netlist: native D/M models reach the builder through
                        // the registry (`new()` alone gives UnknownModel), and the
                        // two decks deliberately carry *different* cards.
                        let mut reg = DeviceRegistry::new();
                        reg.register_builtin_models(&n.models);
                        if variable_step {
                            tran_nr_with_registry_var_opts(n, 5e-9, 300e-9, &reg, &opts).unwrap()
                        } else {
                            tran_nr_with_registry_opts(n, 5e-9, 300e-9, &reg, &opts).unwrap()
                        }
                    };
                    let (rd, rr) = (run(&nd), run(&nr));
                    for &t in &[25e-9, 45e-9, 70e-9, 120e-9, 200e-9, 290e-9] {
                        let (a, b) = (
                            rd.voltage_at("a", t).unwrap(),
                            rr.voltage_at("a", t).unwrap(),
                        );
                        assert!(
                            (a - b).abs() < 1e-9,
                            "{what} {mode:?} variable_step={variable_step} t={t:e}: \
                             device-internal C gives {a:.9}, explicit C gives {b:.9} — \
                             the device is not honouring the integration method"
                        );
                    }
                    // Sanity: the capacitance is actually doing something, so the
                    // assertion above is not passing on two identical flat lines.
                    let (mid, end) = (
                        rd.voltage_at("a", 45e-9).unwrap(),
                        rd.voltage_at("a", 290e-9).unwrap(),
                    );
                    assert!(
                        (mid / end).abs() < 0.6,
                        "{what} {mode:?} variable_step={variable_step}: expected an RC \
                         lag, got {mid} at 45 ns vs {end} settled"
                    );
                }
            }
        }
    }

    /// Any non-zero `K` coupling used to make the variable-step integrator fail
    /// outright — `NoConvergence` after 150 iterations — because its inductor
    /// history update used the standalone `i = G_eq·v + I_hist` and dropped the
    /// mutual term, so the stored current drifted from what was stamped.
    /// `K=0` ran, `K=0.1` did not, and the same netlist was fine fixed-step.
    #[test]
    fn var_step_handles_coupled_inductors() {
        // L/R = 10 µs on each winding.  Primary is pulsed; the secondary sees
        // energy only through the mutual term, so it is a direct probe of it.
        let netlist_for = |k: &str| {
            parse_spice(&format!(
                "* transformer\nV1 a1 0 PULSE(0 1 0 1n 1n 1 2)\n\
                 R1 a1 b1 100\nL1 b1 0 1m\nR2 a2 b2 100\nL2 b2 0 1m\n\
                 K1 L1 L2 {k}\n.options method=be\n.tran 1u 200u\n.end\n"
            ))
            .unwrap()
        };

        for k in ["0.1", "0.5", "0.8", "0.95"] {
            let net = netlist_for(k);
            let reg = DeviceRegistry::new();
            let opts = SimOptions::from_netlist(&net);

            // Must run at all — this is what regressed.
            let var = tran_nr_with_registry_var_opts(&net, 1e-6, 200e-6, &reg, &opts)
                .unwrap_or_else(|e| panic!("k={k} failed under variable-step: {e:?}"));
            // A fine fixed-step run is the reference; BE is first order, so 200 ns
            // on a 10 µs tau is ~20x better resolved than the variable-step start.
            let refr = tran_nr_with_registry_opts(&net, 200e-9, 200e-6, &reg, &opts).unwrap();

            for &t in &[2e-6, 10e-6, 30e-6, 50e-6] {
                let (a, b) = (
                    var.voltage_at("b2", t).unwrap(),
                    refr.voltage_at("b2", t).unwrap(),
                );
                assert!(
                    (a - b).abs() < 0.02,
                    "k={k} t={t:e}: variable-step V(b2)={a:.6}, fine fixed-step={b:.6}"
                );
            }
            // Not vacuous: the mutual term must actually deliver something.
            let peak = var.voltage_at("b2", 2e-6).unwrap();
            let expected = k.parse::<f64>().unwrap() * 0.8; // rough, k-proportional
            assert!(
                peak > 0.5 * expected,
                "k={k}: secondary saw {peak:.6}, expected coupling-proportional transfer"
            );
        }
    }

    #[test]
    fn tran_nr_var_diode_steady_state() {
        // R-D series at DC: var-step should converge to the same OP as fixed-step.
        let netlist_str = "* Diode DC\nVdd a 0 DC 5\nR1 a b 10k\nD1 b 0 myd\n\
             .model myd D (Is=1e-14 N=1)\n.tran 1u 10u\n.end\n";
        let netlist = parse_spice(netlist_str).unwrap();

        let dc = crate::newton::dc_op_nr(&netlist).unwrap();
        let vb_dc = dc.node_voltage("b").unwrap();

        let tr = tran_nr_var(&netlist, 1e-6, 10e-6).unwrap();
        let vb_tran = tr.voltage_at("b", 10e-6).unwrap();
        assert!(
            (vb_tran - vb_dc).abs() < 1e-5,
            "var-step V(b)={vb_tran:.6e}  dc_op={vb_dc:.6e}"
        );
    }

    #[test]
    fn tran_nr_var_gear2_more_accurate_than_be() {
        // RC step response, analytic V(t) = 1 - exp(-t/τ), τ=RC=1ms.
        // Both methods are run with the same LTE tolerance; GEAR-2 is second-
        // order so it should be at least as accurate as BE at the matched τ.
        let src = "* RC\nV1 in 0 PULSE(0 1 0 1n 1n 10m 20m)\n\
                   R1 in out 1k\nC1 out 0 1u\n.tran 200u 4m\n.end\n";
        let net = parse_spice(src).unwrap();
        let mut registry = crate::device_registry::DeviceRegistry::new();
        registry.register_builtin_models(&net.models);

        let opts_be = crate::options::SimOptions {
            method: IntegratorMode::BackwardEuler,
            ..Default::default()
        };
        let r_be = tran_nr_with_registry_var_opts(&net, 200e-6, 4e-3, &registry, &opts_be).unwrap();

        let opts_g = crate::options::SimOptions {
            method: IntegratorMode::Gear,
            ..Default::default()
        };
        let r_g = tran_nr_with_registry_var_opts(&net, 200e-6, 4e-3, &registry, &opts_g).unwrap();

        let tau = 1e-3;
        let exact = |t: f64| 1.0 - (-t / tau).exp();
        let err_at = |r: &TranResult, t: f64| (r.voltage_at("out", t).unwrap() - exact(t)).abs();
        // Sample at 3·τ where both methods are well past the transient.
        let e_be = err_at(&r_be, 3e-3);
        let e_g = err_at(&r_g, 3e-3);
        assert!(
            e_g <= e_be + 1e-4,
            "GEAR-2 error {e_g:.3e} should not exceed BE error {e_be:.3e}"
        );
    }

    #[test]
    fn tran_nr_var_breakpoint_landing() {
        // PULSE with tr=100ns, step=5µs (step >> tr). Without breakpoint insertion the
        // variable-step solver can skip over the rise edge entirely; with it, the solver
        // must land at td+tr and the plateau voltage V(out) at that time must be ≈ v1=5V.
        //
        // Circuit: V1 → R1 → out (no reactive element so V(out) = V1 instantly).
        let netlist = parse_spice(
            "* breakpoint test\nV1 in 0 PULSE(0 5 1u 100n 100n 5u 10u)\nR1 in out 1\n.tran 5u 20u\n.end\n"
        ).unwrap();
        let result = tran_nr_var(&netlist, 5e-6, 20e-6).unwrap();

        // At t = td + tr = 1.1µs the source should be fully risen; V(out) ≈ 5V.
        let v_post_rise = result.voltage_at("out", 1.1e-6).unwrap();
        assert!(
            (v_post_rise - 5.0).abs() < 0.01,
            "V(out) at t=1.1µs should be ≈5V after rise; got {v_post_rise:.4}"
        );
    }

    // ---------- .ic / .nodeset / UIC tests ----------

    #[test]
    fn uic_uses_ic_at_t0() {
        // RC circuit with V1=0 at t=0 (then ramps to 1V via PULSE).  Without
        // UIC, V(out) at t=0 is the DC OP = 0V.  With UIC and `.ic V(out)=0.5`,
        // V(out) at t=0 should be 0.5V.
        let net = fairchild_parser::parse_spice(
            "* uic test\nV1 in 0 PULSE(0 1 1m 1n 1n 100m 200m)\n\
             R1 in out 1k\nC1 out 0 1u\n\
             .ic V(out)=0.5\n\
             .options uic=1\n\
             .tran 10u 100u\n.end\n",
        )
        .unwrap();
        let r = tran_nr_var(&net, 10e-6, 100e-6).unwrap();
        let v0 = r.voltage_at("out", 0.0).unwrap();
        assert!(
            (v0 - 0.5).abs() < 1e-6,
            "UIC: V(out) at t=0 should equal .ic value (0.5), got {v0}"
        );
    }

    #[test]
    fn no_uic_starts_from_dc_op() {
        // Without UIC the t=0 voltage comes from the DC OP, ignoring .ic.
        let net = fairchild_parser::parse_spice(
            "* no uic\nV1 in 0 PULSE(0 1 1m 1n 1n 100m 200m)\n\
             R1 in out 1k\nC1 out 0 1u\n\
             .ic V(out)=0.5\n\
             .tran 10u 100u\n.end\n",
        )
        .unwrap();
        let r = tran_nr_var(&net, 10e-6, 100e-6).unwrap();
        let v0 = r.voltage_at("out", 0.0).unwrap();
        assert!(
            v0.abs() < 1e-6,
            "no UIC: V(out) at t=0 should be DC OP (0V), got {v0}"
        );
    }

    #[test]
    fn tran_nr_fixed_step_errors_on_nonconvergence() {
        // RC driven by a 1V step. With itl4=1 and reltol=1e-30, even a
        // single NR iteration on a diode/RC circuit should fail to converge.
        // We just want to confirm the function returns Err rather than Ok.
        use crate::device_registry::DeviceRegistry;
        use crate::options::SimOptions;
        let netlist = parse_spice(
            "* non-converge test\nVdd a 0 DC 5\nR1 a b 1k\nD1 b 0 myd\n\
             .model myd D (Is=1e-14 N=1)\n.tran 1u 2u\n.end\n",
        )
        .unwrap();
        let registry = {
            let mut r = DeviceRegistry::new();
            r.register_builtin_models(&netlist.models);
            r
        };
        // itl4=1 + reltol=0 is enough to guarantee NR never converges.
        let mut opts = SimOptions::from_netlist(&netlist);
        opts.itl4 = 1;
        opts.reltol = 0.0;
        opts.vntol = 0.0;
        let result = tran_nr_with_registry_opts(&netlist, 1e-6, 2e-6, &registry, &opts);
        assert!(result.is_err(), "expected Err on non-convergence, got Ok");
        let err_str = result.err().unwrap().to_string();
        assert!(
            err_str.contains("did not converge"),
            "unexpected error: {err_str}"
        );
    }

    #[test]
    fn write_nutmeg_tran() {
        let netlist =
            parse_spice("* RC\nV1 in 0 DC 1\nR1 in out 1k\nC1 out 0 1u\n.tran 1u 10u\n.end\n")
                .unwrap();
        let result = tran_be(&netlist, 1e-6, 10e-6);
        let mut buf = Vec::new();
        result.write_nutmeg(&mut buf, "test").unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("Plotname: Transient Analysis"), "plotname: {s}");
        assert!(s.contains("Flags: real"), "flags: {s}");
        assert!(s.contains("time\ttime"), "time var: {s}");
        assert!(s.contains("v(out)\tvoltage"), "v(out): {s}");
        assert!(s.contains("Values:"), "values section: {s}");
    }
}
