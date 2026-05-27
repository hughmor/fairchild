/// Transient integrators: fixed-step BE/TR (linear) and variable-step BE+LTE (nonlinear).
///
/// `run_tran` / `run_tran_tr` — linear-only (R, L, C, V, I).
/// `tran_nr` / `tran_nr_tr`   — fixed-step nonlinear.
/// `tran_nr_var`              — variable-step nonlinear with LTE control.
use indexmap::IndexMap;
use std::collections::HashSet;

use fairchild_parser::{Element, Netlist};

use crate::device::{Device, EvalFlags, ReactiveKind};
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::mna::{
    cap_companion, cap_companion_be_to_tr, cap_companion_gear2, cap_companion_tr_advance,
    ind_companion, ind_companion_be_to_tr, ind_companion_gear2, ind_companion_tr_advance,
    stamp_netlist, CircuitTopology, MnaMatrix,
};
use crate::newton::{build_devices, dc_op_nr_with_registry_opts};
use crate::options::SimOptions;
use crate::solver::lu_solve;

/// Transient integration method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntegratorMode {
    /// Backward Euler (BDF-1): first-order, unconditionally stable, single-history.
    BackwardEuler,
    /// Trapezoidal Rule (TR / BDF-2-like): second-order, A-stable, minimal overhead.
    Trapezoidal,
    /// GEAR / BDF-2: second-order, L-stable, two-step history.  First step and
    /// the step after any rejection demote to BE (order control 1↔2).  Applies
    /// to linear L/C companions; device-internal reactive terms remain BDF-1.
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

/// Run a fixed-step Backward Euler transient simulation (linear circuits only).
///
/// `step` and `stop` come from the `.tran` directive.
/// Initial conditions: V_C = 0, I_L = 0 for all reactive elements.
pub fn run_tran(netlist: &Netlist, step: f64, stop: f64) -> Result<TranResult, SimError> {
    run_tran_mode(netlist, step, stop, IntegratorMode::BackwardEuler)
}

/// Run a fixed-step Trapezoidal Rule transient simulation (linear circuits only).
///
/// Second-order accurate; preferred over `run_tran` for smooth waveforms.
pub fn run_tran_tr(netlist: &Netlist, step: f64, stop: f64) -> Result<TranResult, SimError> {
    run_tran_mode(netlist, step, stop, IntegratorMode::Trapezoidal)
}

fn run_tran_mode(
    netlist: &Netlist,
    step: f64,
    stop: f64,
    mode: IntegratorMode,
) -> Result<TranResult, SimError> {
    let topo = CircuitTopology::build(netlist);
    // True initial conditions: all zeros (IC for sources not yet active).
    let x0 = vec![0.0f64; topo.size];

    let n_steps = ((stop / step).ceil() as usize) + 2;
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

    // Store the true IC at t=0 before any integration.
    push_timepoint(&mut result, 0.0, &topo, &x0);

    // Initialize companion state from IC.
    let (mut cap_state, mut ind_state) = init_companions(netlist, &topo, step, &x0, mode);

    let mut t = step;
    let mut x;
    let mut first_tr = true;
    loop {
        let mat = stamp_netlist(&topo, netlist, t, &cap_state, &ind_state);
        x = lu_solve(&mat.a, &mat.b)?;
        push_timepoint(&mut result, t, &topo, &x);
        if t >= stop {
            break;
        }
        advance_companions(
            netlist,
            &topo,
            step,
            &x,
            &mut cap_state,
            &mut ind_state,
            mode,
            first_tr,
        );
        first_tr = false;
        t = (t + step).min(stop);
    }

    Ok(result)
}

// ---------------------------------------------------------------------------
// Helpers shared by run_tran and tran_nr
// ---------------------------------------------------------------------------

/// Companion state maps for capacitors and inductors: (Geq, Ieq) per element name.
type CompanionStatePair = (IndexMap<String, (f64, f64)>, IndexMap<String, (f64, f64)>);

/// Initialise capacitor and inductor companion state from a solution vector.
///
/// For capacitors: V_C_init = V(pos) - V(neg) from the provided x (typically the DC OP).
/// For inductors: I_L_init = 0 (inductors are open-circuit in the DC solver).
fn init_companions(
    netlist: &Netlist,
    topo: &CircuitTopology,
    step: f64,
    x: &[f64],
    _mode: IntegratorMode,
) -> CompanionStatePair {
    let mut cap_state: IndexMap<String, (f64, f64)> = IndexMap::new();
    let mut ind_state: IndexMap<String, (f64, f64)> = IndexMap::new();
    for el in &netlist.elements {
        match el {
            Element::Capacitor {
                name,
                pos,
                neg,
                capacitance,
            } => {
                let vc = topo.node_voltage(pos, x).unwrap_or(0.0)
                    - topo.node_voltage(neg, x).unwrap_or(0.0);
                // TR always begins with a BE step for stability at discontinuities.
                cap_state.insert(name.clone(), cap_companion(*capacitance, step, vc));
            }
            Element::Inductor {
                name, inductance, ..
            } => {
                // TR always begins with a BE step for stability at discontinuities.
                ind_state.insert(name.clone(), ind_companion(*inductance, step, 0.0));
            }
            _ => {}
        }
    }
    (cap_state, ind_state)
}

/// Advance companion state: read post-step voltages/currents from x, update state maps.
///
/// `first_tr` — true only on the first advance when mode is TR; performs the BE→TR
/// transition that seeds the TR with correct initial capacitor/inductor currents.
fn advance_companions(
    netlist: &Netlist,
    topo: &CircuitTopology,
    step: f64,
    x: &[f64],
    cap_state: &mut IndexMap<String, (f64, f64)>,
    ind_state: &mut IndexMap<String, (f64, f64)>,
    mode: IntegratorMode,
    first_tr: bool,
) {
    // Snapshot i_hist values BEFORE the main loop advances them.
    // Needed by the K-element post-pass which must read old history.
    use std::collections::HashMap;
    let i_hist_snapshot: HashMap<String, f64> = ind_state
        .iter()
        .map(|(k, &(_, ih))| (k.clone(), ih))
        .collect();

    for el in &netlist.elements {
        match el {
            Element::Capacitor {
                name,
                pos,
                neg,
                capacitance,
            } => {
                let vc = topo.node_voltage(pos, x).unwrap_or(0.0)
                    - topo.node_voltage(neg, x).unwrap_or(0.0);
                let next = match mode {
                    IntegratorMode::BackwardEuler | IntegratorMode::Gear => {
                        cap_companion(*capacitance, step, vc)
                    }
                    IntegratorMode::Trapezoidal if first_tr => {
                        let (g_eq_be, i_hist_be) = cap_state[name];
                        cap_companion_be_to_tr(g_eq_be, i_hist_be, vc)
                    }
                    IntegratorMode::Trapezoidal => {
                        let (g_eq, i_hist_old) = cap_state[name];
                        cap_companion_tr_advance(g_eq, i_hist_old, vc)
                    }
                };
                cap_state.insert(name.clone(), next);
            }
            Element::Inductor {
                name,
                pos,
                neg,
                inductance,
            } => {
                let (g_eq, i_hist) = ind_state[name];
                let vl = topo.node_voltage(pos, x).unwrap_or(0.0)
                    - topo.node_voltage(neg, x).unwrap_or(0.0);
                let next = match mode {
                    IntegratorMode::BackwardEuler | IntegratorMode::Gear => {
                        let il = g_eq * vl + i_hist;
                        ind_companion(*inductance, step, il)
                    }
                    IntegratorMode::Trapezoidal if first_tr => {
                        ind_companion_be_to_tr(g_eq, i_hist, vl)
                    }
                    IntegratorMode::Trapezoidal => ind_companion_tr_advance(g_eq, i_hist, vl),
                };
                ind_state.insert(name.clone(), next);
            }
            _ => {}
        }
    }

    // Post-pass: correct ind_state for coupled inductor pairs (K elements).
    //
    // The main loop above used the standalone formula IL = G_eq*VL + I_hist.
    // For coupled pairs we need:
    //   IL1 = G11*VL1 + G12*VL2 + I_hist1_old
    //   IL2 = G12*VL1 + G22*VL2 + I_hist2_old
    //
    // We use i_hist_snapshot (old values) to avoid reading partially-updated state.
    //
    // K correction is only valid for BE/Gear.  In TR mode the companion update
    // uses the trapezoidal formula (ind_companion_tr_advance) whose state is
    // incompatible with the BE-derived G11/G22/G12 stamps.  Silently skipping
    // avoids wrong numbers; a proper TR K correction is a future TODO.
    if matches!(mode, IntegratorMode::BackwardEuler | IntegratorMode::Gear) {
        let mut ind_vals_local: HashMap<String, f64> = HashMap::new();
        for el in &netlist.elements {
            if let Element::Inductor {
                name, inductance, ..
            } = el
            {
                ind_vals_local.insert(name.clone(), *inductance);
            }
        }

        for el in &netlist.elements {
            if let Element::CoupledInductors {
                l1, l2, coupling, ..
            } = el
            {
                let Some(&val1) = ind_vals_local.get(l1) else {
                    continue;
                };
                let Some(&val2) = ind_vals_local.get(l2) else {
                    continue;
                };
                // g_eq1 after the main loop advance = h/val1 (unchanged by the main loop
                // since ind_companion preserves the g_eq formula).
                let Some(&(g_eq1, _)) = ind_state.get(l1) else {
                    continue;
                };
                let Some(&(g_eq2, _)) = ind_state.get(l2) else {
                    continue;
                };
                let i_hist1 = i_hist_snapshot.get(l1).copied().unwrap_or(0.0);
                let i_hist2 = i_hist_snapshot.get(l2).copied().unwrap_or(0.0);

                // Find terminal nodes for voltage lookup.
                let vl1 = {
                    let (pos1, neg1) = find_inductor_terminals_by_name(netlist, l1);
                    topo.node_voltage(pos1, x).unwrap_or(0.0)
                        - topo.node_voltage(neg1, x).unwrap_or(0.0)
                };
                let vl2 = {
                    let (pos2, neg2) = find_inductor_terminals_by_name(netlist, l2);
                    topo.node_voltage(pos2, x).unwrap_or(0.0)
                        - topo.node_voltage(neg2, x).unwrap_or(0.0)
                };

                let k = *coupling;
                let m = k * (val1 * val2).sqrt();
                let det = val1 * val2 - m * m;
                if det.abs() < 1e-40 {
                    continue;
                }

                // h = g_eq * L  (g_eq = h/L)
                let h = g_eq1 * val1;
                let g11 = h * val2 / det;
                let g22 = h * val1 / det;
                let g12 = -h * m / det;

                let il1 = g11 * vl1 + g12 * vl2 + i_hist1;
                let il2 = g12 * vl1 + g22 * vl2 + i_hist2;

                // Overwrite with corrected currents.
                ind_state.insert(l1.clone(), ind_companion(val1, step, il1));
                ind_state.insert(l2.clone(), ind_companion(val2, step, il2));

                // Suppress unused variable warnings for g_eq2.
                let _ = g_eq2;
            }
        }
    } // end if BackwardEuler | Gear
}

/// Find the (pos, neg) terminal names for a named inductor — tran.rs variant.
fn find_inductor_terminals_by_name<'a>(netlist: &'a Netlist, name: &str) -> (&'a str, &'a str) {
    for el in &netlist.elements {
        if let Element::Inductor {
            name: n, pos, neg, ..
        } = el
        {
            if n == name {
                return (pos, neg);
            }
        }
    }
    panic!("coupled inductor '{name}' not found in netlist")
}

/// Seed device-internal reactive-branch companion state from the DC OP.
///
/// For each device, queries `reactive_branches()` once and computes the
/// initial (G_eq, I_hist) using the DC-OP voltage across the branch's
/// (pos, neg) terminals.  Returned as `[dev_idx][branch_idx] -> (G_eq, I_hist)`.
fn init_device_reactive_state(
    devices: &[Box<dyn Device>],
    x: &[f64],
    step: f64,
    _mode: IntegratorMode,
) -> Vec<Vec<(f64, f64)>> {
    let mut out = Vec::with_capacity(devices.len());
    for dev in devices {
        let branches = dev.reactive_branches();
        let mut dev_states = Vec::with_capacity(branches.len());
        for br in &branches {
            let v = match (br.pos, br.neg) {
                (Some(p), Some(n)) => x[p] - x[n],
                (Some(p), None) => x[p],
                (None, Some(n)) => -x[n],
                (None, None) => 0.0,
            };
            // TR always begins with a BE step for stability at discontinuities
            // — mirrors the built-in Element::Capacitor handling.
            let companion = match br.kind {
                ReactiveKind::Capacitor => cap_companion(br.value, step, v),
                ReactiveKind::Inductor => ind_companion(br.value, step, 0.0),
            };
            dev_states.push(companion);
        }
        out.push(dev_states);
    }
    out
}

/// Stamp the companion-model contributions of every device-internal
/// reactive branch into `mat`.  Called inside the NR loop, AFTER each
/// device's `load_jacobian_tran` (so the device's eval cache reflects the
/// current iterate before the integrator re-queries the value).
fn stamp_device_reactive_companions(
    devices: &[Box<dyn Device>],
    state: &[Vec<(f64, f64)>],
    mat: &mut MnaMatrix,
    step: f64,
) {
    for (dev_idx, dev) in devices.iter().enumerate() {
        let branches = dev.reactive_branches();
        for (br_idx, br) in branches.iter().enumerate() {
            // Re-compute companion params using the CURRENT (per-NR-iter)
            // device value, but reuse I_hist from the previous timestep.
            let (_g_old, i_hist) = state[dev_idx][br_idx];
            let g_eq = match br.kind {
                ReactiveKind::Capacitor => br.value / step,
                ReactiveKind::Inductor => step / br.value.max(1e-30),
            };
            // Stamp G_eq between pos and neg (resistor pattern) and inject
            // ±I_hist into the corresponding KCL rows.  For an inductor the
            // current is i = G_eq · v + I_hist (history adds, doesn't subtract).
            let i_sign = match br.kind {
                ReactiveKind::Capacitor => 1.0,
                ReactiveKind::Inductor => -1.0,
            };
            if let Some(p) = br.pos {
                mat.a[p][p] += g_eq;
                if let Some(n) = br.neg {
                    mat.a[p][n] -= g_eq;
                }
                mat.b[p] += i_sign * i_hist;
            }
            if let Some(n) = br.neg {
                mat.a[n][n] += g_eq;
                if let Some(p) = br.pos {
                    mat.a[n][p] -= g_eq;
                }
                mat.b[n] -= i_sign * i_hist;
            }
        }
    }
}

/// Advance the device-internal reactive companion state after a successful
/// timestep.  Reads V_C / I_L from the converged `x` and computes the next
/// (G_eq, I_hist) using the device's reported value at the new operating
/// point (already updated via `commit_timestep` + the next NR loop's eval).
fn advance_device_reactive_state(
    devices: &[Box<dyn Device>],
    x: &[f64],
    state: &mut [Vec<(f64, f64)>],
    step: f64,
) {
    for (dev_idx, dev) in devices.iter().enumerate() {
        let branches = dev.reactive_branches();
        for (br_idx, br) in branches.iter().enumerate() {
            let v = match (br.pos, br.neg) {
                (Some(p), Some(n)) => x[p] - x[n],
                (Some(p), None) => x[p],
                (None, Some(n)) => -x[n],
                (None, None) => 0.0,
            };
            let next = match br.kind {
                ReactiveKind::Capacitor => cap_companion(br.value, step, v),
                ReactiveKind::Inductor => {
                    // i_L = G_eq · v_L + I_hist (history), then form next companion.
                    let (g_eq, i_hist) = state[dev_idx][br_idx];
                    let i_l = g_eq * v + i_hist;
                    ind_companion(br.value, step, i_l)
                }
            };
            state[dev_idx][br_idx] = next;
        }
    }
}

/// Append one time-point to a TranResult.
fn push_timepoint(result: &mut TranResult, t: f64, topo: &CircuitTopology, x: &[f64]) {
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
pub fn tran_nr_with_registry_opts(
    netlist: &Netlist,
    step: f64,
    stop: f64,
    registry: &DeviceRegistry,
    opts: &SimOptions,
) -> Result<TranResult, SimError> {
    // Sanity check fires here when UIC bypasses DC OP; the non-UIC path
    // gets it via the dc_op call below (and the check is cheap enough
    // that we don't bother de-duplicating).
    if opts.sanity_check && opts.uic {
        crate::sanity::check_netlist_sanity(netlist);
    }
    crate::connectivity::check_connectivity(netlist)?;
    let ctx = opts.sim_context();
    let mode = opts.method;

    // With UIC: skip DC OP, seed x from `.ic` (or 0 where unspecified).
    // Without UIC (the default): use DC OP as t=0 condition.
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
        // DC OP already allocated extras; reuse its topology so the row
        // layout (and matrix size) stays consistent through the transient.
        (dc.topo, dc.x)
    };
    let mut topo = topo;

    let mut devices = build_devices(netlist, &mut topo, &ctx, registry)?;
    // After build_devices, topo.size is final.  Pad the initial x vector
    // to match — when we came via UIC, x was sized to the pre-allocation
    // topology; OSDI internal nodes get a zero initial guess.
    x.resize(topo.size, 0.0);
    // Topology is fixed across the whole transient — build the linear
    // solver once and reuse for every NR iter.
    let solver = opts.linear_solver(topo.size);

    // Seed x_tprev from DC OP (or UIC initial conditions) so reactive
    // history is defined before the first step.
    for dev in &mut devices {
        dev.commit_timestep(&x);
    }

    // Honour opts.max_step as an upper bound on the step size.
    let step = step.min(opts.max_step);

    // Reactive companion state seeded from the DC OP.
    let (mut cap_state, mut ind_state) = init_companions(netlist, &topo, step, &x, mode);
    // Device-internal reactive branches (e.g., bias-dependent C_j on a
    // depletion-mode PN-PS).  One companion-state pair per device per
    // declared branch.  Indexed as dev_reactive_state[dev_idx][branch_idx].
    let mut dev_reactive_state: Vec<Vec<(f64, f64)>> =
        init_device_reactive_state(&devices, &x, step, mode);

    let n_steps = ((stop / step).ceil() as usize) + 2;
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
    push_timepoint(&mut result, 0.0, &topo, &x);

    // Cached factorisation: in transient, the sparsity pattern is fixed
    // across the entire run (devices don't appear / disappear, timestep
    // changes only scale values via `α = 1/h`).  One symbolic factor-
    // isation up front, `klu_refactor` (or faer-sparse value-only
    // rebuild) on every NR iteration of every timestep.
    let mut fact: Option<Box<dyn crate::solver::Factorisation>> = None;

    // One MnaMatrix reused across every NR iteration of every timestep —
    // for a 10 k-step transient with ~3 NR iters each this skips 30 k ×
    // (N+1) heap allocations.
    let mut mat = crate::mna::MnaMatrix::zeros(topo.size);

    let mut t = step;
    let mut first_tr = true;
    loop {
        // --- NR loop for this time step ---
        let alpha = 1.0 / step;
        let mut step_converged = false;
        for _iter in 0..opts.itl4 {
            crate::mna::stamp_netlist_in_place(&mut mat, &topo, netlist, t, &cap_state, &ind_state);

            for dev in &mut devices {
                dev.eval(&x, EvalFlags::tran(), &ctx);
                dev.load_residual_tran(&mut mat.b, alpha);
                dev.load_jacobian_tran(&mut mat, alpha);
            }
            // Stamp integrator-managed reactive companions for every
            // device-declared linear reactive branch (uses the device's
            // current bias-dependent value AND the history from the
            // previous accepted timestep).
            stamp_device_reactive_companions(&devices, &dev_reactive_state, &mut mat, step);

            for i in 0..topo.n_nodes() {
                mat.a[i][i] += opts.gmin;
            }
            // gmin on OSDI internal-node rows (see newton.rs for rationale).
            let vsrc_end = topo.n_nodes() + topo.vsrc_index.len();
            for i in vsrc_end..topo.size {
                mat.a[i][i] += opts.gmin;
            }

            let x_new = if let Some(f) = fact.as_mut() {
                f.refactor_and_solve(&mat.a, &mat.b)?
            } else {
                let mut f = solver.factorise(&mat.a)?;
                let r = f.refactor_and_solve(&mat.a, &mat.b)?;
                fact = Some(f);
                r
            };

            let max_dv = x_new
                .iter()
                .zip(x.iter())
                .take(topo.n_nodes())
                .map(|(n, o)| (n - o).abs())
                .fold(0.0f64, f64::max);

            let x_next: Vec<f64> = if max_dv > opts.vmax {
                let scale = opts.vmax / max_dv;
                x.iter()
                    .zip(x_new.iter())
                    .map(|(o, n)| o + scale * (n - o))
                    .collect()
            } else {
                x_new
            };

            let converged = x_next
                .iter()
                .zip(x.iter())
                .all(|(n, o)| (n - o).abs() < opts.vntol + opts.reltol * n.abs());

            x = x_next;
            if converged {
                step_converged = true;
                break;
            }
        }

        if !step_converged {
            return Err(SimError::NoConvergence { iters: opts.itl4 });
        }

        push_timepoint(&mut result, t, &topo, &x);

        for dev in &mut devices {
            dev.commit_timestep(&x);
        }

        if t >= stop {
            break;
        }
        advance_companions(
            netlist,
            &topo,
            step,
            &x,
            &mut cap_state,
            &mut ind_state,
            mode,
            first_tr,
        );
        advance_device_reactive_state(&devices, &x, &mut dev_reactive_state, step);
        first_tr = false;
        t = (t + step).min(stop);
    }

    Ok(result)
}

/// Fixed-step Backward Euler transient using only built-in models from `.model` cards.
pub fn tran_nr(netlist: &Netlist, step: f64, stop: f64) -> Result<TranResult, SimError> {
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_diodes(&netlist.models);
    registry.register_builtin_mosfets(&netlist.models);
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
    registry.register_builtin_diodes(&netlist.models);
    registry.register_builtin_mosfets(&netlist.models);
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
    if opts.sanity_check && opts.uic {
        crate::sanity::check_netlist_sanity(netlist);
    }
    crate::connectivity::check_connectivity(netlist)?;
    let ctx = opts.sim_context();
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
    let mut devices = build_devices(netlist, &mut topo, &ctx, registry)?;
    // Pad x for any OSDI internal-node rows allocated by build_devices.
    x.resize(topo.size, 0.0);
    let solver = opts.linear_solver(topo.size);

    let n_nodes = topo.n_nodes();
    let h_min = step * 1e-6;

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

    // Raw physical state: cap voltage and inductor current seeded from DC OP.
    // `cap_v_prev2` / `ind_i_prev2` hold the state two timesteps back (only
    // populated once at least one step has been accepted) for GEAR-2.
    let mut cap_v: IndexMap<String, f64> = IndexMap::new();
    let mut cap_v_prev2: IndexMap<String, f64> = IndexMap::new();
    let mut ind_i: IndexMap<String, f64> = IndexMap::new();
    let mut ind_i_prev2: IndexMap<String, f64> = IndexMap::new();
    for el in &netlist.elements {
        match el {
            Element::Capacitor { name, pos, neg, .. } => {
                let vc = topo.node_voltage(pos, &x).unwrap_or(0.0)
                    - topo.node_voltage(neg, &x).unwrap_or(0.0);
                cap_v.insert(name.clone(), vc);
            }
            Element::Inductor { name, .. } => {
                ind_i.insert(name.clone(), 0.0);
            }
            _ => {}
        }
    }
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
    let mut mat = crate::mna::MnaMatrix::zeros(topo.size);

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

        // Order control: GEAR-2 needs two accepted steps of history.  Demote
        // to BE on the first two steps, after any rejection (history stale),
        // and on extreme step ratios where BDF-2 would amplify noise.
        let history_ready = h_prev_accepted > 0.0
            && cap_v.keys().all(|n| cap_v_prev2.contains_key(n))
            && ind_i.keys().all(|n| ind_i_prev2.contains_key(n));
        let step_ratio = if h_prev_accepted > 0.0 {
            h_actual / h_prev_accepted
        } else {
            1.0
        };
        let use_gear2 = matches!(opts.method, IntegratorMode::Gear)
            && history_ready
            && consecutive_rejects == 0
            && step_ratio > 0.25
            && step_ratio < 4.0;

        // Build companion maps for h_actual from the stored raw state.
        let mut cap_state: IndexMap<String, (f64, f64)> = IndexMap::new();
        let mut ind_state: IndexMap<String, (f64, f64)> = IndexMap::new();
        for el in &netlist.elements {
            match el {
                Element::Capacitor {
                    name, capacitance, ..
                } => {
                    if let Some(&vc) = cap_v.get(name) {
                        let stamp = if use_gear2 {
                            let vc_prev2 = cap_v_prev2.get(name).copied().unwrap_or(vc);
                            cap_companion_gear2(
                                *capacitance,
                                h_actual,
                                h_prev_accepted,
                                vc,
                                vc_prev2,
                            )
                        } else {
                            cap_companion(*capacitance, h_actual, vc)
                        };
                        cap_state.insert(name.clone(), stamp);
                    }
                }
                Element::Inductor {
                    name, inductance, ..
                } => {
                    if let Some(&il) = ind_i.get(name) {
                        let stamp = if use_gear2 {
                            let il_prev2 = ind_i_prev2.get(name).copied().unwrap_or(il);
                            ind_companion_gear2(
                                *inductance,
                                h_actual,
                                h_prev_accepted,
                                il,
                                il_prev2,
                            )
                        } else {
                            ind_companion(*inductance, h_actual, il)
                        };
                        ind_state.insert(name.clone(), stamp);
                    }
                }
                _ => {}
            }
        }

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
        let mut x_try = x_pred.clone();
        let mut nr_converged = false;

        for _iter in 0..opts.itl4 {
            crate::mna::stamp_netlist_in_place(
                &mut mat, &topo, netlist, t_next, &cap_state, &ind_state,
            );

            for dev in devices.iter_mut() {
                dev.eval(&x_try, EvalFlags::tran(), &ctx);
                dev.load_residual_tran(&mut mat.b, alpha);
                dev.load_jacobian_tran(&mut mat, alpha);
            }

            for i in 0..n_nodes {
                mat.a[i][i] += opts.gmin;
            }
            // gmin on OSDI internal-node rows (see newton.rs for rationale).
            let vsrc_end = n_nodes + topo.vsrc_index.len();
            for i in vsrc_end..topo.size {
                mat.a[i][i] += opts.gmin;
            }

            let x_new = if let Some(f) = fact.as_mut() {
                f.refactor_and_solve(&mat.a, &mat.b)?
            } else {
                let mut f = solver.factorise(&mat.a)?;
                let r = f.refactor_and_solve(&mat.a, &mat.b)?;
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

            let converged = x_next
                .iter()
                .zip(x_try.iter())
                .all(|(n, o)| (n - o).abs() < opts.vntol + opts.reltol * n.abs());

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
                .map(|(_, (xc, xp))| (xc - xp).abs() * 0.5 / (opts.vntol + opts.reltol * xc.abs()))
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

            // Update raw companion state from accepted solution.
            // Shift the BDF-2 history (prev2 ← prev) BEFORE writing the new
            // prev value so the next step has a valid two-step window.
            for el in &netlist.elements {
                match el {
                    Element::Capacitor { name, pos, neg, .. } => {
                        if let Some(&vc_prev) = cap_v.get(name) {
                            cap_v_prev2.insert(name.clone(), vc_prev);
                        }
                        let vc = topo.node_voltage(pos, &x).unwrap_or(0.0)
                            - topo.node_voltage(neg, &x).unwrap_or(0.0);
                        cap_v.insert(name.clone(), vc);
                    }
                    Element::Inductor {
                        name,
                        pos,
                        neg,
                        inductance,
                    } => {
                        let vl = topo.node_voltage(pos, &x).unwrap_or(0.0)
                            - topo.node_voltage(neg, &x).unwrap_or(0.0);
                        let il_prev = ind_i.get(name).copied().unwrap_or(0.0);
                        // Closed-form current update: i_new = G_eq·v_new + I_hist
                        // where (G_eq, I_hist) is the companion we just stamped.
                        let il_new = if let Some(&(g_eq, i_hist)) = ind_state.get(name) {
                            g_eq * vl + i_hist
                        } else {
                            il_prev + (h_actual / inductance) * vl
                        };
                        ind_i_prev2.insert(name.clone(), il_prev);
                        ind_i.insert(name.clone(), il_new);
                    }
                    _ => {}
                }
            }
            h_prev_accepted = h_actual;

            for dev in &mut devices {
                dev.commit_timestep(&x);
            }

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

    Ok(result)
}

/// Variable-step BE + LTE transient using only built-in models from `.model` cards.
pub fn tran_nr_var(netlist: &Netlist, step: f64, stop: f64) -> Result<TranResult, SimError> {
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_diodes(&netlist.models);
    registry.register_builtin_mosfets(&netlist.models);
    tran_nr_with_registry_var(netlist, step, stop, &registry)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

    // ---------- tran_nr tests ----------

    #[test]
    fn tran_nr_matches_linear_for_rc() {
        // Pure RC (no diode): tran_nr and run_tran should agree within 0.1%.
        let netlist = parse_spice(
            "* RC step\nV1 in 0 PULSE(0 1 0 1n 1n 10m 20m)\nR1 in out 1k\nC1 out 0 1u\n.tran 1u 2m\n.end\n"
        ).unwrap();
        let r_linear = run_tran(&netlist, 1e-6, 2e-3).unwrap();
        let r_nr = tran_nr(&netlist, 1e-6, 2e-3).unwrap();

        let v_lin = r_linear.voltage_at("out", 1e-3).unwrap();
        let v_nr = r_nr.voltage_at("out", 1e-3).unwrap();
        assert!(
            (v_lin - v_nr).abs() < 1e-4,
            "tran_nr diverges from run_tran at t=1ms: linear={v_lin:.6}  nr={v_nr:.6}"
        );
    }

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

    // ---------- run_tran regression tests ----------

    #[test]
    fn rc_step_response_shape() {
        // R=1k C=1µF τ=1ms, step to 1V at t=0. V(out) should be ~0.632 at t=τ.
        let netlist = parse_spice(
            "* RC step\nV1 in 0 PULSE(0 1 0 1n 1n 10m 20m)\nR1 in out 1k\nC1 out 0 1u\n.tran 1u 5m\n.end\n"
        ).unwrap();

        let result = run_tran(&netlist, 1e-6, 5e-3).unwrap();

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

        let result = run_tran(&netlist, 1e-6, 5e-3).unwrap();

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
        let netlist = parse_spice(
            "* RC step\nV1 in 0 DC 1\nR1 in out 1k\nC1 out 0 1u\n.tran 200u 5m\n.end\n",
        )
        .unwrap();
        let h = 200e-6; // 5 steps per τ — large enough to show BE error
        let exact = 1.0 - (-1.0_f64).exp(); // ≈ 0.6321

        let r_be = run_tran(&netlist, h, 5e-3).unwrap();
        let r_tr = run_tran_tr(&netlist, h, 5e-3).unwrap();

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

    #[test]
    fn tran_nr_tr_matches_linear_tr() {
        // For a pure RC (no nonlinear device), tran_nr_tr and run_tran_tr should agree.
        let netlist = parse_spice(
            "* RC\nV1 in 0 PULSE(0 1 0 1n 1n 10m 20m)\nR1 in out 1k\nC1 out 0 1u\n.tran 100u 2m\n.end\n"
        ).unwrap();
        let r_lin = run_tran_tr(&netlist, 100e-6, 2e-3).unwrap();
        let r_nr = tran_nr_tr(&netlist, 100e-6, 2e-3).unwrap();

        let v_lin = r_lin.voltage_at("out", 1e-3).unwrap();
        let v_nr = r_nr.voltage_at("out", 1e-3).unwrap();
        assert!(
            (v_lin - v_nr).abs() < 1e-3,
            "tran_nr_tr vs run_tran_tr at t=1ms: lin={v_lin:.6} nr={v_nr:.6}"
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
        registry.register_builtin_diodes(&net.models);

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
            r.register_builtin_diodes(&netlist.models);
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
        let result = run_tran(&netlist, 1e-6, 10e-6).unwrap();
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
