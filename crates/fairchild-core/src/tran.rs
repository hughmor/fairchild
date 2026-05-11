/// Fixed-step transient integrators: Backward Euler (BDF-1) and Trapezoidal Rule (BDF-2 / TR).
///
/// `run_tran` / `run_tran_tr` — linear-only (R, L, C, V, I).
/// `tran_nr` / `tran_nr_tr`   — adds Newton-Raphson for nonlinear devices (diodes, etc.).

use indexmap::IndexMap;

use fairchild_parser::{Element, Netlist};

use crate::device::{EvalFlags, SimContext};
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::mna::{
    cap_companion, cap_companion_be_to_tr, cap_companion_tr_advance,
    ind_companion, ind_companion_be_to_tr, ind_companion_tr_advance,
    stamp_netlist, CircuitTopology,
};
use crate::newton::{build_devices, dc_op_nr_with_registry, GMIN, VMAX, VNTOL, RELTOL, MAX_ITER};
use crate::solver::lu_solve;

/// Transient integration method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntegratorMode {
    /// Backward Euler (BDF-1): first-order, unconditionally stable, single-history.
    BackwardEuler,
    /// Trapezoidal Rule (TR / BDF-2-like): second-order, A-stable, minimal overhead.
    Trapezoidal,
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
        if node == "0" || node == "gnd" { return Some(0.0); }
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
    if xs.is_empty() { return None; }
    if x <= xs[0] { return Some(ys[0]); }
    if x >= *xs.last().unwrap() { return Some(*ys.last().unwrap()); }
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
        node_voltages: topo.node_index.keys()
            .map(|k| (k.clone(), Vec::with_capacity(n_steps)))
            .collect(),
        vsrc_currents: topo.vsrc_index.keys()
            .map(|k| (k.clone(), Vec::with_capacity(n_steps)))
            .collect(),
    };

    // Store the true IC at t=0 before any integration.
    push_timepoint(&mut result, 0.0, &topo, &x0);

    // Initialize companion state from IC.
    let (mut cap_state, mut ind_state) = init_companions(netlist, &topo, step, &x0, mode);

    let mut t = step;
    let mut x = x0;
    let mut first_tr = true;
    loop {
        let mat = stamp_netlist(&topo, netlist, t, &cap_state, &ind_state);
        x = lu_solve(&mat.a, &mat.b)?;
        push_timepoint(&mut result, t, &topo, &x);
        if t >= stop { break; }
        advance_companions(netlist, &topo, step, &x, &mut cap_state, &mut ind_state, mode, first_tr);
        first_tr = false;
        t = (t + step).min(stop);
    }

    Ok(result)
}

// ---------------------------------------------------------------------------
// Helpers shared by run_tran and tran_nr
// ---------------------------------------------------------------------------

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
) -> (IndexMap<String, (f64, f64)>, IndexMap<String, (f64, f64)>) {
    let mut cap_state: IndexMap<String, (f64, f64)> = IndexMap::new();
    let mut ind_state: IndexMap<String, (f64, f64)> = IndexMap::new();
    for el in &netlist.elements {
        match el {
            Element::Capacitor { name, pos, neg, capacitance } => {
                let vc = topo.node_voltage(pos, x).unwrap_or(0.0)
                    - topo.node_voltage(neg, x).unwrap_or(0.0);
                // TR always begins with a BE step for stability at discontinuities.
                cap_state.insert(name.clone(), cap_companion(*capacitance, step, vc));
            }
            Element::Inductor { name, inductance, .. } => {
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
    for el in &netlist.elements {
        match el {
            Element::Capacitor { name, pos, neg, capacitance } => {
                let vc = topo.node_voltage(pos, x).unwrap_or(0.0)
                    - topo.node_voltage(neg, x).unwrap_or(0.0);
                let next = match mode {
                    IntegratorMode::BackwardEuler => cap_companion(*capacitance, step, vc),
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
            Element::Inductor { name, pos, neg, inductance } => {
                let (g_eq, i_hist) = ind_state[name];
                let vl = topo.node_voltage(pos, x).unwrap_or(0.0)
                    - topo.node_voltage(neg, x).unwrap_or(0.0);
                let next = match mode {
                    IntegratorMode::BackwardEuler => {
                        let il = g_eq * vl + i_hist;
                        ind_companion(*inductance, step, il)
                    }
                    IntegratorMode::Trapezoidal if first_tr => {
                        ind_companion_be_to_tr(g_eq, i_hist, vl)
                    }
                    IntegratorMode::Trapezoidal => {
                        ind_companion_tr_advance(g_eq, i_hist, vl)
                    }
                };
                ind_state.insert(name.clone(), next);
            }
            _ => {}
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
pub fn tran_nr_with_registry(
    netlist: &Netlist,
    step: f64,
    stop: f64,
    registry: &DeviceRegistry,
) -> Result<TranResult, SimError> {
    tran_nr_with_registry_mode(netlist, step, stop, registry, IntegratorMode::BackwardEuler)
}

fn tran_nr_with_registry_mode(
    netlist: &Netlist,
    step: f64,
    stop: f64,
    registry: &DeviceRegistry,
    mode: IntegratorMode,
) -> Result<TranResult, SimError> {
    let ctx = SimContext::default();

    let dc = dc_op_nr_with_registry(netlist, registry)?;
    let topo = dc.topo;
    let mut x = dc.x;

    let mut devices = build_devices(netlist, &topo, &ctx, registry)?;

    // Reactive companion state seeded from the DC OP.
    let (mut cap_state, mut ind_state) = init_companions(netlist, &topo, step, &x, mode);

    let n_steps = ((stop / step).ceil() as usize) + 2;
    let mut result = TranResult {
        time: Vec::with_capacity(n_steps),
        node_voltages: topo.node_index.keys()
            .map(|k| (k.clone(), Vec::with_capacity(n_steps)))
            .collect(),
        vsrc_currents: topo.vsrc_index.keys()
            .map(|k| (k.clone(), Vec::with_capacity(n_steps)))
            .collect(),
    };

    // Store t = 0 from DC OP.
    push_timepoint(&mut result, 0.0, &topo, &x);

    let mut t = step;
    let mut first_tr = true;
    loop {
        // --- NR loop for this time step ---
        for _iter in 0..MAX_ITER {
            let mut mat = stamp_netlist(&topo, netlist, t, &cap_state, &ind_state);

            for dev in &mut devices {
                dev.eval(&x, EvalFlags::tran(), &ctx);
                dev.load_residual_tran(&mut mat.b, 1.0);
                dev.load_jacobian_tran(&mut mat, 1.0);
            }

            for i in 0..topo.n_nodes() {
                mat.a[i][i] += GMIN;
            }

            let x_new = lu_solve(&mat.a, &mat.b)?;

            let max_dv = x_new.iter()
                .zip(x.iter())
                .take(topo.n_nodes())
                .map(|(n, o)| (n - o).abs())
                .fold(0.0f64, f64::max);

            let x_next: Vec<f64> = if max_dv > VMAX {
                let scale = VMAX / max_dv;
                x.iter().zip(x_new.iter()).map(|(o, n)| o + scale * (n - o)).collect()
            } else {
                x_new
            };

            let converged = x_next.iter()
                .zip(x.iter())
                .all(|(n, o)| (n - o).abs() < VNTOL + RELTOL * n.abs());

            x = x_next;
            if converged { break; }
        }

        push_timepoint(&mut result, t, &topo, &x);

        if t >= stop { break; }
        advance_companions(netlist, &topo, step, &x, &mut cap_state, &mut ind_state, mode, first_tr);
        first_tr = false;
        t = (t + step).min(stop);
    }

    Ok(result)
}

/// Fixed-step Backward Euler transient using only built-in models from `.model` cards.
pub fn tran_nr(netlist: &Netlist, step: f64, stop: f64) -> Result<TranResult, SimError> {
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_diodes(&netlist.models);
    tran_nr_with_registry(netlist, step, stop, &registry)
}

/// Fixed-step Trapezoidal Rule transient with Newton-Raphson and a pre-built registry.
///
/// Second-order accurate; same interface as `tran_nr_with_registry` but using TR integration.
pub fn tran_nr_with_registry_tr(
    netlist: &Netlist,
    step: f64,
    stop: f64,
    registry: &DeviceRegistry,
) -> Result<TranResult, SimError> {
    tran_nr_with_registry_mode(netlist, step, stop, registry, IntegratorMode::Trapezoidal)
}

/// Fixed-step Trapezoidal Rule transient using only built-in models from `.model` cards.
pub fn tran_nr_tr(netlist: &Netlist, step: f64, stop: f64) -> Result<TranResult, SimError> {
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_diodes(&netlist.models);
    tran_nr_with_registry_tr(netlist, step, stop, &registry)
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
        let v_nr  = r_nr.voltage_at("out", 1e-3).unwrap();
        assert!(
            (v_lin - v_nr).abs() < 1e-4,
            "tran_nr diverges from run_tran at t=1ms: linear={v_lin:.6}  nr={v_nr:.6}"
        );
    }

    #[test]
    fn tran_nr_diode_steady_state() {
        // R-D series, constant V=5V, no reactive elements.
        // tran_nr result at t=1µs must match dc_op_nr within 0.1%.
        let netlist_str =
            "* Diode DC via transient\nVdd a 0 DC 5\nR1 a b 10k\nD1 b 0 myd\n\
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
            "* RC step\nV1 in 0 DC 1\nR1 in out 1k\nC1 out 0 1u\n.tran 200u 5m\n.end\n"
        ).unwrap();
        let h = 200e-6;          // 5 steps per τ — large enough to show BE error
        let exact = 1.0 - (-1.0_f64).exp();  // ≈ 0.6321

        let r_be = run_tran(&netlist, h, 5e-3).unwrap();
        let r_tr = run_tran_tr(&netlist, h, 5e-3).unwrap();

        let v_be = r_be.voltage_at("out", 1e-3).unwrap();
        let v_tr = r_tr.voltage_at("out", 1e-3).unwrap();

        let err_be = (v_be - exact).abs();
        let err_tr = (v_tr - exact).abs();

        assert!(err_tr < err_be, "TR should be more accurate than BE at same step: be_err={err_be:.4e} tr_err={err_tr:.4e}");
        assert!(err_tr < 0.01, "TR error at t=τ should be < 1%: {err_tr:.4e}");
    }

    #[test]
    fn tran_nr_tr_matches_linear_tr() {
        // For a pure RC (no nonlinear device), tran_nr_tr and run_tran_tr should agree.
        let netlist = parse_spice(
            "* RC\nV1 in 0 PULSE(0 1 0 1n 1n 10m 20m)\nR1 in out 1k\nC1 out 0 1u\n.tran 100u 2m\n.end\n"
        ).unwrap();
        let r_lin = run_tran_tr(&netlist, 100e-6, 2e-3).unwrap();
        let r_nr  = tran_nr_tr(&netlist, 100e-6, 2e-3).unwrap();

        let v_lin = r_lin.voltage_at("out", 1e-3).unwrap();
        let v_nr  = r_nr.voltage_at("out", 1e-3).unwrap();
        assert!(
            (v_lin - v_nr).abs() < 1e-3,
            "tran_nr_tr vs run_tran_tr at t=1ms: lin={v_lin:.6} nr={v_nr:.6}"
        );
    }

    #[test]
    fn write_nutmeg_tran() {
        let netlist = parse_spice(
            "* RC\nV1 in 0 DC 1\nR1 in out 1k\nC1 out 0 1u\n.tran 1u 10u\n.end\n",
        ).unwrap();
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
