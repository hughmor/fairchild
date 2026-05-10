/// Fixed-step Backward Euler (BDF-1) transient integrator.
///
/// At each time step tn = t0 + n·h:
///   1. Build MNA with BE companion models for capacitors and inductors.
///   2. Solve A·x = b (faer LU).
///   3. Update companion state (V_C, I_L) for the next step.
///   4. Store the solution in TranResult.

use indexmap::IndexMap;

use fairchild_parser::{Element, Netlist};

use crate::error::SimError;
use crate::mna::{cap_companion, ind_companion, stamp_netlist, CircuitTopology};
use crate::solver::lu_solve;

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
}

fn interp(xs: &[f64], ys: &[f64], x: f64) -> Option<f64> {
    if xs.is_empty() { return None; }
    if x <= xs[0] { return Some(ys[0]); }
    if x >= *xs.last().unwrap() { return Some(*ys.last().unwrap()); }
    let i = xs.partition_point(|&xi| xi <= x).saturating_sub(1);
    let t = (x - xs[i]) / (xs[i + 1] - xs[i]);
    Some(ys[i] + t * (ys[i + 1] - ys[i]))
}

/// Run a fixed-step Backward Euler transient simulation.
///
/// `step` and `stop` come from the `.tran` directive.
/// Initial conditions: V_C = 0, I_L = 0 for all reactive elements.
/// (Matching ngspice behaviour when no `.ic` is specified and UIC is not used.)
pub fn run_tran(netlist: &Netlist, step: f64, stop: f64) -> Result<TranResult, SimError> {
    let topo = CircuitTopology::build(netlist);

    // Companion state: capacitor name → (G_eq, I_hist), inductor name → (G_eq, I_hist)
    let mut cap_state: IndexMap<String, (f64, f64)> = IndexMap::new();
    let mut ind_state: IndexMap<String, (f64, f64)> = IndexMap::new();

    // Initialise companion models with zero initial conditions.
    for el in &netlist.elements {
        match el {
            Element::Capacitor { name, capacitance, .. } => {
                cap_state.insert(name.clone(), cap_companion(*capacitance, step, 0.0));
            }
            Element::Inductor { name, inductance, .. } => {
                ind_state.insert(name.clone(), ind_companion(*inductance, step, 0.0));
            }
            _ => {}
        }
    }

    // Pre-allocate result storage.
    let n_steps = ((stop / step).ceil() as usize) + 1;
    let mut result = TranResult {
        time: Vec::with_capacity(n_steps),
        node_voltages: topo.node_index.keys()
            .map(|k| (k.clone(), Vec::with_capacity(n_steps)))
            .collect(),
        vsrc_currents: topo.vsrc_index.keys()
            .map(|k| (k.clone(), Vec::with_capacity(n_steps)))
            .collect(),
    };

    let mut t = 0.0_f64;
    loop {
        // Stamp and solve.
        let mat = stamp_netlist(&topo, netlist, t, &cap_state, &ind_state);
        let x = lu_solve(&mat.a, &mat.b)?;

        // Store timepoint.
        result.time.push(t);
        for (name, &idx) in &topo.node_index {
            result.node_voltages.get_mut(name).unwrap().push(x[idx]);
        }
        let n_nodes = topo.n_nodes();
        for (name, &idx) in &topo.vsrc_index {
            result.vsrc_currents.get_mut(name).unwrap().push(x[n_nodes + idx]);
        }

        if t >= stop { break; }
        t = (t + step).min(stop);

        // Update companion models for next step.
        for el in &netlist.elements {
            match el {
                Element::Capacitor { name, pos, neg, capacitance } => {
                    let v_p = topo.node_voltage(pos, &x).unwrap_or(0.0);
                    let v_n = topo.node_voltage(neg, &x).unwrap_or(0.0);
                    let v_c = v_p - v_n;
                    cap_state.insert(name.clone(), cap_companion(*capacitance, step, v_c));
                }
                Element::Inductor { name, pos, neg, inductance } => {
                    // Current through inductor = G_eq * V_L + I_hist (Norton companion).
                    let (g_eq, i_hist) = ind_state[name];
                    let v_p = topo.node_voltage(pos, &x).unwrap_or(0.0);
                    let v_n = topo.node_voltage(neg, &x).unwrap_or(0.0);
                    let i_l = g_eq * (v_p - v_n) + i_hist;
                    ind_state.insert(name.clone(), ind_companion(*inductance, step, i_l));
                }
                _ => {}
            }
        }
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

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
}
