/// Modified Nodal Analysis (MNA) matrix assembler.
///
/// Matrix layout (0-based):
///   rows/cols [0 .. n_nodes)          → node voltages (ground excluded)
///   rows/cols [n_nodes .. n_nodes+m)  → voltage-source branch currents
///
/// Ground node "0" is eliminated from the matrix.

use fairchild_parser::{Element, Netlist};
use indexmap::IndexMap;

use crate::error::SimError;

/// The topology/index maps for a circuit — built once, reused across time steps.
#[derive(Clone)]
pub struct CircuitTopology {
    /// Map from node name → matrix row/col (ground excluded).
    pub node_index: IndexMap<String, usize>,
    /// Map from voltage source name → aux current index within A/b.
    pub vsrc_index: IndexMap<String, usize>,
    /// Total matrix dimension = n_nodes + n_vsources.
    pub size: usize,
}

impl CircuitTopology {
    pub fn build(netlist: &Netlist) -> CircuitTopology {
        let (node_index, vsrc_index) = index_circuit(netlist);
        let size = node_index.len() + vsrc_index.len();
        CircuitTopology { node_index, vsrc_index, size }
    }

    pub fn n_nodes(&self) -> usize { self.node_index.len() }

    /// Retrieve a node voltage from a solution vector.
    pub fn node_voltage(&self, node: &str, x: &[f64]) -> Result<f64, SimError> {
        if node == "0" || node == "gnd" { return Ok(0.0); }
        self.node_index.get(node)
            .map(|&i| x[i])
            .ok_or_else(|| SimError::UnknownNode(node.to_string()))
    }

    /// Retrieve a voltage-source branch current from a solution vector.
    pub fn vsrc_current(&self, vsrc_name: &str, x: &[f64]) -> Result<f64, SimError> {
        self.vsrc_index.get(vsrc_name)
            .map(|&i| x[self.n_nodes() + i])
            .ok_or_else(|| SimError::UnknownNode(vsrc_name.to_string()))
    }
}

/// An assembled MNA system: A·x = b.
pub struct MnaMatrix {
    pub a: Vec<Vec<f64>>,
    pub b: Vec<f64>,
}

impl MnaMatrix {
    fn new(size: usize) -> Self {
        MnaMatrix {
            a: vec![vec![0.0f64; size]; size],
            b: vec![0.0f64; size],
        }
    }
}

// ---------------------------------------------------------------------------
// DC stamp (time-independent elements, or companion equivalents for C/L)
// ---------------------------------------------------------------------------

/// Like `stamp_netlist` but scales all independent source amplitudes by `source_scale`.
///
/// Used by source-stepping homotopy: call with source_scale ∈ [0,1] to ramp
/// voltage/current sources from zero to their nominal values.
pub fn stamp_netlist_scaled(
    topo: &CircuitTopology,
    netlist: &Netlist,
    source_scale: f64,
    cap_state: &IndexMap<String, (f64, f64)>,
    ind_state: &IndexMap<String, (f64, f64)>,
) -> MnaMatrix {
    let mut mat = stamp_netlist(topo, netlist, 0.0, cap_state, ind_state);
    if (source_scale - 1.0).abs() < 1e-15 {
        return mat;
    }
    // Re-scale source contributions: voltage sources stamp into b[vi], current sources into b[node].
    // Rather than re-stamping from scratch, rescale: b[vi] already has waveform.at(0) from
    // stamp_netlist (t=0 = DC value). We just scale b entries for source rows.
    let n_nodes = topo.n_nodes();
    for (name, &vi_idx) in &topo.vsrc_index {
        let vi = n_nodes + vi_idx;
        // Find the source's waveform value and rescale.
        if let Some(el) = netlist.elements.iter().find(|e| {
            matches!(e, Element::VoltageSource { name: n, .. } if n == name)
        }) {
            if let Element::VoltageSource { waveform, .. } = el {
                let v_full = waveform.at(0.0);
                // The stamp put v_full into b[vi]; replace with v_full * source_scale.
                mat.b[vi] = mat.b[vi] - v_full + v_full * source_scale;
            }
        }
    }
    for el in &netlist.elements {
        if let Element::CurrentSource { pos, neg, waveform, .. } = el {
            let i_full = waveform.at(0.0);
            if i_full == 0.0 { continue; }
            let delta = i_full * (source_scale - 1.0);
            // stamp_current_source: b[pos] -= i, b[neg] += i. Scale correction = delta.
            if let Some(&p) = topo.node_index.get(pos) { mat.b[p] -= delta; }
            if let Some(&n) = topo.node_index.get(neg) { mat.b[n] += delta; }
        }
    }
    mat
}

/// Stamp the full MNA system for DC operating-point or a single transient step.
///
/// `t`           — current simulation time (used for PULSE sources)
/// `cap_state`   — map from capacitor name → (G_eq, I_hist) BE companion pair
/// `ind_state`   — map from inductor name  → (G_eq, I_hist) BE companion pair
pub fn stamp_netlist(
    topo: &CircuitTopology,
    netlist: &Netlist,
    t: f64,
    cap_state: &IndexMap<String, (f64, f64)>,
    ind_state: &IndexMap<String, (f64, f64)>,
) -> MnaMatrix {
    let mut mat = MnaMatrix::new(topo.size);
    let n_nodes = topo.n_nodes();

    for el in &netlist.elements {
        match el {
            Element::Resistor { pos, neg, resistance, .. } => {
                stamp_conductance(&mut mat.a, &topo.node_index, pos, neg, 1.0 / resistance);
            }
            Element::Capacitor { name, pos, neg, .. } => {
                if let Some(&(g_eq, i_hist)) = cap_state.get(name) {
                    stamp_conductance(&mut mat.a, &topo.node_index, pos, neg, g_eq);
                    // BE companion: KCL at pos gives b[pos] += I_hist.
                    // stamp_current_source(neg, pos, v) adds v to b[pos].
                    stamp_current_source(&mut mat.b, &topo.node_index, neg, pos, i_hist);
                }
                // Capacitor absent from cap_state = open circuit (correct for DC OP).
            }
            Element::Inductor { name, pos, neg, .. } => {
                if let Some(&(g_eq, i_hist)) = ind_state.get(name) {
                    stamp_conductance(&mut mat.a, &topo.node_index, pos, neg, g_eq);
                    // BE companion: KCL at pos gives b[pos] -= I_hist.
                    // stamp_current_source(pos, neg, v) subtracts v from b[pos].
                    stamp_current_source(&mut mat.b, &topo.node_index, pos, neg, i_hist);
                }
                // If not in ind_state (no transient yet), inductor = short circuit at DC.
                // For DC OP we treat L as wire — but that would make a voltage source loop.
                // In practice we skip the inductor in DC OP (treated as ideal short via V=0 source).
                // For now we leave it out of DC (open circuit), which is correct for initial tran.
            }
            Element::VoltageSource { name, pos, neg, waveform } => {
                let vi = n_nodes + topo.vsrc_index[name];
                stamp_vsource(&mut mat.a, &mut mat.b, &topo.node_index, pos, neg, vi, waveform.at(t));
            }
            Element::CurrentSource { pos, neg, waveform, .. } => {
                // SPICE: current flows from n+ through source to n- → subtract from n+, add to n-.
                stamp_current_source(&mut mat.b, &topo.node_index, pos, neg, waveform.at(t));
            }
            Element::Diode { .. } | Element::Mosfet { .. }
            | Element::XOsdi { .. } | Element::Behavioral { .. } => {
                // Nonlinear; stamped by the Device trait inside the Newton-Raphson loop.
            }
        }
    }

    mat
}

// ---------------------------------------------------------------------------
// Stamp primitives
// ---------------------------------------------------------------------------

/// Add a conductance G between pos and neg (same as resistor with R=1/G).
pub fn stamp_conductance(
    a: &mut Vec<Vec<f64>>,
    idx: &IndexMap<String, usize>,
    pos: &str,
    neg: &str,
    g: f64,
) {
    if let Some(&p) = idx.get(pos) {
        a[p][p] += g;
        if let Some(&n) = idx.get(neg) {
            a[p][n] -= g;
            a[n][p] -= g;
        }
    }
    if let Some(&n) = idx.get(neg) {
        a[n][n] += g;
    }
}

/// Stamp a voltage source at aux row `vi`: V(pos) - V(neg) = value.
fn stamp_vsource(
    a: &mut Vec<Vec<f64>>,
    b: &mut Vec<f64>,
    idx: &IndexMap<String, usize>,
    pos: &str,
    neg: &str,
    vi: usize,
    value: f64,
) {
    if let Some(&p) = idx.get(pos) {
        a[p][vi] += 1.0;
        a[vi][p] += 1.0;
    }
    if let Some(&n) = idx.get(neg) {
        a[n][vi] -= 1.0;
        a[vi][n] -= 1.0;
    }
    b[vi] += value;
}

/// Stamp a current source.
/// SPICE convention: positive current flows from pos through source to neg.
///   → b[pos] -= value,  b[neg] += value
fn stamp_current_source(
    b: &mut Vec<f64>,
    idx: &IndexMap<String, usize>,
    pos: &str,
    neg: &str,
    value: f64,
) {
    if let Some(&p) = idx.get(pos) { b[p] -= value; }
    if let Some(&n) = idx.get(neg) { b[n] += value; }
}

// ---------------------------------------------------------------------------
// Index builder
// ---------------------------------------------------------------------------

fn index_circuit(netlist: &Netlist) -> (IndexMap<String, usize>, IndexMap<String, usize>) {
    let mut node_index: IndexMap<String, usize> = IndexMap::new();
    let mut vsrc_index: IndexMap<String, usize> = IndexMap::new();

    for el in &netlist.elements {
        match el {
            Element::Resistor { pos, neg, .. }
            | Element::Capacitor { pos, neg, .. }
            | Element::Inductor { pos, neg, .. }
            | Element::CurrentSource { pos, neg, .. } => {
                add_node(&mut node_index, pos);
                add_node(&mut node_index, neg);
            }
            Element::Diode { anode, cathode, .. } => {
                add_node(&mut node_index, anode);
                add_node(&mut node_index, cathode);
            }
            Element::Mosfet { drain, gate, source, bulk, .. } => {
                add_node(&mut node_index, drain);
                add_node(&mut node_index, gate);
                add_node(&mut node_index, source);
                add_node(&mut node_index, bulk);
            }
            Element::VoltageSource { name, pos, neg, .. } => {
                add_node(&mut node_index, pos);
                add_node(&mut node_index, neg);
                let n = vsrc_index.len();
                vsrc_index.entry(name.clone()).or_insert(n);
            }
            Element::XOsdi { nets, .. } => {
                for net in nets {
                    add_node(&mut node_index, net);
                }
            }
            Element::Behavioral { name, pos, neg, kind, .. } => {
                add_node(&mut node_index, pos);
                add_node(&mut node_index, neg);
                if *kind == fairchild_parser::BehavioralKind::Voltage {
                    let n = vsrc_index.len();
                    vsrc_index.entry(name.clone()).or_insert(n);
                }
            }
        }
    }

    (node_index, vsrc_index)
}

fn add_node(map: &mut IndexMap<String, usize>, node: &str) {
    if node != "0" {
        let n = map.len();
        map.entry(node.to_string()).or_insert(n);
    }
}

/// Backward Euler companion state for a capacitor at one time step.
/// Returns (G_eq, I_hist) where G_eq = C/h and I_hist = G_eq * V_prev.
pub fn cap_companion(capacitance: f64, h: f64, v_prev: f64) -> (f64, f64) {
    let g = capacitance / h;
    (g, g * v_prev)
}

/// Backward Euler companion state for an inductor at one time step.
/// Returns (G_eq, I_hist) where G_eq = h/L and I_hist = I_L_prev.
pub fn ind_companion(inductance: f64, h: f64, i_prev: f64) -> (f64, f64) {
    let g = h / inductance;
    (g, i_prev)
}

/// Trapezoidal Rule (TR) companion initial state for a capacitor.
/// G_eq = 2C/h.  I_hist = G_eq * V_prev (assumes I_C(0) = 0 at DC OP).
pub fn cap_companion_tr(capacitance: f64, h: f64, v_prev: f64) -> (f64, f64) {
    let g = 2.0 * capacitance / h;
    (g, g * v_prev)
}

/// Advance the TR capacitor companion after a solved timestep.
///
/// Given the stored state `(G_eq, I_hist_old)` and the new capacitor voltage `v_c_new`,
/// returns the companion for the *next* step.
/// Derivation: I_hist_new = 2·G_eq·V_C_new − I_hist_old.
pub fn cap_companion_tr_advance(g_eq: f64, i_hist_old: f64, v_c_new: f64) -> (f64, f64) {
    (g_eq, 2.0 * g_eq * v_c_new - i_hist_old)
}

/// Transition a capacitor companion from BE→TR after the first BE step.
///
/// The BE companion `(g_eq_be, i_hist_be)` encodes G_be = C/h, I_hist_be = C/h·V_C_prev.
/// After solving for v_c_new, the capacitor current was I_C = G_be·V_C_new − I_hist_be.
/// Returns the TR companion for the *next* step with G_eq = 2·G_be.
pub fn cap_companion_be_to_tr(g_eq_be: f64, i_hist_be: f64, v_c_new: f64) -> (f64, f64) {
    let i_c = g_eq_be * v_c_new - i_hist_be;
    let g_eq_tr = 2.0 * g_eq_be;
    (g_eq_tr, g_eq_tr * v_c_new + i_c)
}

/// Trapezoidal Rule companion initial state for an inductor.
/// G_eq = h/(2L).  I_hist = I_L_prev (= 0 when starting from DC OP).
pub fn ind_companion_tr(inductance: f64, h: f64, i_prev: f64) -> (f64, f64) {
    let g = h / (2.0 * inductance);
    (g, i_prev)
}

/// Advance the TR inductor companion after a solved timestep.
///
/// Given the stored state `(G_eq, I_hist_old)` and new inductor voltage `v_l_new`,
/// returns the companion for the next step.
/// Derivation: I_hist_new = 2·G_eq·V_L_new + I_hist_old.
pub fn ind_companion_tr_advance(g_eq: f64, i_hist_old: f64, v_l_new: f64) -> (f64, f64) {
    (g_eq, 2.0 * g_eq * v_l_new + i_hist_old)
}

/// Transition an inductor companion from BE→TR after the first BE step.
///
/// The inductor current during the BE step was I_L = G_be·V_L_new + I_hist_be.
/// Returns the TR companion for the next step with G_eq = G_be/2.
pub fn ind_companion_be_to_tr(g_eq_be: f64, i_hist_be: f64, v_l_new: f64) -> (f64, f64) {
    let i_l = g_eq_be * v_l_new + i_hist_be;
    let g_eq_tr = g_eq_be / 2.0;
    (g_eq_tr, i_l + g_eq_tr * v_l_new)
}

/// GEAR / BDF-2 companion for a capacitor with non-uniform step.
///
/// Two-step BDF: \dot{v}_{n+1} ≈ α·v_{n+1} − β₁·v_n + β₂·v_{n-1}
/// with ρ = h_n / h_{n-1}, α = (1+2ρ)/(h(1+ρ)), β₁ = (1+ρ)/h, β₂ = ρ²/(h(1+ρ)).
///
/// Returns the Norton (G_eq, I_hist) such that i_C = G_eq·v − I_hist.
pub fn cap_companion_gear2(
    capacitance: f64,
    h: f64,
    h_prev: f64,
    v_prev: f64,
    v_prev2: f64,
) -> (f64, f64) {
    let rho = h / h_prev;
    let denom = h * (1.0 + rho);
    let g_eq = capacitance * (1.0 + 2.0 * rho) / denom;
    let i_hist = capacitance * ((1.0 + rho) / h * v_prev - (rho * rho) / denom * v_prev2);
    (g_eq, i_hist)
}

/// GEAR / BDF-2 companion for an inductor with non-uniform step.
///
/// For v_L = L·di/dt, the BDF-2 Norton form is:
///   i_{n+1} = G_eq·v_{n+1} + I_hist
/// with G_eq = h(1+ρ)/(L(1+2ρ)),
///      I_hist = ((1+ρ)²/(1+2ρ))·i_n − (ρ²/(1+2ρ))·i_{n-1}.
pub fn ind_companion_gear2(
    inductance: f64,
    h: f64,
    h_prev: f64,
    i_prev: f64,
    i_prev2: f64,
) -> (f64, f64) {
    let rho = h / h_prev;
    let denom_a = 1.0 + 2.0 * rho;
    let g_eq = h * (1.0 + rho) / (inductance * denom_a);
    let i_hist = ((1.0 + rho).powi(2) / denom_a) * i_prev - (rho * rho / denom_a) * i_prev2;
    (g_eq, i_hist)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

    #[test]
    fn voltage_divider_topology_size() {
        let net = parse_spice(
            "* divider\nV1 in 0 1.0\nR1 in mid 1000\nR2 mid 0 1000\n.op\n.end\n",
        )
        .unwrap();
        let topo = CircuitTopology::build(&net);
        // 2 nodes (in, mid) + 1 vsource = 3×3
        assert_eq!(topo.size, 3);
    }

    #[test]
    fn rc_circuit_topology() {
        let net = parse_spice(
            "* RC\nV1 in 0 1.0\nR1 in out 1k\nC1 out 0 1u\n.op\n.end\n",
        )
        .unwrap();
        let topo = CircuitTopology::build(&net);
        // nodes: in, out (2) + 1 vsource = 3
        assert_eq!(topo.size, 3);
    }
}
