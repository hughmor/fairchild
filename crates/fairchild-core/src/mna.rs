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
            Element::Diode { .. } => {
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
            Element::VoltageSource { name, pos, neg, .. } => {
                add_node(&mut node_index, pos);
                add_node(&mut node_index, neg);
                let n = vsrc_index.len();
                vsrc_index.entry(name.clone()).or_insert(n);
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

// ---------------------------------------------------------------------------
// Legacy compatibility wrapper used by DC solver and existing tests
// ---------------------------------------------------------------------------

/// Simple DC-only MNA system (no C/L companion, sources at DC value).
pub struct MnaSystem {
    pub a: Vec<Vec<f64>>,
    pub b: Vec<f64>,
    pub node_index: IndexMap<String, usize>,
    pub vsrc_index: IndexMap<String, usize>,
    pub size: usize,
}

impl MnaSystem {
    pub fn build(netlist: &Netlist) -> Result<MnaSystem, SimError> {
        let topo = CircuitTopology::build(netlist);
        let empty: IndexMap<String, (f64, f64)> = IndexMap::new();
        let mat = stamp_netlist(&topo, netlist, 0.0, &empty, &empty);
        Ok(MnaSystem {
            a: mat.a,
            b: mat.b,
            node_index: topo.node_index,
            vsrc_index: topo.vsrc_index,
            size: topo.size,
        })
    }

    pub fn node_voltage(&self, node: &str, x: &[f64]) -> Result<f64, SimError> {
        if node == "0" || node == "gnd" { return Ok(0.0); }
        self.node_index.get(node)
            .map(|&i| x[i])
            .ok_or_else(|| SimError::UnknownNode(node.to_string()))
    }

    pub fn vsrc_current(&self, vsrc_name: &str, x: &[f64]) -> Result<f64, SimError> {
        let n_nodes = self.node_index.len();
        self.vsrc_index.get(vsrc_name)
            .map(|&i| x[n_nodes + i])
            .ok_or_else(|| SimError::UnknownNode(vsrc_name.to_string()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

    #[test]
    fn voltage_divider_matrix_size() {
        let net = parse_spice(
            "* divider\nV1 in 0 1.0\nR1 in mid 1000\nR2 mid 0 1000\n.op\n.end\n",
        )
        .unwrap();
        let sys = MnaSystem::build(&net).unwrap();
        // 2 nodes (in, mid) + 1 vsource = 3×3
        assert_eq!(sys.size, 3);
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
