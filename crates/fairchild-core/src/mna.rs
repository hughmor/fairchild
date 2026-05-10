/// Modified Nodal Analysis matrix assembler.
///
/// Matrix layout (all indices are 0-based into the MNA matrix):
///   rows/cols [0 .. n_nodes-1)     → node voltages (ground excluded)
///   rows/cols [n_nodes-1 .. n_nodes-1+m) → voltage-source branch currents
///
/// Ground is always node "0" and is eliminated from the matrix.

use fairchild_parser::{Element, Netlist};
use indexmap::IndexMap;

use crate::error::SimError;

/// Completed MNA system: A * x = b.
pub struct MnaSystem {
    /// Matrix A (dense for now; sparse via faer-rs used in solver).
    pub a: Vec<Vec<f64>>,
    /// RHS vector b.
    pub b: Vec<f64>,
    /// Map from node name to matrix row/col index (ground excluded).
    pub node_index: IndexMap<String, usize>,
    /// Map from voltage source name to aux current index (within a/b).
    pub vsrc_index: IndexMap<String, usize>,
    /// Total matrix dimension = (n_nodes - 1) + n_vsources.
    pub size: usize,
}

impl MnaSystem {
    /// Build the MNA system from a parsed netlist.
    pub fn build(netlist: &Netlist) -> Result<MnaSystem, SimError> {
        let (node_index, vsrc_index) = index_circuit(netlist);
        let n_nodes = node_index.len();       // excludes ground
        let n_vsrc = vsrc_index.len();
        let size = n_nodes + n_vsrc;

        let mut a = vec![vec![0.0f64; size]; size];
        let mut b = vec![0.0f64; size];

        for el in &netlist.elements {
            match el {
                Element::Resistor { pos, neg, resistance, .. } => {
                    stamp_resistor(&mut a, &node_index, pos, neg, *resistance);
                }
                Element::VoltageSource { name, pos, neg, dc } => {
                    let vi = n_nodes + vsrc_index[name];
                    stamp_vsource(&mut a, &mut b, &node_index, pos, neg, vi, *dc);
                }
                Element::CurrentSource { pos, neg, dc, .. } => {
                    stamp_isource(&mut b, &node_index, pos, neg, *dc);
                }
            }
        }

        Ok(MnaSystem { a, b, node_index, vsrc_index, size })
    }

    /// Look up a solved node voltage from the solution vector.
    pub fn node_voltage(&self, node: &str, x: &[f64]) -> Result<f64, SimError> {
        if node == "0" || node == "gnd" {
            return Ok(0.0);
        }
        self.node_index.get(node)
            .map(|&i| x[i])
            .ok_or_else(|| SimError::UnknownNode(node.to_string()))
    }

    /// Look up a solved voltage-source branch current from the solution vector.
    pub fn vsrc_current(&self, vsrc_name: &str, x: &[f64]) -> Result<f64, SimError> {
        let n_nodes = self.node_index.len();
        self.vsrc_index.get(vsrc_name)
            .map(|&i| x[n_nodes + i])
            .ok_or_else(|| SimError::UnknownNode(vsrc_name.to_string()))
    }
}

/// Assign integer indices to every non-ground node and every voltage source.
fn index_circuit(netlist: &Netlist) -> (IndexMap<String, usize>, IndexMap<String, usize>) {
    let mut node_index: IndexMap<String, usize> = IndexMap::new();
    let mut vsrc_index: IndexMap<String, usize> = IndexMap::new();

    for el in &netlist.elements {
        match el {
            Element::Resistor { pos, neg, .. }
            | Element::CurrentSource { pos, neg, .. } => {
                add_node(&mut node_index, pos);
                add_node(&mut node_index, neg);
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

/// Stamp a resistor with conductance G = 1/R into the A matrix.
fn stamp_resistor(
    a: &mut Vec<Vec<f64>>,
    idx: &IndexMap<String, usize>,
    pos: &str,
    neg: &str,
    r: f64,
) {
    let g = 1.0 / r;
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

/// Stamp a voltage source (KVL row + KCL coupling).
fn stamp_vsource(
    a: &mut Vec<Vec<f64>>,
    b: &mut Vec<f64>,
    idx: &IndexMap<String, usize>,
    pos: &str,
    neg: &str,
    vi: usize,
    dc: f64,
) {
    if let Some(&p) = idx.get(pos) {
        a[p][vi] += 1.0;
        a[vi][p] += 1.0;
    }
    if let Some(&n) = idx.get(neg) {
        a[n][vi] -= 1.0;
        a[vi][n] -= 1.0;
    }
    b[vi] += dc;
}

/// Stamp an independent current source.
///
/// SPICE convention for `I name n+ n-`: positive current flows from n+ through
/// the source to n-. So current LEAVES n+ and ENTERS n-.
///   b[n+] -= dc    (removes dc from n+ KCL)
///   b[n-] += dc    (injects dc into n- KCL)
fn stamp_isource(
    b: &mut Vec<f64>,
    idx: &IndexMap<String, usize>,
    pos: &str,
    neg: &str,
    dc: f64,
) {
    if let Some(&p) = idx.get(pos) {
        b[p] -= dc;
    }
    if let Some(&n) = idx.get(neg) {
        b[n] += dc;
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
        // 2 nodes (in, mid) + 1 vsource = 3x3 matrix
        assert_eq!(sys.size, 3);
    }
}
