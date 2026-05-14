//! Pre-flight connectivity check: every non-ground node must have a DC path
//! to ground (resistor, inductor, voltage source, or device terminal pair).
//!
//! Catches the most common "circuit looks fine but matrix is singular" class
//! of errors before LU returns NaN.  Floating sub-networks come back with
//! a clean `SimError::FloatingNodes` diagnostic listing the orphan nodes.

use std::collections::{HashMap, HashSet, VecDeque};

use fairchild_parser::{Element, Netlist};

use crate::error::SimError;

/// Verify that every node in the netlist has at least one DC-connected path
/// to ground ("0").  Returns `Ok(())` if the circuit is fully connected, else
/// `Err(SimError::FloatingNodes { nodes })` listing the orphans.
///
/// Connections considered:
///   - R, L:                       both terminals are DC-connected
///   - VoltageSource:              both terminals are DC-connected (short at DC)
///   - CurrentSource, Capacitor:   open at DC — NO connection contributed
///   - Diode, Mosfet, XOsdi:       all device terminals are mutually connected
///     (conservative: assumes the device provides some finite conductance path
///     between every pair of its terminals).  This avoids false positives on
///     active circuits where the only DC path is through a transistor.
///
/// Ground is "0".  If a netlist contains no non-ground nodes, returns Ok.
pub fn check_connectivity(netlist: &Netlist) -> Result<(), SimError> {
    // Collect every node mentioned by any element.
    let mut nodes: HashSet<String> = HashSet::new();
    let mut record = |n: &str, set: &mut HashSet<String>| {
        if !n.is_empty() {
            set.insert(n.to_string());
        }
    };

    // Build union-find by merging connected groups.
    let mut adj: HashMap<String, HashSet<String>> = HashMap::new();
    let mut add_edge = |u: &str, v: &str, m: &mut HashMap<String, HashSet<String>>| {
        m.entry(u.to_string()).or_default().insert(v.to_string());
        m.entry(v.to_string()).or_default().insert(u.to_string());
    };

    for el in &netlist.elements {
        match el {
            Element::Resistor { pos, neg, .. }
            | Element::Inductor { pos, neg, .. }
            | Element::VoltageSource { pos, neg, .. } => {
                record(pos, &mut nodes); record(neg, &mut nodes);
                add_edge(pos, neg, &mut adj);
            }
            // Capacitor and current source are open at DC; their nodes still
            // need to be reachable through some *other* path.  Just record the
            // nodes so they participate in the orphan check.
            Element::Capacitor { pos, neg, .. }
            | Element::CurrentSource { pos, neg, .. } => {
                record(pos, &mut nodes); record(neg, &mut nodes);
            }
            Element::Diode { anode, cathode, .. } => {
                record(anode, &mut nodes); record(cathode, &mut nodes);
                add_edge(anode, cathode, &mut adj);
            }
            Element::Mosfet { drain, gate, source, bulk, .. } => {
                for n in &[drain, gate, source, bulk] { record(n, &mut nodes); }
                // Treat every pair as connected (conservative).
                let terms = [drain.as_str(), gate.as_str(), source.as_str(), bulk.as_str()];
                for i in 0..terms.len() {
                    for j in (i + 1)..terms.len() {
                        add_edge(terms[i], terms[j], &mut adj);
                    }
                }
            }
            Element::XOsdi { nets, .. } => {
                // OSDI device models (including the photonic Norton-equivalent
                // library) self-stamp their own diagonals — every terminal has
                // a finite DC conductance to ground via the model's internals.
                // Treat each XOsdi terminal as connected to ground directly so
                // the check focuses on genuinely orphan R-L-C-V islands.
                for n in nets {
                    record(n, &mut nodes);
                    add_edge(n, "0", &mut adj);
                }
            }
            Element::Behavioral { pos, neg, .. } => {
                // B-element provides a finite DC stamp between (pos, neg):
                // V= form is an aux row (≈ V-source short), I= form is a
                // current source (open at DC), but in either case the
                // expression Jacobian connects every referenced node to
                // (pos, neg).  Treat as a short for connectivity purposes.
                record(pos, &mut nodes); record(neg, &mut nodes);
                add_edge(pos, neg, &mut adj);
            }
        }
    }

    if nodes.is_empty() {
        return Ok(());
    }

    // BFS from ground ("0").
    let mut visited: HashSet<String> = HashSet::new();
    let mut queue: VecDeque<String> = VecDeque::new();
    visited.insert("0".to_string());
    queue.push_back("0".to_string());

    while let Some(u) = queue.pop_front() {
        if let Some(neighbours) = adj.get(&u) {
            for v in neighbours {
                if visited.insert(v.clone()) {
                    queue.push_back(v.clone());
                }
            }
        }
    }

    // Any node mentioned in the circuit but not visited from ground is floating.
    let mut floating: Vec<String> = nodes.into_iter()
        .filter(|n| n != "0" && !visited.contains(n))
        .collect();
    if floating.is_empty() {
        Ok(())
    } else {
        floating.sort();
        Err(SimError::FloatingNodes { nodes: floating })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

    #[test]
    fn connected_divider_passes() {
        let net = parse_spice(
            "* divider\nV1 in 0 DC 1\nR1 in out 1k\nR2 out 0 1k\n.op\n.end\n"
        ).unwrap();
        check_connectivity(&net).unwrap();
    }

    #[test]
    fn isolated_node_caught() {
        // 'floater' is only on a capacitor (open at DC) → unreachable from ground.
        let net = parse_spice(
            "* float\nV1 in 0 DC 1\nR1 in out 1k\nR2 out 0 1k\n\
             C1 floater 0 1u\n.op\n.end\n"
        ).unwrap();
        let err = check_connectivity(&net).unwrap_err();
        match err {
            SimError::FloatingNodes { nodes } => {
                assert!(nodes.contains(&"floater".to_string()), "nodes={nodes:?}");
            }
            other => panic!("expected FloatingNodes, got {other:?}"),
        }
    }

    #[test]
    fn capacitor_island_is_floating() {
        // Two completely disconnected sub-circuits.  The C1→C2 island has
        // no DC path to ground at all.
        let net = parse_spice(
            "* split\nV1 in 0 DC 1\nR1 in 0 1k\n\
             C1 a b 1u\nC2 b 0 1u\n.op\n.end\n"
        ).unwrap();
        let err = check_connectivity(&net).unwrap_err();
        match err {
            SimError::FloatingNodes { nodes } => {
                assert!(nodes.contains(&"a".to_string()) || nodes.contains(&"b".to_string()),
                    "expected 'a' or 'b' as floating, got {nodes:?}");
            }
            other => panic!("expected FloatingNodes, got {other:?}"),
        }
    }

    #[test]
    fn diode_provides_dc_path() {
        // D1 connects b to ground (R1 connects a to b, V1 connects a to ground).
        let net = parse_spice(
            "* rd\nV1 a 0 DC 1\nR1 a b 1k\nD1 b 0 myd\n\
             .model myd D (Is=1e-14 N=1)\n.op\n.end\n"
        ).unwrap();
        check_connectivity(&net).unwrap();
    }
}
