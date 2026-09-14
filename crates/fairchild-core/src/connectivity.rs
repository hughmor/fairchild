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
    let record = |n: &str, set: &mut HashSet<String>| {
        if !n.is_empty() {
            set.insert(n.to_string());
        }
    };

    // Build union-find by merging connected groups.
    let mut adj: HashMap<String, HashSet<String>> = HashMap::new();
    let add_edge = |u: &str, v: &str, m: &mut HashMap<String, HashSet<String>>| {
        m.entry(u.to_string()).or_default().insert(v.to_string());
        m.entry(v.to_string()).or_default().insert(u.to_string());
    };

    for el in &netlist.elements {
        match el {
            Element::Resistor { pos, neg, .. }
            | Element::Inductor { pos, neg, .. }
            | Element::VoltageSource { pos, neg, .. } => {
                record(pos, &mut nodes);
                record(neg, &mut nodes);
                add_edge(pos, neg, &mut adj);
            }
            // Capacitor and current source are open at DC; their nodes still
            // need to be reachable through some *other* path.  Just record the
            // nodes so they participate in the orphan check.
            Element::Capacitor { pos, neg, .. } | Element::CurrentSource { pos, neg, .. } => {
                record(pos, &mut nodes);
                record(neg, &mut nodes);
            }
            Element::Diode { anode, cathode, .. } => {
                record(anode, &mut nodes);
                record(cathode, &mut nodes);
                add_edge(anode, cathode, &mut adj);
            }
            Element::Mosfet {
                drain,
                gate,
                source,
                bulk,
                ..
            } => {
                for n in &[drain, gate, source, bulk] {
                    record(n, &mut nodes);
                }
                // Treat every pair as connected (conservative).
                let terms = [
                    drain.as_str(),
                    gate.as_str(),
                    source.as_str(),
                    bulk.as_str(),
                ];
                for i in 0..terms.len() {
                    for j in (i + 1)..terms.len() {
                        add_edge(terms[i], terms[j], &mut adj);
                    }
                }
            }
            Element::Bjt {
                collector,
                base,
                emitter,
                substrate,
                ..
            } => {
                for n in &[collector, base, emitter, substrate] {
                    record(n, &mut nodes);
                }
                let terms = [
                    collector.as_str(),
                    base.as_str(),
                    emitter.as_str(),
                    substrate.as_str(),
                ];
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
                record(pos, &mut nodes);
                record(neg, &mut nodes);
                add_edge(pos, neg, &mut adj);
            }
            Element::CoupledInductors { .. } => {
                // K element only affects transient; L1/L2 terminals are
                // already handled by their Inductor elements.
            }
            Element::VoltageSwitch {
                pos,
                neg,
                ctrl_pos,
                ctrl_neg,
                ..
            } => {
                // Conducting either way (RON or ROFF), so the switched pair is
                // always joined for the orphan-island check. The control pair
                // is a sense input and joins nothing.
                for n in &[pos, neg, ctrl_pos, ctrl_neg] {
                    record(n, &mut nodes);
                }
                add_edge(pos, neg, &mut adj);
            }
            Element::CurrentSwitch { pos, neg, .. } => {
                record(pos, &mut nodes);
                record(neg, &mut nodes);
                add_edge(pos, neg, &mut adj);
            }
            Element::TransmissionLine {
                a_pos,
                a_neg,
                b_pos,
                b_neg,
                ..
            } => {
                // At DC a lossless line is an ideal through-connection
                // (A+↔B+, A−↔B−), so treat the matching conductors as joined
                // for the orphan-island check.
                for n in &[a_pos, a_neg, b_pos, b_neg] {
                    record(n, &mut nodes);
                }
                add_edge(a_pos, b_pos, &mut adj);
                add_edge(a_neg, b_neg, &mut adj);
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
    let mut floating: Vec<String> = nodes
        .into_iter()
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
        let net = parse_spice("* divider\nV1 in 0 DC 1\nR1 in out 1k\nR2 out 0 1k\n.op\n").unwrap();
        check_connectivity(&net).unwrap();
    }

    #[test]
    fn isolated_node_caught() {
        // 'floater' is only on a capacitor (open at DC) → unreachable from ground.
        let net = parse_spice(
            "* float\nV1 in 0 DC 1\nR1 in out 1k\nR2 out 0 1k\n\
             C1 floater 0 1u\n.op\n",
        )
        .unwrap();
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
             C1 a b 1u\nC2 b 0 1u\n.op\n",
        )
        .unwrap();
        let err = check_connectivity(&net).unwrap_err();
        match err {
            SimError::FloatingNodes { nodes } => {
                assert!(
                    nodes.contains(&"a".to_string()) || nodes.contains(&"b".to_string()),
                    "expected 'a' or 'b' as floating, got {nodes:?}"
                );
            }
            other => panic!("expected FloatingNodes, got {other:?}"),
        }
    }

    #[test]
    fn diode_provides_dc_path() {
        // D1 connects b to ground (R1 connects a to b, V1 connects a to ground).
        let net = parse_spice(
            "* rd\nV1 a 0 DC 1\nR1 a b 1k\nD1 b 0 myd\n\
             .model myd D (Is=1e-14 N=1)\n.op\n",
        )
        .unwrap();
        check_connectivity(&net).unwrap();
    }
}

// ───────────────────────────────────────────────────────────────────────────
// Optical loops
// ───────────────────────────────────────────────────────────────────────────

/// Optical ports that lie on a closed path, if there is one.
///
/// A resonator here is not a device — a ring is two couplers and a waveguide,
/// composed in the netlist, which is the design. So "does this deck contain a
/// cavity" can only be answered from the topology, and the same answer covers a
/// Sagnac loop, a ring-assisted MZI, and a feedback path that goes out through
/// a detector and back in through a driver.
///
/// # Bipartite, deliberately
///
/// The graph is *elements on one side, optical ports on the other*, with an
/// edge wherever an element touches a port. A cycle in that graph means two
/// distinct elements share two distinct ports, which is a loop.
///
/// Connecting an element's own ports to each other instead would report a
/// 4-port coupler as a cycle all by itself — three ports make a triangle — and
/// every deck with a coupler in it would be a resonator. The bipartite form
/// makes a single element a star, which has no cycle whatever its port count.
pub fn optical_loop(netlist: &Netlist) -> Option<Vec<String>> {
    // Port name → the wires it expands to, so an element's net list can be read
    // back as the ports it touches.
    let mut wire_to_port: HashMap<String, String> = HashMap::new();
    for bp in netlist.bundle_ports.iter().filter(|b| b.is_optical()) {
        for w in bp.all_wires() {
            wire_to_port.insert(w.to_lowercase(), bp.name.to_lowercase());
        }
    }
    if wire_to_port.is_empty() {
        return None;
    }

    // One node per element and one per port, in a single index space.
    let mut port_id: HashMap<String, usize> = HashMap::new();
    let mut port_name: Vec<String> = Vec::new();
    let mut edges: Vec<(usize, usize)> = Vec::new();
    let mut elem_count = 0usize;
    for el in &netlist.elements {
        let Element::XOsdi { nets, .. } = el else {
            continue;
        };
        let mut touched: Vec<usize> = nets
            .iter()
            .filter_map(|n| wire_to_port.get(&n.to_lowercase()))
            .map(|p| {
                let next = port_id.len();
                *port_id.entry(p.clone()).or_insert_with(|| {
                    port_name.push(p.clone());
                    next
                })
            })
            .collect();
        touched.sort_unstable();
        touched.dedup();
        if touched.len() < 2 {
            // A laser or a detector is an endpoint; it cannot close anything.
            continue;
        }
        let elem = elem_count;
        elem_count += 1;
        for p in touched {
            edges.push((elem, p));
        }
    }
    if elem_count == 0 {
        return None;
    }

    // Union-find over `elements ∪ ports`: an edge joining two already-joined
    // nodes closes a cycle.
    let n = elem_count + port_id.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    let mut on_loop: Vec<String> = Vec::new();
    for (elem, port) in edges {
        let (a, b) = (
            find(&mut parent, elem),
            find(&mut parent, elem_count + port),
        );
        if a == b {
            // This edge closes a cycle. Naming the port it closes on is enough
            // to find the loop by eye, and cheaper than walking it out.
            on_loop.push(port_name[port].clone());
        } else {
            parent[a] = b;
        }
    }
    on_loop.sort();
    on_loop.dedup();
    (!on_loop.is_empty()).then_some(on_loop)
}

/// Warn when a transient runs an optical cavity with its delays switched off.
///
/// Without a group delay a loop's round trip is instantaneous, so the field
/// takes its steady-state value for whatever phase the drive has *right now*.
/// The wavelength-domain answer is untouched — the propagation phase is in the
/// transfer, so a `.dc` sweep over λ still shows the resonance — but the cavity
/// has no photon lifetime, and for a high-Q ring modulator that lifetime is
/// often the bandwidth limit.
///
/// It warns rather than refusing, and it warns on the *physics* rather than on
/// the option. Delays are off by default, so refusing an explicit
/// `optical_delay=0` while permitting the identical default would be
/// incoherent — the two produce the same circuit. Saying so once, and letting
/// an explicit setting acknowledge it, is the version that means something
/// (#123).
pub fn warn_if_cavity_without_delay(netlist: &Netlist, opts: &crate::options::SimOptions) {
    if !opts.optical_delay.should_warn_about_cavities() {
        return;
    }
    let Some(ports) = optical_loop(netlist) else {
        return;
    };
    crate::warn_user!(
        "this deck has a closed optical path (through {}), and optical group \
         delays are off — so its round trip is instantaneous and the cavity has \
         no photon lifetime. Wavelength sweeps are unaffected; a transient \
         response is not. Set `.options waveguide_delay=1` to model it, or \
         `.options optical_delay=0` to say you meant this",
        ports.join(", ")
    );
}

#[cfg(test)]
mod optical_loop_tests {
    use super::optical_loop;
    use fairchild_parser::parse_spice;

    const RING: &str = "* all-pass ring\n\
        .optical_port lin\n\
        .optical_port thru\n\
        .optical_port ra\n\
        .optical_port rb\n\
        Xlas lin fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
        Xc lin ra thru rb fc_dcoupler kappa_L=0.3\n\
        Xwg rb ra fc_waveguide L_um=100 n_g=4.2\n\
        Xpd thru det 0 fc_photodetector responsivity=0.9\n\
        Rl det 0 1k\n";

    const STRAIGHT: &str = "* straight link\n\
        .optical_port lin\n\
        .optical_port mid\n\
        Xlas lin fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
        Xwg lin mid fc_waveguide L_um=100 n_g=4.2\n\
        Xpd mid det 0 fc_photodetector responsivity=0.9\n\
        Rl det 0 1k\n";

    #[test]
    fn a_ring_is_a_loop_and_a_link_is_not() {
        let ring = parse_spice(RING).expect("parse");
        assert!(
            optical_loop(&ring).is_some(),
            "a coupler and a waveguide sharing two ports is a cavity"
        );
        let link = parse_spice(STRAIGHT).expect("parse");
        assert_eq!(
            optical_loop(&link),
            None,
            "a laser, a waveguide and a detector in a line close nothing"
        );
    }

    /// The reason the graph is bipartite: a multi-port device is a star, not a
    /// cycle.
    ///
    /// Joining an element's own ports to each other instead would make any
    /// device with three or more optical ports a loop by itself, and every deck
    /// with a coupler in it would be reported as a resonator. This is the case
    /// that distinguishes the two constructions.
    #[test]
    fn one_four_port_coupler_is_not_a_loop() {
        let nl = parse_spice(
            "* a coupler, alone\n\
             .optical_port a\n\
             .optical_port b\n\
             .optical_port c\n\
             .optical_port d\n\
             Xlas a fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
             Xc a b c d fc_dcoupler kappa_L=0.3\n\
             Xp1 c p1 0 fc_photodetector responsivity=0.9\n\
             Xp2 d p2 0 fc_photodetector responsivity=0.9\n\
             R1 p1 0 1k\n\
             R2 p2 0 1k\n",
        )
        .expect("parse");
        assert_eq!(
            optical_loop(&nl),
            None,
            "four ports on one device is a star; only two devices sharing two \
             ports close a path"
        );
    }

    /// An MZI is two arms between the same pair of couplers — a loop
    /// topologically, and correctly so: light divides and recombines, which is
    /// what makes its transfer interferometric rather than a simple cascade.
    #[test]
    fn an_mzi_counts_as_a_closed_path() {
        let nl = parse_spice(
            "* MZI\n\
             .optical_port lin\n\
             .optical_port dark\n\
             .optical_port a1\n\
             .optical_port a2\n\
             .optical_port b1\n\
             .optical_port b2\n\
             .optical_port out\n\
             .optical_port unused\n\
             Xlas lin fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
             Xc1 lin dark a1 a2 fc_dcoupler kappa_L=0.785\n\
             Xw1 a1 b1 fc_waveguide L_um=100\n\
             Xw2 a2 b2 fc_waveguide L_um=110\n\
             Xc2 b1 b2 out unused fc_dcoupler kappa_L=0.785\n\
             Xpd out det 0 fc_photodetector responsivity=0.9\n\
             Rl det 0 1k\n",
        )
        .expect("parse");
        assert!(
            optical_loop(&nl).is_some(),
            "two arms between the same couplers is a closed path"
        );
    }

    #[test]
    fn a_deck_with_no_optical_ports_has_nothing_to_say() {
        let nl = parse_spice("* rc\nV1 a 0 DC 1\nR1 a b 1k\nC1 b 0 1p\n").expect("parse");
        assert_eq!(optical_loop(&nl), None);
    }
}
