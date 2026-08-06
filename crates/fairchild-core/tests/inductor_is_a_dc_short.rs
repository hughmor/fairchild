//! An ideal inductor is a SHORT at DC, not an open.
//!
//! The DC assembler skipped any inductor with no companion state, with the
//! comment "inductor = open circuit at DC". That is the capacitor's rule, not
//! the inductor's: the dual of an open capacitor is a shorted inductor. The
//! effect was silent and severe.
//!
//! * A source feeding a load through a choke or a bondwire read 0 V.
//! * A *current* source into one left its node with nothing but gmin on the
//!   diagonal, so the first Newton step demanded I/gmin volts — 1e9 V for 1 mA.
//!   That collapsed the vmax trust region, which scaled every other unknown
//!   (including λ, in metres) into oblivion, and nothing converged.
//!
//! On the giona chip — 105 bondwire inductors — the operating point was
//! unreachable. With the short in place it solves in about a second, and the
//! heater node sits at 1 mA × 200 Ω through its two series heaters, which is
//! what a hand trace of the bias network predicts.

use fairchild_core::dc_op_nr_with_registry;
use fairchild_core::device_registry::DeviceRegistry;
use fairchild_parser::parse_spice;

fn v(deck: &str, node: &str) -> f64 {
    let net = parse_spice(deck).expect("parse");
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    r.node_voltage(node).expect("node")
}

#[test]
fn a_voltage_source_drives_through_an_inductor() {
    // The choke passes DC, so the whole 1 V lands on the load.
    let got = v(
        "* choke\nV1 a 0 DC 1\nL1 a b 1n\nR1 b 0 1k\n.op\n.end\n",
        "b",
    );
    assert!(
        (got - 1.0).abs() < 1e-6,
        "expected 1 V across the load through a DC short, got {got}"
    );
}

#[test]
fn a_current_source_drives_through_an_inductor() {
    // This is the case that used to leave the node on gmin alone: 1 mA into
    // 1e-12 S is 1e9 V, which then wrecked the trust region for every unknown.
    let got = v(
        "* bondwire\nI1 0 a DC 1m\nL1 a b 1n\nR1 b 0 1k\n.op\n.end\n",
        "a",
    );
    assert!(
        (got - 1.0).abs() < 1e-4,
        "expected 1 mA × 1 kΩ = 1 V, got {got}"
    );
}

#[test]
fn inductors_in_series_all_pass_dc() {
    let got = v(
        "* two chokes\nV1 a 0 DC 2\nL1 a b 1n\nL2 b c 2n\nR1 c 0 1k\n.op\n.end\n",
        "c",
    );
    assert!((got - 2.0).abs() < 1e-6, "expected 2 V, got {got}");
}

#[test]
fn an_inductor_to_ground_pulls_its_node_to_ground() {
    // The dual check: a short to ground must actually short, so the divider
    // collapses instead of holding half the supply.
    let got = v(
        "* L to ground\nV1 a 0 DC 1\nR1 a b 1k\nL1 b 0 1n\n.op\n.end\n",
        "b",
    );
    assert!(
        got.abs() < 1e-6,
        "an inductor to ground must pull its node to ~0 V, got {got}"
    );
}

#[test]
fn the_heater_chain_that_broke_giona_now_biases_correctly() {
    // The giona_fc bias path, reduced: a current source reaches ground through
    // two bondwires and two 100 Ω heater resistances, so the driven node sits at
    // 1 mA × 200 Ω. Before the fix the bondwires were open and this node wanted
    // 1e9 V.
    let got = v(
        "* bondwire - heater - heater - bondwire\n\
         I1 0 w11 DC 1m\n\
         L1 w11 h11 1n\n\
         Rh1 h11 mid 100\n\
         Rh2 mid g1 100\n\
         L2 g1 0 1n\n.op\n.end\n",
        "w11",
    );
    assert!(
        (got - 0.2).abs() < 1e-4,
        "expected 1 mA × 200 Ω = 0.2 V at the heater node, got {got}"
    );
}
