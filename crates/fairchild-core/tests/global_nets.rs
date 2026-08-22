//! `.global` — a net that is the same node in every subcircuit scope.
//!
//! The parser-side tests check the names; this checks the *answer*, which is the
//! only thing that proves the flattened nodes actually joined up. A deck whose
//! supplies reach two nested instances without appearing in any port list has to
//! produce the same voltages as the equivalent flat circuit, and if the global
//! nets had been namespaced per instance the dividers would float instead.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

const TOL: f64 = 1e-7;

fn op(deck: &str) -> fairchild_core::NrResult {
    let nl = parse_spice(deck).unwrap_or_else(|e| panic!("parse failed: {e}"));
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&nl.models);
    dc_op_nr_with_registry(&nl, &reg).expect("DC OP failed")
}

/// Two `inv` instances, nested one level inside `chain`, each a 1 k/1 k divider
/// between `vdd` and `vss`. Neither supply is a port of either subcircuit.
///
/// Both midpoints must sit at 1.8/2 = 0.9 V, and the supply must carry both
/// dividers: 2 × 1.8 V / 2 kΩ = 1.8 mA.
#[test]
fn global_supplies_reach_two_levels_of_nesting() {
    let r = op("\
* CDL-style supplies
.global vdd vss
.subckt inv a y
Rpull vdd y 1k
Rdown y vss 1k
.ends
.subckt chain i o
X1 i m inv
X2 m o inv
.ends
Vsup vdd 0 1.8
Vgnd vss 0 0
Xtop in out chain
.op
");
    let mid = r.node_voltage("xtop.m").expect("no node `xtop.m`");
    let out = r.node_voltage("out").expect("no node `out`");
    assert!((mid - 0.9).abs() < TOL, "inner divider: {mid}");
    assert!((out - 0.9).abs() < TOL, "outer divider: {out}");
    let i_sup = r.vsrc_current("vsup").expect("no source `vsup`");
    assert!(
        (i_sup + 1.8e-3).abs() < 1e-9,
        "the supply must carry both dividers, got {i_sup}"
    );
}

/// The anchor: the same deck without the declaration. `vdd` and `vss` become
/// per-instance nodes that nothing drives, and the run fails naming them. Without
/// this, the test above would pass on a build that ignored `.global` entirely and
/// happened to leave the names unqualified.
#[test]
fn the_same_deck_without_global_does_not_connect() {
    let deck = "\
* no .global
.subckt inv a y
Rpull vdd y 1k
Rdown y vss 1k
.ends
.subckt chain i o
X1 i m inv
X2 m o inv
.ends
Vsup vdd 0 1.8
Vgnd vss 0 0
Xtop in out chain
.op
";
    let nl = parse_spice(deck).unwrap();
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&nl.models);
    let msg = match dc_op_nr_with_registry(&nl, &reg) {
        Err(e) => format!("{e:?}"),
        Ok(_) => panic!("per-instance supplies are driven by nothing; this must fail"),
    };
    assert!(
        msg.contains("xtop.x1.vdd") && msg.contains("xtop.x2.vdd"),
        "the failure should name the namespaced supplies: {msg}"
    );
}
