//! A Spectre deck shaped like a foundry wrapper, solved for its answer.
//!
//! The parser tests check that each construct transliterates. This checks the
//! *number*, which is the only thing that proves the four layers actually
//! compose: `parameters` hoisted onto the `.subckt` header, resolved per
//! instance, a body conditional selecting from that instance's values, and an
//! `m=` multiplier scaling what the body flattened to. Every expected value here
//! is worked out by hand from Ohm's law, so nothing agrees with itself.
//!
//! The fixtures are hand-written to mimic the *shape* of a foundry wrapper. No
//! foundry text appears in this repository.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spectre;

const TOL: f64 = 1e-9;

fn op(deck: &str) -> fairchild_core::NrResult {
    let nl = parse_spectre(deck).unwrap_or_else(|e| panic!("parse failed: {e}"));
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&nl.models);
    dc_op_nr_with_registry(&nl, &reg).expect("DC OP failed")
}

/// A wrapper whose sheet resistance is a `parameters` expression over two others,
/// instantiated twice with different geometry — the shape of every parasitic
/// wrapper in a PDK.
///
/// `rval = sheet*l/w`, so the two instances are 100 Ω (l=2u, w=2u) and 400 Ω
/// (l=4u, w=1u) in series across 1 V: 1/500 = 2 mA, and the midpoint sits at
/// 1 − 100/500 = 0.8 V.
#[test]
fn a_wrapper_resolves_its_geometry_per_instance() {
    let r = op("\
// two instances of one wrapper
simulator lang=spectre
inline subckt res_wrap (a b)
parameters w=1u l=1u sheet=100
parameters rval='sheet*l/w'
R1 (a b) resistor r=rval
ends res_wrap
V1 (in 0) vsource dc=1
X1 (in mid) res_wrap w=2u l=2u
X2 (mid 0) res_wrap w=1u l=4u
dc1 dc
");
    let mid = r.node_voltage("mid").expect("no node `mid`");
    assert!(
        (mid - 0.8).abs() < TOL,
        "midpoint {mid} V: the two instances did not resolve their own geometry"
    );
    let i = r.vsrc_current("v1").expect("no source `v1`");
    assert!((i.abs() - 2e-3).abs() < TOL, "supply current {i} A");
}

/// A switch inside the wrapper, on a parameter the caller sets — a corner or a
/// self-heating flag, structurally. Instance `a` takes the branch, `b` does not,
/// so the two arms of the divider differ: 2 kΩ over 2 kΩ + 1 kΩ.
#[test]
fn a_wrapper_switch_selects_per_instance() {
    let r = op("\
// one definition, two switch settings
simulator lang=spectre
subckt rsel (a b)
parameters mode=0
if (mode == 1) {
R1 (a b) resistor r=2k
} else {
R1 (a b) resistor r=1k
}
ends rsel
V1 (in 0) vsource dc=3
Xa (in mid) rsel mode=1
Xb (mid 0) rsel
dc1 dc
");
    let mid = r.node_voltage("mid").expect("no node `mid`");
    // 3 V across 3 kΩ; the 1 kΩ arm is the lower one.
    assert!(
        (mid - 1.0).abs() < TOL,
        "midpoint {mid} V: both instances took the same branch"
    );
}

/// `m=` on the instance, over a body the wrapper computed: four copies of a
/// 4 kΩ wrapper in parallel is 1 kΩ, so 1 V draws 1 mA.
#[test]
fn a_multiplier_scales_a_wrapper() {
    let r = op("\
// four in parallel
simulator lang=spectre
inline subckt res_wrap (a b)
parameters sheet=1000 n=4
parameters rval='sheet*n'
R1 (a b) resistor r=rval
ends res_wrap
V1 (in 0) vsource dc=1
X1 (in 0) res_wrap m=4
dc1 dc
");
    let i = r.vsrc_current("v1").expect("no source `v1`");
    assert!(
        (i.abs() - 1e-3).abs() < TOL,
        "supply current {i} A: want 1 mA (4 kΩ / 4)"
    );
}

/// A function definition and a global net, the two remaining pieces of the
/// dialect a wrapper leans on: `vdd!` is the same node everywhere without being
/// a port, and `half()` is a `.func`. 1.8 V over 2 × 900 Ω = 1 mA.
#[test]
fn a_function_and_a_global_net_reach_the_solve() {
    let r = op("\
// a supply that is not a port, and a function
simulator lang=spectre
real half(real x) { return x/2; }
subckt leg (o)
parameters rtot=1800
R1 (vdd! o) resistor r='half(rtot)'
R2 (o 0) resistor r='half(rtot)'
ends leg
V1 (vdd! 0) vsource dc=1.8
X1 (out) leg
dc1 dc
");
    let out = r
        .node_voltage("x1.out")
        .or_else(|_| r.node_voltage("out"))
        .expect("no node `out`");
    assert!((out - 0.9).abs() < TOL, "midpoint {out} V");
    let i = r.vsrc_current("v1").expect("no source `v1`");
    assert!((i.abs() - 1e-3).abs() < TOL, "supply current {i} A");
}
