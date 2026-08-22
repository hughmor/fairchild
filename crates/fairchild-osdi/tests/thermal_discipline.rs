//! A temperature as a solved unknown: the `thermal` discipline end to end.
//!
//! `tests/models/self_heated_r.va` is a resistor whose value follows its own
//! temperature and whose dissipation sets that temperature — a fixed point, so
//! nothing here can pass by evaluating a formula in the right order.
//!
//! Two claims are pinned, and they need different anchors:
//!
//!   * the fixed point is *found*, against the closed form below;
//!   * the row is known to be kelvin, so `temptol` bounds it and not `vntol`.
//!
//! The second cannot be inferred from the first — a converged answer is
//! converged whichever tolerance was used, since `vntol` on a temperature is
//! tighter than needed rather than looser. It is asserted directly on the
//! topology instead.

use std::path::Path;
use std::sync::Arc;

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_osdi::OsdiLibrary;
use fairchild_parser::parse_spice;

mod common;

/// Closed form for the self-heated resistor's operating point.
///
/// `ΔT = R_th·V²/R(ΔT)` with `R(ΔT) = R0·(1 + α·ΔT)` rearranges to
/// `α·ΔT² + ΔT − K = 0`, `K = R_th·V²/R0`, whose positive root is:
fn analytic_rise(v: f64, r0: f64, alpha: f64, r_th: f64) -> f64 {
    let k = r_th * v * v / r0;
    ((1.0 + 4.0 * alpha * k).sqrt() - 1.0) / (2.0 * alpha)
}

fn registry_with_model(path: &Path) -> DeviceRegistry {
    let lib = Arc::new(unsafe { OsdiLibrary::open(path) }.expect("dlopen"));
    let mut reg = DeviceRegistry::new();
    lib.register_into(&mut reg);
    reg
}

/// The thermal network lives in the deck, in ordinary SPICE primitives — which
/// is the point. `Rth` is 500 K/W because the flow into a thermal node is watts.
///
/// 30 V, not 10, so the rise is 233 K. That is large enough to matter to the
/// `vmax` trust region: a thermal row allowed to set a VOLT-scaled clamp asks
/// for a 233-unit step, scales every electrical unknown by 0.5/233, and Newton
/// then satisfies its step test on the clamped deltas rather than on the
/// residual — 10× the iterations *and* an answer 9e-6 K off. Restore that and
/// `self_heating_reaches_the_fixed_point` fails on the closed form below.
const DECK: &str = "\
* self-heated resistor with its thermal network written in the deck
V1 p 0 30
X1 p 0 h self_heated_r
Rth h 0 500
.op
.end
";

#[test]
fn self_heating_reaches_the_fixed_point() {
    let Some(path) = common::compiled("self_heated_r") else {
        return;
    };
    let reg = registry_with_model(&path);
    let net = parse_spice(DECK).expect("deck parses");
    let res = dc_op_nr_with_registry(&net, &reg).expect("DC OP converges");

    let want = analytic_rise(30.0, 1000.0, 4.0e-3, 500.0);
    let got = res.node_voltage("h").expect("h is a row");
    assert!(
        (got - want).abs() < 1e-6,
        "temperature rise {got} K, closed form says {want} K"
    );

    // The coupling, not just the number: at 42.7 K the resistance is 17 % above
    // its cold value, so a solver that ignored the thermal node would read
    // 10 mA here. Anchored on Ohm's law at the solved temperature rather than
    // on the model, so the two cannot agree by sharing a fault.
    let i = res
        .vsrc_current("v1")
        .expect("V1 branch current is a row")
        .abs();
    let r_hot = 1000.0 * (1.0 + 4.0e-3 * want);
    assert!(
        (i - 30.0 / r_hot).abs() < 1e-9,
        "current {i} A does not match Ohm's law at the solved temperature"
    );
    assert!(
        (i - 30.0 / 1000.0).abs() > 1e-4,
        "current {i} A is the cold-resistance answer — self-heating did nothing"
    );
}

/// `Temp` is a rise above ambient, not an absolute temperature.
///
/// The distinction is invisible in a single solve and produces a converged wrong
/// answer if it is got backwards, so the fixture separates the two: `alpha_a`
/// acts on ambient, `alpha_r` on the rise, and they differ by 4×. At 85 °C the
/// cold resistance moves by `alpha_a·58 K` and the rise is then computed from
/// *that* — a model reading 358.15 K into `Temp` would land nowhere near.
#[test]
fn rise_is_measured_from_ambient_not_absolute() {
    let Some(path) = common::compiled("self_heated_r") else {
        return;
    };
    let reg = registry_with_model(&path);
    let hot = DECK.replace(".op", ".options temp=85\n.op");
    let net = parse_spice(&hot).expect("deck parses");
    let res = dc_op_nr_with_registry(&net, &reg).expect("DC OP converges");

    // 85 °C = 358.15 K, so ambient is 58 K above the 27 °C the parameters
    // are quoted at.
    let r_cold = 1000.0 * (1.0 + 1.0e-3 * 58.0);
    let want = analytic_rise(30.0, r_cold, 4.0e-3, 500.0);
    let got = res.node_voltage("h").expect("h is a row");
    assert!(
        (got - want).abs() < 1e-6,
        "rise {got} K, closed form at 85 °C ambient says {want} K"
    );
    // And it is genuinely a different answer from the 27 °C one, so the
    // assertion above is not satisfied by ignoring `.options temp` entirely.
    let cold = analytic_rise(30.0, 1000.0, 4.0e-3, 500.0);
    assert!(
        (want - cold).abs() > 1e-3,
        "ambient made no difference — this test cannot fail"
    );
}

/// The row is bounded in kelvin, not volts.
///
/// Read off the topology because it is not observable in the answer: `vntol` is
/// *tighter* than a temperature needs, so a wrongly-classified row still
/// converges, to the same value, more slowly. Sabotage check: delete
/// `OsdiDevice::thermal_nodes` and this fails while every other test here passes.
#[test]
fn a_thermal_row_takes_temptol() {
    let Some(path) = common::compiled("self_heated_r") else {
        return;
    };
    let reg = registry_with_model(&path);
    let net = parse_spice(DECK).expect("deck parses");
    let ctx = fairchild_core::device::SimContext::default();
    let mut topo = fairchild_core::mna::CircuitTopology::build_resolved(&net, &ctx, &reg);
    let _devices =
        fairchild_core::newton::build_devices(&net, &mut topo, &ctx, &reg).expect("devices build");

    let opts = fairchild_core::options::SimOptions::default();
    let tol = fairchild_core::tolerance::Tolerances::build(&topo, &opts);

    let h = topo.node_index["h"];
    let p = topo.node_index["p"];
    assert_eq!(
        tol.bound(h, 0.0),
        opts.temptol,
        "the thermal node is bounded by vntol, so the model's `thermal h` \
         never reached the solver"
    );
    assert_eq!(
        tol.bound(p, 0.0),
        opts.vntol,
        "an electrical node picked up temptol — the classification is too broad"
    );
    assert_eq!(
        topo.thermal_rows,
        vec![h],
        "exactly one row carries kelvin in this deck"
    );
}

/// The same fixed point with the thermal node kept inside the model.
///
/// A different route to the same row — `num_extra_nodes` rather than a terminal
/// — and the one an electro-thermal model ported from a vendor library takes,
/// since those keep `r_th`/`c_th` as parameters. The answer must agree with the
/// ported-out version to the last digit, because it is the same physics; if it
/// does not, `push_device` is mapping the extras block wrong.
#[test]
fn an_internal_thermal_node_is_a_thermal_row_too() {
    let Some(path) = common::compiled("self_heated_r") else {
        return;
    };
    let reg = registry_with_model(&path);
    let deck = "\
* thermal network inside the model
V1 p 0 30
X1 p 0 self_heated_r_int
.op
.end
";
    let net = parse_spice(deck).expect("deck parses");
    let ctx = fairchild_core::device::SimContext::default();
    let mut topo = fairchild_core::mna::CircuitTopology::build_resolved(&net, &ctx, &reg);
    let _devices =
        fairchild_core::newton::build_devices(&net, &mut topo, &ctx, &reg).expect("devices build");

    // The internal node is not in `node_index` — it lives above the branch rows,
    // which is exactly why it needs `push_device` to place it.
    assert_eq!(
        topo.thermal_rows.len(),
        1,
        "the model's internal `thermal h` did not become a thermal row"
    );
    let row = topo.thermal_rows[0];
    assert!(
        row >= topo.node_index.len() + topo.vsrc_index.len(),
        "row {row} is not in the device-internal block — the extras offset is wrong"
    );

    let opts = fairchild_core::options::SimOptions::default();
    let tol = fairchild_core::tolerance::Tolerances::build(&topo, &opts);
    assert_eq!(tol.bound(row, 0.0), opts.temptol);

    let res = dc_op_nr_with_registry(&net, &reg).expect("DC OP converges");
    let i = res.vsrc_current("v1").expect("branch current").abs();
    let want = analytic_rise(30.0, 1000.0, 4.0e-3, 500.0);
    assert!(
        (i - 30.0 / (1000.0 * (1.0 + 4.0e-3 * want))).abs() < 1e-9,
        "internal-node answer {i} A disagrees with the ported-out one"
    );
}

/// Two devices heating each other through a plain `R`, which is what a thermal
/// discipline is *for*.
///
/// `X2` dissipates nothing — its own bias is zero — so any rise it shows came
/// through `Rc` from `X1`. There is no parameter that can express this: an
/// `r_th` inside a model is a path to ambient and nothing else, which is why
/// `mrm_wdm.va` disclaimed ring-to-ring crosstalk while its temperature was
/// internal.
///
/// The anchor is the thermal network solved by hand, which is legitimate here
/// because it is a *linear* resistor network with one injection: from `h1`,
/// `Rc + Rth2 = 700` sits in parallel with `Rth1 = 500`, so
/// `R_eff = 3500/12 K/W`, and `h2` divides that down by `500/700`.
#[test]
fn two_devices_couple_through_a_thermal_resistor() {
    let Some(path) = common::compiled("self_heated_r") else {
        return;
    };
    let reg = registry_with_model(&path);
    let deck = "\
* X1 dissipates, X2 does not, and they share a thermal network
V1 p 0 30
X1 p 0 h1 self_heated_r
Rth1 h1 0 500
Rc h1 h2 200
Vq q 0 0
X2 q 0 h2 self_heated_r
Rth2 h2 0 500
.op
.end
";
    let net = parse_spice(deck).expect("deck parses");
    let res = dc_op_nr_with_registry(&net, &reg).expect("DC OP converges");

    let t1 = analytic_rise(30.0, 1000.0, 4.0e-3, 3500.0 / 12.0);
    let t2 = t1 * 500.0 / 700.0;
    let got1 = res.node_voltage("h1").expect("h1 is a row");
    let got2 = res.node_voltage("h2").expect("h2 is a row");
    assert!(
        (got1 - t1).abs() < 1e-6,
        "hot device: {got1} K, network says {t1} K"
    );
    assert!(
        (got2 - t2).abs() < 1e-6,
        "cold device: {got2} K, network says {t2} K — heat did not cross `Rc`"
    );
    assert!(
        got2 > 1.0,
        "the unbiased device is at {got2} K, so nothing coupled at all"
    );
}
