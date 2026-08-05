//! `E`, `F`, `G`, `H` — the four linear controlled sources, against ngspice 46.
//!
//! They are desugared onto the B-element rather than given their own stamps (a
//! VCVS *is* `B… V=gain*(V(cp)-V(cn))`), so what needs pinning is not the
//! algebra but the **conventions**: which pair is the output, which the control,
//! and which direction a current-output source pushes. Every expected value
//! below was read off ngspice on the same deck — getting a sign backwards here
//! would otherwise look perfectly plausible.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};

use fairchild_parser::parse_spice;
/// Absolute tolerance for a node voltage. Loose enough to absorb `gmin`, which
/// fairchild adds to every node diagonal: on a 1 k divider that is worth a few
/// times 1e-10, so a 1e-9 bound would fail on arithmetic rather than on physics.
/// Still four orders tighter than ngspice's printed precision.
const TOL: f64 = 1e-7;

fn v_out(deck: &str) -> f64 {
    let nl = parse_spice(deck).unwrap_or_else(|e| panic!("parse failed: {e}"));
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&nl.models);
    let r = dc_op_nr_with_registry(&nl, &reg).expect("DC OP failed");
    r.node_voltage("out").expect("no node `out`")
}

/// `E<n> p n nc+ nc- gain` — V(out) = gain·V(in). ngspice: 3.0 V.
#[test]
fn vcvs_gain_matches_ngspice() {
    let v = v_out("* E\nVin in 0 DC 1.5\nE1 out 0 in 0 2.0\nRl out 0 1k\n.op\n.end\n");
    assert!((v - 3.0).abs() < TOL, "V(out) = {v}, ngspice gives 3.0");
}

/// `G<n> p n nc+ nc- gm` — a transconductance. The current leaves `n+` and
/// enters `n-`, the same convention the resistor stamp uses (a VCCS whose
/// control nodes *are* its output nodes, with gm = 1/R, must be that resistor),
/// so driving into a grounded 1 k from a +1.5 V control gives a **negative**
/// output. ngspice: −1.5 V.
#[test]
fn vccs_sign_matches_ngspice() {
    let v = v_out("* G\nVin in 0 DC 1.5\nG1 out 0 in 0 1m\nRl out 0 1k\n.op\n.end\n");
    assert!((v + 1.5).abs() < TOL, "V(out) = {v}, ngspice gives -1.5");
}

/// `F<n> p n Vctrl gain` — current gain off another source's branch current.
/// I(Vs) = 1 mA through the 1 k, gain 3, into 1 k. ngspice: 3.0 V.
#[test]
fn cccs_gain_matches_ngspice() {
    let v = v_out("* F\nVs a 0 DC 1\nRs a 0 1k\nF1 out 0 Vs 3.0\nRl out 0 1k\n.op\n.end\n");
    assert!((v - 3.0).abs() < TOL, "V(out) = {v}, ngspice gives 3.0");
}

/// `H<n> p n Vctrl gain` — a transresistance. ngspice: −0.5 V.
#[test]
fn ccvs_gain_matches_ngspice() {
    let v = v_out("* H\nVs a 0 DC 1\nRs a 0 1k\nH1 out 0 Vs 500\nRl out 0 1k\n.op\n.end\n");
    assert!((v + 0.5).abs() < TOL, "V(out) = {v}, ngspice gives -0.5");
}

/// A VCCS whose control nodes are its own output nodes is a resistor. This is
/// the statement that fixes the sign convention rather than merely asserting it:
/// if the stamp were transposed or negated, `gm = 1/R` would not reproduce `R`.
#[test]
fn a_vccs_across_itself_is_a_resistor() {
    let via_g = v_out("* G as R\nVin in 0 DC 1\nR1 in out 1k\nG1 out 0 out 0 1m\n.op\n.end\n");
    let via_r = v_out("* real R\nVin in 0 DC 1\nR1 in out 1k\nR2 out 0 1k\n.op\n.end\n");
    assert!(
        (via_g - via_r).abs() < 1e-12,
        "G1 with gm=1/R gives {via_g}, a real 1k gives {via_r}"
    );
    assert!(
        (via_r - 0.5).abs() < TOL,
        "the divider itself is wrong: {via_r}"
    );
}

/// Controlled sources have to work in AC too — they are stamped every
/// iteration, not folded into the DC operating point.
#[test]
fn vcvs_gain_carries_into_ac() {
    use fairchild_core::ac::ac_analysis;
    let run = |gain: &str| {
        let nl = parse_spice(&format!(
            "* E in ac\nVin in 0 AC 1\nE1 mid 0 in 0 {gain}\n\
             R1 mid out 1k\nC1 out 0 1n\n.end\n"
        ))
        .unwrap();
        let r = ac_analysis(&nl, &[1e3], None, &DeviceRegistry::new()).expect("ac failed");
        let (re, im) = r.voltages.get("out").unwrap()[0];
        (re * re + im * im).sqrt()
    };
    let (one, two) = (run("1.0"), run("2.0"));
    assert!(
        (two - 2.0 * one).abs() < TOL,
        "doubling the VCVS gain should double the AC response: {one} -> {two}"
    );
}
