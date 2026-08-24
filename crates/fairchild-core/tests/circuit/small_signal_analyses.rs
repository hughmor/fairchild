//! `.tf`, `.sens` and `.pz` against closed forms.
//!
//! Every circuit here is one whose answer can be written down: a divider whose
//! gain and both resistances are algebra, an RC whose pole is `−1/RC`, an RLC
//! whose pair is `−ζω₀ ± jω₀√(1−ζ²)`.  That is deliberate — these three
//! analyses have sign and normalisation conventions that are easy to get
//! *self-consistently* wrong, and a test that only checks two of our own
//! subsystems against each other cannot see a fault they share.  The ngspice
//! goldens in `tests/ngspice/` are the second, independent anchor.

use fairchild_core::pz::pole_zero;
use fairchild_core::sens::sensitivity;
use fairchild_core::tf::transfer_function;
use fairchild_core::{DeviceRegistry, SimOptions};
use fairchild_parser::{parse_spice, Analysis, Netlist, PzDrive, PzWant};

fn registry(net: &Netlist) -> DeviceRegistry {
    let mut r = DeviceRegistry::new();
    r.register_builtin_models(&net.models);
    r
}

/// Tight enough that a closed form is the thing under test, not `reltol`.
fn tight() -> SimOptions {
    SimOptions {
        reltol: 1e-12,
        vntol: 1e-14,
        ..SimOptions::default()
    }
}

/// The one `.tf` / `.sens` / `.pz` card a deck declares, as parsed.
fn sole(src: &str) -> (Netlist, Analysis) {
    let net = parse_spice(src).unwrap();
    let a = net
        .analyses
        .first()
        .expect("deck declares no analysis")
        .clone();
    (net, a)
}

/// The adjoint's finite difference on the residual is good to ~1e-9 relative
/// (see `crate::adjoint`), and a resistance is a reciprocal of one, so `1e-7`
/// is the honest band here — tighter would be testing the FD step, not the
/// analysis.
const REL: f64 = 1e-7;

fn assert_close(got: f64, want: f64, rel: f64, what: &str) {
    let tol = rel * want.abs().max(1e-12);
    assert!(
        (got - want).abs() <= tol,
        "{what}: got {got:e}, want {want:e} (tol {tol:e})"
    );
}

// ---------------------------------------------------------------------------
// .tf
// ---------------------------------------------------------------------------

/// Resistive divider: every one of the three numbers is algebra.
///
/// `gain = R2/(R1+R2)`, `r_in = R1+R2` (what the source drives), and
/// `r_out = R1∥R2` (what the output port sees looking back, with the source
/// replaced by its small-signal short).
#[test]
fn tf_divider_matches_algebra() {
    let src = "* divider\n\
               Vin in 0 DC 1\n\
               R1 in out 1k\n\
               R2 out 0 3k\n\
               .tf v(out) Vin\n\
               .end\n";
    let (net, a) = sole(src);
    let Analysis::Tf { out, input_src } = &a else {
        panic!("expected .tf, got {a:?}")
    };
    let r = transfer_function(&net, &registry(&net), &tight(), out, input_src).unwrap();

    assert_close(r.gain, 3000.0 / 4000.0, REL, "gain");
    assert_close(r.r_in, 4000.0, REL, "r_in");
    assert_close(r.r_out, 1000.0 * 3000.0 / 4000.0, REL, "r_out");
    assert_close(r.out_value, 0.75, REL, "out_value");
}

/// The output port is not always ground-referenced.  `v(a,b)` must measure
/// across the pair, and the output resistance must be the one *that port* sees
/// — here R2 alone, since the probe's current returns through R2's own nodes.
#[test]
fn tf_differential_output_port() {
    let src = "* divider\n\
               Vin in 0 DC 1\n\
               R1 in mid 1k\n\
               R2 mid out 2k\n\
               R3 out 0 1k\n\
               .tf v(mid,out) Vin\n\
               .end\n";
    let (net, a) = sole(src);
    let Analysis::Tf { out, input_src } = &a else {
        panic!()
    };
    let r = transfer_function(&net, &registry(&net), &tight(), out, input_src).unwrap();

    // V(mid,out) is R2's share of the 4k chain.
    assert_close(r.gain, 2000.0 / 4000.0, REL, "gain");
    assert_close(r.r_in, 4000.0, REL, "r_in");
    // Looking into (mid,out): R2 in parallel with (R1 + R3), Vin shorted.
    assert_close(r.r_out, 1.0 / (1.0 / 2000.0 + 1.0 / 2000.0), REL, "r_out");
}

/// A current-source input: the transfer is a transresistance and `r_in` is the
/// resistance across the source's own terminals, not a reciprocal.
#[test]
fn tf_current_input_is_a_transresistance() {
    let src = "* shunt\n\
               Iin 0 out DC 1m\n\
               R1 out 0 2k\n\
               R2 out 0 2k\n\
               .tf v(out) Iin\n\
               .end\n";
    let (net, a) = sole(src);
    let Analysis::Tf { out, input_src } = &a else {
        panic!()
    };
    let r = transfer_function(&net, &registry(&net), &tight(), out, input_src).unwrap();

    // `I 0 out` drives current *into* `out`, so V(out) = +I·(R1∥R2).
    assert_close(r.gain, 1000.0, REL, "transresistance");
    assert_close(r.r_in, 1000.0, REL, "r_in");
    assert_close(r.r_out, 1000.0, REL, "r_out");
}

/// A current output through a sense source.  `.tf i(vsense) Vin` on a series
/// loop: the transfer is `1/R_total` and both resistances are `R_total`.
///
/// The sign is the interesting part, and it is *not* the same as the driving
/// source's.  SPICE counts a source's branch current as flowing into its `+`
/// terminal: Vin pushes current out of its own `+`, so `I(Vin)` is negative,
/// while the loop current arrives at Vsense's `+` from R1, so `I(Vsense)` is
/// positive.  Getting these two the same way round is the classic `.tf` sign
/// error.
#[test]
fn tf_current_output() {
    let src = "* loop\n\
               Vin in 0 DC 2\n\
               R1 in mid 1k\n\
               Vsense mid 0 DC 0\n\
               .tf i(Vsense) Vin\n\
               .end\n";
    let (net, a) = sole(src);
    let Analysis::Tf { out, input_src } = &a else {
        panic!()
    };
    let r = transfer_function(&net, &registry(&net), &tight(), out, input_src).unwrap();

    assert_close(r.gain, 1.0 / 1000.0, REL, "transconductance");
    assert_close(r.r_in, 1000.0, REL, "r_in");
    assert_close(r.r_out, 1000.0, REL, "r_out");
}

/// A `.tf` whose named input is not a source is a mistake worth a sentence,
/// not a zero.
#[test]
fn tf_rejects_a_non_source_input() {
    let src = "* divider\n\
               Vin in 0 DC 1\n\
               R1 in out 1k\n\
               R2 out 0 1k\n\
               .tf v(out) R1\n\
               .end\n";
    let (net, a) = sole(src);
    let Analysis::Tf { out, input_src } = &a else {
        panic!()
    };
    let err = transfer_function(&net, &registry(&net), &tight(), out, input_src).unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("not an independent source"), "{msg}");
}

// ---------------------------------------------------------------------------
// .sens
// ---------------------------------------------------------------------------

/// Divider sensitivities have closed forms:
/// `∂V/∂R1 = −Vin·R2/(R1+R2)²`, `∂V/∂R2 = +Vin·R1/(R1+R2)²`,
/// `∂V/∂Vin = R2/(R1+R2)`.
#[test]
fn sens_divider_matches_algebra() {
    let src = "* divider\n\
               Vin in 0 DC 2\n\
               R1 in out 1k\n\
               R2 out 0 3k\n\
               .sens v(out)\n\
               .end\n";
    let (net, a) = sole(src);
    let Analysis::Sens { out, params } = &a else {
        panic!("expected .sens, got {a:?}")
    };
    let r = sensitivity(&net, &registry(&net), &tight(), out, params).unwrap();

    // Bare `.sens` takes every element value, in netlist order.  Names come
    // back as the parser canonicalises them — lowercase, like every other
    // element reference in a result.
    let names: Vec<&str> = r.rows.iter().map(|row| row.name.as_str()).collect();
    assert_eq!(names, vec!["vin.value", "r1.value", "r2.value"]);
    assert!(r.rows.iter().all(|row| row.reached), "{:?}", r.rows);

    let (vin, r1, r2) = (2.0, 1000.0, 3000.0);
    let denom = (r1 + r2) * (r1 + r2);
    assert_close(r.rows[0].sensitivity, r2 / (r1 + r2), 1e-7, "dV/dVin");
    assert_close(r.rows[1].sensitivity, -vin * r2 / denom, 1e-7, "dV/dR1");
    assert_close(r.rows[2].sensitivity, vin * r1 / denom, 1e-7, "dV/dR2");

    // Normalised is the per-100 % figure, which is what makes a 1 kΩ and a 2 V
    // parameter comparable in one table.
    assert_close(
        r.rows[1].normalised,
        -vin * r2 * r1 / denom,
        1e-7,
        "normalised dV/dR1",
    );
}

/// Naming parameters explicitly keeps the card's order, and a parameter the
/// adjoint cannot reach is reported as unreached — never as a zero that reads
/// like a real insensitivity.
#[test]
fn sens_named_params_and_honest_unreached() {
    let src = "* divider\n\
               Vin in 0 DC 2\n\
               R1 in out 1k\n\
               R2 out 0 3k\n\
               .sens v(out) R2 R1\n\
               .end\n";
    let (net, a) = sole(src);
    let Analysis::Sens { out, params } = &a else {
        panic!()
    };
    let r = sensitivity(&net, &registry(&net), &tight(), out, params).unwrap();
    let names: Vec<&str> = r.rows.iter().map(|row| row.name.as_str()).collect();
    assert_eq!(names, vec!["r2.value", "r1.value"]);
    assert!(r.rows[0].sensitivity > 0.0, "dV/dR2 should raise the tap");
    assert!(r.rows[1].sensitivity < 0.0, "dV/dR1 should lower the tap");
    assert!(r.unreached().is_empty());
}

// ---------------------------------------------------------------------------
// .pz
// ---------------------------------------------------------------------------

/// First-order RC: one pole at `−1/RC`, no zeros, and the algebraic modes
/// reported as infinite rather than dropped.
#[test]
fn pz_rc_pole_is_minus_one_over_rc() {
    let src = "* rc\n\
               Vin in 0 DC 0 AC 1\n\
               R1 in out 1k\n\
               C1 out 0 1n\n\
               .pz in 0 out 0 vol pz\n\
               .end\n";
    let (net, a) = sole(src);
    let Analysis::Pz {
        in_pos,
        in_neg,
        out_pos,
        out_neg,
        drive,
        want,
    } = &a
    else {
        panic!("expected .pz, got {a:?}")
    };
    assert_eq!(*drive, PzDrive::Vol);
    assert_eq!(*want, PzWant::Both);

    let r = pole_zero(
        &net,
        &registry(&net),
        &tight(),
        in_pos,
        in_neg,
        out_pos,
        out_neg,
        *drive,
        *want,
    )
    .unwrap();

    assert_eq!(r.poles.len(), 1, "poles: {:?}", r.poles);
    assert_close(r.poles[0].re, -1.0 / (1000.0 * 1e-9), 1e-6, "pole");
    assert!(r.poles[0].im.abs() < 1.0, "pole should be real");
    assert!(r.zeros.is_empty(), "zeros: {:?}", r.zeros);
    assert!(r.infinite_poles > 0, "the vsrc row is algebraic");
}

/// Series RLC, output across C: a complex pair at `−R/2L ± j√(1/LC − (R/2L)²)`.
///
/// This is the case the module's whole branch-current construction exists for
/// — with the inductor left as its `1/(jωL)` admittance the pencil would be
/// quadratic, and the extra roots would have to be told apart from these two.
#[test]
fn pz_series_rlc_complex_pair() {
    let (r, l, c) = (100.0_f64, 1e-3_f64, 1e-9_f64);
    let src = "* rlc\n\
               Vin in 0 DC 0 AC 1\n\
               R1 in a 100\n\
               L1 a out 1m\n\
               C1 out 0 1n\n\
               .pz in 0 out 0 vol pol\n\
               .end\n";
    let (net, a) = sole(src);
    let Analysis::Pz {
        in_pos,
        in_neg,
        out_pos,
        out_neg,
        drive,
        want,
    } = &a
    else {
        panic!()
    };
    let res = pole_zero(
        &net,
        &registry(&net),
        &tight(),
        in_pos,
        in_neg,
        out_pos,
        out_neg,
        *drive,
        *want,
    )
    .unwrap();

    let sigma = -r / (2.0 * l);
    let omega = (1.0 / (l * c) - sigma * sigma).sqrt();
    assert_eq!(res.poles.len(), 2, "poles: {:?}", res.poles);
    for p in &res.poles {
        assert_close(p.re, sigma, 1e-6, "pole real part");
        assert_close(p.im.abs(), omega, 1e-6, "pole imaginary part");
    }
    // Reported as a conjugate pair, both halves, not folded.
    assert!(
        res.poles[0].im * res.poles[1].im < 0.0,
        "expected a conjugate pair: {:?}",
        res.poles
    );
    assert!(res.zeros.is_empty(), "pol asked for poles only");
}

/// A zero, put somewhere it can be checked: a series RC to ground gives a
/// transfer with a zero at `−1/(R2·C)` and a pole further out.
#[test]
fn pz_finds_a_zero() {
    let (r1, r2, c) = (1000.0, 250.0, 1e-9);
    let src = "* zero\n\
               Vin in 0 DC 0 AC 1\n\
               R1 in out 1k\n\
               R2 out mid 250\n\
               C1 mid 0 1n\n\
               .pz in 0 out 0 vol pz\n\
               .end\n";
    let (net, a) = sole(src);
    let Analysis::Pz {
        in_pos,
        in_neg,
        out_pos,
        out_neg,
        drive,
        want,
    } = &a
    else {
        panic!()
    };
    let res = pole_zero(
        &net,
        &registry(&net),
        &tight(),
        in_pos,
        in_neg,
        out_pos,
        out_neg,
        *drive,
        *want,
    )
    .unwrap();

    // H(s) = (R2 + 1/sC) / (R1 + R2 + 1/sC): zero at −1/(R2·C), pole at
    // −1/((R1+R2)·C).
    assert_eq!(res.zeros.len(), 1, "zeros: {:?}", res.zeros);
    assert_close(res.zeros[0].re, -1.0 / (r2 * c), 1e-6, "zero");
    assert_eq!(res.poles.len(), 1, "poles: {:?}", res.poles);
    assert_close(res.poles[0].re, -1.0 / ((r1 + r2) * c), 1e-6, "pole");
}

/// A `cur` drive leaves the input port open where `vol` shorts it, and the two
/// pole sets differ because of it.  Reporting the same numbers for both would
/// mean the drive keyword was being ignored.
#[test]
fn pz_cur_and_vol_see_different_networks() {
    // The input port carries no voltage source, so the drive keyword is what
    // decides whether it is shorted or open.  (With a `Vin` sitting across it,
    // both answers would be the shorted one — the deck's own source pins the
    // node whatever the card says, which is why this deck excites with an `I`.)
    let src = "* twopole\n\
               Iin 0 in DC 0\n\
               R1 in out 1k\n\
               C1 out 0 1n\n\
               R2 in 0 1k\n\
               C2 in 0 2n\n\
               .end\n";
    let net = parse_spice(src).unwrap();
    let reg = registry(&net);
    let go = |d| {
        pole_zero(
            &net,
            &reg,
            &tight(),
            "in",
            "0",
            "out",
            "0",
            d,
            PzWant::Poles,
        )
        .unwrap()
    };

    // `vol` puts a source across (in, 0): node `in` is pinned, so R2 and C2 are
    // shorted out and only the R1–C1 pole survives.
    let vol = go(PzDrive::Vol);
    assert_eq!(vol.poles.len(), 1, "{:?}", vol.poles);
    assert_close(vol.poles[0].re, -1.0 / (1000.0 * 1e-9), 1e-6, "vol pole");

    // `cur` injects into an open port, so both capacitors are in play.
    let cur = go(PzDrive::Cur);
    assert_eq!(cur.poles.len(), 2, "{:?}", cur.poles);
}

/// The size ceiling is a refusal with a number in it, not an hour of silence.
#[test]
fn pz_refuses_past_the_dense_limit() {
    let mut src = String::from("* big\nVin in 0 DC 0 AC 1\n");
    // A long RC ladder: one node and one capacitor per rung.
    let rungs = fairchild_core::MAX_PZ_SIZE + 10;
    src.push_str("R0 in n1 1k\n");
    for i in 1..rungs {
        src.push_str(&format!("R{i} n{i} n{} 1k\n", i + 1));
        src.push_str(&format!("C{i} n{i} 0 1n\n"));
    }
    src.push_str(&format!("Cend n{rungs} 0 1n\n.end\n"));
    let net = parse_spice(&src).unwrap();
    let reg = registry(&net);
    let err = pole_zero(
        &net,
        &reg,
        &SimOptions::default(),
        "in",
        "0",
        "n2",
        "0",
        PzDrive::Vol,
        PzWant::Poles,
    )
    .unwrap_err();
    let msg = err.to_string();
    assert!(
        msg.contains(&fairchild_core::MAX_PZ_SIZE.to_string()),
        "{msg}"
    );
    assert!(msg.contains("dense"), "{msg}");
}

/// `.pz`'s poles must be the same numbers `.ac` bends at.  An independent
/// anchor on the eigensolve: the magnitude response of the RC at its computed
/// pole frequency is 1/√2 of DC.
#[test]
fn pz_pole_agrees_with_the_ac_sweep() {
    let src = "* rc\n\
               Vin in 0 DC 0 AC 1\n\
               R1 in out 3.3k\n\
               C1 out 0 4.7n\n\
               .end\n";
    let net = parse_spice(src).unwrap();
    let reg = registry(&net);
    let pz = pole_zero(
        &net,
        &reg,
        &tight(),
        "in",
        "0",
        "out",
        "0",
        PzDrive::Vol,
        PzWant::Poles,
    )
    .unwrap();
    assert_eq!(pz.poles.len(), 1);
    let f_pole = pz.poles[0].re.abs() / std::f64::consts::TAU;

    let ac = fairchild_core::ac_analysis_opts(
        &net,
        &[f_pole],
        Some("vin"),
        &reg,
        &SimOptions::default(),
    )
    .unwrap();
    let mag = ac.magnitude("out", 0).unwrap();
    assert_close(mag, 1.0 / 2.0_f64.sqrt(), 1e-6, "|H| at the pole frequency");
}
