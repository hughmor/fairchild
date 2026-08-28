//! A parameter is honoured, or it is named. Nothing in between.
//!
//! Three separate failures used to live behind one symptom — the deck runs and
//! you get a number:
//!
//! * an instance parameter that reached the netlist and stopped (`area=2` on a
//!   diode or a BJT: parsed, carried, dropped, silent);
//! * a model-card parameter the model matched explicitly and discarded (`IKF`,
//!   `MJSW`: no diagnostic of any kind);
//! * the audit in `docs/model_status.md` drifting away from either.
//!
//! The positive direction is asserted too, in both places: a test that only
//! checks warnings can be satisfied by warning about everything, and a test that
//! only checks scaling can be satisfied by a model that ignores the audit.

use fairchild_core::dc_op_nr;
use fairchild_core::unmodelled;
use fairchild_parser::parse_spice;

fn current(src: &str, vsrc: &str) -> f64 {
    let netlist = parse_spice(src).expect("parse");
    let r = dc_op_nr(&netlist).expect("solve");
    r.vsrc_current(vsrc).expect("source current")
}

/// `AREA` is how a deck scales a device without writing a second model card —
/// parallel output diodes, ESD structures, a scaled current mirror.
#[test]
fn area_scales_the_diode_exactly() {
    const ONE: &str = "* one\n.model dm D (IS=1e-14 N=1)\nV1 a 0 DC 0.7\nD1 a 0 dm\n.op\n";
    const TWO: &str =
        "* area=2\n.model dm D (IS=1e-14 N=1)\nV1 a 0 DC 0.7\nD1 a 0 dm area=2\n.op\n";
    const PAIR: &str =
        "* two in parallel\n.model dm D (IS=1e-14 N=1)\nV1 a 0 DC 0.7\nD1 a 0 dm\nD2 a 0 dm\n.op\n";

    let (one, two, pair) = (current(ONE, "v1"), current(TWO, "v1"), current(PAIR, "v1"));
    // Three separate Newton solves, so the comparison is limited by the
    // convergence tolerance rather than by the scaling — 1e-7 is six orders
    // tighter than anything a wrong scaling could land inside (the bug was a
    // factor of exactly 1).
    assert!(
        ((two / one) - 2.0).abs() < 1e-7,
        "area=2 must double the current: {two:.9e} vs {one:.9e}"
    );
    // The anchor that matters: AREA=2 *is* two devices, so it has to agree with
    // two devices rather than merely being twice something.
    assert!(
        (two - pair).abs() <= 1e-7 * pair.abs(),
        "area=2 ({two:.9e}) must equal two parallel diodes ({pair:.9e})"
    );
}

/// AREA divides RS: N junctions in parallel each carry their own series
/// resistance, so the pair and the scaled device must still agree once RS makes
/// the current sub-exponential (where a wrong RS shows up at all).
#[test]
fn area_divides_the_diode_series_resistance() {
    const TWO: &str =
        "* area\n.model dm D (IS=1e-14 N=1 RS=10)\nV1 a 0 DC 1.0\nD1 a 0 dm area=2\n.op\n";
    const PAIR: &str = "* pair\n.model dm D (IS=1e-14 N=1 RS=10)\nV1 a 0 DC 1.0\n\
                        D1 a 0 dm\nD2 a 0 dm\n.op\n";
    let (two, pair) = (current(TWO, "v1"), current(PAIR, "v1"));
    assert!(
        (two - pair).abs() < 1e-7 * pair.abs(),
        "with RS in play, area=2 ({two:.9e}) must still equal two diodes ({pair:.9e})"
    );
}

/// The other half of #26: the BJT arm of `build_devices` did not destructure the
/// element's parameters at all, so nothing on a `Q` line could reach the device.
#[test]
fn area_scales_the_bjt_exactly() {
    const SRC: &str = "* BJT area\n\
                       .model qm NPN (IS=1e-16 BF=100 VAF=50 CJE=2p CJC=1p)\n\
                       Vb b 0 DC 0.7\n\
                       Vc c 0 DC 2\n\
                       Vc2 c2 0 DC 2\n\
                       Q1 c b 0 qm area=2\n\
                       Q2 c2 b 0 qm\n\
                       .op\n";
    let netlist = parse_spice(SRC).expect("parse");
    let r = dc_op_nr(&netlist).expect("solve");
    let (i2, i1) = (
        r.vsrc_current("vc").unwrap(),
        r.vsrc_current("vc2").unwrap(),
    );
    assert!(
        ((i2 / i1) - 2.0).abs() < 1e-7,
        "area=2 must double IC: {i2:.9e} vs {i1:.9e}"
    );
}

/// An unusable value used to be dropped by the element parser, so a typo in a
/// swept parameter read as "this simulator ignores AREA".
#[test]
fn an_instance_parameter_that_cannot_be_read_is_an_error() {
    for src in [
        "* bad diode value\n.model dm D (IS=1e-14)\nV1 a 0 DC 0.7\nD1 a 0 dm area=2x\n.op\n",
        "* bad bjt value\n.model qm NPN (IS=1e-16)\nVb b 0 DC 0.7\nQ1 b b 0 qm area=lots\n.op\n",
        "* bad mosfet value\n.model nm NMOS (VTO=0.7)\nV1 a 0 DC 1\nM1 a a 0 0 nm W=2u L=oops\n.op\n",
    ] {
        let err = parse_spice(src).expect_err("a value that cannot be read must not be dropped");
        let msg = err.to_string();
        assert!(
            msg.contains("2x") || msg.contains("lots") || msg.contains("oops"),
            "the error has to quote the value: {msg}"
        );
    }
}

/// `set_instance_params` reports what it could not honour, which is what the
/// registry turns into a diagnostic. Both directions: AREA is consumed, and a
/// parameter this model does not implement comes back named.
#[test]
fn a_device_reports_the_instance_parameters_it_cannot_honour() {
    use fairchild_core::models::GummelPoonBjt;
    let (mut q, _) = GummelPoonBjt::from_model_params(false, &[("is".into(), 1e-16)]);
    let unknown = q.set_instance_params(&[("area".into(), 2.0), ("banana".into(), 3.0)]);
    assert_eq!(unknown, vec!["banana".to_string()]);

    // A negative or zero AREA is not a scaling — it is reported rather than
    // silently producing a dead device.
    let (mut q, _) = GummelPoonBjt::from_model_params(false, &[("is".into(), 1e-16)]);
    assert_eq!(q.set_instance_params(&[("area".into(), 0.0)]), vec!["area"]);
}

/// The classification is a pure function so it can be asserted here rather than
/// by scraping stderr, and it says what the deck loses — `unknown parameter IKF`
/// does not tell a user whether it matters to them.
#[test]
fn the_unmodelled_report_names_the_parameter_and_the_consequence() {
    let card = vec![
        ("is".to_string(), 1e-16),
        ("ikf".to_string(), 1e-3),  // modelled — must NOT be reported
        ("kf".to_string(), 1e-15),  // modelled now too — must NOT be reported
        ("xcjc".to_string(), 0.5),  // not modelled
        ("cjs".to_string(), 1e-12), // modelled now too — must NOT be reported
        ("fcs".to_string(), 0.5),   // not modelled, and ngspice ignores it too
    ];
    let lines = unmodelled::report(unmodelled::BJT, &card);
    assert_eq!(lines.len(), 2, "expected XCJC and FCS only, got {lines:?}");
    assert!(lines.iter().any(|l| l.starts_with("XCJC ignored:")));
    assert!(lines.iter().any(|l| l.starts_with("FCS ignored:")));
    assert!(
        lines.iter().all(|l| l.len() > 30),
        "a diagnostic that does not say what was lost is not much better than \
         silence: {lines:?}"
    );
    // The positive direction: nothing is reported for a card that only sets
    // parameters the model actually stamps.
    assert!(unmodelled::report(
        unmodelled::BJT,
        &[
            ("is".into(), 1e-16),
            ("bf".into(), 100.0),
            ("vaf".into(), 50.0)
        ]
    )
    .is_empty());
}

/// `docs/model_status.md` is an audit, so it goes stale silently — the worst
/// property a contract can have. This is the check that stops it: every
/// parameter on an "accepted, not modelled" row in the document must be on the
/// matching table in `crate::unmodelled`, and vice versa.
#[test]
fn the_unmodelled_tables_match_the_audit_document() {
    const MARKER: &str = "accepted, not modelled";
    let doc = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../docs/model_status.md")
            .canonicalize()
            .expect("docs/model_status.md"),
    )
    .expect("read model_status.md");

    // Section headings are stable ("## 3. Diode …"); the table rows in between
    // are what this reads.
    let sections: Vec<(&str, unmodelled::Unmodelled)> = vec![
        ("## 3. Diode", unmodelled::DIODE),
        ("## 4. BJT", unmodelled::BJT),
        ("## 5. MOSFET", unmodelled::MOSFET),
    ];
    for (heading, table) in sections {
        let start = doc
            .find(heading)
            .unwrap_or_else(|| panic!("{heading} is gone from model_status.md"));
        let rest = &doc[start + heading.len()..];
        let end = rest.find("\n## ").unwrap_or(rest.len());
        let mut documented: Vec<String> = Vec::new();
        for line in rest[..end].lines() {
            if !line.contains(MARKER) {
                continue;
            }
            // Every parameter on the row is backticked; the row's prose is not.
            let cell = line.split('|').nth(1).unwrap_or("");
            for tok in cell.split('`').skip(1).step_by(2) {
                documented.push(tok.trim().to_lowercase());
            }
        }
        documented.sort();
        let mut listed: Vec<String> = table.iter().map(|(k, _)| k.to_string()).collect();
        listed.sort();
        assert_eq!(
            documented, listed,
            "{heading}: docs/model_status.md and crate::unmodelled disagree. \
             A contract that disagrees with the code is worse than none — update \
             both in the same commit."
        );
    }
}

/// MJSW was parsed, stored, and never read: the sidewall junction capacitance
/// used MJ, so a card setting them differently did not get what it asked for.
///
/// Anchored on the closed form rather than on the model's own other branch: with
/// only a sidewall (CJ=0, CJSW≠0) the cap at reverse bias is
/// `CJSW·PS·(1 − V/PB)^−MJSW`, and MJ must not appear in it.
#[test]
fn mjsw_grades_the_sidewall_and_mj_does_not() {
    use fairchild_core::device::{Device, EvalFlags, SimContext};
    use fairchild_core::models::Mosfet1;

    const CJSW: f64 = 1e-9; // F/m
    const PS: f64 = 20e-6; // m
    const PB: f64 = 0.8;
    const MJ: f64 = 0.5;
    const MJSW: f64 = 0.2; // deliberately far from MJ
    const VBS: f64 = -1.0;

    let (mut m, _) = Mosfet1::from_model_params(
        false,
        &[
            ("vto".into(), 0.7),
            ("kp".into(), 1e-4),
            ("cj".into(), 0.0), // bottom removed: only the sidewall is left
            ("cjsw".into(), CJSW),
            ("pb".into(), PB),
            ("mj".into(), MJ),
            ("mjsw".into(), MJSW),
        ],
    );
    m.set_instance_params(&[
        ("w".into(), 2e-6),
        ("l".into(), 1e-6),
        ("ps".into(), PS),
        ("as".into(), 50e-12),
    ]);
    let ctx = SimContext::default();
    m.setup_model(&ctx);
    // D=0, G=1, S=2, B=3 — bulk at −1 V against a grounded source.
    m.setup_instance(&[Some(0), Some(1), Some(2), Some(3)], &ctx);
    m.eval(&[0.0, 0.0, 0.0, VBS], EvalFlags::tran(), &ctx);

    let expect = CJSW * PS * (1.0 - VBS / PB).powf(-MJSW);
    let with_mj = CJSW * PS * (1.0 - VBS / PB).powf(-MJ);
    let got = m.cbs_at_last_eval();
    assert!(
        (got - expect).abs() < 1e-12 * expect,
        "sidewall cap {got:.6e} should be {expect:.6e} (MJSW={MJSW}); \
         with MJ it would be {with_mj:.6e}"
    );
}
