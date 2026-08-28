//! `.temp` against ngspice, for all three native model families.
//!
//! Before this, `.temp` scaled `kT/q` and nothing else. A `.temp 125` run gave
//! the thermal voltage at 125 °C and every device parameter at nominal, so:
//!
//! * a silicon diode leaked its 27 °C current where it should leak ~9e4 times it;
//! * a BJT's collector current at fixed `V_BE` was wrong by the same factor;
//! * a Level 1 MOSFET returned the **bit-identical** 27 °C drain current, because
//!   its DC current does not use `vt` at all.
//!
//! The warnings fired for `TNOM`/`EG`/`XTI`, so none of it was silent — but
//! "temperature sweeps work" was a much stronger claim than what was true (#77
//! §5).
//!
//! # What each test is anchored on
//!
//! ngspice, at runtime, because the SPICE literature carries several variants of
//! each law and a PDK relies on the one the reference simulator ships. Every
//! constant in `crate::temperature` was back-solved from these decks.
//!
//! Two traps this file is shaped around:
//!
//! * **A netlist's first line is a title.** A probe deck whose first line was
//!   `.temp 125` had the card silently eaten, and the "measurement" showed
//!   ngspice ignoring temperature entirely for the BJT and the MOSFET. Every deck
//!   here starts with a comment.
//! * **Isolating one law at a time.** A forward current at fixed bias moves with
//!   `IS(T)` *and* `vt(T)`; a MOSFET's drain current moves with `KP(T)` *and*
//!   `VTO(T)`. Each test below picks a probe that separates them.
//!
//! Requires ngspice on PATH; skipped, not failed, without it.

use std::process::Command;

use fairchild_core::{dc_op_nr_with_registry_opts, options::SimOptions, DeviceRegistry};
use fairchild_parser::parse_spice;

/// Solve one deck in fairchild. `probe` is a source name for a current.
fn fairchild(body: &str, vsrc: &str) -> f64 {
    let src = format!("* temperature\n{body}.op\n");
    let net = parse_spice(&src).expect("parse");
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);
    let opts = SimOptions::from_netlist(&net);
    dc_op_nr_with_registry_opts(&net, &reg, &opts)
        .unwrap_or_else(|e| panic!("fairchild failed on\n{src}\n{e:?}"))
        .vsrc_current(vsrc)
        .expect("source current")
}

fn ngspice(body: &str, vsrc: &str) -> Option<f64> {
    let dir = std::env::temp_dir().join("fc_temp_golden");
    std::fs::create_dir_all(&dir).ok()?;
    let tag: String = body
        .chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .take(48)
        .collect();
    let path = dir.join(format!("t_{tag}.sp"));
    // The leading comment matters: SPICE reads a netlist's first line as its
    // title, so a deck starting with `.temp` loses the card.
    std::fs::write(
        &path,
        format!("* temperature\n{body}.control\nop\nprint i({vsrc})\n.endc\n.end\n"),
    )
    .ok()?;
    let out = Command::new("ngspice").arg("-b").arg(&path).output().ok()?;
    let text = String::from_utf8_lossy(&out.stdout);
    for line in text.lines() {
        if line.trim_start().starts_with(&format!("i({vsrc})")) {
            let (_, rhs) = line.split_once('=')?;
            if let Ok(x) = rhs.split_whitespace().next()?.parse::<f64>() {
                return Some(x);
            }
        }
    }
    None
}

fn agree(what: &str, body: &str, vsrc: &str, rel: f64) {
    let fc = fairchild(body, vsrc);
    let Some(ng) = ngspice(body, vsrc) else {
        eprintln!("ngspice not available — skipping '{what}'");
        return;
    };
    let err = (fc - ng).abs() / ng.abs().max(1e-300);
    assert!(
        err < rel,
        "{what}: fairchild {fc:.9e}, ngspice {ng:.9e} (rel {err:.2e} > {rel:.0e})"
    );
}

const TEMPS: [f64; 4] = [-40.0, 27.0, 75.0, 125.0];

// ---------------------------------------------------------------------------
// Diode
// ---------------------------------------------------------------------------

/// `IS(T)`, read where nothing else can move it.
///
/// A diode at −1 V carries `−IS(T)` with the exponential already dead, so this is
/// the saturation current directly rather than a current that also moved because
/// `vt` did. `gmin=0` keeps the leakage floor out of it.
#[test]
fn diode_saturation_current_follows_temperature() {
    for tc in TEMPS {
        let body = format!(
            ".options gmin=0\n.temp {tc}\n.model dm D (IS=1e-14 N=1)\n\
             V1 a 0 DC 1\nD1 0 a dm\n"
        );
        agree(&format!("diode IS at {tc} C"), &body, "v1", 1e-3);
    }
}

/// `EG` and `XTI` are read from the card, not hardcoded.
///
/// Both defaults are silicon's, so a model that ignored the card entirely would
/// pass the sweep above. Non-default values are the only way to see that.
#[test]
fn diode_eg_and_xti_come_from_the_card() {
    for (eg, xti) in [(0.69, 3.0), (1.11, 2.0), (0.69, 2.0)] {
        let body = format!(
            ".options gmin=0\n.temp 125\n.model dm D (IS=1e-14 N=1 EG={eg} XTI={xti})\n\
             V1 a 0 DC 1\nD1 0 a dm\n"
        );
        agree(&format!("diode EG={eg} XTI={xti}"), &body, "v1", 1e-3);
    }
}

/// `TNOM` re-references the card. A card extracted at 125 °C and run at 125 °C
/// must give its nominal `IS` — the ratio is `T/TNOM`, not `T/300.15`.
#[test]
fn diode_tnom_re_references_the_card() {
    let body = ".options gmin=0\n.temp 125\n.model dm D (IS=1e-14 N=1 TNOM=125)\n\
                V1 a 0 DC 1\nD1 0 a dm\n";
    let fc = fairchild(body, "v1").abs();
    assert!(
        (fc - 1e-14).abs() < 1e-16,
        "run at the card's own TNOM, IS must be the card's 1e-14, not {fc:.6e}"
    );
    agree("diode TNOM=125 at 125 C", body, "v1", 1e-3);
}

/// The emission coefficient divides *both* temperature terms in the diode law,
/// which is the one thing that distinguishes it from the BJT's. At `N = 1` the
/// two forms coincide, so this is the test that separates them.
#[test]
fn the_diode_law_divides_by_n() {
    for n in [1.0, 1.5, 2.0] {
        let body = format!(
            ".options gmin=0\n.temp 125\n.model dm D (IS=1e-14 N={n})\n\
             V1 a 0 DC 1\nD1 0 a dm\n"
        );
        agree(&format!("diode N={n} at 125 C"), &body, "v1", 1e-3);
    }
}

// ---------------------------------------------------------------------------
// BJT
// ---------------------------------------------------------------------------

/// `IS(T)` through the collector current at fixed `V_BE`.
///
/// This moves with `IS(T)` and `vt(T)` together, which is fine against an
/// external anchor: ngspice moves both too, so agreement pins the product. The
/// diode tests above are where each is separated.
#[test]
fn bjt_collector_current_follows_temperature() {
    for tc in TEMPS {
        let body = format!(
            ".temp {tc}\n.model qm NPN (IS=1e-16 BF=100)\n\
             VB b 0 DC 0.65\nVC c 0 DC 2\nQ1 c b 0 qm\n"
        );
        agree(&format!("BJT IC at {tc} C"), &body, "vc", 1e-3);
    }
}

/// `XTB` moves the betas and nothing else.
///
/// Beta is `IC/IB`, and both currents carry `IS(T)`, so the ratio isolates
/// `XTB`. `XTB=0` is the control: the ratio must not move at all.
#[test]
fn bjt_beta_follows_xtb() {
    for xtb in [0.0, 1.5] {
        let mut betas = Vec::new();
        for tc in [27.0, 125.0] {
            let body = format!(
                ".temp {tc}\n.model qm NPN (IS=1e-16 BF=100 XTB={xtb})\n\
                 VB b 0 DC 0.65\nVC c 0 DC 2\nQ1 c b 0 qm\n"
            );
            agree(&format!("BJT XTB={xtb} IC at {tc} C"), &body, "vc", 1e-3);
            agree(&format!("BJT XTB={xtb} IB at {tc} C"), &body, "vb", 2e-3);
            betas.push(fairchild(&body, "vc") / fairchild(&body, "vb"));
        }
        let ratio = betas[1] / betas[0];
        let want = (398.15_f64 / 300.15).powf(xtb);
        assert!(
            (ratio - want).abs() / want < 2e-3,
            "XTB={xtb}: beta(125 C)/beta(27 C) is {ratio:.4} and (T/TNOM)^XTB is \
             {want:.4}. At XTB=0 a ratio away from 1 means the beta picked up a \
             temperature dependence it was not given."
        );
    }
}

// ---------------------------------------------------------------------------
// MOSFET
// ---------------------------------------------------------------------------

/// The mobility law and the threshold shift, separated.
///
/// `sqrt(Id)` is linear in `Vgs − Vth` in saturation, so two gate voltages give
/// both the slope (`KP(T)`) and the intercept (`VTO(T)`). Asserting the current
/// at one `Vgs` cannot tell a mobility error from a threshold error — they trade
/// off against each other, which is exactly how a wrong pair passes.
#[test]
fn mosfet_mobility_and_threshold_follow_temperature() {
    for tc in TEMPS {
        let deck = |vgs: f64| {
            format!(
                ".temp {tc}\n.model nm NMOS (VTO=0.7 KP=100u)\n\
                 VG g 0 DC {vgs}\nVD d 0 DC 3\nM1 d g 0 0 nm W=10u L=1u\n"
            )
        };
        agree(
            &format!("MOSFET Id(Vgs=1.2) at {tc} C"),
            &deck(1.2),
            "vd",
            1e-3,
        );
        agree(
            &format!("MOSFET Id(Vgs=1.8) at {tc} C"),
            &deck(1.8),
            "vd",
            1e-3,
        );

        // And the pair really does separate them: fit and compare to the law.
        let (i12, i18) = (
            fairchild(&deck(1.2), "vd").abs(),
            fairchild(&deck(1.8), "vd").abs(),
        );
        let r = (i18 / i12).sqrt();
        let vth = (1.8 - r * 1.2) / (1.0 - r);
        let t = tc + 273.15;
        let phi = fairchild_core::temperature::scaled_phi(0.6, t, 300.15);
        let want = fairchild_core::temperature::scaled_vto(0.7, 0.0, 0.6, phi, t, 300.15, false);
        assert!(
            (vth - want).abs() < 2e-3,
            "at {tc} C the fitted threshold is {vth:.4} V and the law gives \
             {want:.4} V. A mismatch here with both currents agreeing would mean \
             the mobility and the threshold are compensating for each other."
        );
    }
}

/// `GAMMA` does not shift with temperature, but it multiplies `sqrt(PHI(T))`, so
/// the body effect moves anyway — and `sqrt(PHI(T))` enters the threshold twice,
/// where a sign error cancels at `GAMMA = 0`.
#[test]
fn the_body_effect_moves_with_temperature() {
    for tc in [27.0, 125.0] {
        let body = format!(
            ".temp {tc}\n.model nm NMOS (VTO=0.7 KP=100u GAMMA=0.5 PHI=0.6)\n\
             VG g 0 DC 1.8\nVD d 0 DC 3\nVB bk 0 DC -1\n\
             M1 d g 0 bk nm W=10u L=1u\n"
        );
        agree(
            &format!("MOSFET GAMMA=0.5, Vbs=-1 at {tc} C"),
            &body,
            "vd",
            2e-3,
        );
    }
}

/// A deck with no `.temp` must be untouched. Every existing golden is this case,
/// so a law that is not the identity at nominal moves all of them at once — and
/// that is the one property no cross-simulator comparison can catch, because both
/// simulators would share the offset.
#[test]
fn a_deck_without_temp_is_unchanged() {
    let cases = [
        (
            ".options gmin=0\n.model dm D (IS=1e-14 N=1)\nV1 a 0 DC 1\nD1 0 a dm\n",
            "v1",
            1e-14,
        ),
        (
            ".model nm NMOS (VTO=0.7 KP=100u)\nVG g 0 DC 1.8\nVD d 0 DC 3\n\
             M1 d g 0 0 nm W=10u L=1u\n",
            "vd",
            // 0.5 * 100u * 10 * (1.8 - 0.7)^2
            6.05e-4,
        ),
    ];
    for (body, vsrc, want) in cases {
        let got = fairchild(body, vsrc).abs();
        let rel = (got - want).abs() / want;
        assert!(
            rel < 1e-6,
            "with no .temp the answer must be the nominal closed form {want:.6e}, \
             got {got:.6e} (rel {rel:.2e})"
        );
    }
}
