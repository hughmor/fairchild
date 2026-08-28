//! Junction potential and capacitance against temperature, versus ngspice.
//!
//! #97 §6, and the last of #77 §5. `TNOM` re-referenced `IS`, `BF`/`BR`, `KP`,
//! `PHI` and `VTO`, and not the junction potentials or capacitances — so a DC
//! operating point was fully temperature-corrected while a transient or AC answer
//! still used nominal capacitances.
//!
//! # How the capacitance is measured
//!
//! A reverse-biased diode's own `Cj` as the C of a single-pole RC: `V1` at
//! `DC 5 AC 1` through 1 MΩ into the cathode holds the junction at −5 V, and the
//! junction capacitance shunts that node to ground. Then
//! `|V(a)/V(in)| = 1/sqrt(1 + (2πfRC)²)`, which inverts to `C` with nothing else
//! in the way.
//!
//! The probe frequency sits near the pole on purpose. A first attempt used 1 kHz
//! against a ~450 kHz pole, where the divider attenuates by 6 ppm and six printed
//! digits carry no information about `C` at all.
//!
//! # Why `M` is swept
//!
//! `M` appears three times in the capacitance law — once in each of the two
//! grading corrections and once in the depletion exponent — so a single value can
//! hide a term applied with the wrong sign or in the wrong place.
//!
//! Requires ngspice on PATH; skipped, not failed, without it.

use std::f64::consts::PI;
use std::process::Command;

use fairchild_core::{ac_analysis_opts, options::SimOptions, DeviceRegistry};
use fairchild_parser::parse_spice;

const R: f64 = 1e6;
const F: f64 = 5e5;
const VREV: f64 = 5.0;

fn deck(model: &str, tc: f64) -> String {
    format!(
        "* junction cap vs temperature\n.temp {tc}\n.model dm D ({model})\n\
         V1 in 0 DC {VREV} AC 1\nR1 in a {R:e}\nD1 0 a dm\n"
    )
}

/// `C` implied by the divider's magnitude at `F`.
fn cap_from_mag(mag: f64) -> f64 {
    (1.0 / (mag * mag) - 1.0).sqrt() / (2.0 * PI * F * R)
}

fn fairchild_cap(model: &str, tc: f64) -> f64 {
    let src = deck(model, tc);
    let net = parse_spice(&src).expect("parse");
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);
    let opts = SimOptions::from_netlist(&net);
    let r = ac_analysis_opts(&net, &[F], Some("v1"), &reg, &opts)
        .unwrap_or_else(|e| panic!("ac failed on\n{src}\n{e:?}"));
    let mag = r.magnitude("a", 0).expect("|V(a)| at the single frequency");
    cap_from_mag(mag)
}

fn ngspice_cap(model: &str, tc: f64) -> Option<f64> {
    let dir = std::env::temp_dir().join("fc_cjt_golden");
    std::fs::create_dir_all(&dir).ok()?;
    let tag: String = model
        .chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .take(40)
        .collect();
    let path = dir.join(format!("cjt_{tag}_{tc}.sp"));
    std::fs::write(
        &path,
        format!(
            "{}.control\nac lin 1 {F:e} {F:e}\nprint vm(a)\n.endc\n.end\n",
            deck(model, tc)
        ),
    )
    .ok()?;
    let out = Command::new("ngspice").arg("-b").arg(&path).output().ok()?;
    let text = String::from_utf8_lossy(&out.stdout);
    for line in text.lines() {
        let t = line.trim();
        // A single-point `print` gives `vm(a) = <value>`, not the indexed table a
        // swept `print` produces.
        if t.starts_with("vm(a)") {
            if let Ok(mag) = t.split('=').nth(1)?.trim().parse::<f64>() {
                if mag > 0.0 && mag < 1.0 {
                    return Some(cap_from_mag(mag));
                }
            }
        }
    }
    None
}

const TEMPS: [f64; 4] = [-40.0, 27.0, 75.0, 125.0];

/// `CJO(T)` and `VJ(T)` on a diode, across temperature and grading coefficient.
#[test]
fn diode_junction_capacitance_follows_temperature() {
    for m in [0.33, 0.5] {
        for tc in TEMPS {
            let model = format!("IS=1e-14 N=1 CJO=1p VJ=0.75 M={m}");
            let fc = fairchild_cap(&model, tc);
            let Some(ng) = ngspice_cap(&model, tc) else {
                eprintln!("ngspice not available — skipping");
                return;
            };
            let rel = (fc - ng).abs() / ng;
            assert!(
                rel < 2e-3,
                "M={m}, {tc} C: fairchild Cj={fc:.6e} F, ngspice {ng:.6e} F \
                 (rel {rel:.2e})"
            );
        }
    }
}

/// The capacitance *moves* with temperature, and in the right direction.
///
/// The comparison above would pass if both simulators held `Cj` fixed. A hotter
/// junction has a narrower built-in potential and so a thinner depletion layer,
/// which is more capacitance.
#[test]
fn a_hotter_junction_has_more_capacitance() {
    let model = "IS=1e-14 N=1 CJO=1p VJ=0.75 M=0.5";
    let (cold, nominal, hot) = (
        fairchild_cap(model, -40.0),
        fairchild_cap(model, 27.0),
        fairchild_cap(model, 125.0),
    );
    assert!(
        cold < nominal && nominal < hot,
        "Cj must rise with temperature: -40 C {cold:.6e}, 27 C {nominal:.6e}, \
         125 C {hot:.6e}. All three equal means the law is not applied."
    );
    // And by a material amount, so this cannot pass on rounding.
    let span = (hot - cold) / nominal;
    assert!(
        span > 0.02,
        "the -40…125 C span is {span:.4} of the nominal capacitance, which is too \
         small to be the law rather than noise"
    );
}

/// `TNOM` re-references the card: extracted at 125 °C and run at 125 °C must give
/// the nominal capacitance, not a shifted one.
///
/// This is the property no cross-simulator comparison can check — both would share
/// an offset — and it is what breaks every existing transient golden at once if
/// the two halves of the law fail to cancel.
#[test]
fn tnom_re_references_the_junction() {
    let at_nominal = fairchild_cap("IS=1e-14 N=1 CJO=1p VJ=0.75 M=0.5", 27.0);
    let card_at_125 = fairchild_cap("IS=1e-14 N=1 CJO=1p VJ=0.75 M=0.5 TNOM=125", 125.0);
    let rel = (card_at_125 - at_nominal).abs() / at_nominal;
    assert!(
        rel < 1e-6,
        "a card extracted at TNOM=125 and run at 125 C must give its nominal \
         capacitance {at_nominal:.9e} F, got {card_at_125:.9e} (rel {rel:.2e})"
    );
    // And against ngspice, which agrees about the identity.
    if let Some(ng) = ngspice_cap("IS=1e-14 N=1 CJO=1p VJ=0.75 M=0.5 TNOM=125", 125.0) {
        let rel_ng = (card_at_125 - ng).abs() / ng;
        assert!(
            rel_ng < 2e-3,
            "TNOM=125 at 125 C: {card_at_125:.6e} vs {ng:.6e}"
        );
    }
}

/// A deck with no `.temp` is untouched — every existing transient golden is this
/// case, so a law that is not the identity at nominal moves all of them.
#[test]
fn a_deck_without_temp_keeps_its_nominal_capacitance() {
    let model = "IS=1e-14 N=1 CJO=1p VJ=0.75 M=0.5";
    let got = fairchild_cap(model, 27.0);
    // The closed form at -5 V: CJO·(1 − V/VJ)^(−M) with V = −5.
    let want = 1e-12 * (1.0 + VREV / 0.75_f64).powf(-0.5);
    let rel = (got - want).abs() / want;
    assert!(
        rel < 1e-3,
        "at nominal the capacitance must be the closed form {want:.6e} F, got \
         {got:.6e} (rel {rel:.2e})"
    );
}
