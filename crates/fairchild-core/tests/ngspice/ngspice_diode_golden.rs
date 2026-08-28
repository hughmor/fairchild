//! Golden comparison tests for nonlinear DC (diode circuits).
//! Runs fairchild dc_op_nr and ngspice, asserts agreement within tolerance.
//!
//! Requires ngspice on PATH. Tests are skipped (not failed) when ngspice is absent.

use std::collections::HashMap;
use std::io::Write;
use std::process::Command;

use fairchild_core::{dc_op_nr, SimError};
use fairchild_parser::{parse_spice, Element};

// Slightly relaxed compared to linear tests because ngspice includes gmin
// and temperature-model details not in our Shockley implementation.
const REL_TOL: f64 = 1e-3; // 0.1%
const ABS_TOL_V: f64 = 1e-5; // 10 μV floor

fn find_ngspice() -> Option<std::path::PathBuf> {
    if Command::new("ngspice").arg("--version").output().is_ok() {
        return Some("ngspice".into());
    }
    for candidate in &[
        "/opt/homebrew/bin/ngspice",
        "/usr/local/bin/ngspice",
        "/usr/bin/ngspice",
    ] {
        let p = std::path::Path::new(candidate);
        if p.exists() {
            return Some(p.to_owned());
        }
    }
    None
}

fn strip_control_and_end(netlist: &str) -> String {
    let mut out = String::new();
    let mut in_control = false;
    for line in netlist.lines() {
        let lc = line.trim().to_lowercase();
        if lc.starts_with(".control") {
            in_control = true;
            continue;
        }
        if in_control {
            if lc.starts_with(".endc") {
                in_control = false;
            }
            continue;
        }
        if lc == ".end" {
            continue;
        }
        out.push_str(line);
        out.push('\n');
    }
    out
}

fn parse_ngspice_print(output: &str) -> Option<HashMap<String, f64>> {
    let mut map = HashMap::new();
    for line in output.lines() {
        let line = line.trim();
        if let Some((lhs, rhs)) = line.split_once('=') {
            let key = lhs.trim().to_lowercase();
            if key.starts_with("v(") || key.starts_with("i(") {
                if let Ok(val) = rhs.trim().parse::<f64>() {
                    map.insert(key, val);
                }
            }
        }
    }
    if map.is_empty() {
        None
    } else {
        Some(map)
    }
}

fn ngspice_op(netlist: &str, queries: &[&str]) -> Option<HashMap<String, f64>> {
    let ngspice_bin = find_ngspice()?;
    let mut tmp = tempfile::NamedTempFile::new().ok()?;
    let stripped = strip_control_and_end(netlist);
    let print_vars = queries.join(" ");
    let control_block = format!(".control\nop\nprint {print_vars}\n.endc\n");
    write!(tmp, "{stripped}\n{control_block}").ok()?;
    let output = Command::new(&ngspice_bin)
        .arg("-b")
        .arg(tmp.path())
        .output()
        .ok()?;
    let stdout = String::from_utf8_lossy(&output.stdout);
    parse_ngspice_print(&stdout)
}

fn fairchild_nr_op(netlist_str: &str) -> Result<HashMap<String, f64>, SimError> {
    let netlist = parse_spice(netlist_str).map_err(SimError::Parse)?;
    let result = dc_op_nr(&netlist)?;
    let mut map = HashMap::new();
    for (name, v) in result.all_voltages() {
        map.insert(format!("v({name})"), v);
    }
    for el in &netlist.elements {
        if let Element::VoltageSource { name, .. } = el {
            if let Ok(i) = result.vsrc_current(name) {
                map.insert(format!("i({name})"), i);
            }
        }
    }
    Ok(map)
}

macro_rules! diode_golden_test {
    ($name:ident, $netlist_file:expr, $queries:expr, $tol_rel:expr) => {
        #[test]
        fn $name() {
            let netlist_str = include_str!(concat!("../../../../tests/golden/", $netlist_file));

            let fc = fairchild_nr_op(netlist_str).expect("fairchild NR solve failed");

            let queries: &[&str] = $queries;
            let Some(ng) = ngspice_op(netlist_str, queries) else {
                eprintln!("ngspice not available — skipping golden comparison");
                assert!(!fc.is_empty(), "fairchild produced no results");
                return;
            };

            for key in queries {
                let key_lc = key.to_lowercase();
                let fc_val = fc
                    .get(&key_lc)
                    .copied()
                    .unwrap_or_else(|| panic!("fairchild missing '{key_lc}'; available: {fc:?}"));
                let ng_val = ng
                    .get(&key_lc)
                    .copied()
                    .unwrap_or_else(|| panic!("ngspice missing '{key_lc}'; available: {ng:?}"));
                let tol = f64::max(ABS_TOL_V, $tol_rel * ng_val.abs());
                assert!(
                    (fc_val - ng_val).abs() <= tol,
                    "{key_lc}: fairchild={fc_val:.6e}  ngspice={ng_val:.6e}  \
                     diff={:.2e}  tol={tol:.2e}",
                    (fc_val - ng_val).abs()
                );
            }
        }
    };
}

diode_golden_test!(diode_current_source_bias, "diode_iv.sp", &["v(b)"], REL_TOL);

diode_golden_test!(
    diode_series_rd,
    "diode_series_rd.sp",
    &["v(b)", "i(vdd)"],
    REL_TOL
);

/// `RS` in the model card is a series resistance, so it must equal an explicit
/// resistor of the same value — and both must equal ngspice.
///
/// # What this caught
///
/// `RS` was eliminated by a *lagged* fixed point: `vd_j = vd_terminal − Id·RS`
/// with `Id` from the previous eval, one step per outer Newton iteration. The
/// junction voltage is internal state the outer Newton cannot see, so its
/// convergence test could be satisfied with the lag still wide open — and is,
/// immediately, whenever a voltage source pins the diode's terminals and the
/// visible unknowns stop moving on iteration one. It read 2.7% low at 0.7 V, and
/// 100% wrong at 1.0 V once `gmin` became a real conductance across the junction
/// and gave the reverse-biased branch something to converge onto.
///
/// # Why both halves
///
/// The self-consistency half — internal `RS` against an external resistor —
/// cannot detect a fault common to both forms, so ngspice anchors one side
/// absolutely. ngspice has no lag to hide: it gives the intrinsic junction its
/// own internal node (`dio_posPrime`) and lets the matrix solve for it.
///
/// 1.0 V is in the sweep because that is where the series drop dominates and the
/// lag was largest; 0.4 V is below the knee, where it was invisible.
#[test]
fn rs_in_the_model_equals_an_external_resistor() {
    // 1e-4 is ~30x above the convergence noise between two separate solves and
    // ~270x below the smallest error the lag produced.
    const TOL: f64 = 1e-4;
    for v in [0.4, 0.7, 1.0] {
        let internal = format!(
            "* rs internal\n.model dm D (IS=1e-14 N=1 RS=10)\nV1 a 0 DC {v}\nD1 a 0 dm\n.op\n"
        );
        // Same circuit, drawn instead of parameterised.
        let external = format!(
            "* rs external\n.model dm D (IS=1e-14 N=1)\nV1 a 0 DC {v}\nR1 a m 10\nD1 m 0 dm\n.op\n"
        );

        let pick = |deck: &str| {
            fairchild_nr_op(deck)
                .unwrap_or_else(|e| panic!("V={v}: fairchild failed on\n{deck}\n{e:?}"))["i(v1)"]
        };
        let (int_i, ext_i) = (pick(&internal), pick(&external));
        let rel = (int_i - ext_i).abs() / ext_i.abs();
        assert!(
            rel < TOL,
            "V={v}: RS=10 in the model gives i(v1)={int_i:.9e} but an external \
             10 Ohm gives {ext_i:.9e} (rel {rel:.2e}) — the same circuit. A lag \
             in the RS elimination reads low, and reads low by more the further \
             the series drop dominates."
        );

        let Some(ng) = ngspice_op(&internal, &["i(v1)"]) else {
            eprintln!("ngspice not available — self-consistency checked, anchor skipped");
            continue;
        };
        let ng_i = ng["i(v1)"];
        let rel_ng = (int_i - ng_i).abs() / ng_i.abs();
        assert!(
            rel_ng < TOL,
            "V={v}: RS=10 gives i(v1)={int_i:.9e}, ngspice {ng_i:.9e} \
             (rel {rel_ng:.2e}). Both forms agreeing with each other and not \
             with ngspice would mean the fault is in the junction, not the \
             elimination."
        );
    }
}

// ---------------------------------------------------------------------------
// Recombination, high injection, and RS with temperature (#97 section 5)
// ---------------------------------------------------------------------------

/// A DC solve that honours the deck's `.options` and `.temp`.
fn d5_current(deck: &str) -> f64 {
    let net = fairchild_parser::parse_spice(deck).expect("parse");
    let mut registry = fairchild_core::DeviceRegistry::new();
    registry.register_builtin_models(&net.models);
    let opts = fairchild_core::options::SimOptions::from_netlist(&net);
    fairchild_core::dc_op_nr_with_registry_opts(&net, &registry, &opts)
        .unwrap_or_else(|e| panic!("solve failed on\n{deck}\n{e:?}"))
        .vsrc_current("v1")
        .expect("i(v1)")
        .abs()
}

/// ngspice's `i(v1)` for the same deck.
fn d5_ngspice(deck: &str) -> Option<f64> {
    let ng = find_ngspice()?;
    let mut tmp = tempfile::NamedTempFile::new().ok()?;
    use std::io::Write as _;
    write!(tmp, "{deck}.control\nop\nprint i(v1)\n.endc\n.end\n").ok()?;
    let out = std::process::Command::new(&ng)
        .arg("-b")
        .arg(tmp.path())
        .output()
        .ok()?;
    let text = String::from_utf8_lossy(&out.stdout);
    for line in text.lines() {
        let t = line.trim();
        if t.starts_with("i(v1)") && t.contains('=') {
            if let Ok(v) = t
                .split('=')
                .nth(1)?
                .split_whitespace()
                .next()?
                .parse::<f64>()
            {
                return Some(v.abs());
            }
        }
    }
    None
}

fn d5_deck(model: &str, v: f64, temp: Option<f64>) -> String {
    let t = temp.map(|t| format!(".temp {t}\n")).unwrap_or_default();
    format!("* diode group\n.options gmin=0\n{t}.model dm D ({model})\nV1 a 0 DC {v}\nD1 a 0 dm\n")
}

/// `ISR` adds the recombination current, which dominates below about 0.4 V.
///
/// Without it the low-bias forward current is the ideal exponential and nothing
/// else, so a real diode's low-current ideality — nearer two than one — cannot be
/// described at all. Measured law, with SPICE's generation factor:
///
/// ```text
/// Irec = ISR·(exp(V/(NR·vt)) − 1)·((1 − V/VJ)² + 0.005)^(M/2)
/// ```
///
/// The `VJ`/`M` sweep is what pins the generation factor: at 0.2 V with
/// `ISR = 1e-10` the recombination current is 184× the ideal one, and moving `M` to
/// zero changes it by 12%, so a missing or wrong factor cannot pass.
#[test]
fn isr_adds_the_recombination_current() {
    let mut compared = 0;
    for (vj, m) in [(1.0, 0.5), (0.75, 0.33), (0.6, 0.5), (1.0, 0.0)] {
        for v in [0.2, 0.35, 0.5] {
            let model = format!("IS=1e-14 N=1 ISR=1e-10 VJ={vj} M={m}");
            let deck = d5_deck(&model, v, None);
            let got = d5_current(&deck);
            let Some(ng) = d5_ngspice(&deck) else {
                eprintln!("ngspice not available — skipping");
                return;
            };
            let rel = (got - ng).abs() / ng;
            assert!(
                rel < 2e-3,
                "VJ={vj} M={m} V={v}: fairchild {got:.6e}, ngspice {ng:.6e} \
                 (rel {rel:.2e}). Dropping the generation factor moves this by \
                 10 to 30%; dropping ISR entirely by a factor of 184 at 0.2 V."
            );
            compared += 1;
        }
    }
    assert_eq!(compared, 12, "every VJ/M and bias must have been compared");
}

/// `IKF` bends the forward current from exponential towards `sqrt`, and it acts on
/// the **total** current rather than the ideal part alone.
///
/// ```text
/// Id = Id_total/(1 + sqrt(Id_total/IKF))
/// ```
///
/// Measured to 1e-6 at twelve points. The second half of the test is the part that
/// needed measuring: with `ISR` and `IKF` both on the card, ngspice reads
/// 8.555990e-05 at 0.5 V, against 8.555977e-05 for the knee on the total and
/// 1.143950e-04 for the knee on the ideal current alone. A 34% difference, and the
/// wrong choice looks perfectly reasonable.
#[test]
fn ikf_bends_the_forward_current_and_acts_on_the_total() {
    let mut compared = 0;
    for ikf in ["1e-6", "1e-3", "1e-1"] {
        for v in [0.5, 0.6, 0.75, 0.9] {
            let deck = d5_deck(&format!("IS=1e-14 N=1 IKF={ikf}"), v, None);
            let got = d5_current(&deck);
            let Some(ng) = d5_ngspice(&deck) else {
                eprintln!("ngspice not available — skipping");
                return;
            };
            let rel = (got - ng).abs() / ng;
            assert!(
                rel < 2e-3,
                "IKF={ikf} V={v}: fairchild {got:.6e}, ngspice {ng:.6e} \
                 (rel {rel:.2e}). Ignoring IKF is a factor of 3600 at 0.9 V."
            );
            compared += 1;
        }
    }
    // The knee on the total, not on the ideal alone.
    for v in [0.5, 0.6] {
        let deck = d5_deck("IS=1e-14 N=1 ISR=1e-8 IKF=1e-3", v, None);
        let got = d5_current(&deck);
        let Some(ng) = d5_ngspice(&deck) else { return };
        let rel = (got - ng).abs() / ng;
        assert!(
            rel < 2e-3,
            "V={v} with both ISR and IKF: fairchild {got:.6e}, ngspice {ng:.6e} \
             (rel {rel:.2e}). Applying the knee to the ideal current alone and \
             adding the recombination afterwards reads 1.143950e-04 at 0.5 V \
             where ngspice reads 8.555990e-05."
        );
        compared += 1;
    }
    assert_eq!(compared, 14, "every point must have been compared");
}

/// `TRS1` and `TRS2` move `RS` with temperature.
///
/// `RS(T) = RS·(1 + TRS1·dT + TRS2·dT²)`, `dT` in kelvin from `TNOM`. Read at a
/// hard forward drive where the series drop dominates the answer.
///
/// # The coefficients are sized to the span
///
/// `RS(T) = RS·(1 + TRS1·dT + TRS2·dT²)` goes **negative** for a large enough `|dT|`
/// of the wrong sign, and a negative series resistance is unphysical. `TRS1=1e-2`
/// with `TNOM=75` at −40 °C gives `1 − 1.15 = −0.15`, and ngspice answers "DC
/// solution failed" rather than a number. So the grid below uses `TRS1=1e-3`, which
/// stays positive across the whole 165 K span, and
/// [`a_negative_series_resistance_is_refused_by_name`] covers the other case.
///
/// # `TNOM` is swept, not fixed at 27 °C
///
/// A card with `TNOM=27` cannot tell `dT = T − TNOM` from `dT = T − 300.15 K`,
/// because those are the same temperature. So `TNOM=75` is here too, and the
/// zero-`dT` check runs at both — the coefficients must do nothing at each card's
/// own `TNOM`. Referencing `dT` to 300.15 K passes the whole `TNOM=27` sweep and
/// fails that.
#[test]
fn trs1_and_trs2_move_the_series_resistance_with_temperature() {
    let mut compared = 0;
    for tnom in [27.0, 75.0] {
        for temp in [-40.0, 0.0, 27.0, 75.0, 125.0] {
            for card in ["", "TRS1=1e-3", "TRS2=1e-6", "TRS1=1e-3 TRS2=1e-6"] {
                let deck = d5_deck(
                    &format!("IS=1e-14 N=1 RS=100 TNOM={tnom} {card}"),
                    2.0,
                    Some(temp),
                );
                let got = d5_current(&deck);
                let Some(ng) = d5_ngspice(&deck) else {
                    eprintln!("ngspice not available — skipping");
                    return;
                };
                let rel = (got - ng).abs() / ng;
                assert!(
                    rel < 2e-3,
                    "TNOM={tnom} T={temp} '{card}': fairchild {got:.6e}, \
                     ngspice {ng:.6e} (rel {rel:.2e})"
                );
                compared += 1;
            }
        }
    }
    assert_eq!(
        compared, 40,
        "every TNOM, temperature and card must have been compared"
    );

    // At each card's own TNOM, dT is zero and the coefficients must do nothing.
    // `TNOM=75` is the row that matters: with `TNOM=27` a law referenced to
    // 300.15 K is indistinguishable, because that *is* 27 C.
    for tnom in [27.0, 75.0] {
        let plain = d5_current(&d5_deck(
            &format!("IS=1e-14 N=1 RS=100 TNOM={tnom}"),
            2.0,
            Some(tnom),
        ));
        let with_tc = d5_current(&d5_deck(
            &format!("IS=1e-14 N=1 RS=100 TNOM={tnom} TRS1=1e-3 TRS2=1e-6"),
            2.0,
            Some(tnom),
        ));
        assert!(
            (plain - with_tc).abs() / plain < 1e-12,
            "at TNOM={tnom}, dT is zero and the coefficients must change \
             nothing: {plain:.9e} against {with_tc:.9e}. A law referenced to \
             0 C, or to 300.15 K rather than to TNOM, fails here."
        );
    }
}

/// `NR` is honoured, which is a deliberate divergence from ngspice.
///
/// **ngspice ignores `NR`** — its answer is bit-identical with and without `NR=2`,
/// so it hardcodes the default. Honouring it agrees with ngspice on every card
/// ngspice can represent, and honours one it cannot.
///
/// Both directions are asserted. `NR=2` and an absent `NR` must match ngspice, and
/// `NR=1` must move — by a large factor, because it halves the exponent's divisor.
/// Without the second half this would be the worthless `X_is_accepted` shape.
///
/// Same call as the diode's `AREA` in reverse breakdown, where ngspice disagrees
/// with its own parallel pair.
#[test]
fn nr_is_honoured_although_ngspice_ignores_it() {
    for v in [0.2, 0.35] {
        let default = d5_current(&d5_deck("IS=1e-14 N=1 ISR=1e-10", v, None));
        let explicit = d5_current(&d5_deck("IS=1e-14 N=1 ISR=1e-10 NR=2", v, None));
        assert!(
            (default - explicit).abs() / default < 1e-12,
            "V={v}: NR defaults to 2, so an absent NR and NR=2 must be \
             identical: {default:.9e} against {explicit:.9e}"
        );
        if let Some(ng) = d5_ngspice(&d5_deck("IS=1e-14 N=1 ISR=1e-10 NR=2", v, None)) {
            let rel = (explicit - ng).abs() / ng;
            assert!(
                rel < 2e-3,
                "V={v}: at NR=2 this must agree with ngspice, which hardcodes \
                 that value: {explicit:.6e} against {ng:.6e}"
            );
        }
        // And it has to actually do something, or honouring it is a claim only.
        let sharper = d5_current(&d5_deck("IS=1e-14 N=1 ISR=1e-10 NR=1", v, None));
        assert!(
            sharper / explicit > 2.0,
            "V={v}: NR=1 halves the exponent's divisor, so the recombination \
             current must rise sharply: {sharper:.6e} against {explicit:.6e}. \
             Equal values would mean NR is accepted and dropped."
        );
    }
}

/// A negative `RS(T)` is refused by name, not stamped.
///
/// `1 + TRS1·dT + TRS2·dT²` goes negative for a large enough `|dT|` of the wrong
/// sign. `TNOM=75` with `TRS1=1e-2` at −40 °C gives `1 − 1.15 = −0.15`, so
/// `RS(T) = −15 Ω`, which feeds the junction rather than dropping across it.
///
/// ngspice's answer there is "DC solution failed" — the right outcome reached by
/// the wrong route, because it names neither the parameter nor the temperature. The
/// error here names `RS`, both coefficients, the resulting `RS(T)` and `TNOM`.
///
/// The positive control is the same card at its own `TNOM`, where the factor is
/// exactly 1 and the deck must solve. Without it this would pass on a diode that
/// refuses every card.
#[test]
fn a_negative_series_resistance_is_refused_by_name() {
    let bad = d5_deck("IS=1e-14 N=1 RS=100 TNOM=75 TRS1=1e-2", 2.0, Some(-40.0));
    let net = fairchild_parser::parse_spice(&bad).expect("parse");
    let mut registry = fairchild_core::DeviceRegistry::new();
    registry.register_builtin_models(&net.models);
    let opts = fairchild_core::options::SimOptions::from_netlist(&net);
    let msg = match fairchild_core::dc_op_nr_with_registry_opts(&net, &registry, &opts) {
        Err(e) => format!("{e:?}"),
        Ok(_) => panic!("a negative RS(T) must be refused, not solved"),
    };
    for needle in ["RS", "TRS1", "TRS2", "TNOM"] {
        assert!(
            msg.contains(needle),
            "the refusal must name {needle}, or a user cannot act on it: {msg}"
        );
    }

    // The control: the same card at its own TNOM has a factor of exactly 1.
    let good = d5_deck("IS=1e-14 N=1 RS=100 TNOM=75 TRS1=1e-2", 2.0, Some(75.0));
    let i = d5_current(&good);
    assert!(
        i > 0.0 && i.is_finite(),
        "the control failed: at TNOM the factor is 1 and this deck must solve, \
         or the refusal above is refusing everything. Got {i}"
    );
}
