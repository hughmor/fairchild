//! Golden comparison tests for Gummel-Poon Level 1 BJT circuits.
//! DC tests compare fairchild dc_op_nr against ngspice within 0.5%.
//! Transient test compares BJT CE stage with CJE/CJC caps at one timepoint.
//! Tests are skipped (not failed) when ngspice is absent.

use std::collections::HashMap;
use std::io::Write;
use std::process::Command;

use fairchild_core::{dc_op_nr, options::SimOptions, tran_nr_with_registry_opts, DeviceRegistry};
use fairchild_parser::parse_spice;

// Level 1 GP vs ngspice GP: same equations, same defaults. Allow 0.5 % relative
// to account for minor numerical differences (GMIN, iteration tolerances).
const REL_TOL: f64 = 5e-3;
const ABS_TOL_V: f64 = 1e-4; // 100 µV floor

// ---------------------------------------------------------------------------
// ngspice harness (shared with other golden tests)
// ---------------------------------------------------------------------------

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

fn fairchild_op(netlist_str: &str) -> HashMap<String, f64> {
    let netlist = parse_spice(netlist_str).expect("parse failed");
    let result = dc_op_nr(&netlist).expect("DC solve failed");
    let mut map = HashMap::new();
    for (name, v) in result.all_voltages() {
        map.insert(format!("v({name})"), v);
    }
    for el in &netlist.elements {
        if let fairchild_parser::Element::VoltageSource { name, .. } = el {
            if let Ok(i) = result.vsrc_current(name) {
                map.insert(format!("i({name})"), i);
            }
        }
    }
    map
}

fn assert_close(key: &str, fc: f64, ng: f64) {
    let tol = f64::max(ABS_TOL_V, REL_TOL * ng.abs());
    assert!(
        (fc - ng).abs() <= tol,
        "{key}: fairchild={fc:.6e}  ngspice={ng:.6e}  diff={:.2e}  tol={tol:.2e}",
        (fc - ng).abs()
    );
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn npn_ce_forward_active() {
    // NPN common-emitter amplifier, emitter grounded.
    // VBB=0.8V through RB=10k drives the base; RC=3.3k from VCC=5V to collector.
    // Expected: V(b) ≈ 0.69–0.72V (VBE), V(c) in [1.5, 4.5] (forward active).
    let netlist_str = "\
* NPN CE amplifier — forward active\n\
.model npn1 NPN (IS=1e-15 BF=100 BR=1)\n\
VCC cc 0 DC 5\n\
VBB bb 0 DC 0.8\n\
RB  bb b  10k\n\
RC  cc c  3.3k\n\
Q1  c b 0 0 npn1\n\
.op\n\
.end\n";

    let fc = fairchild_op(netlist_str);
    let vb = fc["v(b)"];
    let vc = fc["v(c)"];

    // Sanity: VBE forward biased, BJT in forward active (VCE > ~0.2V).
    assert!(vb > 0.6 && vb < 0.8, "V(b)={vb:.4}V — expected ~0.7V (VBE)");
    assert!(
        vc > 1.0 && vc < 5.0,
        "V(c)={vc:.4}V — expected forward active [1,5]V"
    );

    let Some(ng) = ngspice_op(netlist_str, &["v(b)", "v(c)"]) else {
        eprintln!("ngspice not available — skipping comparison");
        return;
    };
    assert_close("v(b)", vb, ng["v(b)"]);
    assert_close("v(c)", vc, ng["v(c)"]);
}

#[test]
fn npn_ce_saturation() {
    // Heavy base drive (VBB=2V through RB=10k) forces saturation.
    // Expected: V(c) < 0.5V (BJT is saturated, VCE ≈ VCE_sat ≈ 0.1–0.3V).
    let netlist_str = "\
* NPN CE amplifier — saturation\n\
.model npn1 NPN (IS=1e-15 BF=100 BR=1)\n\
VCC cc 0 DC 5\n\
VBB bb 0 DC 2.0\n\
RB  bb b  10k\n\
RC  cc c  3.3k\n\
Q1  c b 0 0 npn1\n\
.op\n\
.end\n";

    let fc = fairchild_op(netlist_str);
    let vc = fc["v(c)"];

    assert!(vc < 0.5, "V(c)={vc:.4}V — expected saturation (<0.5V)");

    let Some(ng) = ngspice_op(netlist_str, &["v(b)", "v(c)"]) else {
        eprintln!("ngspice not available — skipping comparison");
        return;
    };
    assert_close("v(b)", fc["v(b)"], ng["v(b)"]);
    assert_close("v(c)", vc, ng["v(c)"]);
}

#[test]
fn pnp_ce_forward_active() {
    // PNP common-emitter: emitter tied to VCC=5V, base biased at ~4.3V
    // through RB=10k, collector through RC=3.3k to ground.
    // VEB ≈ 5 − 4.3 = 0.7V → forward active.
    // IC flows from emitter to collector; V(c) = IC * RC ≈ 1–2V.
    let netlist_str = "\
* PNP CE amplifier — forward active\n\
.model pnp1 PNP (IS=1e-15 BF=100 BR=1)\n\
VCC cc 0 DC 5\n\
VBB bb 0 DC 4.3\n\
RB  bb b  10k\n\
RC  c  0  3.3k\n\
Q1  c b cc 0 pnp1\n\
.op\n\
.end\n";

    let fc = fairchild_op(netlist_str);
    let vb = fc["v(b)"];
    let vc = fc["v(c)"];

    // VEB = VCC - VB should be ~0.7V forward biased; VC should be positive.
    let veb = 5.0 - vb;
    assert!(
        veb > 0.5 && veb < 0.85,
        "VEB={veb:.4}V — expected ~0.7V for PNP"
    );
    assert!(vc > 0.5 && vc < 4.5, "V(c)={vc:.4}V — expected [0.5, 4.5]V");

    let Some(ng) = ngspice_op(netlist_str, &["v(b)", "v(c)"]) else {
        eprintln!("ngspice not available — skipping comparison");
        return;
    };
    assert_close("v(b)", vb, ng["v(b)"]);
    assert_close("v(c)", vc, ng["v(c)"]);
}

// ---------------------------------------------------------------------------
// CJE/CJC transient golden test
// ---------------------------------------------------------------------------

/// Parse ngspice batch output for `.meas tran` results.
/// Lines look like:  `v_6n         =  3.456789e+00`
fn ngspice_meas_bjt(netlist: &str) -> Option<HashMap<String, f64>> {
    let ngspice_bin = find_ngspice()?;
    let mut tmp = tempfile::NamedTempFile::new().ok()?;
    write!(tmp, "{}", netlist).ok()?;
    let output = Command::new(&ngspice_bin)
        .arg("-b")
        .arg(tmp.path())
        .output()
        .ok()?;
    let stdout = String::from_utf8_lossy(&output.stdout);
    let mut map = HashMap::new();
    for line in stdout.lines() {
        let line = line.trim();
        if let Some((lhs, rhs)) = line.split_once('=') {
            let key = lhs.trim().to_lowercase();
            if key.starts_with("v_") {
                if let Ok(val) = rhs.split_whitespace().next()?.parse::<f64>() {
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

#[test]
fn npn_ce_cje_cjc_transient() {
    // CE inverter with explicit CJE=4pF / CJC=2pF.  Input steps 0→1V at t=0
    // through RB=1kΩ; collector through RC=1kΩ from VCC=5V.
    //
    // CJE slows base-charge build-up (τ_BE ≈ RB·CJE = 4 ns);
    // CJC adds Miller capacitance, further slowing collector fall.
    // We sample V(c) at t=6ns — mid-transition — and compare vs ngspice.
    //
    // Tolerance: 5% relative (BE truncation error vs ngspice adaptive stepper).
    // Without CJE/CJC the edge would be instantaneous and the test would fail.
    let netlist_str = "\
* BJT CE inverter — CJE/CJC transient golden\n\
.model QNPN NPN (IS=1e-15 BF=100 CJE=4e-12 VJE=0.75 MJE=0.33 CJC=2e-12 VJC=0.75 MJC=0.33 FC=0.5)\n\
VCC vcc 0 DC 5\n\
VIN in 0 PULSE(0 1 0 100p 100p 50n 100n)\n\
RB in b 1k\n\
RC vcc c 1k\n\
Q1 c b 0 0 QNPN\n\
.ic V(c)=5 V(b)=0\n\
.options uic=1\n\
.tran 20p 20n\n\
.meas tran v_14n FIND V(c) AT=14n\n\
.end\n";

    // Fairchild: parse and run Newton-Raphson transient with h=20ps, sample at t=14ns.
    // At t=14ns the collector has dropped from 5V to ~3.6V — clear mid-transition.
    let netlist = parse_spice(netlist_str).expect("parse failed");
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_bjts(&netlist.models);
    let mut opts = SimOptions::from_netlist(&netlist);
    opts.uic = true; // start from .ic values, skip DC OP
    let result = tran_nr_with_registry_opts(&netlist, 20e-12, 20e-9, &registry, &opts)
        .expect("transient failed");
    let fc_vc = result.voltage_at("c", 14e-9).expect("node c not found");

    // Shape check: V(c) must be mid-transition (caps delay the edge).
    // Without CJE/CJC the BJT would saturate instantly — V(c) < 0.5V at t=14ns.
    assert!(
        fc_vc > 1.5 && fc_vc < 4.8,
        "V(c) at 14ns = {fc_vc:.4}V — expected mid-transition [1.5, 4.8]V"
    );

    let Some(ng) = ngspice_meas_bjt(netlist_str) else {
        eprintln!("ngspice not available — shape check passed, skipping accuracy comparison");
        return;
    };
    let ng_vc = ng["v_14n"];
    let tol = f64::max(0.15, 0.05 * ng_vc.abs()); // 5% relative, 150mV absolute floor
    assert!(
        (fc_vc - ng_vc).abs() <= tol,
        "V(c) at 14ns: fairchild={fc_vc:.4e}  ngspice={ng_vc:.4e}  diff={:.2e}  tol={tol:.2e}",
        (fc_vc - ng_vc).abs()
    );
}

/// The Early effect, against ngspice, at four collector voltages.
///
/// This is the test whose absence let `VAF` be applied upside down for the
/// model's whole life (#63): every other golden leaves `VAF` infinite, where the
/// base-charge factor is exactly 1 and the sign cannot be seen.  The in-tree
/// checks could not have caught it either — they assert `IC ≈ IS·exp(VBE/VT)`
/// and `IC/IB ≈ BF`, both of which hold with the factor inverted.
#[test]
fn early_effect_output_conductance_matches_ngspice() {
    // Base held by a source so IC(VCE) is the raw output characteristic.
    let deck = |vc: f64| {
        format!(
            "* Early effect\n\
             .model qm NPN (IS=1e-16 BF=100 VAF=50)\n\
             Vb b 0 DC 0.7\n\
             Vc c 0 DC {vc}\n\
             Q1 c b 0 qm\n\
             .op\n"
        )
    };

    let mut fc_ic = Vec::new();
    let mut ng_ic = Vec::new();
    for vc in [1.0, 2.0, 4.0, 5.0] {
        let src = deck(vc);
        let fc = fairchild_op(&src);
        fc_ic.push(fc["i(vc)"]);
        let Some(ng) = ngspice_op(&src, &["i(vc)"]) else {
            eprintln!("ngspice not available — skipping comparison");
            return;
        };
        ng_ic.push(ng["i(vc)"]);
        assert_close(
            &format!("i(vc) at VC={vc}"),
            *fc_ic.last().unwrap(),
            ng["i(vc)"],
        );
    }

    // The *sign* of the output conductance, stated separately: agreeing with
    // ngspice to 0.5 % already implies it, but a reader of a failure needs to be
    // told which way round the model went.
    assert!(
        fc_ic[3].abs() > fc_ic[0].abs(),
        "|IC| must RISE with VCE — got {:.4e} at 1 V and {:.4e} at 5 V, which is \
         negative output conductance",
        fc_ic[0].abs(),
        fc_ic[3].abs()
    );
}

/// High-injection roll-off (`IKF`/`IKR`) and non-ideal junction leakage
/// (`ISE`/`NE`, `ISC`/`NC`), over seven decades of base current.
///
/// These parameters were matched by the model and discarded, with nothing on
/// stderr (#27).  Sweeping VBE from 0.4 V to 0.9 V walks through the leakage
/// floor (where beta is low because of `ISE`), the ideal region, and the knee
/// (where beta falls because of `IKF`), so a missing term shows up somewhere in
/// the sweep whichever one it is.
#[test]
fn high_injection_and_leakage_match_ngspice() {
    const CARD: &str = ".model qm NPN (IS=1e-16 BF=100 BR=2 VAF=50 \
                        IKF=1e-3 IKR=1e-2 ISE=1e-13 NE=1.6 ISC=1e-13 NC=2)\n";
    for vb in [0.4, 0.6, 0.7, 0.8, 0.9] {
        let src = format!(
            "* Gummel-Poon, all terms\n{CARD}\
             Vb b 0 DC {vb}\n\
             Vc c 0 DC 2\n\
             Q1 c b 0 qm\n\
             .op\n"
        );
        let fc = fairchild_op(&src);
        let Some(ng) = ngspice_op(&src, &["i(vb)", "i(vc)"]) else {
            eprintln!("ngspice not available — skipping comparison");
            return;
        };
        // A floor of 100 µV/100 µA is meaningless at 1 nA, so compare these on
        // relative error alone — the currents span 10⁻⁹ to 10⁻².
        for key in ["i(vb)", "i(vc)"] {
            let (a, b) = (fc[key], ng[key]);
            let rel = (a - b).abs() / b.abs().max(1e-15);
            assert!(
                rel < 5e-3,
                "{key} at VBE={vb}: fairchild={a:.6e} ngspice={b:.6e} rel={rel:.2e}"
            );
        }
        // And the point of IKF: beta must fall off at the top of the sweep.
        if vb >= 0.9 {
            let beta = fc["i(vc)"] / fc["i(vb)"];
            assert!(
                beta < 20.0,
                "beta at VBE=0.9 is {beta:.1}: the high-injection knee is not \
                 rolling it off (BF is 100)"
            );
        }
    }
}

/// `AREA` on the element line, against ngspice: two of the same transistor.
#[test]
fn bjt_area_scales_the_device_like_ngspice() {
    const SRC: &str = "* AREA\n\
                       .model qm NPN (IS=1e-16 BF=100 VAF=50)\n\
                       Vb b 0 DC 0.7\n\
                       Vc c 0 DC 2\n\
                       Vc2 c2 0 DC 2\n\
                       Q1 c b 0 qm area=2\n\
                       Q2 c2 b 0 qm\n\
                       .op\n";
    let fc = fairchild_op(SRC);
    let (i2, i1) = (fc["i(vc)"], fc["i(vc2)"]);
    assert!(
        ((i2 / i1) - 2.0).abs() < 1e-7,
        "area=2 must double IC: {i2:.9e} vs {i1:.9e}"
    );
    let Some(ng) = ngspice_op(SRC, &["i(vc)", "i(vc2)"]) else {
        eprintln!("ngspice not available — skipping comparison");
        return;
    };
    assert_close("i(vc) area=2", i2, ng["i(vc)"]);
    assert_close("i(vc2) area=1", i1, ng["i(vc2)"]);
}

// ---------------------------------------------------------------------------
// Collector-substrate junction (#97 section 3)
// ---------------------------------------------------------------------------

/// A DC solve that honours the deck's `.options`, which `fairchild_op` does not.
fn bjt_op(deck: &str) -> fairchild_core::newton::NrResult {
    let net = parse_spice(deck).expect("parse");
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&net.models);
    let opts = SimOptions::from_netlist(&net);
    fairchild_core::dc_op_nr_with_registry_opts(&net, &registry, &opts)
        .unwrap_or_else(|e| panic!("solve failed on\n{deck}\n{e:?}"))
}

/// ngspice's `.op` answer for one printed quantity on `deck`, verbatim.
///
/// The deck is shared with fairchild exactly, so the two cannot be given
/// different circuits. Only the `.control` block is appended.
fn ngspice_one(deck: &str, query: &str) -> Option<f64> {
    let ng = find_ngspice()?;
    let mut tmp = tempfile::NamedTempFile::new().ok()?;
    write!(tmp, "{deck}.control\nop\nprint {query}\n.endc\n.end\n").ok()?;
    let out = Command::new(&ng).arg("-b").arg(tmp.path()).output().ok()?;
    let text = String::from_utf8_lossy(&out.stdout);
    for line in text.lines() {
        let t = line.trim();
        if t.to_lowercase().starts_with(&query.to_lowercase()) && t.contains('=') {
            if let Ok(v) = t.split('=').nth(1)?.split_whitespace().next()?.parse() {
                return Some(v);
            }
        }
    }
    None
}

/// A reverse-biased BJT leaks **two** `gmin·V`, not one.
///
/// The headline for the substrate junction. ngspice hangs one `gmin` on the
/// collector-substrate junction, and this model had no such junction — so a
/// reverse-biased transistor read exactly half of ngspice's leakage at every
/// `gmin`. The gap was recorded rather than faked, and this closes it.
///
/// The third case is the control that identifies *which* junction the second
/// `gmin` belongs to: tying the substrate to the collector puts zero volts across
/// it, which removes exactly one `gmin·V` and nothing else. If the second `gmin`
/// came from somewhere other than the substrate junction, that case would still
/// read 2.
#[test]
fn a_reverse_biased_bjt_leaks_two_gmin_one_per_junction() {
    const MODEL: &str = ".model qn NPN (IS=1e-16 BF=100)\n";
    for g in [1e-12, 1e-9, 1e-6] {
        for (label, q_line, want_multiple) in [
            ("substrate implicit", "Q1 c b 0 qn", 2.0),
            ("substrate grounded", "Q1 c b 0 0 qn", 2.0),
            ("substrate at the collector", "Q1 c b 0 c qn", 1.0),
        ] {
            let deck = format!(
                "* reverse biased bjt\n.options gmin={g:e}\n{MODEL}\
                 VC c 0 DC 1\nVB b 0 DC 0\n{q_line}\n"
            );
            let got = bjt_op(&deck).vsrc_current("vc").expect("i(vc)").abs();
            // `gmin·1V` per junction; `IS = 1e-16` is at least four orders below
            // the smallest `gmin` here, so this reads the conductance directly.
            let want = want_multiple * g;
            let rel = (got - want).abs() / want;
            assert!(
                rel < 2e-3,
                "gmin={g:e}, {label}: collector leakage {got:.6e} A is \
                 {:.4}·gmin·V, expected {want_multiple}. One means the substrate \
                 junction is missing; three means it was added twice.",
                got / g
            );
            if let Some(ng) = ngspice_one(&deck, "i(vc)") {
                let rel = (got - ng.abs()).abs() / ng.abs();
                assert!(
                    rel < 2e-3,
                    "gmin={g:e}, {label}: fairchild {got:.6e}, ngspice {:.6e}",
                    ng.abs()
                );
            }
        }
    }
}

/// `ISS` gives the substrate junction a DC branch, and it is plain Shockley.
///
/// Not the flat reverse branch ngspice's MOS1 bulk diodes use. Measured, with
/// `gmin = 0` so nothing else conducts: at -0.05 V ngspice reads 8.553040e-16
/// where Shockley gives 8.553119e-16 and a flat `-ISS` would give 1e-15. Five
/// orders of current are covered, so a wrong exponent fails this and not only a
/// wrong prefactor.
#[test]
fn iss_gives_the_substrate_junction_a_shockley_branch() {
    const MODEL: &str = ".model qn NPN (IS=1e-16 BF=100 ISS=1e-15)\n";
    let mut compared = 0;
    for vs in [-2.0, -0.5, -0.05, 0.2, 0.4, 0.5] {
        let deck = format!(
            "* substrate dc branch\n.options gmin=0\n{MODEL}\
             VC c 0 DC 0\nVB b 0 DC 0\nVS s 0 DC {vs}\nQ1 c b 0 s qn\n"
        );
        let got = bjt_op(&deck).vsrc_current("vs").expect("i(vs)");
        let Some(ng) = ngspice_one(&deck, "i(vs)") else {
            eprintln!("ngspice not available — skipping");
            return;
        };
        let rel = (got.abs() - ng.abs()).abs() / ng.abs().max(1e-30);
        assert!(
            rel < 2e-3,
            "Vs={vs}: fairchild I(vs)={got:.6e}, ngspice {ng:.6e} (rel {rel:.2e}). \
             Exactly zero would mean ISS is not stamped; a flat -ISS in reverse \
             would read 1e-15 at -0.05 V where Shockley reads 8.553e-16."
        );
        compared += 1;
    }
    assert_eq!(compared, 6, "every bias point must have been compared");
}

/// The substrate junction is there with no `CJS` and no `ISS` on the card.
///
/// Measured: ngspice's leakage is `2·gmin·V` for a bare `IS`/`BF` card. So the
/// junction is not conditional on being given a capacitance, and neither is this.
#[test]
fn the_substrate_junction_does_not_need_cjs_to_exist() {
    for model in ["IS=1e-16 BF=100", "IS=1e-16 BF=100 CJS=2p"] {
        let deck = format!(
            "* bare card\n.options gmin=1e-9\n.model qn NPN ({model})\n\
             VC c 0 DC 1\nVB b 0 DC 0\nVS s 0 DC 0\nQ1 c b 0 s qn\n"
        );
        let got = bjt_op(&deck).vsrc_current("vs").expect("i(vs)").abs();
        assert!(
            (got / 1e-9 - 1.0).abs() < 2e-3,
            "{model}: the substrate carries {got:.6e} A, expected one gmin·V = \
             1e-9. Zero would mean the junction only exists when CJS is given."
        );
    }
}

/// Series resistance and frequency of the capacitance probe.
///
/// `omega·R·C` lands near 1.26 for a 2 pF junction, which is where the divider is
/// most sensitive to `C`. `AcResult` reports node voltages and not source
/// currents, so the capacitance is read through the divider rather than by
/// dividing a current by `omega` — and ngspice is probed the same way, on the same
/// deck, so neither simulator gets a different circuit.
const CAP_PROBE_R: f64 = 10e3;
const CAP_PROBE_F: f64 = 1e7;

/// The capacitance on `node`, from the RC divider formed with [`CAP_PROBE_R`].
///
/// `|V| = 1/sqrt(1 + (omega·R·C)²)`, so `C = sqrt(1/|V|² − 1)/(omega·R)`.
fn cap_from_divider(mag: f64) -> f64 {
    let w = 2.0 * std::f64::consts::PI * CAP_PROBE_F;
    (1.0 / (mag * mag) - 1.0).max(0.0).sqrt() / (w * CAP_PROBE_R)
}

/// fairchild's `|V(node)|` for the divider deck.
fn ac_mag(deck: &str, node: &str) -> f64 {
    let net = parse_spice(deck).expect("parse");
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&net.models);
    let opts = SimOptions::from_netlist(&net);
    let r = fairchild_core::ac_analysis_opts(&net, &[CAP_PROBE_F], Some("vs"), &registry, &opts)
        .unwrap_or_else(|e| panic!("ac failed on\n{deck}\n{e:?}"));
    r.magnitude(node, 0)
        .unwrap_or_else(|| panic!("no node '{node}' in\n{deck}"))
}

/// The probe deck: the substrate driven through [`CAP_PROBE_R`], collector at
/// `vc`, substrate DC bias `vs`.
fn cap_probe_deck(model: &str, vc: f64, vs: f64) -> String {
    format!(
        "* substrate cap probe\n.options gmin=0\n.model qn NPN ({model})\n\
         VS in 0 DC {vs} AC 1\nRS in s {CAP_PROBE_R:e}\n\
         VC c 0 DC {vc}\nVB b 0 DC 0\nQ1 c b 0 s qn\n"
    )
}

fn ngspice_ac_mag(deck: &str, query: &str, f: f64) -> Option<f64> {
    let ng = find_ngspice()?;
    let mut tmp = tempfile::NamedTempFile::new().ok()?;
    write!(
        tmp,
        "{deck}.control\nac lin 1 {f:e} {f:e}\nprint {query}\n.endc\n.end\n"
    )
    .ok()?;
    let out = Command::new(&ng).arg("-b").arg(tmp.path()).output().ok()?;
    let text = String::from_utf8_lossy(&out.stdout);
    for line in text.lines() {
        let t = line.trim();
        if t.to_lowercase().starts_with(&query.to_lowercase()) && t.contains('=') {
            if let Ok(v) = t
                .split('=')
                .nth(1)?
                .split_whitespace()
                .next()?
                .parse::<f64>()
            {
                return Some(v);
            }
        }
    }
    None
}

/// `CJS`/`VJS`/`MJS` are the collector-substrate depletion capacitance.
///
/// Two laws, and the forward one is not the one `CJE`/`CJC` use. Measured against
/// ngspice, matched to 5e-8 in reverse and 1.7e-7 forward:
///
/// ```text
/// v <= 0   CJS·(1 − v/VJS)^−MJS      the depletion law
/// v >  0   CJS·(1 + MJS·v/VJS)       a straight line from the ZERO-bias value
/// ```
///
/// The forward points past `VJS` are the load-bearing ones: at 0.8, 1.0 and 2.0 V
/// the depletion law is singular or complex, and ngspice is still exactly on the
/// straight line. `FCS` has no part in it — see
/// [`fcs_is_inert_in_ngspice_too`].
#[test]
fn cjs_follows_the_depletion_law_in_reverse_and_a_straight_line_forward() {
    const MODEL: &str = "IS=1e-16 BF=100 CJS=2p VJS=0.75 MJS=0.33";
    let mut compared = 0;
    // `vsub` is the substrate potential minus the collector's. Negative is
    // reverse for an NPN, because the substrate is p and the collector n.
    for (vc, vs) in [
        (3.0, 0.0),
        (1.0, 0.0),
        (0.5, 0.0),
        (0.0, 0.0),
        (0.0, 0.3),
        (0.0, 0.5),
        (0.0, 0.8),
        (0.0, 2.0),
    ] {
        let vsub = vs - vc;
        let deck = cap_probe_deck(MODEL, vc, vs);
        let got = cap_from_divider(ac_mag(&deck, "s"));
        let want = if vsub > 0.0 {
            2e-12 * (1.0 + 0.33 * vsub / 0.75)
        } else {
            2e-12 * (1.0 - vsub / 0.75).powf(-0.33)
        };
        let rel = (got - want).abs() / want;
        assert!(
            rel < 5e-6,
            "Vsub={vsub}: C={got:.6e} F, the law gives {want:.6e} (rel {rel:.2e}). \
             Zero would mean CJS is not stamped. Using `cj_depl` instead would \
             diverge forward, and go complex past VJS."
        );
        if let Some(m) = ngspice_ac_mag(&deck, "mag(v(s))", CAP_PROBE_F) {
            let ng = cap_from_divider(m);
            let rel = (got - ng).abs() / ng;
            assert!(
                rel < 2e-3,
                "Vsub={vsub}: fairchild {got:.6e}, ngspice {ng:.6e}"
            );
            compared += 1;
        }
    }
    if compared > 0 {
        assert_eq!(compared, 8, "every bias must have been compared to ngspice");
    }
}

/// `FCS` is accepted and correctly does nothing, because ngspice ignores it too.
///
/// Normally the worthless test shape. It earns its place the same way the MOSFET's
/// mobility-group test does: the claim is not "we accept it" but "the reference
/// ignores it, so honouring it would be a divergence", and that claim is the only
/// thing standing between `FCS` and someone re-opening it as a to-do.
///
/// Measured: the substrate capacitance at 0.5 V forward is bit-identical for `FCS`
/// of 0.1, 0.5, 0.9 and absent. `FC` is the control — it is a real parameter on
/// the *other* two junctions, so the probe would see it move if the mechanism
/// worked.
#[test]
fn fcs_is_inert_in_ngspice_too() {
    let mut values = Vec::new();
    for fcs in ["", "FCS=0.1", "FCS=0.5", "FCS=0.9"] {
        let model = format!("IS=1e-16 BF=100 CJS=2p VJS=0.75 MJS=0.33 {fcs}");
        let deck = cap_probe_deck(&model, 0.0, 0.5);
        let got = cap_from_divider(ac_mag(&deck, "s"));
        values.push(got);
        if let Some(m) = ngspice_ac_mag(&deck, "mag(v(s))", CAP_PROBE_F) {
            let ng = cap_from_divider(m);
            let rel = (got - ng).abs() / ng;
            assert!(
                rel < 2e-3,
                "FCS='{fcs}': fairchild {got:.6e}, ngspice {ng:.6e}"
            );
        }
    }
    for v in &values[1..] {
        assert!(
            (v - values[0]).abs() / values[0] < 1e-12,
            "FCS moved the substrate capacitance: {values:?}. ngspice's is \
             bit-identical across the same four cards, so honouring FCS here \
             would be a divergence from the reference, not a fix."
        );
    }
    // The control: `MJS` in the same position *does* move it, so the probe works.
    let with_mjs = cap_from_divider(ac_mag(
        &cap_probe_deck("IS=1e-16 BF=100 CJS=2p VJS=0.75 MJS=0.9", 0.0, 0.5),
        "s",
    ));
    assert!(
        (with_mjs - values[0]).abs() / values[0] > 0.1,
        "the control failed: MJS must move the capacitance this probe reads, or \
         the probe cannot tell an ignored parameter from a broken one"
    );
}

/// A BJT's own capacitances reach `.ac`, which they did not.
///
/// The BJT stamps its transient companions itself and overrode neither
/// `reactive_branches` nor `small_signal_reactances`, whose default is empty. So
/// `.ac` and `.noise` saw a transistor with no capacitance at all. Measured before
/// the fix, with `CJE = CJC = CJS = 100p` behind a 1 kOhm resistor into the base:
/// `|V(b)| = 1.000000` at 1 kHz, 1 MHz, 10 MHz and 100 MHz, where the corner is
/// 1.04 MHz.
///
/// Every transient test passed throughout, because transient takes the other path.
/// The structural gate that stops this recurring for any device is
/// `tests/circuit/reactances_reach_ac.rs`.
///
/// ngspice is the anchor, and it agrees to all six printed digits.
#[test]
fn a_bjts_own_capacitance_rolls_off_an_ac_sweep() {
    let deck = "* bjt cap in ac\n\
        .model qn NPN (IS=1e-16 BF=100 CJE=100p CJC=100p CJS=100p VJS=0.75 MJS=0.33)\n\
        V1 in 0 DC 0 AC 1\n\
        R1 in b 1k\n\
        VC c 0 DC 5\n\
        Q1 c b 0 0 qn\n";
    let net = parse_spice(deck).expect("parse");
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&net.models);
    let opts = SimOptions::from_netlist(&net);
    let freqs = [1e3, 1e6, 1e7];
    let r =
        fairchild_core::ac_analysis_opts(&net, &freqs, Some("v1"), &registry, &opts).expect("ac");
    let mut compared = 0;
    for (i, f) in freqs.iter().enumerate() {
        let mag = r.magnitude("b", i).expect("v(b)");
        // Flat at 1.0 is the bug: no capacitance anywhere in the sweep.
        if *f >= 1e6 {
            assert!(
                mag < 0.9,
                "f={f:e}: |V(b)|={mag:.6e}. With 100 pF on every junction behind \
                 1 kOhm the corner is at 1.04 MHz, so a value at 1.0 means the \
                 transistor's capacitance never reached the AC matrix."
            );
        }
        if let Some(ng) = ngspice_ac_mag(deck, "mag(v(b))", *f) {
            let rel = (mag - ng).abs() / ng;
            assert!(
                rel < 2e-3,
                "f={f:e}: fairchild |V(b)|={mag:.6e}, ngspice {ng:.6e} (rel {rel:.2e})"
            );
            compared += 1;
        }
    }
    if compared > 0 {
        assert_eq!(compared, 3, "every frequency must have been compared");
    }
}
