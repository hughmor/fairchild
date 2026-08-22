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
