//! Lossless transmission line (`T`) golden comparison vs ngspice.
//!
//! Validates the Branin companion against ngspice's built-in lossless line on
//! the three canonical cases: a Z0-matched line (clean delayed half-step), an
//! open far end (+1 reflection, far-end doubling at TD, near-end step at 2·TD),
//! and a shorted far end (−1 reflection). Measurements are taken on the signal
//! plateaus (well away from the ~20 ps edges) so the comparison is not
//! sensitive to fixed-step-vs-adaptive edge sampling.
//!
//! `.op` and `.ac` are here too, because a delay behaves differently in each
//! analysis and the difference used to be silent — the DC seed made the line a
//! terminator (#111) and the frequency domain dropped it entirely (#110).

use std::collections::HashMap;
use std::io::Write;
use std::process::Command;

use fairchild_core::tran_nr;
use fairchild_parser::parse_spice;

const ABS_TOL: f64 = 0.01; // 10 mV — plateaus match ngspice to <1 mV in practice

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

/// Run ngspice on a netlist containing `.meas tran v_… FIND v(node) AT=t`.
fn ngspice_meas(netlist: &str) -> Option<HashMap<String, f64>> {
    let bin = find_ngspice()?;
    let mut tmp = tempfile::NamedTempFile::new().ok()?;
    write!(tmp, "{}", netlist).ok()?;
    let out = Command::new(&bin).arg("-b").arg(tmp.path()).output().ok()?;
    let stdout = String::from_utf8_lossy(&out.stdout);
    let mut map = HashMap::new();
    for line in stdout.lines() {
        if let Some((lhs, rhs)) = line.trim().split_once('=') {
            let key = lhs.trim().to_lowercase();
            if key.starts_with("v_") {
                if let Ok(val) = rhs.trim().parse::<f64>() {
                    map.insert(key, val);
                }
            }
        }
    }
    (!map.is_empty()).then_some(map)
}

/// fairchild values at given (node, time) points. Strips `.meas`/`.end` so the
/// parser only sees the circuit + `.tran`.
fn fairchild_at(netlist: &str, points: &[(&str, f64)]) -> Vec<f64> {
    let stripped: String = netlist
        .lines()
        .filter(|l| {
            let lc = l.trim().to_lowercase();
            !lc.starts_with(".meas")
        })
        .collect::<Vec<_>>()
        .join("\n");
    let parsed = parse_spice(&stripped).expect("parse failed");
    // tran_nr builds devices (Newton-Raphson path); run_tran is linear-only.
    let result = tran_nr(&parsed, 20e-12, 4e-9).expect("transient failed");
    points
        .iter()
        .map(|&(n, t)| result.voltage_at(n, t).expect("node not found"))
        .collect()
}

fn netlist(termination: &str, meas: &str) -> String {
    format!(
        "* lossless T-line golden\n\
         Vs s 0 PULSE(0 1 0.5n 10p 10p 100n 200n)\n\
         Rs s a 50\n\
         T1 a 0 b 0 Z0=50 TD=1n\n\
         {termination}\n\
         .tran 20p 4n\n\
         {meas}\
         .end\n"
    )
}

/// Compare fairchild against ngspice (or fall back to absolute expected values
/// when ngspice is unavailable, since these lossless cases have exact answers).
fn check(net: &str, points: &[(&str, f64)], keys: &[&str], expected: &[f64]) {
    let fc = fairchild_at(net, points);
    match ngspice_meas(net) {
        Some(ng) => {
            for (i, k) in keys.iter().enumerate() {
                let ngv = ng[*k];
                assert!(
                    (fc[i] - ngv).abs() <= ABS_TOL,
                    "{k}: fairchild={:.4} ngspice={:.4} diff={:.4}",
                    fc[i],
                    ngv,
                    (fc[i] - ngv).abs()
                );
            }
        }
        None => {
            eprintln!("ngspice unavailable — checking exact lossless expectations");
            for (i, &e) in expected.iter().enumerate() {
                assert!(
                    (fc[i] - e).abs() <= ABS_TOL,
                    "point {i}: fairchild={:.4} expected={:.4}",
                    fc[i],
                    e
                );
            }
        }
    }
}

#[test]
fn matched_line_clean_delayed_step() {
    // Step launched at A (0.5 V into matched Z0), arrives at B after TD=1ns
    // (so at t=1.5ns); both ends matched ⇒ no reflections, both hold 0.5 V.
    let net = netlist(
        "RL b 0 50",
        ".meas tran v_a FIND v(a) AT=2.0n\n.meas tran v_b FIND v(b) AT=2.0n\n",
    );
    check(
        &net,
        &[("a", 2.0e-9), ("b", 2.0e-9)],
        &["v_a", "v_b"],
        &[0.5, 0.5],
    );
}

#[test]
fn open_far_end_reflection() {
    // Open B: +1 reflection. V(b) doubles to 1.0 V at TD; the reflected wave
    // reaches A at 2·TD (t=2.5ns) so V(a): 0.5 → 1.0. Steady state = source.
    let net = netlist(
        "RL b 0 1e12",
        ".meas tran v_b FIND v(b) AT=2.0n\n\
         .meas tran v_a_mid FIND v(a) AT=2.0n\n\
         .meas tran v_a_late FIND v(a) AT=3.2n\n",
    );
    check(
        &net,
        &[("b", 2.0e-9), ("a", 2.0e-9), ("a", 3.2e-9)],
        &["v_b", "v_a_mid", "v_a_late"],
        &[1.0, 0.5, 1.0],
    );
}

#[test]
fn shorted_far_end_reflection() {
    // Short B: −1 reflection. V(b)=0 always; the inverted wave reaches A at
    // 2·TD so V(a): 0.5 → 0. Steady state = shorted line, V(a)=0.
    let net = netlist(
        "RL b 0 1e-6",
        ".meas tran v_b FIND v(b) AT=2.0n\n\
         .meas tran v_a_mid FIND v(a) AT=2.0n\n\
         .meas tran v_a_late FIND v(a) AT=3.2n\n",
    );
    check(
        &net,
        &[("b", 2.0e-9), ("a", 2.0e-9), ("a", 3.2e-9)],
        &["v_b", "v_a_mid", "v_a_late"],
        &[0.0, 0.5, 0.0],
    );
}

/// The DC limit of the Branin relations: adding them gives `i1 + i2 = 0` and
/// subtracting them gives `v1 = v2`, so a lossless line is an ideal
/// through-connection at DC.
///
/// The load is 1 kΩ, deliberately not `Z0`. A line that terminated itself in
/// `Z0` would draw `1/(50+50)` and leave the far end dead, and no test using a
/// 50 Ω load could tell the two apart.
#[test]
fn dc_operating_point_conducts_through_the_line() {
    let net = "* DC bias through a lossless line\n\
               Vs s 0 DC 1\n\
               Rs s a 50\n\
               T1 a 0 b 0 Z0=50 TD=1n\n\
               Rload b 0 1k\n\
               .op\n";
    let parsed = parse_spice(net).expect("parse failed");
    let op = fairchild_core::dc_op_nr(&parsed).expect("dc op failed");
    let v_a = op.node_voltage("a").expect("node a");
    let v_b = op.node_voltage("b").expect("node b");
    let i_vs = op.vsrc_current("vs").expect("I(Vs)");

    // ngspice on the same deck: v(a) = v(b) = 0.952381, i(vs) = −9.52381e−4.
    // An absolute anchor, and also the closed form 1/(50+1000).
    let expect_v = 1000.0 / 1050.0;
    let expect_i = -1.0 / 1050.0;
    assert!(
        (v_a - expect_v).abs() < 1e-9,
        "V(a)={v_a:.9}, expected {expect_v:.9} — a Z0 terminator would give 0.5"
    );
    assert!(
        (v_b - expect_v).abs() < 1e-9,
        "V(b)={v_b:.9}, expected {expect_v:.9} — a dead far end would give 0"
    );
    assert!(
        (i_vs - expect_i).abs() < 1e-12,
        "I(Vs)={i_vs:.9}, expected {expect_i:.9} — 1/(50+50) would give −1e−2"
    );
}

/// A lossless line in `.ac`, against the closed form for a loaded line.
///
/// The anchor is analytic, not another fairchild path:
///
/// ```text
///   Zin  = Z0·(Z_L + j·Z0·tan βl) / (Z0 + j·Z_L·tan βl)
///   V(a) = Vs·Zin/(Rs + Zin)
///   V(b) = V(a) / (cos βl + j·(Z0/Z_L)·sin βl),     βl = ω·TD
/// ```
///
/// The sweep straddles the quarter-wave resonance at `1/(4·TD)` = 250 MHz,
/// where a near-open line pulls `V(a)` down to almost nothing. A model that
/// dropped the delay would report a flat 0.5 at every frequency, which is what
/// this catches. ngspice on the same deck agrees digit for digit at all five
/// points (`vm(a)` = 0.9999303, 0.7026148, 0.006283029, 0.7115001, 0.9999303).
#[test]
fn ac_sweep_matches_the_loaded_line_closed_form() {
    const Z0: f64 = 50.0;
    const RS: f64 = 50.0;
    const TD: f64 = 1e-9;
    const ZL: f64 = 1e6;

    let net = "* AC through a lossless line\n\
               Vs s 0 DC 0 AC 1\n\
               Rs s a 50\n\
               T1 a 0 b 0 Z0=50 TD=1n\n\
               Rterm b 0 1Meg\n";
    let parsed = parse_spice(net).expect("parse failed");
    let freqs: Vec<f64> = (0..5).map(|k| 1e6 + k as f64 * 125e6).collect();
    let ac = fairchild_core::ac_analysis(
        &parsed,
        &freqs,
        Some("vs"),
        &fairchild_core::device_registry::DeviceRegistry::new(),
    )
    .expect("ac failed");

    // Minimal complex helpers — a test is not worth a dependency.
    let mul = |a: (f64, f64), b: (f64, f64)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0);
    let div = |a: (f64, f64), b: (f64, f64)| {
        let d = b.0 * b.0 + b.1 * b.1;
        ((a.0 * b.0 + a.1 * b.1) / d, (a.1 * b.0 - a.0 * b.1) / d)
    };
    let mag = |a: (f64, f64)| (a.0 * a.0 + a.1 * a.1).sqrt();

    let mut saw_resonance = false;
    for (fi, &f) in freqs.iter().enumerate() {
        let bl = 2.0 * std::f64::consts::PI * f * TD;
        let t = bl.tan();
        let zin = mul((Z0, 0.0), div((ZL, Z0 * t), (Z0, ZL * t)));
        let va = div(zin, (zin.0 + RS, zin.1));
        let vb = div(va, (bl.cos(), (Z0 / ZL) * bl.sin()));

        let got_a = ac.magnitude("a", fi).expect("node a");
        let got_b = ac.magnitude("b", fi).expect("node b");
        assert!(
            (got_a - mag(va)).abs() < 1e-6,
            "f={f:.3e}: |V(a)|={got_a:.7}, closed form {:.7}",
            mag(va)
        );
        assert!(
            (got_b - mag(vb)).abs() < 1e-6,
            "f={f:.3e}: |V(b)|={got_b:.7}, closed form {:.7}",
            mag(vb)
        );
        if mag(va) < 0.05 {
            saw_resonance = true;
        }
    }
    assert!(
        saw_resonance,
        "the sweep must cross the quarter-wave resonance, or a delay-free \
         model would pass it"
    );
}

/// `.pz` has no pencil form for `exp(-s·TD)`, so it must refuse rather than
/// return the poles of the circuit with the line deleted.
#[test]
fn pz_refuses_a_circuit_containing_a_delay() {
    let net = "* pole-zero through a lossless line\n\
               Vs s 0 DC 0 AC 1\n\
               Rs s a 50\n\
               T1 a 0 b 0 Z0=50 TD=1n\n\
               Rterm b 0 1k\n\
               C1 b 0 1p\n";
    let parsed = parse_spice(net).expect("parse failed");
    let err = fairchild_core::pz::pole_zero(
        &parsed,
        &fairchild_core::device_registry::DeviceRegistry::new(),
        &fairchild_core::SimOptions::default(),
        "s",
        "0",
        "b",
        "0",
        fairchild_parser::PzDrive::Vol,
        fairchild_parser::PzWant::Both,
    )
    .expect_err("`.pz` must refuse a delay, not truncate it");
    let msg = format!("{err}");
    assert!(
        msg.contains("delay"),
        "the error must name the delay as the reason, got: {msg}"
    );
}
