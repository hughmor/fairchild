//! Voltage- and current-controlled switch (`S` / `W`) golden comparison vs ngspice.
//!
//! The switching law is a hard step with an optional hysteresis band, so the
//! things worth pinning are the *boundaries*: which side of `VT ± VH` flips the
//! switch, and that `RON`/`ROFF` are the resistances either side. Those were
//! read off ngspice 46 directly (see `models/switch.rs`), and this test keeps
//! them honest against whatever ngspice is installed.
//!
//! The transient case is a clocked sample-and-hold — what switches are actually
//! for. Hysteresis-band behaviour is pinned by the unit tests in
//! `models/switch.rs`, which do not need an external simulator.

use std::collections::HashMap;
use std::io::Write;
use std::process::Command;

use fairchild_core::{dc_op_nr, tran_nr};
use fairchild_parser::parse_spice;

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

/// Run ngspice and collect every `name = value` line whose name starts `m_`.
fn ngspice_meas(netlist: &str) -> Option<HashMap<String, f64>> {
    let bin = find_ngspice()?;
    let mut tmp = tempfile::NamedTempFile::new().ok()?;
    write!(tmp, "{netlist}").ok()?;
    let out = Command::new(&bin).arg("-b").arg(tmp.path()).output().ok()?;
    let stdout = String::from_utf8_lossy(&out.stdout);
    let mut map = HashMap::new();
    for line in stdout.lines() {
        if let Some((lhs, rhs)) = line.trim().split_once('=') {
            let key = lhs.trim().to_lowercase();
            if key.starts_with("m_") {
                if let Ok(val) = rhs.trim().parse::<f64>() {
                    map.insert(key, val);
                }
            }
        }
    }
    (!map.is_empty()).then_some(map)
}

fn strip_meas(netlist: &str) -> String {
    netlist
        .lines()
        .filter(|l| !l.trim().to_lowercase().starts_with(".meas"))
        .collect::<Vec<_>>()
        .join("\n")
}

/// A divider whose upper leg is the switch, so `v(b)` reads the state directly:
/// `RON = 10` gives 0.990099 V and `ROFF = 1e6` gives 0.999 mV across `RL = 1k`.
fn s_divider(vctrl: f64, vt: f64, vh: f64, keyword: &str) -> String {
    format!(
        "* S switch divider\n\
         .model swmod SW (VT={vt} VH={vh} RON=10 ROFF=1e6)\n\
         Vc c 0 DC {vctrl}\n\
         V1 a 0 DC 1\n\
         S1 a b c 0 swmod {keyword}\n\
         RL b 0 1k\n\
         .op\n\
         .meas dc m_vb FIND v(b)\n\
         .end\n"
    )
}

#[test]
fn voltage_switch_thresholds_match_ngspice() {
    // 0.9/1.0 are below-and-at VT (both OFF: the law is strictly greater), 1.1
    // is above. This is the boundary ngspice was measured on.
    for vctrl in [0.9, 1.0, 1.1, 5.0] {
        let deck = s_divider(vctrl, 1.0, 0.0, "");
        let got = dc_op_nr(&parse_spice(&strip_meas(&deck)).unwrap())
            .expect("DC OP failed")
            .node_voltage("b")
            .unwrap();
        let want_on = vctrl > 1.0;
        let want = if want_on {
            1000.0 / 1010.0
        } else {
            1000.0 / 1_001_000.0
        };
        assert!(
            (got - want).abs() < 1e-9,
            "Vctrl={vctrl}: got v(b)={got}, want {want} ({})",
            if want_on { "ON" } else { "OFF" }
        );

        let Some(ng) = ngspice_meas(&deck) else {
            eprintln!("ngspice not available; skipping the comparison half");
            continue;
        };
        let ref_v = ng["m_vb"];
        assert!(
            (got - ref_v).abs() < 1e-6,
            "Vctrl={vctrl}: fairchild {got}, ngspice {ref_v}"
        );
    }
}

/// Inside the band the `ON`/`OFF` keyword decides, and it must reach the device.
#[test]
fn the_on_off_keyword_selects_the_state_inside_the_band() {
    for (keyword, want_on) in [("ON", true), ("OFF", false), ("", false)] {
        // VT=1, VH=0.5 → the band is [0.5, 1.5]; 1.0 V sits inside it.
        let deck = s_divider(1.0, 1.0, 0.5, keyword);
        let got = dc_op_nr(&parse_spice(&strip_meas(&deck)).unwrap())
            .expect("DC OP failed")
            .node_voltage("b")
            .unwrap();
        let want = if want_on {
            1000.0 / 1010.0
        } else {
            1000.0 / 1_001_000.0
        };
        assert!(
            (got - want).abs() < 1e-9,
            "keyword '{keyword}': got {got}, want {want}"
        );
    }
}

#[test]
fn current_switch_thresholds_match_ngspice() {
    // `Ic` pushes current through the sense source; IT = 1 mA.
    for i_ctrl in [0.5e-3, 2e-3] {
        let deck = format!(
            "* W switch divider\n\
             .model cswmod CSW (IT=1m IH=0 RON=10 ROFF=1e6)\n\
             Vsense s 0 DC 0\n\
             Ic 0 s DC {i_ctrl}\n\
             V1 a 0 DC 1\n\
             W1 a b Vsense cswmod\n\
             RL b 0 1k\n\
             .op\n\
             .meas dc m_vb FIND v(b)\n\
             .end\n"
        );
        let got = dc_op_nr(&parse_spice(&strip_meas(&deck)).unwrap())
            .expect("DC OP failed")
            .node_voltage("b")
            .unwrap();
        let want_on = i_ctrl > 1e-3;
        let want = if want_on {
            1000.0 / 1010.0
        } else {
            1000.0 / 1_001_000.0
        };
        assert!(
            (got - want).abs() < 1e-9,
            "I={i_ctrl}: got v(b)={got}, want {want}"
        );

        let Some(ng) = ngspice_meas(&deck) else {
            continue;
        };
        assert!(
            (got - ng["m_vb"]).abs() < 1e-6,
            "I={i_ctrl}: fairchild {got}, ngspice {}",
            ng["m_vb"]
        );
    }
}

/// A clocked sample-and-hold — the switch's actual job.
///
/// A ramp is sampled onto a cap while the clock is high and held while it is
/// low, so every hold plateau is a different voltage and a state error shows up
/// as a wrong plateau rather than a wrong slope.
///
/// Deliberately *not* a self-referencing switch (control node == switched
/// node). That circuit chatters under any hard-switch model — ngspice 46 gives
/// up on its DC operating point ("gmin stepping failed", "source stepping
/// failed") and then creeps to the band edge and stalls rather than
/// oscillating. It is not a fair reference, and it is not what switches are
/// used for.
#[test]
fn a_clocked_sample_and_hold_matches_ngspice() {
    let deck = "* switched-capacitor sample and hold\n\
         .model swmod SW (VT=2.5 VH=0 RON=10 ROFF=1e9)\n\
         Vin in 0 PWL(0 0 100u 5)\n\
         Vclk c 0 PULSE(0 5 0 1n 1n 10u 20u)\n\
         S1 in out c 0 swmod OFF\n\
         C1 out 0 1n\n\
         .tran 0.2u 100u\n\
         .meas tran m_h1 FIND v(out) AT=18u\n\
         .meas tran m_h2 FIND v(out) AT=38u\n\
         .meas tran m_h3 FIND v(out) AT=58u\n\
         .meas tran m_t1 FIND v(out) AT=25u\n\
         .end\n";
    let parsed = parse_spice(&strip_meas(deck)).unwrap();
    let r = tran_nr(&parsed, 0.2e-6, 100e-6).expect("transient failed");

    // Hold plateaus: the input at the falling edge (10u, 30u, 50u) is
    // 0.5/1.5/2.5 V on a 5 V / 100 µs ramp. Tolerance covers one step of ramp.
    for (t_hold, want) in [(18e-6, 0.5), (38e-6, 1.5), (58e-6, 2.5)] {
        let got = r.voltage_at("out", t_hold).unwrap();
        assert!(
            (got - want).abs() < 0.03,
            "hold at {t_hold:e}: got {got}, want ~{want}"
        );
    }
    // And it tracks while the clock is high.
    let tracking = r.voltage_at("out", 25e-6).unwrap();
    assert!(
        (tracking - 1.25).abs() < 0.03,
        "tracking at 25 µs: got {tracking}, want ~1.25"
    );

    let Some(ng) = ngspice_meas(deck) else {
        eprintln!("ngspice not available; skipping the comparison half");
        return;
    };
    for (key, t) in [
        ("m_h1", 18e-6),
        ("m_h2", 38e-6),
        ("m_h3", 58e-6),
        ("m_t1", 25e-6),
    ] {
        let got = r.voltage_at("out", t).unwrap();
        assert!(
            (got - ng[key]).abs() < 0.03,
            "{key}: fairchild {got}, ngspice {}",
            ng[key]
        );
    }
}
