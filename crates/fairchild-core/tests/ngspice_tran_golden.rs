//! Transient golden comparison tests: run RC/RL netlists through ngspice (.meas)
//! and compare fairchild Backward Euler results at the same timepoints.
//!
//! Tolerance: 1% relative. This accounts for Backward Euler truncation error
//! (h/τ = 1µs/1ms = 0.1%) vs ngspice's higher-order adaptive integrator.

use std::collections::HashMap;
use std::io::Write;
use std::process::Command;

use fairchild_core::run_tran;
use fairchild_parser::parse_spice;

const REL_TOL: f64 = 0.01; // 1% relative tolerance

// ---------------------------------------------------------------------------
// ngspice harness for transient (.meas)
// ---------------------------------------------------------------------------

fn find_ngspice() -> Option<std::path::PathBuf> {
    if Command::new("ngspice").arg("--version").output().is_ok() {
        return Some("ngspice".into());
    }
    for candidate in &["/opt/homebrew/bin/ngspice", "/usr/local/bin/ngspice", "/usr/bin/ngspice"] {
        let p = std::path::Path::new(candidate);
        if p.exists() { return Some(p.to_owned()); }
    }
    None
}

/// Run ngspice on a netlist that already contains `.meas tran` directives.
/// Parses output lines like `v_1tau              =  6.321204e-01`.
fn ngspice_meas(netlist: &str) -> Option<HashMap<String, f64>> {
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
        // Lines look like: `v_1tau              =  6.321204e-01`
        if let Some((lhs, rhs)) = line.split_once('=') {
            let key = lhs.trim().to_lowercase();
            if key.starts_with("v_") || key.starts_with("i_") {
                if let Ok(val) = rhs.trim().parse::<f64>() {
                    map.insert(key, val);
                }
            }
        }
    }
    if map.is_empty() { None } else { Some(map) }
}

// ---------------------------------------------------------------------------
// Fairchild transient runner
// ---------------------------------------------------------------------------

fn fairchild_tran_at(netlist_str: &str, step: f64, stop: f64, node: &str, at_times: &[f64]) -> Vec<f64> {
    let netlist = parse_spice(netlist_str).expect("parse failed");
    let result = run_tran(&netlist, step, stop).expect("transient failed");
    at_times.iter().map(|&t| result.voltage_at(node, t).expect("node not found")).collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn rc_step_vs_ngspice() {
    let netlist_str = include_str!("../../../tests/golden/rc_step.sp");

    // Fairchild: 1µs fixed step (h/τ = 0.1%, BE error < 0.05%)
    let times = [1e-3_f64, 2e-3, 5e-3];
    let fc_vals = fairchild_tran_at(netlist_str, 1e-6, 5e-3, "out", &times);

    let Some(ng) = ngspice_meas(netlist_str) else {
        eprintln!("ngspice not available — validating fairchild shape only");
        // Shape check: monotone increasing, 0.5 < V(τ) < 0.75
        assert!(fc_vals[0] > 0.5 && fc_vals[0] < 0.75,  "V(out) at 1τ = {:.4}", fc_vals[0]);
        assert!(fc_vals[1] > fc_vals[0], "not monotone");
        assert!(fc_vals[2] > 0.99, "V(out) at 5τ = {:.4}", fc_vals[2]);
        return;
    };

    let ng_vals = [
        ng["v_1tau"],
        ng["v_2tau"],
        ng["v_5tau"],
    ];

    for (i, (&fc, &ng_v)) in fc_vals.iter().zip(ng_vals.iter()).enumerate() {
        let tol = REL_TOL * ng_v.abs();
        assert!(
            (fc - ng_v).abs() <= tol,
            "RC at t[{i}]: fairchild={fc:.5e}  ngspice={ng_v:.5e}  diff={:.2e}  tol={tol:.2e}",
            (fc - ng_v).abs()
        );
    }
}

#[test]
fn rl_step_vs_ngspice() {
    let netlist_str = include_str!("../../../tests/golden/rl_step.sp");

    let times = [1e-3_f64, 2e-3, 5e-3];
    let fc_vals = fairchild_tran_at(netlist_str, 1e-6, 5e-3, "out", &times);

    let Some(ng) = ngspice_meas(netlist_str) else {
        eprintln!("ngspice not available — validating fairchild shape only");
        // Shape check: monotone decreasing, V(τ) ≈ 0.37
        assert!(fc_vals[0] > 0.3 && fc_vals[0] < 0.45, "V(out) at 1τ = {:.4}", fc_vals[0]);
        assert!(fc_vals[1] < fc_vals[0], "not monotone decreasing");
        assert!(fc_vals[2] < 0.02, "V(out) at 5τ = {:.4}", fc_vals[2]);
        return;
    };

    let ng_vals = [
        ng["v_1tau"],
        ng["v_2tau"],
        ng["v_5tau"],
    ];

    for (i, (&fc, &ng_v)) in fc_vals.iter().zip(ng_vals.iter()).enumerate() {
        let tol = REL_TOL * ng_v.abs();
        assert!(
            (fc - ng_v).abs() <= tol,
            "RL at t[{i}]: fairchild={fc:.5e}  ngspice={ng_v:.5e}  diff={:.2e}  tol={tol:.2e}",
            (fc - ng_v).abs()
        );
    }
}
