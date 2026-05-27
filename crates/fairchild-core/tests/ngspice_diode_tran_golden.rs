//! Golden comparison tests for nonlinear transient (diode circuits).
//! Runs fairchild tran_nr and ngspice, asserts agreement within tolerance.
//!
//! Requires ngspice on PATH. Tests are skipped (not failed) when ngspice is absent.

use std::io::Write;
use std::process::Command;

use fairchild_core::{tran_nr, SimError};
use fairchild_parser::parse_spice;

const REL_TOL: f64 = 2e-2; // 2% — ngspice uses finer timestep internally
const ABS_TOL_V: f64 = 1e-3; // 1 mV floor

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

/// Run ngspice, extract the node voltage at a specific simulation time via .meas.
fn ngspice_tran_at(netlist: &str, node: &str, at_time: f64) -> Option<f64> {
    let ngspice_bin = find_ngspice()?;
    let mut tmp = tempfile::NamedTempFile::new().ok()?;

    // Strip any existing .control/.endc and .end, inject our own.
    let body: String = netlist
        .lines()
        .filter(|l| {
            let lc = l.trim().to_lowercase();
            !lc.starts_with(".control") && !lc.starts_with(".endc") && lc != ".end"
        })
        .chain(std::iter::once(".end"))
        .map(|l| format!("{l}\n"))
        .collect();

    // Use .meas tran to extract the value at a specific time.
    let meas_name = format!("v_{node}_at");
    let control =
        format!(".control\ntran\n.endc\n.meas tran {meas_name} FIND v({node}) AT={at_time:.3e}\n");

    write!(tmp, "{body}\n{control}").ok()?;

    let output = Command::new(&ngspice_bin)
        .arg("-b")
        .arg(tmp.path())
        .output()
        .ok()?;
    let stdout = String::from_utf8_lossy(&output.stdout);

    // Parse: "v_cap_at = 5.234567e-01 at ..."
    for line in stdout.lines() {
        let lc = line.trim().to_lowercase();
        if lc.starts_with(&meas_name.to_lowercase()) {
            if let Some(rest) = lc.strip_prefix(&meas_name.to_lowercase()) {
                if let Some(eq_rest) = rest.trim().strip_prefix('=') {
                    let val_str = eq_rest.trim().split_whitespace().next().unwrap_or("");
                    if let Ok(v) = val_str.parse::<f64>() {
                        return Some(v);
                    }
                }
            }
        }
    }
    None
}

fn fairchild_tran_at(
    netlist_str: &str,
    node: &str,
    at_time: f64,
    step: f64,
    stop: f64,
) -> Result<f64, SimError> {
    let netlist = parse_spice(netlist_str).map_err(SimError::Parse)?;
    let result = tran_nr(&netlist, step, stop)?;
    result
        .voltage_at(node, at_time)
        .ok_or_else(|| SimError::UnknownNode(node.to_owned()))
}

// ---------------------------------------------------------------------------
// Test: half-wave rectifier — capacitor voltage after first half-cycle
// ---------------------------------------------------------------------------

#[test]
fn diode_tran_rc_halfwave() {
    let netlist_str = include_str!("../../../tests/golden/diode_tran_rc.sp");
    let step = 10e-9_f64;
    let stop = 600e-6_f64;
    let at = 550e-6_f64; // near end of first positive half-cycle

    let fc =
        fairchild_tran_at(netlist_str, "cap", at, step, stop).expect("fairchild tran_nr failed");

    // Without ngspice, just sanity-check: cap should be charged (>0.4 V) but below 5V.
    assert!(
        fc > 0.4 && fc < 5.0,
        "V(cap) at t=550µs = {fc:.4e}; expected 0.4..5 V"
    );

    let Some(ng) = ngspice_tran_at(netlist_str, "cap", at) else {
        eprintln!(
            "ngspice not available — skipping golden comparison (fairchild V(cap) = {fc:.4e})"
        );
        return;
    };

    let tol = f64::max(ABS_TOL_V, REL_TOL * ng.abs());
    assert!(
        (fc - ng).abs() <= tol,
        "V(cap) at t={at:.3e}: fairchild={fc:.6e}  ngspice={ng:.6e}  diff={:.2e}  tol={tol:.2e}",
        (fc - ng).abs()
    );
}
