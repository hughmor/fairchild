use std::collections::HashMap;
use std::io::Write;
use std::process::Command;

use fairchild_core::{tran_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

const REL_TOL: f64 = 0.02;
const ABS_TOL_V: f64 = 10e-3;

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

fn ngspice_meas_tran(netlist: &str) -> Option<HashMap<String, f64>> {
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
fn transformer_1_1_step_response_vs_ngspice() {
    let netlist_str = "\
* 1:1 transformer step response k=0.9\n\
V1 prim 0 PULSE(0 1 0 1n 1n 200u 500u)\n\
R1 prim n1 100\n\
L1 n1 0 1m\n\
L2 n2 0 1m\n\
R2 n2 0 100\n\
K1 l1 l2 0.9\n\
.tran 10n 20u\n\
.meas tran v_n2_5u FIND V(n2) AT=5u\n\
.meas tran v_n1_5u FIND V(n1) AT=5u\n\
.end\n";

    let netlist = parse_spice(netlist_str).expect("parse");
    let registry = DeviceRegistry::new();
    let result = tran_nr_with_registry(&netlist, 10e-9, 20e-6, &registry).expect("sim");
    let v_n1 = result.voltage_at("n1", 5e-6).expect("n1");
    let v_n2 = result.voltage_at("n2", 5e-6).expect("n2");

    assert!(
        (0.2..=0.9).contains(&v_n1),
        "V(n1) at 5µs = {v_n1:.4}V — expected mid-transient [0.2, 0.9]V"
    );
    assert!(
        (0.1..=0.9).contains(&v_n2),
        "V(n2) at 5µs = {v_n2:.4}V — expected coupled [0.1, 0.9]V"
    );
    let ratio = v_n2 / v_n1;
    assert!(
        (0.3..=1.1).contains(&ratio),
        "V(n2)/V(n1) = {ratio:.4} — expected coupling ratio [0.3, 1.1]"
    );

    let Some(ng) = ngspice_meas_tran(netlist_str) else {
        eprintln!("ngspice not available — shape checks passed, skipping accuracy comparison");
        return;
    };
    let ng_v_n1 = ng["v_n1_5u"];
    let ng_v_n2 = ng["v_n2_5u"];

    let tol_n1 = f64::max(ABS_TOL_V, REL_TOL * ng_v_n1.abs());
    assert!(
        (v_n1 - ng_v_n1).abs() <= tol_n1,
        "V(n1) at 5µs: fairchild={v_n1:.6e}  ngspice={ng_v_n1:.6e}  diff={:.2e}  tol={tol_n1:.2e}",
        (v_n1 - ng_v_n1).abs()
    );

    let tol_n2 = f64::max(ABS_TOL_V, REL_TOL * ng_v_n2.abs());
    assert!(
        (v_n2 - ng_v_n2).abs() <= tol_n2,
        "V(n2) at 5µs: fairchild={v_n2:.6e}  ngspice={ng_v_n2:.6e}  diff={:.2e}  tol={tol_n2:.2e}",
        (v_n2 - ng_v_n2).abs()
    );
}
