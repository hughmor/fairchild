//! Golden comparison tests for nonlinear DC (diode circuits).
//! Runs fairchild dc_op_nr and ngspice, asserts agreement within tolerance.
//!
//! Requires ngspice on PATH. Tests are skipped (not failed) when ngspice is absent.

use std::collections::HashMap;
use std::io::Write;
use std::process::Command;

use fairchild_core::{dc_op_nr, SimError};
use fairchild_parser::{Element, parse_spice};

// Slightly relaxed compared to linear tests because ngspice includes gmin
// and temperature-model details not in our Shockley implementation.
const REL_TOL: f64 = 1e-3;     // 0.1%
const ABS_TOL_V: f64 = 1e-5;   // 10 μV floor

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
        if p.exists() { return Some(p.to_owned()); }
    }
    None
}

fn strip_control_and_end(netlist: &str) -> String {
    let mut out = String::new();
    let mut in_control = false;
    for line in netlist.lines() {
        let lc = line.trim().to_lowercase();
        if lc.starts_with(".control") { in_control = true; continue; }
        if in_control {
            if lc.starts_with(".endc") { in_control = false; }
            continue;
        }
        if lc == ".end" { continue; }
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
    if map.is_empty() { None } else { Some(map) }
}

fn ngspice_op(netlist: &str, queries: &[&str]) -> Option<HashMap<String, f64>> {
    let ngspice_bin = find_ngspice()?;
    let mut tmp = tempfile::NamedTempFile::new().ok()?;
    let stripped = strip_control_and_end(netlist);
    let print_vars = queries.join(" ");
    let control_block = format!(".control\nop\nprint {print_vars}\n.endc\n.end\n");
    write!(tmp, "{stripped}\n{control_block}").ok()?;
    let output = Command::new(&ngspice_bin).arg("-b").arg(tmp.path()).output().ok()?;
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
            let netlist_str =
                include_str!(concat!("../../../tests/golden/", $netlist_file));

            let fc = fairchild_nr_op(netlist_str).expect("fairchild NR solve failed");

            let queries: &[&str] = $queries;
            let Some(ng) = ngspice_op(netlist_str, queries) else {
                eprintln!("ngspice not available — skipping golden comparison");
                assert!(!fc.is_empty(), "fairchild produced no results");
                return;
            };

            for key in queries {
                let key_lc = key.to_lowercase();
                let fc_val = fc.get(&key_lc).copied().unwrap_or_else(|| {
                    panic!("fairchild missing '{key_lc}'; available: {fc:?}")
                });
                let ng_val = ng.get(&key_lc).copied().unwrap_or_else(|| {
                    panic!("ngspice missing '{key_lc}'; available: {ng:?}")
                });
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

diode_golden_test!(
    diode_current_source_bias,
    "diode_iv.sp",
    &["v(b)"],
    REL_TOL
);

diode_golden_test!(
    diode_series_rd,
    "diode_series_rd.sp",
    &["v(b)", "i(vdd)"],
    REL_TOL
);
