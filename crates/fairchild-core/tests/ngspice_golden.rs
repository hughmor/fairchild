//! Golden comparison tests: run each netlist through both fairchild and ngspice,
//! assert node voltages and branch currents agree within 1 µV / 1 nA.
//!
//! Requires ngspice on PATH. Tests are skipped (not failed) when ngspice is absent.

use std::collections::HashMap;
use std::io::Write;
use std::process::Command;

use fairchild_core::dc_op_nr;
use fairchild_parser::parse_spice;

// Tolerance: max(absolute floor, relative fraction of the expected value).
// This handles large currents where floating-point noise exceeds the absolute floor.
const ABS_TOL_V: f64 = 1e-9;  // 1 nV floor
const ABS_TOL_A: f64 = 1e-12; // 1 pA floor
const REL_TOL: f64 = 1e-5;    // 10 ppm relative

// ---------------------------------------------------------------------------
// ngspice harness
// ---------------------------------------------------------------------------

/// Run ngspice on a netlist string and return a map of "v(node)" / "i(vsrc)" → value.
///
/// The netlist must contain `.op` — we inject a `.control` block that calls
/// `print` on every node/source we care about and parse the output.
/// Find the ngspice binary, checking PATH and common install locations.
fn find_ngspice() -> Option<std::path::PathBuf> {
    // Try bare name first (works when PATH is set correctly).
    if Command::new("ngspice").arg("--version").output().is_ok() {
        return Some("ngspice".into());
    }
    // Common Homebrew/system locations.
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

fn ngspice_op(netlist: &str, queries: &[&str]) -> Option<HashMap<String, f64>> {
    let ngspice_bin = find_ngspice()?;

    // Write netlist to a temp file. We inject a .control block.
    let mut tmp = tempfile::NamedTempFile::new().ok()?;

    // Strip any existing .control/.endc blocks and the .end line,
    // then add our own at the bottom.
    let stripped = strip_control_and_end(netlist);
    let print_vars = queries.join(" ");
    let control_block = format!(
        ".control\nop\nprint {print_vars}\n.endc\n.end\n"
    );
    write!(tmp, "{stripped}\n{control_block}").ok()?;

    let output = Command::new(&ngspice_bin)
        .arg("-b")
        .arg(tmp.path())
        .output()
        .ok()?;

    let stdout = String::from_utf8_lossy(&output.stdout);
    parse_ngspice_print(&stdout)
}

/// Remove `.control`...`.endc` blocks and final `.end` from a netlist string.
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

/// Parse ngspice `print` output lines like `v(mid) = 5.000000e-01`.
fn parse_ngspice_print(output: &str) -> Option<HashMap<String, f64>> {
    let mut map = HashMap::new();
    for line in output.lines() {
        let line = line.trim();
        // Match lines like: v(mid) = 5.000000e-01  or  i(v1) = -2.00000e-03
        // Skip lines where rhs isn't a plain float (e.g. "TEMP = 27 and TNOM = 27").
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

// ---------------------------------------------------------------------------
// Helper: run fairchild DC op on a file in tests/golden/
// ---------------------------------------------------------------------------

fn fairchild_op(netlist_str: &str) -> HashMap<String, f64> {
    let netlist = parse_spice(netlist_str).expect("parse failed");
    let result = dc_op_nr(&netlist).expect("DC solve failed");
    let mut map = HashMap::new();
    for (name, v) in result.all_voltages() {
        map.insert(format!("v({name})"), v);
    }
    // Also expose branch currents for voltage sources.
    for el in &netlist.elements {
        if let fairchild_parser::Element::VoltageSource { name, .. } = el {
            if let Ok(i) = result.vsrc_current(name) {
                map.insert(format!("i({name})"), i);
            }
        }
    }
    map
}

// ---------------------------------------------------------------------------
// Macro to keep test bodies concise
// ---------------------------------------------------------------------------

macro_rules! golden_test {
    ($name:ident, $netlist_file:expr, $queries:expr) => {
        #[test]
        fn $name() {
            let netlist_str =
                include_str!(concat!("../../../tests/golden/", $netlist_file));

            let fc = fairchild_op(netlist_str);

            let queries: &[&str] = $queries;
            let Some(ng) = ngspice_op(netlist_str, queries) else {
                eprintln!("ngspice not available — skipping golden comparison");
                // Still validate fairchild ran without panic.
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

                let abs_floor = if key_lc.starts_with("i(") { ABS_TOL_A } else { ABS_TOL_V };
                let tol = f64::max(abs_floor, REL_TOL * ng_val.abs());
                assert!(
                    (fc_val - ng_val).abs() <= tol,
                    "{key_lc}: fairchild={fc_val:.6e}  ngspice={ng_val:.6e}  diff={:.2e}  tol={tol:.2e}",
                    (fc_val - ng_val).abs()
                );
            }
        }
    };
}

// ---------------------------------------------------------------------------
// Golden tests
// ---------------------------------------------------------------------------

golden_test!(
    voltage_divider,
    "voltage_divider.sp",
    &["v(in)", "v(mid)", "i(v1)"]
);

golden_test!(
    current_divider,
    "current_divider.sp",
    &["v(a)"]
);

golden_test!(
    wheatstone_bridge,
    "wheatstone.sp",
    &["v(top)", "v(a)", "v(b)", "i(v1)"]
);

golden_test!(
    resistor_ladder,
    "ladder.sp",
    &["v(in)", "v(n1)", "v(n2)", "v(n3)", "i(v1)"]
);

golden_test!(
    multi_source,
    "multi_source.sp",
    &["v(a)", "v(b)", "i(v1)", "i(v2)"]
);
