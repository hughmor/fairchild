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
