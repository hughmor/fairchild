use std::collections::HashMap;
use std::io::Write;

use fairchild_core::{freq_decade, noise_analysis, options::SimOptions, DeviceRegistry};
use fairchild_parser::parse_spice;

const NETLIST: &str = "\
* RC thermal noise
V1 in 0 DC 1 AC 1
R1 in out 1k
C1 out 0 1u
";

fn find_ngspice() -> Option<std::path::PathBuf> {
    if std::process::Command::new("ngspice")
        .arg("--version")
        .output()
        .is_ok()
    {
        return Some("ngspice".into());
    }
    for candidate in &[
        "/opt/homebrew/bin/ngspice",
        "/usr/local/bin/ngspice",
        "/usr/bin/ngspice",
    ] {
        if std::path::Path::new(candidate).exists() {
            return Some((*candidate).into());
        }
    }
    None
}

/// Run ngspice noise analysis and return parsed PSD values keyed by label.
/// Returns None if ngspice is absent or output can't be parsed.
fn ngspice_noise() -> Option<HashMap<String, f64>> {
    let ngspice_bin = find_ngspice()?;

    let control_netlist = "\
* RC thermal noise
V1 in 0 DC 1 AC 1
R1 in out 1k
C1 out 0 1u
.control
noise V(out) V1 DEC 3 10 1000
setplot noise1
let on_0 = onoise_spectrum[0]
let on_3 = onoise_spectrum[3]
let on_6 = onoise_spectrum[6]
echo \"on_10hz = $&on_0\"
echo \"on_100hz = $&on_3\"
echo \"on_1khz = $&on_6\"
.endc
";

    let mut tmp = tempfile::NamedTempFile::new().ok()?;
    write!(tmp, "{control_netlist}").ok()?;

    let output = std::process::Command::new(&ngspice_bin)
        .arg("-b")
        .arg(tmp.path())
        .output()
        .ok()?;

    let stdout = String::from_utf8_lossy(&output.stdout);
    let mut map = HashMap::new();
    for line in stdout.lines() {
        // Match lines like: on_10hz = 1.234e-17
        if let Some((lhs, rhs)) = line.trim().split_once('=') {
            let key = lhs.trim().to_string();
            if let Ok(val) = rhs.trim().parse::<f64>() {
                map.insert(key, val);
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
fn rc_thermal_noise_vs_analytic() {
    let netlist = parse_spice(NETLIST).expect("parse");
    let registry = DeviceRegistry::new();
    let opts = SimOptions::from_netlist(&netlist);
    let freqs = freq_decade(10.0, 1000.0, 3);
    let result = noise_analysis(&netlist, &freqs, "out", "0", "v1", &registry, &opts)
        .expect("noise analysis");

    // T = 300.15 K (SPICE default), R = 1kΩ, C = 1µF, f_c = 1/(2πRC).
    let r = 1000.0_f64;
    let c = 1e-6_f64;
    let t = opts.temp_k;
    let kb = 1.380649e-23_f64;
    let s_v = 4.0 * kb * t * r;
    let f_c = 1.0 / (2.0 * std::f64::consts::PI * r * c);

    let check_analytic = |idx: usize, freq: f64| {
        let s_expected = s_v / (1.0 + (freq / f_c).powi(2));
        let s = result.onoise_psd[idx];
        let abs_floor = 1e-21_f64;
        let tol = f64::max(abs_floor, 0.005 * s_expected);
        assert!(
            (s - s_expected).abs() <= tol,
            "f={freq:.1} Hz: onoise={s:.4e} expected={s_expected:.4e} diff={:.2e} tol={tol:.2e}",
            (s - s_expected).abs()
        );
    };

    check_analytic(0, freqs[0]); // ~10 Hz
    check_analytic(3, freqs[3]); // ~100 Hz
    check_analytic(6, freqs[6]); // ~1000 Hz

    // ngspice comparison (optional — skip if absent or output unparseable).
    let Some(ng) = ngspice_noise() else {
        eprintln!("ngspice not available — skipping golden comparison");
        return;
    };

    // ngspice onoise_spectrum is V/√Hz; fairchild onoise_psd is V²/Hz.
    let ng_check = |label: &str, idx: usize| {
        let Some(&ng_vrthz) = ng.get(label) else {
            eprintln!("ngspice output missing '{label}' — skipping");
            return;
        };
        let ng_v2hz = ng_vrthz * ng_vrthz;
        let fc_val = result.onoise_psd[idx];
        let abs_floor = 1e-21_f64;
        let tol = f64::max(abs_floor, 0.01 * ng_v2hz.abs());
        assert!(
            (fc_val - ng_v2hz).abs() <= tol,
            "{label}: fairchild={fc_val:.4e} ngspice={ng_v2hz:.4e} diff={:.2e} tol={tol:.2e}",
            (fc_val - ng_v2hz).abs()
        );
    };

    ng_check("on_10hz", 0);
    ng_check("on_100hz", 3);
    ng_check("on_1khz", 6);
}
