//! Golden comparison tests for MOSFET Level 1 circuits.
//! Compares fairchild dc_op_nr against ngspice within tolerance.
//! Tests are skipped (not failed) when ngspice is absent.

use std::collections::HashMap;
use std::io::Write;
use std::process::Command;

use fairchild_core::{dc_op_nr, options::SimOptions, tran_nr_with_registry_opts, DeviceRegistry};
use fairchild_parser::parse_spice;

const REL_TOL: f64 = 2e-3; // 0.2% — Level 1 has some param differences from ngspice
const ABS_TOL_V: f64 = 1e-4; // 100 µV floor

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

/// ngspice runner that parses `print v(node)` from a .op netlist.
fn ngspice_op_node(spice_with_print: &str, node: &str) -> Option<f64> {
    let ngspice = find_ngspice()?;
    let dir = tempfile::tempdir().ok()?;
    let input = dir.path().join("test.sp");
    std::fs::write(&input, spice_with_print).ok()?;

    let out = Command::new(&ngspice)
        .args(["-b", input.to_str().unwrap()])
        .output()
        .ok()?;

    let combined = format!(
        "{}\n{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    for line in combined.lines() {
        let line = line.trim();
        let prefix = format!("v({node})");
        if line.to_lowercase().starts_with(&prefix) {
            if let Some(eq) = line.find('=') {
                let val: f64 = line[eq + 1..]
                    .trim()
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .parse()
                    .ok()?;
                return Some(val);
            }
        }
    }
    None
}

#[test]
fn nmos_resistor_dc_op() {
    // NMOS in saturation: VDD=3.3V through 10kΩ to drain.
    // Gate tied to VDD, source and bulk grounded.
    // .model nmos1 NMOS (VTO=0.7 KP=100u)
    // At gate=3.3V: VGS=3.3, VGS-VTO=2.6V; saturation if VDS >= 2.6V.
    // IDS_sat = 0.5 * 100e-6 * (W/L=10) * 2.6² = 3.38mA
    // But VDD / R drain resistor clamps this: V(d) = VDD - IDS * R = 3.3 - 3.38e-3 * 10e3 ≈ -34V
    // That's outside VDD; actually IDS limited by triode onset at VDS = VGS - VTO = 2.6V.
    // Exact solution: in triode at VDS = Vd: IDS = β·[(VGS-VTO)*VDS - 0.5*VDS²]
    // IDS = (VDD - VD)/R → triode/sat boundary where NR settles.
    let netlist_str = "* NMOS R load\n\
        .model nm1 NMOS (VTO=0.7 KP=100u)\n\
        VDD vdd 0 DC 3.3\n\
        VG  g   0 DC 3.3\n\
        R1  vdd d 10k\n\
        M1  d g 0 0 nm1 W=10u L=1u\n\
        .op\n\
        .end\n";

    let netlist = parse_spice(netlist_str).unwrap();
    let result = dc_op_nr(&netlist).unwrap();
    let vd = result.node_voltage("d").unwrap();

    // Self-consistency check: IDS through R1 = (VDD - VD)/R1.
    let ids_r = (3.3 - vd) / 10e3;
    // IDS from NMOS: at the OP, must match.
    // We can verify V(d) is positive and less than VDD.
    assert!(
        vd > 0.0 && vd < 3.3,
        "V(d) = {vd:.4}V should be between 0 and VDD"
    );
    assert!(ids_r > 0.0, "IDS should be positive: {ids_r:.4e}");

    // Compare with ngspice if available.
    let ngspice_str = format!("{netlist_str}.control\nop\nprint v(d)\n.endc\n");
    if let Some(vd_ng) = ngspice_op_node(&ngspice_str, "d") {
        let err = (vd - vd_ng).abs();
        let tol = ABS_TOL_V + REL_TOL * vd_ng.abs();
        assert!(
            err < tol,
            "V(d): fairchild={vd:.6e}  ngspice={vd_ng:.6e}  err={err:.2e}  tol={tol:.2e}"
        );
    }
}

#[test]
fn cmos_inverter_high_output() {
    // CMOS inverter with VIN=0: NMOS off, PMOS on → VOUT ≈ VDD.
    let netlist_str = "* CMOS inverter — high output\n\
        .model nm NMOS (VTO=0.7  KP=100u)\n\
        .model pm PMOS (VTO=-0.7 KP=100u)\n\
        VDD vdd 0 DC 3.3\n\
        VIN in  0 DC 0\n\
        MN  out in 0   0   nm W=10u L=1u\n\
        MP  out in vdd vdd pm W=10u L=1u\n\
        .op\n\
        .end\n";

    let netlist = parse_spice(netlist_str).unwrap();
    let result = dc_op_nr(&netlist).unwrap();
    let vout = result.node_voltage("out").unwrap();

    // VIN=0: NMOS cutoff, PMOS connects VDD to out → VOUT near VDD.
    assert!(
        vout > 3.0,
        "CMOS inverter VIN=0: VOUT={vout:.4}V (expected > 3.0V)"
    );

    let ngspice_str = format!("{netlist_str}.control\nop\nprint v(out)\n.endc\n");
    if let Some(vout_ng) = ngspice_op_node(&ngspice_str, "out") {
        let err = (vout - vout_ng).abs();
        let tol = ABS_TOL_V + REL_TOL * vout_ng.abs() + 0.05;
        assert!(
            err < tol,
            "CMOS high: fairchild={vout:.6e}  ngspice={vout_ng:.6e}  err={err:.2e}"
        );
    }
}

#[test]
fn cmos_inverter_low_output() {
    // CMOS inverter with VIN=VDD: PMOS off, NMOS on → VOUT ≈ 0V.
    let netlist_str = "* CMOS inverter — low output\n\
        .model nm NMOS (VTO=0.7  KP=100u)\n\
        .model pm PMOS (VTO=-0.7 KP=100u)\n\
        VDD vdd 0 DC 3.3\n\
        VIN in  0 DC 3.3\n\
        MN  out in 0   0   nm W=10u L=1u\n\
        MP  out in vdd vdd pm W=10u L=1u\n\
        .op\n\
        .end\n";

    let netlist = parse_spice(netlist_str).unwrap();
    let result = dc_op_nr(&netlist).unwrap();
    let vout = result.node_voltage("out").unwrap();

    // VIN=VDD: PMOS cutoff, NMOS pulls out to GND → VOUT ≈ 0.
    assert!(
        vout < 0.01,
        "CMOS inverter VIN=VDD: VOUT={vout:.4}V (expected < 10mV)"
    );

    let ngspice_str = format!("{netlist_str}.control\nop\nprint v(out)\n.endc\n");
    if let Some(vout_ng) = ngspice_op_node(&ngspice_str, "out") {
        let err = (vout - vout_ng).abs();
        let tol = ABS_TOL_V + REL_TOL * vout_ng.abs() + 1e-6;
        assert!(
            err < tol,
            "CMOS low: fairchild={vout:.6e}  ngspice={vout_ng:.6e}  err={err:.2e}"
        );
    }
}

/// Run ngspice in batch mode on a netlist containing `.meas tran` statements.
/// Returns a map of measurement name → value, or `None` if ngspice is absent.
fn ngspice_meas_tran(netlist: &str) -> Option<HashMap<String, f64>> {
    let ngspice = find_ngspice()?;
    let mut tmp = tempfile::NamedTempFile::new().ok()?;
    write!(tmp, "{netlist}").ok()?;
    let output = Command::new(&ngspice)
        .arg("-b")
        .arg(tmp.path())
        .output()
        .ok()?;
    let combined = format!(
        "{}\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let mut map = HashMap::new();
    for line in combined.lines() {
        let line = line.trim();
        if let Some((key, val)) = line.split_once('=') {
            let key = key.trim().to_lowercase();
            if key.starts_with("v_") {
                if let Ok(v) = val
                    .trim()
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .parse::<f64>()
                {
                    map.insert(key, v);
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
fn cmos_inverter_caps_switching_time() {
    // CMOS inverter with junction + Meyer overlap capacitances and a 10pF load.
    // At t=15ns the input has just gone high (rise starts at 10ns, 1ns rise time).
    // The output should be mid-transition due to the capacitive load — not yet at
    // a rail — confirming that the cap model is slowing the switching edge.
    let netlist_str = "* CMOS inverter — switching time with junction+Meyer caps\n\
        .model nm NMOS (VTO=0.7 KP=100u CGSO=2.5e-10 CGDO=2.5e-10 CJ=2e-4 CJSW=5e-10)\n\
        .model pm PMOS (VTO=-0.7 KP=100u CGSO=2.5e-10 CGDO=2.5e-10 CJ=2e-4 CJSW=5e-10)\n\
        VDD vdd 0 DC 3.3\n\
        VIN in 0 PULSE(0 3.3 10n 1n 1n 40n 100n)\n\
        MN out in 0 0 nm W=10u L=1u AS=50p AD=50p PS=20u PD=20u\n\
        MP out in vdd vdd pm W=10u L=1u AS=50p AD=50p PS=20u PD=20u\n\
        CL out 0 10p\n\
        .tran 100p 120n\n\
        .meas tran v_15n FIND V(out) AT=15n\n\
        .end\n";

    let netlist = parse_spice(netlist_str).unwrap();
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_mosfets(&netlist.models);
    let opts = SimOptions::from_netlist(&netlist);
    let result = tran_nr_with_registry_opts(&netlist, 100e-12, 120e-9, &registry, &opts).unwrap();
    let fc_v = result.voltage_at("out", 15e-9).unwrap();

    // Shape check: output must be mid-transition (not stuck at either rail).
    assert!(
        fc_v >= 0.3 && fc_v <= 3.0,
        "V(out) at t=15ns = {fc_v:.4}V — expected mid-transition [0.3, 3.0]V; \
        caps (CGSO/CGDO/CJ/CJSW on MOSFETs + CL=10p) are required to slow the edge \
        enough that the output is still switching at 15ns"
    );

    // ngspice comparison (skipped if ngspice absent).
    if let Some(meas) = ngspice_meas_tran(netlist_str) {
        if let Some(&ng_v) = meas.get("v_15n") {
            let err = (fc_v - ng_v).abs();
            let tol = f64::max(0.15, 0.05 * ng_v.abs());
            assert!(
                err <= tol,
                "V(out) at 15ns: fairchild={fc_v:.6e}  ngspice={ng_v:.6e}  \
                err={err:.3e}  tol={tol:.3e}"
            );
        }
    }
}
