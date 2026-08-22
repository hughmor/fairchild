//! Regression tests for fc_splitter's alpha/r parameters.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

/// Default splitter (no params) remains the classic 3 dB lossless
/// equal-power split: both outputs at 1/√2 of the input amplitude.
#[test]
fn splitter_default_is_3db_lossless() {
    let netlist = "\
V_re in_re 0 DC 1.0
V_im in_im 0 DC 0.0
V_wl in_wl 0 DC 1.55e-6
X1 in_re in_im in_wl a_re a_im a_wl b_re b_im b_wl fc_splitter
.op
";
    let net = parse_spice(netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let half = 1.0 / 2f64.sqrt();
    assert!((r.node_voltage("a_re").unwrap() - half).abs() < 1e-6);
    assert!((r.node_voltage("b_re").unwrap() - half).abs() < 1e-6);
}

/// Asymmetric split: r = 0.9 routes 90% intensity to out_a, 10% to out_b.
/// Amplitude coefficients are √0.9 ≈ 0.9487 and √0.1 ≈ 0.3162.
#[test]
fn splitter_asymmetric_routes_by_intensity_fraction() {
    let netlist = "\
V_re in_re 0 DC 1.0
V_im in_im 0 DC 0.0
V_wl in_wl 0 DC 1.55e-6
X1 in_re in_im in_wl a_re a_im a_wl b_re b_im b_wl fc_splitter r=0.9
.op
";
    let net = parse_spice(netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let v_a = r.node_voltage("a_re").unwrap();
    let v_b = r.node_voltage("b_re").unwrap();
    assert!((v_a - 0.9f64.sqrt()).abs() < 1e-6, "V(a_re) = {v_a}");
    assert!((v_b - 0.1f64.sqrt()).abs() < 1e-6, "V(b_re) = {v_b}");
    // Total power adds to 1 (lossless when α default 1.0).
    let p_sum = v_a * v_a + v_b * v_b;
    assert!(
        (p_sum - 1.0).abs() < 1e-6,
        "P_a + P_b = {p_sum}; expected 1.0 (α=1)"
    );
}

/// Insertion loss alpha_dB = 3 → total intensity transmission ≈ 0.5.
/// Symmetric default split (r untouched by alpha_dB change) means each output
/// carries 0.25 intensity → amplitude 0.5.
#[test]
fn splitter_alpha_db_introduces_intensity_loss() {
    let netlist = "\
V_re in_re 0 DC 1.0
V_im in_im 0 DC 0.0
V_wl in_wl 0 DC 1.55e-6
X1 in_re in_im in_wl a_re a_im a_wl b_re b_im b_wl fc_splitter alpha_dB=3.0
.op
";
    let net = parse_spice(netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    // After alpha_dB=3, alpha ≈ 0.501.  r is left at 0.5 (≤ alpha), so out_a
    // gets 0.5 intensity → amplitude √0.5 ≈ 0.707; out_b gets the remaining
    // alpha − r = 0.001 intensity → amplitude ≈ 0.0316.  Verify total
    // intensity matches alpha.
    let v_a = r.node_voltage("a_re").unwrap();
    let v_b = r.node_voltage("b_re").unwrap();
    let p_total = v_a * v_a + v_b * v_b;
    let alpha = 10f64.powf(-3.0 / 10.0);
    assert!(
        (p_total - alpha).abs() < 1e-3,
        "P_a + P_b = {p_total}; expected {alpha} (α from 3 dB)"
    );
}

/// Joint asymmetric + lossy: alpha = 0.8, r = 0.6 → P_a = 0.6, P_b = 0.2.
#[test]
fn splitter_alpha_and_r_compose() {
    let netlist = "\
V_re in_re 0 DC 1.0
V_im in_im 0 DC 0.0
V_wl in_wl 0 DC 1.55e-6
X1 in_re in_im in_wl a_re a_im a_wl b_re b_im b_wl fc_splitter alpha=0.8 r=0.6
.op
";
    let net = parse_spice(netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let v_a = r.node_voltage("a_re").unwrap();
    let v_b = r.node_voltage("b_re").unwrap();
    assert!(
        (v_a * v_a - 0.6).abs() < 1e-6,
        "P_a = {}; expected 0.6",
        v_a * v_a
    );
    assert!(
        (v_b * v_b - 0.2).abs() < 1e-6,
        "P_b = {}; expected 0.2",
        v_b * v_b
    );
}
