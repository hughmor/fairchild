//! Regression tests for fc_pn_th_ps — combined PN-junction + thermal
//! heater phase shifter.  Verify the two physics add linearly.

use fairchild_core::{DeviceRegistry, dc_op_nr_with_registry};
use fairchild_parser::parse_spice;

/// At V_pn = V_h = 0 the combined device must be identical to a passive
/// waveguide (loss-free here) — output amplitude = input amplitude,
/// modulo any propagation phase that's wrapped through.
#[test]
fn pn_th_ps_zero_bias_passes_through() {
    let netlist = "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xpnth ch0 out0 vmod 0 vh 0 fc_pn_th_ps L_um=100 dn_dv=1e-4 g_pn=1e-3 r_heater=1k p_pi=10m
Vmod vmod 0 DC 0.0
Vh vh 0 DC 0.0
.op
.end
";
    let net = parse_spice(netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let v_re = r.node_voltage("out0_re_0").unwrap();
    let v_im = r.node_voltage("out0_im_0").unwrap();
    let amp  = (v_re * v_re + v_im * v_im).sqrt();
    let p_in = 1e-3_f64;
    assert!((amp - p_in.sqrt()).abs() < 1e-6,
        "zero-bias output amp = {amp:.6}, expected {:.6}", p_in.sqrt());
}

/// Driving the heater alone produces the same phase shift as fc_thermal_ps
/// (P = V²/R, φ = π·P/P_pi).  Driving the PN alone produces the same shift
/// as fc_pn_ps.  Driving BOTH should produce the sum.  Test by comparing
/// power at the output via a photodetector — at φ = π/2 both PN and heater
/// individually rotate the field to pure imaginary; together they rotate
/// to φ = π → output = -input (real, negated).
#[test]
fn pn_th_ps_phases_add() {
    // V_h chosen so heater alone gives φ_th = π/2 → P_h = P_pi/2 → V_h = √(R·P_pi/2).
    // With R = 1k, P_pi = 10m → V_h = √(5) ≈ 2.236 V.
    let v_h_for_pi_over_2: f64 = (1000.0 * 10e-3 / 2.0_f64).sqrt();
    // V_pn chosen so PN alone gives φ_eo = π/2.
    //   φ_eo = 2π L dn/dV V_pn / λ → V_pn = λ / (4 L dn/dV).
    // L = 100 µm, dn/dV = 1e-4, λ = 1550 nm → V_pn = 1.55e-6 / (4 · 1e-4 · 1e-4) = 38.75 V.
    let v_pn_for_pi_over_2: f64 = 1.55e-6 / (4.0 * 100e-6 * 1e-4);
    let netlist = format!("\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xpnth ch0 out0 vmod 0 vh 0 fc_pn_th_ps L_um=100 dn_dv=1e-4 g_pn=1e-3 r_heater=1k p_pi=10m
Vmod vmod 0 DC {v_pn_for_pi_over_2}
Vh vh 0 DC {v_h_for_pi_over_2}
.op
.end
");
    let net = parse_spice(&netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let v_re = r.node_voltage("out0_re_0").unwrap();
    let v_im = r.node_voltage("out0_im_0").unwrap();
    // Total φ = π (mod 2π): re = -A_in, im ≈ 0.
    let p_in: f64 = 1e-3;
    let amp_in = p_in.sqrt();
    assert!((v_re - (-amp_in)).abs() < 1e-4,
        "out.re = {v_re}; expected ≈ {} for combined φ = π", -amp_in);
    assert!(v_im.abs() < 1e-4, "out.im = {v_im}; expected ≈ 0");
}

/// The two electrical interfaces must be independent — driving the heater
/// shouldn't pull current through the PN junction and vice versa.  Sanity
/// check via Vsrc current sensing.
#[test]
fn pn_th_ps_electrical_interfaces_independent() {
    let netlist = "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xpnth ch0 out0 vmod 0 vh 0 fc_pn_th_ps L_um=100 dn_dv=1e-4 g_pn=1e-3 r_heater=1k p_pi=10m
Vmod vmod 0 DC 1.0
Vh vh 0 DC 0.0
.op
.end
";
    let net = parse_spice(netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let i_pn = r.vsrc_current("vmod").unwrap().abs();
    let i_h  = r.vsrc_current("vh").unwrap().abs();
    // PN draws 1 V · 1 mS = 1 mA.
    assert!((i_pn - 1e-3).abs() < 1e-9, "I(Vmod) = {i_pn}; expected 1 mA");
    // Heater draws nothing at V_h = 0.
    assert!(i_h < 1e-9, "I(Vh) = {i_h}; expected ~0");
}
