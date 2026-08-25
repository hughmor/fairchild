//! Regression tests for fc_mzm — idealised lab-bench Mach-Zehnder modulator.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

/// V_sig = 0 → output at maximum (no insertion loss when alpha = 1) →
/// amplitude = input amplitude.
#[test]
fn mzm_at_zero_bias_passes_through() {
    let netlist = "\
.optical_port in0
.optical_port out0
Xl0 in0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xmzm in0 out0 vsig 0 fc_mzm V_pi=3.0 alpha=1.0 e_r=1k
Vsig vsig 0 DC 0.0
.op
";
    let net = parse_spice(netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let v_re = r.node_voltage("out0_re_0").unwrap();
    let p_out = v_re * v_re;
    let p_in = 1e-3_f64;
    // T(0) = α · (1 − 1/E_r) · (1 + 1)/2 + α/E_r = α (independent of E_r).
    assert!(
        (p_out - p_in).abs() < 1e-9,
        "P_out(V=0) = {p_out}; expected {p_in}"
    );
}

/// At V_sig = V_pi the output drops to α/E_r of the input intensity.
#[test]
fn mzm_at_v_pi_reaches_extinction() {
    let netlist = "\
.optical_port in0
.optical_port out0
Xl0 in0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xmzm in0 out0 vsig 0 fc_mzm V_pi=3.0 alpha=1.0 e_r=100
Vsig vsig 0 DC 3.0
.op
";
    let net = parse_spice(netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let v_re = r.node_voltage("out0_re_0").unwrap();
    let v_im = r.node_voltage("out0_im_0").unwrap();
    let p_out = v_re * v_re + v_im * v_im;
    let p_expected = 1e-3 / 100.0;
    assert!(
        (p_out - p_expected).abs() < 1e-6,
        "P_out(Vπ) = {p_out}; expected {p_expected} (1/E_r of input)"
    );
}

/// Insertion-loss `alpha_dB` scales the WHOLE transmission curve.  At V_sig
/// = 0 with `alpha_dB = 3` (intensity transmission ≈ 0.501), expect output
/// power = 0.501 mW.
#[test]
fn mzm_alpha_db_scales_max_transmission() {
    let netlist = "\
.optical_port in0
.optical_port out0
Xl0 in0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xmzm in0 out0 vsig 0 fc_mzm V_pi=3.0 alpha_dB=3.0 e_r=10k
Vsig vsig 0 DC 0.0
.op
";
    let net = parse_spice(netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let v_re = r.node_voltage("out0_re_0").unwrap();
    let p_out = v_re * v_re;
    let expected = 1e-3 * 10f64.powf(-3.0 / 10.0);
    assert!(
        (p_out - expected).abs() < 1e-6,
        "P_out = {p_out}; expected {expected} (3 dB IL on 1 mW input)"
    );
}

/// `e_r_dB = 20` ⇒ E_r = 100.  At V_sig = V_pi the floor is 1/100 of α·P_in.
#[test]
fn mzm_e_r_db_keyword_works() {
    let netlist = "\
.optical_port in0
.optical_port out0
Xl0 in0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xmzm in0 out0 vsig 0 fc_mzm V_pi=3.0 alpha=1.0 e_r_dB=20
Vsig vsig 0 DC 3.0
.op
";
    let net = parse_spice(netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let v_re = r.node_voltage("out0_re_0").unwrap();
    let v_im = r.node_voltage("out0_im_0").unwrap();
    let p_out = v_re * v_re + v_im * v_im;
    let p_expected = 1e-3 / 100.0;
    assert!(
        (p_out - p_expected).abs() < 1e-6,
        "P_out at Vπ with 20 dB E_r = {p_out}; expected {p_expected}"
    );
}

/// `f_c` is accepted as a forward-compat keyword and is **not wired into the
/// response**, and this pins both halves of that.
///
/// It used to assert only that a deck naming `f_c` still solved — a test that
/// passes whether `f_c` is implemented, dropped, or deleted from the model. The
/// no-op is the actual claim (`mzm.rs`: "accepted but not yet wired into the
/// response"), so the no-op is what gets asserted. When someone implements the
/// EO cutoff, this fails and points at the doc comment that has to change with
/// it.
///
/// Two settings three decades apart: if `f_c` ever reaches the response, a DC
/// point may still not move, so the assertion is deliberately about *this*
/// analysis and says so in its name.
#[test]
fn mzm_f_c_is_accepted_and_does_not_move_the_dc_point() {
    let deck = |extra: &str| {
        format!(
            ".optical_port in0\n.optical_port out0\n\
             Xl0 in0 fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
             Xmzm in0 out0 vsig 0 fc_mzm V_pi=3.0 {extra}\n\
             Vsig vsig 0 DC 0.0\n.op\n"
        )
    };
    let solve = |src: String| {
        let net = parse_spice(&src).unwrap();
        let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
        r.node_voltage("out0_re_0").expect("out0_re_0")
    };
    let none = solve(deck(""));
    let slow = solve(deck("f_c=1G"));
    let fast = solve(deck("f_c=1000G"));
    assert_eq!(
        none, slow,
        "f_c is not wired into the response; it moved V(out0_re_0)"
    );
    assert_eq!(
        slow, fast,
        "f_c=1G and f_c=1000G must agree while f_c is a no-op"
    );

    // …and the keyword is genuinely recognised rather than ignored as a typo.
    use fairchild_core::device::Device;
    let mut mzm = fairchild_core::models::NativeMzm::default();
    assert!(
        mzm.set_real_param("f_c", 1e10),
        "f_c must be recognised, not fall through as unknown"
    );
}
