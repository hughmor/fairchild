//! Characterization tests for fc_pn_ps_full / fc_pn_th_ps_full (L3 "full" PN
//! modulator: depletion + injection + TPA + static self-heating + r_series).
//! There is no external golden for this model, so these pin the as-shipped
//! behaviour at several operating points to guard the OpticalSegment +
//! PhotonicActiveModel migration. The reference values were captured from the
//! original monolithic implementation via the CLI.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

fn solve(class: &str, vbias: f64, power_mw: f64, extra: &str) -> fairchild_core::NrResult {
    let netlist = format!(
        "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW={power_mw} wavelength_nm=1550
Xpn ch0 out0 a 0 {class} L_um=500 pin_at_ref=1 {extra}
Vb a 0 DC {vbias}
.op
.end
"
    );
    let net = parse_spice(&netlist).unwrap();
    dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP converges")
}

/// Reverse bias −2 V: depletion Δn + reverse FCA loss + propagation phase.
#[test]
fn pn_ps_full_reverse_bias_depletion() {
    let r = solve("fc_pn_ps_full", -2.0, 1.0, "");
    let re = r.node_voltage("out0_re_0").unwrap();
    let im = r.node_voltage("out0_im_0").unwrap();
    assert!((re - 0.030_759_10).abs() < 1e-7, "out_re={re:.8}");
    assert!((im - 0.006_352_350).abs() < 1e-8, "out_im={im:.9}");
}

/// Mild forward bias 0.1 V: injection Δn (negative) + forward FCA loss.
#[test]
fn pn_ps_full_forward_injection() {
    let r = solve("fc_pn_ps_full", 0.1, 1.0, "");
    let re = r.node_voltage("out0_re_0").unwrap();
    let im = r.node_voltage("out0_im_0").unwrap();
    assert!((re - (-0.004_833_139)).abs() < 1e-8, "out_re={re:.9}");
    assert!((im - (-0.005_575_206)).abs() < 1e-8, "out_im={im:.9}");
}

/// Forward bias 0.7 V through R_series = 100 Ω: the junction voltage is solved
/// implicitly each iterate, so the source current reflects the V_pn ≠ V_junc
/// split. Pins that implicit solve.
#[test]
fn pn_ps_full_series_resistance_current() {
    let r = solve("fc_pn_ps_full", 0.7, 1.0, "r_series=100");
    let i = r.vsrc_current("vb").unwrap();
    assert!((i - (-1.300_568e-3)).abs() < 1e-8, "I(vb)={i:.7e}");
}

/// fc_pn_th_ps_full = full + heater. The heater (P = V_h²/R, φ_th = π·P/P_π)
/// adds a wavelength-independent rotation on top of the full PN physics.
#[test]
fn pn_th_ps_full_heater_adds_phase() {
    let base = "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xpn ch0 out0 a 0 hp 0 fc_pn_th_ps_full L_um=500 pin_at_ref=1
Vb a 0 DC -1.0
Vh hp 0 DC {vh}
.op
.end
";
    let phase = |vh: f64| {
        let net = parse_spice(&base.replace("{vh}", &format!("{vh}"))).unwrap();
        let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).unwrap();
        let re = r.node_voltage("out0_re_0").unwrap();
        let im = r.node_voltage("out0_im_0").unwrap();
        im.atan2(re)
    };
    let phase_off = phase(0.0);
    // P = V_h²/1000 = P_π/2 = 5 mW → V_h = √5 → φ_th = π/2 (observed sign −π/2).
    let phase_on = phase(5.0_f64.sqrt());
    let dphi = phase_on - phase_off;
    assert!(
        (dphi - (-std::f64::consts::FRAC_PI_2)).abs() < 1e-3,
        "heater should rotate by −π/2; got Δφ={dphi:.6}"
    );
}
