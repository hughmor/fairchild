//! Native grating coupler regression tests.

use fairchild_core::{DeviceRegistry, dc_op_nr_with_registry};
use fairchild_parser::parse_spice;

/// Drive a 1 V amplitude into a 3 dB grating coupler.  Expected
/// V(out_re) = 10^(-3/20) ≈ 0.7079.
#[test]
fn grating_coupler_attenuates_by_alpha_db() {
    let netlist = "\
* grating coupler 3 dB
V_re in_re 0 DC 1.0
V_im in_im 0 DC 0.0
V_wl in_wl 0 DC 1.55e-6
X1 in_re in_im in_wl out_re out_im out_wl fc_grating_coupler alpha_dB=3.0
.op
.end
";
    let net = parse_spice(netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let v_re = r.node_voltage("out_re").unwrap();
    let v_im = r.node_voltage("out_im").unwrap();
    let v_wl = r.node_voltage("out_wl").unwrap();
    let t_expected = 10f64.powf(-3.0 / 20.0);
    assert!((v_re - t_expected).abs() < 1e-6,
        "V(out_re) = {v_re}; expected {t_expected}");
    assert!(v_im.abs() < 1e-6, "V(out_im) = {v_im}; expected 0");
    assert!((v_wl - 1.55e-6).abs() < 1e-12, "V(out_wl) should pass through");
}

/// The `alpha` keyword takes an amplitude transmission (linear scale).
/// alpha = 0.5 → IL = -20·log10(0.5) ≈ 6.02 dB → V(out_re) = 0.5.
#[test]
fn grating_coupler_alpha_keyword_is_linear_amplitude() {
    let netlist = "\
V_re in_re 0 DC 1.0
V_im in_im 0 DC 0.0
V_wl in_wl 0 DC 1.55e-6
X1 in_re in_im in_wl out_re out_im out_wl fc_grating_coupler alpha=0.5
.op
.end
";
    let net = parse_spice(netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let v_re = r.node_voltage("out_re").unwrap();
    assert!((v_re - 0.5).abs() < 1e-6, "V(out_re) = {v_re}; expected 0.5");
}

/// Grating coupler on a 2-channel WDM bundle: parser replicates per channel,
/// each channel sees the same attenuation.
#[test]
fn grating_coupler_wdm_replicates_per_channel() {
    let netlist = "\
.optical_port ch0
.optical_port ch1
.optical_port bus_in 2
.optical_port bus_out 2
.optical_port o0
.optical_port o1
Vmod0 ch0_re_0 0 DC 1.0
Vmod0i ch0_im_0 0 DC 0.0
Vmod0w ch0_wl_0 0 DC 1.55e-6
Vmod1 ch1_re_0 0 DC 0.5
Vmod1i ch1_im_0 0 DC 0.0
Vmod1w ch1_wl_0 0 DC 1.551e-6
Xmux  bus_in ch0 ch1 fc_mux
Xgc   bus_in bus_out fc_grating_coupler alpha_dB=6.02
Xdmx  bus_out o0 o1 fc_demux
.op
.end
";
    let net = parse_spice(netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    // 6.02 dB → t ≈ 0.5.  Channel 0: 1.0 → 0.5; channel 1: 0.5 → 0.25.
    let o0 = r.node_voltage("o0_re_0").unwrap();
    let o1 = r.node_voltage("o1_re_0").unwrap();
    let t = 10f64.powf(-6.02 / 20.0);
    assert!((o0 - 1.0 * t).abs() < 1e-3, "V(o0.re) = {o0}; expected ≈ 0.5");
    assert!((o1 - 0.5 * t).abs() < 1e-3, "V(o1.re) = {o1}; expected ≈ 0.25");
}
