//! Native WDM MUX / DEMUX regression tests.
//!
//! Two lasers at different powers and wavelengths share one waveguide via
//! `fc_mux` (combine) and `fc_demux` (split).  Each laser is on its own
//! single-channel bundle; the bus between MUX and DEMUX is a 2-channel
//! bundle.  Verifies that:
//!
//!   - the parser doesn't replicate fc_mux / fc_demux per channel,
//!   - the bus-side bundle width matches the channel count of the MUX inputs,
//!   - channels are independent (asymmetric input powers produce asymmetric
//!     outputs),
//!   - the parser's `bundle ports must have matching channel counts` check
//!     is bypassed for the bridging devices.

use fairchild_core::{DeviceRegistry, dc_op_nr_with_registry};
use fairchild_parser::parse_spice;

fn wdm_via_mux(p1_mw: f64, p2_mw: f64) -> String {
    format!("\
* WDM via fc_mux / fc_demux
.optical_port ch0
.optical_port ch1
.optical_port wdm_bus 2
.optical_port wg_out 2
.optical_port d0
.optical_port d1

Xl1 ch0 fc_cw_laser power_mW={p1_mw} wavelength_nm=1549.9
Xl2 ch1 fc_cw_laser power_mW={p2_mw} wavelength_nm=1550.1

Xmux  wdm_bus ch0 ch1 fc_mux
Xwg   wdm_bus wg_out fc_waveguide L_um=100 n_g=4.2 alpha_dB_cm=2.0 wavelength_nm=1550
Xdemux wg_out d0 d1 fc_demux

Xpd1 d0 v_pd1 0 fc_photodetector responsivity=0.8
Xpd2 d1 v_pd2 0 fc_photodetector responsivity=0.8
Vbias bias 0 DC 1.0
Rl1 v_pd1 bias 1k
Rl2 v_pd2 bias 1k
.op
.end
")
}

/// Equal power on both channels → equal PD outputs.
#[test]
fn wdm_mux_demux_symmetric_at_equal_power() {
    let net = parse_spice(&wdm_via_mux(1.0, 1.0)).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let v1 = r.node_voltage("v_pd1").unwrap();
    let v2 = r.node_voltage("v_pd2").unwrap();
    assert!((v1 - v2).abs() < 1e-6,
        "equal-power channels should produce equal PD voltages: v1={v1:.6} v2={v2:.6}");
    // Expected: V_bias + P · R_load · responsivity · waveguide_loss
    //   = 1 + 1e-3 · 0.8 · 1000 · exp(-α·L/2) ≈ 1 + 0.8 · 0.9988 ≈ 1.799 V.
    assert!((v1 - 1.798).abs() < 0.01,
        "v_pd1 = {v1:.4}, expected ≈ 1.798 V (1 V bias + 0.8 V photocurrent through 1k after waveguide loss)");
}

/// Asymmetric powers (2 mW / 0.5 mW) → linearly asymmetric outputs.  This
/// pins channel independence: any cross-talk through the MUX/DEMUX/bus path
/// would mix the channels and break the linear scaling.
#[test]
fn wdm_mux_demux_channels_are_independent() {
    let net = parse_spice(&wdm_via_mux(2.0, 0.5)).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let v1 = r.node_voltage("v_pd1").unwrap();
    let v2 = r.node_voltage("v_pd2").unwrap();
    // Channel 0 carries 2 mW → V_pd1 ≈ 1 + 2 · 0.8 · 0.9988 ≈ 2.598 V.
    assert!((v1 - 2.598).abs() < 0.01,
        "v_pd1 = {v1:.4}, expected ≈ 2.598 V at 2 mW input");
    // Channel 1 carries 0.5 mW → V_pd2 ≈ 1 + 0.5 · 0.8 · 0.9988 ≈ 1.399 V.
    assert!((v2 - 1.399).abs() < 0.01,
        "v_pd2 = {v2:.4}, expected ≈ 1.399 V at 0.5 mW input");
    // The photocurrents should differ by exactly the laser-power ratio.
    let i1 = (v1 - 1.0) / 1000.0;
    let i2 = (v2 - 1.0) / 1000.0;
    let ratio = i1 / i2;
    assert!((ratio - 4.0).abs() < 0.01,
        "photocurrent ratio {ratio:.4}, expected ≈ 4.0 (2 mW / 0.5 mW)");
}

/// 4-channel MUX/DEMUX — verify the inferred bundle width scales beyond N=2.
#[test]
fn wdm_mux_demux_n4_routes_four_channels() {
    let netlist_str = "\
* 4-channel WDM
.optical_port c0
.optical_port c1
.optical_port c2
.optical_port c3
.optical_port bus 4
.optical_port out_bus 4
.optical_port d0
.optical_port d1
.optical_port d2
.optical_port d3

Xl0 c0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xl1 c1 fc_cw_laser power_mW=2.0 wavelength_nm=1551
Xl2 c2 fc_cw_laser power_mW=3.0 wavelength_nm=1552
Xl3 c3 fc_cw_laser power_mW=4.0 wavelength_nm=1553

Xmux   bus c0 c1 c2 c3 fc_mux
Xwg    bus out_bus fc_waveguide L_um=100 n_g=4.2 alpha_dB_cm=2.0 wavelength_nm=1550
Xdemux out_bus d0 d1 d2 d3 fc_demux

Xpd0 d0 v0 0 fc_photodetector responsivity=0.8
Xpd1 d1 v1 0 fc_photodetector responsivity=0.8
Xpd2 d2 v2 0 fc_photodetector responsivity=0.8
Xpd3 d3 v3 0 fc_photodetector responsivity=0.8
Vbias bias 0 DC 1.0
R0 v0 bias 1k
R1 v1 bias 1k
R2 v2 bias 1k
R3 v3 bias 1k
.op
.end
";
    let net = parse_spice(netlist_str).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let v0 = r.node_voltage("v0").unwrap();
    let v1 = r.node_voltage("v1").unwrap();
    let v2 = r.node_voltage("v2").unwrap();
    let v3 = r.node_voltage("v3").unwrap();
    // Each PD voltage = 1 + P_mW · 0.8 · 1k · 1e-3 · waveguide_loss(≈0.9988).
    let expected = |p_mw: f64| 1.0 + p_mw * 0.8 * 0.9988;
    for (probe, p_mw, v) in [("v0", 1.0, v0), ("v1", 2.0, v1), ("v2", 3.0, v2), ("v3", 4.0, v3)] {
        let exp = expected(p_mw);
        assert!((v - exp).abs() < 0.01,
            "{probe} = {v:.4}, expected ≈ {exp:.4} at {p_mw} mW input");
    }
}
