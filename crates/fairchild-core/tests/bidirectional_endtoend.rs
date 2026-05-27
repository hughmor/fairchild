//! End-to-end bidirectional propagation tests using the C.2-refactored
//! devices that are bidir-aware so far: fc_cw_laser, fc_waveguide, and
//! fc_photodetector.  These are the minimal trio that can drive a forward
//! signal chain — and they must also produce the same answer as the
//! unidirectional baseline because no backward light enters the chain.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

/// Laser → waveguide → PD chain.  Compare the PD voltage with bidir off
/// vs. on — the answer should be identical: the laser drives fw only,
/// the waveguide propagates fw, the PD absorbs fw.  No bw light exists.
#[test]
fn bidir_chain_matches_unidir_for_forward_only_drive() {
    let unidir = "\
.optical_port src
.optical_port wg_out
Xl src fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xwg src wg_out fc_waveguide L_um=100 n_g=4.2 alpha_dB_cm=2.0
Xpd wg_out pd_a 0 fc_photodetector responsivity=0.8 i_dark_a=1e-12 r_shunt=1Meg
Vb bias 0 DC 1.0
Rload pd_a bias 1k
.op
.end
";
    let bidir = "\
.options enable_bidirectional=1
.optical_port src
.optical_port wg_out
Xl src fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xwg src wg_out fc_waveguide L_um=100 n_g=4.2 alpha_dB_cm=2.0
Xpd wg_out pd_a 0 fc_photodetector responsivity=0.8 i_dark_a=1e-12 r_shunt=1Meg
Vb bias 0 DC 1.0
Rload pd_a bias 1k
.op
.end
";
    let r_uni = dc_op_nr_with_registry(&parse_spice(unidir).unwrap(), &DeviceRegistry::new())
        .expect("DC OP (unidir)");
    let r_bi = dc_op_nr_with_registry(&parse_spice(bidir).unwrap(), &DeviceRegistry::new())
        .expect("DC OP (bidir)");
    let v_uni = r_uni.node_voltage("pd_a").unwrap();
    let v_bi = r_bi.node_voltage("pd_a").unwrap();
    assert!(
        (v_uni - v_bi).abs() < 1e-6,
        "PD anode should match between unidir and bidir for fw-only drive: \
         uni={v_uni:.6} bi={v_bi:.6}"
    );
}

/// In bidir mode the laser explicitly zeroes the bw wires.  Verify by
/// probing the source bundle's `re_bw_0` net directly.
#[test]
fn bidir_laser_zeroes_backward_wires() {
    let netlist = "\
.options enable_bidirectional=1
.optical_port src
Xl src fc_cw_laser power_mW=1.0 wavelength_nm=1550
.op
.end
";
    let net = parse_spice(netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let re_fw = r.node_voltage("src_re_fw_0").unwrap();
    let re_bw = r.node_voltage("src_re_bw_0").unwrap();
    let im_bw = r.node_voltage("src_im_bw_0").unwrap();
    let p_in = 1e-3_f64;
    assert!(
        (re_fw - p_in.sqrt()).abs() < 1e-9,
        "fw amplitude should be √(1 mW); got {re_fw}"
    );
    assert!(re_bw.abs() < 1e-12, "bw real should be 0; got {re_bw}");
    assert!(im_bw.abs() < 1e-12, "bw imag should be 0; got {im_bw}");
}

/// Drive the bw direction directly through an external voltage source on
/// `wg_out_re_bw_0` and verify the waveguide propagates it back to the
/// `src_re_bw_0` wire — the bidir reciprocity of the waveguide.
#[test]
fn bidir_waveguide_propagates_backward_path() {
    let netlist = "\
.options enable_bidirectional=1
.optical_port src
.optical_port wg_out
* Tie fw inputs to 0 — only backward drive.
Vsrc_fw_re src_re_fw_0 0 DC 0
Vsrc_fw_im src_im_fw_0 0 DC 0
Vsrc_wl    src_wl_0    0 DC 1.55e-6
* Drive bw from the OUTPUT side.
Vbw_re wg_out_re_bw_0 0 DC 1.0
Vbw_im wg_out_im_bw_0 0 DC 0.0
Vout_fw_re wg_out_re_fw_0 0 DC 0
Vout_fw_im wg_out_im_fw_0 0 DC 0
Vout_wl    wg_out_wl_0    0 DC 1.55e-6
Xwg src wg_out fc_waveguide L_um=100 n_g=4.2 alpha_dB_cm=2.0
.op
.end
";
    let net = parse_spice(netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let src_re_bw = r.node_voltage("src_re_bw_0").unwrap();
    let src_im_bw = r.node_voltage("src_im_bw_0").unwrap();
    let amp = (src_re_bw * src_re_bw + src_im_bw * src_im_bw).sqrt();
    // Expected amplitude: t_amp · 1.0 = exp(-α·L/2) ≈ exp(-2·100·ln10/20·100e-6/2)
    //                            = exp(-1.151e-3) ≈ 0.998850.
    let t_expected = (-2.0_f64 * 100.0 * std::f64::consts::LN_10 / 20.0 * 100e-6 / 2.0).exp();
    assert!(
        (amp - t_expected).abs() < 1e-5,
        "bw amplitude through waveguide = {amp:.6}; expected {t_expected:.6}"
    );
}
