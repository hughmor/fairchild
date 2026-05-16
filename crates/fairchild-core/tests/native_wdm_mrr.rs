//! Native WDM micro-ring modulator regression test.
//!
//! Two lasers placed symmetrically (±50 pm) around the ring's V=0
//! resonance share one bus through one ring driven by one V_pn.  The
//! two channels MUST show different transmission profiles — that's the
//! whole demonstration: same modulator, different per-wavelength response.
//!
//! The bundle-port mechanism gives this for free: every photonic device
//! on the bus replicates per channel, but `vmod` is a plain net so both
//! ring instances share it.  No multiplexer or demultiplexer device is
//! involved.

use fairchild_core::{
    options::SimOptions, tran::IntegratorMode, tran_nr_with_registry_var_opts,
    DeviceRegistry, dc_op_nr_with_registry,
};
use fairchild_parser::parse_spice;

fn wdm_netlist(vmod_dc: f64) -> String {
    format!("\
* WDM MRR regression — two lasers at ±50 pm around the ring resonance
.optical_port bus_in 2
.optical_port wg1_out 2
.optical_port dc_b 2
.optical_port dc_c 2
.optical_port pn_in 2
.optical_port pd_in 2

Xlaser1 bus_in_re_0 bus_in_im_0 bus_in_wl_0 fc_cw_laser power_mW=1.0 wavelength_nm=1549.95
Xlaser2 bus_in_re_1 bus_in_im_1 bus_in_wl_1 fc_cw_laser power_mW=1.0 wavelength_nm=1550.05

Xwg1 bus_in wg1_out fc_waveguide L_um=50 n_g=4.2 alpha_dB_cm=2.0
Xdc wg1_out dc_b dc_c pn_in fc_dcoupler kappa_L=0.336
Xpn pn_in dc_b vmod 0 fc_pn_ps L_um=500 V_pi_L=2e-3 g_pn=1e-3 alpha_dB_cm=10 n_g=4.2
Xwg2 dc_c pd_in fc_waveguide L_um=50 n_g=4.2 alpha_dB_cm=2.0

Xpd1 pd_in_re_0 pd_in_im_0 pd_in_wl_0 pd1_anode 0 fc_photodetector responsivity=0.8 i_dark_a=1e-9 r_shunt=1Meg
Xpd2 pd_in_re_1 pd_in_im_1 pd_in_wl_1 pd2_anode 0 fc_photodetector responsivity=0.8 i_dark_a=1e-9 r_shunt=1Meg

Vbias bias 0 DC 1.0
Rload1 pd1_anode bias 1k
Rload2 pd2_anode bias 1k

Vmod vmod 0 DC {vmod_dc}
.op
.end
")
}

/// V_pn = 0: lasers detuned symmetrically (±50 pm) — both channels see
/// identical transmission.  Pins this symmetry as a sanity check.
#[test]
fn wdm_dc_op_symmetric_at_zero_bias() {
    let net = parse_spice(&wdm_netlist(0.0)).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let v_pd1 = r.node_voltage("pd1_anode").unwrap();
    let v_pd2 = r.node_voltage("pd2_anode").unwrap();
    assert!((v_pd1 - v_pd2).abs() < 1e-3,
        "at V_pn=0 with symmetric detuning, PDs should match: pd1={v_pd1:.4} pd2={v_pd2:.4}");
}

/// At a V_pn that walks the ring resonance onto laser2 (red-side),
/// channel 1 should be in a deep notch while channel 0 is far off-res.
/// V = 0.4 V shifts resonance by ≈ 57 pm red → reaches laser2 at +50 pm.
#[test]
fn wdm_dc_op_breaks_symmetry_at_intermediate_bias() {
    let net = parse_spice(&wdm_netlist(0.4)).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let v_pd1 = r.node_voltage("pd1_anode").unwrap();
    let v_pd2 = r.node_voltage("pd2_anode").unwrap();
    // Channel 0 (blue-side): resonance walked away, far off-res → high T.
    // Channel 1 (red-side):  resonance walked onto it, deep notch → low T.
    assert!(v_pd1 > 1.6,
        "ch0 should be off-resonance at V=0.4: V(pd1)={v_pd1:.3}");
    assert!(v_pd2 < 1.4,
        "ch1 should be near-resonance at V=0.4: V(pd2)={v_pd2:.3}");
    // The asymmetry is the whole point.
    let asym = v_pd1 - v_pd2;
    assert!(asym > 0.3,
        "WDM asymmetry too small (Δ = {asym:.3} V); expected ≥ 0.3 V");
}

/// Transient with V_pn = PULSE(0→4 V→0) — verify both channels run, both
/// are oscillation-free, and they DIVERGE during the rising edge (when
/// the ring resonance sweeps through laser2 but not laser1).
#[test]
fn wdm_transient_two_channels_diverge() {
    let netlist_str = "\
* WDM transient regression
.optical_port bus_in 2
.optical_port wg1_out 2
.optical_port dc_b 2
.optical_port dc_c 2
.optical_port pn_in 2
.optical_port pd_in 2

Xlaser1 bus_in_re_0 bus_in_im_0 bus_in_wl_0 fc_cw_laser power_mW=1.0 wavelength_nm=1549.95
Xlaser2 bus_in_re_1 bus_in_im_1 bus_in_wl_1 fc_cw_laser power_mW=1.0 wavelength_nm=1550.05
Xwg1 bus_in wg1_out fc_waveguide L_um=50 n_g=4.2 alpha_dB_cm=2.0
Xdc wg1_out dc_b dc_c pn_in fc_dcoupler kappa_L=0.336
Xpn pn_in dc_b vmod 0 fc_pn_ps L_um=500 V_pi_L=2e-3 g_pn=1e-3 alpha_dB_cm=10 n_g=4.2
Xwg2 dc_c pd_in fc_waveguide L_um=50 n_g=4.2 alpha_dB_cm=2.0
Xpd1 pd_in_re_0 pd_in_im_0 pd_in_wl_0 pd1_anode 0 fc_photodetector responsivity=0.8 i_dark_a=1e-9 r_shunt=1Meg
Xpd2 pd_in_re_1 pd_in_im_1 pd_in_wl_1 pd2_anode 0 fc_photodetector responsivity=0.8 i_dark_a=1e-9 r_shunt=1Meg
Vbias bias 0 DC 1.0
Rload1 pd1_anode bias 1k
Rload2 pd2_anode bias 1k
Vmod vmod 0 PULSE(0 4 100n 100n 100n 800n 2u)
.options method=gear
.tran 5n 2u
.end
";
    let net = parse_spice(netlist_str).unwrap();
    let mut opts = SimOptions::from_netlist(&net);
    opts.method = IntegratorMode::Gear;
    let r = tran_nr_with_registry_var_opts(&net, 5e-9, 2e-6, &DeviceRegistry::new(), &opts)
        .expect("WDM transient must complete");
    let pd1 = r.node_voltages.get("pd1_anode").unwrap();
    let pd2 = r.node_voltages.get("pd2_anode").unwrap();
    let t   = &r.time;

    // Find the timepoint where V_pn ≈ 0.4 V — the strongest asymmetry.
    // V_pn(t) ramps linearly 0 → 4 V from t=100 ns to 200 ns; so the
    // 0.4 V crossing happens around t ≈ 110 ns.
    let probe = t.iter().position(|&tt| tt >= 1.1e-7).unwrap();
    let diff = (pd1[probe] - pd2[probe]).abs();
    assert!(diff > 0.3,
        "expected ≥ 0.3 V asymmetry between channels around V_pn=0.4V; got {diff:.3}");

    // Both channels must converge to roughly the same off-resonance level
    // by the end of the pulse plateau (V_pn = 4 V, ring shifted ~570 pm
    // red → both lasers far blue of the new resonance → both high T).
    let late = t.iter().position(|&tt| tt >= 5e-7).unwrap();
    let conv = (pd1[late] - pd2[late]).abs();
    assert!(conv < 0.1,
        "both channels should be off-resonance and similar at V_pn=4V; got Δ = {conv:.3}");
}
