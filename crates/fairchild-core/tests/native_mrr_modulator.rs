//! Native micro-ring modulator regression test.
//!
//! Mirrors the user-facing example at `examples/photonic/native_mrr_modulator.sp`
//! to keep the published example honest as the photonic codebase evolves.
//!
//! Topology:
//!   fc_cw_laser → fc_waveguide → fc_dcoupler → fc_waveguide → fc_photodetector
//!                                  ↑   ↓                          │
//!                                  └─ fc_pn_ps (ring) ─┘           │
//!                                                                  ↓
//!                                                         R_load to V_bias
//!
//! Asserts:
//!   - DC OP converges with V_pn = 0 (on-resonance, deep notch).
//!   - Transient with V_pn = 0 → 4 V → 0 shows the photodetector output
//!     swinging from the notch level to ~unity transmission and back.
//!
//! Uses native Rust photonic devices only — no `.osdi`, no Verilog-A.

use fairchild_core::{
    dc_op_nr_with_registry, options::SimOptions, tran::IntegratorMode,
    tran_nr_with_registry_var_opts, DeviceRegistry,
};
use fairchild_parser::parse_spice;

/// Build the netlist programmatically so the test doesn't rely on the
/// example file's exact path or formatting.
fn mrr_netlist(vmod_dc: f64) -> String {
    format!(
        "\
* Native MRR for regression
.optical_port laser_out
.optical_port wg1_out
.optical_port dc_b
.optical_port dc_c
.optical_port pn_in
.optical_port pd_in

Xlaser laser_out fc_cw_laser power_mW=1.0 wavelength_nm=1550

Xwg1 laser_out wg1_out fc_waveguide L_um=50 n_g=4.2 alpha_dB_cm=2.0

Xdc wg1_out dc_b dc_c pn_in fc_dcoupler kappa_L=0.336

Xpn pn_in dc_b vmod 0 fc_pn_ps L_um=500 V_pi_L=2e-3 g_pn=1e-3 alpha_dB_cm=10 pin_at_ref=1

Xwg2 dc_c pd_in fc_waveguide L_um=50 n_g=4.2 alpha_dB_cm=2.0

Xpd pd_in pd_anode 0 fc_photodetector responsivity=0.8 i_dark_a=1e-9 r_shunt=1Meg

Vbias bias 0 DC 1.0
Rload pd_anode bias 1k
Vmod vmod 0 DC {vmod_dc}

.op
.end
"
    )
}

/// V_pn = 0: ring is on-resonance, transmission low, PD output near V_bias.
#[test]
fn mrr_dc_op_at_resonance_deep_notch() {
    let net = parse_spice(&mrr_netlist(0.0)).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let v_pd = r.node_voltage("pd_anode").unwrap();
    // On-resonance: should be much closer to V_bias (1.0 V) than to the
    // off-resonance asymptote (~1.8 V).
    assert!(
        v_pd < 1.3,
        "V(pd_anode) = {v_pd:.3} should be near V_bias=1V at resonance, far below 1.8V off-res"
    );
}

/// V_pn = 4 V (= Vπ): ring is fully off-resonance, near-unity transmission.
#[test]
fn mrr_dc_op_off_resonance_full_transmission() {
    let net = parse_spice(&mrr_netlist(4.0)).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let v_pd = r.node_voltage("pd_anode").unwrap();
    assert!(
        v_pd > 1.6,
        "V(pd_anode) = {v_pd:.3} should be ≥ 1.6 V at V_π (full transmission)"
    );
}

/// Wavelength sweep through resonance.  The PN-PS's propagation phase
/// must depend on λ (not just V_pn) — otherwise sweeping λ produces a
/// flat response and the ring isn't really a resonator.  Sample three
/// wavelengths: on-resonance (1550 nm = the device's reference), and
/// well off-resonance on either side.  Assert that the off-resonance
/// values are close to the high-transmission asymptote (~1.8 V) and the
/// on-resonance value is in the notch (~1.1 V).
#[test]
fn mrr_wavelength_sweep_resolves_resonance() {
    let build = |wl_nm: f64| {
        let s =
            mrr_netlist(0.0).replace("wavelength_nm=1550\n", &format!("wavelength_nm={wl_nm}\n"));
        parse_spice(&s).unwrap()
    };
    let on_res = build(1550.0);
    let off_res_lo = build(1549.5);
    let off_res_hi = build(1550.5);

    let r_on = dc_op_nr_with_registry(&on_res, &DeviceRegistry::new()).unwrap();
    let r_lo = dc_op_nr_with_registry(&off_res_lo, &DeviceRegistry::new()).unwrap();
    let r_hi = dc_op_nr_with_registry(&off_res_hi, &DeviceRegistry::new()).unwrap();

    let v_on = r_on.node_voltage("pd_anode").unwrap();
    let v_lo = r_lo.node_voltage("pd_anode").unwrap();
    let v_hi = r_hi.node_voltage("pd_anode").unwrap();

    // Off-resonance: PD output should be near the full-transmission level.
    assert!(
        v_lo > 1.6,
        "off-res low (λ=1549.5): V(pd_anode)={v_lo}; expected ≳ 1.6 V"
    );
    assert!(
        v_hi > 1.6,
        "off-res high (λ=1550.5): V(pd_anode)={v_hi}; expected ≳ 1.6 V"
    );
    // On-resonance: deep notch — must be visibly different from off-res.
    assert!(
        v_on < v_lo - 0.3,
        "on-res V(pd_anode)={v_on:.3} should be ≥ 0.3 V below off-res-low {v_lo:.3}"
    );
    assert!(
        v_on < v_hi - 0.3,
        "on-res V(pd_anode)={v_on:.3} should be ≥ 0.3 V below off-res-high {v_hi:.3}"
    );
}

/// Full transient with a 100 ns rise / 100 ns fall pulse — verify the PD
/// output traces the pulse and the extinction-ratio swing.
#[test]
fn mrr_transient_traces_pn_pulse() {
    let netlist_str = "\
* MRR transient
.optical_port laser_out
.optical_port wg1_out
.optical_port dc_b
.optical_port dc_c
.optical_port pn_in
.optical_port pd_in

Xlaser laser_out fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xwg1 laser_out wg1_out fc_waveguide L_um=50 n_g=4.2 alpha_dB_cm=2.0
Xdc wg1_out dc_b dc_c pn_in fc_dcoupler kappa_L=0.336
Xpn pn_in dc_b vmod 0 fc_pn_ps L_um=500 V_pi_L=2e-3 g_pn=1e-3 alpha_dB_cm=10 pin_at_ref=1
Xwg2 dc_c pd_in fc_waveguide L_um=50 n_g=4.2 alpha_dB_cm=2.0
Xpd pd_in pd_anode 0 fc_photodetector responsivity=0.8 i_dark_a=1e-9 r_shunt=1Meg

Vbias bias 0 DC 1.0
Rload pd_anode bias 1k
Vmod vmod 0 PULSE(0 4 100n 100n 100n 800n 2u)

.options method=gear
.tran 5n 2u
.end
";
    let net = parse_spice(netlist_str).unwrap();
    let mut opts = SimOptions::from_netlist(&net);
    opts.method = IntegratorMode::Gear;
    let r = tran_nr_with_registry_var_opts(&net, 5e-9, 2e-6, &DeviceRegistry::new(), &opts)
        .expect("tran must complete");
    let pd = r.node_voltages.get("pd_anode").unwrap();
    let v_min = pd.iter().cloned().fold(f64::INFINITY, f64::min);
    let v_max = pd.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    // V_pn = 0: V_pd ≈ 1.09 V (notch);  V_pn = 4V: V_pd ≈ 1.80 V (full).
    // The pulse spends meaningful time at both extremes, so v_min and v_max
    // should bracket the modulation swing closely.
    assert!(
        v_min < 1.2,
        "min(V(pd_anode)) = {v_min:.3}, expected ≈ 1.09 V at resonance"
    );
    assert!(
        v_max > 1.7,
        "max(V(pd_anode)) = {v_max:.3}, expected ≈ 1.80 V at full transmission"
    );
    // Modulation depth ≥ 500 mV — a clean, visible swing.
    assert!(
        v_max - v_min > 0.5,
        "swing {:.3} V too small; should be ≥ 0.5 V",
        v_max - v_min
    );
}
