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
    assert!((re - 0.030_670_69).abs() < 1e-7, "out_re={re:.8}");
    assert!((im - 0.006_334_093).abs() < 1e-8, "out_im={im:.9}");
}

/// Mild forward bias 0.1 V: injection Δn (negative) + forward FCA loss.
#[test]
fn pn_ps_full_forward_injection() {
    let r = solve("fc_pn_ps_full", 0.1, 1.0, "");
    let re = r.node_voltage("out0_re_0").unwrap();
    let im = r.node_voltage("out0_im_0").unwrap();
    assert!((re - (-0.004_819_248)).abs() < 1e-8, "out_re={re:.9}");
    assert!((im - (-0.005_559_183)).abs() < 1e-8, "out_im={im:.9}");
}

/// Current-parametrized injection (dn_di/da_di) is the same physics as the
/// (e−1)-factor form: I_fwd = i_sat·(e−1), so dn_di = dn_dv_inj/i_sat must
/// reproduce the forward-injection golden bit-for-bit (r_series = 0 here).
#[test]
fn pn_ps_full_dn_di_reparametrization_equivalent() {
    let r = solve(
        "fc_pn_ps_full",
        0.1,
        1.0,
        // defaults: dn_dv_inj=1.311e-4, da_dv_inj=150, i_sat=1e-12
        "dn_dv_inj=0 da_dv_inj=0 dn_di=1.311e8 da_di=1.5e14",
    );
    let re = r.node_voltage("out0_re_0").unwrap();
    let im = r.node_voltage("out0_im_0").unwrap();
    assert!((re - (-0.004_819_248)).abs() < 1e-8, "out_re={re:.9}");
    assert!((im - (-0.005_559_183)).abs() < 1e-8, "out_im={im:.9}");
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

// ── WDM back-action: the whole bus heats, not channel 0 (#51) ──────────────
//
// `FullPnDrive` used to read its optical back-action from `intensity_w.first()`
// — channel 0 alone — and then apply the resulting Δn/Δα to every channel. So a
// shared effect was driven by a single slot: 1/N of the truth with all channels
// lit, and *exactly zero* whenever channel 0 alone was dark.
//
// Both tests below use one wavelength on all four lasers, which makes the
// channels physically identical and lets slot-independence be asserted directly.
// `r_th` defaults to 0 (self-heating off), which is why nothing else caught this.

/// Four single-channel lasers muxed onto one 4-channel bus, into one junction.
fn wdm4(p: [f64; 4]) -> fairchild_core::NrResult {
    let netlist = format!(
        "\
.optical_port c0
.optical_port c1
.optical_port c2
.optical_port c3
.optical_port bus 4
.optical_port out 4
Xl0 c0 fc_cw_laser power_mW={} wavelength_nm=1550
Xl1 c1 fc_cw_laser power_mW={} wavelength_nm=1550
Xl2 c2 fc_cw_laser power_mW={} wavelength_nm=1550
Xl3 c3 fc_cw_laser power_mW={} wavelength_nm=1550
Xmux bus c0 c1 c2 c3 fc_mux
Xpn bus out a 0 fc_pn_ps_full L_um=500 alpha_dB_cm=2.0 dn_dt=1e-4 r_th=5000 beta_tpa=0
Vb a 0 DC 0
.op
.end
",
        p[0], p[1], p[2], p[3]
    );
    let net = parse_spice(&netlist).unwrap();
    dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP converges")
}

fn phase_of(r: &fairchild_core::NrResult, ch: usize) -> f64 {
    let re = r.node_voltage(&format!("out_re_{ch}")).unwrap();
    let im = r.node_voltage(&format!("out_im_{ch}")).unwrap();
    im.atan2(re)
}

/// ABSOLUTE anchor: the self-heating phase shift against a hand-computed power
/// budget. Tripling the lit power must triple Δn, and the increment is a closed
/// form, so this cannot be satisfied by the segment and the drive merely
/// agreeing with each other (`#32`'s shared-fault trap).
///
///   α      = 2.0 dB/cm · ln(10)/10 · 100 = 46.0517 Np/m
///   Δn     = dn_dt · R_th · α · L · ΣI = 1e-4 · 5000 · 46.0517 · 500e-6 · ΣI
///   Δφ     = 2π · L · Δ(Δn) / λ
///
/// With 10 mW per channel, going from one lit channel (ΣI = 10 mW) to three
/// (ΣI = 30 mW) gives Δ(Δn) = 2 · 1.15129e-4 and Δφ = 0.466757 rad.
/// Against the bug both runs have channel 0 dark, so both see zero heating and
/// Δφ is exactly 0.
#[test]
fn wdm_self_heating_sums_the_bus_and_matches_a_hand_power_budget() {
    let alpha = 2.0 * (10.0f64).ln() / 10.0 * 100.0;
    let l_m = 500e-6;
    let dn_per_watt = 1e-4 * 5000.0 * alpha * l_m;
    let expect = 2.0 * std::f64::consts::PI * l_m * (2.0 * dn_per_watt * 10e-3) / 1550e-9;

    let one = wdm4([0.0, 10.0, 0.0, 0.0]);
    let three = wdm4([0.0, 10.0, 10.0, 10.0]);
    let got = phase_of(&one, 1) - phase_of(&three, 1);

    assert!(
        (got.abs() - expect).abs() < 1e-4,
        "Δφ on channel 1 from tripling the lit power: got {got:.6} rad, \
         hand-computed {expect:.6} rad (α={alpha:.4} Np/m, Δn/W={dn_per_watt:.6})"
    );
}

/// Slot independence: with identical wavelengths, one lit channel must heat the
/// junction the same amount whichever slot it occupies. Against the bug, slot 0
/// is the only one that heats at all, so the two outputs differ.
#[test]
fn wdm_back_action_does_not_depend_on_which_slot_is_lit() {
    let ch0 = wdm4([10.0, 0.0, 0.0, 0.0]);
    let ch2 = wdm4([0.0, 0.0, 10.0, 0.0]);
    let (re0, im0) = (
        ch0.node_voltage("out_re_0").unwrap(),
        ch0.node_voltage("out_im_0").unwrap(),
    );
    let (re2, im2) = (
        ch2.node_voltage("out_re_2").unwrap(),
        ch2.node_voltage("out_im_2").unwrap(),
    );
    assert!(
        (re0 - re2).abs() < 1e-9 && (im0 - im2).abs() < 1e-9,
        "lit slot 0 gives ({re0:.9}, {im0:.9}), lit slot 2 gives ({re2:.9}, {im2:.9}) \
         — identical wavelengths and powers, so the back-action must not care which slot"
    );
}

/// The parser's `bundle_arity_for` is a second, hand-maintained list of a fact
/// the registry already knows, because `fairchild-core` depends on
/// `fairchild-parser` and the dispatch therefore cannot be unified. Two lists
/// are two chances to disagree: six LEVEL leaf names (`fc_pn_ps_full`,
/// `fc_pn_ps_inj`, `fc_pn_th_ps_{cap,inj,full}`, `fc_phase_shifter_expr`) were
/// registered as devices but absent from the arity list, so each was silently
/// `Scalar` and refused a WDM bus that the identical device accepted under its
/// family name. A test can see both lists even though the parser cannot.
#[test]
fn every_registered_photonic_model_declares_its_arity() {
    use fairchild_parser::{bundle_arity_for, BundleArity};
    // Deliberately single-channel: one laser emits one wavelength. Combine them
    // with `fc_mux` for WDM. Anything else `fc_*` is expected bundle-aware.
    const SCALAR_BY_DESIGN: &[&str] = &["fc_cw_laser", "fc_driven_laser"];

    let reg = DeviceRegistry::new();
    let mut wrong: Vec<&str> = reg
        .registered_names()
        .filter(|n| n.starts_with("fc_"))
        .filter(|n| !SCALAR_BY_DESIGN.contains(n))
        .filter(|n| bundle_arity_for(n) == BundleArity::Scalar)
        .collect();
    wrong.sort_unstable();
    assert!(
        wrong.is_empty(),
        "these photonic models are registered but fall through to BundleArity::Scalar, \
         so a multi-channel bundle is refused on them: {wrong:?}. Add them to \
         bundle_arity_for in fairchild-parser, or to SCALAR_BY_DESIGN here if being \
         single-channel is the intent."
    );
}

// ── cross-TPA is per channel, which one Δα cannot express (#54) ────────────

/// Two-photon absorption between *distinct* frequencies is twice as strong as
/// self-TPA, so on a loaded bus each channel sees its own loss:
///
///   α_TPA,j = β/A_eff · (I_j + 2·Σ_{k≠j} I_k) = β/A_eff · (2·Σ_k I_k − I_j)
///
/// With unequal channel powers those are different numbers, so this cannot be
/// satisfied by any single bus-wide Δα — which is exactly what
/// `OpticalPerturbation` used to carry. The previous best was the
/// no-cross-enhancement bound `β/A_eff·Σ_k I_k`, identical on every channel.
///
/// Everything else is switched off so the anchor is closed-form: V = 0 kills
/// depletion and injection, `r_th` defaults to 0 so nothing self-heats, and
/// `alpha_dB_cm=0` leaves TPA as the only loss.
#[test]
fn cross_tpa_gives_each_channel_its_own_loss() {
    const BETA: f64 = 7.9e-10;
    const A_EFF: f64 = 1.257e-13;
    const L_M: f64 = 500e-6;
    const P0: f64 = 20e-3;
    const P1: f64 = 5e-3;

    let netlist = format!(
        "\
.optical_port c0
.optical_port c1
.optical_port bus 2
.optical_port out 2
Xl0 c0 fc_cw_laser power_mW={p0} wavelength_nm=1550
Xl1 c1 fc_cw_laser power_mW={p1} wavelength_nm=1550
Xmux bus c0 c1 fc_mux
Xpn bus out a 0 fc_pn_ps_full L_um=500 alpha_dB_cm=0 beta_tpa={BETA} a_eff_m2={A_EFF}
Vb a 0 DC 0
.op
.end
",
        p0 = P0 * 1e3,
        p1 = P1 * 1e3
    );
    let net = parse_spice(&netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP converges");

    let out = |ch: usize| {
        let re = r.node_voltage(&format!("out_re_{ch}")).unwrap();
        let im = r.node_voltage(&format!("out_im_{ch}")).unwrap();
        re * re + im * im
    };
    // α_TPA,j = β/A_eff·(2·total − I_j); power survives exp(−α·L).
    let total = P0 + P1;
    let want = |i_j: f64| i_j * (-(BETA / A_EFF * (2.0 * total - i_j)) * L_M).exp();

    for (ch, i_j) in [(0usize, P0), (1usize, P1)] {
        let (got, expect) = (out(ch), want(i_j));
        assert!(
            (got - expect).abs() < 1e-9,
            "channel {ch}: got {got:.9} W, hand-computed {expect:.9} W"
        );
    }

    // And the discriminator: a single bus-wide Δα would attenuate both channels
    // by the same factor. Here the factors differ, which is the whole point.
    let f0 = out(0) / P0;
    let f1 = out(1) / P1;
    assert!(
        (f0 - f1).abs() > 1e-3,
        "the two channels are attenuated by {f0:.6} and {f1:.6} — a single Δα \
         would make these identical, so this test would not be testing anything"
    );
}
