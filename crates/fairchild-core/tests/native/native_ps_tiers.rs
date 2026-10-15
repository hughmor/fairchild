//! A tier is a set of parameters, not a model.
//!
//! Every `fc_pn_ps*` card is `fc_pn_ps_full` with terms switched off. This file
//! is the check on that claim, and it is the reason the three separate drive
//! models could be deleted: if a tier really is its parameters, then driving
//! the full card with a tier's own parameters has to reproduce that tier
//! *exactly* — the same code reached by a different route, not a close match.
//!
//! The tier parameters come from `Preset::params`, so this cannot drift into a
//! second copy of the table. Add a row to `PRESETS` and it is checked here.

use fairchild_core::models::PRESETS;
use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

/// Solve one bias point and return `(out_re, out_im, I(vb))`.
///
/// `r_series=40` is on every deck: it is the one parameter that changes the
/// device's *structure* rather than a coefficient, so a tier that reproduced
/// the coefficients but allocated its rows differently would still fail here.
fn probe(card: &str, extra: &str, vbias: f64) -> (f64, f64, f64) {
    let netlist = format!(
        "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xpn ch0 out0 a 0 {card} L_um=500 r_series=40 {extra}
Vb a 0 DC {vbias}
.op
"
    );
    let net = parse_spice(&netlist).unwrap_or_else(|e| panic!("{card}: parse: {e}"));
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new())
        .unwrap_or_else(|e| panic!("{card} at {vbias} V: {e}"));
    (
        r.node_voltage("out0_re_0").expect("out0_re_0"),
        r.node_voltage("out0_im_0").expect("out0_im_0"),
        r.vsrc_current("vb").expect("I(vb)"),
    )
}

/// `fc_pn_ps_full` plus a tier's parameters is that tier, bit for bit.
///
/// This pins the *mechanism*: there is no per-tier code left, only a parameter
/// list, so a card is reachable from the full drive and nothing else.
///
/// **What it cannot catch** is a wrong list, because both sides read the same
/// one — delete an entry from `T_INJ` and this still passes. The tiers are
/// anchored against their physics in
/// `each_tier_models_the_half_of_the_junction_it_claims` below, which is where
/// an absolute reference lives. The pair is deliberate.
///
/// Both bias signs, because the tiers differ in which half of the junction they
/// model: at −1.5 V the depletion terms carry the answer and at +0.6 V the
/// injection ones do.
#[test]
fn every_pn_tier_is_the_full_drive_with_terms_switched_off() {
    for p in PRESETS {
        // The heater cards have different terminals, so they are not reachable
        // from `fc_pn_ps_full` and are covered by the composite test below.
        if !p.name.starts_with("fc_pn_ps") || p.name == "fc_pn_ps_full" {
            continue;
        }
        let extra = p
            .params()
            .iter()
            .map(|(k, v)| format!("{k}={v:e}"))
            .collect::<Vec<_>>()
            .join(" ");

        for vbias in [-1.5, 0.6] {
            let tier = probe(p.name, "", vbias);
            let full = probe("fc_pn_ps_full", &extra, vbias);
            assert_eq!(
                tier, full,
                "{} at {vbias} V is not 'fc_pn_ps_full {extra}'. The tier is \
                 supposed to be those parameters and nothing else.",
                p.name
            );
        }
    }
}

/// ABSOLUTE anchor on each tier: the phase it produces, against the closed form
/// of the physics that tier claims to model.
///
/// The reference here is `2π·L·Δn/λ` computed in the test from the extracted
/// coefficients, not from the preset table, so a wrong row fails. That is the
/// half the equivalence test above cannot do.
///
/// What each tier claims:
///
/// * linear and depletion — `Δn = dn_dv·V`, the *same* number at +0.6 V as at
///   −0.6 V reflected. An injection term left switched on is `exp(V/n·V_T) − 1`,
///   about 4e9 at +0.6 V, so it is not a subtle failure.
/// * injection — no reverse term at all, so reverse bias moves the phase by
///   nothing. Leaving `dn_dv_rev` on gives 1.5 V × 5.024e-5 over 500 µm, which
///   is 0.153 rad and 24000x the tolerance here.
#[test]
fn each_tier_models_the_half_of_the_junction_it_claims() {
    // `r_series=0` here, unlike the equivalence test above. With it on, the
    // junction sees `V/(1 + g_pn·r_s)` rather than `V`, which is a true
    // statement about the divider and a distraction from the closed form being
    // anchored. Leaving it in cost exactly the factor 1/1.04.
    let phase = |card: &str, vbias: f64| {
        let (re, im, _) = probe(card, "r_series=0", vbias);
        im.atan2(re)
    };

    // Reverse bias moves the depletion tiers by dn_dv·ΔV and nothing else.
    const DN_DV: f64 = 5.024e-5; // the one extracted reverse coefficient
    const L_M: f64 = 500e-6;
    const LAMBDA: f64 = 1550e-9;
    let want = 2.0 * std::f64::consts::PI * L_M * DN_DV * 1.0 / LAMBDA;
    for card in ["fc_pn_ps", "fc_pn_ps_cap"] {
        let got = phase(card, -0.5) - phase(card, -1.5);
        assert!(
            (got.abs() - want).abs() < 1e-6,
            "{card}: 1 V of reverse bias should rotate by {want:.6} rad \
             (2π·L·dn_dv·ΔV/λ); got {got:.6}"
        );
    }

    // And the same tiers are linear *through* zero: no injection term, so a
    // forward volt is worth exactly what a reverse volt is.
    for card in ["fc_pn_ps", "fc_pn_ps_cap"] {
        let forward = phase(card, 0.6) - phase(card, 0.1);
        let reverse = phase(card, -0.1) - phase(card, -0.6);
        assert!(
            (forward - reverse).abs() < 1e-6,
            "{card} must be linear through 0 V: 0.1→0.6 V rotates by \
             {forward:.6} rad and −0.6→−0.1 V by {reverse:.6}. An injection \
             term left on would make the forward side enormous."
        );
    }

    // The injection tier has no depletion index, so reverse bias does nothing.
    let inj = phase("fc_pn_ps_inj", -0.5) - phase("fc_pn_ps_inj", -1.5);
    assert!(
        inj.abs() < 1e-9,
        "fc_pn_ps_inj models forward injection only, so 1 V of reverse bias \
         must rotate by nothing; got {inj:.9} rad"
    );
}

/// The heater axis is orthogonal to the tier axis, so adding a heater must not
/// disturb the PN answer when the heater is off.
///
/// `fc_pn_th_ps_*` is `fc_pn_ps_*` plus a heater, and an unpowered heater is a
/// resistor to ground that the optics never see. This is what lets one tier
/// table serve both families.
#[test]
fn a_heater_at_rest_does_not_change_the_pn_answer() {
    let with_heater = |card: &str, vbias: f64| {
        let netlist = format!(
            "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xpn ch0 out0 a 0 hp 0 {card} L_um=500 r_series=40
Vb a 0 DC {vbias}
Vh hp 0 DC 0
.op
"
        );
        let net = parse_spice(&netlist).unwrap();
        let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).unwrap();
        (
            r.node_voltage("out0_re_0").unwrap(),
            r.node_voltage("out0_im_0").unwrap(),
            r.vsrc_current("vb").unwrap(),
        )
    };

    for tier in ["", "_cap", "_inj", "_full"] {
        for vbias in [-1.5, 0.6] {
            let bare = probe(&format!("fc_pn_ps{tier}"), "", vbias);
            let heated = with_heater(&format!("fc_pn_th_ps{tier}"), vbias);
            assert_eq!(
                bare, heated,
                "fc_pn_th_ps{tier} with its heater at 0 V should equal \
                 fc_pn_ps{tier} at {vbias} V"
            );
        }
    }
}

/// Every PN card answers to `LEVEL` on a `.model` card, and the two routes have
/// to reach the same preset.
///
/// Two lists of tier names would be two chances to disagree, so `LEVEL` builds
/// a card name and looks it up. This checks the arithmetic of that.
#[test]
fn level_selects_the_same_preset_as_the_suffixed_card_name() {
    for (level, suffix) in [(1, ""), (2, "_cap"), (3, "_inj"), (4, "_full")] {
        let netlist = format!(
            "\
.optical_port ch0
.optical_port out0
.model myps fc_pn_ps LEVEL={level}
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xpn ch0 out0 a 0 myps L_um=500 r_series=40
Vb a 0 DC 0.6
.op
"
        );
        let net = parse_spice(&netlist).unwrap();
        let mut reg = DeviceRegistry::new();
        reg.register_builtin_models(&net.models);
        let r = dc_op_nr_with_registry(&net, &reg).unwrap();
        let carded = (
            r.node_voltage("out0_re_0").unwrap(),
            r.node_voltage("out0_im_0").unwrap(),
            r.vsrc_current("vb").unwrap(),
        );
        assert_eq!(
            carded,
            probe(&format!("fc_pn_ps{suffix}"), "", 0.6),
            "LEVEL={level} should be fc_pn_ps{suffix}"
        );
    }
}

/// `tau_th` is an ordinary parameter now, so a PN composite can have a thermal
/// lag. It could not before: `fc_pn_th_ps_*` was always built with the
/// instantaneous heater, and `fc_thermal_ps_rc` had no junction.
///
/// The lag is what the row buys, so check the row appears and that DC is
/// unchanged by it — at steady state the filter's state equation reduces to
/// `T = P`, so an operating point must not depend on `tau_th` at all.
#[test]
fn a_pn_composite_can_have_a_thermal_lag() {
    let op = |extra: &str| {
        let netlist = format!(
            "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xpn ch0 out0 a 0 hp 0 fc_pn_th_ps_cap L_um=500 {extra}
Vb a 0 DC -1.5
Vh hp 0 DC 2.0
.op
"
        );
        let net = parse_spice(&netlist).unwrap();
        let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new())
            .unwrap_or_else(|e| panic!("fc_pn_th_ps_cap {extra}: {e}"));
        (
            r.node_voltage("out0_re_0").unwrap(),
            r.node_voltage("out0_im_0").unwrap(),
        )
    };
    let instant = op("");
    let lagged = op("tau_th=20u");
    // Not `assert_eq!` like the tiers above, and the difference is the point:
    // those reach identical code, while this adds a row and changes the order
    // the solve sums in. 1e-15 is the arithmetic, not the physics.
    assert!(
        (instant.0 - lagged.0).abs() < 1e-15 && (instant.1 - lagged.1).abs() < 1e-15,
        "tau_th must not move the operating point — the state equation reduces \
         to T = P at DC, whichever side of it you come from. Got {instant:?} \
         without the lag and {lagged:?} with it"
    );
}
