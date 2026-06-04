//! The active-photonic LEVEL system: a `.model <name> <family> LEVEL=<n>` card
//! selects a phase-shifter variant (à la MOSFET LEVEL), equivalent to the
//! corresponding direct `fc_*` device name. Verifies the LEVEL dispatch maps to
//! the same physics as the dedicated device classes.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

/// Build a registry with the netlist's `.model` cards registered (mirrors what
/// the analysis entry points do).
fn registry_for(net: &fairchild_parser::Netlist) -> DeviceRegistry {
    let mut r = DeviceRegistry::new();
    r.register_builtin_models(&net.models);
    r
}

fn out_re_im(netlist: &str) -> (f64, f64) {
    let net = parse_spice(netlist).unwrap();
    let reg = registry_for(&net);
    let r = dc_op_nr_with_registry(&net, &reg).expect("DC OP converges");
    (
        r.node_voltage("out0_re_0").unwrap(),
        r.node_voltage("out0_im_0").unwrap(),
    )
}

/// `.model … fc_pn_ps LEVEL=2` ≡ the direct `fc_pn_ps_cap` device.
#[test]
fn level2_pn_ps_matches_cap() {
    let leveled = out_re_im(
        "\
.optical_port ch0
.optical_port out0
.model myps fc_pn_ps LEVEL=2
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xpn ch0 out0 a 0 myps pin_at_ref=1
Vb a 0 DC -1.0
.op
.end
",
    );
    let direct = out_re_im(
        "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xpn ch0 out0 a 0 fc_pn_ps_cap pin_at_ref=1
Vb a 0 DC -1.0
.op
.end
",
    );
    assert!(
        (leveled.0 - direct.0).abs() < 1e-12 && (leveled.1 - direct.1).abs() < 1e-12,
        "LEVEL=2 {leveled:?} should equal fc_pn_ps_cap {direct:?}"
    );
}

/// A bare `.model … fc_pn_ps` (no LEVEL) ≡ the plain `fc_pn_ps`.
#[test]
fn level_absent_defaults_to_level1() {
    let leveled = out_re_im(
        "\
.optical_port ch0
.optical_port out0
.model myps fc_pn_ps
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xpn ch0 out0 a 0 myps pin_at_ref=1
Vb a 0 DC -1.0
.op
.end
",
    );
    let direct = out_re_im(
        "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xpn ch0 out0 a 0 fc_pn_ps pin_at_ref=1
Vb a 0 DC -1.0
.op
.end
",
    );
    assert!(
        (leveled.0 - direct.0).abs() < 1e-12 && (leveled.1 - direct.1).abs() < 1e-12,
        "LEVEL-absent {leveled:?} should equal fc_pn_ps {direct:?}"
    );
}

/// `.model … fc_thermal_ps LEVEL=2` ≡ `fc_thermal_ps_rc`, and the model card's
/// own params (here `p_pi`) are honoured.
#[test]
fn level2_thermal_matches_rc_and_card_params_apply() {
    let leveled = out_re_im(
        "\
.optical_port ch0
.optical_port out0
.model myth fc_thermal_ps LEVEL=2 p_pi=20m
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xth ch0 out0 hp 0 myth
Vh hp 0 DC 2.0
.op
.end
",
    );
    let direct = out_re_im(
        "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xth ch0 out0 hp 0 fc_thermal_ps_rc p_pi=20m
Vh hp 0 DC 2.0
.op
.end
",
    );
    assert!(
        (leveled.0 - direct.0).abs() < 1e-12 && (leveled.1 - direct.1).abs() < 1e-12,
        "LEVEL=2 thermal {leveled:?} should equal fc_thermal_ps_rc {direct:?}"
    );
}

/// `.model … fc_pn_th_ps LEVEL=4` ≡ the direct `fc_pn_th_ps_full`.
#[test]
fn level4_pn_thermal_matches_full() {
    let leveled = out_re_im(
        "\
.optical_port ch0
.optical_port out0
.model myfull fc_pn_th_ps LEVEL=4
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xpn ch0 out0 a 0 hp 0 myfull pin_at_ref=1
Vb a 0 DC -1.0
Vh hp 0 DC 1.0
.op
.end
",
    );
    let direct = out_re_im(
        "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xpn ch0 out0 a 0 hp 0 fc_pn_th_ps_full pin_at_ref=1
Vb a 0 DC -1.0
Vh hp 0 DC 1.0
.op
.end
",
    );
    assert!(
        (leveled.0 - direct.0).abs() < 1e-12 && (leveled.1 - direct.1).abs() < 1e-12,
        "LEVEL=4 PN+thermal {leveled:?} should equal fc_pn_th_ps_full {direct:?}"
    );
}

/// Param precedence (the ModelFactory/ParamSet plumbing): a model-card param is
/// overridden by the same param on the instance line.
#[test]
fn instance_param_overrides_model_card() {
    // Card dn_dv=1e-5; instance dn_dv=9e-5 → result must match a direct device
    // with dn_dv=9e-5 (instance wins), not 1e-5.
    let leveled = out_re_im(
        "\
.optical_port ch0
.optical_port out0
.model myps fc_pn_ps dn_dv=1e-5
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xps ch0 out0 a 0 myps pin_at_ref=1 dn_dv=9e-5
Vb a 0 DC -1.0
.op
.end
",
    );
    let direct = out_re_im(
        "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xps ch0 out0 a 0 fc_pn_ps pin_at_ref=1 dn_dv=9e-5
Vb a 0 DC -1.0
.op
.end
",
    );
    assert!(
        (leveled.0 - direct.0).abs() < 1e-12 && (leveled.1 - direct.1).abs() < 1e-12,
        "instance dn_dv should override card dn_dv: leveled={leveled:?} direct={direct:?}"
    );
    // And confirm it actually differs from the card-only value (sanity).
    let card_only = out_re_im(
        "\
.optical_port ch0
.optical_port out0
.model myps fc_pn_ps dn_dv=1e-5
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xps ch0 out0 a 0 myps pin_at_ref=1
Vb a 0 DC -1.0
.op
.end
",
    );
    assert!(
        (leveled.1 - card_only.1).abs() > 1e-4,
        "override should change the result vs card-only: {leveled:?} vs {card_only:?}"
    );
}
