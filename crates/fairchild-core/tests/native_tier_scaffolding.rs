//! Sanity tests for the `level=` tier-selection parameter on fc_pn_ps and
//! fc_thermal_ps.  At this commit (B.1), L2/L3 are scaffolded but fall back
//! to L1 stamping — these tests verify the keyword is accepted, in range
//! validation works, and L1 behaviour is unchanged.

use fairchild_core::{DeviceRegistry, dc_op_nr_with_registry};
use fairchild_parser::parse_spice;

/// Default `level=1` and explicit `level=1` must produce identical results.
#[test]
fn pn_ps_default_and_explicit_l1_agree() {
    let nl_default = "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xpn ch0 out0 vmod 0 fc_pn_ps L_um=100 g_pn=1e-3
Vmod vmod 0 DC 1.0
.op
.end
";
    let nl_l1 = nl_default.replace("g_pn=1e-3", "g_pn=1e-3 level=1");
    let r1 = dc_op_nr_with_registry(&parse_spice(nl_default).unwrap(), &DeviceRegistry::new()).unwrap();
    let r2 = dc_op_nr_with_registry(&parse_spice(&nl_l1).unwrap(), &DeviceRegistry::new()).unwrap();
    let v1 = r1.vsrc_current("vmod").unwrap().abs();
    let v2 = r2.vsrc_current("vmod").unwrap().abs();
    assert!((v1 - v2).abs() < 1e-12, "level=1 must match default; v1={v1} v2={v2}");
}

/// `level=2` is accepted (with a stderr warning) but currently falls back
/// to L1 behaviour.  This pins the contract so users can ship netlists
/// that opt into L2/L3 before the implementation lands.
#[test]
fn pn_ps_level2_falls_back_to_l1_for_now() {
    let nl_l1 = "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xpn ch0 out0 vmod 0 fc_pn_ps L_um=100 g_pn=1e-3 level=1
Vmod vmod 0 DC 1.0
.op
.end
";
    let nl_l2 = nl_l1.replace("level=1", "level=2 c_j0=20f tau_carrier=1n");
    let r1 = dc_op_nr_with_registry(&parse_spice(nl_l1).unwrap(), &DeviceRegistry::new()).unwrap();
    let r2 = dc_op_nr_with_registry(&parse_spice(&nl_l2).unwrap(), &DeviceRegistry::new()).unwrap();
    let v1 = r1.vsrc_current("vmod").unwrap().abs();
    let v2 = r2.vsrc_current("vmod").unwrap().abs();
    assert!((v1 - v2).abs() < 1e-9,
        "level=2 should fall back to L1 (no L2 impl yet); v1={v1} v2={v2}");
}

/// `level=` is accepted on fc_thermal_ps too.
#[test]
fn thermal_ps_level_keyword_is_accepted() {
    let netlist = "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xth ch0 out0 heat_p 0 fc_thermal_ps r_heater=1k p_pi=10m level=2 tau_th=10u
Vh heat_p 0 DC 1.0
.op
.end
";
    let net = parse_spice(netlist).unwrap();
    let _r = dc_op_nr_with_registry(&net, &DeviceRegistry::new())
        .expect("thermal_ps should accept level= and still solve");
}
