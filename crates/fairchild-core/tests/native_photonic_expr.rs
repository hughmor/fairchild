//! Tier-1 runtime-loadable optical models: a `.model … fc_phase_shifter_expr`
//! card with declarative constitutive expressions (`dneff`, `dalpha`) over the
//! bias `V`, evaluated per NR-iterate — no recompile. Verifies a linear
//! expression map reproduces the equivalent hard-coded `fc_pn_ps`, and that a
//! nonlinear map takes effect.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

fn out_re_im(netlist: &str) -> (f64, f64) {
    let net = parse_spice(netlist).unwrap();
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);
    let r = dc_op_nr_with_registry(&net, &reg).expect("DC OP converges");
    (
        r.node_voltage("out0_re_0").unwrap(),
        r.node_voltage("out0_im_0").unwrap(),
    )
}

/// A declarative linear map `dneff = "5e-5*V"` reproduces `fc_pn_ps` with
/// `dn_dv=5e-5` (same optics) — proving the expression path matches the
/// hard-coded drive.
#[test]
fn expr_linear_map_matches_fc_pn_ps() {
    let expr = out_re_im(
        "\
.optical_port ch0
.optical_port out0
.model myps fc_phase_shifter_expr dneff=\"5.0e-5*V\" g_pn=1e-3 alpha_db_cm=20
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xps ch0 out0 a 0 myps pin_at_ref=1
Vb a 0 DC -1.5
.op
.end
",
    );
    let direct = out_re_im(
        "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xps ch0 out0 a 0 fc_pn_ps dn_dv=5.0e-5 g_pn=1e-3 alpha_db_cm=20 pin_at_ref=1
Vb a 0 DC -1.5
.op
.end
",
    );
    assert!(
        (expr.0 - direct.0).abs() < 1e-10 && (expr.1 - direct.1).abs() < 1e-10,
        "expr map {expr:?} should match fc_pn_ps {direct:?}"
    );
}

/// A nonlinear declarative map (with a SPICE-suffixed coefficient and a `V²`
/// term) parses and takes effect — the output differs from the linear-only map
/// at the same bias.
#[test]
fn expr_nonlinear_map_takes_effect() {
    let common = |dneff: &str| {
        format!(
            "\
.optical_port ch0
.optical_port out0
.model myps fc_phase_shifter_expr dneff={dneff} g_pn=1e-3 alpha_db_cm=0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xps ch0 out0 a 0 myps pin_at_ref=1
Vb a 0 DC 2.0
.op
.end
"
        )
    };
    let linear = out_re_im(&common("\"5.0e-5*V\""));
    let quad = out_re_im(&common("\"5.0e-5*V + 1.0e-5*V*V\""));
    // Same |A| (alpha=0 ⇒ lossless), different phase ⇒ different (re, im).
    let amp_lin = (linear.0 * linear.0 + linear.1 * linear.1).sqrt();
    let amp_quad = (quad.0 * quad.0 + quad.1 * quad.1).sqrt();
    assert!(
        (amp_lin - amp_quad).abs() < 1e-9,
        "lossless: amplitudes equal"
    );
    assert!(
        (linear.0 - quad.0).abs() > 1e-4 || (linear.1 - quad.1).abs() > 1e-4,
        "the V² term must change the phase: linear={linear:?} quad={quad:?}"
    );
}
