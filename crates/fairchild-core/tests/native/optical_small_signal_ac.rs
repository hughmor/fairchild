//! An optical field must have a small-signal response in `.ac`.
//!
//! A device that stamps a coefficient frozen at the previous iterate leaves its
//! control column structurally empty and says so through
//! `Device::frozen_jacobian_columns`. The operating point repairs those columns
//! before differentiating, so `.sens` and `.tf` see through such a device.
//! `.ac` and `.pz` did not, and the column stayed empty: the small-signal
//! optical path through `fc_mzm`, `fc_xfer` and the laser was **exactly zero**,
//! with no diagnostic, because `gmin` and the passive network keep the matrix
//! non-singular and the solve succeeds (#113).
//!
//! Zero is the shape of the bug, so a test asserting "not zero" would be enough
//! to catch this one. It is not enough to keep the repair honest, so these
//! assert the value against the modulator's own closed form.

use fairchild_core::{ac_analysis, dc_op_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

const V_PI: f64 = 3.0;
const E_R: f64 = 1e3;

/// `fc_mzm`'s documented transfer, in the test rather than from the model:
/// `T(V) = α·[(1 − 1/E_r)·(1 + cos(πV/V_π))/2 + 1/E_r]`, amplitude `√T`.
fn t_amp(v: f64) -> f64 {
    let phi = std::f64::consts::PI * v / V_PI;
    (((1.0 - 1.0 / E_R) * (1.0 + phi.cos()) / 2.0) + 1.0 / E_R).sqrt()
}

/// `d(√T)/dV`, differentiated by hand from the same expression.
fn dt_amp_dv(v: f64) -> f64 {
    let phi = std::f64::consts::PI * v / V_PI;
    let dt_dv = (1.0 - 1.0 / E_R) * (-std::f64::consts::PI / V_PI) * phi.sin() / 2.0;
    dt_dv / (2.0 * t_amp(v))
}

fn deck(bias: f64, ac: &str) -> String {
    format!(
        "* MZM small-signal optical response\n\
         .optical_port in0\n\
         .optical_port out0\n\
         Xl0 in0 fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
         Xmzm in0 out0 vsig 0 fc_mzm V_pi=3.0 alpha=1.0 e_r=1k\n\
         Vsig vsig 0 DC {bias} {ac}\n"
    )
}

/// The small-signal optical response of a modulator is the slope of its own
/// transfer curve. Asserted as a ratio to the operating-point field, so the
/// laser power and any global phase cancel and what is left is `t'/t`, a pure
/// closed form.
#[test]
fn an_mzm_has_the_small_signal_response_its_transfer_curve_implies() {
    // Quadrature, where the slope is steepest and an error is most visible.
    for bias in [0.75, 1.5, 2.25] {
        let reg = DeviceRegistry::new();
        let op = dc_op_nr_with_registry(&parse_spice(&deck(bias, "")).unwrap(), &reg)
            .expect("operating point");
        let field_dc = op.node_voltage("out0_re_0").unwrap();

        let ac = ac_analysis(
            &parse_spice(&deck(bias, "AC 1")).unwrap(),
            &[1e6],
            Some("vsig"),
            &reg,
        )
        .expect("ac sweep");
        let field_ac = ac.magnitude("out0_re_0", 0).expect("optical output");

        let want = (dt_amp_dv(bias) / t_amp(bias)).abs();
        let got = field_ac / field_dc.abs();
        assert!(
            (got - want).abs() < 1e-6 * want.max(1.0),
            "bias {bias} V: |dE/dV|/|E| = {got:.6}, closed form {want:.6}. \
             Exactly zero means the frozen control column was never repaired"
        );
    }
}

/// The bias where the closed form says the response vanishes must also read
/// zero. Without it, "nonzero everywhere" would pass a stamp of the wrong
/// magnitude, and the null is the one point a wrong slope cannot fake.
#[test]
fn an_mzm_has_no_response_at_its_transfer_peak() {
    let reg = DeviceRegistry::new();
    let ac = ac_analysis(
        &parse_spice(&deck(0.0, "AC 1")).unwrap(),
        &[1e6],
        Some("vsig"),
        &reg,
    )
    .expect("ac sweep");
    let field_ac = ac.magnitude("out0_re_0", 0).expect("optical output");
    assert!(
        field_ac < 1e-9,
        "at V = 0 the transfer is flat (sin(0) = 0), so the response must \
         vanish; got {field_ac:.3e}"
    );
}

/// The photodetector converts that field response to a current, and only an
/// *intensity* change reaches it. This is the end-to-end statement: a link's
/// AC response is not zero, and it is zero at the transfer peak for the same
/// reason the field response is.
#[test]
fn the_detected_response_follows_the_field_response() {
    let link = |bias: f64| {
        let net = format!(
            "* MZM link, detected\n\
             .optical_port in0\n\
             .optical_port out0\n\
             Xl0 in0 fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
             Xmzm in0 out0 vsig 0 fc_mzm V_pi=3.0 alpha=1.0 e_r=1k\n\
             Xpd out0 pd 0 fc_photodetector responsivity=0.8 r_shunt=1Meg\n\
             Vb b 0 DC 1.0\n\
             Rl pd b 1k\n\
             Vsig vsig 0 DC {bias} AC 1\n"
        );
        let ac = ac_analysis(
            &parse_spice(&net).unwrap(),
            &[1e6],
            Some("vsig"),
            &DeviceRegistry::new(),
        )
        .expect("ac sweep");
        ac.magnitude("pd", 0).expect("detector node")
    };
    assert!(
        link(1.5) > 1e-3,
        "a modulator at quadrature must produce a detectable response"
    );
    assert!(
        link(0.0) < 1e-9,
        "at the transfer peak there is no intensity change to detect"
    );
}
