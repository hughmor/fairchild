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

/// The envelope group delay is a phase slope, and the slope is `−τ_g`.
///
/// `.ac` on an optical path used to be flat in phase whatever the geometry,
/// because the delay lived only in the residual and no frequency-domain
/// assembly reads one (#114). The anchor is the definition of group delay:
/// `out(Ω) = H·exp(−jΩτ_g)·in(Ω)`, so `dφ/df = −360°·τ_g` exactly, and the
/// magnitude does not move at all on a lossless segment.
///
/// The on/off pair matters as much as the slope. With the option off the
/// transmission is documented as instantaneous, so a phase slope there would be
/// a different bug — the delay leaking into a run that did not ask for it.
#[test]
fn the_envelope_group_delay_is_a_phase_slope_in_ac() {
    const L_M: f64 = 1e-2;
    const N_G: f64 = 4.2;
    let tau_g = L_M * N_G / 299_792_458.0;

    let sweep = |delay: bool| {
        let net = format!(
            "* envelope delay in .ac\n\
             .optical_port l_out\n\
             .optical_port m_out\n\
             .optical_port w_out\n\
             Xlaser l_out fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
             Xmzm l_out m_out vsig 0 fc_mzm V_pi=3.0 alpha=1.0 e_r=1k\n\
             Xwg m_out w_out fc_waveguide L_um=10000 n_g=4.2 alpha_dB_cm=0.0\n\
             Vsig vsig 0 DC 1.5 AC 1\n\
             .options waveguide_delay={}\n",
            if delay { 1 } else { 0 }
        );
        let freqs = [1e9, 3e9, 5e9];
        let ac = ac_analysis(
            &parse_spice(&net).unwrap(),
            &freqs,
            Some("vsig"),
            &DeviceRegistry::new(),
        )
        .expect("ac sweep");
        let phases: Vec<f64> = (0..freqs.len())
            .map(|i| ac.phase_deg("w_out_re_0", i).expect("optical output"))
            .collect();
        let mags: Vec<f64> = (0..freqs.len())
            .map(|i| ac.magnitude("w_out_re_0", i).expect("optical output"))
            .collect();
        (freqs, phases, mags)
    };

    let (freqs, phases, mags) = sweep(true);
    let want = -360.0 * tau_g * (freqs[1] - freqs[0]);
    for w in phases.windows(2) {
        // Unwrap: the step is under 180 degrees, so the shorter arc is the
        // real one.
        let mut d = w[1] - w[0];
        while d > 180.0 {
            d -= 360.0;
        }
        while d < -180.0 {
            d += 360.0;
        }
        assert!(
            (d - want).abs() < 1e-6,
            "phase step {d:.4} deg per {:.1} GHz, expected {want:.4} \
             (tau_g = {:.2} ps)",
            (freqs[1] - freqs[0]) / 1e9,
            tau_g * 1e12
        );
    }
    assert!(
        mags.windows(2).all(|w| (w[1] - w[0]).abs() < 1e-12),
        "a lossless delay must not change the magnitude: {mags:?}"
    );

    let (_, flat, _) = sweep(false);
    assert!(
        flat.windows(2).all(|w| (w[1] - w[0]).abs() < 1e-9),
        "with the option off the transmission is instantaneous, so the phase \
         must not move: {flat:?}"
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
