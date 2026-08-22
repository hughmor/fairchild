//! The Jacobian the devices stamp must equal `∂f/∂x`, or be declared as not.
//!
//! `dL/dp = −λᵀ·∂f/∂p` is the total derivative only if `Jᵀλ = ∂L/∂x` was solved
//! with the true `J`.  Newton has no such requirement — any contracting
//! iteration matrix converges to the same fixed point, so a device may freeze a
//! coefficient at the previous iterate and the forward answer is unaffected.
//! The gradient is not: a missing block silently contributes **zero** to every
//! path through it, which reads as a genuine insensitivity.
//!
//! So every frozen block has to be either found generically (λ wires) or
//! declared by its device (`Device::frozen_jacobian_columns`), and this file is
//! the check.  It is also the tool to reach for when adding a photonic device:
//! if a new model freezes something and forgets to say so, the gradient through
//! it is zero and nothing else in the suite notices.

use fairchild_core::adjoint::jacobian_check;
use fairchild_core::{dc_op_nr_with_registry_opts, DeviceRegistry, SimOptions};
use fairchild_parser::parse_spice;

fn opts() -> SimOptions {
    SimOptions {
        reltol: 1e-12,
        vntol: 1e-14,
        ..SimOptions::default()
    }
}

/// Assert that every disagreement between the stamped Jacobian and `∂f/∂x` sits
/// in a column the adjoint already knows to re-derive.
fn assert_complete(tag: &str, src: &str) {
    let net = parse_spice(src).unwrap_or_else(|e| panic!("{tag}: parse: {e}"));
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);
    let o = opts();
    let op = dc_op_nr_with_registry_opts(&net, &reg, &o).unwrap_or_else(|e| panic!("{tag}: {e}"));

    // 1e-6 relative.  This started at 1e-4, on the reasoning that the reference
    // side is itself a finite difference and the failure mode being hunted — a
    // wholly missing block — is off by 100 %.  That reasoning was wrong, and it
    // cost a real bug: `fc_photodetector` stamped a 1 µS shunt it then cancelled
    // out of the residual, which against the 50 Ω load below is 5e-5 relative
    // and passed.  It was still a genuine `∂f/∂x` error, and it made every
    // adjoint gradient through a detector wrong by `R_load/r_shunt`.  A stamp
    // that is right to four figures and no further is not right.
    let bad = jacobian_check(&net, &reg, &o, &op.x, 1e-6, 1e-9).unwrap();
    let undeclared: Vec<_> = bad.iter().filter(|m| !m.frozen).collect();
    assert!(
        undeclared.is_empty(),
        "{tag}: {} undeclared Jacobian mismatch(es) — a gradient through any of \
         these columns is silently zero. First few: {:?}",
        undeclared.len(),
        undeclared.iter().take(4).collect::<Vec<_>>()
    );
}

// ---------------------------------------------------------------------------
// Electrical — these have no excuse; the stamps are analytic derivatives
// ---------------------------------------------------------------------------

#[test]
fn electrical_devices_stamp_a_complete_jacobian() {
    assert_complete(
        "divider",
        "* r\nV1 in 0 DC 1\nR1 in out 1k\nR2 out 0 2k\n.op\n",
    );
    assert_complete(
        "diode",
        "* d\n.model dmod D (IS=1e-14 N=1.0)\nV1 in 0 DC 2\nR1 in mid 1k\nD1 mid 0 dmod\n.op\n",
    );
    assert_complete(
        "mosfet",
        "* m\n.model nm NMOS (VTO=0.7 KP=100u)\nVdd d 0 DC 3\nVg g 0 DC 2\nRd d dr 1k\nM1 dr g 0 0 nm W=10u L=1u\n.op\n",
    );
    assert_complete(
        "bjt",
        "* q\n.model qn NPN (IS=1e-16 BF=100)\nVcc c 0 DC 5\nVb b 0 DC 0.7\nRc c cc 1k\nQ1 cc b 0 0 qn\n.op\n",
    );
}

// ---------------------------------------------------------------------------
// Photonic — passive first, then the electro-optic families
// ---------------------------------------------------------------------------

const LASER: &str = "Xl0 in0 fc_cw_laser power_mW=1.0 wavelength_nm=1550\n";

fn optical(body: &str) -> String {
    format!(".optical_port in0\n.optical_port out0\n{LASER}{body}.op\n")
}

#[test]
fn passive_photonic_devices_stamp_a_complete_jacobian() {
    assert_complete(
        "waveguide",
        &optical("Xwg in0 out0 fc_waveguide L_um=250 n_g=4.2 alpha_dB_cm=2.0\n"),
    );
    assert_complete("grating", &optical("Xgc in0 out0 fc_grating_coupler\n"));
    assert_complete(
        "splitter",
        &format!(
            ".optical_port in0\n.optical_port o1\n.optical_port o2\n{LASER}\
             Xsp in0 o1 o2 fc_splitter\n.op\n"
        ),
    );
    assert_complete(
        "dcoupler",
        &format!(
            ".optical_port in0\n.optical_port in1\n.optical_port o1\n.optical_port o2\n{LASER}\
             Xdc in0 in1 o1 o2 fc_dcoupler kappa_L=0.336\n.op\n"
        ),
    );
}

/// The electro-optic families: every one of these freezes its optical
/// coefficient at the previous Newton iterate, so every one has to declare the
/// electrical column it reads.
#[test]
fn electro_optic_devices_declare_what_they_freeze() {
    assert_complete(
        "mzm",
        &format!(
            ".optical_port in0\n.optical_port out0\n{LASER}\
             Xmzm in0 out0 vsig 0 fc_mzm V_pi=3.0 alpha=1.0 e_r=1000\n\
             Vsig vsig 0 DC 1.5\n.op\n"
        ),
    );
    // `dn_dv` is what makes these electro-optic at all.  Without it the device
    // has no voltage dependence to freeze and the check passes vacuously —
    // which it did, on the first draft of this test.
    assert_complete(
        "pn phase shifter",
        &format!(
            ".optical_port in0\n.optical_port out0\n{LASER}\
             Xpn in0 out0 vmod 0 fc_pn_ps L_um=500 dn_dv=5e-5 g_pn=1e-3 alpha_dB_cm=10 pin_at_ref=1\n\
             Vm vmod 0 DC -1.0\n.op\n"
        ),
    );
    assert_complete(
        "pn phase shifter with junction cap",
        &format!(
            ".optical_port in0\n.optical_port out0\n{LASER}\
             Xpn in0 out0 vmod 0 fc_pn_ps_cap L_um=100 dn_dv=5e-5 g_pn=1e-3 c_j0=100f v_bi=0.7 m_j=0.5\n\
             Vm vmod 0 DC -1.0\n.op\n"
        ),
    );
    assert_complete(
        "thermal phase shifter",
        &format!(
            ".optical_port in0\n.optical_port out0\n{LASER}\
             Xth in0 out0 vh 0 fc_thermal_ps L_um=200 p_pi_th=20m r_heater=500\n\
             Vh vh 0 DC 2.0\n.op\n"
        ),
    );
    assert_complete(
        "optical 2x2 weight bank",
        ".optical_port bus\n.optical_port dark\n.optical_port thru\n.optical_port drop\n\
         Xl1 bus fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
         Xwb bus dark thru drop wctl 0 fc_optical_2x2 w=0 dw_dv_0=0.5\n\
         Vw wctl 0 DC 0.4\n.op\n",
    );
}

/// A full electro-optic link: modulator drives light, photodiode turns it back
/// into current, and the gradient has to traverse both directions.
///
/// Two load resistances, because the size of a detector-stamp error relative to
/// its own row scales with the load: the 1 µS shunt bug was 5e-5 against 50 Ω
/// and 1e-3 against 1 kΩ, and only one of those was ever going to be noticed.
#[test]
fn a_full_eo_link_stamps_a_complete_jacobian() {
    for r_load in ["50", "1k"] {
        assert_complete(
            &format!("link with a {r_load} load"),
            &format!(
                ".optical_port in0\n.optical_port out0\n{LASER}\
                 Xmzm in0 out0 vsig 0 fc_mzm V_pi=3.0 alpha=1.0 e_r=1000\n\
                 Xpd out0 pout 0 fc_photodetector responsivity=0.8\n\
                 Rl pout 0 {r_load}\n\
                 Vsig vsig 0 DC 1.5\n.op\n"
            ),
        );
    }
}
