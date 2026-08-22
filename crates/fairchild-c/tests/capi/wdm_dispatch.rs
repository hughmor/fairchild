//! WDM dispatch through the C ABI (#52).
//!
//! `fc_load_file` / `fc_load_string` used to call the no-oracle parse, so the C
//! ABI was the one front end still deciding bundle arity from the parser's
//! static name list. That list is keyed on the name written on the X-line, so a
//! `.model`-card-named instance could never be found in it, and every deck
//! below was refused before a device was built.
//!
//! Nothing loaded a bundle deck through these entry points, which is why the
//! gap outlived the fix everywhere else. This file is that missing coverage.

use std::ffi::{CStr, CString};

use fairchild_c::{
    fc_error, fc_load_string, fc_op_node, fc_run_op, fc_sim_free, fc_sim_new, FC_OK,
};

/// Two channels at *different* powers through a card-named waveguide. The
/// powers differ so that taking channel 0's answer and applying it bus-wide —
/// the shape of the original bug — cannot pass.
const CARD_WAVEGUIDE: &str = "\
.optical_port c0
.optical_port c1
.optical_port bus 2
.optical_port out 2
Xl0 c0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xl1 c1 fc_cw_laser power_mW=2.5 wavelength_nm=1551
Xmux bus c0 c1 fc_mux
.model mywg fc_waveguide L_um=1000 alpha_dB_cm=3.0 n_g=4.2
Xw bus out mywg
.op
";

const L_CM: f64 = 0.1;
const ALPHA_DB_CM: f64 = 3.0;

struct Sim(*mut fairchild_c::FcSim);

impl Sim {
    fn new() -> Self {
        let p = fc_sim_new();
        assert!(!p.is_null(), "fc_sim_new returned NULL");
        Sim(p)
    }

    fn load(&self, deck: &str) -> i32 {
        let c = CString::new(deck).unwrap();
        unsafe { fc_load_string(self.0, c.as_ptr()) }
    }

    fn err(&self) -> String {
        let p = unsafe { fc_error(self.0) };
        if p.is_null() {
            String::new()
        } else {
            unsafe { CStr::from_ptr(p) }.to_string_lossy().into_owned()
        }
    }

    fn node(&self, name: &str) -> f64 {
        let c = CString::new(name).unwrap();
        let mut v = 0.0;
        let rc = unsafe { fc_op_node(self.0, c.as_ptr(), &mut v) };
        assert_eq!(rc, FC_OK, "fc_op_node({name}) failed: {}", self.err());
        v
    }

    /// Optical power on one channel of a bundle, in the field's own units.
    fn power(&self, port: &str, ch: usize) -> f64 {
        let re = self.node(&format!("{port}_re_{ch}"));
        let im = self.node(&format!("{port}_im_{ch}"));
        re * re + im * im
    }
}

impl Drop for Sim {
    fn drop(&mut self) {
        unsafe { fc_sim_free(self.0) }
    }
}

/// The absolute anchor: a waveguide's power transmission is `10^(-αL/10)`, and
/// it applies to each channel separately. Checked as a *ratio* across the
/// device, so whatever the mux does to the input does not enter the assertion —
/// and checked on both channels, which is what a bus-wide broadcast of channel
/// 0's result fails.
#[test]
fn a_card_named_photonic_device_carries_a_bundle_through_the_c_abi() {
    let sim = Sim::new();
    assert_eq!(
        sim.load(CARD_WAVEGUIDE),
        FC_OK,
        "a carded fc_waveguide on a 2-channel bundle must load: {}",
        sim.err()
    );
    assert_eq!(
        unsafe { fc_run_op(sim.0) },
        FC_OK,
        "the op should solve: {}",
        sim.err()
    );

    let expected = 10f64.powf(-ALPHA_DB_CM * L_CM / 10.0);
    let (p_in_0, p_in_1) = (sim.power("bus", 0), sim.power("bus", 1));

    // A degenerate deck would make the per-channel check vacuous.
    assert!(
        (p_in_1 / p_in_0 - 2.5).abs() < 1e-6,
        "the two channels must arrive at different powers for this test to mean \
         anything, got {p_in_0} and {p_in_1}"
    );

    for ch in 0..2 {
        let (p_in, p_out) = (sim.power("bus", ch), sim.power("out", ch));
        let t = p_out / p_in;
        assert!(
            (t - expected).abs() < 1e-9,
            "channel {ch}: transmission {t} should be {expected} \
             (α={ALPHA_DB_CM} dB/cm over {L_CM} cm); in={p_in} out={p_out}"
        );
    }
}

/// Dispatch must not have become a free-for-all: one laser still cannot serve
/// four channels, and the refusal still reaches the caller through `fc_error`
/// with the fix in it.
#[test]
fn a_scalar_device_on_a_wide_bundle_is_still_refused_and_says_why() {
    let sim = Sim::new();
    let rc = sim.load(
        "\
.optical_port bus 4
Xl bus fc_cw_laser power_mW=1 wavelength_nm=1550
.op
",
    );
    assert_ne!(rc, FC_OK, "one laser cannot serve 4 channels");
    let msg = sim.err();
    assert!(
        msg.contains("no WDM semantics"),
        "the refusal should still explain itself, got: {msg}"
    );
}
