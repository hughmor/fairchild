//! Every parameter a deck writes reaches the device before it decides anything.
//!
//! #31 replaced four spellings of "set the device up, then tell it its
//! parameters" with one, and made construction fallible so a device can refuse a
//! configuration instead of asserting from inside `eval`. The named risk of that
//! refactor is precise: *a device silently receiving different parameters than
//! before*. A device that gets a default where the deck gave a value usually
//! still produces a plausible number, and no existing test would notice.
//!
//! So these are the assertions that would catch it. Each one is a closed form or
//! an exact equality, and each was checked by breaking the thing it covers.
//!
//! The `.model`-card cases matter most: card parameters used to be written onto
//! the device *after* the factory returned, and are now merged into the
//! `ParamSet` before construction. That is the change with the most room to drop
//! a value on the floor.

use fairchild_core::{dc_op_nr, dc_op_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

/// Reflected power off a one-port facet, in watts. `re = 1 V` is 1 W, so a
/// launched field of 1 makes the reflected power numerically equal to `R`.
///
/// `fc_facet` is the probe of choice because its parameter *is* the answer:
/// `P_out = R · P_in` exactly, with no fitting constant in between, so a
/// parameter that failed to arrive shows up as a different number rather than a
/// slightly different one.
fn reflected_power(deck_body: &str) -> f64 {
    let src = format!(
        ".options enable_bidirectional=1\n{deck_body}\
         V_re p_re_fw 0 DC 1\nV_im p_im_fw 0 DC 0\n\
         V_wl p_wl 0 DC 1.55e-6\n.op\n"
    );
    let net = parse_spice(&src).expect("parse");
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&net.models);
    registry.register_loaded_model_cards(&net.models);
    let r = dc_op_nr_with_registry(&net, &registry).expect("solve");
    let re = r.node_voltage("p_re_bw").unwrap();
    let im = r.node_voltage("p_im_bw").unwrap();
    re * re + im * im
}

const PORT: &str = "p_re_fw p_im_fw p_re_bw p_im_bw p_wl";

/// The baseline: a parameter on the element line arrives.
#[test]
fn an_instance_parameter_reaches_the_device() {
    let p = reflected_power(&format!("Xf {PORT} fc_facet reflectance=0.25\n"));
    assert!((p - 0.25).abs() < 1e-15, "{p:.12e}");
}

/// A `.model` card parameter arrives, with nothing on the element line.
///
/// This is the case the refactor moved. Card parameters used to be applied by a
/// wrapper closure *after* the target factory had already built, set up, and (as
/// of #31) validated the device.
#[test]
fn a_model_card_parameter_reaches_the_device() {
    let p = reflected_power(&format!(
        ".model mirror fc_facet (reflectance=0.25)\nXf {PORT} mirror\n"
    ));
    assert!((p - 0.25).abs() < 1e-15, "{p:.12e}");
}

/// Both present: the element line wins, and it wins *completely* — the card's
/// value must not be blended, averaged, or applied afterwards on top.
#[test]
fn the_element_line_beats_the_card() {
    let p = reflected_power(&format!(
        ".model mirror fc_facet (reflectance=0.25)\nXf {PORT} mirror reflectance=0.64\n"
    ));
    assert!((p - 0.64).abs() < 1e-15, "{p:.12e}");
    // …and the two values are far enough apart that the wrong one cannot be
    // mistaken for a tolerance: 0.25 and 0.64 differ by more than a factor of 2.
}

/// A card parameter and an element parameter that are *different parameters* both
/// arrive. Merging one set into another is where a key gets lost.
#[test]
fn a_card_and_an_element_parameter_both_arrive() {
    // Card says transmittance, element says loss; R is the remainder, 0.3.
    let p = reflected_power(&format!(
        ".model part fc_facet (transmittance=0.5)\nXf {PORT} part loss=0.2\n"
    ));
    assert!((p - 0.3).abs() < 1e-15, "{p:.12e}");
}

/// Construction refuses a configuration only the two halves together make
/// illegal — the card is fine on its own, the element line is fine on its own.
///
/// This is what fallible construction buys: neither `setup_instance` (which has
/// no parameters) nor `set_real_param` (which has one at a time) can see this.
#[test]
fn a_budget_broken_only_by_the_combination_is_refused() {
    let src = format!(
        ".options enable_bidirectional=1\n\
         .model part fc_facet (transmittance=0.6)\n\
         Xf {PORT} part reflectance=0.7\n\
         V_re p_re_fw 0 DC 1\nV_im p_im_fw 0 DC 0\n\
         V_wl p_wl 0 DC 1.55e-6\n.op\n"
    );
    let net = parse_spice(&src).expect("parse");
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&net.models);
    registry.register_loaded_model_cards(&net.models);
    let e = match dc_op_nr_with_registry(&net, &registry) {
        Ok(_) => panic!("0.7 + 0.6 > 1 must be refused"),
        Err(e) => e.to_string(),
    };
    assert!(e.contains("no power left"), "{e}");
    // The element is named. A device's own `Err` cannot know this; the caller
    // adds it, and that is the only reason the message is usable in a deck with
    // more than one facet.
    assert!(e.contains("xf"), "the refusal must name the element: {e}");
}

/// The same, for the parameter path that has nothing to do with photonics: a
/// diode's `.model` card and its instance `AREA`.
///
/// Included because `fc_facet` reaches the registry through
/// `register_native_photonics` and a diode through `register_builtin_diodes` —
/// two different registration sites, and #31 rewrote both.
#[test]
fn a_diode_card_and_its_instance_area_both_arrive() {
    let one = current("* base\n.model dm D (IS=1e-14 N=1)\nV1 a 0 DC 0.7\nD1 a 0 dm\n.op\n");
    let two = current("* area\n.model dm D (IS=1e-14 N=1)\nV1 a 0 DC 0.7\nD1 a 0 dm area=2\n.op\n");
    // AREA=2 is exactly two of them in parallel.
    assert!(
        ((two / one) - 2.0).abs() < 1e-7,
        "area=2 gave {:.9} of the unit current",
        two / one
    );
    // And the card's own IS reaches the device: a decade of IS is a decade of
    // current at fixed bias.
    let ten = current("* IS×10\n.model dm D (IS=1e-13 N=1)\nV1 a 0 DC 0.7\nD1 a 0 dm\n.op\n");
    assert!(
        ((ten / one) - 10.0).abs() < 1e-5,
        "IS×10 gave {:.9}× the current",
        ten / one
    );
}

fn current(src: &str) -> f64 {
    let net = parse_spice(src).expect("parse");
    dc_op_nr(&net)
        .expect("solve")
        .vsrc_current("v1")
        .expect("source current")
        .abs()
}
