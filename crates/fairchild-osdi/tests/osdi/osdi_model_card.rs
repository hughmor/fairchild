//! `.model <card> <module> (params)` must reach an OSDI device.
//!
//! This is the idiom every foundry PDK ships — `.osdi bsim4.osdi` + a `.model`
//! card of a few hundred parameters + `M1 d g s b <card> W= L=` — and it is
//! what the user guide documents. It used to fail outright with
//! "unknown model", because `.osdi` registered descriptors under the *module*
//! name and nothing ever bound a card name to one.
//!
//! Pre-condition: legacy/va-models/build/diode_shockley.osdi must exist.
//! Build it with:
//!   DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib \
//!   openvaf-r legacy/va-models/electronic/diode_shockley.va \
//!     --output legacy/va-models/build/diode_shockley.osdi

use std::path::PathBuf;
use std::sync::Arc;

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_osdi::OsdiLibrary;
use fairchild_parser::parse_spice;

/// 1 mA into a diode: V = N·Vt·ln(I/Is), so Is and N are both visible in V(b).
fn v_at(diode_lines: &str) -> Option<f64> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../legacy/va-models/build/diode_shockley.osdi");
    if !path.exists() {
        eprintln!("Skipping: {} not found — compile it first.", path.display());
        return None;
    }

    let netlist = parse_spice(&format!("* card test\nIb 0 b 1m\n{diode_lines}.op\n")).unwrap();
    let lib = Arc::new(unsafe { OsdiLibrary::open(&path) }.expect("dlopen failed"));
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&netlist.models);
    lib.register_into(&mut registry);
    registry.register_loaded_model_cards(&netlist.models);

    let result = dc_op_nr_with_registry(&netlist, &registry).expect("DC OP failed");
    Some(result.node_voltage("b").unwrap())
}

#[test]
fn model_card_matches_the_same_params_inline() {
    let Some(inline) = v_at("Xd b 0 diode_shockley Is=1e-12 N=1.4\n") else {
        return;
    };
    let card = v_at("Xd b 0 mydiode\n.model mydiode diode_shockley (Is=1e-12 N=1.4)\n").unwrap();
    assert!(
        (inline - card).abs() < 1e-12,
        "card {card:.9} != inline {inline:.9}"
    );

    // …and the params genuinely bit, or the comparison above is vacuous.
    let defaults = v_at("Xd b 0 diode_shockley\n").unwrap();
    assert!(
        (defaults - card).abs() > 0.1,
        "Is/N had no effect: defaults {defaults:.6} vs card {card:.6}"
    );
}

#[test]
fn instance_param_overrides_the_card() {
    let Some(overridden) =
        v_at("Xd b 0 mydiode N=1.4\n.model mydiode diode_shockley (Is=1e-12 N=1.0)\n")
    else {
        return;
    };
    let expected = v_at("Xd b 0 diode_shockley Is=1e-12 N=1.4\n").unwrap();
    assert!(
        (overridden - expected).abs() < 1e-12,
        "instance N=1.4 did not win over the card's N=1.0: {overridden:.9} vs {expected:.9}"
    );
}

/// A two-terminal OSDI model instantiated as `D` must behave exactly as the
/// same model instantiated as `X`. The `D` parser used to drop instance
/// params, so the two forms silently disagreed.
#[test]
fn d_form_and_x_form_agree() {
    let Some(d_form) = v_at("D1 b 0 diode_shockley Is=1e-12 N=1.4\n") else {
        return;
    };
    let x_form = v_at("Xd b 0 diode_shockley Is=1e-12 N=1.4\n").unwrap();
    assert!(
        (d_form - x_form).abs() < 1e-12,
        "D-form {d_form:.9} != X-form {x_form:.9}"
    );
}

#[test]
fn native_cards_keep_their_own_handler() {
    // A `.model … d(…)` card must still go through register_builtin_diodes,
    // which does construction-time work register_loaded_model_cards cannot.
    // The OSDI library is loaded too, so both passes run over the same cards.
    let Some(native) = v_at("D1 b 0 dnat\n.model dnat d(is=1e-12 n=1.4)\n") else {
        return;
    };
    let osdi = v_at("Xd b 0 diode_shockley Is=1e-12 N=1.4\n").unwrap();
    assert!(
        (native - osdi).abs() < 5e-3,
        "native diode {native:.6} and OSDI diode {osdi:.6} disagree"
    );
}
