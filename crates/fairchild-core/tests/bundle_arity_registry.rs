//! WDM dispatch is the registry's answer, not a name list in the parser (#52).
//!
//! `bundle_arity_for` could only ever match a name written literally on an
//! X-line. A `.model`-card-named instance is looked up under the *card's* name,
//! so no card-based device was ever found there — and every photonic device
//! whose only route is a card (`fc_awgr` in table mode, `fc_phase_shifter_expr`)
//! was therefore refused on the very bundles it exists to carry, as was the
//! whole documented `LEVEL` idiom.
//!
//! These decks are the ones that used to be rejected. Each is a parse-level
//! check: the failure mode was a `ParseError`, before any device was built.

use fairchild_core::DeviceRegistry;
use fairchild_parser::{parse_spice_with_arity, PermissiveArity};

/// Parse the way a real load does: permissive pass to harvest cards, build the
/// registry, then parse again with the registry deciding dispatch.
fn two_pass(src: &str) -> Result<fairchild_parser::Netlist, fairchild_parser::ParseError> {
    let probe = parse_spice_with_arity(src, &PermissiveArity)?;
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&probe.models);
    parse_spice_with_arity(src, &reg)
}

const AWGR_TABLE: &str = "\
.optical_port i0 2
.optical_port i1 2
.optical_port o0 2
.optical_port o1 2
.model awg2 fc_awgr fwhm_ghz=40
Xr i0 i1 o0 o1 awg2
.op
.end
";

const LEVEL4_ON_A_BUS: &str = "\
.optical_port c0
.optical_port c1
.optical_port bus 2
.optical_port out 2
.model myfull fc_pn_th_ps LEVEL=4
Xl0 c0 fc_cw_laser power_mW=1 wavelength_nm=1550
Xl1 c1 fc_cw_laser power_mW=1 wavelength_nm=1551
Xmux bus c0 c1 fc_mux
Xpn bus out a 0 hp 0 myfull
Vb a 0 DC -1.0
Vh hp 0 DC 0
.op
.end
";

const EXPR_ON_A_BUS: &str = "\
.optical_port a 2
.optical_port b 2
.model myexpr fc_phase_shifter_expr dneff=\"1e-4*v\"
Xp a b c 0 myexpr
Vb c 0 DC 1
.op
.end
";

const CARD_WAVEGUIDE: &str = "\
.optical_port a 2
.optical_port b 2
.model mywg fc_waveguide L_um=100
Xw a b mywg
.op
.end
";

/// A card inherits its kind's dispatch. Each of these is refused outright by
/// the static table, because the name on the X-line is the card's.
#[test]
fn a_card_named_photonic_device_can_carry_a_bundle() {
    for (label, src) in [
        ("fc_awgr in table mode", AWGR_TABLE),
        (".model … fc_pn_th_ps LEVEL=4", LEVEL4_ON_A_BUS),
        ("fc_phase_shifter_expr", EXPR_ON_A_BUS),
        ("a carded fc_waveguide", CARD_WAVEGUIDE),
    ] {
        two_pass(src).unwrap_or_else(|e| panic!("{label} should parse on a 2-channel bundle: {e}"));
    }
}

/// The same decks under the static table, which is what the registry replaces
/// as the authority. If this ever starts passing, the fix above has stopped
/// being load-bearing and these tests no longer prove anything.
#[test]
fn the_static_table_alone_still_refuses_every_one_of_them() {
    for (label, src) in [
        ("fc_awgr in table mode", AWGR_TABLE),
        (".model … fc_pn_th_ps LEVEL=4", LEVEL4_ON_A_BUS),
        ("fc_phase_shifter_expr", EXPR_ON_A_BUS),
        ("a carded fc_waveguide", CARD_WAVEGUIDE),
    ] {
        assert!(
            fairchild_parser::parse_spice(src).is_err(),
            "{label}: the static table is expected to refuse this — that is the \
             bug the registry oracle exists to fix"
        );
    }
}

/// Dispatch must not become a free-for-all: a genuinely single-channel device
/// is still refused on a wide bundle, and the message still names the fix.
#[test]
fn a_scalar_device_is_still_refused_on_a_wide_bundle() {
    let src = "\
.optical_port bus 4
Xl bus fc_cw_laser power_mW=1 wavelength_nm=1550
.op
.end
";
    let err = two_pass(src).expect_err("one laser cannot serve 4 channels");
    let msg = format!("{err}");
    assert!(
        msg.contains("no WDM semantics"),
        "the refusal should still explain itself, got: {msg}"
    );
}
