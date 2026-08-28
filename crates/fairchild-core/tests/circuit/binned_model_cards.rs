//! A binned PDK card loads, and the geometry picks the bin.
//!
//! This was the last thing standing between fairchild and a production PDK
//! (#77 §1). `.model nch.1 nmos (LMIN=… )` used to register a model *called*
//! `nch.1`, so the element line asking for `nch` failed as `unknown model`.
//! Essentially every foundry PDK bins by W/L, so this blocked loading on its own.
//!
//! # What each test would catch
//!
//! The interesting failures here are not "it errored". They are:
//!
//! * the deck loads and every instance quietly gets the *same* bin;
//! * geometry outside every window silently takes the nearest bin;
//! * the four selectors reach the device as unknown parameters;
//! * a dotted model name that was never meant as a bin gets split.
//!
//! So each test compares an *observable* — drain current, which is
//! monotonic in `VTO` — rather than asserting that a solve happened.

use fairchild_core::models::mosfet1::Mosfet1;
use fairchild_core::{dc_op_nr, DeviceRegistry, SimError};
use fairchild_parser::parse_spice;

/// Two bins that differ only in `VTO`, so the selected bin is readable off the
/// drain current. `VTO` is the right knob: a wrong bin cannot hide in it.
const PDK: &str = "\
* a binned PDK card, cut down to the two things that matter
.model nch.1 nmos (LMIN=0.18u LMAX=0.30u WMIN=0.22u WMAX=1u  VTO=0.40 KP=200u)
.model nch.2 nmos (LMIN=0.30u LMAX=1.00u WMIN=0.22u WMAX=1u  VTO=0.70 KP=200u)
VDD d 0 DC 1.8
VG  g 0 DC 1.2
";

fn id_for(l_um: f64, w_um: f64) -> Result<f64, SimError> {
    let deck = format!("{PDK}M1 d g 0 0 nch W={w_um}u L={l_um}u\n.op\n");
    let net = parse_spice(&deck).expect("parse");
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);
    let r = fairchild_core::dc_op_nr_with_registry(&net, &reg)?;
    Ok(r.vsrc_current("vdd").expect("I(vdd)").abs())
}

/// The gate: an instance names the *base* model and gets the bin its geometry
/// falls in.
///
/// # Why the geometry is chosen the way it is
///
/// Both instances have **W/L = 1**. A Level-1 saturation current is
/// `0.5·KP·(W/L)·(Vgs−Vth)²`, so with W/L held equal the only thing left that can
/// move the current is `VTO` — which is to say, the bin. The first version of
/// this test compared `L=0.25u` against `L=0.5u` at fixed `W`, and so measured
/// the geometry change rather than the bin change: sabotaging `select` to return
/// bin one for every geometry left it green.
///
/// The expected ratio is a closed form, not a previous run:
/// `((1.2−0.40)/(1.2−0.70))² = 2.56`.
#[test]
fn geometry_selects_the_bin_and_not_merely_the_geometry() {
    let bin1 = id_for(0.25, 0.25).expect("bin 1 must solve");
    let bin2 = id_for(0.50, 0.50).expect("bin 2 must solve");
    let ratio = bin1 / bin2;
    let want = ((1.2 - 0.40) / (1.2 - 0.70_f64)).powi(2);
    assert!(
        (ratio - want).abs() / want < 0.02,
        "W/L is 1 in both, so I(bin1)/I(bin2) is fixed by VTO alone and must be \
         {want:.4} — got {ratio:.4} (bin1={bin1:.6e}, bin2={bin2:.6e}). A ratio of \
         1 means both instances got the same card."
    );
}

/// Which bin, not merely a different one: the selected card must give the same
/// answer as a deck in which that card's parameters are the whole model.
///
/// The ratio test above pins that the two bins differ. This pins that bin 2 is
/// *bin 2*, against an anchor outside the binning code — swapping the two cards'
/// parameters would satisfy the ratio and fail this.
#[test]
fn the_selected_bin_matches_that_card_used_alone() {
    for (l_um, vto) in [(0.25, 0.40), (0.50, 0.70)] {
        let binned = id_for(l_um, l_um).expect("binned deck must solve");
        let alone = format!(
            "* one card, no bins\n.model nch nmos (VTO={vto} KP=200u)\n\
             VDD d 0 DC 1.8\nVG g 0 DC 1.2\n\
             M1 d g 0 0 nch W={l_um}u L={l_um}u\n.op\n"
        );
        let net = parse_spice(&alone).expect("parse");
        let unbinned = dc_op_nr(&net)
            .expect("solve")
            .vsrc_current("vdd")
            .expect("I(vdd)")
            .abs();
        let rel = (binned - unbinned).abs() / unbinned;
        assert!(
            rel < 1e-9,
            "L={l_um}u selected a bin whose current is {binned:.9e}, but VTO={vto} \
             used as the whole model gives {unbinned:.9e} (rel {rel:.2e}). The \
             geometry picked a card, and it was not this one."
        );
    }
}

/// Outside every window is a hard error naming the geometry and the windows.
/// The alternative — nearest bin — is a wrong answer with nothing to read.
#[test]
fn geometry_outside_every_bin_is_refused_by_name() {
    let err = id_for(5.0, 0.5).expect_err("L=5um is past every LMAX");
    let msg = err.to_string();
    assert!(
        matches!(err, SimError::NoMatchingBin { .. }),
        "should be NoMatchingBin, not UnknownModel: {msg}"
    );
    for needle in ["nch", "nch.1", "nch.2"] {
        assert!(msg.contains(needle), "error must name {needle}: {msg}");
    }
    // And it must not read as a missing model, which is what the old code said.
    assert!(
        !msg.contains("unknown model"),
        "a binned card that exists must not report as missing: {msg}"
    );
}

/// The selectors choose a bin. They are not Level 1 parameters, so a binned deck
/// must not warn about four unknown parameters per card — a PDK has hundreds of
/// cards and that is the difference between a usable log and an unreadable one.
#[test]
fn the_geometry_selectors_are_not_reported_as_unknown_parameters() {
    let net = parse_spice(&format!("{PDK}M1 d g 0 0 nch W=0.5u L=0.25u\n.op\n")).expect("parse");
    // `classify` strips the selectors, and the device is what decides what it
    // does not recognise — so ask both, the way registration does.
    for card in &net.models {
        let (_, bin) = fairchild_core::binning::classify(&card.name, &card.params);
        let (_, unknown) = Mosfet1::from_model_params(false, &bin.params);
        assert!(
            unknown.is_empty(),
            "card '{}' left {unknown:?} unrecognised — the four geometry \
             selectors choose a bin and are not model parameters, and a PDK has \
             hundreds of cards, so four spurious warnings each is the difference \
             between a usable log and an unreadable one",
            card.name
        );
    }
}

/// A dot in a model name is legal and predates binning. Reinterpreting
/// `my.model` as bin `model` of `my` would break a deck that never asked.
#[test]
fn a_dotted_name_without_selectors_stays_one_model() {
    let deck = "* dotted, not binned\n\
                .model my.nmos nmos (VTO=0.4 KP=200u)\n\
                VDD d 0 DC 1.8\nVG g 0 DC 1.2\n\
                M1 d g 0 0 my.nmos W=0.5u L=0.25u\n.op\n";
    let net = parse_spice(deck).expect("parse");
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);
    let r = fairchild_core::dc_op_nr_with_registry(&net, &reg)
        .expect("a dotted model name must still resolve");
    assert!(
        r.vsrc_current("vdd").expect("I(vdd)").abs() > 0.0,
        "the device conducts, so the card reached it"
    );
}

/// Spectre writes binning as a braced body of numbered sections. It must reach
/// the same selection, or the two dialects disagree about the same PDK.
#[test]
fn spectre_braced_bins_reach_the_same_selection() {
    let spectre = "\
simulator lang=spectre
model nch bsim4 {
  1: lmin=0.18u lmax=0.30u wmin=0.22u wmax=1u vto=0.40 kp=200u
  2: lmin=0.30u lmax=1.00u wmin=0.22u wmax=1u vto=0.70 kp=200u
}
";
    let spice = fairchild_parser::spectre::to_spice(spectre).expect("translate");
    assert!(
        spice.contains(".model nch.1") && spice.contains(".model nch.2"),
        "each section becomes its own card: {spice}"
    );
    assert!(
        spice.contains("lmin=0.18u") || spice.contains("lmin=1.8e-7"),
        "the window survives the translation: {spice}"
    );
}

/// A braced body that is not numbered sections is still refused. Guessing which
/// card an instance gets is the wrong answer this whole path exists to avoid.
#[test]
fn an_unreadable_braced_body_is_still_refused() {
    let spectre = "simulator lang=spectre\nmodel nch bsim4 { vto=0.4 kp=200u }\n";
    let err = fairchild_parser::spectre::to_spice(spectre)
        .expect_err("a braced body with no sections cannot be read as bins");
    let msg = err.to_string();
    assert!(
        msg.contains("numbered sections"),
        "the refusal should say what it wanted: {msg}"
    );
}

/// Unbinned cards must keep working, and by the same code path — an unbinned
/// card is a group of one with an unbounded window. A regression here would show
/// up as every non-PDK deck in the suite failing, but this states the intent.
#[test]
fn an_unbinned_card_still_resolves() {
    let deck = "* plain\n.model nm nmos (VTO=0.4 KP=200u)\n\
                VDD d 0 DC 1.8\nVG g 0 DC 1.2\n\
                M1 d g 0 0 nm W=10u L=1u\n.op\n";
    let net = parse_spice(deck).expect("parse");
    let r = dc_op_nr(&net).expect("an unbinned card must solve");
    assert!(r.vsrc_current("vdd").expect("I(vdd)").abs() > 0.0);
}
