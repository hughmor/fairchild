//! Resolving λ from declared routing must agree with solving for it.
//!
//! The oracle is the solver itself: λ is currently an MNA unknown, so a deck's
//! solved λ values are the answer that resolution has to reproduce. Getting them
//! to agree on real decks is what makes it safe to stop solving for λ at all —
//! and this test keeps working afterwards, because the solved side is what is
//! being replaced and the declared side is what replaces it.
//!
//! Routing is declared per device (`Device::lambda_routing`) rather than read
//! off the matrix. An earlier attempt inferred it from the assembled Jacobian
//! and could not survive removing the rows it inferred from — and it was not
//! even where the constraint lives, since `V(out_λ) − V(in_λ) = 0` is stamped
//! into a branch row while the λ node rows carry only KCL over those branches.

use fairchild_core::{dc_op_nr_with_registry, lambda, DeviceRegistry};
use fairchild_parser::parse_spice;

/// Every λ net that resolution reached must match the solved value.
fn assert_resolution_matches_the_solve(label: &str, deck: &str) -> usize {
    let net = parse_spice(deck).expect("deck parses");
    let reg = DeviceRegistry::new();
    let r = dc_op_nr_with_registry(&net, &reg).expect("DC OP converges");

    let opts = fairchild_core::SimOptions::default();
    let ctx = opts.sim_context();
    let map = lambda::resolve(&net, &ctx, &reg);

    let mut checked = 0;
    for name in net.optical_nets.iter() {
        if !fairchild_parser::is_lambda_wire(name) {
            continue;
        }
        let Ok(solved) = r.node_voltage(name) else {
            continue;
        };
        let Some(resolved) = map.get(name) else {
            panic!("{label}: {name} has no resolved λ at all");
        };
        // An undriven wire is solved as 0 and resolved to the band centre; that
        // is the documented bootstrap, not a disagreement.
        if solved == 0.0 {
            continue;
        }
        checked += 1;
        assert!(
            (solved - resolved).abs() < 1e-15,
            "{label}: {name} solved to {solved:e} m but resolution says \
             {resolved:e} m — declared routing disagrees with the matrix"
        );
    }
    assert!(checked > 0, "{label}: nothing was compared");
    checked
}

#[test]
fn a_chain_of_segments_carries_its_source_wavelength() {
    let n = assert_resolution_matches_the_solve(
        "waveguide chain",
        "\
.optical_port a
.optical_port b
.optical_port c
Xl a fc_cw_laser power_mW=2 wavelength_nm=1531.5
Xw1 a b fc_waveguide L_um=100 n_g=4.2
Xw2 b c fc_waveguide L_um=250 n_g=4.2
.op
",
    );
    assert!(n >= 3, "expected a λ per port, compared {n}");
}

/// An active device between two passives: the same declaration path, but through
/// `ActiveOpticalDevice`, where the optical terminals are followed by electrical
/// ones and the routing indices must still line up.
#[test]
fn an_active_device_passes_the_label_through_its_electrical_terminals() {
    assert_resolution_matches_the_solve(
        "phase shifter in a chain",
        "\
.optical_port a
.optical_port b
.optical_port c
Xl a fc_cw_laser power_mW=1 wavelength_nm=1310.0
Xps a b d 0 fc_pn_ps L_um=200
Vb d 0 DC -1.0
Xw b c fc_waveguide L_um=50 n_g=4.2
.op
",
    );
}

/// Two different wavelengths in one deck, so a single global answer cannot pass.
#[test]
fn two_sources_at_different_wavelengths_stay_apart() {
    assert_resolution_matches_the_solve(
        "two independent paths",
        "\
.optical_port a1
.optical_port b1
.optical_port a2
.optical_port b2
Xl1 a1 fc_cw_laser power_mW=1 wavelength_nm=1270.0
Xl2 a2 fc_cw_laser power_mW=1 wavelength_nm=1610.0
Xw1 a1 b1 fc_waveguide L_um=100 n_g=4.2
Xw2 a2 b2 fc_waveguide L_um=100 n_g=4.2
.op
",
    );
}

// ── the three devices that are not slot-for-slot ───────────────────────────
//
// `OpticalSegment`'s declaration covers everything that passes a label straight
// through, channel k in to channel k out. These three do not, so each declares
// its own — and each is a different kind of not-straight-through, which is why
// all three are here rather than one standing in for the others.

/// A mux moves port `k`'s label onto bus slot `k`; a demux moves it back. The
/// wavelengths are deliberately unequal, so a routing that crossed two channels
/// would put the wrong number on the wrong wire rather than being invisible.
#[test]
fn a_mux_and_demux_move_labels_between_ports_and_slots() {
    assert_resolution_matches_the_solve(
        "mux then demux",
        "\
.optical_port c0
.optical_port c1
.optical_port c2
.optical_port bus 3
.optical_port mid 3
.optical_port d0
.optical_port d1
.optical_port d2
Xl0 c0 fc_cw_laser power_mW=1 wavelength_nm=1546.12
Xl1 c1 fc_cw_laser power_mW=2 wavelength_nm=1550.00
Xl2 c2 fc_cw_laser power_mW=3 wavelength_nm=1553.88
Xmux bus c0 c1 c2 fc_mux
Xwg bus mid fc_waveguide L_um=500 n_g=4.2
Xdemux mid d0 d1 d2 fc_demux
.op
",
    );
}

/// The router, and the one that corrected my assumption. A label does NOT
/// follow the field's cyclic permutation: every output slot `k` mirrors one
/// chosen input port's slot `k` tag (`lambda_src`, port 0 by default), because a
/// slot *is* a wavelength across the whole router. Which port a photon arrived
/// on decides where its energy goes, not what colour it is.
///
/// Sabotage-checked against the permutation reading, which puts channel 0's
/// 1546.12 nm on the wrong output and shows up as an unreached net.
#[test]
fn an_awgr_output_mirrors_one_input_ports_tag_not_the_field_route() {
    assert_resolution_matches_the_solve(
        "2x2 awgr",
        "\
.optical_port a0
.optical_port a1
.optical_port i0 2
.optical_port i1 2
.optical_port o0 2
.optical_port o1 2
Xl0 a0 fc_cw_laser power_mW=1 wavelength_nm=1546.12
Xl1 a1 fc_cw_laser power_mW=1 wavelength_nm=1550.00
Xmux i0 a0 a1 fc_mux
Xr i0 i1 o0 o1 fc_awgr
.op
",
    );
}

/// The giona bus: eight lasers muxed onto one bundle through a ring bank. The
/// shape the whole exercise is for, and the one where an eight-way routing
/// mistake would be least visible.
#[test]
fn an_eight_channel_bus_through_a_ring_resolves_every_label() {
    let mut deck = String::new();
    for k in 0..8 {
        deck += &format!(".optical_port c{k}\n");
    }
    deck += ".optical_port bus 8\n.optical_port out 8\n";
    for k in 0..8 {
        deck += &format!(
            "Xl{k} c{k} fc_cw_laser power_mW=4 wavelength_nm={:.2}\n",
            1546.12 + 0.8 * k as f64
        );
    }
    deck += "Xmux bus";
    for k in 0..8 {
        deck += &format!(" c{k}");
    }
    deck += " fc_mux\nXpn bus out a 0 fc_pn_ps_full L_um=50\nVb a 0 DC 0\n.op\n";
    assert_resolution_matches_the_solve("giona-shaped bus", &deck);
}

// ── a model that declares nothing must say so ──────────────────────────────
//
// λ is resolved before the solve from what devices declare, so a model that
// declares nothing is invisible to resolution: everything downstream of it takes
// the band centre instead of the wavelength actually present. For a ring or a
// filter that means evaluating a passband at a colour nowhere in the circuit.
//
// A hand-written fixed-port Verilog-A model is exactly this case and cannot fix
// itself — its routing is inside compiled Verilog-A that fairchild does not
// parse, whereas a bundle-dialect model escapes only because we generate its
// header and therefore know its layout. Guessing from port names was ruled out:
// pairing λ ports in declaration order is right until it silently is not, and a
// PCell that hand-wires its bundle (`source_bank.sp` names them `a1w`) defeats
// name-based reasoning outright.
//
// So it warns. It cannot error without refusing every optical Verilog-A model
// written before the dialect existed, and a model whose λ ports pass a tag to
// nothing downstream is harmless. These two tests pin *when* it speaks, because
// a diagnostic that fires on the harmless case trains people to ignore it.

/// Resolution over a deck, returning the λ nets it could not reach.
fn unreached_nets(deck: &str) -> Vec<String> {
    let net = parse_spice(deck).expect("deck parses");
    let reg = DeviceRegistry::new();
    let opts = fairchild_core::SimOptions::default();
    let map = lambda::resolve(&net, &opts.sim_context(), &reg);
    map.unreached().to_vec()
}

/// A device with no λ declarations breaks the chain: the net past it is
/// unreached even though light is flowing through.
///
/// `fc_cw_laser` → `fc_waveguide` both declare, so `a` and `b` resolve. There is
/// no native device that carries light and declares nothing — every one of them
/// is declared now — so the condition is asserted the only way it can be
/// without a compiler: the chain resolves end to end, and a *dark* branch does
/// not. If a future device forgets to declare, the first half of this test is
/// what changes.
#[test]
fn a_declared_chain_resolves_end_to_end_and_a_dark_branch_does_not() {
    let dark = unreached_nets(
        "\
.optical_port a
.optical_port b
.optical_port lonely
Xl a fc_cw_laser power_mW=1 wavelength_nm=1543.0
Xw a b fc_waveguide L_um=100 n_g=4.2
Xw2 lonely b fc_waveguide L_um=100 n_g=4.2
.op
",
    );
    // `lonely` is an input nothing drives — legitimately dark, and reported as
    // such rather than silently taken for a wavelength.
    assert!(
        dark.iter().any(|n| n.starts_with("lonely")),
        "an undriven optical input should be reported unreached, got {dark:?}"
    );
    assert!(
        !dark.iter().any(|n| n.starts_with("a_") || n == "a"),
        "the lit path must resolve, got {dark:?}"
    );
}

/// A wavelength named by hand on the wire seeds resolution. This is the idiom
/// every hand-wired bundle in the tree uses, and the only way to label light
/// arriving from outside the deck — substituting the band centre for a wire
/// someone explicitly drove to 1551 nm would be the wrong answer, quietly.
#[test]
fn a_voltage_source_on_a_lambda_net_names_the_wavelength() {
    let deck = "\
.optical_port a
.optical_port b
Vre a_re_0 0 DC 0.03162
Vwl a_wl_0 0 DC 1.551e-6
Xw a b fc_waveguide L_um=100 n_g=4.2
.op
";
    let net = parse_spice(deck).expect("deck parses");
    let reg = DeviceRegistry::new();
    let opts = fairchild_core::SimOptions::default();
    let map = lambda::resolve(&net, &opts.sim_context(), &reg);
    let got = map.get("b_wl_0").expect("the output λ net resolves");
    assert!(
        (got - 1.551e-6).abs() < 1e-15,
        "a hand-driven λ wire must propagate: got {got:e} m, expected 1.551e-6 m"
    );
}
