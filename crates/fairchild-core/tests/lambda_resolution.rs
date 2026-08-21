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

    // Same devices the solve used, built the same way.
    let opts = fairchild_core::SimOptions::default();
    let ctx = opts.sim_context();
    let mut topo = fairchild_core::CircuitTopology::build(&net);
    let devices = fairchild_core::build_devices(&net, &mut topo, &ctx, &reg)
        .expect("devices build for the resolution pass");
    let map = lambda::resolve(&net, &devices, ctx.lambda_center_m);

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
.end
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
.end
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
.end
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
.end
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
.end
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
    deck += " fc_mux\nXpn bus out a 0 fc_pn_ps_full L_um=50\nVb a 0 DC 0\n.op\n.end\n";
    assert_resolution_matches_the_solve("giona-shaped bus", &deck);
}
