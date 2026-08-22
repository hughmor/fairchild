//! What happens when a device's "unused" optical port turns out to carry light.
//!
//! Every device in this file is correct on its own. The defect this suite hunts
//! only appears in composition: a device that drives a wire it does not own
//! meets a neighbour driving the same wire, and the block goes rank-deficient
//! without the solve reporting anything — `gmin` on the two branch diagonals
//! breaks the exact dependence, so LU succeeds and returns a weighted average
//! of the two assertions.
//!
//! So every case here closes a **power budget** rather than checking that two
//! stamps agree with each other. An agreement invariant cannot see a fault the
//! two sides share, and "the optical stamps agree" is exactly such an
//! invariant. The anchor on the other side is always a number computed by hand
//! from the deck: for a lossless chain ended in a perfect mirror it is
//! "everything launched comes back", which no averaged answer satisfies.
//!
//! Wire order under `enable_bidirectional=1` is `[re_fw, im_fw, re_bw, im_bw, λ]`.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

fn run(body: &str) -> fairchild_core::newton::NrResult {
    let src = format!(".options enable_bidirectional=1\n{body}.op\n");
    let net = parse_spice(&src).expect("netlist should parse");
    dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP should converge")
}

/// Optical power on one bundle wire pair, W. 1 V of `re` is 1 W.
fn power(r: &fairchild_core::newton::NrResult, port: &str, dir: &str, ch: usize) -> f64 {
    let v = |w: &str| {
        r.node_voltage(&format!("{port}_{w}_{dir}_{ch}"))
            .unwrap_or_else(|e| panic!("{port}_{w}_{dir}_{ch}: {e}"))
    };
    let (re, im) = (v("re"), v("im"));
    re * re + im * im
}

/// Power transmission of `l_um` of waveguide at `alpha_db_cm`, one way.
fn one_way(l_um: f64, alpha_db_cm: f64) -> f64 {
    10f64.powf(-alpha_db_cm * (l_um * 1e-4) / 10.0)
}

// ── fc_mux / fc_demux ──────────────────────────────────────────────────────

/// Two lasers, a mux, a lossless waveguide and a perfect mirror. Everything
/// launched has to come back out of the channel port it was launched on.
///
/// The anchor is the launched power itself, so nothing internal to the optical
/// stamps can satisfy it by agreeing with itself. Asymmetric powers make the
/// two channels distinguishable: if the return path mixed them, both would come
/// back at the mean instead of at their own value.
///
/// A mux used to stamp its backward pair in the forward direction — driving the
/// bus's backward wires, which the waveguide downstream already drives, and
/// leaving the channel ports' backward wires driven by nobody. The reflection
/// then read exactly zero at both channel ports while the deck converged
/// happily.
#[test]
fn a_mirror_behind_a_mux_returns_every_launched_watt() {
    let r = run("\
.optical_port ch0\n\
.optical_port ch1\n\
.optical_port bus 2\n\
.optical_port far 2\n\
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
Xl1 ch1 fc_cw_laser power_mW=0.25 wavelength_nm=1551\n\
Xmux bus ch0 ch1 fc_mux\n\
Xwg  bus far fc_waveguide L_um=500 n_eff=2.445 n_g=4.2 alpha_dB_cm=0\n\
Xf   far fc_facet reflectance=1.0\n");
    for (port, launched) in [("ch0", 1.0e-3), ("ch1", 0.25e-3)] {
        let back = power(&r, port, "bw", 0);
        assert!(
            (back - launched).abs() / launched < 1e-9,
            "{port}: launched {launched:.6e} W, {back:.6e} W came back"
        );
    }
}

/// The same chain with real numbers in it: the round trip pays the propagation
/// loss twice and the reflectance once, per channel, and the mux's own
/// insertion loss twice.
///
/// `il_db` is a field-symmetric per-channel loss, so a round trip through the
/// mux costs it twice — which also pins that the backward path goes through the
/// filter at all rather than around it.
#[test]
fn the_round_trip_through_a_mux_pays_every_loss_twice() {
    const L_UM: f64 = 1000.0;
    const ALPHA: f64 = 2.0;
    const R_POWER: f64 = 0.3;
    const IL_DB: f64 = 1.5;
    let r = run(&format!(
        "\
.optical_port ch0\n\
.optical_port ch1\n\
.optical_port bus 2\n\
.optical_port far 2\n\
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
Xl1 ch1 fc_cw_laser power_mW=0.25 wavelength_nm=1551\n\
Xmux bus ch0 ch1 fc_mux il_db={IL_DB}\n\
Xwg  bus far fc_waveguide L_um={L_UM} n_eff=2.445 n_g=4.2 alpha_dB_cm={ALPHA}\n\
Xf   far fc_facet reflectance={R_POWER}\n"
    ));
    let mux_pass = 10f64.powf(-IL_DB / 10.0);
    let trip = mux_pass * one_way(L_UM, ALPHA);
    for (port, launched) in [("ch0", 1.0e-3), ("ch1", 0.25e-3)] {
        let expect = launched * trip * R_POWER * trip;
        let back = power(&r, port, "bw", 0);
        assert!(
            (back - expect).abs() / expect < 1e-9,
            "{port}: {back:.6e} W came back, expected {expect:.6e} W"
        );
    }
}

/// A demux with a mirror on one output port and a terminator on the other.
/// Only the mirrored channel comes back, and it comes back whole.
///
/// A demux had the mirror-image of the mux's fault: it drove its channel ports'
/// backward wires, which every device wired onto a channel port already drives,
/// and never drove the bus's. `fc_facet` is that device here, so before the fix
/// this deck was two drivers on `d0_re_bw_0`.
#[test]
fn a_mirror_behind_a_demux_returns_only_its_own_channel() {
    let r = run("\
.optical_port ch0\n\
.optical_port ch1\n\
.optical_port bus 2\n\
.optical_port d0\n\
.optical_port d1\n\
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
Xl1 ch1 fc_cw_laser power_mW=0.25 wavelength_nm=1551\n\
Xmux  bus ch0 ch1 fc_mux\n\
Xdmx  bus d0 d1 fc_demux\n\
Xmir  d0 fc_facet reflectance=1.0\n\
Xterm d1 fc_facet\n");
    let back0 = power(&r, "ch0", "bw", 0);
    let back1 = power(&r, "ch1", "bw", 0);
    assert!(
        (back0 - 1.0e-3).abs() / 1.0e-3 < 1e-9,
        "the mirrored channel should return all 1 mW; got {back0:.6e} W"
    );
    assert_eq!(
        back1, 0.0,
        "the terminated channel should return nothing; got {back1:.6e} W"
    );
}

// ── fc_circulator ──────────────────────────────────────────────────────────

/// The deck a circulator exists for: launch at port 1, stimulate a device under
/// test out of port 2, and read what comes back at port 3. With a lossless
/// waveguide and a perfect mirror as the DUT, every launched watt has to arrive
/// at port 3 — and none of it may leak back toward the laser.
///
/// Both ends of the circulator are reached through a waveguide rather than
/// wired to it directly, which is the part that used to be impossible: the
/// convention was port-relative (`fw` meant "into me" at all three ports), so
/// every port behaved like an `in` port and colliding with the neighbouring
/// device's `in` port was unavoidable.
#[test]
fn a_lossless_circulator_loop_delivers_every_watt_to_port_three() {
    let r = run("\
.optical_port src\n\
.optical_port p1\n\
.optical_port p2\n\
.optical_port p3\n\
.optical_port dut\n\
.optical_port out\n\
Xl   src fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
Xin  src p1 fc_waveguide L_um=250 n_eff=2.445 n_g=4.2 alpha_dB_cm=0\n\
Xc   p1 p2 p3 fc_circulator\n\
Xdw  p2 dut fc_waveguide L_um=500 n_eff=2.445 n_g=4.2 alpha_dB_cm=0\n\
Xm   dut fc_facet reflectance=1.0\n\
Xow  p3 out fc_waveguide L_um=750 n_eff=2.445 n_g=4.2 alpha_dB_cm=0\n\
Xt   out fc_facet\n");
    let delivered = power(&r, "out", "fw", 0);
    assert!(
        (delivered - 1.0e-3).abs() / 1.0e-3 < 1e-9,
        "port 3 should receive all 1 mW; got {delivered:.6e} W"
    );
    // Isolation: light entering port 2 leaves at port 3, never back at port 1.
    // Nothing drives port 3's return, so the laser's port sees exactly zero.
    assert_eq!(
        power(&r, "src", "bw", 0),
        0.0,
        "a circulator must not send the DUT's reflection back to the laser"
    );
}

/// The same loop with loss in it. The reflection pays the DUT waveguide twice
/// and the reflectance once; the launch and read waveguides once each.
#[test]
fn the_circulator_return_path_pays_the_dut_round_trip() {
    const R_POWER: f64 = 0.3;
    const ALPHA: f64 = 2.0;
    let r = run(&format!(
        "\
.optical_port src\n\
.optical_port p1\n\
.optical_port p2\n\
.optical_port p3\n\
.optical_port dut\n\
.optical_port out\n\
Xl   src fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
Xin  src p1 fc_waveguide L_um=250 n_eff=2.445 n_g=4.2 alpha_dB_cm={ALPHA}\n\
Xc   p1 p2 p3 fc_circulator\n\
Xdw  p2 dut fc_waveguide L_um=500 n_eff=2.445 n_g=4.2 alpha_dB_cm={ALPHA}\n\
Xm   dut fc_facet reflectance={R_POWER}\n\
Xow  p3 out fc_waveguide L_um=750 n_eff=2.445 n_g=4.2 alpha_dB_cm={ALPHA}\n\
Xt   out fc_facet\n"
    ));
    let expect = 1.0e-3
        * one_way(250.0, ALPHA)
        * one_way(500.0, ALPHA)
        * R_POWER
        * one_way(500.0, ALPHA)
        * one_way(750.0, ALPHA);
    let got = power(&r, "out", "fw", 0);
    assert!(
        (got - expect).abs() / expect < 1e-9,
        "port 3 delivered {got:.6e} W, expected {expect:.6e} W"
    );
}

// ── fc_photodetector ───────────────────────────────────────────────────────

/// A photodetector terminates an optical path without claiming the backward
/// wire: it binds no auxiliary row on any optical wire at all, so it absorbs
/// whatever arrives and asserts nothing about what does not.
///
/// Saying that needs a deck where the PD's backward wire is driven by a real
/// device rather than a `V` source — a source would win the argument outright,
/// because its branch row carries no `gmin` while a device's does. The deck
/// would then pass whether or not the PD claimed the wire, which is a test that
/// cannot fail. A circulator supplies a real driver: the loop below sends the
/// launch out of port 2, through the DUT to a mirror, back in at port 2, out at
/// port 3, off a second mirror, and back in at port 3 — from where it leaves at
/// port 1, on the very wire the PD is sitting on.
///
/// The PD therefore sees both directions at once, and the budget is exact: the
/// laser's outgoing field plus the round trip's return, one photocurrent.
#[test]
fn a_photodetector_absorbs_backward_light_without_owning_the_wire() {
    const P_MW: f64 = 1.0;
    const L_UM: f64 = 1000.0;
    const ALPHA: f64 = 2.0;
    const RESP: f64 = 0.8;
    let r = run(&format!(
        "\
.optical_port p1\n\
.optical_port p2\n\
.optical_port p3\n\
.optical_port dut\n\
Xl   p1 fc_cw_laser power_mW={P_MW} wavelength_nm=1550\n\
Xc   p1 p2 p3 fc_circulator\n\
Xdw  p2 dut fc_waveguide L_um={L_UM} n_eff=2.445 n_g=4.2 alpha_dB_cm={ALPHA}\n\
Xm2  dut fc_facet reflectance=1.0\n\
Xm3  p3 fc_facet reflectance=1.0\n\
Xpd  p1 a 0 fc_photodetector responsivity={RESP} i_dark_a=0 r_shunt=1e12\n\
Rl   a 0 1k\n"
    ));
    // Out at port 2, twice through the DUT waveguide, in at port 2, out at
    // port 3, straight back in at port 3, out at port 1.
    let p_launch = P_MW * 1e-3;
    let p_return = p_launch * one_way(L_UM, ALPHA).powi(2);
    let back = power(&r, "p1", "bw", 0);
    assert!(
        (back - p_return).abs() / p_return < 1e-9,
        "the loop returned {back:.6e} W to port 1, expected {p_return:.6e} W"
    );
    // One junction, both directions, one photocurrent.
    let expect_v = RESP * (p_launch + p_return) * 1e3;
    let got_v = r.node_voltage("a").unwrap().abs();
    assert!(
        (got_v - expect_v).abs() / expect_v < 1e-6,
        "PD load read {got_v:.6e} V, expected {expect_v:.6e} V"
    );
}

// ── the structural guard ───────────────────────────────────────────────────

/// Two devices pinning one node is refused by name, not averaged.
///
/// This is the backstop for every device this suite does not cover and every
/// one not written yet. Two `out` ports on the same bundle both drive its
/// forward wires — the same shape as the faults above, with nothing optical
/// about it.
#[test]
fn two_devices_driving_one_wire_is_an_error_naming_both() {
    let net = parse_spice(
        ".options enable_bidirectional=1\n\
         .optical_port a\n.optical_port b\n.optical_port c\n\
         Xw1 a b fc_waveguide L_um=100\n\
         Xw2 c b fc_waveguide L_um=100\n\
         .op\n",
    )
    .unwrap();
    let Err(err) = dc_op_nr_with_registry(&net, &DeviceRegistry::new()) else {
        panic!("two waveguide outputs on one bundle must not build");
    };
    let msg = format!("{err}");
    assert!(
        msg.contains("xw1") && msg.contains("xw2") && msg.contains("b_re_fw_0"),
        "the error must name both devices and the wire: {msg}"
    );
}

/// The unidirectional path is untouched: with three wires per channel there is
/// no backward pair to get the direction of, and the same decks read the same
/// numbers they always did.
#[test]
fn the_unidirectional_route_is_unchanged() {
    let net = parse_spice(
        "\
.optical_port ch0\n.optical_port ch1\n.optical_port bus 2\n\
.optical_port d0\n.optical_port d1\n\
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
Xl1 ch1 fc_cw_laser power_mW=0.25 wavelength_nm=1551\n\
Xmux bus ch0 ch1 fc_mux\n\
Xdmx bus d0 d1 fc_demux\n\
Xp0 d0 a0 0 fc_photodetector responsivity=1.0 i_dark_a=0 r_shunt=1e12\n\
Xp1 d1 a1 0 fc_photodetector responsivity=1.0 i_dark_a=0 r_shunt=1e12\n\
R0 a0 0 1k\nR1 a1 0 1k\n.op\n",
    )
    .unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    // The photocurrent flows cathode → anode, so the load pulls the anode
    // below ground; the magnitude is the budget.
    for (node, p) in [("a0", 1.0e-3), ("a1", 0.25e-3)] {
        let v = r.node_voltage(node).unwrap().abs();
        // 1e-6, not tighter: the 1e12 Ω shunt across the 1 kΩ load costs 2 ppb.
        assert!(
            (v - p * 1e3).abs() / (p * 1e3) < 1e-6,
            "{node}: {v:.6e} V, expected {:.6e} V",
            p * 1e3
        );
    }
}

// ── absorbed power does not care which way a photon travels (#51) ──────────

/// `channel_intensities` read the forward wires only, so under
/// `enable_bidirectional=1` the returning field heated nothing: a junction in
/// front of a mirror saw half the light it was actually absorbing.
///
/// The anchor is a closed form fed by the backward field the solve itself
/// reports, so it is not an agreement between the segment and the drive — the
/// prediction is arithmetic on a measured amplitude:
///
///   Δn      = dn_dt · R_th · α · L · (I_fw_in + I_bw_in)
///   Δφ(R)   = 2π · L · (Δn(R) − Δn(0)) / λ = 2π · L · C · I_bw_in / λ
///
/// The forward input is the laser and does not move with `reflectance`, so the
/// whole difference between the two runs is the backward field. Against the old
/// code the difference is exactly zero.
#[test]
fn backward_light_heats_the_junction_it_passes_through() {
    let deck = |refl: f64| {
        format!(
            ".optical_port a\n\
             .optical_port b\n\
             Xl a fc_cw_laser power_mW=10 wavelength_nm=1550\n\
             Xpn a b c 0 fc_pn_ps_full L_um=500 alpha_dB_cm=2.0 dn_dt=1e-4 r_th=5000 beta_tpa=0\n\
             Vb c 0 DC 0\n\
             Xf b fc_facet reflectance={refl}\n"
        )
    };
    let phase = |r: &fairchild_core::newton::NrResult| {
        let re = r.node_voltage("b_re_fw_0").unwrap();
        let im = r.node_voltage("b_im_fw_0").unwrap();
        im.atan2(re)
    };

    let dark = run(&deck(0.0));
    let lit = run(&deck(0.9));

    let i_bw = power(&lit, "b", "bw", 0);
    assert!(
        i_bw > 1e-4,
        "the mirror must actually return light: I_bw={i_bw:.6} W"
    );

    let alpha = 2.0 * (10.0f64).ln() / 10.0 * 100.0; // dB/cm → Np/m
    let l_m = 500e-6;
    let c = 1e-4 * 5000.0 * alpha * l_m; // Δn per watt of absorbed-power input
    let expect = 2.0 * std::f64::consts::PI * l_m * (c * i_bw) / 1550e-9;
    let got = (phase(&dark) - phase(&lit)).abs();

    assert!(
        (got - expect).abs() < 2e-4,
        "backward light contributed a Δφ of {got:.6} rad; the field the solve \
         reports (I_bw={i_bw:.6} W) implies {expect:.6} rad"
    );
}
