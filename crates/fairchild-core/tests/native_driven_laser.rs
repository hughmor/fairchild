//! `fc_driven_laser` — optical power that follows an electrical input, so one
//! SPICE source produces a modulated optical waveform with no external
//! modulator in the deck.

use fairchild_core::{dc_op_nr_with_registry, tran_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

const SLOPE: f64 = 1e-3; // W/V
const V_TH: f64 = 0.5;
const P_FLOOR: f64 = 1e-12;

fn deck(drive: &str, extra: &str) -> String {
    format!(
        "{drive}\n\
         Xlas o_re o_im o_wl drv 0 fc_driven_laser \
           slope_w_v={SLOPE} v_th={V_TH} p_floor_w={P_FLOOR} r_in=1e12 wavelength_nm=1550\n\
         {extra}"
    )
}

/// Optical power at the laser's port for a fixed drive voltage.
fn power_at(v: f64) -> f64 {
    let net = parse_spice(&deck(&format!("Vd drv 0 DC {v}"), ".op\n.end\n")).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let re = r.node_voltage("o_re").unwrap();
    let im = r.node_voltage("o_im").unwrap();
    re * re + im * im
}

/// Above threshold the L–V curve is exactly the straight line it claims to be.
#[test]
fn power_follows_the_drive_above_threshold() {
    for v in [0.6, 1.0, 2.0, 5.0] {
        let want = P_FLOOR + SLOPE * (v - V_TH);
        let got = power_at(v);
        assert!(
            (got - want).abs() / want < 1e-12,
            "V={v}: {got:.9e} W, expected {want:.9e} W"
        );
    }
}

/// Below threshold the output sits on the spontaneous-emission floor, not at
/// zero — and not at a negative power, which is what an unclamped straight
/// line would ask `sqrt` for.
///
/// "Below" means by more than a few `v_knee`, the width over which the L-V
/// curve bends. At 50 knee widths down the softplus contributes `e^-50`, which
/// is 10 orders below the floor itself. Exactly *at* threshold it does not —
/// see `at_threshold_the_knee_carries_half_the_slope`.
#[test]
fn below_threshold_the_output_is_the_floor() {
    for v in [-2.0, 0.0, V_TH - 0.05] {
        let got = power_at(v);
        assert!(
            (got - P_FLOOR).abs() < 1e-18,
            "V={v}: {got:.6e} W, expected the {P_FLOOR:.0e} W floor"
        );
    }
}

/// At threshold the curve bends rather than corners, so `dP/dV` is half the
/// slope and the power sits `slope·v_knee·ln2` above the floor.
///
/// This is the whole fix for the falling-edge convergence failure, so it is
/// worth pinning as a number rather than as a vibe: a hard `max(0, ·)` would
/// put exactly `P_FLOOR` here and hand Newton a `dA/dV` of `slope/(2√P_FLOOR)`
/// = 500, discontinuously.
#[test]
fn at_threshold_the_knee_carries_half_the_slope() {
    const V_KNEE: f64 = 1e-3; // the default
    let want = P_FLOOR + SLOPE * V_KNEE * std::f64::consts::LN_2;
    let got = power_at(V_TH);
    assert!(
        (got - want).abs() / want < 1e-9,
        "at threshold: {got:.6e} W, expected {want:.6e} W"
    );
    // And the resulting Jacobian entry is small, which is the point.
    let da_dv = 0.5 * SLOPE / (2.0 * got.sqrt());
    assert!(
        da_dv < 1.0,
        "dA/dV at the knee is {da_dv:.1}; the hard corner gave 500"
    );
}

/// The floor is what keeps `dA/dV = slope/(2√P)` finite at threshold. Without
/// it the Jacobian entry diverges exactly where a modulated laser spends every
/// falling edge.
#[test]
fn the_floor_bounds_the_jacobian_at_threshold() {
    // Sweep across the corner in small steps; each one has to converge, and
    // the field has to stay monotone through it.
    let mut prev = 0.0;
    let mut v = V_TH - 1e-6;
    while v < V_TH + 1e-5 {
        let a = power_at(v).sqrt();
        assert!(
            a >= prev,
            "field went backwards crossing threshold at V={v}"
        );
        prev = a;
        v += 1e-6;
    }
    // The steepest slope the model can present, at the floor itself.
    let max_da_dv = SLOPE / (2.0 * P_FLOOR.sqrt());
    assert!(max_da_dv < 1e3, "dA/dV bound {max_da_dv:.3e} is too stiff");
}

/// The drive derivative is stamped, not frozen: a closed opto-electronic loop
/// converges in a handful of Newton iterations.
///
/// The laser feeds a photodiode whose photocurrent develops the very voltage
/// that drives the laser — loop gain `responsivity · slope · R_load`, set to
/// 0.95 here. With `∂A/∂V` stamped, Newton solves it. Frozen (the coefficient
/// refreshed from the previous iterate, which is what most photonic devices in
/// this tree do for their optical coefficients) it degenerates into successive
/// substitution and needs ~135 iterations to contract 0.95^n below `reltol`,
/// which is why the iteration limit is part of the assertion.
#[test]
fn an_optoelectronic_loop_converges_because_the_derivative_is_stamped() {
    const RESPONSIVITY: f64 = 0.8;
    const R_LOAD: f64 = 1e3;
    const I_SEED: f64 = 1e-6;
    const LOOP_GAIN: f64 = 0.95;
    let slope = LOOP_GAIN / (RESPONSIVITY * R_LOAD);

    let src = format!(
        "* laser → PD → the laser's own drive node\n\
         .options itl1=15\n\
         Xlas o_re o_im o_wl drv 0 fc_driven_laser \
           slope_w_v={slope} v_th=0 p_floor_w=1e-15 r_in=1e12 wavelength_nm=1550\n\
         Xpd  o_re o_im o_wl drv 0 fc_photodetector \
           responsivity={RESPONSIVITY} r_shunt=1e12 i_dark_a=0\n\
         Rl   drv 0 {R_LOAD}\n\
         Iseed 0 drv {I_SEED}\n\
         .op\n.end\n"
    );
    let net = parse_spice(&src).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new())
        .expect("the loop must converge inside itl1=15");
    let v = r.node_voltage("drv").unwrap();
    // V = (I_seed + responsivity·slope·V)·R_load  ⇒  V = I_seed·R_load / (1 − g)
    let want = I_SEED * R_LOAD / (1.0 - LOOP_GAIN);
    assert!(
        (v - want).abs() / want < 1e-6,
        "loop settled at {v:.9e} V, expected {want:.9e} V"
    );
}

/// The point of the device: one voltage source, a modulated optical output,
/// and a detector that sees it — no external modulator anywhere in the deck.
#[test]
fn a_pulsed_drive_produces_a_modulated_optical_waveform() {
    let src = deck(
        "Vd drv 0 PULSE(0 2 0 10p 10p 490p 1n)",
        "Xpd o_re o_im o_wl pa 0 fc_photodetector responsivity=0.8 r_shunt=1Meg i_dark_a=0\n\
         Rl pa 0 1k\n\
         .tran 5p 3n\n.end\n",
    );
    let net = parse_spice(&src).unwrap();
    let r = tran_nr_with_registry(&net, 5e-12, 3e-9, &DeviceRegistry::new())
        .expect("transient should run");
    let v = &r.node_voltages["pa"];

    let hi = v.iter().cloned().fold(f64::MIN, f64::max);
    let lo = v.iter().cloned().fold(f64::MAX, f64::min);
    // Drive high = 2 V ⇒ P = 1.5 mW ⇒ I = 1.2 mA ⇒ 1.2 V across 1 kΩ.
    let want_hi = 0.8 * (P_FLOOR + SLOPE * (2.0 - V_TH)) * 1e3;
    assert!(
        (hi - want_hi).abs() / want_hi < 1e-3,
        "peak {hi:.6} V, expected {want_hi:.6} V"
    );
    // Drive low = 0 V ⇒ the floor, 90 dB down. An extinction ratio that big is
    // the floor's whole design brief: numerically safe, optically invisible.
    assert!(lo < want_hi * 1e-6, "off level {lo:.3e} V is not dark");
}

/// **Modulating through threshold must converge with nothing reactive in the
/// circuit.**
///
/// The L-V curve used to be `P = p_floor + max(0, slope·(V−V_th))`, so
/// `dA/dV = slope/(2√P)` jumped from 0 below threshold to `slope/(2√p_floor)`
/// above it — 500 here, three orders above anything else in the Jacobian, and
/// switching on and off as Newton stepped across the corner. The iterate
/// ping-ponged and never converged.
///
/// It looked fine for months because every deck in the tree had a load
/// capacitor, whose `C/h` on the detector diagonal damped the oscillation.
/// Remove the capacitor — which is the only thing this deck does differently
/// from `a_pulsed_drive_produces_a_modulated_optical_waveform` — and the
/// falling edge failed with `NoConvergence` and no hint as to why.
///
/// `v_knee` smooths the knee over a few mV, which is what a real laser's L-I
/// curve does anyway. Set `v_knee=0` to restore the corner: this test then
/// fails, which is how the fix was verified.
#[test]
fn modulating_through_threshold_converges_with_no_reactive_element() {
    let src = deck(
        // Falls from 2 V to 0 V, i.e. straight through V_TH = 0.5.
        "Vd drv 0 PULSE(0 2 0 30p 30p 70p 200p)",
        "Xpd o_re o_im o_wl pa 0 fc_photodetector responsivity=0.8 r_shunt=1Meg i_dark_a=0\n\
         Rl pa 0 1k\n\
         .tran 1p 400p\n.end\n",
    );
    let net = parse_spice(&src).unwrap();
    let r = tran_nr_with_registry(&net, 1e-12, 400e-12, &DeviceRegistry::new())
        .expect("a purely resistive link must still converge through threshold");
    let v = &r.node_voltages["pa"];

    // And it must converge to the right waveform, not merely to something.
    let hi = v.iter().cloned().fold(f64::MIN, f64::max);
    let want_hi = 0.8 * (P_FLOOR + SLOPE * (2.0 - V_TH)) * 1e3;
    assert!(
        (hi - want_hi).abs() / want_hi < 1e-3,
        "peak {hi:.6} V, expected {want_hi:.6} V"
    );
    assert!(
        v.iter().cloned().fold(f64::MAX, f64::min) < want_hi * 1e-6,
        "the off level must still be dark"
    );
}

/// The smooth knee must not move the L-V curve anywhere anyone uses it.
///
/// `softplus(u)` and `max(0,u)` differ by less than `e^-u`, so at more than a
/// few tens of `v_knee` above threshold they agree to machine precision. This
/// pins that: the same drive, with and without the smoothing, to 1e-9.
#[test]
fn the_threshold_knee_does_not_move_the_curve_above_threshold() {
    // Read the field wires rather than a detector's load voltage: the branch
    // rows enforce them exactly, so this compares the model and not Newton's
    // stopping tolerance.
    let power = |v: f64, knee: &str| {
        let src = format!(
            "Vd drv 0 DC {v}\n\
             Xlas o_re o_im o_wl drv 0 fc_driven_laser \
               slope_w_v={SLOPE} v_th={V_TH} p_floor_w={P_FLOOR} r_in=1e12 {knee}\n\
             .op\n.end\n"
        );
        let net = parse_spice(&src).unwrap();
        let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("op");
        let re = r.node_voltage("o_re").unwrap();
        let im = r.node_voltage("o_im").unwrap();
        re * re + im * im
    };
    for v in [0.6, 1.0, 2.0, 5.0] {
        let smooth = power(v, ""); // default v_knee = 1 mV
        let corner = power(v, "v_knee=0"); // the old hard max(0, ·)
        assert!(
            (smooth - corner).abs() / corner < 1e-9,
            "V={v} is {} knee widths above threshold; the two forms must agree: \
             {smooth:.9e} vs {corner:.9e}",
            (v - V_TH) / 1e-3
        );
    }
}
