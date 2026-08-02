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
#[test]
fn below_threshold_the_output_is_the_floor() {
    for v in [-2.0, 0.0, V_TH] {
        let got = power_at(v);
        assert!(
            (got - P_FLOOR).abs() < 1e-18,
            "V={v}: {got:.6e} W, expected the {P_FLOOR:.0e} W floor"
        );
    }
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
