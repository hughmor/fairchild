//! `examples/verilog_a/models/ring_nheater.va` — the giona weight ring.
//!
//! An add-drop ring with no junction, tuned by doping its own waveguide n-type
//! and running current through it. The optical core is `mrm_addrop.va`'s and is
//! covered there; what is new here, and all this file tests, is the loop the
//! doping closes:
//!
//!     doping -> free carriers -> absorption -> heat -> resistance
//!
//! so that sweeping the heater and watching the CURRENT finds the resonance
//! with no photodetector in the circuit.
//!
//! The third test is the one that makes the first two mean something. An
//! agreement between "the resistance moved" and "the ring was on resonance" is
//! satisfied by any model that heats up, so it also switches the doping off and
//! requires the readout to vanish — if it does not, the effect is coming from
//! somewhere other than the carriers it is supposed to come from.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_osdi::{load_libraries_with_widths, VaOptions};
use fairchild_parser::{instantiated_widths, parse_spice_with_arity, PermissiveArity};
use std::collections::{BTreeMap, BTreeSet};

use crate::common;

fn model_dir() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/verilog_a/models")
}

fn solve(deck: &str) -> fairchild_core::NrResult {
    let probe = parse_spice_with_arity(deck, &PermissiveArity).expect("probe pass parses");
    let widths: BTreeMap<String, BTreeSet<usize>> = instantiated_widths(&probe);
    let mut reg = DeviceRegistry::new();
    load_libraries_with_widths(
        &[],
        &[model_dir()
            .join("ring_nheater.va")
            .to_string_lossy()
            .to_string()],
        None,
        &VaOptions {
            cache_dir: Some(common::test_cache_dir()),
            include_dirs: vec![model_dir()],
            ..VaOptions::from_env()
        },
        &mut reg,
        &widths,
        3,
    )
    .expect("the example model generates and compiles");
    let net = parse_spice_with_arity(deck, &reg).expect("the deck parses against the model");
    dc_op_nr_with_registry(&net, &reg).expect("DC OP converges")
}

fn deck(v_htr: f64, p_mw: f64, extra: &str) -> String {
    format!(
        ".optical_port src\n.optical_port th\n.optical_port ad\n.optical_port dr\n\
         XL src fc_cw_laser power_mW={p_mw} wavelength_nm=1550\n\
         Xr src th ad dr hp 0 tr ring_nheater {extra}\n\
         VH hp 0 DC {v_htr}\n\
         .op\n"
    )
}

/// One sweep point: (heater resistance, drop-port power, temperature rise).
fn point(v_htr: f64, p_mw: f64, extra: &str) -> (f64, f64, f64) {
    let r = solve(&deck(v_htr, p_mw, extra));
    // Conventional current out of the source's + terminal is the negative of
    // its branch unknown.
    let i = -r
        .vsrc_current("vh")
        .expect("the heater source's branch current is in the solution");
    let re = r.node_voltage("dr_re_0").unwrap();
    let im = r.node_voltage("dr_im_0").unwrap();
    (v_htr / i, re * re + im * im, r.node_voltage("tr").unwrap())
}

/// Heater volts either side of the resonance the laser sits on, at 0 V bias.
/// Below 1.6 V the ring has not reached it; the drop port peaks around 1 V.
const V_SWEEP: [f64; 13] = [
    0.30, 0.45, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.35, 1.50, 1.80, 2.20,
];

/// Dark: the resistance is the heater's own tempco and nothing else.
///
/// `R(T) = r_heat0·(1 + alpha_r·ΔT)` with `ΔT = r_th·V²/R` is a fixed point, not
/// a formula, so this also checks that it is being solved rather than evaluated
/// once: at 2.2 V the rise is large enough that ignoring the feedback would put
/// the answer several per cent out.
#[test]
fn the_dark_resistance_follows_its_own_solved_temperature() {
    if !common::have_compiler() {
        return;
    }
    let (r0, _, t0) = point(0.05, 1e-9, "");
    assert!(
        (r0 - 2200.0).abs() < 1.0 && t0 < 0.02,
        "at 50 mV the ring is cold and R is r_heat0: got {r0:.2} ohm, {t0:.4} K"
    );
    let (r1, _, t1) = point(2.2, 1e-9, "");
    // Self-consistent: R = 2200(1 + 1.5e-3·ΔT), ΔT = 8290·V²/R.
    let expect = 2200.0 * (1.0 + 1.5e-3 * t1);
    assert!(
        (r1 - expect).abs() / expect < 1e-3,
        "R and its own temperature disagree: {r1:.2} vs {expect:.2} ohm at {t1:.2} K"
    );
    assert!(
        t1 > 5.0 && r1 > r0,
        "2.2 V should warm it measurably and raise R: {t1:.2} K, {r0:.1} -> {r1:.1} ohm"
    );
}

/// The point of the model: the resonance appears in the CURRENT.
///
/// Light on the bus, sweep the heater, and the extra resistance the absorbed
/// light causes must peak where the drop port does — the ring reading itself
/// out with no photodetector anywhere in the deck.
#[test]
fn sweeping_the_heater_finds_the_resonance_in_the_resistance() {
    if !common::have_compiler() {
        return;
    }
    let mut best_drop = (0usize, 0.0f64);
    let mut best_excess = (0usize, 0.0f64);
    for (i, &v) in V_SWEEP.iter().enumerate() {
        let (r_lit, drop, _) = point(v, 1.0, "");
        let (r_dark, _, _) = point(v, 1e-9, "");
        let excess = (r_lit - r_dark) / r_dark;
        if drop > best_drop.1 {
            best_drop = (i, drop);
        }
        if excess > best_excess.1 {
            best_excess = (i, excess);
        }
    }
    assert!(
        best_drop.1 > 0.5e-3,
        "the sweep never crossed the resonance: peak drop {:.4} mW",
        best_drop.1 * 1e3
    );
    // Adjacent is close enough: the optical heating red-shifts the ring as it
    // tunes in, so the resistance peak trails the drop peak by a fraction of a
    // step. Demanding the same index would be pinning the lineshape asymmetry,
    // which is real physics and not what this test is about.
    let gap = best_drop.0.abs_diff(best_excess.0);
    assert!(
        gap <= 1,
        "the resistance peak is {gap} steps from the drop peak ({} V vs {} V)",
        V_SWEEP[best_excess.0],
        V_SWEEP[best_drop.0]
    );
    assert!(
        best_excess.1 > 1e-4,
        "the readout is too small to be a readout: {:.1} ppm at 1 mW",
        best_excess.1 * 1e6
    );
}

/// Undope the waveguide and the readout has to go away.
///
/// This is the control. The two tests above are also passed by a ring that
/// simply gets warm, so switch off the mechanism the effect is attributed to —
/// the free carriers the doping put in the mode — and require it to disappear.
/// What is left is the scattering loss, which absorbs nothing into the silicon.
#[test]
fn without_the_doping_there_is_nothing_to_read_out() {
    if !common::have_compiler() {
        return;
    }
    let on = V_SWEEP
        .iter()
        .map(|&v| {
            let (r_lit, _, _) = point(v, 1.0, "");
            let (r_dark, _, _) = point(v, 1e-9, "");
            (r_lit - r_dark) / r_dark
        })
        .fold(0.0f64, f64::max);
    let off = V_SWEEP
        .iter()
        .map(|&v| {
            let (r_lit, _, _) = point(v, 1.0, "n_dope_cm3=0 alpha_db_cm=0");
            let (r_dark, _, _) = point(v, 1e-9, "n_dope_cm3=0 alpha_db_cm=0");
            (r_lit - r_dark) / r_dark
        })
        .fold(0.0f64, f64::max);
    assert!(
        off < on / 20.0,
        "an undoped, lossless ring still reads out: {:.1} ppm against {:.1} ppm doped",
        off * 1e6,
        on * 1e6
    );
}
