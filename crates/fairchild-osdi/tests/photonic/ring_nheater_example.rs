//! `examples/verilog_a/models/ring_nheater.va` — the giona weight ring.
//!
//! An add-drop ring with no junction, tuned by doping its own waveguide n-type
//! and running current through it. The optical core is `mrm_addrop.va`'s and is
//! covered there; what is new here, and all this file tests, is the loop the
//! doping closes:
//!
//!     doping -> free carriers -> absorption -> heat -> resistance
//!
//!     absorbed photons -> carriers -> conductance up -> RESISTANCE DOWN
//!
//! and the competing one, which is weaker and pushes the other way:
//!
//!     absorbed photons -> heat -> positive tempco -> resistance up
//!
//! so that ramping the heater current and watching the voltage finds the
//! resonance as a DIP in resistance, with no photodetector in the circuit.
//!
//! Everything here drives CURRENT, not voltage, because the self-heating loop
//! has opposite signs in the two: at constant current a falling R dissipates
//! less and falls further, at constant voltage it dissipates more and is pulled
//! back. The last test pins that difference rather than leaving it asserted.
//!
//! The doping control is what makes the rest mean anything. "The resistance
//! moved when the ring was on resonance" is satisfied by any model that heats
//! up — so it switches off the surface-state absorption, leaving only the
//! thermal half, and requires the sign to INVERT.

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

/// Current drive, which is how these are characterised — see the header.
fn deck(i_ma: f64, p_mw: f64, extra: &str) -> String {
    let i_a = i_ma * 1e-3;
    format!(
        ".optical_port src\n.optical_port th\n.optical_port ad\n.optical_port dr\n\
         XL src fc_cw_laser power_mW={p_mw} wavelength_nm=1550\n\
         Xr src th ad dr hp 0 tr ring_nheater {extra}\n\
         IH 0 hp DC {i_a}\n\
         .op\n"
    )
}

/// Voltage drive, for the one test that compares the two.
fn deck_v(v_htr: f64, p_mw: f64, extra: &str) -> String {
    format!(
        ".optical_port src\n.optical_port th\n.optical_port ad\n.optical_port dr\n\
         XL src fc_cw_laser power_mW={p_mw} wavelength_nm=1550\n\
         Xr src th ad dr hp 0 tr ring_nheater {extra}\n\
         VH hp 0 DC {v_htr}\n\
         .op\n"
    )
}

/// One sweep point: (heater resistance, drop-port power, temperature rise).
fn point(i_ma: f64, p_mw: f64, extra: &str) -> (f64, f64, f64) {
    let r = solve(&deck(i_ma, p_mw, extra));
    let v = r.node_voltage("hp").unwrap();
    let re = r.node_voltage("dr_re_0").unwrap();
    let im = r.node_voltage("dr_im_0").unwrap();
    (
        v / (i_ma * 1e-3),
        re * re + im * im,
        r.node_voltage("tr").unwrap(),
    )
}

/// Heater milliamps, spanning the published 0-1.25 mA. The cold ring sits about
/// a third of a nanometre blue of the laser, so the drop port peaks near
/// 0.77 mA and the sweep straddles it.
const I_SWEEP: [f64; 11] = [
    0.05, 0.20, 0.35, 0.50, 0.62, 0.77, 0.90, 1.00, 1.10, 1.18, 1.25,
];

/// Dark: the resistance is the heater's own tempco and nothing else.
///
/// `R(T) = r_heat0·(1 + alpha_r·ΔT)` with `ΔT = r_th·V²/R` is a fixed point, not
/// a formula, so this also checks that it is being solved rather than evaluated
/// once: at 2.2 V the rise is large enough that ignoring the feedback would put
/// the answer several per cent out.
/// The dark curve, against the published characterisation it is matched to.
///
/// The electro-thermal parameters come from a paper quoting 134 ohm/mW of
/// self-heating on a 2.2 kohm heater; its dark sweep runs to about 2800 ohm and
/// just over 3 V at 1.25 mA. `R = r_heat0/(1 - alpha_r*r_th*r_heat0*I^2)` is a
/// fixed point, not a formula, so reproducing the endpoint also checks that it
/// is being solved: the closed form is 2785 ohm, and ignoring the feedback
/// would give 2200 + 134*4.35 = 2783 only by coincidence of this operating
/// point — at the halfway mark the two differ by 5 %.
#[test]
fn the_dark_current_sweep_matches_the_published_device() {
    if !common::have_compiler() {
        return;
    }
    // 20 uA still dissipates 0.9 uW, which is 0.03 K on a 38 400 K/W ring — so
    // "cold" here means a tenth of an ohm, not zero.
    let (r0, _, t0) = point(0.02, 1e-7, "");
    assert!(
        (r0 - 2200.0).abs() < 1.0 && t0 < 0.05,
        "at 20 uA the ring is all but cold and R is r_heat0: got {r0:.2} ohm, {t0:.4} K"
    );
    let (r1, _, t1) = point(1.25, 1e-7, "");
    let v1 = r1 * 1.25e-3;
    assert!(
        (r1 - 2785.0).abs() < 40.0 && (v1 - 3.48).abs() < 0.06,
        "1.25 mA dark should give ~2785 ohm and ~3.48 V (paper: ~2800, just over \
         3 V); got {r1:.1} ohm, {v1:.3} V"
    );
    // R and the temperature it is a function of have to agree with each other.
    let expect = 2200.0 * (1.0 + 1.5909e-3 * t1);
    assert!(
        (r1 - expect).abs() / expect < 1e-3,
        "R and its own temperature disagree: {r1:.2} vs {expect:.2} ohm at {t1:.2} K"
    );
}

/// The point of the model: the resonance appears in the CURRENT.
///
/// Light on the bus, sweep the heater, and the extra resistance the absorbed
/// light causes must peak where the drop port does — the ring reading itself
/// out with no photodetector anywhere in the deck.
/// The point of the model: the resonance appears in the resistance, as a DIP.
///
/// The published device shifts its whole R-vs-I curve DOWN under illumination,
/// by roughly 300 ohm on resonance. Down, because the carriers absorbed light
/// creates conduct, and that beats the heating those same carriers cause.
#[test]
fn the_resonance_shows_up_as_a_drop_in_resistance() {
    if !common::have_compiler() {
        return;
    }
    let mut best_drop = (0usize, 0.0f64);
    let mut deepest = (0usize, 0.0f64);
    for (i, &ma) in I_SWEEP.iter().enumerate() {
        let (r_lit, drop, _) = point(ma, 1.0, "");
        let (r_dark, _, _) = point(ma, 1e-7, "");
        let d_r = r_lit - r_dark;
        if drop > best_drop.1 {
            best_drop = (i, drop);
        }
        if d_r < deepest.1 {
            deepest = (i, d_r);
        }
    }
    assert!(
        best_drop.1 > 0.5e-3,
        "the sweep never crossed the resonance: peak drop {:.4} mW",
        best_drop.1 * 1e3
    );
    // Adjacent is close enough: the light's own heating red-shifts the ring as
    // it tunes in, so the electrical feature trails the optical one by a
    // fraction of a step. That asymmetry is real physics, not the thing under
    // test here.
    let gap = best_drop.0.abs_diff(deepest.0);
    assert!(
        gap <= 1,
        "the resistance dip is {gap} steps from the drop peak ({} mA vs {} mA)",
        I_SWEEP[deepest.0],
        I_SWEEP[best_drop.0]
    );
    assert!(
        (-deepest.1 - 300.0).abs() < 80.0,
        "on resonance at 1 mW the paper sees about -300 ohm; got {:.1}",
        deepest.1
    );
}

/// Undope the waveguide and the readout has to go away.
///
/// This is the control. The two tests above are also passed by a ring that
/// simply gets warm, so switch off the mechanism the effect is attributed to —
/// the free carriers the doping put in the mode — and require it to disappear.
/// What is left is the scattering loss, which absorbs nothing into the silicon.
/// Switch off the absorption that makes carriers and the sign has to INVERT.
///
/// This is the control, and it is stronger than "the effect gets smaller". With
/// `alpha_ssa_db_cm = 0` the free-carrier absorption and the scattering are
/// still there, so the ring still absorbs and still warms — but nothing frees a
/// carrier, so only the thermal half of the loop survives and the resistance
/// must go UP. If the dip persists, it is not coming from the carriers.
#[test]
fn with_no_surface_absorption_the_light_raises_the_resistance_instead() {
    if !common::have_compiler() {
        return;
    }
    let at_res = |extra: &str| {
        let (r_lit, _, _) = point(0.77, 1.0, extra);
        let (r_dark, _, _) = point(0.77, 1e-7, extra);
        r_lit - r_dark
    };
    let with = at_res("");
    let without = at_res("alpha_ssa_db_cm=0");
    assert!(
        with < -100.0,
        "with surface absorption the light should pull R well down; got {with:.1} ohm"
    );
    assert!(
        without > 0.0,
        "with the carrier channel off only heating is left, so R must RISE; got \
         {without:.1} ohm"
    );
    assert!(
        without < -with / 20.0,
        "the thermal half should be far smaller than the carrier half: {without:.1} \
         against {with:.1} ohm"
    );
}

/// Current drive amplifies the readout; voltage drive fights it.
///
/// Same ring, same light, same operating point — only the source type differs.
/// At constant current a falling R dissipates less, cools, and falls further;
/// at constant voltage it dissipates more, warms, and is pulled back. So the
/// fractional signal has to be bigger in current mode, and this is the reason
/// the rest of this file drives current.
#[test]
fn current_drive_amplifies_the_readout_and_voltage_drive_suppresses_it() {
    if !common::have_compiler() {
        return;
    }
    // Current mode at the resonance, and the voltage the DARK ring develops
    // there — so the voltage-mode run starts from the same operating point.
    let (r_dark_i, _, _) = point(0.77, 1e-7, "");
    let (r_lit_i, _, _) = point(0.77, 1.0, "");
    let v_op = r_dark_i * 0.77e-3;

    let r_of = |p_mw: f64| {
        let r = solve(&deck_v(v_op, p_mw, ""));
        let i = -r.vsrc_current("vh").unwrap();
        v_op / i
    };
    let (r_dark_v, r_lit_v) = (r_of(1e-7), r_of(1.0));

    let frac_i = (r_lit_i - r_dark_i) / r_dark_i;
    let frac_v = (r_lit_v - r_dark_v) / r_dark_v;
    assert!(
        frac_i < 0.0 && frac_v < 0.0,
        "both drives should still show a drop: {frac_i:.4} current, {frac_v:.4} voltage"
    );
    assert!(
        frac_i < frac_v * 1.15,
        "current drive should give the larger drop by a clear margin: {:.2} % against \
         {:.2} %",
        frac_i * 100.0,
        frac_v * 100.0
    );
}

/// Q and the FSR, because between them they set two parameters.
///
/// `n_g` comes from the free spectral range at a fixed 8 um radius —
/// `FSR = lambda^2/(n_g*L)` leaves nothing else to give — and `kappa_l` comes
/// from the loaded Q, because the linewidth of this ring is 31:1
/// coupling-dominated and the loss channels cannot reach it. Neither was fitted
/// alongside anything else, so if either drifts, the parameter that was solved
/// from it is wrong rather than merely stale.
#[test]
fn the_linewidth_and_free_spectral_range_are_the_ones_they_were_solved_from() {
    if !common::have_compiler() {
        return;
    }
    // Dark and cold: a nanowatt moves nothing, so this is the passive cavity.
    // The wavelength belongs to the source, not the ring, so it is swept by
    // rewriting the deck rather than through the instance-parameter string.
    let drop = |nm: f64| {
        let d = deck(1e-3, 1e-6, "").replace("wavelength_nm=1550", &format!("wavelength_nm={nm}"));
        let r = solve(&d);
        let re = r.node_voltage("dr_re_0").unwrap();
        let im = r.node_voltage("dr_im_0").unwrap();
        re * re + im * im
    };
    /// Resonance and FWHM (nm) of the strongest notch in `[lo, hi]`.
    fn peak(f: &dyn Fn(f64) -> f64, lo: f64, hi: f64, n: usize) -> (f64, f64) {
        let step = (hi - lo) / (n - 1) as f64;
        let y: Vec<f64> = (0..n).map(|i| f(lo + step * i as f64)).collect();
        let i = y
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .unwrap()
            .0;
        let half = y[i] / 2.0;
        let mut a = i;
        let mut b = i;
        while a > 0 && y[a] > half {
            a -= 1;
        }
        while b < n - 1 && y[b] > half {
            b += 1;
        }
        (lo + step * i as f64, step * (b - a) as f64)
    }

    let (res, fwhm) = peak(&drop, 1549.2, 1550.1, 46);
    let q = res / fwhm;
    assert!(
        (q - 5900.0).abs() / 5900.0 < 0.12,
        "loaded Q should be the 5900 kappa_l was solved from; got {q:.0} \
         (resonance {res:.3} nm, FWHM {:.0} pm)",
        fwhm * 1e3
    );
    let (next, _) = peak(&drop, 1560.5, 1563.5, 31);
    let fsr = next - res;
    assert!(
        (fsr - 12.1).abs() < 0.35,
        "the FSR should be the 12.1 nm n_g was solved from; got {fsr:.2} nm"
    );
}
