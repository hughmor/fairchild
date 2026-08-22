//! The worked WDM example, exercised: `examples/verilog_a/models/mrm_wdm.va`.
//!
//! A silicon microring modulator written in the bundle-port dialect, with the
//! channel count nowhere in the file. This is the model a user is pointed at, so
//! it is worth proving it does what its header claims rather than merely
//! compiling — a broken example teaches the wrong thing more effectively than no
//! example at all.
//!
//! What each test pins is a *coupling*, because coupling is what makes the model
//! more than a transfer function: the ring's resonance moves with carriers the
//! junction injects, and with heat the heater puts in. Both are solved, not
//! assumed — carrier density and temperature are internal nodes with real MNA
//! rows.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_osdi::{load_libraries_with_widths, VaOptions};
use fairchild_parser::{instantiated_widths, parse_spice_with_arity, PermissiveArity};
use std::collections::{BTreeMap, BTreeSet};

mod common;

fn model_dir() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/verilog_a/models")
}

/// Build the deck's registry, generating `mrm_wdm` at whatever widths it uses.
fn solve(deck: &str) -> fairchild_core::NrResult {
    let probe = parse_spice_with_arity(deck, &PermissiveArity).expect("probe pass parses");
    let widths: BTreeMap<String, BTreeSet<usize>> = instantiated_widths(&probe);
    let mut reg = DeviceRegistry::new();
    load_libraries_with_widths(
        &[],
        &[model_dir().join("mrm_wdm.va").to_string_lossy().to_string()],
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

/// A resonance of this ring, so the deck actually exercises one.
///
/// `λ_res = n_eff·L/m`. With the defaults — n_eff 2.3022, r = 8 µm — the C-band
/// order is m = 75, giving 1542.93 nm. Off resonance an all-pass ring is
/// all-pass: unity transmission, phase only, and every test would read 0.4999 mW
/// whatever the bias did. The first version of this file made exactly that
/// mistake.
const LAMBDA_RES_NM: f64 = 1542.93;

/// Channel spacing, ~0.7 of a linewidth. FSR = λ²/(n_g·L) = 11.3 nm and the
/// finesse is ≈ 69, so a linewidth is ≈ 0.16 nm: this spreads the comb across
/// the resonance instead of bunching it on one flank.
const SPACING_NM: f64 = 0.12;

/// `n` lasers muxed onto one bus through the ring, with the given bias and
/// heater drive. `extra` goes on the instance line, for switching parts of the
/// physics off.
fn deck_with(n: usize, v_bias: f64, v_htr: f64, extra: &str) -> String {
    let mut s = String::new();
    for k in 0..n {
        s += &format!(".optical_port c{k}\n");
    }
    s += &format!(".optical_port bus {n}\n.optical_port out {n}\n");
    for k in 0..n {
        // Straddle the resonance: channel 0 sits on it, the rest walk off.
        s += &format!(
            "Xl{k} c{k} fc_cw_laser power_mW=0.5 wavelength_nm={:.4}\n",
            LAMBDA_RES_NM + SPACING_NM * k as f64
        );
    }
    s += "Xmux bus";
    for k in 0..n {
        s += &format!(" c{k}");
    }
    s += " fc_mux\n";
    s += &format!(
        "Xr bus out a 0 hp 0 tr mrm_wdm {extra}\n\
         Vb a 0 DC {v_bias}\n\
         Vh hp 0 DC {v_htr}\n\
         .op\n"
    );
    s
}

fn deck(n: usize, v_bias: f64, v_htr: f64) -> String {
    deck_with(n, v_bias, v_htr, "")
}

/// The same deck with everything that couples channels together switched off:
/// no free-carrier absorption, no TPA, and a thermal resistance small enough
/// that absorbed power raises no temperature worth having.
fn deck_uncoupled(n: usize) -> String {
    deck_with(n, 0.0, 0.0, "dalpha_dnc=0 beta_tpa=0 r_th=1e-9")
}

fn power_mw(r: &fairchild_core::NrResult, ch: usize) -> f64 {
    let re = r.node_voltage(&format!("out_re_{ch}")).unwrap();
    let im = r.node_voltage(&format!("out_im_{ch}")).unwrap();
    (re * re + im * im) * 1e3
}

/// One source file, two channel counts, with the channels made independent.
///
/// With free-carrier absorption, TPA and self-heating switched off, nothing
/// couples one channel to another, so channel 0's answer must be identical
/// whether it runs alone or beside three more. This is what catches a
/// per-channel indexing error — the failure mode a single-channel test cannot
/// see, because it needs two channels to be confused with each other.
#[test]
fn an_uncoupled_channel_is_unchanged_by_its_neighbours() {
    if !common::have_compiler() {
        return;
    }
    let one = power_mw(&solve(&deck_uncoupled(1)), 0);
    let four = power_mw(&solve(&deck_uncoupled(4)), 0);
    assert!(
        (one - four).abs() < 1e-12,
        "with the shared physics off, channel 0 gives {one:.12} mW alone and \
         {four:.12} mW beside three others — they should be identical"
    );
}

/// And with the shared physics back on, the neighbours DO matter.
///
/// One ring, one carrier pool, one temperature: light in the other channels is
/// absorbed by the same waveguide and heats the same resonator, which moves the
/// resonance channel 0 is sitting on. A model where this made no difference
/// would not be modelling a device, only three independent filters.
#[test]
fn a_shared_ring_couples_its_channels_through_heat_and_carriers() {
    if !common::have_compiler() {
        return;
    }
    let one = power_mw(&solve(&deck(1, 0.0, 0.0)), 0);
    let four = power_mw(&solve(&deck(4, 0.0, 0.0)), 0);
    assert!(
        (one - four).abs() > 1e-9,
        "channel 0 reads {one:.9} mW alone and {four:.9} mW with three neighbours \
         — absorbed power from the others shares one ring and should move it"
    );
}

/// The four channels sit at four detunings across the resonance, so they must
/// come out at four different powers. A model that ignored `LAMBDA(bus_in, k)`
/// and used one wavelength for the bundle would return four identical numbers.
#[test]
fn each_channel_sees_its_own_detuning() {
    if !common::have_compiler() {
        return;
    }
    let r = solve(&deck(4, 0.0, 0.0));
    let p: Vec<f64> = (0..4).map(|k| power_mw(&r, k)).collect();
    let spread =
        p.iter().cloned().fold(f64::MIN, f64::max) - p.iter().cloned().fold(f64::MAX, f64::min);
    assert!(
        spread > 1e-3,
        "four wavelengths across one resonance should give four transmissions, \
         got {p:?}"
    );
}

/// Forward bias injects carriers, carriers lower the index, the resonance moves
/// — so the transmitted power changes. This is the modulator working, and it is
/// the whole point of the junction being a transport model rather than a
/// constant.
#[test]
fn forward_bias_injects_carriers_and_moves_the_resonance() {
    if !common::have_compiler() {
        return;
    }
    let off = power_mw(&solve(&deck(2, 0.0, 0.0)), 0);
    let on = power_mw(&solve(&deck(2, 0.40, 0.0)), 0);
    assert!(
        (on - off).abs() > 1e-4,
        "0.40 V of forward bias should move the resonance and change the \
         transmission: {off:.6} mW unbiased vs {on:.6} mW biased"
    );
}

/// The heater raises the solved temperature, which shifts the index the other
/// way. That the *heater* and the *junction* both move it, through different
/// physics onto one shared resonance, is what makes this a device rather than
/// two independent effects.
#[test]
fn the_heater_shifts_the_resonance_through_a_solved_temperature() {
    if !common::have_compiler() {
        return;
    }
    let cold = power_mw(&solve(&deck(2, 0.0, 0.0)), 0);
    let warm = power_mw(&solve(&deck(2, 0.0, 0.20)), 0);
    assert!(
        (warm - cold).abs() > 1e-4,
        "0.20 V across the heater should warm the ring and move its resonance: \
         {cold:.6} mW cold vs {warm:.6} mW warm"
    );
}
