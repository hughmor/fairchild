//! One Verilog-A source, three channel counts, checked against native (#55).
//!
//! `wg_bundle.va` never mentions a channel count. These runs compile it at
//! N = 1, 2 and 4 from that one file and check each against `fc_waveguide`,
//! which is the independent anchor — comparing a generated model to a second
//! generated model would share any expansion fault, and comparing it to
//! `wg_wdm2.va` would only cover N = 2.
//!
//! The N = 1 vs N = 4 pair is the case that matters most. A per-channel indexing
//! error shows up as one channel carrying another channel's field, and a
//! single-channel test cannot see it — so every channel is launched at a
//! different power and checked against its own.

use std::collections::{BTreeMap, BTreeSet};

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_osdi::{load_libraries_with_widths, VaOptions};
use fairchild_parser::{instantiated_widths, parse_spice_with_arity, PermissiveArity};

mod common;

const L_UM: f64 = 1000.0;
const N_G: f64 = 4.2;
const ALPHA: f64 = 3.0;
/// Deliberately *not* the 1550 nm band centre. `wl_default` in the generated
/// source is 1550 nm, and an unreached λ port falls back to the band centre, so
/// a wavelength that never reached the model would land on 1550 either way. At
/// 1310 nm the phase through 1 mm differs by thousands of radians, which shows
/// up as a different `(re, im)` split against the native anchor.
const WL_NM: f64 = 1310.0;

/// A deck putting `model` on an n-channel bundle beside a native waveguide fed
/// from the same lasers, so the two can be compared channel by channel.
fn deck(model: &str, powers: &[f64]) -> String {
    let n = powers.len();
    let mut s = String::new();
    for k in 0..n {
        s += &format!(".optical_port c{k}\n");
    }
    s += &format!(".optical_port bus {n}\n.optical_port dut {n}\n.optical_port ref {n}\n");
    for (k, p) in powers.iter().enumerate() {
        s += &format!("Xl{k} c{k} fc_cw_laser power_mW={p} wavelength_nm={WL_NM}\n");
    }
    s += "Xmux bus";
    for k in 0..n {
        s += &format!(" c{k}");
    }
    s += " fc_mux\n";
    // The model under test, and the native anchor, from the same bus.
    s += &format!(
        "Xdut bus dut {model} l_um={L_UM} n_g={N_G} alpha_dB_cm={ALPHA}\n\
         Xref bus ref fc_waveguide L_um={L_UM} n_g={N_G} alpha_dB_cm={ALPHA} wavelength_nm={WL_NM}\n\
         .op\n.end\n"
    );
    s
}

fn power_of(r: &fairchild_core::NrResult, port: &str, ch: usize) -> f64 {
    let re = r.node_voltage(&format!("{port}_re_{ch}")).unwrap();
    let im = r.node_voltage(&format!("{port}_im_{ch}")).unwrap();
    re * re + im * im
}

fn field_of(r: &fairchild_core::NrResult, port: &str, ch: usize) -> (f64, f64) {
    (
        r.node_voltage(&format!("{port}_re_{ch}")).unwrap(),
        r.node_voltage(&format!("{port}_im_{ch}")).unwrap(),
    )
}

/// Build a registry with `wg_bundle.va` generated for exactly the widths this
/// deck asks for — which is what the two-pass load does for real.
fn registry_for(src: &str) -> (DeviceRegistry, fairchild_parser::Netlist) {
    let probe = parse_spice_with_arity(src, &PermissiveArity).expect("probe pass parses");
    let widths: BTreeMap<String, BTreeSet<usize>> = instantiated_widths(&probe);

    let va = PathBuf_of("wg_bundle.va");
    let mut reg = DeviceRegistry::new();
    load_libraries_with_widths(
        &[],
        &[va.to_string_lossy().to_string()],
        None,
        &VaOptions {
            cache_dir: Some(common::test_cache_dir()),
            ..VaOptions::from_env()
        },
        &mut reg,
        &widths,
        3,
    )
    .expect("the dialect source generates and compiles");
    let net = parse_spice_with_arity(src, &reg).expect("a bundle model takes the declared width");
    (reg, net)
}

#[allow(non_snake_case)]
fn PathBuf_of(name: &str) -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/models")
        .join(name)
}

fn run_at(powers: &[f64]) {
    let n = powers.len();
    let src = deck("wg_bundle", powers);
    let (reg, net) = registry_for(&src);
    let r = dc_op_nr_with_registry(&net, &reg).expect("DC OP converges");

    for (ch, &launched_mw) in powers.iter().enumerate() {
        let dut = power_of(&r, "dut", ch);
        let anchor = power_of(&r, "ref", ch);
        assert!(
            (dut - anchor).abs() < 1e-12,
            "N={n} channel {ch}: generated model gives {dut:.12} W, native \
             fc_waveguide gives {anchor:.12} W"
        );
        // And the channel really is its own: the launched powers differ, so a
        // crossed index would show up as another channel's number.
        let launched = launched_mw * 1e-3;
        let expect = launched * 10f64.powf(-ALPHA * (L_UM * 1e-4) / 10.0);
        assert!(
            (dut - expect).abs() < 1e-12,
            "N={n} channel {ch}: got {dut:.12} W, hand-computed {expect:.12} W \
             from its own launched {launched:.6} W"
        );
        // Power alone cannot see the wavelength — loss does not depend on it —
        // and the native anchor cannot arbitrate the phase either, because
        // `wg_bundle.va` propagates on `n_g` while `fc_waveguide` uses a
        // dispersion-corrected `n_eff`. So the phase is checked against the
        // model's own law, `A_out = A_in·amp·exp(−j·2π·n_g·L/λ)`, computed here
        // from the laser's parameter. This is the half that fails if `wl_<k>` is
        // not filled from resolution: at the 1550 nm `wl_default` the field
        // comes out (−0.0356, +0.0726) instead of (+0.0633, −0.0503) on the
        // 7 mW channel — 3000 rad of accumulated phase apart.
        let ph = 2.0 * std::f64::consts::PI * N_G * (L_UM * 1e-6) / (WL_NM * 1e-9);
        let amp = 10f64.powf(-ALPHA * (L_UM * 1e-4) / 20.0);
        let a_in = launched.sqrt();
        let (want_re, want_im) = (amp * a_in * ph.cos(), -amp * a_in * ph.sin());
        let (dut_re, dut_im) = field_of(&r, "dut", ch);
        assert!(
            (dut_re - want_re).abs() < 1e-9 && (dut_im - want_im).abs() < 1e-9,
            "N={n} channel {ch}: generated model gives ({dut_re:.9}, {dut_im:.9}), its \
             own law at the laser's {WL_NM} nm gives ({want_re:.9}, {want_im:.9}) — the \
             wavelength it evaluated its phase at is not the one the deck emits"
        );
    }
}

/// A label has to reach what comes *after* a generated model.
///
/// The dialect used to make this the author's job: `OWL(WL(b,k)) <+ OWL(WL(a,k))`
/// carried the tag through the matrix. λ is not a matrix quantity any more, so
/// a bundle model instead *declares* its λ terminals and its slot-for-slot
/// routing, and resolution walks through it. Left undeclared, everything
/// downstream resolves to the band centre — the wrong wavelength with no
/// diagnostic, which is the failure this checks for.
///
/// The anchor is the native waveguide's own dispersion-corrected law at the
/// laser's 1310 nm, evaluated here. At the 1550 nm band centre the phase
/// increment differs by hundreds of radians.
#[test]
fn a_label_reaches_the_native_device_downstream_of_a_generated_one() {
    if !common::have_compiler() {
        return;
    }
    const N_EFF: f64 = 2.445;
    const TAIL_UM: f64 = 40.0;

    let mut src = deck("wg_bundle", &[7.0]);
    src = src.replace(
        ".op\n.end\n",
        &format!(
            ".optical_port tail\n\
             Xtail dut tail fc_waveguide L_um={TAIL_UM} n_eff={N_EFF} n_g={N_G} alpha_dB_cm=0\n\
             .op\n.end\n"
        ),
    );
    let (reg, net) = registry_for(&src);
    let r = dc_op_nr_with_registry(&net, &reg).expect("DC OP converges");

    let (a_re, a_im) = field_of(&r, "dut", 0);
    let (b_re, b_im) = field_of(&r, "tail", 0);
    // The native segment linearises dispersion about `wl_ref_m`, which defaults
    // to the band centre.
    let lambda = WL_NM * 1e-9;
    let wl_ref = 1.55e-6;
    let n_eff_lam = N_EFF + (lambda - wl_ref) * (N_EFF - N_G) / wl_ref;
    let want = 2.0 * std::f64::consts::PI * n_eff_lam * (TAIL_UM * 1e-6) / lambda;
    let got = -((b_im * a_re - b_re * a_im).atan2(b_re * a_re + b_im * a_im));
    let two_pi = 2.0 * std::f64::consts::PI;
    let diff = (got - want).rem_euclid(two_pi);
    let diff = diff.min(two_pi - diff);
    assert!(
        diff < 1e-6,
        "the waveguide after the generated model turned the field by {got:.6} rad; \
         the laser's {WL_NM} nm implies {want:.6} rad (mod 2π) — the label did not \
         propagate through the Verilog-A bundle model"
    );
}

#[test]
fn one_source_serves_one_channel() {
    if !common::have_compiler() {
        return;
    }
    run_at(&[7.0]);
}

#[test]
fn the_same_source_serves_two_channels() {
    if !common::have_compiler() {
        return;
    }
    run_at(&[7.0, 3.0]);
}

/// Four channels at four different powers — the case a single-channel test
/// cannot cover, because a per-channel indexing error is invisible until the
/// channels differ.
#[test]
fn the_same_source_serves_four_channels_and_keeps_them_apart() {
    if !common::have_compiler() {
        return;
    }
    run_at(&[9.0, 5.0, 2.0, 11.0]);
}
