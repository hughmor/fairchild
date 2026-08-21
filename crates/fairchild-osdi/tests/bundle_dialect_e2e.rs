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
const WL_NM: f64 = 1550.0;

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
    }
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
