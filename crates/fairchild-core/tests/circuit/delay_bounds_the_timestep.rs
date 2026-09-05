//! Every device that owns a delay must bound the timestep, and the list of
//! those devices must not be able to grow silently.
//!
//! `DelayLine::sample` clamps above its newest sample, so a step longer than
//! the delay reconstructs it from the previous accepted point: the effective
//! delay becomes `max(tau, h)` and tracks the step size rather than the
//! geometry. Nothing in the LTE norm can see that — it measures the error of a
//! step already taken, and this is an error in what the circuit *is*.
//!
//! The expensive failure here is an absence: a new delay device that forgets
//! `requested_max_timestep` looks exactly like a correct one, and every
//! existing test still passes. So this file carries a completeness gate as well
//! as the behaviour, and the gate reads the source rather than a hand-kept
//! list.

use std::collections::BTreeSet;
use std::path::Path;

use fairchild_core::device::{Device, SimContext};
use fairchild_core::models::{NativeTLine, NativeWaveguide};

/// Devices known to own a `DelayLine`, by the file that declares the field.
///
/// Adding a delay to a device means adding it here *and* giving it a
/// `requested_max_timestep`. The gate below fails on a file that has one and is
/// not listed, which is the only way to notice the omission.
const KNOWN_DELAY_OWNERS: &[&str] = &[
    "models/tline.rs",            // NativeTLine — TD
    "models/photonic/segment.rs", // OpticalSegment — tau_g, under waveguide_delay
    "models/photonic/xfer.rs",    // NativeOptical2x2 — tau_s
];

/// Walk `src/` and report every file declaring a `DelayLine` field.
fn files_owning_a_delay_line() -> BTreeSet<String> {
    let src = Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
    let mut found = BTreeSet::new();
    let mut stack = vec![src.clone()];
    while let Some(dir) = stack.pop() {
        for entry in std::fs::read_dir(&dir).expect("read src") {
            let path = entry.expect("dir entry").path();
            if path.is_dir() {
                stack.push(path);
                continue;
            }
            if path.extension().is_none_or(|e| e != "rs") {
                continue;
            }
            let text = std::fs::read_to_string(&path).expect("read source");
            // The field declaration, not the `use` line or a doc mention.
            if text.contains(": DelayLine,") {
                let rel = path
                    .strip_prefix(&src)
                    .expect("under src")
                    .to_string_lossy()
                    .replace('\\', "/");
                found.insert(rel);
            }
        }
    }
    found
}

#[test]
fn every_delay_owner_is_accounted_for() {
    let found = files_owning_a_delay_line();
    let known: BTreeSet<String> = KNOWN_DELAY_OWNERS.iter().map(|s| s.to_string()).collect();
    let unlisted: Vec<_> = found.difference(&known).collect();
    assert!(
        unlisted.is_empty(),
        "these files declare a DelayLine and are not in KNOWN_DELAY_OWNERS: {unlisted:?}. \
         Add each one here and give its device a `requested_max_timestep`, or a step \
         longer than its delay will silently become its delay."
    );
    let vanished: Vec<_> = known.difference(&found).collect();
    assert!(
        vanished.is_empty(),
        "KNOWN_DELAY_OWNERS lists files with no DelayLine any more: {vanished:?}. \
         Remove them so the gate keeps meaning something."
    );
}

#[test]
fn a_transmission_line_bounds_the_step_to_half_its_delay() {
    let t = NativeTLine::new(50.0, 1e-9);
    assert_eq!(
        t.requested_max_timestep(),
        Some(0.5e-9),
        "a T element must ask for TD/2"
    );
    // TD = 0 is a wire, not a delay, and must not pin the step to zero.
    let wire = NativeTLine::new(50.0, 0.0);
    assert_eq!(wire.requested_max_timestep(), None);
}

/// A waveguide asks only when the option that engages its delay is on. This is
/// the on/off pair `CLAUDE.md` asks for: a bound that appeared unconditionally
/// would throttle every photonic run for a delay it is not modelling.
#[test]
fn a_waveguide_bounds_the_step_only_when_the_delay_is_engaged() {
    let build = |delay: bool| {
        let ctx = SimContext {
            waveguide_delay: delay,
            ..Default::default()
        };
        let mut wg = NativeWaveguide::new();
        wg.setup_model(&ctx);
        // 1 cm at n_g = 4.19 (the default strip) => tau_g ~ 140 ps.
        wg.set_real_param("l_m", 1e-2);
        wg.setup_model(&ctx);
        wg
    };
    assert_eq!(
        build(false).requested_max_timestep(),
        None,
        "delay off: no bound, or every photonic run pays for a delay it is not modelling"
    );
    let bound = build(true)
        .requested_max_timestep()
        .expect("delay on: a bound");
    let tau_g = 1e-2 * 4.19 / 299_792_458.0;
    assert!(
        (bound - tau_g / 2.0).abs() < 1e-18,
        "expected tau_g/2 = {:.4e}, got {bound:.4e}",
        tau_g / 2.0
    );
}

/// A fixed-step run promises a sample grid, so it cannot honour a smaller
/// requested step by sub-stepping. It must refuse instead: the answer at the
/// requested size is a circuit with a different delay in it.
#[test]
fn a_fixed_step_run_refuses_a_step_longer_than_a_delay() {
    let net = "* coarse step on a 1 ns line\n\
               Vs s 0 PULSE(0 1 0.5n 10p 10p 100n 200n)\n\
               Rs s a 50\n\
               T1 a 0 b 0 Z0=50 TD=1n\n\
               Rterm b 0 1Meg\n";
    let parsed = fairchild_parser::parse_spice(net).expect("parse");
    let msg = match fairchild_core::tran_nr(&parsed, 2e-9, 12e-9) {
        Err(e) => format!("{e}"),
        Ok(_) => panic!("a 2 ns step on a 1 ns line must not be answered"),
    };
    assert!(
        msg.contains("too large for a delay"),
        "the error must name the delay as the reason, got: {msg}"
    );
    // And the same deck at a step the delay can carry runs.
    fairchild_core::tran_nr(&parsed, 20e-12, 12e-9).expect("20 ps step is fine");
}

/// The variable-step controller honours the bound instead of refusing, because
/// there the step is its to choose. The `.tran` card asks for 2 ns; every
/// accepted step must come in under TD/2.
#[test]
fn the_variable_step_controller_honours_the_bound_from_the_first_step() {
    let net = "* coarse card, adaptive stepping, 1 ns line\n\
               Vs s 0 PULSE(0 1 0.5n 10p 10p 100n 200n)\n\
               Rs s a 50\n\
               T1 a 0 b 0 Z0=50 TD=1n\n\
               Rterm b 0 1Meg\n\
               .options variable_step=1\n";
    let parsed = fairchild_parser::parse_spice(net).expect("parse");
    let res = fairchild_core::tran_nr_var(&parsed, 2e-9, 6e-9).expect("transient");
    let worst = res
        .time
        .windows(2)
        .map(|w| w[1] - w[0])
        .fold(0.0f64, f64::max);
    assert!(
        worst <= 0.5e-9 + 1e-15,
        "largest accepted step {worst:.3e} s exceeds TD/2; the first step is the \
         one that escapes if the bound is only applied after an acceptance"
    );
    // The physics the bound exists to protect: the far end doubles one TD after
    // the launch, not one timestep after it.
    let v_b = res.voltage_at("b", 1.6e-9).expect("node b");
    assert!(
        v_b > 0.9,
        "far end should have doubled to ~1 V by t = 1.6 ns (launch 0.5 ns + TD), got {v_b:.4}"
    );
    let v_b_early = res.voltage_at("b", 1.2e-9).expect("node b");
    assert!(
        v_b_early < 0.1,
        "far end must still be quiet at t = 1.2 ns, got {v_b_early:.4}"
    );
}
