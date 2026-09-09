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

/// A fixed-step run keeps its output grid **and** resolves the delay, by taking
/// an integer number of internal steps per requested point.
///
/// The grid is a promise: shrinking the step outright would move every output
/// time. An integer factor moves none of them — each requested time is still
/// landed on exactly, with no interpolation — and the delay gets the resolution
/// it asked for. The cost is reported once rather than hidden.
///
/// This used to refuse, which was worse: it handed the user an arithmetic
/// problem the device had already solved.
#[test]
fn a_fixed_step_run_substeps_rather_than_losing_the_delay() {
    let net = "* coarse step on a 1 ns line\n\
               Vs s 0 PULSE(0 1 0.5n 10p 10p 100n 200n)\n\
               Rs s a 50\n\
               T1 a 0 b 0 Z0=50 TD=1n\n\
               Rterm b 0 1Meg\n";
    let parsed = fairchild_parser::parse_spice(net).expect("parse");
    // 2 ns steps on a 1 ns line: four internal steps per output point.
    let coarse = fairchild_core::tran_nr(&parsed, 2e-9, 12e-9).expect("must not refuse");

    // The output grid is exactly the one asked for.
    for (k, t) in coarse.time.iter().enumerate().take(6) {
        let want = k as f64 * 2e-9;
        assert!(
            (t - want).abs() < 1e-15,
            "output point {k} landed at {t:.4e}, the card asked for {want:.4e}"
        );
    }

    // And the physics is the physics: the far end doubles one TD after the
    // launch at 0.5 ns, so by the 2 ns sample it is up. Before sub-stepping
    // this read 0 here and rose at 4 ns — the line behaving as TD = 2 ns.
    let v_b = coarse.voltage_at("b", 2e-9).expect("node b");
    assert!(
        v_b > 0.9,
        "far end should be up by t = 2 ns (launch 0.5 ns + TD = 1 ns), got {v_b:.4}"
    );
    // The near-end reflection arrives at 2·TD after the launch, so it is up by
    // the 4 ns sample and not by the 2 ns one.
    let v_a_early = coarse.voltage_at("a", 2e-9).expect("node a");
    let v_a_late = coarse.voltage_at("a", 4e-9).expect("node a");
    assert!(
        v_a_early < 0.6 && v_a_late > 0.9,
        "near end should still be at the launched half-step at 2 ns and up at 4 ns, \
         got {v_a_early:.4} then {v_a_late:.4}"
    );

    // A step the delay can carry needs no sub-stepping and must be unchanged.
    let fine = fairchild_core::tran_nr(&parsed, 20e-12, 12e-9).expect("20 ps step");
    let fine_b = fine.voltage_at("b", 2e-9).expect("node b");
    assert!(
        (v_b - fine_b).abs() < 0.02,
        "the sub-stepped answer must agree with a natively fine run: \
         {v_b:.4} vs {fine_b:.4}"
    );
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

// ── the part of the error a tolerance decides (#112 part 2) ──────────────────

/// A line between two ideal sources: every node is forced, so the LTE estimate
/// has **nothing** to look at.
///
/// The estimate covers node rows and excludes the ones a source pins. Here that
/// is all of them, and the only unknowns left are the two port currents, which
/// the norm never reaches. Without a bound the controller sees an error of
/// identically zero and runs at the card's step, whatever the delayed wave is
/// doing.
///
/// The line is lossy so the two sources are not a DC short of each other, which
/// would be a voltage-source loop and is correctly refused.
fn two_source_line(freq: f64) -> String {
    format!(
        "* line between two ideal sources\n\
         Vs s 0 SIN(0 1 {freq:e} 0 0 0)\n\
         Vt b 0 DC 0\n\
         T1 s 0 b 0 Z0=50 TD=0.7n loss_db=3\n\
         .options variable_step=1\n"
    )
}

#[test]
fn a_delay_bounds_the_step_where_the_lte_estimate_is_blind() {
    let parsed = fairchild_parser::parse_spice(&two_source_line(1e9)).expect("parse");
    // The card asks for 100 ps and TD/2 allows 350, so nothing else in the
    // controller would go below 100 ps here.
    let res = fairchild_core::tran_nr_var(&parsed, 100e-12, 4e-9).expect("transient");
    let steps: Vec<f64> = res.time.windows(2).map(|w| w[1] - w[0]).collect();
    // Skip the opening steps: the bound needs three samples before it can
    // estimate a curvature at all, so the first two are unconstrained (#120).
    let settled = &steps[3..];
    let worst = settled.iter().copied().fold(0.0f64, f64::max);
    assert!(
        worst < 20e-12,
        "with every node forced, only the delay's own curvature can bound the \
         step; largest settled step was {:.2} ps against a 100 ps card",
        worst * 1e12
    );
}

/// And the bound is the formula, not a constant.
///
/// `h ≤ √(8·tol/|y''|)` with `y'' ∝ ω²`, so the step goes as `1/ω`: doubling
/// the source frequency roughly halves it. That scaling is a property of the
/// expression rather than of this circuit, which is what makes it an anchor —
/// a hard-coded step limit would pass the test above and fail this one, and a
/// `1/ω²` law would overshoot it by the same margin.
///
/// "Roughly", because `tol = vntol + reltol·|x|` moves with the amplitude at
/// the port, and the reflection pattern of a fixed `TD` is not the same at the
/// two frequencies. Measured 2.1 against an ideal 2.0; the window below admits
/// nothing that is not a `1/ω` law.
#[test]
fn the_curvature_bound_scales_as_one_over_frequency() {
    let settled_max = |freq: f64| {
        let parsed = fairchild_parser::parse_spice(&two_source_line(freq)).expect("parse");
        let res = fairchild_core::tran_nr_var(&parsed, 100e-12, 3e-9).expect("transient");
        let mut steps: Vec<f64> = res.time[3..].windows(2).map(|w| w[1] - w[0]).collect();
        // The median, not the extremes. A sine's curvature passes through zero
        // twice a cycle, where this bound goes to infinity and something else
        // limits the step, so the largest step says more about the other
        // limiter than about this one.
        steps.sort_by(|a, b| a.partial_cmp(b).unwrap());
        steps[steps.len() / 2]
    };
    let (lo, hi) = (settled_max(1e9), settled_max(2e9));
    let ratio = lo / hi;
    assert!(
        (1.5..2.6).contains(&ratio),
        "doubling the frequency should roughly halve the bound: {:.3} ps at \
         1 GHz and {:.3} ps at 2 GHz is a ratio of {ratio:.2}. A fixed limit \
         would give 1 and a 1/omega^2 law would give 4",
        lo * 1e12,
        hi * 1e12
    );
}
