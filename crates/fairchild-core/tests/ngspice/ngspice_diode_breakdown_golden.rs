//! Reverse breakdown (`BV`/`IBV`) against ngspice, at runtime.
//!
//! #77 §3 called this the gap most likely to answer a real question wrongly, and
//! that is the shape of it: a Zener and an ESD clamp are *specified* by their
//! knee. Without one, both simulate as ordinary diodes that block whatever you
//! put across them, the deck runs, and the number is wrong by orders of
//! magnitude with no diagnostic.
//!
//! # Why runtime comparison and not a stored golden
//!
//! The knee is an exponential, so every point on it is a statement about the
//! *offset* of that exponential. A stored number would pin one point and let the
//! shape drift. Running ngspice pins the shape, and the sweep below crosses the
//! knee — flat leakage, the exponential tail, the knee itself, and runaway.
//!
//! Requires ngspice on PATH; skipped, not failed, without it.

use std::process::Command;

use fairchild_core::{dc_op_nr_with_registry_opts, options::SimOptions, DeviceRegistry};
use fairchild_parser::parse_spice;

/// `D1 0 a` with `V1 a 0 DC v` puts the diode in **reverse** for positive `v`,
/// so the sweep below reads as a reverse-bias sweep.
fn deck(model: &str, v: f64) -> String {
    format!("* breakdown\n.options gmin=0\n.model dm D ({model})\nV1 a 0 DC {v}\nD1 0 a dm\n.op\n")
}

fn fairchild(model: &str, v: f64) -> f64 {
    let src = deck(model, v);
    let net = parse_spice(&src).expect("parse");
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);
    // `from_netlist`, not `default()`: the deck sets `gmin=0` so the comparison
    // is about the junction and not about a leakage floor either side.
    let opts = SimOptions::from_netlist(&net);
    dc_op_nr_with_registry_opts(&net, &reg, &opts)
        .unwrap_or_else(|e| panic!("fairchild failed on\n{src}\n{e:?}"))
        .vsrc_current("v1")
        .expect("I(v1)")
}

fn ngspice(model: &str, v: f64) -> Option<f64> {
    let dir = std::env::temp_dir().join("fc_bv_golden");
    std::fs::create_dir_all(&dir).ok()?;
    // The whole model in the name, not its length: `IBV=1e-3` and `IBV=1e-6` are
    // the same length, and cargo runs these tests in parallel, so keying on the
    // length had one test reading another test's deck.
    let tag: String = model
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect();
    let path = dir.join(format!("bv_{tag}_{v}.sp"));
    let body = deck(model, v).replace(".op\n", "");
    std::fs::write(
        &path,
        format!("{body}.control\nop\nprint i(v1)\n.endc\n.end\n"),
    )
    .ok()?;
    let out = Command::new("ngspice").arg("-b").arg(&path).output().ok()?;
    let text = String::from_utf8_lossy(&out.stdout);
    for line in text.lines() {
        if line.trim_start().starts_with("i(v1)") {
            let (_, rhs) = line.split_once('=')?;
            if let Ok(x) = rhs.split_whitespace().next()?.parse::<f64>() {
                return Some(x);
            }
        }
    }
    None
}

fn agree(what: &str, model: &str, v: f64, rel: f64) {
    let fc = fairchild(model, v);
    let Some(ng) = ngspice(model, v) else {
        eprintln!("ngspice not available — skipping '{what}'");
        return;
    };
    let denom = ng.abs().max(1e-300);
    let err = (fc - ng).abs() / denom;
    assert!(
        err < rel,
        "{what} at V_rev={v}: fairchild {fc:.9e}, ngspice {ng:.9e} (rel {err:.2e} > {rel:.0e})"
    );
}

const ZENER: &str = "IS=1e-14 N=1 BV=5 IBV=1e-3";

/// The knee, crossed: the exponential tail approaching it, the knee current
/// itself, and runaway past it. Each point is a separate statement about the
/// offset of the same exponential, which is why the sweep and not one number.
///
/// Mild reverse is deliberately not here — see
/// `mild_reverse_is_exact_shockley_and_ngspice_is_fitted`.
///
/// 1e-4 rather than tighter because ngspice prints six digits and the runaway
/// points are steep enough that the last one moves.
#[test]
fn the_breakdown_knee_agrees_with_ngspice_across_it() {
    for v in [4.5, 4.9, 5.0, 5.1, 5.2, 5.5] {
        agree("Zener knee", ZENER, v, 1e-4);
    }
}

/// Below the knee the two simulators differ, and fairchild is the exact one.
///
/// ngspice smooths its reverse saturation with a cubic fit,
/// `I = -IS·(1 + (3·vte/(vd·e))³)`, which at -0.5 V reads 1.86e-4 low against the
/// Shockley law it is fitting. fairchild evaluates the law. Matching ngspice here
/// would be *less* correct, so the divergence stands — and it is asserted in both
/// directions so it cannot drift into an accident:
///
/// * fairchild against the closed form, to 1e-12;
/// * the gap to ngspice against ngspice's own fit formula, to 2%.
///
/// The second half is what makes this a test rather than a loosened tolerance. It
/// says the difference is the one term we know about, and nothing else.
#[test]
fn mild_reverse_is_exact_shockley_and_ngspice_is_fitted() {
    let vt = SimOptions::default().sim_context().vt();
    for v in [0.5, 1.0, 3.0] {
        let vd = -v;
        let fc = fairchild(ZENER, v);
        let exact = 1e-14 * ((vd / vt).exp() - 1.0);
        assert!(
            (fc - exact).abs() <= 1e-12 * exact.abs(),
            "V_rev={v}: fairchild {fc:.9e} should be the Shockley law's              {exact:.9e} exactly"
        );

        let Some(ng) = ngspice(ZENER, v) else {
            eprintln!("ngspice not available — closed form checked, gap skipped");
            continue;
        };
        // ngspice's own smoothing term, from its source's shape and confirmed by
        // this arithmetic reproducing its printed current.
        let arg = (3.0 * vt / (vd * std::f64::consts::E)).powi(3);
        let ng_predicted = -1e-14 * (1.0 + arg);
        let rel = (ng - ng_predicted).abs() / ng_predicted.abs();
        assert!(
            rel < 2e-2,
            "V_rev={v}: ngspice reads {ng:.9e} and its cubic fit predicts              {ng_predicted:.9e} (rel {rel:.2e}). If this fails, the difference              from fairchild is no longer the term we think it is, and the              divergence needs re-measuring rather than tolerating."
        );
    }
}

/// `IBV` is the current *at* `-BV`, which only holds if the breakdown voltage is
/// adjusted for it. This is the test that says the adjustment exists: with the
/// exponential placed at `-BV` unshifted, `I(-BV)` would be `-IS` and this fails
/// by eleven orders of magnitude.
#[test]
fn the_current_at_minus_bv_is_ibv() {
    for ibv in [1e-3, 1e-2, 1e-6] {
        let model = format!("IS=1e-14 N=1 BV=5 IBV={ibv:e}");
        let got = fairchild(&model, 5.0).abs();
        let rel = (got - ibv).abs() / ibv;
        assert!(
            rel < 1e-3,
            "IBV={ibv:e}: I(-BV) is {got:.6e}, and IBV means it should be {ibv:e} \
             (rel {rel:.2e}). `-IS` here would mean the knee was placed at -BV \
             without solving for the offset."
        );
        agree("I(-BV) = IBV", &model, 5.0, 1e-4);
    }
}

/// The exponential's slope is `1/(N·vt)`, not `1/vt`. Fitted from two points
/// rather than asserted, so a wrong `N` in the breakdown branch shows up as the
/// wrong slope even where the knee itself still lands correctly.
#[test]
fn the_breakdown_slope_follows_n() {
    let vt = SimOptions::default().sim_context().vt();
    for n in [1.0, 1.5, 2.0] {
        let model = format!("IS=1e-14 N={n} BV=5 IBV=1e-3");
        let (i1, i2) = (fairchild(&model, 4.6).abs(), fairchild(&model, 4.8).abs());
        let fitted = (4.8 - 4.6) / (i2 / i1).ln();
        let rel = (fitted - n * vt).abs() / (n * vt);
        assert!(
            rel < 2e-3,
            "N={n}: the breakdown exponential's fitted vte is {fitted:.6} and \
             N·vt is {:.6} (rel {rel:.2e}). Using vt instead of N·vt would fit \
             {vt:.6} for every N.",
            n * vt
        );
    }
}

/// `IBV` below the leakage the card already has at `-BV` has no offset to solve
/// for. ngspice leaves `BV` unshifted there and reads `-IS` at the knee, so the
/// clamp is a real branch and not a guard against a number nobody writes.
#[test]
fn an_ibv_under_the_leakage_floor_leaves_bv_unshifted() {
    // IS·BV/vt = 1.93e-12 is the threshold, so 1e-12 is under it and 1e-9 over.
    for ibv in [1e-12, 1e-14] {
        let model = format!("IS=1e-14 N=1 BV=5 IBV={ibv:e}");
        agree("clamped IBV, at the knee", &model, 5.0, 1e-4);
        agree("clamped IBV, past the knee", &model, 5.2, 1e-4);
    }
    agree(
        "unclamped IBV just above the floor",
        "IS=1e-14 N=1 BV=5 IBV=1e-9",
        5.0,
        1e-3,
    );
}

/// A card with no `BV` must not acquire a knee. This is the regression direction:
/// every existing diode deck is this case, and a breakdown branch that triggers
/// without `BV` would break all of them at once.
#[test]
fn a_card_without_bv_never_breaks_down() {
    for v in [5.0, 20.0, 100.0] {
        let i = fairchild("IS=1e-14 N=1", v);
        assert!(
            (i.abs() - 1e-14).abs() < 1e-16,
            "no BV means flat reverse saturation, but V_rev={v} gives {i:.6e}"
        );
        agree("no BV", "IS=1e-14 N=1", v, 1e-4);
    }
}

/// Breakdown survives `RS`, which makes the junction voltage a solve rather
/// than the terminal voltage.
#[test]
fn breakdown_survives_rs() {
    for v in [5.1, 5.3] {
        agree(
            "breakdown with RS",
            "IS=1e-14 N=1 BV=5 IBV=1e-3 RS=10",
            v,
            1e-3,
        );
    }
}

/// `AREA=2` in breakdown equals **two diodes in parallel**, and this is a
/// deliberate divergence from ngspice's `area=2`.
///
/// # What was measured
///
/// ngspice's breakdown branch is exactly independent of `area` — 4.8, 5.0, 5.1
/// and 5.3 V all return the identical current at `area=1` and `area=2`, ratio
/// 1.0000, while the forward current doubles correctly. The cause is structural:
/// deriving the knee offset from `IS·AREA` doubles the prefactor and lifts the
/// offset by `vte·ln 2`, and the two cancel to the last bit.
///
/// But ngspice does not agree with itself. Two diodes in parallel give exactly
/// twice the breakdown current of one `area=2` diode there.
///
/// # Why the pair is the anchor
///
/// `area_scales_the_diode_exactly` already states this tree's rule: "AREA=2 *is*
/// two devices, so it has to agree with two devices rather than merely being
/// twice something." An `area=10` Zener silently carrying a tenth of its knee
/// current is the failure this codebase refuses. So the comparison is against the
/// pair — which is also ngspice's own answer for the pair — and not against
/// ngspice's `area=N`.
#[test]
fn area_in_breakdown_equals_devices_in_parallel() {
    let solve = |instances: &str| {
        let src = format!(
            "* area vs parallel\n.options gmin=0\n\
             .model dm D (IS=1e-14 N=1 BV=5 IBV=1e-3)\n\
             V1 a 0 DC 5.1\n{instances}.op\n"
        );
        let net = parse_spice(&src).expect("parse");
        let mut reg = DeviceRegistry::new();
        reg.register_builtin_models(&net.models);
        let opts = SimOptions::from_netlist(&net);
        dc_op_nr_with_registry_opts(&net, &reg, &opts)
            .unwrap_or_else(|e| panic!("solve failed on\n{src}\n{e:?}"))
            .vsrc_current("v1")
            .expect("I(v1)")
    };
    let area2 = solve("D1 0 a dm area=2\n");
    let pair = solve("D1 0 a dm\nD2 0 a dm\n");
    let one = solve("D1 0 a dm\n");
    let rel = (area2 - pair).abs() / pair.abs();
    assert!(
        rel < 1e-9,
        "in breakdown, area=2 gives {area2:.9e} and two in parallel give \
         {pair:.9e} (rel {rel:.2e}). Equal to the single-device {one:.9e} instead \
         would mean AREA cancelled out of the knee — which is what ngspice does, \
         and what disagrees with its own parallel pair."
    );
    let doubling = area2.abs() / one.abs();
    assert!(
        (doubling - 2.0).abs() < 1e-6,
        "and it should be a clean doubling, not merely equal to the pair: got \
         {doubling:.6}"
    );
}

/// A Zener shunt regulator, which is the circuit a `BV` card exists for.
///
/// 12 V through 1 kΩ into a 5 V Zener. Three separate things have to work, and
/// each fails differently:
///
/// * **breakdown modelled at all** — without it the diode blocks and `out` sits
///   at the 12 V supply instead of regulating near 5 V. That is #77 §3's "runs and
///   answers the wrong question" in one number, and it is a 2.4× error.
/// * **the knee in the right place** — `out` is set by where the exponential
///   crosses the load line, so it reads the offset, not just the shape.
/// * **limiting on the way in** — unlike every other test here the node is
///   *solved for*, so Newton walks into the knee from 0 V. The breakdown
///   exponential is as steep as the forward one, and one unlimited step into it
///   overflows. This is the only test in the file that exercises that: dropping
///   the mirrored `pnjlim` leaves the other eight green.
#[test]
fn a_zener_shunt_regulator_regulates() {
    const DECK: &str = "* zener shunt regulator\n\
                        .model dz D (IS=1e-14 N=1 BV=5 IBV=1e-3)\n\
                        V1 in 0 DC 12\nR1 in out 1k\nDz 0 out dz\n.op\n";
    let net = parse_spice(DECK).expect("parse");
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);
    let opts = SimOptions::from_netlist(&net);
    let r = dc_op_nr_with_registry_opts(&net, &reg, &opts).expect(
        "a Zener regulator must converge — an unlimited step into the knee \
                 overflows the exponential",
    );
    let out = r.node_voltage("out").expect("V(out)");

    assert!(
        (4.5..5.5).contains(&out),
        "V(out) is {out:.4} V. A 5 V Zener fed 12 V through 1 k regulates just \
         above its knee; 12 V means breakdown is not modelled and the diode is \
         blocking, which is the whole failure this parameter exists to fix."
    );

    // Anchor: ngspice on the same deck. The regulated voltage is where the
    // exponential meets the load line, so agreeing here is agreeing about the
    // offset and the slope at once.
    let dir = std::env::temp_dir().join("fc_bv_golden");
    std::fs::create_dir_all(&dir).ok();
    let path = dir.join("zener_regulator.sp");
    let body = DECK.replace(".op\n", "");
    if std::fs::write(
        &path,
        format!("{body}.control\nop\nprint v(out)\n.endc\n.end\n"),
    )
    .is_err()
    {
        return;
    }
    let Ok(run) = Command::new("ngspice").arg("-b").arg(&path).output() else {
        eprintln!("ngspice not available — regulation checked, anchor skipped");
        return;
    };
    let text = String::from_utf8_lossy(&run.stdout);
    let ng = text.lines().find_map(|l| {
        let l = l.trim_start();
        if l.starts_with("v(out)") {
            l.split_once('=')?
                .1
                .split_whitespace()
                .next()?
                .parse::<f64>()
                .ok()
        } else {
            None
        }
    });
    let Some(ng) = ng else {
        eprintln!("ngspice printed no v(out) — regulation checked, anchor skipped");
        return;
    };
    let rel = (out - ng).abs() / ng.abs();
    assert!(
        rel < 1e-4,
        "V(out): fairchild {out:.9} V, ngspice {ng:.9} V (rel {rel:.2e})"
    );
}

/// The regulator must still regulate with the trust region loosened.
///
/// # What this caught
///
/// A mirrored `pnjlim` for the breakdown exponential — which looks obviously
/// right, since the knee is as steep as forward conduction. At the default
/// `vmax=0.5` it changed no answer on any deck, including 1 kV into 10 Ω, so
/// every other test in this file stayed green with it and without it.
///
/// At `vmax=1e6` it read `out = 12 V` with 1.2e-11 A through the Zener — the
/// diode blocking, a 2.4× error, silent. `vd_prev` is state the outer Newton
/// cannot see: the mirror compressed the walk into the knee while the free node
/// jumped straight to the supply, so the stamp kept saying "barely conducting",
/// and in that region the terminal current is under `abstol` so the visible
/// unknowns stopped moving and Newton reported success.
///
/// So this is not a convergence test with a wide tolerance. It is a *value* test
/// under a setting that removes the thing accidentally covering for a model
/// fault, which is the category `options_take_effect.rs` exists for: a feature
/// tested only at its default cannot tell you what the default was hiding.
#[test]
fn a_zener_regulator_regulates_with_a_loosened_trust_region() {
    const DECK: &str = "* zener shunt regulator\n\
                        .model dz D (IS=1e-14 N=1 BV=5 IBV=1e-3)\n\
                        V1 in 0 DC 12\nR1 in out 1k\nDz 0 out dz\n.op\n";
    let net = parse_spice(DECK).expect("parse");
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);

    // The default is the control: both settings must give the same answer, since
    // `vmax` bounds the path and not the fixed point.
    let base = SimOptions::from_netlist(&net);
    let tight = dc_op_nr_with_registry_opts(&net, &reg, &base)
        .expect("solve at the default vmax")
        .node_voltage("out")
        .expect("V(out)");

    for vmax in [1e3, 1e6, 1e9] {
        let opts = SimOptions {
            vmax,
            ..base.clone()
        };
        let out = dc_op_nr_with_registry_opts(&net, &reg, &opts)
            .unwrap_or_else(|e| panic!("vmax={vmax:e}: {e:?}"))
            .node_voltage("out")
            .expect("V(out)");
        let rel = (out - tight).abs() / tight;
        assert!(
            rel < 1e-6,
            "vmax={vmax:e}: V(out) is {out:.6} V against {tight:.6} V at the \
             default. `vmax` bounds the Newton path, not the answer, so these must \
             agree. 12 V here means the Zener converged to a blocking state — a \
             limiter compressing the walk into the knee while the node ran ahead."
        );
    }
}
