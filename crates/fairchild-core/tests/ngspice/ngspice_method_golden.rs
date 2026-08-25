//! `.options method` against ngspice, one integrator at a time.
//!
//! The existing transient goldens all pin Backward Euler, deliberately — the
//! tolerances are derived from BE's truncation error at the step they use. That
//! makes them silent about the other two methods, and #93 lived in that silence:
//! `method=gear` ran Backward Euler on the fixed-step path and every golden
//! still passed, because none of them ever asked for anything else.
//!
//! ngspice honours `.options method=gear|trap` on the same deck, so the fix has
//! an external anchor and not just an internal one. That matters here more than
//! usual: the failure mode was two *fairchild* code paths agreeing with each
//! other about the wrong method, which no amount of internal cross-checking can
//! see.
//!
//! # What each method is being held to
//!
//! Not "matches ngspice to 1e-9" — the two simulators do not have to take the
//! same step sequence or make the same first-step choice. What is pinned is the
//! thing a wrong method gets wrong by an order of magnitude:
//!
//! * each method lands on the analytic RC step response within *its own* error
//!   bound at this step size, and
//! * the three methods do not agree with each other, because if they did one of
//!   them would not be running.
//!
//! The analytic answer is the absolute anchor; ngspice is the second opinion.

use std::collections::HashMap;
use std::io::Write;
use std::process::Command;

use fairchild_core::options::SimOptions;
use fairchild_core::{tran_nr_configured, DeviceRegistry, IntegratorMode};
use fairchild_parser::parse_spice;

use super::ngspice_golden::find_ngspice;

/// τ = 1 kΩ · 1 nF = 1 µs, stepped at 100 ns — ten points per time constant, so
/// the difference between a first- and second-order method is percent-level and
/// cannot be mistaken for round-off.
const RC: &str = "* rc step\n\
                  V1 in 0 PULSE(0 1 0 1p 1p 1 2)\n\
                  R1 in out 1k\n\
                  C1 out 0 1n\n\
                  .tran 100n 3u\n";

const TAU: f64 = 1e-6;
const PROBE_TIMES: &[f64] = &[5e-7, 1e-6, 2e-6, 3e-6];

/// The exact single-pole step response, which is what both simulators are
/// approximating.
fn exact(t: f64) -> f64 {
    1.0 - (-t / TAU).exp()
}

fn fairchild_at(method: IntegratorMode) -> Vec<f64> {
    let netlist = parse_spice(RC).expect("parse");
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&netlist.models);
    let opts = SimOptions {
        method,
        ..SimOptions::from_netlist(&netlist)
    };
    let r = tran_nr_configured(&netlist, 100e-9, 3e-6, &registry, &opts)
        .unwrap_or_else(|e| panic!("{method:?} failed: {e:?}"));
    PROBE_TIMES
        .iter()
        .map(|&t| r.voltage_at("out", t).expect("V(out)"))
        .collect()
}

/// ngspice's `V(out)` at the same times, under the same method.
///
/// `None` when ngspice is not installed — the same skip every other golden here
/// takes, and CI installs it so the skip cannot go quiet.
fn ngspice_at(method: &str) -> Option<Vec<f64>> {
    let ng = find_ngspice()?;
    // `.meas tran … FIND … AT=` rather than vector indexing: the index form
    // depends on ngspice's own step sequence, which under an adaptive default is
    // not the one this deck asks for, and it silently returned 0.
    let mut deck =
        String::from("* rc step\nV1 in 0 PULSE(0 1 0 1p 1p 1 2)\nR1 in out 1k\nC1 out 0 1n\n");
    deck.push_str(&format!(".options method={method}\n"));
    deck.push_str(".tran 100n 3u\n");
    for (i, t) in PROBE_TIMES.iter().enumerate() {
        deck.push_str(&format!(".meas tran v_{i} FIND V(out) AT={t:e}\n"));
    }
    deck.push_str(".end\n");

    let mut tmp = tempfile::NamedTempFile::new().ok()?;
    tmp.write_all(deck.as_bytes()).ok()?;
    let out = Command::new(ng).arg("-b").arg(tmp.path()).output().ok()?;
    let text = String::from_utf8_lossy(&out.stdout);
    let mut map = HashMap::new();
    for line in text.lines() {
        if let Some((lhs, rhs)) = line.trim().split_once('=') {
            let key = lhs.trim().to_lowercase();
            if key.starts_with("v_") {
                // `.meas` prints `v_0 = 3.93e-01 at= 5.00e-07`; take the first
                // number only.
                let first = rhs.split_whitespace().next()?;
                if let Ok(v) = first.parse::<f64>() {
                    map.insert(key, v);
                }
            }
        }
    }
    // A parse failure must not read as "ngspice is not installed". `find_ngspice`
    // already answered that question; anything missing from here on is ngspice
    // running and not saying what we asked, which is a broken test rather than an
    // absent dependency.
    let vals: Vec<f64> = (0..PROBE_TIMES.len())
        .map(|i| {
            *map.get(&format!("v_{i}")).unwrap_or_else(|| {
                panic!(
                    "ngspice ran but did not report v_{i} for method={method}. \
                     Deck:\n{deck}\nOutput:\n{text}"
                )
            })
        })
        .collect();
    Some(vals)
}

/// Each method lands within its own error bound on the analytic answer.
///
/// The bounds are per-method on purpose. A single loose bound that all three
/// pass is the shape of test that let #93 through: `gear` running Backward Euler
/// would satisfy any tolerance wide enough for Backward Euler.
#[test]
fn each_method_meets_its_own_error_bound_on_the_analytic_answer() {
    // Worst |error| at h/τ = 0.1. BE is first order and lands near 1.8e-2; TR
    // and BDF-2 are second order and land an order below that. The ceilings are
    // set just above what each method achieves, so a method that quietly became
    // a different one fails.
    for (method, ceiling, floor) in [
        (IntegratorMode::BackwardEuler, 2.5e-2, 5e-3),
        (IntegratorMode::Trapezoidal, 8e-3, 0.0),
        (IntegratorMode::Gear, 1.0e-2, 0.0),
    ] {
        let got = fairchild_at(method);
        let worst = PROBE_TIMES
            .iter()
            .zip(&got)
            .map(|(&t, &v)| (v - exact(t)).abs())
            .fold(0.0f64, f64::max);
        assert!(
            worst < ceiling,
            "{method:?}: worst error {worst:.3e} exceeds {ceiling:.1e}"
        );
        // …and BE must *not* be as good as the second-order methods, or the
        // ceiling above is not distinguishing anything.
        assert!(
            worst > floor,
            "{method:?}: worst error {worst:.3e} is below {floor:.1e}, which \
             means it is not the method it says it is"
        );
    }
}

/// The three methods disagree with each other.
///
/// This is the assertion #93 fails. `be` and `gear` were byte-identical.
#[test]
fn the_three_methods_do_not_agree_with_each_other() {
    let be = fairchild_at(IntegratorMode::BackwardEuler);
    let tr = fairchild_at(IntegratorMode::Trapezoidal);
    let gear = fairchild_at(IntegratorMode::Gear);
    assert_ne!(
        be, tr,
        "be and tr agree exactly — one of them is not running"
    );
    assert_ne!(
        be, gear,
        "be and gear agree exactly — this is #93: BDF-2 demoted to Backward \
         Euler because the fixed-step path never supplied the previous step size"
    );
    assert_ne!(tr, gear, "tr and gear agree exactly");
}

/// And ngspice, on the same deck under the same `.options method`, is close.
///
/// Loose on purpose — 2% of full scale. The two simulators need not take the
/// same first step or make the same startup choice, and pinning them tighter
/// would be pinning coincidence. What it does catch is a method that is
/// wholesale the wrong one, which is a percent-level error at this step size.
#[test]
fn each_method_tracks_ngspice_on_the_same_deck() {
    for (method, ng_name) in [
        (IntegratorMode::BackwardEuler, "gear"), // ngspice's `gear` at order 1
        (IntegratorMode::Trapezoidal, "trap"),
    ] {
        let Some(want) = ngspice_at(ng_name) else {
            eprintln!("skipping: ngspice not found");
            return;
        };
        let got = fairchild_at(method);
        for ((&t, &g), &w) in PROBE_TIMES.iter().zip(&got).zip(&want) {
            assert!(
                (g - w).abs() < 2e-2,
                "{method:?} vs ngspice {ng_name} at t={t:.2e}: {g:.6} vs {w:.6}"
            );
        }
    }
}
