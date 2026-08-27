//! A stalled Armijo line search must not be mistaken for convergence.
//!
//! When no α satisfies the Armijo condition the search falls through at
//! `ALPHA_MIN = 1/16`, so the accepted step becomes a *constant* `vmax/16`,
//! independent of how far the Newton direction wanted to go. The iterate then
//! marches at fixed velocity, and the relative convergence test
//! `|Δx| < abstol + reltol·|x|` is satisfied the moment `|x|` grows past
//! `(vmax/16)/reltol ≈ 31 V` — which is the tolerance catching up with a walk,
//! not the residual going to zero.
//!
//! Before the guard, this circuit reported success at 31.25 V where the answer
//! is 19.985 V: a +56% error, at *every* iteration limit, because the stopping
//! point was set by `reltol` rather than by the circuit. Refusing to converge on
//! a fallback step turns that into a failure the caller can see.
//!
//! Since #90 the trust region is `vmax + 1e-3·|v|` per row rather than a flat
//! `vmax`, and this circuit now *reaches* its operating point instead of failing
//! to. That is a better outcome and not a weaker test: the rule under test is
//! "never a wrong number", so the assertions moved from "it must fail" to "if it
//! answers, the answer is the closed form".

use fairchild_core::error::SimError;
use fairchild_core::{dc_op_nr_with_registry_opts, options::SimOptions, DeviceRegistry};
use fairchild_parser::parse_spice;

/// 1 mA forced through a photodetector's shunt, in parallel with `gmin`.
///
/// `r_shunt=1MEG` rather than the original 20k, and that is load-bearing. Until
/// `0c24510` the shunt was stamped into the Jacobian but cancelled out of the
/// residual, so `p` carried no conductance at all and sat at `I/gmin`; 20 kΩ
/// reproduced the stall only because it was inert. Now that it conducts, 20 kΩ
/// converges in a handful of iterations and this file would test nothing.
///
/// **Not 1 TΩ**, which this used while the circuit was unreachable. `gmin` is
/// 1e-12 S, so `r_shunt = 1 TΩ` is *exactly* `1/gmin`: the solver's own nodal
/// leakage becomes an equal second path and takes half the current, and the
/// "right answer" turns into a statement about gmin rather than about the walk.
/// 1 MΩ puts the operating point 1000 V from the seed — two thousand `vmax`
/// steps, far beyond any of the iteration limits below, so the walk is still
/// what is being tested — while the leakage drops to 1 ppm of the shunt. The
/// discrimination is unaffected: before #90 this circuit read 501 V.
const STALLS: &str = "* current source into a photodetector shunt\n\
     .optical_port a\n\
     XPD a p n fc_photodetector r_shunt=1MEG\n\
     I1 0 p DC 1m\n\
     Rn n 0 1\n.op\n";

/// `1 mA × 1 MΩ`. Ohm's law and nothing else — see `STALLS` for why the shunt
/// is 1 MΩ and not `1/gmin`.
const STALLS_ANSWER: f64 = 1e3;
/// `V(n)`: the same 1 mA through the 1 Ω sense resistor. A second, independent
/// statement — if only `STALLS_ANSWER` were checked, a current divided the wrong
/// way between the shunt and the solver's leakage could still satisfy it.
const STALLS_VN: f64 = 1e-3;

/// Same circuit with an ordinary 20 kΩ resistor added in parallel, which the
/// solver handles without stalling. Since `0c24510` made the shunt conduct as
/// well as stamp, the two are genuinely in parallel: 10 kΩ, not 20.
const REFERENCE: &str = "* same, plus an explicit parallel resistor\n\
     .optical_port a\n\
     XPD a p n fc_photodetector r_shunt=20k\n\
     Rx p n 20k\n\
     I1 0 p DC 1m\n\
     Rn n 0 1\n.op\n";

fn solve(deck: &str, itl1: usize) -> Result<f64, SimError> {
    solve_both(deck, itl1).map(|(p, _)| p)
}

/// `(V(p), V(n))`, so the current split can be checked as well as its size.
fn solve_both(deck: &str, itl1: usize) -> Result<(f64, f64), SimError> {
    let net = parse_spice(deck)?;
    let registry = DeviceRegistry::new();
    let mut opts = SimOptions::from_netlist(&net);
    opts.itl1 = itl1;
    let r = dc_op_nr_with_registry_opts(&net, &registry, &opts)?;
    Ok((r.node_voltage("p")?, r.node_voltage("n")?))
}

#[test]
fn the_reference_circuit_gives_the_expected_operating_point() {
    let v = solve(REFERENCE, 150).expect("explicit parallel resistor must solve");
    assert!(
        (v - 10.0).abs() < 0.5,
        "expected ~10 V for 1 mA through 20 kΩ ∥ 20 kΩ, got {v}. Reading ~20 V \
         means the detector's r_shunt is inert again — it must conduct in the \
         residual, not only appear in the Jacobian (0c24510)."
    );
}

#[test]
fn a_stalled_line_search_reports_failure_not_a_wrong_answer() {
    // The point is not that this circuit is unsolvable — it is that if the
    // solver cannot get there, it must say so rather than stopping wherever
    // reltol happens to catch up with a fixed-size step.
    //
    // It *does* get there now. The trust region is `vmax + reltol·|v|` per row
    // rather than a flat `vmax`, so a node whose operating point is 1000 V climbs
    // to it instead of walking 0.5 V at a time; #90's fix made refusing to stop
    // early correct, and that made the walk's length matter. The `Err` arm is
    // kept because this test's subject is the *rule*, not this circuit: whatever
    // comes back, it must not be a wrong number.
    for itl1 in [150usize, 1000, 3000] {
        match solve_both(STALLS, itl1) {
            Err(SimError::NoConvergence { .. }) => {}
            Err(other) => panic!("itl1={itl1}: unexpected error {other:?}"),
            Ok((v, vn)) => {
                assert!(
                    (v - STALLS_ANSWER).abs() < 1e-3 * STALLS_ANSWER,
                    "itl1={itl1}: converged to {v} V, but the operating point is \
                     {STALLS_ANSWER} V. A step of vmax/16 marching until \
                     reltol*|x| catches up lands near 31.25 V, which is what \
                     this test exists to catch."
                );
                // …and all of the current went through the shunt, not some of
                // it through the solver's own leakage.
                assert!(
                    (vn - STALLS_VN).abs() < 1e-3 * STALLS_VN,
                    "itl1={itl1}: V(n) = {vn}, expected {STALLS_VN} — the current \
                     divided differently than the conductances say it should"
                );
            }
        }
    }
}

#[test]
fn the_stopping_point_must_not_depend_on_the_iteration_limit() {
    // The tell-tale of the bug: identical wrong answers at 1000 and 3000
    // iterations, because the walk stops when reltol*|x| overtakes the step
    // rather than when the residual falls. Two limits, and if both converge
    // they must agree.
    let a = solve(STALLS, 1000);
    let b = solve(STALLS, 3000);
    if let (Ok(va), Ok(vb)) = (&a, &b) {
        assert!(
            (va - vb).abs() < 1e-6 * va.abs().max(1.0),
            "converged answers differ with the iteration limit: {va} vs {vb}"
        );
        assert!(
            (va - STALLS_ANSWER).abs() < 1e-3 * STALLS_ANSWER,
            "converged to the wrong value {va}"
        );
    }
}
