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

use fairchild_core::error::SimError;
use fairchild_core::{dc_op_nr_with_registry_opts, options::SimOptions, DeviceRegistry};
use fairchild_parser::parse_spice;

/// 1 mA forced through a photodetector's shunt. There is a well-defined answer
/// at 1 mA × `r_shunt` — the solver just cannot walk to it while the line
/// search is stalling.
///
/// `r_shunt=1T` rather than the original 20k, and that is load-bearing. Until
/// `0c24510` the shunt was stamped into the Jacobian but cancelled out of the
/// residual, so `p` carried no conductance at all and sat at `I/gmin`; 20 kΩ
/// reproduced the stall only because it was inert. Now that it conducts, 20 kΩ
/// converges in a handful of iterations and this file would test nothing. 1 TΩ
/// puts the target back at 1e9 V, out of reach of a `vmax/16` walk, which is
/// the condition the guard exists for.
const STALLS: &str = "* current source into a photodetector shunt\n\
     .optical_port a\n\
     XPD a p n fc_photodetector r_shunt=1T\n\
     I1 0 p DC 1m\n\
     Rn n 0 1\n.op\n.end\n";

/// What `STALLS` would settle at if it ever did: 1 mA × 1 TΩ.
const STALLS_ANSWER: f64 = 1e9;

/// Same circuit with an ordinary 20 kΩ resistor added in parallel, which the
/// solver handles without stalling. Since `0c24510` made the shunt conduct as
/// well as stamp, the two are genuinely in parallel: 10 kΩ, not 20.
const REFERENCE: &str = "* same, plus an explicit parallel resistor\n\
     .optical_port a\n\
     XPD a p n fc_photodetector r_shunt=20k\n\
     Rx p n 20k\n\
     I1 0 p DC 1m\n\
     Rn n 0 1\n.op\n.end\n";

fn solve(deck: &str, itl1: usize) -> Result<f64, SimError> {
    let net = parse_spice(deck)?;
    let registry = DeviceRegistry::new();
    let mut opts = SimOptions::from_netlist(&net);
    opts.itl1 = itl1;
    let r = dc_op_nr_with_registry_opts(&net, &registry, &opts)?;
    r.node_voltage("p")
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
    for itl1 in [150usize, 1000, 3000] {
        match solve(STALLS, itl1) {
            Err(SimError::NoConvergence { .. }) => {}
            Err(other) => panic!("itl1={itl1}: unexpected error {other:?}"),
            Ok(v) => {
                // If a future change makes this converge, that is welcome — but
                // only to the right answer.
                assert!(
                    (v - STALLS_ANSWER).abs() < 1e-3 * STALLS_ANSWER,
                    "itl1={itl1}: converged to {v} V, but the operating point is \
                     {STALLS_ANSWER} V. A step of vmax/16 marching until \
                     reltol*|x| catches up lands near 31.25 V, which is what \
                     this test exists to catch."
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
