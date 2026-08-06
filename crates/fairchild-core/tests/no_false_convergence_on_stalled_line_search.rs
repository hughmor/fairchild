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

/// 1 mA forced through a photodetector's shunt. The shunt is stamped correctly,
/// so there is a well-defined answer near 1 mA × 20 kΩ — the solver just cannot
/// walk to it while the line search is stalling.
const STALLS: &str = "* current source into a photodetector shunt\n\
     .optical_port a\n\
     XPD a p n fc_photodetector r_shunt=20k\n\
     I1 0 p DC 1m\n\
     Rn n 0 1\n.op\n.end\n";

/// Same circuit with the shunt duplicated as an ordinary resistor, which the
/// solver handles without stalling. This is the reference value.
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
        (v - 20.0).abs() < 0.5,
        "expected ~20 V for 1 mA through ~20 kΩ, got {v}"
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
                    (v - 20.0).abs() < 0.5,
                    "itl1={itl1}: converged to {v} V, but the operating point is \
                     ~20 V. A step of vmax/16 marching until reltol*|x| catches \
                     up lands near 31.25 V, which is what this test exists to \
                     catch."
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
        assert!((va - 20.0).abs() < 0.5, "converged to the wrong value {va}");
    }
}
