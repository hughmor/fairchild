//! A circuit with no unique DC solution must say so, not blame Newton.
//!
//! The homotopy tries direct NR, then source-stepping, then gmin-stepping. Each
//! stage detected `SingularMatrix` and logged it under `verbose`, but the final
//! `Err` was hard-coded to `NoConvergence { iters }` — so an impossible topology
//! reported "Newton-Raphson did not converge after 150 iterations" and sent the
//! reader looking for a bias point that cannot exist. Continuation cannot rescue
//! a circuit with no solution, so the real diagnosis has to survive the stages.

use fairchild_core::error::SimError;
use fairchild_core::{dc_op_nr_with_registry_opts, options::SimOptions, DeviceRegistry};
use fairchild_parser::parse_spice;

fn op(deck: &str) -> Result<(), SimError> {
    let net = parse_spice(deck)?;
    let registry = DeviceRegistry::new();
    let opts = SimOptions::from_netlist(&net);
    dc_op_nr_with_registry_opts(&net, &registry, &opts).map(|_| ())
}

#[test]
fn parallel_voltage_sources_that_disagree_report_singular() {
    // V1 and V2 both fix node `a`, to different values: no solution exists.
    let err = op("* conflict\nV1 a 0 DC 1\nV2 a 0 DC 2\nR1 a 0 1k\n.op\n")
        .expect_err("two sources fixing one node to different values has no solution");
    assert!(
        matches!(err, SimError::SingularMatrix),
        "expected SingularMatrix, got {err:?} — a structurally impossible circuit \
         must not be reported as a convergence failure"
    );
}

#[test]
fn a_loop_of_voltage_sources_reports_singular() {
    // V1 + V2 fix V(b) at 2, V3 fixes it at 3: the loop is over-determined.
    let err = op("* loop\nV1 a 0 DC 1\nV2 b a DC 1\nV3 b 0 DC 3\nR1 a 0 1k\n.op\n")
        .expect_err("an over-determined voltage-source loop has no solution");
    assert!(matches!(err, SimError::SingularMatrix), "got {err:?}");
}

#[test]
fn a_solvable_circuit_is_unaffected() {
    op("* divider\nV1 a 0 DC 1\nR1 a m 1k\nR2 m 0 1k\n.op\n")
        .expect("a plain divider must still solve");
}

#[test]
fn a_genuinely_hard_nonlinear_circuit_still_reports_nonconvergence() {
    // The distinction has to cut both ways: when the matrix is fine and Newton
    // simply cannot get there, NoConvergence is still the honest answer. One NR
    // iteration is not enough for a diode, so this fails for the right reason.
    let net = parse_spice(
        "* laser into a photodetector: nonlinear, well-posed, solvable\n\
         .optical_port a\n\
         XL1 a fc_cw_laser power_mW=0.1 wavelength_nm=1550\n\
         XPD1 a p 0 fc_photodetector responsivity=0.7 r_shunt=20k\n\
         Rl p 0 2k\n.op\n",
    )
    .unwrap();
    let registry = DeviceRegistry::new();
    let mut opts = SimOptions::from_netlist(&net);
    opts.itl1 = 1;
    opts.srcsteps = 0;
    opts.gmin_max = 0.0;
    match dc_op_nr_with_registry_opts(&net, &registry, &opts) {
        Err(SimError::NoConvergence { .. }) => {}
        Err(other) => {
            panic!("expected NoConvergence for a well-posed but hard circuit, got {other:?}")
        }
        Ok(_) => panic!("expected NoConvergence: one NR iteration cannot solve a diode"),
    }
}
