//! Pseudo-transient continuation: the fourth homotopy stage.
//!
//! Hangs a fictitious capacitor on every KCL row and a fictitious inductor in
//! series with every voltage source, and integrates to steady state. It is
//! reached only when direct Newton, source stepping and gmin stepping have all
//! failed, so the tests here have two jobs and they pull in opposite
//! directions:
//!
//! 1. It must **reach** a circuit the other three cannot.
//! 2. It must **not touch** a circuit they can, and must not return a point
//!    that solves its own fictitious problem rather than the user's.
//!
//! (2) is the one that could pass by accident, so `iters` records which stage
//! answered and the tests assert on it, rather than asserting that some stage
//! did.
//!
//! ## What these do not cover
//!
//! The final plain-Newton solve at the end of `pseudo_transient` — the one that
//! makes the returned point solve the user's equations rather than the
//! fictitious ones — is **not** observable from here, and saying so is better
//! than counting it as covered. Deleting it and returning the trajectory's last
//! point leaves every test in this file passing, because a ramp that runs to
//! completion has already driven the real residual to ~1e-12; the guard only
//! shows on a trajectory that stops early, which needs a step budget these
//! decks do not reach. It is kept as the guarantee it is, not as tested code.
//!
//! The stamp's own sign is pinned in `mna.rs`, by contraction toward a
//! divider's arithmetic solution. It cannot be pinned by a fixed-point
//! argument: the pseudo term is `s·(x − x_ref)` and vanishes at `x_ref`
//! whichever sign it carries, so a test that seeds at the solution and gets it
//! back passes with the sign inverted — the first version of this file did
//! exactly that and could not fail. What the sign decides is whether the
//! iteration *contracts*, and `ptc_solves_a_latch_every_other_stage_fails` is
//! the end-to-end half of that: with the branch rows inverted the latch's
//! residual grew from 3.1 to 2.9e3 over 600 steps and this file failed.
//!
//! Verified by sabotage, three for three: inverting the branch sign, dropping
//! the KCL capacitors, and inverting the sign end to end each make this file
//! fail.

use fairchild_core::mna::CircuitTopology;
use fairchild_core::options::SimOptions;
use fairchild_core::{dc_op_nr_with_registry_opts, DeviceRegistry, SimError};
use fairchild_parser::parse_spice;

/// A latch whose feedback loop is algebraic: two steep behavioural inverters,
/// cross-coupled, with no capacitance anywhere in the loop itself.
///
/// Every earlier stage fails on this, and so does a real transient — the
/// capacitors the deck does carry sit behind 1 kΩ from the sources' outputs, so
/// they damp nothing in the loop. Only a capacitor on *every* row does.
const LATCH: &str = "\
* cross-coupled steep inverters
Vdd vdd 0 1.8
B1 q  0 V = 0.9 - 0.9*tanh(1e4*(V(qb)-0.9))
B2 qb 0 V = 0.9 - 0.9*tanh(1e4*(V(q)-0.9))
Rq  q  qi 1k
Rqb qb qbi 1k
Cq  qi 0 1f
Cqb qbi 0 1f
.op
.end
";

fn solve(deck: &str) -> Result<(CircuitTopology, Vec<f64>, usize), SimError> {
    let net = parse_spice(deck).unwrap();
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);
    let opts = SimOptions::from_netlist(&net);
    let r = dc_op_nr_with_registry_opts(&net, &reg, &opts)?;
    Ok((r.topo, r.x, r.iters))
}

fn node(topo: &CircuitTopology, x: &[f64], name: &str) -> f64 {
    x[topo.node_index[name]]
}

/// The subject: a circuit no earlier stage reaches.
///
/// `iters` records which stage answered — 1 direct, 2 source stepping, 3 gmin
/// stepping, 4 pseudo-transient — so this asserts *which* stage solved it and
/// not merely that something did. Delete the fourth rung and it fails.
#[test]
fn ptc_solves_a_latch_every_other_stage_fails() {
    let (topo, x, stage) = solve(LATCH).expect("PTC must solve the latch");
    assert_eq!(stage, 4, "the latch must be answered by PTC, not earlier");

    // And the answer must be a real solution, not wherever the fictitious
    // trajectory happened to stop. Each behavioural source states its output
    // exactly; check the deck's own equations at the returned point.
    let q = node(&topo, &x, "q");
    let qb = node(&topo, &x, "qb");
    let expect_q = 0.9 - 0.9 * (1e4 * (qb - 0.9)).tanh();
    let expect_qb = 0.9 - 0.9 * (1e4 * (q - 0.9)).tanh();
    assert!(
        (q - expect_q).abs() < 1e-6 && (qb - expect_qb).abs() < 1e-6,
        "PTC returned q={q}, qb={qb}, which do not satisfy the deck: \
         expected q={expect_q}, qb={expect_qb}"
    );
}

/// A symmetric seed reaches the symmetric equilibrium, and a `.nodeset` that
/// breaks the symmetry reaches a latched state.
///
/// Both are correct answers — the circuit genuinely has three — and this pins
/// that PTC starts from the seed it was given rather than from wherever the
/// previous stage left off. Without the nodeset reaching PTC, the second case
/// would return the metastable point too.
#[test]
fn ptc_starts_from_the_nodeset_it_was_given() {
    let (topo, x, _) = solve(LATCH).unwrap();
    let (q, qb) = (node(&topo, &x, "q"), node(&topo, &x, "qb"));
    assert!(
        (q - 0.9).abs() < 1e-3 && (qb - 0.9).abs() < 1e-3,
        "a symmetric seed must reach the symmetric equilibrium, got q={q} qb={qb}"
    );

    let seeded = LATCH.replace(".op", ".nodeset V(qi)=1.5 V(qbi)=0.2\n.op");
    let (topo, x, stage) = solve(&seeded).expect("PTC must solve the seeded latch");
    assert_eq!(stage, 4);
    let (q, qb) = (node(&topo, &x, "q"), node(&topo, &x, "qb"));
    assert!(
        (q - 1.8).abs() < 1e-3 && qb.abs() < 1e-3,
        "a nodeset toward q high must reach the latched state, got q={q} qb={qb}"
    );
}

/// A circuit that converges must never reach the fourth stage.
///
/// This is the guard on "adding PTC cannot move an answer that already
/// existed": if the ladder ever fell through to it on an easy circuit, the
/// answer would come from a different code path and the goldens would be
/// asserting something else. `iters` makes that observable rather than assumed.
///
/// A table with a completeness gate would be better here, but the population is
/// the whole `benchmarks/` and `examples/` tree, which the golden runs already
/// cover; these are the shapes with more than one way to fail.
#[test]
fn easy_circuits_never_reach_the_fourth_stage() {
    let decks: &[(&str, &str)] = &[
        (
            "resistive divider",
            "* divider\nV1 in 0 1.0\nR1 in mid 1k\nR2 mid 0 1k\n.op\n.end\n",
        ),
        (
            "diode rectifier",
            "* diode\n.model d1 D (is=1e-16 n=1)\nV1 a 0 1.0\nD1 a b d1\nR1 b 0 1k\n.op\n.end\n",
        ),
        (
            "cmos latch, no ic",
            "* latch\n.model nm NMOS (vto=0.5 kp=200u)\n.model pm PMOS (vto=-0.5 kp=100u)\n\
             Vdd vdd 0 1.8\nM1 q qb vdd vdd pm w=2u l=0.18u\nM2 q qb 0 0 nm w=1u l=0.18u\n\
             M3 qb q vdd vdd pm w=2u l=0.18u\nM4 qb q 0 0 nm w=1u l=0.18u\n.op\n.end\n",
        ),
        (
            "bjt amplifier",
            "* bjt\n.model npn NPN (is=1e-16 bf=100 vaf=100)\nVcc vcc 0 5\n\
             Rc vcc c 2k\nRb vcc b 200k\nQ1 c b 0 npn\n.op\n.end\n",
        ),
    ];
    for (name, deck) in decks {
        let (_, _, stage) = solve(deck).unwrap_or_else(|e| panic!("{name} failed to solve: {e}"));
        assert!(
            stage < 4,
            "{name} was answered by stage {stage}; it must not reach pseudo-transient \
             continuation, or PTC has started changing answers that already existed"
        );
    }
}

/// A topology with no solution must still be reported as one, not spent on a
/// fourth stage and then reported as non-convergence.
///
/// gmin stepping recognises a singular matrix and the ladder returns there, so
/// PTC never sees it. Getting this wrong turns a precise diagnosis — "two
/// voltage sources in parallel with different values" — into "did not converge
/// after 150 iterations", which sends the user hunting for a bias point that
/// cannot exist.
#[test]
fn a_singular_topology_keeps_its_diagnosis() {
    let deck = "* two sources, one node\nV1 a 0 1.0\nV2 a 0 2.0\nR1 a 0 1k\n.op\n.end\n";
    match solve(deck) {
        Err(SimError::SingularMatrix) => {}
        Err(e) => panic!("expected SingularMatrix, got {e}"),
        Ok(_) => panic!("a circuit with no solution must not converge"),
    }
}
