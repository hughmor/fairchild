//! The answer must not depend on which factorisation the sparse solver picks.
//!
//! `FaerSparseFactorisation` pins faer's choice to **simplicial**
//! (`SUPERNODAL_THRESHOLD`). Before that it left the choice to faer's flop-ratio
//! heuristic, and a device with an internal node flipped it to supernodal — whose
//! dense-block kernel is the wrong shape of work for a circuit matrix. That cost
//! 10× per timestep on a MOSFET ladder and 5× on a BJT one, with the step count
//! unchanged (#99).
//!
//! Performance lives in `benchmarks/`, not here — a timing assertion in
//! `cargo test` is a flake waiting for a loaded CI runner. What belongs here is
//! the **correctness** invariant the change rests on: pinning the factorisation
//! strategy must not move an answer, on exactly the matrix structure that used to
//! flip it.
//!
//! # Why device internal nodes specifically
//!
//! They are what changes the sparsity structure enough to matter: an internal row
//! is allocated after the voltage-source rows, at the end of the matrix, and
//! couples back to its device's terminals near the front. A test built from plain
//! `R`/`C` ladders cannot reach that structure, which is why the pathology went
//! unnoticed — every existing scaling benchmark is a ladder.

use fairchild_core::device_registry::DeviceRegistry;
use fairchild_core::newton::dc_op_nr_with_registry_opts;
use fairchild_core::options::SimOptions;
use fairchild_core::solver::SolverKind;
use fairchild_core::tran_nr_with_registry_opts;
use fairchild_parser::parse_spice;

/// A ladder of devices whose series resistances give each one internal nodes.
///
/// Both families, because they allocate different numbers (`RD`/`RS` is two per
/// MOSFET, `RB`/`RC`/`RE` three per BJT) and the pathology was worse for the
/// *smaller* count — so a test with only one of them would leave the more
/// surprising case uncovered.
fn deck(n: usize, analysis: &str) -> String {
    let mut s = String::from(
        "* internal nodes, both families\n\
         .model nm NMOS (VTO=0.7 KP=200u RD=50 RS=50)\n\
         .model qm NPN (IS=1e-16 BF=100 RB=100 RC=20 RE=5)\n\
         VDD vdd 0 DC 3\n\
         VG g 0 PULSE(0 2 0 1n 1n 1u 2u)\n\
         VB b 0 DC 0.7\n",
    );
    for i in 0..n {
        s.push_str(&format!("RL{i} vdd d{i} 2k\n"));
        s.push_str(&format!("M{i} d{i} g 0 0 nm W=10u L=1u\n"));
        s.push_str(&format!("RC{i} vdd c{i} 2k\n"));
        s.push_str(&format!("Q{i} c{i} b 0 qm\n"));
    }
    s.push_str(analysis);
    s
}

fn solve_dc(src: &str, kind: SolverKind) -> Vec<(String, f64)> {
    let net = parse_spice(src).expect("parse");
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);
    let opts = SimOptions {
        solver: kind,
        ..SimOptions::from_netlist(&net)
    };
    let r = dc_op_nr_with_registry_opts(&net, &reg, &opts)
        .unwrap_or_else(|e| panic!("{kind:?} failed: {e:?}"));
    let mut v: Vec<(String, f64)> = r.all_voltages().map(|(n, x)| (n.to_string(), x)).collect();
    v.sort_by(|a, b| a.0.cmp(&b.0));
    v
}

/// Dense and sparse must agree on a matrix full of device internal nodes.
///
/// Dense LU is the anchor: it has no ordering, no supernodes and no strategy to
/// pick, so it cannot share a fault with the sparse path. That is what makes this
/// an absolute comparison rather than two subsystems agreeing.
#[test]
fn dense_and_sparse_agree_with_internal_nodes_everywhere() {
    let src = deck(24, ".op\n");
    let dense = solve_dc(&src, SolverKind::Dense);
    let sparse = solve_dc(&src, SolverKind::Sparse);
    assert_eq!(
        dense.len(),
        sparse.len(),
        "the two backends disagree about how many unknowns there are"
    );
    let opts = SimOptions::default();
    for ((na, a), (nb, b)) in dense.iter().zip(&sparse) {
        assert_eq!(na, nb, "node order differs between backends");
        // The solver's own convergence bound: below it the difference is where
        // each Newton stopped, not what the circuit is.
        let bound = opts.vntol + opts.reltol * a.abs();
        assert!(
            (a - b).abs() < bound,
            "V({na}): dense {a:.9}, sparse {b:.9} — they differ by more than the \
             convergence bound {bound:.2e}, so the factorisation strategy changed \
             the answer"
        );
    }
}

/// And in transient, which is where the cost showed up.
///
/// The pathology was per-timestep, so the DC comparison above would have passed
/// throughout. This is the shape of run that pays for it.
#[test]
fn dense_and_sparse_agree_in_transient_with_internal_nodes() {
    let src = deck(12, ".tran 20n 400n\n");
    let net = parse_spice(&src).expect("parse");
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);

    let run = |kind: SolverKind| {
        let opts = SimOptions {
            solver: kind,
            ..SimOptions::from_netlist(&net)
        };
        tran_nr_with_registry_opts(&net, 20e-9, 400e-9, &reg, &opts)
            .unwrap_or_else(|e| panic!("{kind:?} transient failed: {e:?}"))
    };
    let (dense, sparse) = (run(SolverKind::Dense), run(SolverKind::Sparse));
    assert_eq!(
        dense.time.len(),
        sparse.time.len(),
        "the two backends took different numbers of steps, so the factorisation \
         is influencing the step controller"
    );
    let opts = SimOptions::default();
    for probe in ["d0", "c0"] {
        for &t in &[100e-9, 250e-9, 390e-9] {
            let a = dense.voltage_at(probe, t).expect("dense sample");
            let b = sparse.voltage_at(probe, t).expect("sparse sample");
            let bound = opts.vntol + opts.reltol * a.abs();
            assert!(
                (a - b).abs() < bound,
                "V({probe}) at {t:e}: dense {a:.9}, sparse {b:.9}, bound {bound:.2e}"
            );
        }
    }
}
