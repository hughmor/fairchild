//! End-to-end test: `.options solver=klu` directive routes through to
//! the SuiteSparse KLU backend and produces the same DC operating point
//! as the dense / faer-sparse paths.
//!
//! Gated on the `klu` cargo feature.

#![cfg(feature = "klu")]

use fairchild_core::device_registry::DeviceRegistry;
use fairchild_core::newton::dc_op_nr_with_registry_opts;
use fairchild_core::options::SimOptions;
use fairchild_core::solver::SolverKind;
use fairchild_parser::parse_spice;

#[test]
fn klu_solves_resistor_divider() {
    // V1 1V → 1k → out → 1k → 0   ⇒  V(out) = 0.5 V
    let net = parse_spice("* divider\nV1 in 0 DC 1\nR1 in out 1k\nR2 out 0 1k\n.op\n").unwrap();
    let reg = DeviceRegistry::new();
    let mut opts = SimOptions::from_netlist(&net);
    opts.solver = SolverKind::Klu;
    let r = dc_op_nr_with_registry_opts(&net, &reg, &opts).unwrap();
    let v_out = r.node_voltage("out").unwrap();
    assert!(
        (v_out - 0.5).abs() < 1e-9,
        "V(out) = {v_out} (expected 0.5)"
    );
}

#[test]
fn klu_options_directive_selects_klu_backend() {
    // Parser routes `.options solver=klu` into `SimOptions.solver`.
    let net = parse_spice(
        "* divider with directive\nV1 in 0 DC 2\nR1 in out 1k\nR2 out 0 1k\n\
         .options solver=klu\n.op\n",
    )
    .unwrap();
    let opts = SimOptions::from_netlist(&net);
    assert_eq!(opts.solver, SolverKind::Klu);

    let reg = DeviceRegistry::new();
    let r = dc_op_nr_with_registry_opts(&net, &reg, &opts).unwrap();
    let v_out = r.node_voltage("out").unwrap();
    assert!(
        (v_out - 1.0).abs() < 1e-9,
        "V(out) = {v_out} (expected 1.0)"
    );
}

#[test]
fn klu_transient_diode_rc_matches_dense() {
    // Exercises the symbolic/numeric split: one DC OP solve + hundreds
    // of NR iterations across timesteps, all reusing the same KLU
    // symbolic factorisation through `klu_refactor`.  Compare the final
    // V(out) at t=stop against the dense path.
    use fairchild_core::tran::tran_nr_with_registry_var_opts;

    let net = parse_spice(
        "* RC + diode rectifier\n\
         Vin in 0 PULSE(0 1 0 1n 1n 10n 20n)\n\
         R1 in n1 1k\nD1 n1 out myd\nR2 out 0 10k\nC1 out 0 1n\n\
         .model myd D (Is=1e-14 N=1)\n.tran 0.5n 30n\n",
    )
    .unwrap();
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_diodes(&net.models);

    let mut opts_klu = SimOptions::from_netlist(&net);
    opts_klu.solver = SolverKind::Klu;
    let r_klu = tran_nr_with_registry_var_opts(&net, 0.5e-9, 30e-9, &reg, &opts_klu).unwrap();

    let mut opts_dense = SimOptions::from_netlist(&net);
    opts_dense.solver = SolverKind::Dense;
    let r_dense = tran_nr_with_registry_var_opts(&net, 0.5e-9, 30e-9, &reg, &opts_dense).unwrap();

    let v_klu_last = r_klu
        .node_voltages
        .get("out")
        .unwrap()
        .last()
        .copied()
        .unwrap();
    let v_dense_last = r_dense
        .node_voltages
        .get("out")
        .unwrap()
        .last()
        .copied()
        .unwrap();
    assert!(
        (v_klu_last - v_dense_last).abs() < 1e-6,
        "KLU transient diverged from Dense: V(out)_klu={v_klu_last} V(out)_dense={v_dense_last}"
    );
}

#[test]
fn klu_matches_dense_on_nonlinear_diode_circuit() {
    // Forward-biased diode through 1 kΩ — KLU and Dense must agree.
    let net = parse_spice(
        "* diode\nV1 in 0 DC 1.0\nR1 in n1 1k\nD1 n1 0 myd\n\
         .model myd D (Is=1e-14 N=1)\n.op\n",
    )
    .unwrap();
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_diodes(&net.models);

    let mut opts_klu = SimOptions::from_netlist(&net);
    opts_klu.solver = SolverKind::Klu;
    let r_klu = dc_op_nr_with_registry_opts(&net, &reg, &opts_klu).unwrap();

    let mut opts_dense = SimOptions::from_netlist(&net);
    opts_dense.solver = SolverKind::Dense;
    let r_dense = dc_op_nr_with_registry_opts(&net, &reg, &opts_dense).unwrap();

    let v_klu = r_klu.node_voltage("n1").unwrap();
    let v_dense = r_dense.node_voltage("n1").unwrap();
    assert!(
        (v_klu - v_dense).abs() < 1e-9,
        "KLU vs Dense: V(n1) klu={v_klu} dense={v_dense}"
    );
}

// ---------------------------------------------------------------------------
// Cross-backend agreement — the paths the divider tests above never reach.
//
// The KLU factorisation caches a CSC structure and a slot map derived from the
// structural pattern (mirroring `FaerSparseFactorisation`), so the interesting
// cases are the ones that exercise the cache rather than build it once:
//
//   * a nonlinear solve, where `klu_refactor` runs on iterations 2..n;
//   * a device that stamps a cell which is *zero at x = 0*, so the refill walk
//     has to notice the pattern growing and trigger exactly one rebuild;
//   * a transient, where the same handle is reused across every timestep;
//   * `.noise`, which goes through the transpose solve.
//
// Each asserts KLU against faer-sparse on the identical netlist: two spellings
// of one answer, which is the shape that catches a stale or mis-slotted cache.
// ---------------------------------------------------------------------------

use fairchild_core::noise::noise_analysis;
use fairchild_core::tran::tran_nr_with_registry_opts;

fn both_backends<T, F>(net_src: &str, run: F) -> (T, T)
where
    F: Fn(&fairchild_parser::Netlist, &DeviceRegistry, &SimOptions) -> T,
{
    let net = parse_spice(net_src).unwrap();
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);
    let base = SimOptions::from_netlist(&net);

    let mut sparse = base.clone();
    sparse.solver = SolverKind::Sparse;
    let mut klu = base;
    klu.solver = SolverKind::Klu;
    (run(&net, &reg, &sparse), run(&net, &reg, &klu))
}

/// Nonlinear DC: several NR iterations, so most solves take the
/// `klu_refactor` fast path rather than a fresh analyze.
#[test]
fn klu_matches_sparse_on_a_nonlinear_dc_solve() {
    let src = "* diode ladder\n\
               .model dm D (IS=1e-14 N=1)\n\
               V1 in 0 DC 3\n\
               R1 in a 1k\nD1 a b dm\nR2 b c 1k\nD2 c d dm\nR3 d 0 1k\n\
               .op\n";
    let (sp, kl) = both_backends(src, |n, r, o| {
        let res = dc_op_nr_with_registry_opts(n, r, o).expect("DC OP");
        ["a", "b", "c", "d"]
            .iter()
            .map(|k| res.node_voltage(k).unwrap())
            .collect::<Vec<_>>()
    });
    for (i, (s, k)) in sp.iter().zip(kl.iter()).enumerate() {
        assert!((s - k).abs() < 1e-9, "node {i}: sparse {s}, klu {k}");
    }
}

/// A cell that is **zero at the x = 0 initial guess and non-zero at the
/// solution** — the case the cached slot map has to detect and rebuild for.
///
/// `I = V(c)·V(a)` has `∂I/∂V(a) = V(c)`, which is exactly zero on the first
/// iterate and non-zero once `V(c)` moves, so the active set genuinely grows.
/// Verified by instrumenting `KluFactorisation::rebuild`: this deck rebuilds
/// twice (build, then growth) where the linear decks above rebuild once.
///
/// What this pins is that the two backends agree *while* the pattern grows.
/// It is deliberately not claimed as a regression guard for growth detection
/// itself: disabling that detection leaves this test passing, because a
/// Jacobian missing an entry costs Newton iterations rather than accuracy —
/// the residual is computed from the real matrix either way. Catching a
/// dropped rebuild needs an iteration-count or convergence-failure assertion,
/// which is a sharper instrument than this file has.
#[test]
fn klu_matches_sparse_when_the_active_pattern_grows() {
    let src = "* product term: dI/dV(a) is zero at x=0\n\
               Vc c 0 DC 2\n\
               V1 in 0 DC 1\n\
               R1 in a 1k\n\
               B1 a 0 I=V(c)*V(a)*1m\n\
               R2 a 0 10k\n\
               .op\n";
    let (sp, kl) = both_backends(src, |n, r, o| {
        let res = dc_op_nr_with_registry_opts(n, r, o).expect("DC OP");
        res.node_voltage("a").unwrap()
    });
    assert!(sp.abs() > 1e-6, "expected a non-trivial solution, got {sp}");
    assert!((sp - kl).abs() < 1e-9, "sparse {sp}, klu {kl}");
}

/// One factorisation handle reused across every timestep.
#[test]
fn klu_matches_sparse_across_a_transient() {
    let src = "* rc with a diode clamp\n\
               .model dm D (IS=1e-14 N=1)\n\
               V1 in 0 PULSE(0 2 1n 1n 1n 20n 40n)\n\
               R1 in out 1k\nC1 out 0 1n\nD1 out 0 dm\n\
               .tran 1n 60n\n";
    let (sp, kl) = both_backends(src, |n, r, o| {
        let res = tran_nr_with_registry_opts(n, 1e-9, 60e-9, r, o).expect("transient");
        [10e-9, 25e-9, 45e-9]
            .iter()
            .map(|&t| res.voltage_at("out", t).unwrap())
            .collect::<Vec<_>>()
    });
    for (i, (s, k)) in sp.iter().zip(kl.iter()).enumerate() {
        assert!((s - k).abs() < 1e-9, "sample {i}: sparse {s}, klu {k}");
    }
}

/// `.ac` answers come from the one-shot `solve` — since #75 a direct CSC
/// build instead of a dense n×n round-trip. The sweep solves the 2n×2n
/// `[G −B; B G]` embedding, and the deck carries an inductor so branch-current
/// rows sit in the block too: the shape that catches a block placed one row
/// off. Compared against faer-sparse at every frequency, real and imaginary
/// parts separately — two spellings of one factorisation must agree to
/// round-off, not to a tolerance that would hide a misplaced cell.
#[test]
fn klu_matches_sparse_across_an_ac_sweep() {
    use fairchild_core::ac::ac_analysis_opts;
    let src = "* diode-biased rlc: nonlinear op, then ac over the linearisation\n\
               .model dm D (IS=1e-14 N=1)\n\
               V1 in 0 DC 0.7 AC 1\n\
               R1 in mid 1k\n\
               D1 mid 0 dm\n\
               L1 mid out 1m\n\
               C1 out 0 1n\n\
               R2 out 0 10k\n";
    let freqs: Vec<f64> = (0..31).map(|i| 10f64.powf(1.0 + 0.2 * i as f64)).collect();
    let (sp, kl) = both_backends(src, |n, r, o| {
        let res = ac_analysis_opts(n, &freqs, None, r, o).expect("ac");
        res.voltages.get("out").expect("V(out)").clone()
    });
    assert_eq!(sp.len(), freqs.len());
    for (i, ((sre, sim), (kre, kim))) in sp.iter().zip(kl.iter()).enumerate() {
        let scale = (sre * sre + sim * sim).sqrt().max(1e-30);
        assert!(
            (sre - kre).abs() <= 1e-9 * scale && (sim - kim).abs() <= 1e-9 * scale,
            "freq {i}: sparse ({sre:e}, {sim:e}), klu ({kre:e}, {kim:e})"
        );
    }
}

/// `.noise` is the only caller of the transpose solve. KLU now answers it with
/// `klu_tsolve` on the existing factorisation instead of materialising a dense
/// transpose, so this pins that the two agree.
#[test]
fn klu_matches_sparse_through_the_transpose_solve() {
    let src = "* rc thermal noise\nV1 in 0 DC 0\nR1 in out 10k\nC1 out 0 1n\n";
    let freqs = [1e3, 1e4, 1e5];
    let (sp, kl) = both_backends(src, |n, r, o| {
        let res = noise_analysis(n, &freqs, "out", "0", "v1", r, o).expect("noise");
        res.onoise_psd.clone()
    });
    assert_eq!(sp.len(), freqs.len());
    for (i, (s, k)) in sp.iter().zip(kl.iter()).enumerate() {
        assert!(
            (s - k).abs() <= 1e-12 * s.abs().max(1e-30),
            "freq {i}: sparse {s:e}, klu {k:e}"
        );
    }
}
