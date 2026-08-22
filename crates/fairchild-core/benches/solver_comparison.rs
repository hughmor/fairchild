/// Benchmark: dense vs. sparse vs. KLU (optional) linear-system backends.
///
/// # What is tested
///
/// DC operating-point solves on synthetic resistor-mesh circuits of increasing
/// size.  A resistor mesh (N×N grid of 1 kΩ resistors with a corner voltage
/// source) produces an MNA matrix whose sparsity is ~4/N — a good proxy for
/// the kind of nets seen in WDM photonic bus simulations.
///
/// # Why DC OP?
///
/// Each Newton-Raphson iteration is dominated by a single LU factorisation +
/// back-substitution of the MNA matrix.  DC OP is the purest measure of this
/// cost; transient adds integration-companion stamping overhead that is small
/// relative to the solver.
///
/// # Running
///
/// ```bash
/// # Default (dense + sparse):
/// cargo bench -p fairchild-core
///
/// # With KLU (requires `brew install suite-sparse`):
/// cargo bench -p fairchild-core --features klu
///
/// # Quick smoke-check (no warm-up, 1 sample):
/// cargo bench -p fairchild-core -- --quick
/// ```
///
/// HTML report: `target/criterion/report/index.html`
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use fairchild_core::{dc_op_nr_with_registry_opts, DeviceRegistry, SimOptions, SolverKind};
use fairchild_parser::parse_spice;

// ---------------------------------------------------------------------------
// Circuit generators
// ---------------------------------------------------------------------------

/// Build an N×N resistor mesh netlist.
///
/// Topology: corner node `v00` driven to 1 V by `V1`; diagonally opposite
/// corner `vNN` grounded through `R_load`.  Interior nodes connected by 1 kΩ
/// resistors in both horizontal and vertical directions.
///
/// Node count ≈ N² + 2 (including V1 branch row in MNA).
/// Non-zeros per row ≈ 4 (most interior nodes touch 4 resistors) → very sparse.
fn build_resistor_mesh(n: usize) -> String {
    assert!(n >= 2, "mesh must be at least 2×2");
    let mut s = format!("* {}x{} resistor mesh\n", n, n);

    // Voltage source drives the (0,0) corner to 1 V.
    s.push_str("V1 v_0_0 0 DC 1\n");

    // Horizontal resistors: v_(i,j) → v_(i,j+1)
    for i in 0..n {
        for j in 0..(n - 1) {
            s.push_str(&format!("Rh_{i}_{j} v_{i}_{j} v_{i}_{} 1k\n", j + 1));
        }
    }

    // Vertical resistors: v_(i,j) → v_(i+1,j)
    for i in 0..(n - 1) {
        for j in 0..n {
            let i1 = i + 1;
            s.push_str(&format!("Rv_{i}_{j} v_{i}_{j} v_{i1}_{j} 1k\n"));
        }
    }

    // Load resistor grounds the far corner so the circuit has a unique solution.
    s.push_str(&format!("Rload v_{n1}_{n1} 0 1k\n", n1 = n - 1));
    s.push_str(".op\n");
    s
}

/// Build a diode-ladder netlist with N stages.
///
/// Each stage: a series resistor (10 kΩ) and a shockley diode to ground.
/// Introduces nonlinear Jacobian entries — the DC OP needs multiple NR
/// iterations and exercises the full NR+solver path (not just the linear
/// factorisation).  Sparsity is similar to the resistor mesh.
fn build_diode_ladder(n: usize) -> String {
    let mut s = format!("* {n}-stage diode ladder\n");
    s.push_str(".model d1 D (Is=1e-14 N=1)\n");
    s.push_str(&format!("V1 n0 0 DC {:.1}\n", n as f64 * 0.7 + 1.0));
    for i in 0..n {
        s.push_str(&format!("R{i} n{i} n{} 10k\n", i + 1));
        s.push_str(&format!("D{i} n{} 0 d1\n", i + 1));
    }
    s.push_str(".op\n");
    s
}

// ---------------------------------------------------------------------------
// Benchmark helpers
// ---------------------------------------------------------------------------

fn run_dc_op(netlist_src: &str, kind: SolverKind) {
    let netlist = parse_spice(netlist_src).expect("parse");
    let registry = {
        let mut r = DeviceRegistry::new();
        r.register_builtin_diodes(&netlist.models);
        r
    };
    let opts = SimOptions {
        solver: kind,
        ..Default::default()
    };
    dc_op_nr_with_registry_opts(&netlist, &registry, &opts).expect("dc op");
}

// ---------------------------------------------------------------------------
// Benchmark groups
// ---------------------------------------------------------------------------

/// Dense vs. sparse on the resistor mesh across multiple scales.
fn bench_mesh_solver(c: &mut Criterion) {
    let mut group = c.benchmark_group("dc_op/resistor_mesh");

    // Sizes chosen to span the crossover from dense-wins to sparse-wins.
    // N×N mesh → N² nodes.  At N=5 → 25 nodes (dense ok); N=20 → 400 nodes.
    for &n in &[3usize, 5, 8, 12, 16, 20] {
        let src = build_resistor_mesh(n);
        let n_nodes = n * n;
        group.throughput(Throughput::Elements(n_nodes as u64));

        group.bench_with_input(BenchmarkId::new("dense", n_nodes), &src, |b, s| {
            b.iter(|| run_dc_op(s, SolverKind::Dense));
        });
        group.bench_with_input(BenchmarkId::new("sparse", n_nodes), &src, |b, s| {
            b.iter(|| run_dc_op(s, SolverKind::Sparse));
        });
        #[cfg(feature = "klu")]
        group.bench_with_input(BenchmarkId::new("klu", n_nodes), &src, |b, s| {
            b.iter(|| run_dc_op(s, SolverKind::Klu));
        });
    }
    group.finish();
}

/// Dense vs. sparse on the nonlinear diode ladder (full NR iteration count).
fn bench_diode_ladder_solver(c: &mut Criterion) {
    let mut group = c.benchmark_group("dc_op/diode_ladder");

    for &n in &[10usize, 30, 60, 100, 150] {
        let src = build_diode_ladder(n);
        group.throughput(Throughput::Elements(n as u64));

        group.bench_with_input(BenchmarkId::new("dense", n), &src, |b, s| {
            b.iter(|| run_dc_op(s, SolverKind::Dense));
        });
        group.bench_with_input(BenchmarkId::new("sparse", n), &src, |b, s| {
            b.iter(|| run_dc_op(s, SolverKind::Sparse));
        });
        #[cfg(feature = "klu")]
        group.bench_with_input(BenchmarkId::new("klu", n), &src, |b, s| {
            b.iter(|| run_dc_op(s, SolverKind::Klu));
        });
    }
    group.finish();
}

criterion_group!(benches, bench_mesh_solver, bench_diode_ladder_solver);
criterion_main!(benches);
