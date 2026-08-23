/// Benchmark: matrix assembly rate, isolated from the linear solve.
///
/// # Why this exists separately from `solver_comparison`
///
/// A stamp happens per matrix entry, per device, per Newton iteration, per
/// timestep — hundreds of millions of them over a transient. Any change to how
/// `MnaMatrix` finds the cell a stamper is writing shows up here and nowhere
/// else: measured end-to-end it disappears into factorisation noise, and then
/// nobody can tell whether the change was worth making.
///
/// So this benchmark does exactly what one Newton iteration does to the matrix
/// and then stops: clear, stamp the netlist elements, stamp every device's
/// Jacobian, floor the diagonal with gmin. No solve.
///
/// # Circuits
///
/// Ring oscillators (`benchmarks/circuits/ring_osc_*.sp`) are the natural
/// input — a MOSFET stamps a 4×4 block, so the assembly is heavy and the solve
/// on a chain topology is cheap, which is the ratio that makes assembly worth
/// optimising. The resistor mesh is the opposite extreme (2 terminals, 4 cells)
/// and is here as the floor: whatever a stamp costs, the mesh pays it with the
/// least device work hiding it.
///
/// # Running
///
/// ```bash
/// cargo bench -p fairchild-core --bench stamp_rate
/// cargo bench -p fairchild-core --bench stamp_rate -- --quick   # smoke check
/// ```
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use fairchild_core::device::{Device, EvalFlags, SimContext};
use fairchild_core::mna::{
    stamp_netlist_scaled_in_place, CircuitTopology, InductorDc, MnaMatrix, RowFloor, StampPlan,
};
use fairchild_core::newton::build_devices_with_footprints;
use fairchild_core::{DeviceRegistry, SimOptions};
use fairchild_parser::{parse_spice, Netlist};
use indexmap::IndexMap;

/// Everything one Newton iteration needs in order to stamp, built once.
struct Bed {
    netlist: Netlist,
    topo: CircuitTopology,
    devices: Vec<Box<dyn Device>>,
    plan: StampPlan,
    ctx: SimContext,
    opts: SimOptions,
    mat: MnaMatrix,
    x: Vec<f64>,
}

impl Bed {
    fn new(src: &str) -> Bed {
        let netlist = parse_spice(src).expect("parse");
        let registry = {
            let mut r = DeviceRegistry::new();
            r.register_builtin_models(&netlist.models);
            r
        };
        let opts = SimOptions::from_netlist(&netlist);
        let ctx = opts.sim_context();
        let mut topo = CircuitTopology::build_resolved(&netlist, &ctx, &registry);
        let (devices, foot) =
            build_devices_with_footprints(&netlist, &mut topo, &ctx, &registry).expect("devices");
        let plan = StampPlan::new(&topo, &netlist, &foot);
        let mut devices = devices;
        plan.resolve_device_cells(&mut devices);
        let mat = MnaMatrix::with_pattern(topo.size, plan.pattern.clone());
        // A nonzero operating point, so device Jacobians are not all evaluated
        // at the trivial x = 0 where a stamp can be skipped as zero.
        let x = (0..topo.size)
            .map(|i| 0.1 + 0.01 * (i % 7) as f64)
            .collect();
        Bed {
            netlist,
            topo,
            devices,
            plan,
            ctx,
            opts,
            mat,
            x,
        }
    }

    /// One full assembly pass: exactly what `residual_l2` / `nr_inner` do to
    /// the matrix before handing it to the solver.
    fn assemble(&mut self) {
        for dev in self.devices.iter_mut() {
            dev.eval(&self.x, EvalFlags::dc(), &self.ctx);
        }
        self.stamp();
    }

    /// The same pass with device *evaluation* hoisted out — only the writes
    /// into the matrix remain. `eval` is transcendental model math and has
    /// nothing to do with how a cell is addressed, so a change to the
    /// addressing shows up as a fraction of `assemble` and as the whole of
    /// this. The gap between the two is the ceiling on any such change.
    fn stamp(&mut self) {
        let empty: IndexMap<String, (f64, f64)> = IndexMap::new();
        stamp_netlist_scaled_in_place(
            &mut self.mat,
            &self.topo,
            &self.netlist,
            1.0,
            &empty,
            &empty,
            Some(&self.plan),
            InductorDc::Short,
        );
        for dev in self.devices.iter_mut() {
            dev.load_residual(&mut self.mat.b);
            dev.load_jacobian(&mut self.mat);
        }
        self.topo
            .stamp_gmin(&mut self.mat.a, self.opts.gmin, RowFloor::PinEmptyRows);
    }
}

/// N×N resistor mesh — the same generator `solver_comparison` uses, kept here
/// rather than shared because criterion compiles each bench as its own binary.
fn resistor_mesh(n: usize) -> String {
    let mut s = format!("* {n}x{n} resistor mesh\nV1 v_0_0 0 DC 1\n");
    for i in 0..n {
        for j in 0..(n - 1) {
            s.push_str(&format!("Rh_{i}_{j} v_{i}_{j} v_{i}_{} 1k\n", j + 1));
        }
    }
    for i in 0..(n - 1) {
        for j in 0..n {
            s.push_str(&format!("Rv_{i}_{j} v_{i}_{j} v_{}_{j} 1k\n", i + 1));
        }
    }
    s.push_str(&format!("Rload v_{}_{} 0 1k\n.op\n.end\n", n - 1, n - 1));
    s
}

fn bench_ring_osc(c: &mut Criterion) {
    let mut g = c.benchmark_group("stamp/ring_osc");
    for stages in [21usize, 101, 499] {
        let path = format!(
            "{}/../../benchmarks/circuits/ring_osc_{stages}.sp",
            env!("CARGO_MANIFEST_DIR")
        );
        let Ok(src) = std::fs::read_to_string(&path) else {
            eprintln!("skipping ring_osc_{stages}: {path} not found");
            continue;
        };
        let mut bed = Bed::new(&src);
        // One element per stamped cell, so criterion reports cells/second —
        // the number that is actually comparable across circuit sizes.
        g.throughput(Throughput::Elements(bed.plan.pattern.nnz as u64));
        g.bench_with_input(BenchmarkId::new("assemble", stages), &stages, |b, _| {
            b.iter(|| {
                bed.assemble();
                black_box(&bed.mat.b[0]);
            })
        });
        g.bench_with_input(BenchmarkId::new("stamp", stages), &stages, |b, _| {
            b.iter(|| {
                bed.stamp();
                black_box(&bed.mat.b[0]);
            })
        });
    }
    g.finish();
}

fn bench_mesh(c: &mut Criterion) {
    let mut g = c.benchmark_group("stamp/resistor_mesh");
    for n in [20usize, 40, 60] {
        let mut bed = Bed::new(&resistor_mesh(n));
        g.throughput(Throughput::Elements(bed.plan.pattern.nnz as u64));
        g.bench_with_input(BenchmarkId::new("stamp", n), &n, |b, _| {
            b.iter(|| {
                bed.stamp();
                black_box(&bed.mat.b[0]);
            })
        });
    }
    g.finish();
}

criterion_group!(benches, bench_ring_osc, bench_mesh);
criterion_main!(benches);
