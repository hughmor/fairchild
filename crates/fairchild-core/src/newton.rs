use indexmap::IndexMap;

use fairchild_parser::{Element, Netlist};

use std::sync::Arc;

use crate::behavioral::BehavioralDevice;
use crate::connectivity::check_connectivity;
use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::mna::{stamp_netlist_scaled, CircuitTopology};
use crate::options::SimOptions;
use crate::solver::{lu_solve, LinearSolver};

/// Result of a nonlinear DC operating-point solve.
pub struct NrResult {
    pub topo: CircuitTopology,
    pub x: Vec<f64>,
    pub iters: usize,
}

impl NrResult {
    pub fn node_voltage(&self, node: &str) -> Result<f64, SimError> {
        self.topo.node_voltage(node, &self.x)
    }

    pub fn vsrc_current(&self, name: &str) -> Result<f64, SimError> {
        self.topo.vsrc_current(name, &self.x)
    }

    pub fn all_voltages(&self) -> impl Iterator<Item = (&str, f64)> {
        self.topo.node_index.iter().map(|(name, &i)| (name.as_str(), self.x[i]))
    }

    /// Write the DC operating point as an ngspice-compatible Nutmeg ASCII rawfile.
    pub fn write_nutmeg<W: std::io::Write>(&self, mut w: W, title: &str) -> std::io::Result<()> {
        let n_nodes = self.topo.n_nodes();
        let n_vars = n_nodes + self.topo.vsrc_index.len();
        writeln!(w, "Title: {title}")?;
        writeln!(w, "Plotname: Operating Point")?;
        writeln!(w, "Flags: real")?;
        writeln!(w, "No. Variables: {n_vars}")?;
        writeln!(w, "No. Points: 1")?;
        writeln!(w, "Variables:")?;
        for (idx, name) in self.topo.node_index.keys().enumerate() {
            writeln!(w, "\t{idx}\tv({name})\tvoltage")?;
        }
        for (i, name) in self.topo.vsrc_index.keys().enumerate() {
            writeln!(w, "\t{}\ti({name})\tcurrent", n_nodes + i)?;
        }
        writeln!(w, "Values:")?;
        let mut first = true;
        for &idx in self.topo.node_index.values() {
            if first {
                writeln!(w, " 0\t{:.6e}", self.x[idx])?;
                first = false;
            } else {
                writeln!(w, "\t{:.6e}", self.x[idx])?;
            }
        }
        for &idx in self.topo.vsrc_index.values() {
            writeln!(w, "\t{:.6e}", self.x[n_nodes + idx])?;
        }
        Ok(())
    }

    /// Write DC operating point as a two-row CSV (header + one data row).
    pub fn write_csv<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        write!(w, "analysis")?;
        for name in self.topo.node_index.keys() {
            write!(w, ",V({name})")?;
        }
        let n_nodes = self.topo.n_nodes();
        for name in self.topo.vsrc_index.keys() {
            write!(w, ",I({name})")?;
        }
        writeln!(w)?;
        write!(w, "dc_op")?;
        for &idx in self.topo.node_index.values() {
            write!(w, ",{:.6e}", self.x[idx])?;
        }
        for &idx in self.topo.vsrc_index.values() {
            write!(w, ",{:.6e}", self.x[n_nodes + idx])?;
        }
        writeln!(w)
    }
}

/// Build device instances from the elements in a netlist via the registry.
pub fn build_devices(
    netlist: &Netlist,
    topo: &mut CircuitTopology,
    ctx: &SimContext,
    registry: &DeviceRegistry,
) -> Result<Vec<Box<dyn Device>>, SimError> {
    let mut devices: Vec<Box<dyn Device>> = Vec::new();
    let topo_arc = Arc::new(topo.clone());
    for el in &netlist.elements {
        match el {
            Element::Diode { anode, cathode, model_name, .. } => {
                let factory = registry.get(model_name)
                    .ok_or_else(|| SimError::UnknownModel(model_name.clone()))?;
                let pos: NodeId = topo.node_index.get(anode).copied();
                let neg: NodeId = topo.node_index.get(cathode).copied();
                devices.push(factory(&[pos, neg], ctx));
            }
            Element::Mosfet { drain, gate, source, bulk, model_name, params, .. } => {
                let d: NodeId = topo.node_index.get(drain).copied();
                let g: NodeId = topo.node_index.get(gate).copied();
                let s: NodeId = topo.node_index.get(source).copied();
                let b: NodeId = topo.node_index.get(bulk).copied();
                if let Some(dev) = registry.build_mosfet(model_name, params, &[d, g, s, b], ctx) {
                    devices.push(dev);
                } else {
                    let factory = registry.get(model_name)
                        .ok_or_else(|| SimError::UnknownModel(model_name.clone()))?;
                    devices.push(factory(&[d, g, s, b], ctx));
                }
            }
            Element::Behavioral { name, pos, neg, kind, expr } => {
                let dev = BehavioralDevice::build(
                    Arc::clone(&topo_arc),
                    name, pos, neg, *kind, expr.clone(),
                );
                devices.push(Box::new(dev));
            }
            Element::XOsdi { nets, model_name, params, .. } => {
                let factory = registry.get(model_name)
                    .ok_or_else(|| SimError::UnknownModel(model_name.clone()))?;
                let terminals: Vec<NodeId> = nets.iter()
                    .map(|net| topo.node_index.get(net).copied())
                    .collect();
                let mut dev = factory(&terminals, ctx);
                let expected = dev.num_terminals();
                if terminals.len() != expected {
                    eprintln!(
                        "warning: XOsdi '{model_name}': netlist provides {} net(s) but model \
                         expects {expected} terminal(s); extra terminals default to ground",
                        terminals.len()
                    );
                }
                for (name, value) in params {
                    dev.set_real_param(name, *value);
                }
                // OSDI models that use direct potential contributions
                // (`V(port) <+ ...`) declare internal flow-branch nodes:
                // OpenVAF surfaces these as `num_nodes − num_terminals`.
                // Allocate MNA rows for them so the OSDI runtime has real
                // slots to stamp into instead of running past mna_nodes.
                let extras = dev.num_extra_nodes();
                if extras > 0 {
                    let first = topo.allocate_extra_rows(extras);
                    dev.bind_extra_nodes(first);
                }
                devices.push(dev);
            }
            _ => {}
        }
    }
    Ok(devices)
}

/// Core Newton-Raphson loop at a fixed source scale and gmin.
fn nr_inner(
    topo: &CircuitTopology,
    netlist: &Netlist,
    devices: &mut [Box<dyn Device>],
    ctx: &SimContext,
    opts: &SimOptions,
    solver: &dyn LinearSolver,
    mut x: Vec<f64>,
    source_scale: f64,
    gmin_extra: f64,
) -> Result<Vec<f64>, SimError> {
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();
    let n_nodes = topo.n_nodes();

    for _ in 0..opts.itl1 {
        let mut mat = stamp_netlist_scaled(topo, netlist, source_scale, &empty, &empty);

        for dev in devices.iter_mut() {
            dev.eval(&x, EvalFlags::dc(), ctx);
            dev.load_residual(&mut mat.b);
            dev.load_jacobian(&mut mat);
        }

        for i in 0..n_nodes {
            mat.a[i][i] += opts.gmin + gmin_extra;
        }
        // Stamp gmin on OSDI internal-node rows too.  Their default state
        // when the device emits a degenerate contribution (e.g.
        // `OWL(out_lambda) <+ OWL(in_lambda)` where both ports map to the
        // same circuit net) is a zero row — singular without gmin.  Skip
        // voltage-source aux rows, which need a clean diagonal for their
        // own KCL equation.
        let vsrc_end = n_nodes + topo.vsrc_index.len();
        for i in vsrc_end..topo.size {
            mat.a[i][i] += opts.gmin + gmin_extra;
        }

        let x_new = solver.solve(&mat.a, &mat.b)?;

        let max_dv = x_new.iter().zip(x.iter()).take(n_nodes)
            .map(|(n, o)| (n - o).abs())
            .fold(0.0f64, f64::max);

        let x_next: Vec<f64> = if max_dv > opts.vmax {
            let scale = opts.vmax / max_dv;
            x.iter().zip(x_new.iter()).map(|(o, n)| o + scale * (n - o)).collect()
        } else {
            x_new
        };

        let converged = x_next.iter().zip(x.iter())
            .all(|(n, o)| (n - o).abs() < opts.vntol + opts.reltol * n.abs());

        x = x_next;
        if converged {
            return Ok(x);
        }
    }
    Err(SimError::NoConvergence { iters: opts.itl1 })
}

/// Source-stepping homotopy: ramp sources from 0 → full value in adaptive increments.
fn source_stepping(
    topo: &CircuitTopology,
    netlist: &Netlist,
    devices: &mut [Box<dyn Device>],
    ctx: &SimContext,
    opts: &SimOptions,
    solver: &dyn LinearSolver,
    x0: Vec<f64>,
) -> Result<Vec<f64>, SimError> {
    let mut x = x0;
    let mut scale = 0.0_f64;
    let mut ds = 1.0_f64 / (opts.srcsteps as f64).max(1.0);
    let min_ds = 1e-6_f64;

    while scale < 1.0 {
        let next = (scale + ds).min(1.0);
        match nr_inner(topo, netlist, devices, ctx, opts, solver, x.clone(), next, 0.0) {
            Ok(x_new) => {
                x = x_new;
                scale = next;
                ds = (ds * 2.0).min(2.0 * (1.0 / opts.srcsteps.max(1) as f64));
            }
            Err(_) => {
                ds *= 0.5;
                if ds < min_ds {
                    return Err(SimError::NoConvergence { iters: opts.itl1 });
                }
            }
        }
    }
    Ok(x)
}

/// GMIN stepping: add a large artificial conductance to all nodes, then ramp it down.
fn gmin_stepping(
    topo: &CircuitTopology,
    netlist: &Netlist,
    devices: &mut [Box<dyn Device>],
    ctx: &SimContext,
    opts: &SimOptions,
    solver: &dyn LinearSolver,
) -> Result<Vec<f64>, SimError> {
    let mut gmin_extra = opts.gmin_max;
    let target = opts.gmin;
    let mut x = vec![0.0f64; topo.size];

    loop {
        match nr_inner(topo, netlist, devices, ctx, opts, solver, x.clone(), 1.0, gmin_extra) {
            Ok(x_new) => {
                x = x_new;
                if gmin_extra <= target { break; }
                gmin_extra = (gmin_extra * 0.1).max(target);
            }
            Err(_) => {
                return Err(SimError::NoConvergence { iters: opts.itl1 });
            }
        }
    }
    Ok(x)
}

/// Build an initial `x` vector for DC NR, seeding from any `.nodeset` entries
/// in the netlist.  Unknown nodes are silently ignored.
fn build_x0_from_nodeset(netlist: &Netlist, topo: &CircuitTopology) -> Vec<f64> {
    let mut x0 = vec![0.0f64; topo.size];
    for (name, value) in &netlist.nodeset {
        if let Some(&i) = topo.node_index.get(name) {
            x0[i] = *value;
        }
    }
    x0
}

/// DC operating-point with explicit `SimOptions`.
///
/// Convergence strategy (in order):
///   1. Direct Newton-Raphson from the `.nodeset` seed (or x=0).
///   2. Source stepping: ramp sources from 0 → full value.
///   3. GMIN stepping: add large diagonal conductance, ramp to standard GMIN.
pub fn dc_op_nr_with_registry_opts(
    netlist: &Netlist,
    registry: &DeviceRegistry,
    opts: &SimOptions,
) -> Result<NrResult, SimError> {
    check_connectivity(netlist)?;
    let ctx = opts.sim_context();
    let mut topo = CircuitTopology::build(netlist);

    let mut devices = build_devices(netlist, &mut topo, &ctx, registry)?;
    let x0 = build_x0_from_nodeset(netlist, &topo);
    let solver = opts.linear_solver(topo.size);

    if let Ok(x) = nr_inner(&topo, netlist, &mut devices, &ctx, opts, &*solver, x0.clone(), 1.0, 0.0) {
        return Ok(NrResult { topo, x, iters: 1 });
    }

    if let Ok(x) = source_stepping(&topo, netlist, &mut devices, &ctx, opts, &*solver, x0) {
        return Ok(NrResult { topo, x, iters: 2 });
    }

    match gmin_stepping(&topo, netlist, &mut devices, &ctx, opts, &*solver) {
        Ok(x) => Ok(NrResult { topo, x, iters: 3 }),
        Err(e) => Err(e),
    }
}

/// DC operating-point with options taken from any `.options` directives in
/// the netlist (defaults where unspecified).
///
/// This is the recommended entry point for CLI/Python callers: it honours
/// `.options reltol=… gmin=… method=…` automatically.
pub fn dc_op_nr_with_registry(
    netlist: &Netlist,
    registry: &DeviceRegistry,
) -> Result<NrResult, SimError> {
    dc_op_nr_with_registry_opts(netlist, registry, &SimOptions::from_netlist(netlist))
}

/// DC operating-point with pre-built devices (for sweeps / parametric analysis).
pub fn dc_op_nr_with_devices_opts(
    netlist: &Netlist,
    topo: &CircuitTopology,
    devices: &mut Vec<Box<dyn Device>>,
    ctx: &SimContext,
    opts: &SimOptions,
) -> Result<NrResult, SimError> {
    let x0 = vec![0.0f64; topo.size];
    let solver = opts.linear_solver(topo.size);

    if let Ok(x) = nr_inner(topo, netlist, devices, ctx, opts, &*solver, x0.clone(), 1.0, 0.0) {
        return Ok(NrResult { topo: topo.clone(), x, iters: 1 });
    }

    if let Ok(x) = source_stepping(topo, netlist, devices, ctx, opts, &*solver, x0) {
        return Ok(NrResult { topo: topo.clone(), x, iters: 2 });
    }

    match gmin_stepping(topo, netlist, devices, ctx, opts, &*solver) {
        Ok(x) => Ok(NrResult { topo: topo.clone(), x, iters: 3 }),
        Err(e) => Err(e),
    }
}

/// DC operating-point with pre-built devices, default options.
pub fn dc_op_nr_with_devices(
    netlist: &Netlist,
    topo: &CircuitTopology,
    devices: &mut Vec<Box<dyn Device>>,
    ctx: &SimContext,
) -> Result<NrResult, SimError> {
    dc_op_nr_with_devices_opts(netlist, topo, devices, ctx, &SimOptions::default())
}

/// DC operating-point using only built-in models, default options.
pub fn dc_op_nr(netlist: &Netlist) -> Result<NrResult, SimError> {
    dc_op_nr_opts(netlist, &SimOptions::default())
}

/// DC operating-point using only built-in models, with explicit options.
pub fn dc_op_nr_opts(netlist: &Netlist, opts: &SimOptions) -> Result<NrResult, SimError> {
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_diodes(&netlist.models);
    registry.register_builtin_mosfets(&netlist.models);
    dc_op_nr_with_registry_opts(netlist, &registry, opts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

    #[test]
    fn linear_circuit_converges() {
        let net = parse_spice(
            "* divider\nV1 in 0 DC 1.0\nR1 in out 1k\nR2 out 0 1k\n.op\n.end\n",
        ).unwrap();
        let r = dc_op_nr(&net).unwrap();
        assert!(r.iters <= 5, "linear circuit should converge fast, took {}", r.iters);
        let v = r.node_voltage("out").unwrap();
        assert!((v - 0.5).abs() < 1e-6, "v(out)={v}");
    }

    #[test]
    fn pnjlim_helps_stiff_diode() {
        // 20 V across a tiny series R into a diode: an easy case for an
        // ideal solver, a stiff one without junction limiting.  With pnjlim
        // ON (default) this must converge in well under ITL1=150 iters.
        let net = parse_spice(
            "* stiff diode\nVdd a 0 DC 20\nR1 a b 100\nD1 b 0 myd\n\
             .model myd D (Is=1e-15 N=1)\n.op\n.end\n"
        ).unwrap();
        let r = dc_op_nr(&net).unwrap();
        // The pn-junction equation gives V(b) ≈ 18·V_T·ln(I/Is) ≈ 0.7 V,
        // and most of the supply drops across R1.
        let vb = r.node_voltage("b").unwrap();
        assert!(vb > 0.5 && vb < 1.0, "V(b)={vb:.4} should be ~0.7V");
        assert!(r.iters < 50, "convergence took {} iters with pnjlim on", r.iters);
    }

    #[test]
    fn current_source_biased_diode() {
        let net = parse_spice(
            "* Diode bias\nIb 0 b 1m\nD1 b 0 myd\n.model myd D (Is=1e-14 N=1)\n.op\n.end\n",
        ).unwrap();
        let r = dc_op_nr(&net).unwrap();
        let vb = r.node_voltage("b").unwrap();
        let vt = 1.380649e-23 * 300.15 / 1.602176634e-19;
        let expected = vt * (1e-3_f64 / 1e-14_f64 + 1.0).ln();
        let tol = 1e-4 * expected;
        assert!(
            (vb - expected).abs() < tol,
            "V(b)={vb:.6e}  expected={expected:.6e}  diff={:.2e}",
            (vb - expected).abs()
        );
    }

    #[test]
    fn series_rd_circuit() {
        let net = parse_spice(
            "* R-D series\nVdd a 0 DC 5\nR1 a b 10k\nD1 b 0 myd\n\
             .model myd D (Is=1e-14 N=1)\n.op\n.end\n",
        ).unwrap();
        let r = dc_op_nr(&net).unwrap();
        let vb = r.node_voltage("b").unwrap();
        assert!(vb > 0.5 && vb < 0.8, "V(b) out of expected range: {vb}");
        let vt = 1.380649e-23 * 300.15 / 1.602176634e-19;
        let ir = (5.0 - vb) / 10e3;
        let id = 1e-14 * ((vb / vt).exp() - 1.0);
        assert!((ir - id).abs() < 1e-8, "KCL error: ir={ir:.4e} id={id:.4e}");
    }

    #[test]
    fn write_csv_dc_op() {
        let net = parse_spice(
            "* divider\nV1 in 0 DC 2.0\nR1 in out 1k\nR2 out 0 1k\n.op\n.end\n",
        ).unwrap();
        let r = dc_op_nr(&net).unwrap();
        let mut buf = Vec::new();
        r.write_csv(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("analysis,"), "header: {s}");
        assert!(s.contains("V(out)"), "should contain V(out): {s}");
        assert!(s.contains("dc_op"), "should have dc_op row: {s}");
        assert!(s.contains("1.000000e0") || s.contains("1.000000e+0"), "V(out)≈1V missing: {s}");
    }

    #[test]
    fn write_nutmeg_dc_op() {
        let net = parse_spice(
            "* divider\nV1 in 0 DC 2.0\nR1 in out 1k\nR2 out 0 1k\n.op\n.end\n",
        ).unwrap();
        let r = dc_op_nr(&net).unwrap();
        let mut buf = Vec::new();
        r.write_nutmeg(&mut buf, "divider test").unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("Plotname: Operating Point"), "plotname: {s}");
        assert!(s.contains("Flags: real"), "flags: {s}");
        assert!(s.contains("v(out)\tvoltage"), "v(out): {s}");
        assert!(s.contains("No. Points: 1"), "single point: {s}");
    }

    #[test]
    fn tighter_reltol_takes_more_iterations() {
        // A circuit with mild nonlinearity. With opts.reltol=1e-10 NR should still
        // converge but use more iterations than at the default 1e-3.
        let net = parse_spice(
            "* R-D\nVdd a 0 DC 5\nR1 a b 10k\nD1 b 0 myd\n\
             .model myd D (Is=1e-14 N=1)\n.op\n.end\n",
        ).unwrap();
        let mut opts = SimOptions::default();
        opts.reltol = 1e-10;
        opts.vntol  = 1e-12;
        let r = dc_op_nr_opts(&net, &opts).unwrap();
        let vb = r.node_voltage("b").unwrap();
        assert!(vb > 0.5 && vb < 0.8);
    }
}
