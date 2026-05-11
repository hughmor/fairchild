use indexmap::IndexMap;

use fairchild_parser::{Element, Netlist};

use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::mna::{stamp_netlist_scaled, CircuitTopology};
use crate::solver::lu_solve;

/// SPICE-standard convergence tolerances.
pub(crate) const VNTOL: f64 = 1e-6;     // 1 μV absolute floor for voltages
pub(crate) const RELTOL: f64 = 1e-3;    // 0.1 % relative tolerance
/// Backup global step damping: max allowed Δ(node voltage) per iteration.
pub(crate) const VMAX: f64 = 0.5;
pub(crate) const MAX_ITER: usize = 150;
/// Minimum conductance added to the diagonal to prevent near-singular matrices.
/// Standard SPICE heuristic (ngspice default: 1e-12 S).
pub(crate) const GMIN: f64 = 1e-12;

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

/// Build device instances from the Diode elements in a netlist via the registry.
///
/// Called by both `dc_op_nr` and `tran_nr`.
pub(crate) fn build_devices(
    netlist: &Netlist,
    topo: &CircuitTopology,
    ctx: &SimContext,
    registry: &DeviceRegistry,
) -> Result<Vec<Box<dyn Device>>, SimError> {
    let mut devices: Vec<Box<dyn Device>> = Vec::new();
    for el in &netlist.elements {
        match el {
            Element::Diode { anode, cathode, model_name, .. } => {
                let factory = registry.get(model_name)
                    .ok_or_else(|| SimError::UnknownModel(model_name.clone()))?;
                let pos: NodeId = topo.node_index.get(anode).copied();
                let neg: NodeId = topo.node_index.get(cathode).copied();
                devices.push(factory(&[pos, neg], ctx));
            }
            Element::Mosfet { drain, gate, source, bulk, model_name, .. } => {
                let factory = registry.get(model_name)
                    .ok_or_else(|| SimError::UnknownModel(model_name.clone()))?;
                let d: NodeId = topo.node_index.get(drain).copied();
                let g: NodeId = topo.node_index.get(gate).copied();
                let s: NodeId = topo.node_index.get(source).copied();
                let b: NodeId = topo.node_index.get(bulk).copied();
                devices.push(factory(&[d, g, s, b], ctx));
            }
            _ => {}
        }
    }
    Ok(devices)
}

/// Core Newton-Raphson loop at a fixed source scale and gmin.
///
/// `source_scale` ∈ [0,1]: scales all independent source amplitudes.
/// `gmin_extra`: extra diagonal conductance added to every node (for GMIN stepping).
/// Returns Ok(x) if converged within MAX_ITER, Err if not.
fn nr_inner(
    topo: &CircuitTopology,
    netlist: &Netlist,
    devices: &mut [Box<dyn Device>],
    ctx: &SimContext,
    mut x: Vec<f64>,
    source_scale: f64,
    gmin_extra: f64,
) -> Result<Vec<f64>, SimError> {
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();
    let n_nodes = topo.n_nodes();

    for _ in 0..MAX_ITER {
        let mut mat = stamp_netlist_scaled(topo, netlist, source_scale, &empty, &empty);

        for dev in devices.iter_mut() {
            dev.eval(&x, EvalFlags::dc(), ctx);
            dev.load_residual(&mut mat.b);
            dev.load_jacobian(&mut mat);
        }

        for i in 0..n_nodes {
            mat.a[i][i] += GMIN + gmin_extra;
        }

        let x_new = lu_solve(&mat.a, &mat.b)?;

        let max_dv = x_new.iter().zip(x.iter()).take(n_nodes)
            .map(|(n, o)| (n - o).abs())
            .fold(0.0f64, f64::max);

        let x_next: Vec<f64> = if max_dv > VMAX {
            let scale = VMAX / max_dv;
            x.iter().zip(x_new.iter()).map(|(o, n)| o + scale * (n - o)).collect()
        } else {
            x_new
        };

        let converged = x_next.iter().zip(x.iter())
            .all(|(n, o)| (n - o).abs() < VNTOL + RELTOL * n.abs());

        x = x_next;
        if converged {
            return Ok(x);
        }
    }
    Err(SimError::NoConvergence { iters: MAX_ITER })
}

/// Source-stepping homotopy: ramp sources from 0 → full value in at most `n_steps` increments.
///
/// Returns Ok(x) at full source values if the homotopy converges; Err if it can't even
/// converge with tiny source steps.
fn source_stepping(
    topo: &CircuitTopology,
    netlist: &Netlist,
    devices: &mut [Box<dyn Device>],
    ctx: &SimContext,
    x0: Vec<f64>,
) -> Result<Vec<f64>, SimError> {
    let mut x = x0;
    let mut scale = 0.0_f64;
    let mut ds = 0.1_f64;
    let min_ds = 1e-6_f64;

    while scale < 1.0 {
        let next = (scale + ds).min(1.0);
        match nr_inner(topo, netlist, devices, ctx, x.clone(), next, 0.0) {
            Ok(x_new) => {
                x = x_new;
                scale = next;
                ds = (ds * 2.0).min(0.2);
            }
            Err(_) => {
                ds *= 0.5;
                if ds < min_ds {
                    return Err(SimError::NoConvergence { iters: MAX_ITER });
                }
            }
        }
    }
    Ok(x)
}

/// GMIN stepping: add a large artificial conductance to all nodes, solve, then ramp it
/// down to the standard GMIN over logarithmic steps.
fn gmin_stepping(
    topo: &CircuitTopology,
    netlist: &Netlist,
    devices: &mut [Box<dyn Device>],
    ctx: &SimContext,
) -> Result<Vec<f64>, SimError> {
    let mut gmin_extra = 1.0_f64;  // Start at 1 S (≈ 1 Ω across every node)
    let target = GMIN;
    let mut x = vec![0.0f64; topo.size];

    // Ramp down GMIN from 1 S to GMIN over ~12 decades in steps of ÷10.
    loop {
        match nr_inner(topo, netlist, devices, ctx, x.clone(), 1.0, gmin_extra) {
            Ok(x_new) => {
                x = x_new;
                if gmin_extra <= target { break; }
                gmin_extra = (gmin_extra * 0.1).max(target);
            }
            Err(_) => {
                return Err(SimError::NoConvergence { iters: MAX_ITER });
            }
        }
    }
    Ok(x)
}

/// DC operating-point with a pre-built registry (supports OSDI and built-in models).
///
/// Convergence strategy (in order):
///   1. Direct Newton-Raphson from x=0.
///   2. Source stepping: ramp sources from 0 → full value.
///   3. GMIN stepping: add large diagonal conductance, ramp to standard GMIN.
pub fn dc_op_nr_with_registry(
    netlist: &Netlist,
    registry: &DeviceRegistry,
) -> Result<NrResult, SimError> {
    let ctx = SimContext::default();
    let topo = CircuitTopology::build(netlist);

    let mut devices = build_devices(netlist, &topo, &ctx, registry)?;
    let x0 = vec![0.0f64; topo.size];

    // Strategy 1: direct NR.
    if let Ok(x) = nr_inner(&topo, netlist, &mut devices, &ctx, x0.clone(), 1.0, 0.0) {
        return Ok(NrResult { topo, x, iters: 1 });
    }

    // Strategy 2: source stepping.
    if let Ok(x) = source_stepping(&topo, netlist, &mut devices, &ctx, x0) {
        return Ok(NrResult { topo, x, iters: 2 });
    }

    // Strategy 3: GMIN stepping.
    match gmin_stepping(&topo, netlist, &mut devices, &ctx) {
        Ok(x) => Ok(NrResult { topo, x, iters: 3 }),
        Err(e) => Err(e),
    }
}

/// DC operating-point using only built-in models from `.model` cards.
pub fn dc_op_nr(netlist: &Netlist) -> Result<NrResult, SimError> {
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_diodes(&netlist.models);
    dc_op_nr_with_registry(netlist, &registry)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

    /// A purely linear circuit should converge quickly (no nonlinear device).
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

    /// Current-source biased diode: V(b) = N·Vt·ln(Ib/Is).
    #[test]
    fn current_source_biased_diode() {
        // 1 mA into D1 b 0; Is=1e-14, N=1.
        // Expected: V(b) = Vt * ln(1e-3/1e-14) ≈ 0.025852 * 25.328 ≈ 0.6549 V
        let net = parse_spice(
            "* Diode bias\nIb 0 b 1m\nD1 b 0 myd\n.model myd D (Is=1e-14 N=1)\n.op\n.end\n",
        ).unwrap();
        let r = dc_op_nr(&net).unwrap();
        let vb = r.node_voltage("b").unwrap();
        let vt = 1.380649e-23 * 300.15 / 1.602176634e-19;
        let expected = vt * (1e-3_f64 / 1e-14_f64 + 1.0).ln();
        let tol = 1e-4 * expected;  // 0.01 % relative
        assert!(
            (vb - expected).abs() < tol,
            "V(b)={vb:.6e}  expected={expected:.6e}  diff={:.2e}",
            (vb - expected).abs()
        );
    }

    /// Series R-D circuit: Vdd=5V, R=10k, D biased by resistive divider.
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
        // V(out) should be ~1.0 V (voltage divider).
        assert!(s.contains("1.000000e0") || s.contains("1.000000e+0"), "V(out)≈1V missing: {s}");
    }
}
