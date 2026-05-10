use indexmap::IndexMap;

use fairchild_parser::{Element, Netlist};

use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::mna::{stamp_netlist, CircuitTopology};
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

/// DC operating-point with a pre-built registry (supports OSDI and built-in models).
pub fn dc_op_nr_with_registry(
    netlist: &Netlist,
    registry: &DeviceRegistry,
) -> Result<NrResult, SimError> {
    let ctx = SimContext::default();
    let topo = CircuitTopology::build(netlist);
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();

    let mut devices = build_devices(netlist, &topo, &ctx, registry)?;
    let mut x = vec![0.0f64; topo.size];

    for iter in 0..MAX_ITER {
        let mut mat = stamp_netlist(&topo, netlist, 0.0, &empty, &empty);

        for dev in &mut devices {
            dev.eval(&x, EvalFlags::dc(), &ctx);
            dev.load_residual(&mut mat.b);
            dev.load_jacobian(&mut mat);
        }

        // GMIN: add a small conductance to each node diagonal to prevent
        // near-singular matrices when all devices have tiny conductance (e.g.
        // a diode at Vd=0 has gd ≈ Is/Vt ≈ 4e-13 S ≪ GMIN).
        for i in 0..topo.n_nodes() {
            mat.a[i][i] += GMIN;
        }

        let x_new = lu_solve(&mat.a, &mat.b)?;

        let max_dv = x_new.iter()
            .zip(x.iter())
            .take(topo.n_nodes())
            .map(|(n, o)| (n - o).abs())
            .fold(0.0f64, f64::max);

        let x_next: Vec<f64> = if max_dv > VMAX {
            let scale = VMAX / max_dv;
            x.iter().zip(x_new.iter()).map(|(o, n)| o + scale * (n - o)).collect()
        } else {
            x_new
        };

        let converged = x_next.iter()
            .zip(x.iter())
            .all(|(n, o)| (n - o).abs() < VNTOL + RELTOL * n.abs());

        x = x_next;

        if converged {
            return Ok(NrResult { topo, x, iters: iter + 1 });
        }
    }

    Err(SimError::NoConvergence { iters: MAX_ITER })
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
        // Sanity: voltage at diode anode should be in a forward-bias range.
        assert!(vb > 0.5 && vb < 0.8, "V(b) out of expected range: {vb}");
        // KCL: (5 - vb)/10k ≈ Is*(exp(vb/Vt)-1)
        let vt = 1.380649e-23 * 300.15 / 1.602176634e-19;
        let ir = (5.0 - vb) / 10e3;
        let id = 1e-14 * ((vb / vt).exp() - 1.0);
        assert!((ir - id).abs() < 1e-8, "KCL error: ir={ir:.4e} id={id:.4e}");
    }
}
