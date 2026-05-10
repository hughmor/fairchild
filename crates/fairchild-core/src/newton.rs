use indexmap::IndexMap;

use fairchild_parser::{Element, Netlist};

use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::error::SimError;
use crate::mna::{stamp_netlist, CircuitTopology};
use crate::models::ShockleyDiode;
use crate::solver::lu_solve;

/// SPICE-standard convergence tolerances.
const VNTOL: f64 = 1e-6;     // 1 μV absolute floor for voltages
const RELTOL: f64 = 1e-3;    // 0.1 % relative tolerance
/// Backup global step damping: max allowed Δ(node voltage) per iteration.
const VMAX: f64 = 0.5;
const MAX_ITER: usize = 150;

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

/// DC operating-point with Newton-Raphson for nonlinear devices (diodes, etc.).
///
/// Falls back gracefully for purely linear circuits (no Diode elements):
/// runs one NR iteration, which converges in a single step.
pub fn dc_op_nr(netlist: &Netlist) -> Result<NrResult, SimError> {
    let ctx = SimContext::default();
    let topo = CircuitTopology::build(netlist);
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();

    // Build device list from Diode elements.
    let mut devices: Vec<Box<dyn Device>> = Vec::new();
    for el in &netlist.elements {
        if let Element::Diode { anode, cathode, model_name, .. } = el {
            let card = netlist.models.iter()
                .find(|m| m.name == *model_name && m.kind.starts_with('d'))
                .ok_or_else(|| SimError::UnknownModel(model_name.clone()))?;

            let mut dev = Box::new(ShockleyDiode::from_params(&card.params));
            dev.setup_model(&ctx);

            let pos: NodeId = topo.node_index.get(anode).copied();
            let neg: NodeId = topo.node_index.get(cathode).copied();
            dev.setup_instance(&[pos, neg], &ctx);
            devices.push(dev);
        }
    }

    let mut x = vec![0.0f64; topo.size];

    for iter in 0..MAX_ITER {
        let mut mat = stamp_netlist(&topo, netlist, 0.0, &empty, &empty);

        // Evaluate and stamp nonlinear devices.
        for dev in &mut devices {
            dev.eval(&x, EvalFlags::dc(), &ctx);
            dev.load_residual(&mut mat.b);
            dev.load_jacobian(&mut mat);
        }

        let x_new = lu_solve(&mat.a, &mat.b)?;

        // Global voltage step damping (backup — pnjlim handles most cases).
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

        // Convergence: |Δx[i]| < VNTOL + RELTOL·|x_next[i]| for all i.
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
