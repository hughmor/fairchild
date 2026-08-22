//! Behavioural (`B`-element) device: a nonlinear V- or I-source whose value
//! is a user-supplied expression evaluated on each Newton iteration.
//!
//! The Jacobian is computed by numerical differentiation: for each node /
//! branch the expression references, we perturb the corresponding entry of
//! `x`, re-evaluate, and take the central difference.  This is correct but
//! O(K · cost(expr)) per NR iteration where K = number of distinct refs.
//! Acceptable for the typical B-element with 1–3 references.

use std::sync::Arc;

use fairchild_parser::{BehavioralKind, EvalContext, Expr};

use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::{CircuitTopology, MnaMatrix};

/// One expression reference, resolved into MNA-vector coordinates.
#[derive(Clone)]
enum RefKind {
    /// `V(node)` — read `x[idx]` directly.  None means ground (always 0).
    NodeV(NodeId),
    /// `I(vsrc_or_bsrc)` — read `x[n_nodes + idx]`.
    BranchI(usize),
    /// `TIME` — read from the device context (no x-dependence; constant
    /// Jacobian wrt x).  Reserved for the future `time` expression
    /// reference; the match arms below already handle it correctly when
    /// the parser begins emitting it.
    #[allow(dead_code)]
    Time,
}

/// Evaluation context that reads V/I from a snapshot of `x`.
struct XContext<'a> {
    x: &'a [f64],
    topo: &'a CircuitTopology,
    time: f64,
}

impl<'a> EvalContext for XContext<'a> {
    fn node_voltage(&self, node: &str) -> f64 {
        if node == "0" || node == "gnd" {
            return 0.0;
        }
        self.topo
            .node_index
            .get(node)
            .map(|&i| self.x[i])
            .unwrap_or(0.0)
    }
    fn branch_current(&self, vsrc: &str) -> f64 {
        let n = self.topo.n_nodes();
        self.topo
            .vsrc_index
            .get(vsrc)
            .map(|&i| self.x[n + i])
            .unwrap_or(0.0)
    }
    fn time(&self) -> f64 {
        self.time
    }
}

/// A B-element device.
pub struct BehavioralDevice {
    kind: BehavioralKind,
    pos: NodeId,
    neg: NodeId,
    /// Aux MNA row for V= form (None for I=).
    vi: Option<usize>,
    /// Parsed expression.
    expr: Arc<Expr>,
    /// Topology snapshot (cheap: just an Arc-shared reference for name → index lookups).
    topo: Arc<CircuitTopology>,
    /// Distinct x-coordinates the expression depends on.  Used for numerical
    /// Jacobian: perturb each of these and re-evaluate.
    refs: Vec<RefKind>,
    /// Cached `expr(x)` and `d(expr)/d(x[k])` from the last `eval`.
    value: f64,
    grads: Vec<f64>,
    /// Norton-equivalent: `i_eq = f(x_old) − Σ_k g_k · x_old[col_k]`.
    /// Stamped into b at residual time; combined with `grads` in the Jacobian
    /// this gives the correct linearisation for `f(x_new)` at `x_old`.
    i_eq: f64,
    /// Cached time (set at eval time).
    t: f64,
}

impl BehavioralDevice {
    /// Build a B-element device from a parsed element + topology.
    pub fn build(
        topo: Arc<CircuitTopology>,
        name: &str,
        pos: &str,
        neg: &str,
        kind: BehavioralKind,
        expr: Expr,
    ) -> Self {
        let pos_id = if pos == "0" || pos == "gnd" {
            None
        } else {
            topo.node_index.get(pos).copied()
        };
        let neg_id = if neg == "0" || neg == "gnd" {
            None
        } else {
            topo.node_index.get(neg).copied()
        };
        let vi = if kind == BehavioralKind::Voltage {
            topo.vsrc_index.get(name).copied()
        } else {
            None
        };

        // Collect refs from the expression and resolve them.
        let mut v_nodes = Vec::new();
        let mut i_srcs = Vec::new();
        expr.collect_refs(&mut v_nodes, &mut i_srcs);
        let mut refs: Vec<RefKind> = Vec::new();
        for n in &v_nodes {
            let id = if n == "0" || n == "gnd" {
                None
            } else {
                topo.node_index.get(n).copied()
            };
            refs.push(RefKind::NodeV(id));
        }
        for s in &i_srcs {
            if let Some(idx) = topo.vsrc_index.get(s).copied() {
                refs.push(RefKind::BranchI(idx));
            }
        }
        // De-dup refs by their (kind, coord) key.
        let key = |r: &RefKind| match r {
            RefKind::NodeV(Some(i)) => (0u8, *i),
            RefKind::NodeV(None) => (0, usize::MAX),
            RefKind::BranchI(i) => (1, *i),
            RefKind::Time => (2, 0),
        };
        refs.sort_by_key(|r| key(r));
        refs.dedup_by(|a, b| key(a) == key(b));

        let grads = vec![0.0; refs.len()];
        BehavioralDevice {
            kind,
            pos: pos_id,
            neg: neg_id,
            vi,
            expr: Arc::new(expr),
            topo,
            refs,
            value: 0.0,
            grads,
            i_eq: 0.0,
            t: 0.0,
        }
    }

    fn eval_value(&self, x: &[f64], t: f64) -> f64 {
        let ctx = XContext {
            x,
            topo: &self.topo,
            time: t,
        };
        self.expr.eval(&ctx)
    }
}

impl Device for BehavioralDevice {
    fn num_terminals(&self) -> usize {
        2
    }
    fn setup_model(&mut self, _ctx: &SimContext) {}
    fn setup_instance(&mut self, _terminals: &[NodeId], _ctx: &SimContext) {}

    /// A B-element is the one device that stamps outside its own terminals:
    /// the Jacobian row picks up a column for every `V(node)` / `I(vsrc)` the
    /// expression references.  `build_devices` never sees those, so report the
    /// whole footprint here — terminals, aux row, and every referenced coord.
    fn extra_stamp_rows(&self) -> Vec<usize> {
        let n_nodes = self.topo.n_nodes();
        let mut rows: Vec<usize> = self.pos.into_iter().chain(self.neg).collect();
        if let Some(vi) = self.vi {
            rows.push(n_nodes + vi);
        }
        for r in &self.refs {
            match r {
                RefKind::NodeV(Some(i)) => rows.push(*i),
                RefKind::BranchI(i) => rows.push(n_nodes + *i),
                RefKind::NodeV(None) | RefKind::Time => {}
            }
        }
        rows
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        self.value = self.eval_value(x, self.t);

        // Numerical Jacobian: central difference for each referenced coord.
        let n_nodes = self.topo.n_nodes();
        let mut x_pert = x.to_vec();
        let mut sum_g_x = 0.0;
        for (k, r) in self.refs.iter().enumerate() {
            let idx_opt: Option<usize> = match r {
                RefKind::NodeV(Some(i)) => Some(*i),
                RefKind::NodeV(None) => None,
                RefKind::BranchI(i) => Some(n_nodes + *i),
                RefKind::Time => None,
            };
            let Some(i) = idx_opt else {
                self.grads[k] = 0.0;
                continue;
            };

            let xi = x_pert[i];
            let h = 1e-6 * (xi.abs().max(1.0));
            x_pert[i] = xi + h;
            let v_plus = self.eval_value(&x_pert, self.t);
            x_pert[i] = xi - h;
            let v_minus = self.eval_value(&x_pert, self.t);
            x_pert[i] = xi;
            let g = (v_plus - v_minus) / (2.0 * h);
            self.grads[k] = g;
            sum_g_x += g * xi;
        }
        self.i_eq = self.value - sum_g_x;
    }

    fn load_residual(&self, b: &mut [f64]) {
        let n_nodes = self.topo.n_nodes();
        match self.kind {
            BehavioralKind::Voltage => {
                // Aux row: V(pos) − V(neg) − f(x_new) = 0.
                // Linearise: V(pos) − V(neg) − Σ_k g_k·x_new[k] = i_eq
                // So stamp b[vi] += i_eq.
                if let Some(vi) = self.vi {
                    b[n_nodes + vi] += self.i_eq;
                }
            }
            BehavioralKind::Current => {
                // KCL: current OUT of pos is f(x_new) ≈ i_eq + Σ_k g_k·x_new[k].
                // Move Σ g x to A; i_eq goes to b with the "OUT" sign:
                //   b[pos] -= i_eq,  b[neg] += i_eq.
                if let Some(p) = self.pos {
                    b[p] -= self.i_eq;
                }
                if let Some(n) = self.neg {
                    b[n] += self.i_eq;
                }
            }
        }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n_nodes = self.topo.n_nodes();
        match self.kind {
            BehavioralKind::Voltage => {
                if let Some(vi) = self.vi {
                    let row = n_nodes + vi;
                    // Incidence (just like a real voltage source):
                    if let Some(p) = self.pos {
                        mat.a[p][row] += 1.0;
                        mat.a[row][p] += 1.0;
                    }
                    if let Some(n) = self.neg {
                        mat.a[n][row] -= 1.0;
                        mat.a[row][n] -= 1.0;
                    }
                    // The aux row equation V(pos) - V(neg) - Σ g_k·x_new[k] = i_eq
                    // gives -g_k entries at row=row, col=col_k.
                    for (k, r) in self.refs.iter().enumerate() {
                        let col_opt = match r {
                            RefKind::NodeV(Some(i)) => Some(*i),
                            RefKind::NodeV(None) => None,
                            RefKind::BranchI(i) => Some(n_nodes + *i),
                            RefKind::Time => None,
                        };
                        if let Some(col) = col_opt {
                            mat.a[row][col] -= self.grads[k];
                        }
                    }
                }
            }
            BehavioralKind::Current => {
                // KCL_pos coefficient of x_new[col_k] is +g_k (current OUT of pos);
                // KCL_neg coefficient is -g_k (current IN to neg).
                for (k, r) in self.refs.iter().enumerate() {
                    let col_opt = match r {
                        RefKind::NodeV(Some(i)) => Some(*i),
                        RefKind::NodeV(None) => None,
                        RefKind::BranchI(i) => Some(n_nodes + *i),
                        RefKind::Time => None,
                    };
                    if let Some(col) = col_opt {
                        if let Some(p) = self.pos {
                            mat.a[p][col] += self.grads[k];
                        }
                        if let Some(n) = self.neg {
                            mat.a[n][col] -= self.grads[k];
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use fairchild_parser::parse_spice;

    #[test]
    fn behavioral_current_source_linear_expression() {
        // B1 with pos=out, neg=0 pulls I = V(in)·1mA out of `out`.  SPICE
        // convention: positive I from + node through source to − node, so the
        // external network must source that current back into `out` via R2;
        // for R2=1k that pins V(out) = −I·R = −1V.
        let net = parse_spice(
            "* b\nV1 in 0 DC 1\nR1 in nin 1\n\
             B1 out 0 I=V(in)*1m\nR2 out 0 1k\n.op\n",
        )
        .unwrap();
        let r = crate::newton::dc_op_nr(&net).unwrap();
        let v_out = r.node_voltage("out").unwrap();
        assert!((v_out + 1.0).abs() < 1e-4, "V(out)={v_out}");
    }

    #[test]
    fn behavioral_voltage_source_simple() {
        // B1 V = V(in) * 2.  V(in) = 1, so V(out) = 2.
        let net = parse_spice(
            "* bv\nV1 in 0 DC 1\nR1 in 0 1k\n\
             B1 out 0 V=V(in)*2\nR2 out 0 1k\n.op\n",
        )
        .unwrap();
        let r = crate::newton::dc_op_nr(&net).unwrap();
        let v_out = r.node_voltage("out").unwrap();
        assert!((v_out - 2.0).abs() < 1e-4, "V(out)={v_out}");
    }

    #[test]
    fn behavioral_nonlinear_square() {
        // B1 V = V(in) ^ 2.  V(in) = 0.7, so V(out) should converge to 0.49.
        let net = parse_spice(
            "* bv\nV1 in 0 DC 0.7\nR1 in 0 1k\n\
             B1 out 0 V=V(in)^2\nR2 out 0 1k\n.op\n",
        )
        .unwrap();
        let r = crate::newton::dc_op_nr(&net).unwrap();
        let v_out = r.node_voltage("out").unwrap();
        assert!((v_out - 0.49).abs() < 1e-3, "V(out)={v_out}");
    }
}
