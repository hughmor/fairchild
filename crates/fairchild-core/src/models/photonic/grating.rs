use super::stamp_potential_eq;
use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

// ────────────────────────────────────────────────────────────────────────
// Native grating coupler — flat insertion loss, zero length
// ────────────────────────────────────────────────────────────────────────

/// Grating coupler (fibre ↔ chip).  Zero-length waveguide with a flat
/// amplitude attenuation set by `alpha_db` (insertion loss).  Variable-arity
/// bundle-aware: 6·N terminals.
pub struct NativeGratingCoupler {
    alpha_db: f64,
    n_channels: usize,
    wpc: usize,
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
}

impl Default for NativeGratingCoupler {
    fn default() -> Self {
        Self::new()
    }
}

impl NativeGratingCoupler {
    pub fn new() -> Self {
        Self {
            alpha_db: 3.0,
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
        }
    }
}

impl Device for NativeGratingCoupler {
    fn num_terminals(&self) -> usize {
        self.nodes.len()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 2 * wpc;
        assert!(
            !terminals.is_empty() && terminals.len().is_multiple_of(stride),
            "fc_grating_coupler: terminal count must be {stride}·N (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = terminals.len() / stride;
        self.n_channels = n;
        self.nodes = terminals.to_vec();
        let bpc = if wpc == 5 { 4 } else { 2 };
        self.branches = vec![None; bpc * n];
    }

    fn num_extra_nodes(&self) -> usize {
        self.branches.len()
    }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..self.branches.len() {
            self.branches[i] = Some(first_idx + i);
        }
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "alpha_db" | "alpha_db_il" | "il_db" => {
                self.alpha_db = value;
                true
            }
            "alpha" => {
                let t = value.max(1e-30);
                self.alpha_db = -20.0 * t.log10();
                true
            }
            _ => false,
        }
    }

    /// Straight through, channel for channel.
    fn lambda_routing(&self) -> Vec<(usize, usize)> {
        let (wpc, n) = (self.wpc, self.n_channels);
        let lam = wpc - 1;
        (0..n)
            .map(|k| (wpc * k + lam, wpc * n + wpc * k + lam))
            .collect()
    }

    fn eval(&mut self, _x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {}

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let bpc = if wpc == 5 { 4 } else { 2 };
        let out_base = wpc * n;
        let t = 10f64.powf(-self.alpha_db / 20.0);
        for k in 0..n {
            let in_re_fw = self.nodes[wpc * k];
            let in_im_fw = self.nodes[wpc * k + 1];
            let out_re_fw = self.nodes[out_base + wpc * k];
            let out_im_fw = self.nodes[out_base + wpc * k + 1];
            stamp_potential_eq(mat, &self.branches, bpc * k, out_re_fw, &[(in_re_fw, -t)]);
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k + 1,
                out_im_fw,
                &[(in_im_fw, -t)],
            );
            if wpc == 5 {
                let in_re_bw = self.nodes[wpc * k + 2];
                let in_im_bw = self.nodes[wpc * k + 3];
                let out_re_bw = self.nodes[out_base + wpc * k + 2];
                let out_im_bw = self.nodes[out_base + wpc * k + 3];
                stamp_potential_eq(
                    mat,
                    &self.branches,
                    bpc * k + 2,
                    in_re_bw,
                    &[(out_re_bw, -t)],
                );
                stamp_potential_eq(
                    mat,
                    &self.branches,
                    bpc * k + 3,
                    in_im_bw,
                    &[(out_im_bw, -t)],
                );
            }
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }
}
