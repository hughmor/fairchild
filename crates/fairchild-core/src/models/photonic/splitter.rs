use super::stamp_potential_eq;
use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

// ────────────────────────────────────────────────────────────────────────
// Native 1×2 Y-junction splitter (3 dB lossless)
// ────────────────────────────────────────────────────────────────────────

/// 1×2 splitter with optional insertion loss and asymmetric power split.
///
/// Total power transmission `α` (intensity, default 1.0 = lossless) is split
/// across the two outputs with fraction `r` going to `out_a` and `α − r` going
/// to `out_b`.  Defaults (α = 1.0, r = 0.5) reproduce the original 3 dB
/// lossless splitter.  Amplitude coefficients: `k_a = √r`, `k_b = √(α − r)`.
/// Wavelength duplicated to both outputs.
/// Variable-arity bundle-aware 1×2 splitter.  Terminal layout for N channels
/// (9·N total):
///   [a.0.re, a.0.im, a.0.λ, ..., a.{N-1}.λ,
///    c.0.re, c.0.im, c.0.λ, ..., c.{N-1}.λ,
///    d.0.re, d.0.im, d.0.λ, ..., d.{N-1}.λ]
pub struct NativeSplitter {
    alpha: f64,
    r: f64,
    n_channels: usize,
    wpc: usize,
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
}

impl Default for NativeSplitter {
    fn default() -> Self {
        Self::new()
    }
}

impl NativeSplitter {
    pub fn new() -> Self {
        Self {
            alpha: 1.0,
            r: 0.5,
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
        }
    }
}

impl Device for NativeSplitter {
    fn num_terminals(&self) -> usize {
        self.nodes.len()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 3 * wpc; // 3 ports per channel
        assert!(
            !terminals.is_empty() && terminals.len().is_multiple_of(stride),
            "fc_splitter: terminal count must be {stride}·N (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = terminals.len() / stride;
        self.n_channels = n;
        self.nodes = terminals.to_vec();
        // Per channel: 4 fw branches (re, im for out_a + out_b) + 2 λ +
        // (if bidir) 2 bw branches (re, im for in_a).  Note: under bidir the
        // splitter behaves like a combiner in reverse — bw light from out_a
        // and out_b combine back into in.  re_bw_in = k_a · re_bw_a + k_b · re_bw_b.
        let bpc = if wpc == 5 { 8 } else { 6 };
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
            "alpha" => {
                // Intensity transmission, must be ≤ 1.  Re-anchor r at the
                // symmetric midpoint (alpha/2) iff the user hasn't already
                // skewed it — but simplest: if r exceeds the new alpha,
                // clamp it to alpha/2.  Otherwise leave r alone.
                self.alpha = value.clamp(0.0, 1.0);
                if self.r > self.alpha {
                    self.r = self.alpha * 0.5;
                }
                true
            }
            "alpha_db" | "il_db" => {
                self.alpha = 10f64.powf(-value / 10.0);
                if self.r > self.alpha {
                    self.r = self.alpha * 0.5;
                }
                true
            }
            "r" | "split_ratio" => {
                // r = fraction of intensity to out_a; clamped to [0, alpha].
                self.r = value.clamp(0.0, self.alpha);
                true
            }
            _ => false,
        }
    }

    /// One input, two outputs: both copies carry the input's label.
    fn lambda_routing(&self) -> Vec<(usize, usize)> {
        let (wpc, n) = (self.wpc, self.n_channels);
        let lam = wpc - 1;
        (0..n)
            .flat_map(|k| {
                [
                    (wpc * k + lam, wpc * n + wpc * k + lam),
                    (wpc * k + lam, 2 * wpc * n + wpc * k + lam),
                ]
            })
            .collect()
    }

    fn eval(&mut self, _x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {}

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let bpc = if wpc == 5 { 8 } else { 6 };
        let lam = wpc - 1;
        let port_c = wpc * n;
        let port_d = 2 * wpc * n;
        let k_a = self.r.max(0.0).sqrt();
        let k_b = (self.alpha - self.r).max(0.0).sqrt();
        for k in 0..n {
            let a_re_fw = self.nodes[wpc * k];
            let a_im_fw = self.nodes[wpc * k + 1];
            let a_l = self.nodes[wpc * k + lam];
            let c_re_fw = self.nodes[port_c + wpc * k];
            let c_im_fw = self.nodes[port_c + wpc * k + 1];
            let c_l = self.nodes[port_c + wpc * k + lam];
            let d_re_fw = self.nodes[port_d + wpc * k];
            let d_im_fw = self.nodes[port_d + wpc * k + 1];
            let d_l = self.nodes[port_d + wpc * k + lam];
            // Forward: c, d are scaled outputs of a.
            stamp_potential_eq(mat, &self.branches, bpc * k, c_re_fw, &[(a_re_fw, -k_a)]);
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k + 1,
                c_im_fw,
                &[(a_im_fw, -k_a)],
            );
            stamp_potential_eq(mat, &self.branches, bpc * k + 2, c_l, &[(a_l, -1.0)]);
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k + 3,
                d_re_fw,
                &[(a_re_fw, -k_b)],
            );
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k + 4,
                d_im_fw,
                &[(a_im_fw, -k_b)],
            );
            stamp_potential_eq(mat, &self.branches, bpc * k + 5, d_l, &[(a_l, -1.0)]);
            if wpc == 5 {
                // Backward: in.bw is the (reciprocal) combination of out_a.bw and out_b.bw.
                // The same coupling matrix applies on the way back:
                //   a.re_bw = k_a · c.re_bw + k_b · d.re_bw
                //   a.im_bw = k_a · c.im_bw + k_b · d.im_bw
                let a_re_bw = self.nodes[wpc * k + 2];
                let a_im_bw = self.nodes[wpc * k + 3];
                let c_re_bw = self.nodes[port_c + wpc * k + 2];
                let c_im_bw = self.nodes[port_c + wpc * k + 3];
                let d_re_bw = self.nodes[port_d + wpc * k + 2];
                let d_im_bw = self.nodes[port_d + wpc * k + 3];
                stamp_potential_eq(
                    mat,
                    &self.branches,
                    bpc * k + 6,
                    a_re_bw,
                    &[(c_re_bw, -k_a), (d_re_bw, -k_b)],
                );
                stamp_potential_eq(
                    mat,
                    &self.branches,
                    bpc * k + 7,
                    a_im_bw,
                    &[(c_im_bw, -k_a), (d_im_bw, -k_b)],
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
