use crate::device::{Device, EvalFlags, NodeId, ReactiveBranchSpec, ReactiveKind, SimContext};
use crate::mna::MnaMatrix;
use super::stamp_potential_eq;

// ────────────────────────────────────────────────────────────────────────
// Native directional coupler (2×2)
// ────────────────────────────────────────────────────────────────────────

/// 2×2 directional coupler with length-coupled cross-coefficient.
///
/// Lossless coupling matrix:
///   [c]   [ t   k] [a]    with t = cos(κL), k = -j sin(κL),
///   [d] = [ k   t] [b]    so |t|² + |k|² = 1.
///
/// In SVEA real/imag form:
///   c_re = t·a_re + s·b_im     d_re = t·b_re + s·a_im
///   c_im = t·a_im − s·b_re     d_im = t·b_im − s·a_re
/// (with t = cos(κL), s = sin(κL)).  Wavelength passes through unchanged
/// to both outputs from the corresponding input.
/// Variable-arity bundle-aware directional coupler.  Terminal layout for N
/// channels (12·N total terminals):
///   [a.0.re, a.0.im, a.0.λ, ..., a.{N-1}.λ,
///    b.0.re, b.0.im, b.0.λ, ..., b.{N-1}.λ,
///    c.0.re, c.0.im, c.0.λ, ..., c.{N-1}.λ,
///    d.0.re, d.0.im, d.0.λ, ..., d.{N-1}.λ]
/// Each channel uses the same κ·L (wavelength-independent at this tier).
/// The wavelength tag on each output channel mirrors port a's λ wire — for
/// closed-loop topologies that prevents a missing-driver bind-loop on d_λ.
pub struct NativeDirectionalCoupler {
    kappa_per_m: f64,
    length_m:    f64,
    n_channels:  usize,
    wpc:         usize,            // 3 or 5
    nodes:       Vec<NodeId>,
    branches:    Vec<Option<usize>>,
}

impl NativeDirectionalCoupler {
    /// Defaults model a circular MRR add-drop coupling region:
    ///   ring radius R = 8 µm, minimum bus-to-ring gap g₀ = 300 nm,
    ///   rib waveguide (500 nm × 220 nm core on 90 nm slab, SOI).
    /// Femwell supermode sweep + integration along the curved approach
    /// gives κ·L_total = 0.0769 rad (sin²(κL) ≈ 0.59 % power cross-coupling).
    /// `length_m` is an effective coupling length (where κ has dropped to
    /// 1 % of peak); the physical observable is `kappa_per_m · length_m =
    /// κ·L = 0.0769`.
    /// See `scripts/waveguide_simulations/coupler_sim.py` for the derivation.
    pub fn new() -> Self {
        let kappa_l = 0.0769;
        let length_m = 6.4e-6;
        Self {
            kappa_per_m: kappa_l / length_m,
            length_m,
            n_channels:  0,
            wpc:         3,
            nodes:       Vec::new(),
            branches:    Vec::new(),
        }
    }
}

impl Device for NativeDirectionalCoupler {
    fn num_terminals(&self) -> usize { self.nodes.len() }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 4 * wpc; // 4 ports per channel
        assert!(
            !terminals.is_empty() && terminals.len() % stride == 0,
            "fc_dcoupler: terminal count must be {stride}·N (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = terminals.len() / stride;
        self.n_channels = n;
        self.nodes      = terminals.to_vec();
        // Per channel: 4 fw branches (re, im for c and d) + (if bidir) 4 bw
        // branches (re, im for a and b) + 2 λ branches (c_λ, d_λ).
        let bpc = if wpc == 5 { 10 } else { 6 };
        self.branches   = vec![None; bpc * n];
    }

    fn num_extra_nodes(&self) -> usize { self.branches.len() }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..self.branches.len() { self.branches[i] = Some(first_idx + i); }
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "kappa_per_m" | "kappa" => { self.kappa_per_m = value; true }
            "l_um"   => { self.length_m = value * 1e-6; true }
            "l_m" | "length" => { self.length_m = value; true }
            "kappa_l" | "kappal" => {
                self.kappa_per_m = if self.length_m > 0.0 { value / self.length_m } else { 0.0 };
                true
            }
            _ => false,
        }
    }

    fn eval(&mut self, _x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {}

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        let kl  = self.kappa_per_m * self.length_m;
        let t   = kl.cos();
        let s   = kl.sin();
        let bpc = if wpc == 5 { 10 } else { 6 };
        let lam = wpc - 1;
        let port_b = wpc * n;
        let port_c = 2 * wpc * n;
        let port_d = 3 * wpc * n;
        for k in 0..n {
            let a_re_fw = self.nodes[wpc * k];
            let a_im_fw = self.nodes[wpc * k + 1];
            let a_l     = self.nodes[wpc * k + lam];
            let b_re_fw = self.nodes[port_b + wpc * k];
            let b_im_fw = self.nodes[port_b + wpc * k + 1];
            let c_re_fw = self.nodes[port_c + wpc * k];
            let c_im_fw = self.nodes[port_c + wpc * k + 1];
            let c_l     = self.nodes[port_c + wpc * k + lam];
            let d_re_fw = self.nodes[port_d + wpc * k];
            let d_im_fw = self.nodes[port_d + wpc * k + 1];
            let d_l     = self.nodes[port_d + wpc * k + lam];
            // Forward: c, d are outputs computed from a, b inputs.
            stamp_potential_eq(mat, &self.branches, bpc * k,     c_re_fw,
                &[(a_re_fw, -t), (b_im_fw, -s)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + 1, c_im_fw,
                &[(a_im_fw, -t), (b_re_fw,  s)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + 2, d_re_fw,
                &[(b_re_fw, -t), (a_im_fw, -s)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + 3, d_im_fw,
                &[(b_im_fw, -t), (a_re_fw,  s)]);
            // λ wires: c_λ = a_λ, d_λ = a_λ (closed-loop bind safety).
            stamp_potential_eq(mat, &self.branches, bpc * k + 4, c_l,
                &[(a_l, -1.0)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + 5, d_l,
                &[(a_l, -1.0)]);
            if wpc == 5 {
                // Bw: a, b are outputs computed from c, d inputs (same matrix,
                // reciprocal coupler).
                let a_re_bw = self.nodes[wpc * k + 2];
                let a_im_bw = self.nodes[wpc * k + 3];
                let b_re_bw = self.nodes[port_b + wpc * k + 2];
                let b_im_bw = self.nodes[port_b + wpc * k + 3];
                let c_re_bw = self.nodes[port_c + wpc * k + 2];
                let c_im_bw = self.nodes[port_c + wpc * k + 3];
                let d_re_bw = self.nodes[port_d + wpc * k + 2];
                let d_im_bw = self.nodes[port_d + wpc * k + 3];
                stamp_potential_eq(mat, &self.branches, bpc * k + 6, a_re_bw,
                    &[(c_re_bw, -t), (d_im_bw, -s)]);
                stamp_potential_eq(mat, &self.branches, bpc * k + 7, a_im_bw,
                    &[(c_im_bw, -t), (d_re_bw,  s)]);
                stamp_potential_eq(mat, &self.branches, bpc * k + 8, b_re_bw,
                    &[(d_re_bw, -t), (c_im_bw, -s)]);
                stamp_potential_eq(mat, &self.branches, bpc * k + 9, b_im_bw,
                    &[(d_im_bw, -t), (c_re_bw,  s)]);
            }
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) { self.load_residual(b); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) { self.load_jacobian(mat); }
}

