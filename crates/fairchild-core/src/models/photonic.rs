//! Native Rust photonic passive devices (B3).
//!
//! Each device implements `Device` directly — no Verilog-A round-trip, no
//! OSDI shared library, no Norton hack.  Outputs are stamped via the
//! direct-potential pattern (one auxiliary MNA row per output potential
//! equation, requested through `Device::num_extra_nodes` and bound via
//! `bind_extra_nodes`).
//!
//! Port convention (matches B1 discipline scheme, single-channel
//! forward-only): each optical port is a 3-wire bundle (re, im, λ).  The
//! .optical_port bundle directive (B2) maps a user-visible port name to
//! the underlying wires.
//!
//! All models in this module are equivalent-physics implementations from
//! public textbook formulas (Saleh-Teich, Heuck-Englund, Pozar) and carry
//! no PDK-specific calibration.

use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

// ────────────────────────────────────────────────────────────────────────
// Native straight waveguide
// ────────────────────────────────────────────────────────────────────────

/// Straight optical waveguide — propagation loss + accumulated phase.
///
/// Physics: `A_out = A_in · exp(-α·L/2) · exp(-j·β·L)` with `β = 2π·n_g/λ`.
///
/// Variable-arity bundle-aware device.  Terminal layout for N channels:
///   [in.0.re, in.0.im, in.0.λ,  ..., in.{N-1}.λ,
///    out.0.re, out.0.im, out.0.λ, ..., out.{N-1}.λ]   (6·N terminals)
/// Each channel runs independent re/im/λ propagation using its own input
/// wavelength wire.  No per-channel state is shared — this is a pure-optical
/// device — but having one instance for the whole bundle keeps WDM the rule
/// rather than the exception and simplifies stamping when channel count grows.
pub struct NativeWaveguide {
    length_m:        f64,
    n_g:             f64,
    alpha_neper_m:   f64,
    // Bootstrap λ for the first NR iterate (x = 0).  Sourced from
    // `SimContext::lambda_center_m` in `setup_model`.
    lambda_bootstrap_m: f64,
    n_channels:      usize,
    wpc:             usize,        // wires_per_channel: 3 (unidir) or 5 (bidir)
    nodes:    Vec<NodeId>,
    branches: Vec<Option<usize>>,  // wpc per channel
    c_cached: Vec<f64>,
    s_cached: Vec<f64>,
}

impl NativeWaveguide {
    pub fn new() -> Self {
        NativeWaveguide {
            length_m:           100e-6,
            n_g:                4.2,
            alpha_neper_m:      dB_per_cm_to_neper_per_m(2.0),
            lambda_bootstrap_m: 1.55e-6,
            n_channels:         0,
            wpc:                3,
            nodes:              Vec::new(),
            branches:           Vec::new(),
            c_cached:           Vec::new(),
            s_cached:           Vec::new(),
        }
    }
}

impl Device for NativeWaveguide {
    fn num_terminals(&self) -> usize { self.nodes.len() }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.lambda_bootstrap_m = ctx.lambda_center_m;
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 2 * wpc; // in + out
        assert!(
            !terminals.is_empty() && terminals.len() % stride == 0,
            "fc_waveguide: terminal count must be {stride}·N for N ≥ 1 channels (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = terminals.len() / stride;
        self.n_channels = n;
        self.nodes      = terminals.to_vec();
        self.branches   = vec![None; wpc * n];
        self.c_cached   = vec![1.0; n];
        self.s_cached   = vec![0.0; n];
    }

    fn num_extra_nodes(&self) -> usize { self.branches.len() }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..self.branches.len() { self.branches[i] = Some(first_idx + i); }
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "l_um"          => { self.length_m       = value * 1e-6;                       true }
            "l_m" | "length"=> { self.length_m       = value;                              true }
            "n_g"           => { self.n_g            = value;                              true }
            "alpha_db_cm"   => { self.alpha_neper_m  = dB_per_cm_to_neper_per_m(value);    true }
            _ => false,
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let boot = self.lambda_bootstrap_m;
        let t_amp = (-self.alpha_neper_m * self.length_m / 2.0).exp();
        let two_pi = 2.0 * std::f64::consts::PI;
        // λ wire position within each input channel block: last wire.
        let lambda_off = wpc - 1;
        for k in 0..n {
            let lambda = match self.nodes[wpc * k + lambda_off] {
                Some(i) => {
                    let v = x[i];
                    if v.abs() > boot * 0.5 { v } else { boot }
                }
                None => boot,
            };
            let phi = two_pi * self.n_g * self.length_m / lambda;
            self.c_cached[k] = t_amp * phi.cos();
            self.s_cached[k] = t_amp * phi.sin();
        }
    }

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        let in_block_base  = 0;
        let out_block_base = wpc * n;
        for k in 0..n {
            let c = self.c_cached[k];
            let s = self.s_cached[k];
            // Per-channel input / output wires.
            let in_re_fw  = self.nodes[in_block_base  + wpc * k];
            let in_im_fw  = self.nodes[in_block_base  + wpc * k + 1];
            let in_l      = self.nodes[in_block_base  + wpc * k + (wpc - 1)];
            let out_re_fw = self.nodes[out_block_base + wpc * k];
            let out_im_fw = self.nodes[out_block_base + wpc * k + 1];
            let out_l     = self.nodes[out_block_base + wpc * k + (wpc - 1)];
            // Forward path: out.re_fw = c·in.re_fw + s·in.im_fw etc.
            stamp_potential_eq(mat, &self.branches, wpc * k,     out_re_fw,
                &[(in_re_fw, -c), (in_im_fw, -s)]);
            stamp_potential_eq(mat, &self.branches, wpc * k + 1, out_im_fw,
                &[(in_re_fw,  s), (in_im_fw, -c)]);
            // λ passes through unchanged.
            stamp_potential_eq(mat, &self.branches, wpc * k + (wpc - 1), out_l,
                &[(in_l, -1.0)]);
            if wpc == 5 {
                // Backward path mirrors the forward physics with reversed
                // direction: light at out_bw propagates back to in_bw with
                // the same c, s.  in.re_bw = c·out.re_bw + s·out.im_bw,
                // in.im_bw = -s·out.re_bw + c·out.im_bw.
                let in_re_bw  = self.nodes[in_block_base  + wpc * k + 2];
                let in_im_bw  = self.nodes[in_block_base  + wpc * k + 3];
                let out_re_bw = self.nodes[out_block_base + wpc * k + 2];
                let out_im_bw = self.nodes[out_block_base + wpc * k + 3];
                stamp_potential_eq(mat, &self.branches, wpc * k + 2, in_re_bw,
                    &[(out_re_bw, -c), (out_im_bw, -s)]);
                stamp_potential_eq(mat, &self.branches, wpc * k + 3, in_im_bw,
                    &[(out_re_bw,  s), (out_im_bw, -c)]);
            }
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) { self.load_residual(b); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) { self.load_jacobian(mat); }
}

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
    pub fn new() -> Self {
        Self {
            kappa_per_m: 100.0,
            length_m:    5e-3,
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
    alpha:      f64,
    r:          f64,
    n_channels: usize,
    wpc:        usize,
    nodes:      Vec<NodeId>,
    branches:   Vec<Option<usize>>,
}

impl NativeSplitter {
    pub fn new() -> Self {
        Self {
            alpha:      1.0,
            r:          0.5,
            n_channels: 0,
            wpc:        3,
            nodes:      Vec::new(),
            branches:   Vec::new(),
        }
    }
}

impl Device for NativeSplitter {
    fn num_terminals(&self) -> usize { self.nodes.len() }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 3 * wpc; // 3 ports per channel
        assert!(
            !terminals.is_empty() && terminals.len() % stride == 0,
            "fc_splitter: terminal count must be {stride}·N (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = terminals.len() / stride;
        self.n_channels = n;
        self.nodes      = terminals.to_vec();
        // Per channel: 4 fw branches (re, im for out_a + out_b) + 2 λ +
        // (if bidir) 2 bw branches (re, im for in_a).  Note: under bidir the
        // splitter behaves like a combiner in reverse — bw light from out_a
        // and out_b combine back into in.  re_bw_in = k_a · re_bw_a + k_b · re_bw_b.
        let bpc = if wpc == 5 { 8 } else { 6 };
        self.branches   = vec![None; bpc * n];
    }

    fn num_extra_nodes(&self) -> usize { self.branches.len() }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..self.branches.len() { self.branches[i] = Some(first_idx + i); }
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "alpha" => {
                // Intensity transmission, must be ≤ 1.  Re-anchor r at the
                // symmetric midpoint (alpha/2) iff the user hasn't already
                // skewed it — but simplest: if r exceeds the new alpha,
                // clamp it to alpha/2.  Otherwise leave r alone.
                self.alpha = value.clamp(0.0, 1.0);
                if self.r > self.alpha { self.r = self.alpha * 0.5; }
                true
            }
            "alpha_db" | "il_db" => {
                self.alpha = 10f64.powf(-value / 10.0);
                if self.r > self.alpha { self.r = self.alpha * 0.5; }
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

    fn eval(&mut self, _x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {}

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n   = self.n_channels;
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
            let a_l     = self.nodes[wpc * k + lam];
            let c_re_fw = self.nodes[port_c + wpc * k];
            let c_im_fw = self.nodes[port_c + wpc * k + 1];
            let c_l     = self.nodes[port_c + wpc * k + lam];
            let d_re_fw = self.nodes[port_d + wpc * k];
            let d_im_fw = self.nodes[port_d + wpc * k + 1];
            let d_l     = self.nodes[port_d + wpc * k + lam];
            // Forward: c, d are scaled outputs of a.
            stamp_potential_eq(mat, &self.branches, bpc * k,     c_re_fw, &[(a_re_fw, -k_a)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + 1, c_im_fw, &[(a_im_fw, -k_a)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + 2, c_l,     &[(a_l, -1.0)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + 3, d_re_fw, &[(a_re_fw, -k_b)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + 4, d_im_fw, &[(a_im_fw, -k_b)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + 5, d_l,     &[(a_l, -1.0)]);
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
                stamp_potential_eq(mat, &self.branches, bpc * k + 6, a_re_bw,
                    &[(c_re_bw, -k_a), (d_re_bw, -k_b)]);
                stamp_potential_eq(mat, &self.branches, bpc * k + 7, a_im_bw,
                    &[(c_im_bw, -k_a), (d_im_bw, -k_b)]);
            }
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) { self.load_residual(b); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) { self.load_jacobian(mat); }
}

// ────────────────────────────────────────────────────────────────────────
// Native grating coupler — flat insertion loss, zero length
// ────────────────────────────────────────────────────────────────────────

/// Grating coupler (fibre ↔ chip).  Zero-length waveguide with a flat
/// amplitude attenuation set by `alpha_db` (insertion loss).  Variable-arity
/// bundle-aware: 6·N terminals.
pub struct NativeGratingCoupler {
    alpha_db:   f64,
    n_channels: usize,
    wpc:        usize,
    nodes:      Vec<NodeId>,
    branches:   Vec<Option<usize>>,
}

impl NativeGratingCoupler {
    pub fn new() -> Self {
        Self { alpha_db: 3.0, n_channels: 0, wpc: 3, nodes: Vec::new(), branches: Vec::new() }
    }
}

impl Device for NativeGratingCoupler {
    fn num_terminals(&self) -> usize { self.nodes.len() }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 2 * wpc;
        assert!(
            !terminals.is_empty() && terminals.len() % stride == 0,
            "fc_grating_coupler: terminal count must be {stride}·N (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = terminals.len() / stride;
        self.n_channels = n;
        self.nodes      = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 };
        self.branches   = vec![None; bpc * n];
    }

    fn num_extra_nodes(&self) -> usize { self.branches.len() }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..self.branches.len() { self.branches[i] = Some(first_idx + i); }
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "alpha_db" | "alpha_db_il" | "il_db" => { self.alpha_db = value; true }
            "alpha" => {
                let t = value.max(1e-30);
                self.alpha_db = -20.0 * t.log10();
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
        let bpc = if wpc == 5 { 5 } else { 3 };
        let lam = wpc - 1;
        let out_base = wpc * n;
        let t = 10f64.powf(-self.alpha_db / 20.0);
        for k in 0..n {
            let in_re_fw  = self.nodes[wpc * k];
            let in_im_fw  = self.nodes[wpc * k + 1];
            let in_l      = self.nodes[wpc * k + lam];
            let out_re_fw = self.nodes[out_base + wpc * k];
            let out_im_fw = self.nodes[out_base + wpc * k + 1];
            let out_l     = self.nodes[out_base + wpc * k + lam];
            stamp_potential_eq(mat, &self.branches, bpc * k,     out_re_fw, &[(in_re_fw, -t)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + 1, out_im_fw, &[(in_im_fw, -t)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + (bpc - 1), out_l, &[(in_l, -1.0)]);
            if wpc == 5 {
                let in_re_bw  = self.nodes[wpc * k + 2];
                let in_im_bw  = self.nodes[wpc * k + 3];
                let out_re_bw = self.nodes[out_base + wpc * k + 2];
                let out_im_bw = self.nodes[out_base + wpc * k + 3];
                stamp_potential_eq(mat, &self.branches, bpc * k + 2, in_re_bw, &[(out_re_bw, -t)]);
                stamp_potential_eq(mat, &self.branches, bpc * k + 3, in_im_bw, &[(out_im_bw, -t)]);
            }
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) { self.load_residual(b); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) { self.load_jacobian(mat); }
}

// ────────────────────────────────────────────────────────────────────────
// Native CW laser source
// ────────────────────────────────────────────────────────────────────────

/// Constant-amplitude SVEA source.  Drives the three output wires of a
/// single optical-port bundle to a fixed (re, im, λ) value via direct
/// potential contributions — no electrical input.
///
/// `A_re = √P · cos(φ₀)`, `A_im = √P · sin(φ₀)` where `P = power_mW · 1e−3`.
/// CW laser source.  Drives ONE optical channel's bundle wires: 3 wires
/// (re, im, λ) under unidirectional propagation, or 5 wires (re_fw, im_fw,
/// re_bw, im_bw, λ) under bidirectional — the bw wires are forced to 0
/// because a laser emits in one direction only.
pub struct NativeCwLaser {
    re_amp:     f64,
    im_amp:     f64,
    wavelen_m:  f64,
    wpc:        usize,             // 3 (unidir) or 5 (bidir)
    nodes:      Vec<NodeId>,
    branches:   Vec<Option<usize>>,
}

impl NativeCwLaser {
    pub fn new() -> Self {
        // Defaults: 1 mW, 0° phase, 1550 nm.
        let p = 1e-3_f64;
        Self {
            re_amp:    p.sqrt(),
            im_amp:    0.0,
            wavelen_m: 1550e-9,
            wpc:       3,
            nodes:     Vec::new(),
            branches:  Vec::new(),
        }
    }
}

impl Device for NativeCwLaser {
    fn num_terminals(&self) -> usize { self.nodes.len() }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wavelen_m = ctx.lambda_center_m;
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        debug_assert_eq!(terminals.len(), wpc,
            "fc_cw_laser: expected {wpc} terminals (one channel × wpc); got {}",
            terminals.len());
        self.nodes    = terminals.to_vec();
        self.branches = vec![None; wpc];
    }

    fn num_extra_nodes(&self) -> usize { self.branches.len() }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..self.branches.len() { self.branches[i] = Some(first_idx + i); }
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "power_mw" => {
                let p = (value * 1e-3).max(0.0);
                let phi = (self.im_amp / self.re_amp.max(1e-30)).atan();
                let mag = p.sqrt();
                self.re_amp = mag * phi.cos();
                self.im_amp = mag * phi.sin();
                true
            }
            "power_w" => {
                let p = value.max(0.0);
                let phi = (self.im_amp / self.re_amp.max(1e-30)).atan();
                let mag = p.sqrt();
                self.re_amp = mag * phi.cos();
                self.im_amp = mag * phi.sin();
                true
            }
            "phi_0_deg" | "phase_deg" => {
                let mag = (self.re_amp * self.re_amp + self.im_amp * self.im_amp).sqrt();
                let phi = value * std::f64::consts::PI / 180.0;
                self.re_amp = mag * phi.cos();
                self.im_amp = mag * phi.sin();
                true
            }
            "wavelength_nm" => { self.wavelen_m = value * 1e-9; true }
            "wavelength_m"  => { self.wavelen_m = value; true }
            "re_amp" => { self.re_amp = value; true }
            "im_amp" => { self.im_amp = value; true }
            _ => false,
        }
    }

    fn eval(&mut self, _x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {}

    fn load_residual(&self, b: &mut [f64]) {
        // Inhomogeneous branch equations: V(out_re_fw) = re_amp, ...
        // Wire order: [re_fw, im_fw, re_bw, im_bw, λ] (5-wire bidir) or
        //             [re,    im,    λ]               (3-wire unidir).
        if self.wpc == 5 {
            if let Some(j) = self.branches[0] { b[j] += self.re_amp; }
            if let Some(j) = self.branches[1] { b[j] += self.im_amp; }
            // bw wires forced to 0 — no contribution from RHS (branch row
            // already enforces V = 0 because rhs is 0).
            if let Some(j) = self.branches[4] { b[j] += self.wavelen_m; }
        } else {
            if let Some(j) = self.branches[0] { b[j] += self.re_amp; }
            if let Some(j) = self.branches[1] { b[j] += self.im_amp; }
            if let Some(j) = self.branches[2] { b[j] += self.wavelen_m; }
        }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        // Branch rows: V(out_wire) - target = 0.  Stamp +1 at (J, out) and
        // +1 at (out, J).  RHS handled in load_residual.  bw wires use the
        // same shape, with target = 0 (no explicit RHS contribution).
        for (i, out_node) in self.nodes.iter().enumerate() {
            if let (Some(out), Some(j)) = (*out_node, self.branches[i]) {
                mat.a[j][out] += 1.0;
                mat.a[out][j] += 1.0;
            }
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) { self.load_residual(b); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) { self.load_jacobian(mat); }
}

// ────────────────────────────────────────────────────────────────────────
// Native photodetector (PIN — instantaneous responsivity)
// ────────────────────────────────────────────────────────────────────────

/// PIN photodetector with linear responsivity and a shunt resistance.
///
/// Physics:  I_ph = R · (re² + im²) + I_dark
/// Flows cathode → anode (reverse-biased junction convention).  A shunt
/// resistance models the junction impedance.
///
/// The current is a nonlinear function of the optical inputs, so this
/// device contributes both a residual (I_ph) and the linearised Jacobian
/// terms ∂I_ph/∂V(in_re), ∂I_ph/∂V(in_im) and 1/R_shunt for the V/R
/// shunt.  No internal nodes are needed — the photocurrent stamps
/// directly between the electrical terminals.
/// Variable-arity to support WDM: ONE physical photodetector with one shared
/// dark current and shunt, summing photocurrents across N optical channels.
/// Terminal layout (3·N + 2 terminals for N channels):
///
///   [in.0.re, in.0.im, in.0.λ,  in.1.re, in.1.im, in.1.λ,  ...,  in.{N-1}.λ,
///    anode, cathode]
///
/// The parser does NOT replicate this device (see fairchild-parser's
/// `expand_optical_ports` exception list) — instead the device flattens all
/// channels into one terminal block.  Photocurrent
/// `I_ph = responsivity · Σ_k (re_k² + im_k²) + i_dark` is computed in one
/// place; the shunt `1/r_shunt` is stamped once between anode and cathode.
pub struct NativePhotodetector {
    responsivity:  f64,
    i_dark:        f64,
    r_shunt:       f64,
    r_series:      f64,
    n_channels:    usize,
    wpc:           usize,            // 3 (unidir) or 5 (bidir)
    nodes:         Vec<NodeId>,
    v_int_idx:     Option<usize>,
    has_internal:  bool,
    i_ph:          f64,
    // Linearisation coefficients per (channel, direction).  In unidir mode
    // only g[k].0/.1 are used; in bidir mode the .2/.3 entries cover the bw
    // wires.  A real PIN absorbs every photon — both fw and bw light heat
    // the same junction and produce one summed photocurrent.
    g_re_fw:       Vec<f64>,
    g_im_fw:       Vec<f64>,
    g_re_bw:       Vec<f64>,
    g_im_bw:       Vec<f64>,
    v_re_fw_op:    Vec<f64>,
    v_im_fw_op:    Vec<f64>,
    v_re_bw_op:    Vec<f64>,
    v_im_bw_op:    Vec<f64>,
    v_j_op:        f64,
}

impl NativePhotodetector {
    pub fn new() -> Self {
        Self {
            responsivity: 1.0,
            i_dark:       1e-9,
            r_shunt:      1e6,
            r_series:     0.0,
            n_channels:   0,
            wpc:          3,
            nodes:        Vec::new(),
            v_int_idx:    None,
            has_internal: false,
            i_ph:         0.0,
            g_re_fw:      Vec::new(),
            g_im_fw:      Vec::new(),
            g_re_bw:      Vec::new(),
            g_im_bw:      Vec::new(),
            v_re_fw_op:   Vec::new(),
            v_im_fw_op:   Vec::new(),
            v_re_bw_op:   Vec::new(),
            v_im_bw_op:   Vec::new(),
            v_j_op:       0.0,
        }
    }
}

impl Device for NativePhotodetector {
    fn num_terminals(&self) -> usize { self.nodes.len() }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        // Layout: wpc·N (bundle inputs) + 2 (anode, cathode).
        assert!(
            terminals.len() >= wpc + 2 && (terminals.len() - 2) % wpc == 0,
            "fc_photodetector: terminal count must be {wpc}·N + 2 for N ≥ 1 channels; got {}",
            terminals.len()
        );
        let n = (terminals.len() - 2) / wpc;
        self.n_channels = n;
        self.nodes      = terminals.to_vec();
        self.g_re_fw    = vec![0.0; n];
        self.g_im_fw    = vec![0.0; n];
        self.g_re_bw    = vec![0.0; n];
        self.g_im_bw    = vec![0.0; n];
        self.v_re_fw_op = vec![0.0; n];
        self.v_im_fw_op = vec![0.0; n];
        self.v_re_bw_op = vec![0.0; n];
        self.v_im_bw_op = vec![0.0; n];
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "responsivity" => { self.responsivity = value; true }
            "i_dark" | "i_dark_a" => { self.i_dark = value; true }
            "r_shunt" => { self.r_shunt = value; true }
            "r_series" | "r_s" => { self.r_series = value.max(0.0); true }
            // `c_par` (and a bias-dependent C_j(V)) lands with the L2 PD
            // model — it needs device-internal reactive state which the
            // current Device trait doesn't expose.  Accept the keyword to
            // keep schematics forward-compatible; document as no-op.
            "c_par" | "c_j0" | "c_par_f" => true,
            _ => false,
        }
    }

    fn num_extra_nodes(&self) -> usize {
        // Allocate one internal MNA row (the "junction node" between r_series
        // and the parallel R_sh || I_ph stack) when r_series is non-zero.
        if self.r_series > 0.0 { 1 } else { 0 }
    }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        if self.r_series > 0.0 {
            self.v_int_idx    = Some(first_idx);
            self.has_internal = true;
        } else {
            self.v_int_idx    = None;
            self.has_internal = false;
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        let elec_base = wpc * n;
        let v_j_node = match self.v_int_idx {
            Some(i) => x[i],
            None    => self.nodes[elec_base].map_or(0.0, |i| x[i]),
        };
        let v_c = self.nodes[elec_base + 1].map_or(0.0, |i| x[i]);
        let mut p_total = 0.0;
        for k in 0..n {
            let v_re_fw = self.nodes[wpc * k    ].map_or(0.0, |i| x[i]);
            let v_im_fw = self.nodes[wpc * k + 1].map_or(0.0, |i| x[i]);
            p_total += v_re_fw * v_re_fw + v_im_fw * v_im_fw;
            self.g_re_fw[k] = 2.0 * self.responsivity * v_re_fw;
            self.g_im_fw[k] = 2.0 * self.responsivity * v_im_fw;
            self.v_re_fw_op[k] = v_re_fw;
            self.v_im_fw_op[k] = v_im_fw;
            if wpc == 5 {
                let v_re_bw = self.nodes[wpc * k + 2].map_or(0.0, |i| x[i]);
                let v_im_bw = self.nodes[wpc * k + 3].map_or(0.0, |i| x[i]);
                p_total += v_re_bw * v_re_bw + v_im_bw * v_im_bw;
                self.g_re_bw[k] = 2.0 * self.responsivity * v_re_bw;
                self.g_im_bw[k] = 2.0 * self.responsivity * v_im_bw;
                self.v_re_bw_op[k] = v_re_bw;
                self.v_im_bw_op[k] = v_im_bw;
            }
        }
        self.i_ph   = self.responsivity * p_total + self.i_dark;
        self.v_j_op = v_j_node - v_c;
    }

    fn load_residual(&self, b: &mut [f64]) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        let mut nonlin_remainder = self.i_ph;
        for k in 0..n {
            nonlin_remainder -= self.g_re_fw[k] * self.v_re_fw_op[k]
                              + self.g_im_fw[k] * self.v_im_fw_op[k];
            if wpc == 5 {
                nonlin_remainder -= self.g_re_bw[k] * self.v_re_bw_op[k]
                                  + self.g_im_bw[k] * self.v_im_bw_op[k];
            }
        }
        let i_eq = -nonlin_remainder - (self.v_j_op / self.r_shunt);
        let elec_base = wpc * n;
        let v_j_node = self.v_int_idx.or(self.nodes[elec_base]);
        if let Some(j) = v_j_node              { b[j] -= i_eq; }
        if let Some(c) = self.nodes[elec_base + 1] { b[c] += i_eq; }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        let elec_base = wpc * n;
        let a_idx = self.nodes[elec_base];
        let c_idx = self.nodes[elec_base + 1];
        let g_sh  = 1.0 / self.r_shunt;
        let j_idx = self.v_int_idx.or(a_idx);
        if let Some(j) = self.v_int_idx {
            let g_s = 1.0 / self.r_series;
            if let Some(a) = a_idx {
                mat.a[a][a] += g_s;
                mat.a[a][j] -= g_s;
                mat.a[j][a] -= g_s;
            }
            mat.a[j][j] += g_s;
        }
        if let Some(j) = j_idx {
            mat.a[j][j] += g_sh;
            if let Some(c) = c_idx { mat.a[j][c] -= g_sh; }
        }
        if let Some(c) = c_idx {
            mat.a[c][c] += g_sh;
            if let Some(j) = j_idx { mat.a[c][j] -= g_sh; }
        }
        for k in 0..n {
            let r_fw = self.nodes[wpc * k    ];
            let i_fw = self.nodes[wpc * k + 1];
            if let Some(j) = j_idx {
                if let Some(r) = r_fw { mat.a[j][r] -= self.g_re_fw[k]; }
                if let Some(i) = i_fw { mat.a[j][i] -= self.g_im_fw[k]; }
            }
            if let Some(c) = c_idx {
                if let Some(r) = r_fw { mat.a[c][r] += self.g_re_fw[k]; }
                if let Some(i) = i_fw { mat.a[c][i] += self.g_im_fw[k]; }
            }
            if wpc == 5 {
                let r_bw = self.nodes[wpc * k + 2];
                let i_bw = self.nodes[wpc * k + 3];
                if let Some(j) = j_idx {
                    if let Some(r) = r_bw { mat.a[j][r] -= self.g_re_bw[k]; }
                    if let Some(i) = i_bw { mat.a[j][i] -= self.g_im_bw[k]; }
                }
                if let Some(c) = c_idx {
                    if let Some(r) = r_bw { mat.a[c][r] += self.g_re_bw[k]; }
                    if let Some(i) = i_bw { mat.a[c][i] += self.g_im_bw[k]; }
                }
            }
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) { self.load_residual(b); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) { self.load_jacobian(mat); }
}

// ────────────────────────────────────────────────────────────────────────
// Native thermal phase shifter
// ────────────────────────────────────────────────────────────────────────

/// Thermal phase shifter (heater).
///
/// Electrical side: resistive heater with conductance `1/R_heater` between
/// anode and cathode.  Joule power `P = V²/R` is converted to an optical
/// phase shift `φ = π · P / P_pi`, where `P_pi` is the heater power
/// required for a π phase shift.
///
/// Optical side: 3-wire bundle in → 3-wire bundle out, applies `exp(-jφ)`.
/// Wavelength passes through unchanged.
///
/// 8 terminals: [in_re, in_im, in_λ, out_re, out_im, out_λ, anode, cathode]
/// 3 internal branch rows for the three direct-potential outputs.
/// Variable-arity: the parser does NOT replicate this device per channel.
/// One instance handles all N optical channels with one shared heater
/// resistor.  Terminal layout: [in.0.re,...,in.{N-1}.λ, out.0.re,...,
/// out.{N-1}.λ, heat_p, heat_n] = 6N + 2.
pub struct NativeThermalPhaseShifter {
    r_heater: f64,
    p_pi:     f64,
    n_channels: usize,
    wpc:      usize,
    nodes:    Vec<NodeId>,
    branches: Vec<Option<usize>>,
    c_cached: f64,
    s_cached: f64,
}

impl NativeThermalPhaseShifter {
    pub fn new() -> Self {
        Self {
            r_heater: 1000.0,
            p_pi:     10e-3,
            n_channels: 0,
            wpc:      3,
            nodes:    Vec::new(),
            branches: Vec::new(),
            c_cached: 1.0,
            s_cached: 0.0,
        }
    }
}

impl Device for NativeThermalPhaseShifter {
    fn num_terminals(&self) -> usize { self.nodes.len() }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 2 * wpc;
        assert!(
            terminals.len() >= stride + 2 && (terminals.len() - 2) % stride == 0,
            "fc_thermal_ps: terminal count must be {stride}·N + 2 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 2) / stride;
        self.n_channels = n;
        self.nodes      = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 };
        self.branches   = vec![None; bpc * n];
    }

    fn num_extra_nodes(&self) -> usize { self.branches.len() }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..self.branches.len() {
            self.branches[i] = Some(first_idx + i);
        }
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "r_heater" | "r" => { self.r_heater = value; true }
            "p_pi" | "p_pi_w" => { self.p_pi = value; true }
            _ => false,
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        let elec_base = 2 * wpc * n;
        let v_a = self.nodes[elec_base    ].map_or(0.0, |i| x[i]);
        let v_c = self.nodes[elec_base + 1].map_or(0.0, |i| x[i]);
        let v   = v_a - v_c;
        let p   = v * v / self.r_heater;
        let phi = std::f64::consts::PI * p / self.p_pi;
        self.c_cached = phi.cos();
        self.s_cached = phi.sin();
    }

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        let bpc = if wpc == 5 { 5 } else { 3 };
        let lam = wpc - 1;
        let out_base  = wpc * n;
        let elec_base = 2 * wpc * n;
        // Electrical: ONE shared heater resistor.
        let g = 1.0 / self.r_heater;
        let p = self.nodes[elec_base];
        let m = self.nodes[elec_base + 1];
        if let Some(a) = p {
            mat.a[a][a] += g;
            if let Some(c) = m { mat.a[a][c] -= g; }
        }
        if let Some(c) = m {
            mat.a[c][c] += g;
            if let Some(a) = p { mat.a[c][a] -= g; }
        }
        let c_cos = self.c_cached;
        let s_sin = self.s_cached;
        for k in 0..n {
            let in_re_fw  = self.nodes[wpc * k];
            let in_im_fw  = self.nodes[wpc * k + 1];
            let in_l      = self.nodes[wpc * k + lam];
            let out_re_fw = self.nodes[out_base + wpc * k];
            let out_im_fw = self.nodes[out_base + wpc * k + 1];
            let out_l     = self.nodes[out_base + wpc * k + lam];
            stamp_potential_eq(mat, &self.branches, bpc * k,     out_re_fw,
                &[(in_re_fw, -c_cos), (in_im_fw, -s_sin)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + 1, out_im_fw,
                &[(in_re_fw,  s_sin), (in_im_fw, -c_cos)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + (bpc - 1), out_l,
                &[(in_l, -1.0)]);
            if wpc == 5 {
                let in_re_bw  = self.nodes[wpc * k + 2];
                let in_im_bw  = self.nodes[wpc * k + 3];
                let out_re_bw = self.nodes[out_base + wpc * k + 2];
                let out_im_bw = self.nodes[out_base + wpc * k + 3];
                stamp_potential_eq(mat, &self.branches, bpc * k + 2, in_re_bw,
                    &[(out_re_bw, -c_cos), (out_im_bw, -s_sin)]);
                stamp_potential_eq(mat, &self.branches, bpc * k + 3, in_im_bw,
                    &[(out_re_bw,  s_sin), (out_im_bw, -c_cos)]);
            }
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) { self.load_residual(b); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) { self.load_jacobian(mat); }
}

// ────────────────────────────────────────────────────────────────────────
// Native PN-junction phase shifter
// ────────────────────────────────────────────────────────────────────────

/// PN-junction phase shifter (carrier-depletion or carrier-injection).
///
/// Electrical side: linearised PN junction.  Modelled as a parallel
/// combination of a small ohmic resistance (1 / G_pn) and a linear
/// capacitance (used in transient via the standard Norton C stamp; in DC
/// this contributes nothing).  Voltage dependence of the cap is ignored at
/// this first-pass level.
///
/// Optical side: phase shift `φ = 2π · L · Δn_eff / λ`, where the
/// effective-index change is linearised as `Δn_eff = δn_dV · V_pn`.  This
/// reproduces the small-signal behaviour of either depletion or carrier
/// injection modulators when calibrated to measurements (parameter
/// `dn_dv`).  Wavelength passes through unchanged.
/// Variable-arity to support WDM: the parser does NOT replicate this device
/// per channel (see `BUNDLE_AWARE_MODELS` in fairchild-parser).  One instance
/// handles all N optical channels with one shared PN junction.  Terminal
/// layout for N channels (6·N + 2 total terminals):
///
///   [in.0.re, in.0.im, in.0.λ,  in.1.re, in.1.im, in.1.λ,  ...,  in.{N-1}.λ,
///    out.0.re, out.0.im, out.0.λ,  ...,  out.{N-1}.λ,
///    anode, cathode]
///
/// The wavelength wires read independently per channel (a WDM laser can
/// drive each channel at a different λ), but the electrical conductance,
/// the EO Δn_eff, and the loss factor are all shared — the single physical
/// device sees one V_pn across one junction regardless of how many
/// wavelengths pass through.
pub struct NativePnPhaseShifter {
    length_m: f64,
    n_g:      f64,
    wl_ref_m: f64,
    dn_dv:    f64,
    g_pn:     f64,
    alpha_neper_m: f64,
    n_channels: usize,
    wpc:      usize,                 // 3 or 5
    nodes:    Vec<NodeId>,
    branches: Vec<Option<usize>>,
    c_cached: Vec<f64>,
    s_cached: Vec<f64>,
}

impl NativePnPhaseShifter {
    pub fn new() -> Self {
        Self {
            length_m: 1e-3,
            n_g:      4.2,
            wl_ref_m: 1.55e-6,
            dn_dv:    1e-4,
            g_pn:     1e-3,
            alpha_neper_m: 0.0,
            n_channels: 0,
            wpc:      3,
            nodes:    Vec::new(),
            branches: Vec::new(),
            c_cached: Vec::new(),
            s_cached: Vec::new(),
        }
    }
}

impl Device for NativePnPhaseShifter {
    fn num_terminals(&self) -> usize { self.nodes.len() }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wl_ref_m = ctx.lambda_center_m;
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 2 * wpc; // in + out bundle
        // Layout: wpc·N (in) + wpc·N (out) + 2 (anode, cathode).
        assert!(
            terminals.len() >= stride + 2 && (terminals.len() - 2) % stride == 0,
            "fc_pn_ps: terminal count must be {stride}·N + 2 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 2) / stride;
        self.n_channels = n;
        self.nodes      = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 }; // branches per channel
        self.branches   = vec![None; bpc * n];
        self.c_cached   = vec![1.0; n];
        self.s_cached   = vec![0.0; n];
    }

    fn num_extra_nodes(&self) -> usize { self.branches.len() }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..self.branches.len() {
            self.branches[i] = Some(first_idx + i);
        }
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "l_um"   => { self.length_m = value * 1e-6; true }
            "l_m" | "length" => { self.length_m = value; true }
            "n_g"    => { self.n_g = value; true }
            "dn_dv"  => { self.dn_dv = value; true }
            "g_pn"   => { self.g_pn  = value; true }
            "v_pi_l" => {
                // Vπ·L (V·m): solve for dn_dv such that the EO phase shift
                // is π at V = Vπ.  2π·L·dn_dv·Vπ/λ_ref = π →
                // dn_dv = λ_ref / (2·L·Vπ).
                if value > 0.0 {
                    self.dn_dv = self.wl_ref_m / (2.0 * value);
                }
                true
            }
            "alpha_db_cm" => {
                self.alpha_neper_m = dB_per_cm_to_neper_per_m(value);
                true
            }
            _ => false,
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        let elec_base = 2 * wpc * n;
        let v_a = self.nodes[elec_base    ].map_or(0.0, |i| x[i]);
        let v_c = self.nodes[elec_base + 1].map_or(0.0, |i| x[i]);
        let v_pn = v_a - v_c;
        let two_pi = 2.0 * std::f64::consts::PI;
        let t_amp = (-self.alpha_neper_m * self.length_m / 2.0).exp();
        let lam = wpc - 1;
        for k in 0..n {
            let lambda = match self.nodes[wpc * k + lam] {
                Some(i) => {
                    let v = x[i];
                    if v.abs() > 1e-9 { v } else { self.wl_ref_m }
                }
                None => self.wl_ref_m,
            };
            let phi_prop = two_pi * self.n_g * self.length_m
                           * (1.0 / lambda - 1.0 / self.wl_ref_m);
            let phi_eo   = two_pi * self.length_m * self.dn_dv * v_pn / lambda;
            let phi      = phi_prop + phi_eo;
            self.c_cached[k] = t_amp * phi.cos();
            self.s_cached[k] = t_amp * phi.sin();
        }
    }

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        let bpc = if wpc == 5 { 5 } else { 3 };
        let lam = wpc - 1;
        let out_base  = wpc * n;
        let elec_base = 2 * wpc * n;
        // Electrical: ONE shared PN-junction conductance.
        let g = self.g_pn;
        let anode = self.nodes[elec_base];
        let cath  = self.nodes[elec_base + 1];
        if let Some(a) = anode {
            mat.a[a][a] += g;
            if let Some(c) = cath { mat.a[a][c] -= g; }
        }
        if let Some(c) = cath {
            mat.a[c][c] += g;
            if let Some(a) = anode { mat.a[c][a] -= g; }
        }
        // Optical: per-channel branch equations.
        for k in 0..n {
            let in_re_fw  = self.nodes[wpc * k];
            let in_im_fw  = self.nodes[wpc * k + 1];
            let in_l      = self.nodes[wpc * k + lam];
            let out_re_fw = self.nodes[out_base + wpc * k];
            let out_im_fw = self.nodes[out_base + wpc * k + 1];
            let out_l     = self.nodes[out_base + wpc * k + lam];
            let c = self.c_cached[k];
            let s = self.s_cached[k];
            stamp_potential_eq(mat, &self.branches, bpc * k,     out_re_fw,
                &[(in_re_fw, -c), (in_im_fw, -s)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + 1, out_im_fw,
                &[(in_re_fw,  s), (in_im_fw, -c)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + (bpc - 1), out_l,
                &[(in_l, -1.0)]);
            if wpc == 5 {
                // Backward path mirrors fw: bw entering at out exits at in
                // with same propagation + EO phase shift (one physical
                // junction; same Δn applies to either direction).
                let in_re_bw  = self.nodes[wpc * k + 2];
                let in_im_bw  = self.nodes[wpc * k + 3];
                let out_re_bw = self.nodes[out_base + wpc * k + 2];
                let out_im_bw = self.nodes[out_base + wpc * k + 3];
                stamp_potential_eq(mat, &self.branches, bpc * k + 2, in_re_bw,
                    &[(out_re_bw, -c), (out_im_bw, -s)]);
                stamp_potential_eq(mat, &self.branches, bpc * k + 3, in_im_bw,
                    &[(out_re_bw,  s), (out_im_bw, -c)]);
            }
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) { self.load_residual(b); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) { self.load_jacobian(mat); }
}

// ────────────────────────────────────────────────────────────────────────
// Native combined PN + thermal phase shifter
// ────────────────────────────────────────────────────────────────────────

/// Waveguide segment with both a PN junction (electro-optic) AND a thermal
/// heater driving it.  Δn contributions sum.  Two electrical interfaces
/// (anode/cathode for the PN, heat_p/heat_n for the heater) are independent
/// — driving either alone produces only that physics's phase shift; driving
/// both produces the sum.
///
/// Variable-arity bundle-aware.  Terminal layout for N channels
/// (6·N + 4 total):
///   [in.0.re..in.{N-1}.λ,  out.0.re..out.{N-1}.λ,
///    anode, cathode, heat_p, heat_n]
///
/// At L1 the physics is the linear sum of fc_pn_ps (small-signal Δn_eff =
/// dn/dV · V_pn) and fc_thermal_ps (instantaneous Joule heating →
/// φ_th = π · P / P_pi_th).  L2 will add bias-dependent C_j, tau_th, and
/// distinct reverse/forward EO coefficients; L3 will add carrier dynamics
/// and self-heating from optical absorption.
pub struct NativePnThermalPhaseShifter {
    length_m:        f64,
    n_g:             f64,
    wl_ref_m:        f64,
    // PN-side params (small-signal electro-optic).
    dn_dv:           f64,
    g_pn:            f64,
    // Heater-side params (Joule → phase).
    r_heater:        f64,
    p_pi_th:         f64,
    // Shared loss along the segment.
    alpha_neper_m:   f64,
    n_channels:      usize,
    wpc:             usize,
    nodes:           Vec<NodeId>,
    branches:        Vec<Option<usize>>,
    c_cached:        Vec<f64>,
    s_cached:        Vec<f64>,
}

impl NativePnThermalPhaseShifter {
    pub fn new() -> Self {
        Self {
            length_m:        1e-3,
            n_g:             4.2,
            wl_ref_m:        1.55e-6,
            dn_dv:           1e-4,
            g_pn:            1e-3,
            r_heater:        1000.0,
            p_pi_th:         10e-3,
            alpha_neper_m:   0.0,
            n_channels:      0,
            wpc:             3,
            nodes:           Vec::new(),
            branches:        Vec::new(),
            c_cached:        Vec::new(),
            s_cached:        Vec::new(),
        }
    }
}

impl Device for NativePnThermalPhaseShifter {
    fn num_terminals(&self) -> usize { self.nodes.len() }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wl_ref_m = ctx.lambda_center_m;
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 2 * wpc;
        assert!(
            terminals.len() >= stride + 4 && (terminals.len() - 4) % stride == 0,
            "fc_pn_th_ps: terminal count must be {stride}·N + 4 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 4) / stride;
        self.n_channels = n;
        self.nodes      = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 };
        self.branches   = vec![None; bpc * n];
        self.c_cached   = vec![1.0; n];
        self.s_cached   = vec![0.0; n];
    }

    fn num_extra_nodes(&self) -> usize { self.branches.len() }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..self.branches.len() { self.branches[i] = Some(first_idx + i); }
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "l_um"            => { self.length_m = value * 1e-6; true }
            "l_m" | "length"  => { self.length_m = value;        true }
            "n_g"             => { self.n_g      = value;        true }
            "dn_dv"           => { self.dn_dv    = value;        true }
            "g_pn"            => { self.g_pn     = value;        true }
            "v_pi_l"          => {
                if value > 0.0 { self.dn_dv = self.wl_ref_m / (2.0 * value); }
                true
            }
            "r_heater" | "r"  => { self.r_heater = value;        true }
            "p_pi" | "p_pi_w" | "p_pi_th" => { self.p_pi_th = value; true }
            "alpha_db_cm"     => {
                self.alpha_neper_m = dB_per_cm_to_neper_per_m(value);
                true
            }
            _ => false,
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        let elec = 2 * wpc * n;
        let v_a    = self.nodes[elec    ].map_or(0.0, |i| x[i]);
        let v_c    = self.nodes[elec + 1].map_or(0.0, |i| x[i]);
        let v_hp   = self.nodes[elec + 2].map_or(0.0, |i| x[i]);
        let v_hn   = self.nodes[elec + 3].map_or(0.0, |i| x[i]);
        let v_pn   = v_a  - v_c;
        let v_h    = v_hp - v_hn;
        let p_heat = v_h * v_h / self.r_heater;
        let phi_th = std::f64::consts::PI * p_heat / self.p_pi_th;
        let two_pi = 2.0 * std::f64::consts::PI;
        let t_amp  = (-self.alpha_neper_m * self.length_m / 2.0).exp();
        let lam = wpc - 1;
        for k in 0..n {
            let lambda = match self.nodes[wpc * k + lam] {
                Some(i) => {
                    let v = x[i];
                    if v.abs() > 1e-9 { v } else { self.wl_ref_m }
                }
                None => self.wl_ref_m,
            };
            let phi_prop = two_pi * self.n_g * self.length_m
                           * (1.0 / lambda - 1.0 / self.wl_ref_m);
            let phi_eo   = two_pi * self.length_m * self.dn_dv * v_pn / lambda;
            let phi      = phi_prop + phi_eo + phi_th;
            self.c_cached[k] = t_amp * phi.cos();
            self.s_cached[k] = t_amp * phi.sin();
        }
    }

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        let bpc = if wpc == 5 { 5 } else { 3 };
        let lam = wpc - 1;
        let out_base = wpc * n;
        let elec = 2 * wpc * n;
        // Electrical (PN side): ONE g_pn between anode and cathode.
        let anode = self.nodes[elec];
        let cath  = self.nodes[elec + 1];
        let g_pn  = self.g_pn;
        if let Some(a) = anode {
            mat.a[a][a] += g_pn;
            if let Some(c) = cath { mat.a[a][c] -= g_pn; }
        }
        if let Some(c) = cath {
            mat.a[c][c] += g_pn;
            if let Some(a) = anode { mat.a[c][a] -= g_pn; }
        }
        // Electrical (heater side): ONE g_heater between heat_p and heat_n.
        let hp = self.nodes[elec + 2];
        let hn = self.nodes[elec + 3];
        let g_h = 1.0 / self.r_heater;
        if let Some(a) = hp {
            mat.a[a][a] += g_h;
            if let Some(c) = hn { mat.a[a][c] -= g_h; }
        }
        if let Some(c) = hn {
            mat.a[c][c] += g_h;
            if let Some(a) = hp { mat.a[c][a] -= g_h; }
        }
        // Optical: per-channel branch equations using cached rotation.
        for k in 0..n {
            let in_re_fw  = self.nodes[wpc * k];
            let in_im_fw  = self.nodes[wpc * k + 1];
            let in_l      = self.nodes[wpc * k + lam];
            let out_re_fw = self.nodes[out_base + wpc * k];
            let out_im_fw = self.nodes[out_base + wpc * k + 1];
            let out_l     = self.nodes[out_base + wpc * k + lam];
            let c = self.c_cached[k];
            let s = self.s_cached[k];
            stamp_potential_eq(mat, &self.branches, bpc * k,     out_re_fw,
                &[(in_re_fw, -c), (in_im_fw, -s)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + 1, out_im_fw,
                &[(in_re_fw,  s), (in_im_fw, -c)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + (bpc - 1), out_l,
                &[(in_l, -1.0)]);
            if wpc == 5 {
                let in_re_bw  = self.nodes[wpc * k + 2];
                let in_im_bw  = self.nodes[wpc * k + 3];
                let out_re_bw = self.nodes[out_base + wpc * k + 2];
                let out_im_bw = self.nodes[out_base + wpc * k + 3];
                stamp_potential_eq(mat, &self.branches, bpc * k + 2, in_re_bw,
                    &[(out_re_bw, -c), (out_im_bw, -s)]);
                stamp_potential_eq(mat, &self.branches, bpc * k + 3, in_im_bw,
                    &[(out_re_bw,  s), (out_im_bw, -c)]);
            }
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) { self.load_residual(b); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) { self.load_jacobian(mat); }
}

// ────────────────────────────────────────────────────────────────────────
// Native testbench MZM (idealised Mach-Zehnder modulator)
// ────────────────────────────────────────────────────────────────────────

/// Testbench MZM — represents an idealised lab-bench Mach-Zehnder modulator
/// (free-space or fibre-pigtailed), NOT a real on-chip MZI.  Use this to
/// model the modulator instrument feeding a chip, or as a reference against
/// which to compare a chip-level MZI built from `fc_splitter` +
/// `fc_pn_th_ps` + `fc_splitter` (the "real" MZI lives in the schematic).
///
/// Pins (4 + 6·N for N optical channels):
///   1  in     — optical bundle input
///   2  out    — optical bundle output
///   3  sig    — electrical drive (sig − gnd is the modulation voltage)
///   4  gnd    — modulation return
///
/// Intensity transmission as a function of V_sig:
///   T(V) = α · [ (1 − 1/E_r) · (1 + cos(π V_sig / V_π)) / 2  +  1/E_r ]
/// with α = 10^(−`alpha_dB`/10) (intensity) and E_r the extinction ratio.
/// T(V) ranges from α (V=0) to α/E_r (V = V_π).  Amplitude transmission
/// `t_amp = √T`.  Wavelength passes through unchanged.
///
/// `f_c` (first-order EO cutoff) is accepted but not yet wired into the
/// signal path — it requires device-internal reactive state which lands
/// with the L2 framework.  At this commit `V_sig` reaches the modulator
/// instantaneously regardless of f_c.
pub struct NativeMzm {
    v_pi:         f64,
    alpha_int:    f64,
    e_r:          f64,
    f_c:          f64,
    n_channels:   usize,
    wpc:          usize,
    nodes:        Vec<NodeId>,
    branches:     Vec<Option<usize>>,
    t_amp_cached: f64,
}

impl NativeMzm {
    pub fn new() -> Self {
        Self {
            v_pi:         3.0,
            alpha_int:    1.0,
            e_r:          1.0e3,
            f_c:          1.0e10,
            n_channels:   0,
            wpc:          3,
            nodes:        Vec::new(),
            branches:     Vec::new(),
            t_amp_cached: 1.0,
        }
    }
}

impl Device for NativeMzm {
    fn num_terminals(&self) -> usize { self.nodes.len() }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 2 * wpc;
        assert!(
            terminals.len() >= stride + 2 && (terminals.len() - 2) % stride == 0,
            "fc_mzm: terminal count must be {stride}·N + 2 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 2) / stride;
        self.n_channels = n;
        self.nodes      = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 };
        self.branches   = vec![None; bpc * n];
    }

    fn num_extra_nodes(&self) -> usize { self.branches.len() }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..self.branches.len() { self.branches[i] = Some(first_idx + i); }
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "v_pi" | "vpi" => { self.v_pi = value.max(1e-30); true }
            "alpha" => {
                self.alpha_int = value.clamp(0.0, 1.0);
                true
            }
            "alpha_db" | "il_db" => {
                self.alpha_int = 10f64.powf(-value / 10.0);
                true
            }
            "e_r" | "er" => { self.e_r = value.max(1.0); true }
            "e_r_db" | "er_db" => {
                self.e_r = 10f64.powf(value / 10.0).max(1.0);
                true
            }
            "f_c" | "fc" => { self.f_c = value.max(0.0); true }
            _ => false,
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        let elec_base = 2 * wpc * n;
        let v_sig = self.nodes[elec_base    ].map_or(0.0, |i| x[i]);
        let v_gnd = self.nodes[elec_base + 1].map_or(0.0, |i| x[i]);
        let v_mod = v_sig - v_gnd;
        let cos_term = (std::f64::consts::PI * v_mod / self.v_pi).cos();
        let inv_er   = 1.0 / self.e_r;
        let t_int    = self.alpha_int * ((1.0 - inv_er) * 0.5 * (1.0 + cos_term) + inv_er);
        self.t_amp_cached = t_int.max(0.0).sqrt();
    }

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        let bpc = if wpc == 5 { 5 } else { 3 };
        let lam = wpc - 1;
        let out_base = wpc * n;
        let t = self.t_amp_cached;
        for k in 0..n {
            let in_re_fw  = self.nodes[wpc * k];
            let in_im_fw  = self.nodes[wpc * k + 1];
            let in_l      = self.nodes[wpc * k + lam];
            let out_re_fw = self.nodes[out_base + wpc * k];
            let out_im_fw = self.nodes[out_base + wpc * k + 1];
            let out_l     = self.nodes[out_base + wpc * k + lam];
            stamp_potential_eq(mat, &self.branches, bpc * k,     out_re_fw, &[(in_re_fw, -t)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + 1, out_im_fw, &[(in_im_fw, -t)]);
            stamp_potential_eq(mat, &self.branches, bpc * k + (bpc - 1), out_l, &[(in_l, -1.0)]);
            if wpc == 5 {
                // MZM is reciprocal: same t_amp applies to backward path.
                let in_re_bw  = self.nodes[wpc * k + 2];
                let in_im_bw  = self.nodes[wpc * k + 3];
                let out_re_bw = self.nodes[out_base + wpc * k + 2];
                let out_im_bw = self.nodes[out_base + wpc * k + 3];
                stamp_potential_eq(mat, &self.branches, bpc * k + 2, in_re_bw, &[(out_re_bw, -t)]);
                stamp_potential_eq(mat, &self.branches, bpc * k + 3, in_im_bw, &[(out_im_bw, -t)]);
            }
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) { self.load_residual(b); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) { self.load_jacobian(mat); }
}

// ────────────────────────────────────────────────────────────────────────
// Native WDM multiplexer / demultiplexer
// ────────────────────────────────────────────────────────────────────────
//
// `fc_mux` / `fc_demux` bridge between N single-channel optical bundles and
// one N-channel optical bundle.  They are TOPOLOGY MARKERS, not signal
// processors: each device is identity-routing channel-by-channel
// (`bus[k].* = ch_k.*`).  The point is to give the schematic a single place
// where bundle widths change, so users can wire a wavelength-diverse circuit
// without dealing with KiCad's bus syntax (which can't connect directly to
// single symbol pins).
//
// Terminal layout (variable arity, derived in `setup_instance`):
//
//   fc_mux  N=4 has 6·N = 24 terminals.  The first 3·N are the bus output
//           wires interleaved per channel: [bus.0.re, bus.0.im, bus.0.λ,
//           bus.1.re, ..., bus.{N-1}.λ].  The next 3·N are the N single-
//           channel inputs in instance order: [ch0.re, ch0.im, ch0.λ,
//           ch1.re, ..., ch{N-1}.λ].
//   fc_demux same layout — bus first (now input), single channels next
//           (now outputs).
//
// The parser knows these two model names are "bundle-bridging" and must
// (a) skip the channel-count matching check and (b) emit a single instance
// with every bundle flattened to its underlying wires.  See
// `expand_optical_ports` in fairchild-parser.

/// Identity-routing combiner: N single-channel optical bundles → 1 N-channel
/// bundle.  Pin 1 (and the first bundle wire block) is the bus output.
pub struct NativeMux {
    n_channels: usize,
    wpc:        usize,
    nodes:      Vec<NodeId>,
    branches:   Vec<Option<usize>>,
}

impl NativeMux {
    pub fn new() -> Self {
        Self { n_channels: 0, wpc: 3, nodes: Vec::new(), branches: Vec::new() }
    }
}

impl Device for NativeMux {
    fn num_terminals(&self) -> usize { self.nodes.len() }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 2 * wpc; // bus channel (wpc) + per-channel wires (wpc)
        assert!(
            !terminals.is_empty() && terminals.len() % stride == 0,
            "fc_mux: terminal count must be a positive multiple of {stride} \
             (wpc={wpc}: bus wires + per-channel wires); got {}",
            terminals.len()
        );
        let n = terminals.len() / stride;
        self.n_channels = n;
        self.nodes      = terminals.to_vec();
        self.branches   = vec![None; wpc * n];
    }

    fn num_extra_nodes(&self) -> usize { self.branches.len() }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..self.branches.len() {
            self.branches[i] = Some(first_idx + i);
        }
    }

    fn eval(&mut self, _x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {}

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        for k in 0..n {
            for w in 0..wpc {
                let bus_w = self.nodes[wpc * k + w];
                let ch_w  = self.nodes[wpc * (n + k) + w];
                // Identity-route every wire (fw, bw, λ) — bus reads from channel.
                stamp_potential_eq(mat, &self.branches, wpc * k + w, bus_w,
                    &[(ch_w, -1.0)]);
            }
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) { self.load_residual(b); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) { self.load_jacobian(mat); }
}

/// Identity-routing splitter: 1 N-channel optical bundle → N single-channel
/// bundles.  Pin 1 (and the first bundle wire block) is the bus input.
pub struct NativeDemux {
    n_channels: usize,
    wpc:        usize,
    nodes:      Vec<NodeId>,
    branches:   Vec<Option<usize>>,
}

impl NativeDemux {
    pub fn new() -> Self {
        Self { n_channels: 0, wpc: 3, nodes: Vec::new(), branches: Vec::new() }
    }
}

impl Device for NativeDemux {
    fn num_terminals(&self) -> usize { self.nodes.len() }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 2 * wpc;
        assert!(
            !terminals.is_empty() && terminals.len() % stride == 0,
            "fc_demux: terminal count must be a positive multiple of {stride} \
             (wpc={wpc}: bus wires + per-channel wires); got {}",
            terminals.len()
        );
        let n = terminals.len() / stride;
        self.n_channels = n;
        self.nodes      = terminals.to_vec();
        self.branches   = vec![None; wpc * n];
    }

    fn num_extra_nodes(&self) -> usize { self.branches.len() }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..self.branches.len() {
            self.branches[i] = Some(first_idx + i);
        }
    }

    fn eval(&mut self, _x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {}

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        for k in 0..n {
            for w in 0..wpc {
                let bus_w = self.nodes[wpc * k + w];
                let ch_w  = self.nodes[wpc * (n + k) + w];
                // Channels drive FROM bus.
                stamp_potential_eq(mat, &self.branches, wpc * k + w, ch_w,
                    &[(bus_w, -1.0)]);
            }
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) { self.load_residual(b); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) { self.load_jacobian(mat); }
}

/// Stamp `V(out) = Σ k_i · V(in_i)` into one auxiliary branch row.
fn stamp_potential_eq(
    mat: &mut MnaMatrix,
    branches: &[Option<usize>],
    branch_idx: usize,
    out_node: NodeId,
    ins: &[(NodeId, f64)],
) {
    let (Some(out), Some(j)) = (out_node, branches[branch_idx]) else { return };
    mat.a[j][out] += 1.0;
    for &(in_node, k) in ins {
        if let Some(in_i) = in_node { mat.a[j][in_i] += k; }
    }
    mat.a[out][j] += 1.0;
}

// ────────────────────────────────────────────────────────────────────────
// Shared utilities
// ────────────────────────────────────────────────────────────────────────

#[allow(non_snake_case)]
fn dB_per_cm_to_neper_per_m(alpha_db_cm: f64) -> f64 {
    // 1 dB = ln(10)/20 Np ≈ 0.1151 Np; 1 cm = 0.01 m → multiply by 100/cm.
    alpha_db_cm * 100.0 * std::f64::consts::LN_10 / 20.0
}

#[cfg(test)]
mod tests {
    use crate::newton::dc_op_nr_with_registry;
    use crate::device_registry::DeviceRegistry;
    use fairchild_parser::parse_spice;

    /// Drive a native waveguide directly through voltage sources on its
    /// underlying re/im/λ wires and verify the output amplitude matches the
    /// closed-form `exp(-α·L/2) · exp(-jβL)` formula.
    ///
    /// `.optical` declarations are skipped here — we use XOsdi (which the
    /// discipline check exempts) for the native device and treat the wires
    /// as plain electrical nets for the voltage sources.
    #[test]
    fn native_waveguide_amplitude_matches_closed_form() {
        // L = 100 µm, n_g = 4.2, α = 2 dB/cm, λ = 1.55 µm.
        // T = exp(-α_Np/m · L / 2), where α_Np/m = 2·100·ln(10)/20 ≈ 23.026 Np/m.
        //   → T = exp(-23.026·100e-6/2) = exp(-1.1513e-3) ≈ 0.998850
        // φ = 2π·n_g·L/λ ≈ 2π·4.2·100/1550 wraps; |A_out| = T regardless of φ.
        let netlist = parse_spice(
            "* native waveguide test\n\
             V_re in_re 0 DC 1.0\n\
             V_im in_im 0 DC 0.0\n\
             V_wl in_wl 0 DC 1.55e-6\n\
             X1 in_re in_im in_wl out_re out_im out_wl fc_waveguide \
                L_um=100 n_g=4.2 alpha_db_cm=2.0\n\
             .op\n.end\n"
        ).unwrap();
        let registry = DeviceRegistry::new();
        let r = dc_op_nr_with_registry(&netlist, &registry)
            .expect("DC OP should converge");
        let v_re = r.node_voltage("out_re").unwrap();
        let v_im = r.node_voltage("out_im").unwrap();
        let amp  = (v_re * v_re + v_im * v_im).sqrt();
        let expected = (-23.0258509_f64 * 100e-6 / 2.0).exp(); // 0.99885
        assert!((amp - expected).abs() < 1e-5,
            "|A_out|={amp:.6} expected={expected:.6}");
        // Output wavelength must equal input wavelength.
        let v_wl = r.node_voltage("out_wl").unwrap();
        assert!((v_wl - 1.55e-6).abs() < 1e-15);
    }

    /// 3 dB splitter — equal-power split, in-phase.  |c|² + |d|² = |a|²,
    /// and c = d (both halves of the input).
    #[test]
    fn native_splitter_equal_power_split() {
        let netlist = parse_spice(
            "* splitter test\n\
             V_re a_re 0 DC 1.0\n\
             V_im a_im 0 DC 0.0\n\
             V_wl a_wl 0 DC 1.55e-6\n\
             X1 a_re a_im a_wl c_re c_im c_wl d_re d_im d_wl fc_splitter\n\
             .op\n.end\n"
        ).unwrap();
        let r = dc_op_nr_with_registry(&netlist, &DeviceRegistry::new())
            .expect("DC OP should converge");
        let c_re = r.node_voltage("c_re").unwrap();
        let c_im = r.node_voltage("c_im").unwrap();
        let d_re = r.node_voltage("d_re").unwrap();
        let d_im = r.node_voltage("d_im").unwrap();
        // Each output: a / √2 ≈ 0.7071
        let expected_re = 1.0 / 2.0_f64.sqrt();
        assert!((c_re - expected_re).abs() < 1e-9);
        assert!((d_re - expected_re).abs() < 1e-9);
        assert!(c_im.abs() < 1e-9);
        assert!(d_im.abs() < 1e-9);
        // Power conservation: c² + d² ≈ 1 (input was 1.0)
        let p_total = c_re * c_re + c_im * c_im + d_re * d_re + d_im * d_im;
        assert!((p_total - 1.0).abs() < 1e-9);
    }

    #[test]
    fn native_cw_laser_drives_output_potentials() {
        let netlist = parse_spice(
            "* laser test\n\
             X1 out_re out_im out_wl fc_cw_laser \
                power_mW=4.0 phi_0_deg=0.0 wavelength_nm=1550\n\
             .op\n.end\n"
        ).unwrap();
        let r = dc_op_nr_with_registry(&netlist, &DeviceRegistry::new()).unwrap();
        // P = 4 mW → A = √(4e-3) ≈ 0.06325 V/m equivalent.
        let v_re = r.node_voltage("out_re").unwrap();
        let v_im = r.node_voltage("out_im").unwrap();
        let v_wl = r.node_voltage("out_wl").unwrap();
        let expected_amp = 4e-3_f64.sqrt();
        assert!((v_re - expected_amp).abs() < 1e-9, "v_re={v_re}");
        assert!(v_im.abs() < 1e-9);
        assert!((v_wl - 1.55e-6).abs() < 1e-15);
    }

    /// PIN photodetector under a reverse bias.  Optical input drives V(in_re) =
    /// 1, V(in_im) = 0 → P = 1 W; responsivity = 0.8 A/W; expected
    /// photocurrent ≈ 0.8 A flowing cathode → anode.  Verifies by reading
    /// the anode voltage through a load resistor to ground.
    #[test]
    fn native_photodetector_produces_responsivity_current() {
        let netlist = parse_spice(
            "* PD test\n\
             V_re in_re 0 DC 1.0\n\
             V_im in_im 0 DC 0.0\n\
             V_wl in_wl 0 DC 1.55e-6\n\
             V_bias bias 0 DC 1.0\n\
             R_load anode bias 1k\n\
             X1 in_re in_im in_wl anode 0 fc_photodetector \
                responsivity=0.8 i_dark_a=1e-12 r_shunt=1e6\n\
             .op\n.end\n"
        ).unwrap();
        let r = dc_op_nr_with_registry(&netlist, &DeviceRegistry::new())
            .expect("DC OP should converge");
        // P_opt = 1 W; I_ph = 0.8 A flowing cathode→anode in the device frame.
        // Through R_load = 1k from anode to bias=1V, V(anode) settles so that
        // (V(anode) − 1) / 1k + (small shunt) ≈ I_ph.  V(anode) ≈ 1 + 800 V
        // for I_ph = 0.8 A — that's a clipped value in real circuits but
        // mathematically the linear stamp produces it.  Use a tiny power
        // instead for sane numbers:
        let v_anode = r.node_voltage("anode").unwrap();
        // For now just assert: the anode is significantly above bias
        // (current was pushed into the load).
        assert!(v_anode > 1.5, "v_anode = {v_anode} should be > bias (1V)");
    }

    /// Thermal phase shifter at V = 0 has zero phase shift → output = input.
    /// At V = Vπ (chosen so that V²/R = P_pi), phase shift = π → output = −input.
    #[test]
    fn native_thermal_ps_zero_voltage_passthrough() {
        let netlist = parse_spice(
            "* thermal PS at V=0\n\
             V_re in_re 0 DC 1.0\n\
             V_im in_im 0 DC 0.0\n\
             V_wl in_wl 0 DC 1.55e-6\n\
             V_heat heat 0 DC 0.0\n\
             X1 in_re in_im in_wl out_re out_im out_wl heat 0 fc_thermal_ps \
                r_heater=1k p_pi=10m\n\
             .op\n.end\n"
        ).unwrap();
        let r = dc_op_nr_with_registry(&netlist, &DeviceRegistry::new()).unwrap();
        let v_re = r.node_voltage("out_re").unwrap();
        let v_im = r.node_voltage("out_im").unwrap();
        assert!((v_re - 1.0).abs() < 1e-9, "zero-V should pass input through: out_re={v_re}");
        assert!(v_im.abs() < 1e-9);
    }

    #[test]
    fn native_thermal_ps_at_v_pi_inverts() {
        // V_pi = sqrt(P_pi · R) for V²/R = P_pi.  P_pi=10m, R=1k → V_pi=√10 ≈ 3.162.
        let v_pi = (10e-3 * 1000.0_f64).sqrt();
        let netlist = parse_spice(&format!(
            "* thermal PS at V_pi\n\
             V_re in_re 0 DC 1.0\n\
             V_im in_im 0 DC 0.0\n\
             V_wl in_wl 0 DC 1.55e-6\n\
             V_heat heat 0 DC {v_pi}\n\
             X1 in_re in_im in_wl out_re out_im out_wl heat 0 fc_thermal_ps \
                r_heater=1k p_pi=10m\n\
             .op\n.end\n"
        )).unwrap();
        let r = dc_op_nr_with_registry(&netlist, &DeviceRegistry::new()).unwrap();
        let v_re = r.node_voltage("out_re").unwrap();
        let v_im = r.node_voltage("out_im").unwrap();
        // φ = π → exp(-jπ)·(1+0j) = -1 → out_re = -1, out_im = 0.
        assert!((v_re + 1.0).abs() < 1e-6, "at Vπ out_re should be ≈ -1: got {v_re}");
        assert!(v_im.abs() < 1e-6, "at Vπ out_im should be ≈ 0: got {v_im}");
    }

    /// PN phase shifter: zero bias → identity passthrough.
    #[test]
    fn native_pn_ps_zero_bias_passthrough() {
        let netlist = parse_spice(
            "* PN PS at V=0\n\
             V_re in_re 0 DC 1.0\n\
             V_im in_im 0 DC 0.0\n\
             V_wl in_wl 0 DC 1.55e-6\n\
             V_bias bias 0 DC 0.0\n\
             X1 in_re in_im in_wl out_re out_im out_wl bias 0 fc_pn_ps \
                L_um=1000 V_pi_L=2e-3\n\
             .op\n.end\n"
        ).unwrap();
        let r = dc_op_nr_with_registry(&netlist, &DeviceRegistry::new()).unwrap();
        let v_re = r.node_voltage("out_re").unwrap();
        let v_im = r.node_voltage("out_im").unwrap();
        assert!((v_re - 1.0).abs() < 1e-9);
        assert!(v_im.abs() < 1e-9);
    }

    /// Directional coupler: at κL = π/4, transmission and coupling are equal
    /// (cos²(π/4) = sin²(π/4) = 0.5).  With input only on port a, the cross
    /// port d gets power transferred via the imaginary cross-coupling.
    #[test]
    fn native_dcoupler_half_power_at_kl_pi_over_4() {
        let netlist = parse_spice(
            "* dcoupler test\n\
             V_ar a_re 0 DC 1.0\n\
             V_ai a_im 0 DC 0.0\n\
             V_awl a_wl 0 DC 1.55e-6\n\
             V_br b_re 0 DC 0.0\n\
             V_bi b_im 0 DC 0.0\n\
             V_bwl b_wl 0 DC 1.55e-6\n\
             X1 a_re a_im a_wl b_re b_im b_wl \
                c_re c_im c_wl d_re d_im d_wl fc_dcoupler kappa_L=0.7853981633974483\n\
             .op\n.end\n"
        ).unwrap();
        let r = dc_op_nr_with_registry(&netlist, &DeviceRegistry::new())
            .expect("DC OP should converge");
        // With a=(1,0), b=(0,0), κL=π/4 → t=s=1/√2:
        //   c_re = t·a_re = 1/√2,    c_im = -s·b_re = 0
        //   d_re = t·b_re + s·a_im = 0,   d_im = -s·a_re = -1/√2
        let half = 1.0 / 2.0_f64.sqrt();
        let c_re = r.node_voltage("c_re").unwrap();
        let c_im = r.node_voltage("c_im").unwrap();
        let d_re = r.node_voltage("d_re").unwrap();
        let d_im = r.node_voltage("d_im").unwrap();
        assert!((c_re - half).abs() < 1e-9, "c_re={c_re} expected {half}");
        assert!(c_im.abs() < 1e-9);
        assert!(d_re.abs() < 1e-9);
        assert!((d_im + half).abs() < 1e-9, "d_im={d_im} expected {}", -half);
        // Power: |c|² + |d|² = 1 (lossless).
        let p = c_re * c_re + c_im * c_im + d_re * d_re + d_im * d_im;
        assert!((p - 1.0).abs() < 1e-9, "p={p}");
    }
}
