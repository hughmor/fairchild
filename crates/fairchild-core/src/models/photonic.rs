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

use crate::device::{Device, EvalFlags, NodeId, ReactiveBranchSpec, ReactiveKind, SimContext};
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
    n_eff:           f64,    // n_eff at wl_ref_m (default 2.445 for silicon)
    n_g:             f64,    // n_g    at wl_ref_m (default 4.2)
    /// Reference wavelength at which `n_eff` and `n_g` are evaluated.  The
    /// dispersion-corrected index `n_eff(λ)` is linearised around this point.
    wl_ref_m:        f64,
    alpha_neper_m:   f64,
    /// Group delay `τ_g = L·n_g/c` (s).  Computed and exposed for transient
    /// post-processing; this device does not yet implement a true delay line,
    /// so the parameter is informational only at this tier (DC OP and steady-
    /// state spectra are unaffected — τ matters only at modulation
    /// bandwidths comparable to 1/τ).
    tau_g_s:         f64,
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
        // Defaults: classic 500 × 220 nm SOI strip waveguide, straight.
        //  n_eff / n_g at 1550 nm extracted from femwell (see
        //  `scripts/waveguide_simulations/cband_sweep.csv`, strip column).
        //  Phase-shifter device classes (fc_pn_ps, fc_pn_th_ps, fc_pn_ps_cap)
        //  use bent-rib values appropriate to a ring section instead.
        let length_m = 100e-6;
        let n_g      = 4.19;
        NativeWaveguide {
            length_m,
            n_eff:              2.445,
            n_g,
            wl_ref_m:           1.55e-6,
            alpha_neper_m:      dB_per_cm_to_neper_per_m(2.0),
            tau_g_s:            length_m * n_g / C0,
            lambda_bootstrap_m: 1.55e-6,
            n_channels:         0,
            wpc:                3,
            nodes:              Vec::new(),
            branches:           Vec::new(),
            c_cached:           Vec::new(),
            s_cached:           Vec::new(),
        }
    }

    fn refresh_tau(&mut self) { self.tau_g_s = self.length_m * self.n_g / C0; }
}

impl Device for NativeWaveguide {
    fn num_terminals(&self) -> usize { self.nodes.len() }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.lambda_bootstrap_m = ctx.lambda_center_m;
        // Default the dispersion reference to the band centre too; the user
        // can still override via `wl_ref_nm=…`.
        if (self.wl_ref_m - 1.55e-6).abs() < 1e-12 {
            self.wl_ref_m = ctx.lambda_center_m;
        }
        self.wpc = ctx.wires_per_channel();
        self.refresh_tau();
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
            "l_um"          => { self.length_m       = value * 1e-6; self.refresh_tau();  true }
            "l_m" | "length"=> { self.length_m       = value;        self.refresh_tau();  true }
            "n_g"           => { self.n_g            = value;        self.refresh_tau();  true }
            "n_eff"         => { self.n_eff          = value;                              true }
            "wl_ref_m" | "lambda_ref_m" => { self.wl_ref_m = value;                        true }
            "wl_ref_nm" | "lambda_ref_nm" => { self.wl_ref_m = value * 1e-9;               true }
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
            // Dispersion-corrected effective index for accumulated phase.
            let n_eff_lam = n_eff_at_lambda(self.n_eff, self.n_g, self.wl_ref_m, lambda);
            let phi = two_pi * n_eff_lam * self.length_m / lambda;
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
// Native thermal phase shifter with thermal time constant (L2 — path B)
// ────────────────────────────────────────────────────────────────────────

/// Thermal phase shifter with a first-order thermal RC.  Same pin layout
/// as `fc_thermal_ps`.  Adds:
///   - `tau_th` — thermal time constant (s).  The optical phase shift
///     tracks the FILTERED heater power rather than the instantaneous
///     Joule dissipation: `dT/dt = (P − T) / tau_th`, with T in
///     normalised "power-equivalent" units so the steady-state phase
///     equals the L1 model's φ = π · P / P_pi.
///
/// Implementation: T is a *state variable* on the MNA matrix — the
/// device allocates one extra row (above its optical branches) and
/// stamps the BE-discretised state equation directly.  This is the
/// "path B" pattern for nonlinear / nonlinear-coupled state, complementary
/// to the linear-companion path A (used by `fc_pn_ps_cap` for C_j(V)).
/// The previous-timestep T is captured via `commit_timestep` after each
/// successful NR convergence.
pub struct NativeThermalPhaseShifterRc {
    r_heater:   f64,
    p_pi:       f64,
    tau_th:     f64,
    n_channels: usize,
    wpc:        usize,
    nodes:      Vec<NodeId>,
    /// Optical branch rows (re/im/λ per channel).  Same as L1.
    branches:   Vec<Option<usize>>,
    /// State row for T(t) (single MNA index allocated alongside branches).
    t_state_idx: Option<usize>,
    /// Previous-timestep value of T, captured by `commit_timestep`.
    t_old:      f64,
    /// Cached operating-point quantities (per NR iteration).
    t_op:       f64,
    v_h_op:     f64,
    c_cached:   f64,
    s_cached:   f64,
}

impl NativeThermalPhaseShifterRc {
    pub fn new() -> Self {
        Self {
            r_heater:    1000.0,
            p_pi:        10e-3,
            tau_th:      10e-6,    // 10 µs typical waveguide heater
            n_channels:  0,
            wpc:         3,
            nodes:       Vec::new(),
            branches:    Vec::new(),
            t_state_idx: None,
            t_old:       0.0,
            t_op:        0.0,
            v_h_op:      0.0,
            c_cached:    1.0,
            s_cached:    0.0,
        }
    }
}

impl Device for NativeThermalPhaseShifterRc {
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
            "fc_thermal_ps_rc: terminal count must be {stride}·N + 2 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 2) / stride;
        self.n_channels = n;
        self.nodes      = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 };
        self.branches   = vec![None; bpc * n];
    }

    fn num_extra_nodes(&self) -> usize {
        // Optical branches + 1 state row for T(t).
        self.branches.len() + 1
    }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        let n = self.branches.len();
        for i in 0..n { self.branches[i] = Some(first_idx + i); }
        self.t_state_idx = Some(first_idx + n);
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "r_heater" | "r" => { self.r_heater = value; true }
            "p_pi" | "p_pi_w" => { self.p_pi = value; true }
            "tau_th" | "tau" => { self.tau_th = value.max(1e-30); true }
            _ => false,
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        let elec_base = 2 * wpc * n;
        let v_a = self.nodes[elec_base    ].map_or(0.0, |i| x[i]);
        let v_c = self.nodes[elec_base + 1].map_or(0.0, |i| x[i]);
        self.v_h_op = v_a - v_c;
        // Read T from the state row.
        self.t_op = self.t_state_idx.map_or(0.0, |i| x[i]);
        // Phase from T (filtered power).
        let phi = std::f64::consts::PI * self.t_op / self.p_pi;
        self.c_cached = phi.cos();
        self.s_cached = phi.sin();
    }

    fn load_residual(&self, b: &mut [f64]) {
        // DC: T = P (steady state).  Linearised: T − P_lin(V_h) = 0
        // where P_lin = 2·V_h_op·V_h/R − V_h_op²/R, so
        //   T − 2·V_h_op·V_h/R + V_h_op²/R = 0
        //   ⇒ rearranged with V_h_op² constant on RHS:  b[t_idx] = −V_h_op²/R.
        if let Some(t_idx) = self.t_state_idx {
            let p_op = self.v_h_op * self.v_h_op / self.r_heater;
            b[t_idx] -= p_op;
        }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        let bpc = if wpc == 5 { 5 } else { 3 };
        let lam = wpc - 1;
        let out_base  = wpc * n;
        let elec_base = 2 * wpc * n;
        // Electrical: shared heater conductance.
        let g = 1.0 / self.r_heater;
        let p = self.nodes[elec_base];
        let m = self.nodes[elec_base + 1];
        if let Some(a) = p { mat.a[a][a] += g; if let Some(c) = m { mat.a[a][c] -= g; } }
        if let Some(c) = m { mat.a[c][c] += g; if let Some(a) = p { mat.a[c][a] -= g; } }
        // State row for T (DC: T = P, i.e., T - P_linearised = 0).
        // Stamp: row = +1·T - (2·V_h_op/R)·V_hp + (2·V_h_op/R)·V_hn = +V_h_op²/R
        if let Some(t_idx) = self.t_state_idx {
            mat.a[t_idx][t_idx] += 1.0;
            let two_vop_over_r = 2.0 * self.v_h_op / self.r_heater;
            if let Some(hp) = p { mat.a[t_idx][hp] -= two_vop_over_r; }
            if let Some(hn) = m { mat.a[t_idx][hn] += two_vop_over_r; }
        }
        // Optical branches: identical structure to fc_thermal_ps but using
        // c_cached/s_cached derived from T (state).
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

    fn load_residual_tran(&self, b: &mut [f64], alpha: f64) {
        // BE-discretised state-row RHS:  b[t_idx] = T_old·α − V_h_op²/(R·τ).
        // (The −V_h_op² term is the linearisation remainder of P_lin =
        //  2·V_h_op·V_h/R − V_h_op²/R; the 2·V_h_op·V_h part is in the
        //  Jacobian.)  Optical branch rows are homogeneous (no residual).
        if let Some(t_idx) = self.t_state_idx {
            let inv_tau = 1.0 / self.tau_th;
            let p_op = self.v_h_op * self.v_h_op / self.r_heater;
            b[t_idx] += self.t_old * alpha - p_op * inv_tau;
        }
    }

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, alpha: f64) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        let bpc = if wpc == 5 { 5 } else { 3 };
        let lam = wpc - 1;
        let out_base  = wpc * n;
        let elec_base = 2 * wpc * n;
        // Electrical: shared heater conductance.
        let g = 1.0 / self.r_heater;
        let p = self.nodes[elec_base];
        let m = self.nodes[elec_base + 1];
        if let Some(a) = p { mat.a[a][a] += g; if let Some(c) = m { mat.a[a][c] -= g; } }
        if let Some(c) = m { mat.a[c][c] += g; if let Some(a) = p { mat.a[c][a] -= g; } }
        // BE state-row Jacobian: T_new·(α + 1/τ) − 2·V_h_op·V_h/(R·τ) = …
        if let Some(t_idx) = self.t_state_idx {
            let inv_tau = 1.0 / self.tau_th;
            mat.a[t_idx][t_idx] += alpha + inv_tau;
            let two_vop_over_r = 2.0 * self.v_h_op / self.r_heater;
            if let Some(hp) = p { mat.a[t_idx][hp] -= two_vop_over_r * inv_tau; }
            if let Some(hn) = m { mat.a[t_idx][hn] += two_vop_over_r * inv_tau; }
        }
        // Optical branches (same as DC, c/s from cached T).
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

    fn commit_timestep(&mut self, x: &[f64]) {
        // Snapshot T for use as T_old on the next timestep's BE stamp.
        if let Some(idx) = self.t_state_idx {
            self.t_old = x[idx];
        }
    }
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
    n_eff:    f64,
    n_g:      f64,
    wl_ref_m: f64,
    dn_dv:    f64,
    g_pn:     f64,
    alpha_neper_m: f64,
    /// When true (default), subtract the absolute propagation phase at
    /// `wl_ref_m` so the device is "transparent" at λ = λ_ref.  Convenient
    /// for testbench rings where the user wants the ring on-resonance at the
    /// laser wavelength by construction.  Set to false for multi-ring designs
    /// (rings of different L) where you want each ring's natural absolute
    /// resonance position — otherwise all rings cluster at λ_ref regardless
    /// of length.  Set via `pin_at_ref=0|1` SPICE parameter.
    pin_at_ref: bool,
    n_channels: usize,
    wpc:      usize,                 // 3 or 5
    nodes:    Vec<NodeId>,
    branches: Vec<Option<usize>>,
    c_cached: Vec<f64>,
    s_cached: Vec<f64>,
}

impl NativePnPhaseShifter {
    pub fn new() -> Self {
        // Defaults: SOI rib waveguide PN modulator section, R = 8 µm bent.
        //  n_eff / n_g from `scripts/waveguide_simulations/cband_sweep.csv`
        //   (rib_R8 column at 1550 nm).
        //  alpha = 20 dB/cm is dominated by free-carrier absorption from the
        //   typical 5e17 cm⁻³ slab doping; replace with whatever the
        //   `pn_modulator/` sims report for your specific doping profile.
        //  V_pi_L default → dn_dv = wl_ref/(2·V_pi_L).  Set V_pi_L = 0.015
        //   (V·m) so V_pi = 0.015 / L_um·1e-6.  At the typical 1-mm PN
        //   length, V_pi ≈ 15 V (reverse-bias depletion-mode).
        //  pin_at_ref = false: use physical absolute propagation phase so
        //   ring resonances depend on L.  Pin to a ref-wavelength using
        //   `pin_at_ref=1` for testbench rings designed on-resonance.
        Self {
            length_m: 1e-3,
            n_eff:    2.7654,
            n_g:      4.02,
            wl_ref_m: 1.55e-6,
            dn_dv:    1.55e-6 / (2.0 * 0.015),
            g_pn:     1e-3,
            alpha_neper_m: dB_per_cm_to_neper_per_m(20.0),
            pin_at_ref: false,
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
            "n_eff"  => { self.n_eff = value; true }
            "wl_ref_m" | "lambda_ref_m" => { self.wl_ref_m = value; true }
            "wl_ref_nm" | "lambda_ref_nm" => { self.wl_ref_m = value * 1e-9; true }
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
            "pin_at_ref" => {
                self.pin_at_ref = value != 0.0;
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
        // Reference absolute propagation phase at λ_ref.  When `pin_at_ref`
        // is on we subtract it so the device is "transparent" at λ = λ_ref;
        // otherwise we use the full absolute phase so each ring's natural
        // resonance position depends on L (correct physics for multi-ring
        // designs).
        let phi_ref = if self.pin_at_ref {
            two_pi * self.n_eff * self.length_m / self.wl_ref_m
        } else { 0.0 };
        for k in 0..n {
            let lambda = match self.nodes[wpc * k + lam] {
                Some(i) => {
                    let v = x[i];
                    if v.abs() > 1e-9 { v } else { self.wl_ref_m }
                }
                None => self.wl_ref_m,
            };
            let n_eff_lam = n_eff_at_lambda(self.n_eff, self.n_g, self.wl_ref_m, lambda);
            let phi_abs  = two_pi * n_eff_lam * self.length_m / lambda;
            let phi_prop = phi_abs - phi_ref;
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
// Native PN-junction phase shifter — depletion (L2 with bias-dependent C_j)
// ────────────────────────────────────────────────────────────────────────

/// Depletion-mode PN-junction phase shifter with bias-dependent junction
/// capacitance.  Same optical/electrical port layout and stamping as
/// `fc_pn_ps`, plus:
///   - `C_j(V) = C_j0 / (1 − V_pn/V_bi)^m_j` for V_pn ≤ V_bi/2 (depletion),
///     linearly continued above the singularity to keep NR stable when
///     the user wanders into forward bias.  The integrator owns the
///     companion-model state for this capacitance via
///     `Device::reactive_branches`.
///   - `da/dV` — linear loss-vs-bias coefficient (free-carrier absorption
///     in moderate forward bias).
///
/// Tier convention: this is a separate device class, not a `level=` switch
/// on `fc_pn_ps`.  Forward-injection physics (higher dn/dV, large da/dV,
/// carrier-injection time constants) belongs in a future `fc_pn_ps_inj`.
pub struct NativePnPhaseShifterCap {
    length_m:      f64,
    n_eff:         f64,
    n_g:           f64,
    wl_ref_m:      f64,
    /// See `NativePnPhaseShifter::pin_at_ref`.  Default true.
    pin_at_ref:    bool,
    dn_dv:         f64,
    g_pn:          f64,
    alpha_neper_m: f64,
    // Bias-dependent C_j parameters.
    c_j0:          f64,    // F at V_pn = 0
    v_bi:          f64,    // V — built-in voltage (knee)
    m_j:           f64,    // grading coefficient
    // Linear da/dV loss-vs-bias (Np/m per V).  Default 0.
    da_dv:         f64,
    n_channels:    usize,
    wpc:           usize,
    nodes:         Vec<NodeId>,
    branches:      Vec<Option<usize>>,
    c_cached:      Vec<f64>,
    s_cached:      Vec<f64>,
    // Cached per-NR-iteration values used to re-feed the integrator's
    // companion model AND the per-channel optical loss factor.
    c_j_cached:    f64,
    alpha_eff_neper_m: f64,
}

impl NativePnPhaseShifterCap {
    pub fn new() -> Self {
        // Same baseline defaults as `fc_pn_ps` (bent rib SOI PN modulator).
        // Adds: depletion-mode C_j(V) (`c_j0`, `v_bi`, `m_j`) and a
        // linear reverse-bias loss-vs-bias coefficient (`da_dv`).
        Self {
            length_m:      1e-3,
            n_eff:         2.7654,
            n_g:           4.02,
            wl_ref_m:      1.55e-6,
            pin_at_ref:    false,
            dn_dv:         1.55e-6 / (2.0 * 0.015),
            g_pn:          1e-3,
            alpha_neper_m: dB_per_cm_to_neper_per_m(20.0),
            c_j0:          20e-15,
            v_bi:          0.7,
            m_j:           0.5,
            da_dv:         0.0,
            n_channels:    0,
            wpc:           3,
            nodes:         Vec::new(),
            branches:      Vec::new(),
            c_cached:      Vec::new(),
            s_cached:      Vec::new(),
            c_j_cached:    20e-15,
            alpha_eff_neper_m: 0.0,
        }
    }
}

impl Device for NativePnPhaseShifterCap {
    fn num_terminals(&self) -> usize { self.nodes.len() }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wl_ref_m = ctx.lambda_center_m;
        self.wpc      = ctx.wires_per_channel();
        self.alpha_eff_neper_m = self.alpha_neper_m;
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 2 * wpc;
        assert!(
            terminals.len() >= stride + 2 && (terminals.len() - 2) % stride == 0,
            "fc_pn_ps_cap: terminal count must be {stride}·N + 2 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 2) / stride;
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
            "l_um"   => { self.length_m = value * 1e-6; true }
            "l_m" | "length" => { self.length_m = value; true }
            "n_g"    => { self.n_g = value; true }
            "n_eff"  => { self.n_eff = value; true }
            "wl_ref_m" | "lambda_ref_m" => { self.wl_ref_m = value; true }
            "wl_ref_nm" | "lambda_ref_nm" => { self.wl_ref_m = value * 1e-9; true }
            "dn_dv"  => { self.dn_dv = value; true }
            "g_pn"   => { self.g_pn  = value; true }
            "v_pi_l" => {
                if value > 0.0 { self.dn_dv = self.wl_ref_m / (2.0 * value); }
                true
            }
            "alpha_db_cm" => {
                self.alpha_neper_m = dB_per_cm_to_neper_per_m(value);
                true
            }
            "c_j0"   => { self.c_j0 = value.max(0.0); true }
            "v_bi"   => { self.v_bi = value.max(1e-3); true }
            "m_j"    => { self.m_j  = value.clamp(0.0, 0.99); true }
            "da_dv"  => { self.da_dv = value; true }
            "pin_at_ref" => { self.pin_at_ref = value != 0.0; true }
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

        // C_j(V_pn) with linear continuation above the depletion singularity.
        // Knee chosen at V_bi/2 — matches SPICE diode convention.
        let v_knee = 0.5 * self.v_bi;
        self.c_j_cached = if v_pn < v_knee {
            self.c_j0 / (1.0 - v_pn / self.v_bi).powf(self.m_j)
        } else {
            // Linear extrapolation: c(V_knee) + (dc/dV at knee) · (V_pn − V_knee).
            let c_knee = self.c_j0 / (1.0 - v_knee / self.v_bi).powf(self.m_j);
            let dc_dv  = c_knee * self.m_j / (self.v_bi - v_knee);
            c_knee + dc_dv * (v_pn - v_knee)
        };
        // Bias-dependent loss: α(V) = α_0 + (da/dV) · max(0, −V_pn) — only
        // adds extra absorption in reverse bias (free-carrier-like).
        let v_rev = (-v_pn).max(0.0);
        self.alpha_eff_neper_m = self.alpha_neper_m + self.da_dv * v_rev;
        let t_amp  = (-self.alpha_eff_neper_m * self.length_m / 2.0).exp();

        let two_pi = 2.0 * std::f64::consts::PI;
        let lam = wpc - 1;
        let phi_ref = if self.pin_at_ref {
            two_pi * self.n_eff * self.length_m / self.wl_ref_m
        } else { 0.0 };
        for k in 0..n {
            let lambda = match self.nodes[wpc * k + lam] {
                Some(i) => {
                    let v = x[i];
                    if v.abs() > 1e-9 { v } else { self.wl_ref_m }
                }
                None => self.wl_ref_m,
            };
            let n_eff_lam = n_eff_at_lambda(self.n_eff, self.n_g, self.wl_ref_m, lambda);
            let phi_abs  = two_pi * n_eff_lam * self.length_m / lambda;
            let phi_prop = phi_abs - phi_ref;
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

    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> {
        // ONE shared depletion capacitance between anode and cathode (the
        // single physical junction).  Bias-dependent value re-queried per
        // NR iteration; the integrator owns the companion-model state.
        let elec_base = 2 * self.wpc * self.n_channels;
        let anode = self.nodes.get(elec_base    ).copied().flatten();
        let cath  = self.nodes.get(elec_base + 1).copied().flatten();
        vec![ReactiveBranchSpec {
            kind:  ReactiveKind::Capacitor,
            pos:   anode,
            neg:   cath,
            value: self.c_j_cached,
        }]
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
    n_eff:           f64,
    n_g:             f64,
    wl_ref_m:        f64,
    /// See `NativePnPhaseShifter::pin_at_ref`.  Default true.
    pin_at_ref:      bool,
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
        // Same baseline as `fc_pn_ps`; adds a linear thermal heater pair.
        Self {
            length_m:        1e-3,
            n_eff:           2.7654,
            n_g:             4.02,
            wl_ref_m:        1.55e-6,
            pin_at_ref:      false,
            dn_dv:           1.55e-6 / (2.0 * 0.015),
            g_pn:            1e-3,
            r_heater:        1000.0,
            p_pi_th:         10e-3,
            alpha_neper_m:   dB_per_cm_to_neper_per_m(20.0),
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
            "n_eff"           => { self.n_eff    = value;        true }
            "wl_ref_m" | "lambda_ref_m"  => { self.wl_ref_m = value;        true }
            "wl_ref_nm" | "lambda_ref_nm" => { self.wl_ref_m = value * 1e-9; true }
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
            "pin_at_ref"      => { self.pin_at_ref = value != 0.0; true }
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
        let phi_ref = if self.pin_at_ref {
            two_pi * self.n_eff * self.length_m / self.wl_ref_m
        } else { 0.0 };
        for k in 0..n {
            let lambda = match self.nodes[wpc * k + lam] {
                Some(i) => {
                    let v = x[i];
                    if v.abs() > 1e-9 { v } else { self.wl_ref_m }
                }
                None => self.wl_ref_m,
            };
            let n_eff_lam = n_eff_at_lambda(self.n_eff, self.n_g, self.wl_ref_m, lambda);
            let phi_abs  = two_pi * n_eff_lam * self.length_m / lambda;
            let phi_prop = phi_abs - phi_ref;
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
// Native PN+thermal phase shifter, depletion-mode (L2a-with-heater)
// Combines fc_pn_ps_cap (C_j(V), reverse-bias FCA via da_dv) with a heater
// pair (heat_p, heat_n) for thermal trim/biasing on top of the EO shift.
// ────────────────────────────────────────────────────────────────────────

/// Lateral PN modulator + thermal trim heater.  Operates in the reverse-bias
/// regime (V_pn ≤ 0); has depletion C_j(V) and a linear da_dv loss-vs-bias
/// like `fc_pn_ps_cap`, plus heater terminals identical to `fc_pn_th_ps`.
///
/// Terminal layout (N optical channels): 2·wpc·N + 4
///   in.0.re … in.{N-1}.λ, out.0.re … out.{N-1}.λ, anode, cathode, heat_p, heat_n
pub struct NativePnThermalPhaseShifterCap {
    length_m:      f64,
    n_eff:         f64,
    n_g:           f64,
    wl_ref_m:      f64,
    pin_at_ref:    bool,
    // PN side
    dn_dv:         f64,
    g_pn:          f64,
    alpha_neper_m: f64,
    c_j0:          f64,
    v_bi:          f64,
    m_j:           f64,
    da_dv:         f64,
    // Heater side
    r_heater:      f64,
    p_pi_th:       f64,
    n_channels:    usize,
    wpc:           usize,
    nodes:         Vec<NodeId>,
    branches:      Vec<Option<usize>>,
    c_cached:      Vec<f64>,
    s_cached:      Vec<f64>,
    c_j_cached:    f64,
    alpha_eff_neper_m: f64,
}

impl NativePnThermalPhaseShifterCap {
    pub fn new() -> Self {
        // Defaults from pn_modulator.py extraction (5e17/5e17 lateral PN,
        // 100 nm offset, 300 K).
        Self {
            length_m:      1e-3,
            n_eff:         2.7654,
            n_g:           4.02,
            wl_ref_m:      1.55e-6,
            pin_at_ref:    false,
            dn_dv:         5.024e-5,                     // depletion-mode linear coeff
            g_pn:          1e-3,
            alpha_neper_m: 29.78,                         // 2.59 dB/cm at V=0 (FCA, sim)
            c_j0:          1.375e-16 * 1e3,               // F per µm × µm → F-ish per 1 mm L (= 1.375e-13 F/m); see `length_m` default
            v_bi:          0.917,                          // V_bi @ N_A=N_D=5e17, 300 K
            m_j:           0.5,
            da_dv:         7.83,                           // slope: Δα/ΔV ≈ 7.8 Np/m per V
            r_heater:      1000.0,
            p_pi_th:       10e-3,
            n_channels:    0,
            wpc:           3,
            nodes:         Vec::new(),
            branches:      Vec::new(),
            c_cached:      Vec::new(),
            s_cached:      Vec::new(),
            c_j_cached:    1.375e-16,
            alpha_eff_neper_m: 29.78,
        }
    }
}

impl Device for NativePnThermalPhaseShifterCap {
    fn num_terminals(&self) -> usize { self.nodes.len() }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
        self.alpha_eff_neper_m = self.alpha_neper_m;
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel(); self.wpc = wpc;
        let stride = 2 * wpc;
        assert!(
            terminals.len() >= stride + 4 && (terminals.len() - 4) % stride == 0,
            "fc_pn_th_ps_cap: terminal count must be {stride}·N + 4 (wpc={wpc}); got {}",
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
            "l_um"   => { self.length_m = value * 1e-6; true }
            "l_m" | "length" => { self.length_m = value; true }
            "n_g"    => { self.n_g = value; true }
            "n_eff"  => { self.n_eff = value; true }
            "wl_ref_m" | "lambda_ref_m" => { self.wl_ref_m = value; true }
            "wl_ref_nm" | "lambda_ref_nm" => { self.wl_ref_m = value * 1e-9; true }
            "dn_dv"  => { self.dn_dv = value; true }
            "g_pn"   => { self.g_pn = value; true }
            "v_pi_l" => {
                if value > 0.0 { self.dn_dv = self.wl_ref_m / (2.0 * value); }
                true
            }
            "alpha_db_cm" => { self.alpha_neper_m = dB_per_cm_to_neper_per_m(value); true }
            "c_j0"   => { self.c_j0 = value.max(0.0); true }
            "v_bi"   => { self.v_bi = value.max(1e-3); true }
            "m_j"    => { self.m_j  = value.clamp(0.0, 0.99); true }
            "da_dv"  => { self.da_dv = value; true }
            "r_heater" | "r" => { self.r_heater = value; true }
            "p_pi" | "p_pi_w" | "p_pi_th" => { self.p_pi_th = value; true }
            "pin_at_ref" => { self.pin_at_ref = value != 0.0; true }
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

        // C_j(V) with linear continuation past the V_bi/2 knee.
        let v_knee = 0.5 * self.v_bi;
        self.c_j_cached = if v_pn < v_knee {
            self.c_j0 / (1.0 - v_pn / self.v_bi).powf(self.m_j)
        } else {
            let c_knee = self.c_j0 / (1.0 - v_knee / self.v_bi).powf(self.m_j);
            let dc_dv  = c_knee * self.m_j / (self.v_bi - v_knee);
            c_knee + dc_dv * (v_pn - v_knee)
        };
        // Bias-dependent loss: α(V) = α_0 + da_dv·max(0, -V).
        let v_rev = (-v_pn).max(0.0);
        self.alpha_eff_neper_m = self.alpha_neper_m + self.da_dv * v_rev;
        let t_amp  = (-self.alpha_eff_neper_m * self.length_m / 2.0).exp();

        let two_pi = 2.0 * std::f64::consts::PI;
        let lam = wpc - 1;
        let phi_ref = if self.pin_at_ref {
            two_pi * self.n_eff * self.length_m / self.wl_ref_m
        } else { 0.0 };
        for k in 0..n {
            let lambda = match self.nodes[wpc * k + lam] {
                Some(i) => { let v = x[i]; if v.abs() > 1e-9 { v } else { self.wl_ref_m } }
                None    => self.wl_ref_m,
            };
            let n_eff_lam = n_eff_at_lambda(self.n_eff, self.n_g, self.wl_ref_m, lambda);
            let phi_abs  = two_pi * n_eff_lam * self.length_m / lambda;
            let phi_prop = phi_abs - phi_ref;
            let phi_eo   = two_pi * self.length_m * self.dn_dv * v_pn / lambda;
            let phi      = phi_prop + phi_eo + phi_th;
            self.c_cached[k] = t_amp * phi.cos();
            self.s_cached[k] = t_amp * phi.sin();
        }
    }

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        stamp_pn_ths_jacobian(
            mat, &self.nodes, &self.branches, self.n_channels, self.wpc,
            self.g_pn, self.r_heater, &self.c_cached, &self.s_cached,
        );
    }

    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> {
        let elec_base = 2 * self.wpc * self.n_channels;
        let anode = self.nodes.get(elec_base    ).copied().flatten();
        let cath  = self.nodes.get(elec_base + 1).copied().flatten();
        vec![ReactiveBranchSpec {
            kind: ReactiveKind::Capacitor, pos: anode, neg: cath,
            value: self.c_j_cached,
        }]
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) { self.load_residual(b); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) { self.load_jacobian(mat); }
}

// ────────────────────────────────────────────────────────────────────────
// Native PN phase shifter, forward-bias injection (L2b)
// ────────────────────────────────────────────────────────────────────────

/// Forward-bias carrier-injection PN phase shifter.  Suitable for V_pn ∈
/// [0, ~0.8 V] only — does not model depletion physics.  Use this when the
/// design relies on injected free-carrier index change (typical VOA / slow
/// thermal-trim modulator with carrier dynamics on top).
///
/// Physics differences from fc_pn_ps_cap:
///   - Shockley forward I-V:   I = I_s · (exp(V/(n·V_T)) − 1), linearised at
///     the operating point with g_d = (I + I_s)/(n·V_T).
///   - Diffusion capacitance:  C_d = τ_carrier · g_d  (replaces depletion).
///   - Exponential injection Δn(V):  Δn_inj = K_inj · (exp(V/(n·V_T)) − 1),
///     where K_inj is the linearised forward-bias dn_dv coefficient scaled
///     to give the requested fractional dn at V = ~3·V_T.
///   - Exponential injection loss Δα(V):  same exp shape with K_alpha.
///
/// Terminal layout: 2·wpc·N + 2  (same as `fc_pn_ps`).
pub struct NativePnPhaseShifterInj {
    length_m:      f64,
    n_eff:         f64,
    n_g:           f64,
    wl_ref_m:      f64,
    pin_at_ref:    bool,
    // Forward-bias diode
    i_sat:         f64,
    n_diode:       f64,
    tau_carrier:   f64,
    // EO + FCA injection coefficients
    dn_dv_inj:     f64,    // K_inj = Δn_eff(V→Vt) / 1 (slope-equivalent)
    da_dv_inj:     f64,    // Np/m per (exp(V/Vt)-1)
    alpha_neper_m: f64,    // background propagation loss
    n_channels:    usize,
    wpc:           usize,
    nodes:         Vec<NodeId>,
    branches:      Vec<Option<usize>>,
    c_cached:      Vec<f64>,
    s_cached:      Vec<f64>,
    g_d_cached:    f64,
    i_eq_cached:   f64,
    c_d_cached:    f64,
    v_pn_op:       f64,
}

impl NativePnPhaseShifterInj {
    pub fn new() -> Self {
        Self {
            length_m:      1e-3,
            n_eff:         2.7654,
            n_g:           4.02,
            wl_ref_m:      1.55e-6,
            pin_at_ref:    false,
            i_sat:         1e-12,
            n_diode:       1.05,
            tau_carrier:   10e-9,
            dn_dv_inj:     1.311e-4,                       // forward small-signal coeff (sim)
            da_dv_inj:     150.0,                           // exp prefactor for FCA injection (Np/m)
            alpha_neper_m: dB_per_cm_to_neper_per_m(1.0),
            n_channels:    0,
            wpc:           3,
            nodes:         Vec::new(),
            branches:      Vec::new(),
            c_cached:      Vec::new(),
            s_cached:      Vec::new(),
            g_d_cached:    1e-9,
            i_eq_cached:   0.0,
            c_d_cached:    0.0,
            v_pn_op:       0.0,
        }
    }
}

impl Device for NativePnPhaseShifterInj {
    fn num_terminals(&self) -> usize { self.nodes.len() }
    fn setup_model(&mut self, ctx: &SimContext) { self.wpc = ctx.wires_per_channel(); }
    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel(); self.wpc = wpc;
        let stride = 2 * wpc;
        assert!(
            terminals.len() >= stride + 2 && (terminals.len() - 2) % stride == 0,
            "fc_pn_ps_inj: terminal count must be {stride}·N + 2 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 2) / stride;
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
            "l_um"   => { self.length_m = value * 1e-6; true }
            "l_m" | "length" => { self.length_m = value; true }
            "n_g"    => { self.n_g = value; true }
            "n_eff"  => { self.n_eff = value; true }
            "wl_ref_nm" => { self.wl_ref_m = value * 1e-9; true }
            "wl_ref_m" => { self.wl_ref_m = value; true }
            "i_sat" | "is" => { self.i_sat = value.max(0.0); true }
            "n_diode" | "n" => { self.n_diode = value.max(0.5); true }
            "tau_carrier" | "tau" => { self.tau_carrier = value.max(0.0); true }
            "dn_dv_inj" | "dn_dv" => { self.dn_dv_inj = value; true }
            "da_dv_inj" | "da_dv" => { self.da_dv_inj = value; true }
            "alpha_db_cm" => { self.alpha_neper_m = dB_per_cm_to_neper_per_m(value); true }
            "pin_at_ref"  => { self.pin_at_ref = value != 0.0; true }
            _ => false,
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, ctx: &SimContext) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        let elec = 2 * wpc * n;
        let v_a = self.nodes[elec    ].map_or(0.0, |i| x[i]);
        let v_c = self.nodes[elec + 1].map_or(0.0, |i| x[i]);
        let v_pn = v_a - v_c;
        self.v_pn_op = v_pn;

        // Shockley I-V (clamp exp argument for NR stability)
        let vt = ctx.vt() * self.n_diode;
        let arg = (v_pn / vt).min(40.0).max(-40.0);
        let e = arg.exp();
        let i_diode = self.i_sat * (e - 1.0);
        self.g_d_cached = self.i_sat * e / vt;
        // Norton equivalent: I_eq = I(V_op) - g · V_op
        self.i_eq_cached = i_diode - self.g_d_cached * v_pn;
        // Diffusion capacitance (forward bias only — small under reverse)
        self.c_d_cached = self.tau_carrier * self.g_d_cached;

        // Optical: injection-driven Δn and Δα, only forward bias contributes.
        let inj = (e - 1.0).max(0.0);
        let alpha_eff = self.alpha_neper_m + self.da_dv_inj * inj;
        let t_amp = (-alpha_eff * self.length_m / 2.0).exp();

        let two_pi = 2.0 * std::f64::consts::PI;
        let lam = wpc - 1;
        let phi_ref = if self.pin_at_ref {
            two_pi * self.n_eff * self.length_m / self.wl_ref_m
        } else { 0.0 };
        for k in 0..n {
            let lambda = match self.nodes[wpc * k + lam] {
                Some(i) => { let v = x[i]; if v.abs() > 1e-9 { v } else { self.wl_ref_m } }
                None    => self.wl_ref_m,
            };
            let n_eff_lam = n_eff_at_lambda(self.n_eff, self.n_g, self.wl_ref_m, lambda);
            let phi_abs  = two_pi * n_eff_lam * self.length_m / lambda;
            let phi_prop = phi_abs - phi_ref;
            // For injection, we add a Soref-Bennett-shaped Δn (negative for
            // more carriers → less n).  Linear coefficient `dn_dv_inj` is the
            // slope at V=0; full exponential form scales as (e-1)/1V.
            let phi_eo = -two_pi * self.length_m * self.dn_dv_inj * inj / lambda;
            let phi    = phi_prop + phi_eo;
            self.c_cached[k] = t_amp * phi.cos();
            self.s_cached[k] = t_amp * phi.sin();
        }
    }

    fn load_residual(&self, b: &mut [f64]) {
        let elec = 2 * self.wpc * self.n_channels;
        if let Some(a) = self.nodes[elec]     { b[a] -= self.i_eq_cached; }
        if let Some(c) = self.nodes[elec + 1] { b[c] += self.i_eq_cached; }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        stamp_pn_optical(
            mat, &self.nodes, &self.branches, self.n_channels, self.wpc,
            &self.c_cached, &self.s_cached,
        );
        // Electrical: shared diode small-signal conductance
        let elec = 2 * self.wpc * self.n_channels;
        stamp_resistor(mat, self.nodes[elec], self.nodes[elec + 1], self.g_d_cached);
    }

    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> {
        let elec = 2 * self.wpc * self.n_channels;
        let a = self.nodes.get(elec).copied().flatten();
        let c = self.nodes.get(elec + 1).copied().flatten();
        vec![ReactiveBranchSpec {
            kind: ReactiveKind::Capacitor, pos: a, neg: c,
            value: self.c_d_cached,
        }]
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) { self.load_residual(b); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) { self.load_jacobian(mat); }
}

// ────────────────────────────────────────────────────────────────────────
// Native PN+thermal phase shifter, forward injection (L2b-with-heater)
// ────────────────────────────────────────────────────────────────────────

pub struct NativePnThermalPhaseShifterInj {
    inj: NativePnPhaseShifterInj,
    r_heater: f64,
    p_pi_th:  f64,
}

impl NativePnThermalPhaseShifterInj {
    pub fn new() -> Self {
        Self { inj: NativePnPhaseShifterInj::new(),
               r_heater: 1000.0, p_pi_th: 10e-3 }
    }
}

impl Device for NativePnThermalPhaseShifterInj {
    fn num_terminals(&self) -> usize { self.inj.nodes.len() }
    fn setup_model(&mut self, ctx: &SimContext) { self.inj.setup_model(ctx); }
    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        // Same as `_inj` but expects 4 extra electrical pins.
        let wpc = ctx.wires_per_channel(); self.inj.wpc = wpc;
        let stride = 2 * wpc;
        assert!(
            terminals.len() >= stride + 4 && (terminals.len() - 4) % stride == 0,
            "fc_pn_th_ps_inj: terminal count must be {stride}·N + 4 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 4) / stride;
        self.inj.n_channels = n;
        self.inj.nodes      = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 };
        self.inj.branches   = vec![None; bpc * n];
        self.inj.c_cached   = vec![1.0; n];
        self.inj.s_cached   = vec![0.0; n];
    }
    fn num_extra_nodes(&self) -> usize { self.inj.branches.len() }
    fn bind_extra_nodes(&mut self, idx: usize) { self.inj.bind_extra_nodes(idx); }
    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "r_heater" | "r" => { self.r_heater = value; true }
            "p_pi" | "p_pi_w" | "p_pi_th" => { self.p_pi_th = value; true }
            _ => self.inj.set_real_param(name, value),
        }
    }
    fn eval(&mut self, x: &[f64], flags: EvalFlags, ctx: &SimContext) {
        self.inj.eval(x, flags, ctx);
        // Add thermal phase from heater (heat_p, heat_n live at the END of nodes).
        let n   = self.inj.n_channels;
        let wpc = self.inj.wpc;
        let elec = 2 * wpc * n;
        let v_hp = self.inj.nodes[elec + 2].map_or(0.0, |i| x[i]);
        let v_hn = self.inj.nodes[elec + 3].map_or(0.0, |i| x[i]);
        let v_h  = v_hp - v_hn;
        let phi_th = std::f64::consts::PI * (v_h * v_h / self.r_heater) / self.p_pi_th;
        // Rotate the cached c, s by phi_th (no re-eval of EO part needed):
        let cth = phi_th.cos(); let sth = phi_th.sin();
        for k in 0..n {
            let c = self.inj.c_cached[k]; let s = self.inj.s_cached[k];
            self.inj.c_cached[k] = c*cth - s*sth;
            self.inj.s_cached[k] = c*sth + s*cth;
        }
    }
    fn load_residual(&self, b: &mut [f64]) { self.inj.load_residual(b); }
    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        self.inj.load_jacobian(mat);
        // Heater resistor between heat_p and heat_n.
        let elec = 2 * self.inj.wpc * self.inj.n_channels;
        stamp_resistor(mat, self.inj.nodes[elec + 2], self.inj.nodes[elec + 3],
                       1.0 / self.r_heater);
    }
    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> { self.inj.reactive_branches() }
    fn load_residual_tran(&self, b: &mut [f64], a: f64) { self.inj.load_residual_tran(b, a); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, a: f64) {
        self.inj.load_jacobian_tran(mat, a);
        let elec = 2 * self.inj.wpc * self.inj.n_channels;
        stamp_resistor(mat, self.inj.nodes[elec + 2], self.inj.nodes[elec + 3],
                       1.0 / self.r_heater);
    }
}

// ────────────────────────────────────────────────────────────────────────
// Native PN phase shifter, combined depletion + injection + TPA + static
// self-heating (L3).
// ────────────────────────────────────────────────────────────────────────

/// "Full" PN modulator — captures both bias regimes, two-photon absorption,
/// and static thermal self-heating from absorbed optical power.  Intended
/// for steady-state and slow transient analysis of high-performance
/// modulators; carrier and thermal dynamics belong in `_carrier` (L4).
///
/// Physics included:
///   - Reverse-bias depletion Δn_rev(V) = dn_dv_rev · V  (V ≤ 0)
///   - Forward-bias injection Δn_fwd(V) = -dn_dv_inj · (exp(V/V_T)-1)  (V ≥ 0)
///   - Depletion C_j(V) for V ≤ 0; diffusion C_d for V > 0
///   - Reverse FCA loss α_rev(V) = da_dv_rev · max(0,-V)
///   - Forward FCA loss α_fwd(V) = da_dv_inj · (exp(V/V_T)-1)
///   - TPA loss α_TPA = β_TPA · (|A|²/A_eff)
///   - Static self-heating: ΔT_ss = R_th · α_total · L · |A|²;  Δn_th = dn_dT · ΔT
pub struct NativePnPhaseShifterFull {
    length_m:      f64,
    n_eff:         f64,
    n_g:           f64,
    wl_ref_m:      f64,
    pin_at_ref:    bool,
    // Reverse
    dn_dv_rev:     f64,
    da_dv_rev:     f64,
    c_j0:          f64,
    v_bi:          f64,
    m_j:           f64,
    // Forward
    i_sat:         f64,
    n_diode:       f64,
    tau_carrier:   f64,
    dn_dv_inj:     f64,
    da_dv_inj:     f64,
    // Common
    alpha_neper_m: f64,
    // TPA + thermal
    beta_tpa_m_per_w: f64,
    a_eff_m2:      f64,
    r_th_k_per_w:  f64,
    dn_dt:         f64,
    n_channels:    usize,
    wpc:           usize,
    nodes:         Vec<NodeId>,
    branches:      Vec<Option<usize>>,
    c_cached:      Vec<f64>,
    s_cached:      Vec<f64>,
    g_pn_cached:   f64,
    i_eq_cached:   f64,
    c_eff_cached:  f64,
}

impl NativePnPhaseShifterFull {
    pub fn new() -> Self {
        Self {
            length_m:      1e-3,
            n_eff:         2.7654,
            n_g:           4.02,
            wl_ref_m:      1.55e-6,
            pin_at_ref:    false,
            dn_dv_rev:     5.024e-5,
            da_dv_rev:     7.83,
            c_j0:          1.375e-13,                       // F at default L
            v_bi:          0.917,
            m_j:           0.5,
            i_sat:         1e-12,
            n_diode:       1.05,
            tau_carrier:   10e-9,
            dn_dv_inj:     1.311e-4,
            da_dv_inj:     150.0,
            alpha_neper_m: dB_per_cm_to_neper_per_m(1.0),
            beta_tpa_m_per_w: 7.9e-12,
            a_eff_m2:      1.257e-13,
            r_th_k_per_w:  0.0,                              // user must set
            dn_dt:         1.86e-4,                          // crystalline Si
            n_channels:    0,
            wpc:           3,
            nodes:         Vec::new(),
            branches:      Vec::new(),
            c_cached:      Vec::new(),
            s_cached:      Vec::new(),
            g_pn_cached:   1e-9,
            i_eq_cached:   0.0,
            c_eff_cached:  1.375e-13,
        }
    }
}

impl Device for NativePnPhaseShifterFull {
    fn num_terminals(&self) -> usize { self.nodes.len() }
    fn setup_model(&mut self, ctx: &SimContext) { self.wpc = ctx.wires_per_channel(); }
    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel(); self.wpc = wpc;
        let stride = 2 * wpc;
        assert!(
            terminals.len() >= stride + 2 && (terminals.len() - 2) % stride == 0,
            "fc_pn_ps_full: terminal count must be {stride}·N + 2 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 2) / stride;
        self.n_channels = n;
        self.nodes      = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 };
        self.branches   = vec![None; bpc * n];
        self.c_cached   = vec![1.0; n];
        self.s_cached   = vec![0.0; n];
    }
    fn num_extra_nodes(&self) -> usize { self.branches.len() }
    fn bind_extra_nodes(&mut self, idx: usize) {
        for i in 0..self.branches.len() { self.branches[i] = Some(idx + i); }
    }
    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "l_um"   => { self.length_m = value * 1e-6; true }
            "l_m" | "length" => { self.length_m = value; true }
            "n_g" => { self.n_g = value; true }
            "n_eff" => { self.n_eff = value; true }
            "wl_ref_nm" => { self.wl_ref_m = value * 1e-9; true }
            "wl_ref_m" => { self.wl_ref_m = value; true }
            "dn_dv_rev" | "dn_dv" => { self.dn_dv_rev = value; true }
            "da_dv_rev" | "da_dv" => { self.da_dv_rev = value; true }
            "c_j0" => { self.c_j0 = value.max(0.0); true }
            "v_bi" => { self.v_bi = value.max(1e-3); true }
            "m_j"  => { self.m_j  = value.clamp(0.0, 0.99); true }
            "i_sat" | "is" => { self.i_sat = value.max(0.0); true }
            "n_diode" | "n" => { self.n_diode = value.max(0.5); true }
            "tau_carrier" | "tau" => { self.tau_carrier = value.max(0.0); true }
            "dn_dv_inj" => { self.dn_dv_inj = value; true }
            "da_dv_inj" => { self.da_dv_inj = value; true }
            "alpha_db_cm" => { self.alpha_neper_m = dB_per_cm_to_neper_per_m(value); true }
            "beta_tpa" | "beta_tpa_m_per_w" => { self.beta_tpa_m_per_w = value; true }
            "a_eff_m2" => { self.a_eff_m2 = value.max(1e-20); true }
            "a_eff_um2" => { self.a_eff_m2 = value.max(1e-8) * 1e-12; true }
            "r_th" | "r_th_k_per_w" => { self.r_th_k_per_w = value.max(0.0); true }
            "dn_dt" => { self.dn_dt = value; true }
            "pin_at_ref"  => { self.pin_at_ref = value != 0.0; true }
            _ => false,
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, ctx: &SimContext) {
        let n   = self.n_channels;
        let wpc = self.wpc;
        let elec = 2 * wpc * n;
        let v_a = self.nodes[elec    ].map_or(0.0, |i| x[i]);
        let v_c = self.nodes[elec + 1].map_or(0.0, |i| x[i]);
        let v_pn = v_a - v_c;

        // Electrical: Shockley I-V across the whole regime.
        let vt = ctx.vt() * self.n_diode;
        let e = (v_pn / vt).min(40.0).max(-40.0).exp();
        let i_diode = self.i_sat * (e - 1.0);
        let g_d = self.i_sat * e / vt;
        self.g_pn_cached = g_d.max(1e-15);
        self.i_eq_cached = i_diode - g_d * v_pn;

        // Capacitance: piecewise C_j (reverse) vs C_d (forward).
        let c_j_v = {
            let v_knee = 0.5 * self.v_bi;
            if v_pn < v_knee {
                self.c_j0 / (1.0 - v_pn / self.v_bi).powf(self.m_j)
            } else {
                let c_knee = self.c_j0 / (1.0 - v_knee / self.v_bi).powf(self.m_j);
                let dc_dv  = c_knee * self.m_j / (self.v_bi - v_knee);
                c_knee + dc_dv * (v_pn - v_knee)
            }
        };
        let c_d_v = self.tau_carrier * g_d;
        self.c_eff_cached = c_j_v + c_d_v;

        // Optical: sum of contributions.  Compute |A|² at the input for TPA
        // and self-heating (use channel-0; aggregating across WDM is left to L4).
        let v_re_0 = self.nodes[0].map_or(0.0, |i| x[i]);
        let v_im_0 = self.nodes[1].map_or(0.0, |i| x[i]);
        let intensity_w = (v_re_0*v_re_0 + v_im_0*v_im_0).max(0.0);
        let alpha_tpa = self.beta_tpa_m_per_w * intensity_w / self.a_eff_m2;
        let inj = (e - 1.0).max(0.0);
        let v_rev = (-v_pn).max(0.0);
        let alpha_fca = self.alpha_neper_m + self.da_dv_rev * v_rev
                      + self.da_dv_inj * inj;
        let alpha_total = alpha_fca + alpha_tpa;
        let t_amp = (-alpha_total * self.length_m / 2.0).exp();

        // Self-heating Δn (static): ΔT = R_th · P_abs;  P_abs ≈ α · L · |A|².
        let p_abs = alpha_total * self.length_m * intensity_w;
        let dn_self = self.dn_dt * self.r_th_k_per_w * p_abs;

        let two_pi = 2.0 * std::f64::consts::PI;
        let lam = wpc - 1;
        let phi_ref = if self.pin_at_ref {
            two_pi * self.n_eff * self.length_m / self.wl_ref_m
        } else { 0.0 };
        for k in 0..n {
            let lambda = match self.nodes[wpc * k + lam] {
                Some(i) => { let v = x[i]; if v.abs() > 1e-9 { v } else { self.wl_ref_m } }
                None    => self.wl_ref_m,
            };
            let n_eff_lam = n_eff_at_lambda(self.n_eff, self.n_g, self.wl_ref_m, lambda);
            let phi_abs  = two_pi * n_eff_lam * self.length_m / lambda;
            let phi_prop = phi_abs - phi_ref;
            let phi_eo_rev = two_pi * self.length_m * self.dn_dv_rev * v_pn / lambda;
            let phi_eo_inj = -two_pi * self.length_m * self.dn_dv_inj * inj / lambda;
            let phi_self   = two_pi * self.length_m * dn_self / lambda;
            let phi = phi_prop + phi_eo_rev + phi_eo_inj + phi_self;
            self.c_cached[k] = t_amp * phi.cos();
            self.s_cached[k] = t_amp * phi.sin();
        }
    }

    fn load_residual(&self, b: &mut [f64]) {
        let elec = 2 * self.wpc * self.n_channels;
        if let Some(a) = self.nodes[elec]     { b[a] -= self.i_eq_cached; }
        if let Some(c) = self.nodes[elec + 1] { b[c] += self.i_eq_cached; }
    }
    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        stamp_pn_optical(
            mat, &self.nodes, &self.branches, self.n_channels, self.wpc,
            &self.c_cached, &self.s_cached,
        );
        let elec = 2 * self.wpc * self.n_channels;
        stamp_resistor(mat, self.nodes[elec], self.nodes[elec + 1], self.g_pn_cached);
    }
    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> {
        let elec = 2 * self.wpc * self.n_channels;
        let a = self.nodes.get(elec).copied().flatten();
        let c = self.nodes.get(elec + 1).copied().flatten();
        vec![ReactiveBranchSpec { kind: ReactiveKind::Capacitor, pos: a, neg: c,
            value: self.c_eff_cached }]
    }
    fn load_residual_tran(&self, b: &mut [f64], _a: f64) { self.load_residual(b); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _a: f64) { self.load_jacobian(mat); }
}

// ────────────────────────────────────────────────────────────────────────
// Native PN+thermal phase shifter, full (L3-with-heater)
// ────────────────────────────────────────────────────────────────────────

pub struct NativePnThermalPhaseShifterFull {
    full: NativePnPhaseShifterFull,
    r_heater: f64,
    p_pi_th:  f64,
}

impl NativePnThermalPhaseShifterFull {
    pub fn new() -> Self {
        Self { full: NativePnPhaseShifterFull::new(),
               r_heater: 1000.0, p_pi_th: 10e-3 }
    }
}

impl Device for NativePnThermalPhaseShifterFull {
    fn num_terminals(&self) -> usize { self.full.nodes.len() }
    fn setup_model(&mut self, ctx: &SimContext) { self.full.setup_model(ctx); }
    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel(); self.full.wpc = wpc;
        let stride = 2 * wpc;
        assert!(
            terminals.len() >= stride + 4 && (terminals.len() - 4) % stride == 0,
            "fc_pn_th_ps_full: terminal count must be {stride}·N + 4 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 4) / stride;
        self.full.n_channels = n;
        self.full.nodes      = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 };
        self.full.branches   = vec![None; bpc * n];
        self.full.c_cached   = vec![1.0; n];
        self.full.s_cached   = vec![0.0; n];
    }
    fn num_extra_nodes(&self) -> usize { self.full.branches.len() }
    fn bind_extra_nodes(&mut self, idx: usize) { self.full.bind_extra_nodes(idx); }
    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "r_heater" | "r" => { self.r_heater = value; true }
            "p_pi" | "p_pi_w" | "p_pi_th" => { self.p_pi_th = value; true }
            _ => self.full.set_real_param(name, value),
        }
    }
    fn eval(&mut self, x: &[f64], flags: EvalFlags, ctx: &SimContext) {
        self.full.eval(x, flags, ctx);
        let n   = self.full.n_channels;
        let wpc = self.full.wpc;
        let elec = 2 * wpc * n;
        let v_hp = self.full.nodes[elec + 2].map_or(0.0, |i| x[i]);
        let v_hn = self.full.nodes[elec + 3].map_or(0.0, |i| x[i]);
        let v_h  = v_hp - v_hn;
        let phi_th = std::f64::consts::PI * (v_h * v_h / self.r_heater) / self.p_pi_th;
        let cth = phi_th.cos(); let sth = phi_th.sin();
        for k in 0..n {
            let c = self.full.c_cached[k]; let s = self.full.s_cached[k];
            self.full.c_cached[k] = c*cth - s*sth;
            self.full.s_cached[k] = c*sth + s*cth;
        }
    }
    fn load_residual(&self, b: &mut [f64]) { self.full.load_residual(b); }
    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        self.full.load_jacobian(mat);
        let elec = 2 * self.full.wpc * self.full.n_channels;
        stamp_resistor(mat, self.full.nodes[elec + 2], self.full.nodes[elec + 3],
                       1.0 / self.r_heater);
    }
    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> { self.full.reactive_branches() }
    fn load_residual_tran(&self, b: &mut [f64], a: f64) { self.full.load_residual_tran(b, a); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, a: f64) {
        self.full.load_jacobian_tran(mat, a);
        let elec = 2 * self.full.wpc * self.full.n_channels;
        stamp_resistor(mat, self.full.nodes[elec + 2], self.full.nodes[elec + 3],
                       1.0 / self.r_heater);
    }
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

// ────────────────────────────────────────────────────────────────────────
// Native 3-port circulator (bidir-only)
// ────────────────────────────────────────────────────────────────────────

/// 3-port circulator.  Routes light cyclically: light entering port 1
/// exits port 2; entering port 2 exits port 3; entering port 3 exits
/// port 1.  Requires `enable_bidirectional=1` because the routing
/// fundamentally needs each port to support both an incoming wave (fw,
/// inward to the circulator) and an outgoing wave (bw, outward from the
/// circulator).  Errors at setup_instance if bidir is off.
///
/// Wire convention (consistent across the circulator): at every port,
/// `re_fw`/`im_fw` represent the wave flowing INWARD (toward the device)
/// and `re_bw`/`im_bw` represent the wave flowing OUTWARD.  Internal
/// routing:
///   port_p.bw = port_((p+2) mod 3).fw   — for re and im, every channel
/// (light entering port (p-1) leaves at port p, mod 3).
///
/// λ is tied across all three ports: `port_1.λ = port_0.λ`, `port_2.λ =
/// port_0.λ`.  This works whether the laser drives port 0, 1, or 2 —
/// SPICE branch equations resolve the cycle consistently.
///
/// Bundle-aware: 3·wpc·N terminals for N WDM channels.  Per channel
/// branch count: 6 re/im routing + 2 λ ties = 8.
pub struct NativeCirculator {
    n_channels: usize,
    wpc:        usize,
    nodes:      Vec<NodeId>,
    branches:   Vec<Option<usize>>,
}

impl NativeCirculator {
    pub fn new() -> Self {
        Self {
            n_channels: 0,
            wpc:        5,
            nodes:      Vec::new(),
            branches:   Vec::new(),
        }
    }
}

impl Device for NativeCirculator {
    fn num_terminals(&self) -> usize { self.nodes.len() }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        if ctx.wires_per_channel() != 5 {
            panic!("fc_circulator requires bidirectional propagation; \
                    set `.options enable_bidirectional=1` (or via CLI / Python)");
        }
        self.wpc = 5;
        let stride = 3 * 5; // 3 ports × 5 wires per channel
        assert!(
            !terminals.is_empty() && terminals.len() % stride == 0,
            "fc_circulator: terminal count must be {stride}·N for N ≥ 1 channels; got {}",
            terminals.len()
        );
        let n = terminals.len() / stride;
        self.n_channels = n;
        self.nodes      = terminals.to_vec();
        self.branches   = vec![None; 8 * n];
    }

    fn num_extra_nodes(&self) -> usize { self.branches.len() }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..self.branches.len() { self.branches[i] = Some(first_idx + i); }
    }

    fn set_real_param(&mut self, _name: &str, _value: f64) -> bool { false }

    fn eval(&mut self, _x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {}

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n   = self.n_channels;
        let wpc = 5;
        // Per channel: stride = 3 ports × 5 wires = 15.
        for k in 0..n {
            let base = 15 * k;
            // Wires per port: [re_fw, im_fw, re_bw, im_bw, λ].
            let port_wires = |p: usize| -> (NodeId, NodeId, NodeId, NodeId, NodeId) {
                let pb = base + wpc * p;
                (
                    self.nodes[pb],     // re_fw
                    self.nodes[pb + 1], // im_fw
                    self.nodes[pb + 2], // re_bw
                    self.nodes[pb + 3], // im_bw
                    self.nodes[pb + 4], // λ
                )
            };
            let (p0_re_fw, p0_im_fw, p0_re_bw, p0_im_bw, p0_lam) = port_wires(0);
            let (p1_re_fw, p1_im_fw, p1_re_bw, p1_im_bw, p1_lam) = port_wires(1);
            let (p2_re_fw, p2_im_fw, p2_re_bw, p2_im_bw, p2_lam) = port_wires(2);
            let b = 8 * k;
            // port_p.bw = port_((p+2) mod 3).fw
            // port_0.bw = port_2.fw
            stamp_potential_eq(mat, &self.branches, b,     p0_re_bw, &[(p2_re_fw, -1.0)]);
            stamp_potential_eq(mat, &self.branches, b + 1, p0_im_bw, &[(p2_im_fw, -1.0)]);
            // port_1.bw = port_0.fw
            stamp_potential_eq(mat, &self.branches, b + 2, p1_re_bw, &[(p0_re_fw, -1.0)]);
            stamp_potential_eq(mat, &self.branches, b + 3, p1_im_bw, &[(p0_im_fw, -1.0)]);
            // port_2.bw = port_1.fw
            stamp_potential_eq(mat, &self.branches, b + 4, p2_re_bw, &[(p1_re_fw, -1.0)]);
            stamp_potential_eq(mat, &self.branches, b + 5, p2_im_bw, &[(p1_im_fw, -1.0)]);
            // λ ties: port_1.λ = port_0.λ, port_2.λ = port_0.λ.
            stamp_potential_eq(mat, &self.branches, b + 6, p1_lam, &[(p0_lam, -1.0)]);
            stamp_potential_eq(mat, &self.branches, b + 7, p2_lam, &[(p0_lam, -1.0)]);
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

/// Stamp a linear resistor of conductance `g` between two MNA nodes.  No-op
/// for ground terminals.  Symmetric standard 2×2 stamp.
#[inline]
fn stamp_resistor(mat: &mut MnaMatrix, a: NodeId, b: NodeId, g: f64) {
    if let Some(i) = a {
        mat.a[i][i] += g;
        if let Some(j) = b { mat.a[i][j] -= g; }
    }
    if let Some(j) = b {
        mat.a[j][j] += g;
        if let Some(i) = a { mat.a[j][i] -= g; }
    }
}

/// Stamp the per-channel optical-branch equations for a PN-style phase
/// shifter that has already computed `c_cached[k]`, `s_cached[k]` per
/// channel.  Layout assumed: input bundle (wpc·n), output bundle (wpc·n).
/// Used by `fc_pn_ps_inj` and `fc_pn_ps_full` whose electrical side is
/// stamped separately.
fn stamp_pn_optical(
    mat: &mut MnaMatrix, nodes: &[NodeId], branches: &[Option<usize>],
    n: usize, wpc: usize, c_cached: &[f64], s_cached: &[f64],
) {
    let bpc = if wpc == 5 { 5 } else { 3 };
    let lam = wpc - 1;
    let out_base = wpc * n;
    for k in 0..n {
        let in_re_fw  = nodes[wpc * k];
        let in_im_fw  = nodes[wpc * k + 1];
        let in_l      = nodes[wpc * k + lam];
        let out_re_fw = nodes[out_base + wpc * k];
        let out_im_fw = nodes[out_base + wpc * k + 1];
        let out_l     = nodes[out_base + wpc * k + lam];
        let c = c_cached[k]; let s = s_cached[k];
        stamp_potential_eq(mat, branches, bpc * k,     out_re_fw,
            &[(in_re_fw, -c), (in_im_fw, -s)]);
        stamp_potential_eq(mat, branches, bpc * k + 1, out_im_fw,
            &[(in_re_fw,  s), (in_im_fw, -c)]);
        stamp_potential_eq(mat, branches, bpc * k + (bpc - 1), out_l,
            &[(in_l, -1.0)]);
        if wpc == 5 {
            let in_re_bw  = nodes[wpc * k + 2];
            let in_im_bw  = nodes[wpc * k + 3];
            let out_re_bw = nodes[out_base + wpc * k + 2];
            let out_im_bw = nodes[out_base + wpc * k + 3];
            stamp_potential_eq(mat, branches, bpc * k + 2, in_re_bw,
                &[(out_re_bw, -c), (out_im_bw, -s)]);
            stamp_potential_eq(mat, branches, bpc * k + 3, in_im_bw,
                &[(out_re_bw,  s), (out_im_bw, -c)]);
        }
    }
}

/// Full Jacobian stamp for an fc_pn_th_ps_cap-style device (PN + heater +
/// per-channel optics).  Stamps PN g_pn, heater g_h, and the optical
/// branches.  Skips the electrical residual entry (caller stamps if needed).
fn stamp_pn_ths_jacobian(
    mat: &mut MnaMatrix, nodes: &[NodeId], branches: &[Option<usize>],
    n: usize, wpc: usize, g_pn: f64, r_heater: f64,
    c_cached: &[f64], s_cached: &[f64],
) {
    let elec = 2 * wpc * n;
    stamp_resistor(mat, nodes[elec],     nodes[elec + 1], g_pn);
    stamp_resistor(mat, nodes[elec + 2], nodes[elec + 3], 1.0 / r_heater);
    stamp_pn_optical(mat, nodes, branches, n, wpc, c_cached, s_cached);
}

#[allow(non_snake_case)]
fn dB_per_cm_to_neper_per_m(alpha_db_cm: f64) -> f64 {
    // 1 dB = ln(10)/20 Np ≈ 0.1151 Np; 1 cm = 0.01 m → multiply by 100/cm.
    alpha_db_cm * 100.0 * std::f64::consts::LN_10 / 20.0
}

/// Speed of light in vacuum (m/s).
pub const C0: f64 = 299_792_458.0;

/// First-order dispersion-corrected effective index at wavelength `lambda`,
/// given `n_eff_0` and `n_g_0` evaluated at reference wavelength `wl_ref_m`.
///
/// Physics: n_g(λ) = n_eff(λ) − λ·dn_eff/dλ, so the linear Taylor
/// expansion of n_eff around λ_0 is
///     n_eff(λ) ≈ n_eff_0 + (λ − λ_0) · (n_eff_0 − n_g_0) / λ_0
/// (slope `(n_eff_0 − n_g_0)/λ_0` chosen to reproduce `n_g_0` at λ_0).
/// Use this for accumulated propagation phase `φ = 2π·n_eff(λ)·L/λ`.
#[inline]
fn n_eff_at_lambda(n_eff_0: f64, n_g_0: f64, wl_ref_m: f64, lambda: f64) -> f64 {
    if wl_ref_m.abs() < 1e-30 { return n_eff_0; }
    n_eff_0 + (lambda - wl_ref_m) * (n_eff_0 - n_g_0) / wl_ref_m
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
                L_um=1000 V_pi_L=2e-3 pin_at_ref=1 alpha_dB_cm=0\n\
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
