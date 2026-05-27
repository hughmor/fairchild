use crate::device::{Device, EvalFlags, NodeId, ReactiveBranchSpec, ReactiveKind, SimContext};
use crate::mna::MnaMatrix;
use super::{C0, dB_per_cm_to_neper_per_m, n_eff_at_lambda, stamp_potential_eq};

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

