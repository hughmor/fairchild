use super::{dB_per_cm_to_neper_per_m, n_eff_at_lambda, stamp_potential_eq, C0};
use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

/// Input couplings for a `stamp_potential_eq` row: `(node, coefficient)` pairs.
/// Empty in delay mode (the output is driven by a history source on the RHS).
type Couplings<'a> = &'a [(NodeId, f64)];

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
    length_m: f64,
    n_eff: f64, // n_eff at wl_ref_m (default 2.445 for silicon)
    n_g: f64,   // n_g    at wl_ref_m (default 4.2)
    /// Reference wavelength at which `n_eff` and `n_g` are evaluated.  The
    /// dispersion-corrected index `n_eff(λ)` is linearised around this point.
    wl_ref_m: f64,
    alpha_neper_m: f64,
    /// Group delay `τ_g = L·n_g/c` (s).  Computed and exposed for transient
    /// post-processing; this device does not yet implement a true delay line,
    /// so the parameter is informational only at this tier (DC OP and steady-
    /// state spectra are unaffected — τ matters only at modulation
    /// bandwidths comparable to 1/τ).
    tau_g_s: f64,
    // Bootstrap λ for the first NR iterate (x = 0).  Sourced from
    // `SimContext::lambda_center_m` in `setup_model`.
    lambda_bootstrap_m: f64,
    n_channels: usize,
    wpc: usize, // wires_per_channel: 3 (unidir) or 5 (bidir)
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>, // wpc per channel
    c_cached: Vec<f64>,
    s_cached: Vec<f64>,

    // ── Group-delay line state (active only when ctx.waveguide_delay) ────────
    /// True for the current eval iff the delay model is engaged (transient +
    /// `waveguide_delay` option + τ_g > 0).  Governs whether the output port
    /// equations couple to the live input nodes (instantaneous) or to a
    /// history-reconstructed delayed source (delay line).
    delay_active: bool,
    /// Absolute time of the current eval; stamped into the history at commit.
    last_time_s: f64,
    /// Committed timestamps (monotonic increasing).
    hist_t: Vec<f64>,
    /// Committed source amplitudes per timestep.  Each entry is a flat vector
    /// laid out per channel as `[in_re_fw, in_im_fw]` (unidir) or
    /// `[in_re_fw, in_im_fw, out_re_bw, out_im_bw]` (bidir) — the amplitudes
    /// that, delayed by τ_g, drive the opposite port.
    hist_vals: Vec<Vec<f64>>,
    /// Delayed source amplitudes for the current step (same layout as a
    /// `hist_vals` entry), reconstructed at `t − τ_g` in `eval`.
    delayed: Vec<f64>,
}

impl Default for NativeWaveguide {
    fn default() -> Self {
        Self::new()
    }
}

impl NativeWaveguide {
    pub fn new() -> Self {
        // Defaults: classic 500 × 220 nm SOI strip waveguide, straight.
        //  n_eff / n_g at 1550 nm extracted from femwell (see
        //  `scripts/waveguide_simulations/cband_sweep.csv`, strip column).
        //  Phase-shifter device classes (fc_pn_ps, fc_pn_th_ps, fc_pn_ps_cap)
        //  use bent-rib values appropriate to a ring section instead.
        let length_m = 100e-6;
        let n_g = 4.19;
        NativeWaveguide {
            length_m,
            n_eff: 2.445,
            n_g,
            wl_ref_m: 1.55e-6,
            alpha_neper_m: dB_per_cm_to_neper_per_m(2.0),
            tau_g_s: length_m * n_g / C0,
            lambda_bootstrap_m: 1.55e-6,
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
            c_cached: Vec::new(),
            s_cached: Vec::new(),
            delay_active: false,
            last_time_s: 0.0,
            hist_t: Vec::new(),
            hist_vals: Vec::new(),
            delayed: Vec::new(),
        }
    }

    fn refresh_tau(&mut self) {
        self.tau_g_s = self.length_m * self.n_g / C0;
    }

    /// Number of delayed source amplitudes per channel: 2 (forward re/im) for
    /// unidirectional bundles, 4 (+ backward re/im) for bidirectional.
    fn vals_per_channel(&self) -> usize {
        if self.wpc == 5 {
            4
        } else {
            2
        }
    }

    /// Read the source amplitudes that feed delayed ports from the solution
    /// vector: the forward input (and, for bidirectional bundles, the backward
    /// output).  Layout matches `hist_vals` entries.
    fn gather_sources(&self, x: &[f64]) -> Vec<f64> {
        let n = self.n_channels;
        let wpc = self.wpc;
        let out_base = wpc * n;
        let per = self.vals_per_channel();
        let read = |nid: NodeId| nid.map_or(0.0, |i| x[i]);
        let mut v = vec![0.0; n * per];
        for k in 0..n {
            v[per * k] = read(self.nodes[wpc * k]); // in_re_fw
            v[per * k + 1] = read(self.nodes[wpc * k + 1]); // in_im_fw
            if wpc == 5 {
                v[per * k + 2] = read(self.nodes[out_base + wpc * k + 2]); // out_re_bw
                v[per * k + 3] = read(self.nodes[out_base + wpc * k + 3]); // out_im_bw
            }
        }
        v
    }

    /// Reconstruct the source amplitudes at time `tq` by linear interpolation of
    /// the committed history.  Clamps to the endpoints: before the first sample
    /// the input envelope is taken as the initial (DC) value; after the last it
    /// holds the most recent value (delays shorter than the timestep degrade
    /// gracefully to "no delay").
    fn interp_history(&self, tq: f64) -> Vec<f64> {
        let width = self.n_channels * self.vals_per_channel();
        if self.hist_t.is_empty() {
            return vec![0.0; width];
        }
        if tq <= self.hist_t[0] {
            return self.hist_vals[0].clone();
        }
        let last = self.hist_t.len() - 1;
        if tq >= self.hist_t[last] {
            return self.hist_vals[last].clone();
        }
        // Binary search for the bracketing interval [i, i+1].
        let i = match self
            .hist_t
            .binary_search_by(|t| t.partial_cmp(&tq).unwrap_or(std::cmp::Ordering::Less))
        {
            Ok(j) => return self.hist_vals[j].clone(),
            Err(j) => j - 1,
        };
        let (t0, t1) = (self.hist_t[i], self.hist_t[i + 1]);
        let f = if t1 > t0 { (tq - t0) / (t1 - t0) } else { 0.0 };
        let (a, b) = (&self.hist_vals[i], &self.hist_vals[i + 1]);
        (0..width).map(|j| a[j] + f * (b[j] - a[j])).collect()
    }
}

impl Device for NativeWaveguide {
    fn num_terminals(&self) -> usize {
        self.nodes.len()
    }

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
            !terminals.is_empty() && terminals.len().is_multiple_of(stride),
            "fc_waveguide: terminal count must be {stride}·N for N ≥ 1 channels (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = terminals.len() / stride;
        self.n_channels = n;
        self.nodes = terminals.to_vec();
        self.branches = vec![None; wpc * n];
        self.c_cached = vec![1.0; n];
        self.s_cached = vec![0.0; n];
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
            "l_um" => {
                self.length_m = value * 1e-6;
                self.refresh_tau();
                true
            }
            "l_m" | "length" => {
                self.length_m = value;
                self.refresh_tau();
                true
            }
            "n_g" => {
                self.n_g = value;
                self.refresh_tau();
                true
            }
            "n_eff" => {
                self.n_eff = value;
                true
            }
            "wl_ref_m" | "lambda_ref_m" => {
                self.wl_ref_m = value;
                true
            }
            "wl_ref_nm" | "lambda_ref_nm" => {
                self.wl_ref_m = value * 1e-9;
                true
            }
            "alpha_db_cm" => {
                self.alpha_neper_m = dB_per_cm_to_neper_per_m(value);
                true
            }
            _ => false,
        }
    }

    fn eval(&mut self, x: &[f64], flags: EvalFlags, ctx: &SimContext) {
        // Engage the delay line only in transient runs, when the option is on,
        // and when there is a finite group delay to model.  DC/AC and the
        // default (instantaneous) path are unaffected.
        self.last_time_s = ctx.time_s;
        self.delay_active = flags.transient && ctx.waveguide_delay && self.tau_g_s > 0.0;
        if self.delay_active {
            self.delayed = self.interp_history(ctx.time_s - self.tau_g_s);
        }

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
                    if v.abs() > boot * 0.5 {
                        v
                    } else {
                        boot
                    }
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

    fn load_residual(&self, b: &mut [f64]) {
        // Only the delay line contributes to the RHS: each delayed output port
        // equation becomes `V_out = (transmission)·(delayed input)`, a known
        // constant reconstructed from history.  The instantaneous path is fully
        // homogeneous (handled entirely in the Jacobian), so nothing to stamp.
        if !self.delay_active {
            return;
        }
        let n = self.n_channels;
        let per = self.vals_per_channel();
        for k in 0..n {
            let c = self.c_cached[k];
            let s = self.s_cached[k];
            let dly_fw_re = self.delayed[per * k];
            let dly_fw_im = self.delayed[per * k + 1];
            // Forward: out.re = c·in_re(t−τ) + s·in_im(t−τ); out.im = −s·… + c·…
            if let Some(j) = self.branches[self.wpc * k] {
                b[j] += c * dly_fw_re + s * dly_fw_im;
            }
            if let Some(j) = self.branches[self.wpc * k + 1] {
                b[j] += -s * dly_fw_re + c * dly_fw_im;
            }
            if self.wpc == 5 {
                let dly_bw_re = self.delayed[per * k + 2];
                let dly_bw_im = self.delayed[per * k + 3];
                if let Some(j) = self.branches[self.wpc * k + 2] {
                    b[j] += c * dly_bw_re + s * dly_bw_im;
                }
                if let Some(j) = self.branches[self.wpc * k + 3] {
                    b[j] += -s * dly_bw_re + c * dly_bw_im;
                }
            }
        }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let in_block_base = 0;
        let out_block_base = wpc * n;
        // In delay mode the output ports are driven by a history-reconstructed
        // source (stamped into the residual), so their branch equations carry
        // no coupling to the live input nodes — only the identity `V_out = b`.
        let delay = self.delay_active;
        for k in 0..n {
            let c = self.c_cached[k];
            let s = self.s_cached[k];
            // Per-channel input / output wires.
            let in_re_fw = self.nodes[in_block_base + wpc * k];
            let in_im_fw = self.nodes[in_block_base + wpc * k + 1];
            let in_l = self.nodes[in_block_base + wpc * k + (wpc - 1)];
            let out_re_fw = self.nodes[out_block_base + wpc * k];
            let out_im_fw = self.nodes[out_block_base + wpc * k + 1];
            let out_l = self.nodes[out_block_base + wpc * k + (wpc - 1)];
            // Forward path: out.re_fw = c·in.re_fw + s·in.im_fw etc.
            let (re_ins, im_ins): (Couplings, Couplings) = if delay {
                (&[], &[])
            } else {
                (
                    &[(in_re_fw, -c), (in_im_fw, -s)],
                    &[(in_re_fw, s), (in_im_fw, -c)],
                )
            };
            stamp_potential_eq(mat, &self.branches, wpc * k, out_re_fw, re_ins);
            stamp_potential_eq(mat, &self.branches, wpc * k + 1, out_im_fw, im_ins);
            // λ passes through unchanged (a wavelength label is not delayed).
            stamp_potential_eq(
                mat,
                &self.branches,
                wpc * k + (wpc - 1),
                out_l,
                &[(in_l, -1.0)],
            );
            if wpc == 5 {
                // Backward path mirrors the forward physics with reversed
                // direction: light at out_bw propagates back to in_bw with
                // the same c, s.  in.re_bw = c·out.re_bw + s·out.im_bw,
                // in.im_bw = -s·out.re_bw + c·out.im_bw.
                let in_re_bw = self.nodes[in_block_base + wpc * k + 2];
                let in_im_bw = self.nodes[in_block_base + wpc * k + 3];
                let out_re_bw = self.nodes[out_block_base + wpc * k + 2];
                let out_im_bw = self.nodes[out_block_base + wpc * k + 3];
                let (bw_re_ins, bw_im_ins): (Couplings, Couplings) = if delay {
                    (&[], &[])
                } else {
                    (
                        &[(out_re_bw, -c), (out_im_bw, -s)],
                        &[(out_re_bw, s), (out_im_bw, -c)],
                    )
                };
                stamp_potential_eq(mat, &self.branches, wpc * k + 2, in_re_bw, bw_re_ins);
                stamp_potential_eq(mat, &self.branches, wpc * k + 3, in_im_bw, bw_im_ins);
            }
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }

    #[allow(clippy::needless_range_loop)]
    fn commit_timestep(&mut self, x: &[f64]) {
        // Record the port amplitudes for this accepted step so future steps can
        // read them back delayed by τ_g.  Only maintained when the delay line
        // is engaged; otherwise this is a no-op (no memory, no cost).
        if !self.delay_active {
            return;
        }
        self.hist_t.push(self.last_time_s);
        self.hist_vals.push(self.gather_sources(x));
        // Trim history older than one full delay window (keep one sample before
        // the window so interpolation still brackets `t − τ_g`).
        let cutoff = self.last_time_s - self.tau_g_s;
        let mut drop = 0;
        while drop + 1 < self.hist_t.len() && self.hist_t[drop + 1] <= cutoff {
            drop += 1;
        }
        if drop > 0 {
            self.hist_t.drain(0..drop);
            self.hist_vals.drain(0..drop);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::device::EvalFlags;

    /// 1-channel unidirectional waveguide with branch nodes bound at index 6.
    /// Terminal layout: in_re=0, in_im=1, in_λ=2, out_re=3, out_im=4, out_λ=5.
    fn lossless_wg(tau_s: f64) -> NativeWaveguide {
        let ctx = SimContext::default();
        let mut wg = NativeWaveguide::new();
        wg.alpha_neper_m = 0.0; // lossless ⇒ |transmission| = 1
        wg.setup_instance(
            &[Some(0), Some(1), Some(2), Some(3), Some(4), Some(5)],
            &ctx,
        );
        wg.bind_extra_nodes(6);
        wg.tau_g_s = tau_s; // override the geometry-derived value for a clean test
        wg
    }

    fn ctx_delay(t: f64) -> SimContext {
        SimContext {
            waveguide_delay: true,
            time_s: t,
            ..Default::default()
        }
    }

    #[test]
    fn delay_disabled_by_default_and_no_history() {
        let mut wg = lossless_wg(2.0);
        let x = vec![0.7, 0.0, 0.0, 0.0, 0.0, 0.0];
        // Default ctx has waveguide_delay = false → instantaneous path.
        wg.eval(&x, EvalFlags::tran(), &SimContext::default());
        assert!(!wg.delay_active, "delay must be off without the option");
        wg.commit_timestep(&x);
        assert!(
            wg.hist_t.is_empty(),
            "no history should accumulate when off"
        );
        // Residual is homogeneous in the instantaneous path.
        let mut b = vec![0.0; 9];
        wg.load_residual(&mut b);
        assert!(b.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn delay_line_reproduces_delayed_input() {
        // Build a ramp history in_re(t) = t for t = 0,1,2,3, with τ = 2.
        let mut wg = lossless_wg(2.0);
        for step in 0..=3 {
            let t = step as f64;
            let x = vec![t, 0.0, 0.0, 0.0, 0.0, 0.0]; // in_re = t
            wg.eval(&x, EvalFlags::tran(), &ctx_delay(t));
            assert!(wg.delay_active);
            wg.commit_timestep(&x);
        }
        // Query at t = 3 ⇒ delayed source is the input at t − τ = 1 ⇒ in_re = 1.
        let xq = vec![3.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        wg.eval(&xq, EvalFlags::tran(), &ctx_delay(3.0));
        assert!(
            (wg.delayed[0] - 1.0).abs() < 1e-12,
            "delayed in_re should be 1.0 (input at t−τ), got {}",
            wg.delayed[0]
        );
        // Interpolation between samples: query at t − τ = 1.5 ⇒ 1.5.
        wg.eval(&xq, EvalFlags::tran(), &ctx_delay(3.5));
        assert!(
            (wg.delayed[0] - 1.5).abs() < 1e-12,
            "linear interpolation should give 1.5, got {}",
            wg.delayed[0]
        );
    }

    #[test]
    fn delay_preserves_energy_and_stamps_residual() {
        // A lossless waveguide must conserve |amplitude|² through the delay.
        let mut wg = lossless_wg(1.0);
        // History: a complex input (re, im) = (0.6, 0.8) at t = 0, magnitude 1.
        let x0 = vec![0.6, 0.8, 0.0, 0.0, 0.0, 0.0];
        wg.eval(&x0, EvalFlags::tran(), &ctx_delay(0.0));
        wg.commit_timestep(&x0);
        // Query at t = 1 ⇒ delayed source = input at t − τ = 0 = (0.6, 0.8).
        let xq = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        wg.eval(&xq, EvalFlags::tran(), &ctx_delay(1.0));
        let mut b = vec![0.0; 9];
        wg.load_residual(&mut b);
        // Branch rows 6,7 carry the delayed output re/im. Their magnitude is the
        // transmitted (lossless ⇒ unit) magnitude of the delayed input.
        let out_mag2 = b[6] * b[6] + b[7] * b[7];
        assert!(
            (out_mag2 - 1.0).abs() < 1e-9,
            "lossless delay must conserve |A|²=1, got {out_mag2}"
        );
    }
}
