//! `OpticalSegment` — the shared optical core of every active/passive
//! waveguide-like photonic device.
//!
//! A segment owns the per-channel re/im/λ bundle propagation
//! (`A_out = A_in · exp(−(α+Δα)·L/2) · exp(−j·(β·L + Δφ))`), the auxiliary
//! branch-row stamping, the group-delay line, and the λ bootstrap. It knows
//! *nothing* about why the mode is perturbed — callers hand it a per-segment
//! effective-index change `Δn_eff` and excess loss `Δα`, computed by whatever
//! physics drives the device (free-carrier dispersion, thermo-optic, Pockels,
//! a photoconductive back-action, or an external model). That keeps the optical
//! math in one place and lets new device physics be a new perturbation source
//! rather than a new copy of the stamp loop.
//!
//! See `_notes/optical_abstraction_design.md` for the layering rationale and
//! the future device classes this is designed to admit.

use super::{dB_per_cm_to_neper_per_m, n_eff_at_lambda, stamp_potential_eq, C0};
use crate::delay::DelayLine;
use crate::device::{NodeId, SimContext};
use crate::mna::MnaMatrix;

/// Input couplings for a `stamp_potential_eq` row: `(node, coefficient)` pairs.
/// Empty in delay mode (the output is driven by a history source on the RHS).
type Couplings<'a> = &'a [(NodeId, f64)];

/// The optical core: bundle propagation + perturbation + delay + stamping,
/// shared by the waveguide and every active phase-shifter / modulator class.
pub struct OpticalSegment {
    pub length_m: f64,
    pub n_eff: f64,
    pub n_g: f64,
    /// Reference wavelength at which `n_eff`/`n_g` are quoted; dispersion is
    /// linearised around it.
    pub wl_ref_m: f64,
    pub alpha_neper_m: f64,
    /// When true, subtract the absolute propagation phase at `wl_ref_m` so the
    /// segment is "transparent" at λ_ref (testbench-ring convenience). The
    /// passive waveguide leaves this false (physical absolute phase).
    pub pin_at_ref: bool,
    /// Group delay `τ_g = L·n_g/c` (s).
    tau_g_s: f64,
    /// Bootstrap λ for the first NR iterate (x = 0): `wl_ref_m`.
    lambda_bootstrap_m: f64,
    n_channels: usize,
    wpc: usize, // wires per channel: 3 (unidir) or 5 (bidir)
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
    c_cached: Vec<f64>,
    s_cached: Vec<f64>,
    // ── Group-delay line (engaged only when the owning device asks) ──────────
    delay: DelayLine,
    delayed: Vec<f64>,
}

impl OpticalSegment {
    /// Construct a segment with the given geometry. `n_eff`/`n_g`/`wl_ref_m`/
    /// `alpha_neper_m` are typically overridden by the owning device's defaults
    /// and `.model`/instance params.
    pub fn new(length_m: f64, n_eff: f64, n_g: f64, alpha_neper_m: f64) -> Self {
        OpticalSegment {
            length_m,
            n_eff,
            n_g,
            wl_ref_m: 1.55e-6,
            alpha_neper_m,
            pin_at_ref: false,
            tau_g_s: length_m * n_g / C0,
            lambda_bootstrap_m: 1.55e-6,
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
            c_cached: Vec::new(),
            s_cached: Vec::new(),
            delay: DelayLine::new(),
            delayed: Vec::new(),
        }
    }

    pub fn tau_g_s(&self) -> f64 {
        self.tau_g_s
    }

    pub fn n_channels(&self) -> usize {
        self.n_channels
    }

    pub fn wpc(&self) -> usize {
        self.wpc
    }

    pub fn refresh_tau(&mut self) {
        self.tau_g_s = self.length_m * self.n_g / C0;
    }

    /// Pick up the band-centre defaults from the context (call from
    /// `Device::setup_model`). Defaults the dispersion reference to the band
    /// centre too, unless the user already overrode `wl_ref_m`.
    pub fn setup_model(&mut self, ctx: &SimContext) {
        self.lambda_bootstrap_m = ctx.lambda_center_m;
        if (self.wl_ref_m - 1.55e-6).abs() < 1e-12 {
            self.wl_ref_m = ctx.lambda_center_m;
        }
        self.wpc = ctx.wires_per_channel();
        self.refresh_tau();
    }

    /// Bind the segment to its `2·wpc·N` optical bundle wires (in block then out
    /// block). Allocates the per-channel branch slots and the c/s caches.
    pub fn setup_instance(&mut self, optical_terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 2 * wpc; // in + out
        assert!(
            !optical_terminals.is_empty() && optical_terminals.len().is_multiple_of(stride),
            "OpticalSegment: optical wire count must be {stride}·N for N ≥ 1 channels \
             (wpc={wpc}); got {}",
            optical_terminals.len()
        );
        let n = optical_terminals.len() / stride;
        self.n_channels = n;
        self.nodes = optical_terminals.to_vec();
        self.branches = vec![None; wpc * n];
        self.c_cached = vec![1.0; n];
        self.s_cached = vec![0.0; n];
    }

    /// Number of optical bundle wires this segment occupies (`2·wpc·N`).
    pub fn num_optical_wires(&self) -> usize {
        2 * self.wpc * self.n_channels
    }

    /// Number of auxiliary branch rows this segment needs.
    pub fn num_aux_branches(&self) -> usize {
        self.branches.len()
    }

    /// Bind the auxiliary branch rows to a contiguous block starting at
    /// `first_idx`.
    pub fn bind_branches(&mut self, first_idx: usize) {
        for (i, b) in self.branches.iter_mut().enumerate() {
            *b = Some(first_idx + i);
        }
    }

    /// Try to set a geometry parameter shared by all segment-based devices.
    /// Returns `true` if recognised. Device-specific (electrical) params are
    /// handled by the owning device / active model.
    pub fn set_param(&mut self, name: &str, value: f64) -> bool {
        match name {
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
            "pin_at_ref" => {
                self.pin_at_ref = value != 0.0;
                true
            }
            _ => false,
        }
    }

    /// Read the per-channel λ wire (bootstrapped to `wl_ref_m` when undriven).
    fn lambda_of(&self, x: &[f64], k: usize) -> f64 {
        let lam = self.wpc - 1;
        match self.nodes[self.wpc * k + lam] {
            Some(i) => {
                let v = x[i];
                if v.abs() > 1e-9 {
                    v
                } else {
                    self.wl_ref_m
                }
            }
            None => self.wl_ref_m,
        }
    }

    /// Per-channel input optical intensity |A_in|² (W), the input to any
    /// photoconductive / detection back-action. Forward channel only.
    pub fn channel_intensities(&self, x: &[f64]) -> Vec<f64> {
        let read = |nid: NodeId| nid.map_or(0.0, |i| x[i]);
        (0..self.n_channels)
            .map(|k| {
                let re = read(self.nodes[self.wpc * k]);
                let im = read(self.nodes[self.wpc * k + 1]);
                re * re + im * im
            })
            .collect()
    }

    /// Recompute the cached transmission coefficients for every channel from
    /// the segment's own propagation plus the supplied per-segment perturbation
    /// `Δn_eff` (added to the effective index) and `Δα` (added to the loss, in
    /// Neper/m). A passive segment passes `dn_eff = dalpha = 0`.
    ///
    /// `delay_active` engages the group-delay line (the owning device decides
    /// when — e.g. `flags.transient && ctx.waveguide_delay && τ>0`).
    pub fn refresh(
        &mut self,
        x: &[f64],
        dn_eff: f64,
        dalpha_neper_m: f64,
        delay_active: bool,
        ctx: &SimContext,
    ) {
        self.delay.set_state(delay_active, ctx.time_s);
        if delay_active {
            let width = self.n_channels * self.vals_per_channel();
            self.delayed = self.delay.sample(self.tau_g_s, width);
        }

        let two_pi = 2.0 * std::f64::consts::PI;
        let t_amp = (-(self.alpha_neper_m + dalpha_neper_m) * self.length_m / 2.0).exp();
        let phi_ref = if self.pin_at_ref {
            two_pi * self.n_eff * self.length_m / self.wl_ref_m
        } else {
            0.0
        };
        for k in 0..self.n_channels {
            let lambda = self.lambda_of(x, k);
            let n_eff_lam = n_eff_at_lambda(self.n_eff, self.n_g, self.wl_ref_m, lambda);
            // Δn_eff folds into the index, exactly as φ_eo = 2π·Δn_eff·L/λ.
            let phi = two_pi * (n_eff_lam + dn_eff) * self.length_m / lambda - phi_ref;
            self.c_cached[k] = t_amp * phi.cos();
            self.s_cached[k] = t_amp * phi.sin();
        }
    }

    /// Stamp the per-channel optical branch equations. In delay mode the output
    /// ports are driven by a history-reconstructed source (see
    /// [`stamp_residual`](Self::stamp_residual)), so their couplings to the live
    /// input nodes are dropped. The delay state is the one set by the most
    /// recent [`refresh`](Self::refresh).
    pub fn stamp(&self, mat: &mut MnaMatrix) {
        let delay_active = self.delay.is_active();
        let n = self.n_channels;
        let wpc = self.wpc;
        let out_base = wpc * n;
        let lam = wpc - 1;
        for k in 0..n {
            let c = self.c_cached[k];
            let s = self.s_cached[k];
            let in_re_fw = self.nodes[wpc * k];
            let in_im_fw = self.nodes[wpc * k + 1];
            let in_l = self.nodes[wpc * k + lam];
            let out_re_fw = self.nodes[out_base + wpc * k];
            let out_im_fw = self.nodes[out_base + wpc * k + 1];
            let out_l = self.nodes[out_base + wpc * k + lam];
            let (re_ins, im_ins): (Couplings, Couplings) = if delay_active {
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
            stamp_potential_eq(mat, &self.branches, wpc * k + lam, out_l, &[(in_l, -1.0)]);
            if wpc == 5 {
                let in_re_bw = self.nodes[wpc * k + 2];
                let in_im_bw = self.nodes[wpc * k + 3];
                let out_re_bw = self.nodes[out_base + wpc * k + 2];
                let out_im_bw = self.nodes[out_base + wpc * k + 3];
                let (bw_re_ins, bw_im_ins): (Couplings, Couplings) = if delay_active {
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

    /// Stamp the delay-line history source onto the RHS. No-op when the delay
    /// line is inactive (the instantaneous path is fully homogeneous).
    pub fn stamp_residual(&self, b: &mut [f64]) {
        if !self.delay.is_active() {
            return;
        }
        let n = self.n_channels;
        let per = self.vals_per_channel();
        for k in 0..n {
            let c = self.c_cached[k];
            let s = self.s_cached[k];
            let dly_fw_re = self.delayed[per * k];
            let dly_fw_im = self.delayed[per * k + 1];
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

    /// Record the current port amplitudes so future steps can read them back
    /// delayed by `τ_g`. No-op when the delay line is inactive.
    pub fn commit(&mut self, x: &[f64]) {
        if self.delay.is_active() {
            let snapshot = self.gather_sources(x);
            self.delay.record(snapshot, self.tau_g_s);
        }
    }

    /// Whether the group-delay line is currently engaged.
    pub fn delay_active(&self) -> bool {
        self.delay.is_active()
    }

    /// Whether the delay history is empty (no accepted timesteps recorded).
    pub fn delay_is_empty(&self) -> bool {
        self.delay.is_empty()
    }

    /// The delayed source amplitudes reconstructed by the most recent
    /// [`refresh`](Self::refresh) (empty unless the delay line is active).
    pub fn delayed(&self) -> &[f64] {
        &self.delayed
    }

    /// Override the geometry-derived group delay (used by tests that want a
    /// clean integer τ independent of L/n_g).
    pub fn set_tau_g_s(&mut self, t: f64) {
        self.tau_g_s = t;
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

    /// Read the source amplitudes feeding delayed ports (forward input, plus
    /// backward output for bidirectional bundles). Layout matches the `delayed`
    /// buffer.
    fn gather_sources(&self, x: &[f64]) -> Vec<f64> {
        let n = self.n_channels;
        let wpc = self.wpc;
        let out_base = wpc * n;
        let per = self.vals_per_channel();
        let read = |nid: NodeId| nid.map_or(0.0, |i| x[i]);
        let mut v = vec![0.0; n * per];
        for k in 0..n {
            v[per * k] = read(self.nodes[wpc * k]);
            v[per * k + 1] = read(self.nodes[wpc * k + 1]);
            if wpc == 5 {
                v[per * k + 2] = read(self.nodes[out_base + wpc * k + 2]);
                v[per * k + 3] = read(self.nodes[out_base + wpc * k + 3]);
            }
        }
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 1-channel unidirectional lossless segment, branches bound at index 6.
    /// Optical layout: in_re=0, in_im=1, in_λ=2, out_re=3, out_im=4, out_λ=5.
    fn lossless_segment(tau_s: f64) -> OpticalSegment {
        let ctx = SimContext::default();
        let mut seg = OpticalSegment::new(100e-6, 2.445, 4.19, 0.0); // α=0 ⇒ |t|=1
        seg.setup_instance(&[Some(0), Some(1), Some(2), Some(3), Some(4), Some(5)], &ctx);
        seg.bind_branches(6);
        seg.set_tau_g_s(tau_s);
        seg
    }

    fn ctx_delay(t: f64) -> SimContext {
        SimContext {
            waveguide_delay: true,
            time_s: t,
            ..Default::default()
        }
    }

    #[test]
    fn delay_off_by_default_no_history() {
        let mut seg = lossless_segment(2.0);
        let x = vec![0.7, 0.0, 0.0, 0.0, 0.0, 0.0];
        // delay_active=false ⇒ instantaneous path.
        seg.refresh(&x, 0.0, 0.0, false, &SimContext::default());
        assert!(!seg.delay_active());
        seg.commit(&x);
        assert!(seg.delay_is_empty(), "no history accumulates when off");
        let mut b = vec![0.0; 9];
        seg.stamp_residual(&mut b);
        assert!(b.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn delay_line_reproduces_delayed_input() {
        let mut seg = lossless_segment(2.0);
        for step in 0..=3 {
            let t = step as f64;
            let x = vec![t, 0.0, 0.0, 0.0, 0.0, 0.0]; // in_re = t
            seg.refresh(&x, 0.0, 0.0, true, &ctx_delay(t));
            assert!(seg.delay_active());
            seg.commit(&x);
        }
        // Query at t = 3 ⇒ delayed source is the input at t − τ = 1.
        let xq = vec![3.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        seg.refresh(&xq, 0.0, 0.0, true, &ctx_delay(3.0));
        assert!((seg.delayed()[0] - 1.0).abs() < 1e-12, "got {}", seg.delayed()[0]);
        // Linear interpolation between samples: t − τ = 1.5 ⇒ 1.5.
        seg.refresh(&xq, 0.0, 0.0, true, &ctx_delay(3.5));
        assert!((seg.delayed()[0] - 1.5).abs() < 1e-12, "got {}", seg.delayed()[0]);
    }

    #[test]
    fn delay_preserves_energy_and_stamps_residual() {
        let mut seg = lossless_segment(1.0);
        let x0 = vec![0.6, 0.8, 0.0, 0.0, 0.0, 0.0]; // |A| = 1
        seg.refresh(&x0, 0.0, 0.0, true, &ctx_delay(0.0));
        seg.commit(&x0);
        let xq = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        seg.refresh(&xq, 0.0, 0.0, true, &ctx_delay(1.0));
        let mut b = vec![0.0; 9];
        seg.stamp_residual(&mut b);
        let out_mag2 = b[6] * b[6] + b[7] * b[7];
        assert!(
            (out_mag2 - 1.0).abs() < 1e-9,
            "lossless delay must conserve |A|²=1, got {out_mag2}"
        );
    }
}
