use super::dB_per_cm_to_neper_per_m;
use super::segment::{OpticalSegment, PerChannel};
use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

// ────────────────────────────────────────────────────────────────────────
// Native straight waveguide
// ────────────────────────────────────────────────────────────────────────

/// Straight optical waveguide — propagation loss + accumulated phase.
///
/// Physics: `A_out = A_in · exp(-α·L/2) · exp(-j·β·L)` with `β = 2π·n_eff(λ)/λ`.
///
/// This is the simplest [`OpticalSegment`]-based device: a passive segment with
/// **no** electrical drive (zero `Δn_eff` / `Δα` perturbation). The optical
/// propagation, λ bootstrap, branch stamping, and (opt-in) group delay all live
/// in the shared segment; this struct is a thin `Device` adapter over it. The
/// active phase-shifter / modulator classes reuse the same segment and add a
/// perturbation source on top.
///
/// Variable-arity bundle-aware device. Terminal layout for N channels:
///   [in.0.re, in.0.im, in.0.λ, …, out.{N-1}.λ]  (6·N terminals, wpc=3).
pub struct NativeWaveguide {
    seg: OpticalSegment,
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
        NativeWaveguide {
            seg: OpticalSegment::new(100e-6, 2.445, 4.19, dB_per_cm_to_neper_per_m(2.0)),
        }
    }
}

impl Device for NativeWaveguide {
    fn num_terminals(&self) -> usize {
        self.seg.num_optical_wires()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.seg.setup_model(ctx);
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        // The waveguide is pure-optical: every terminal is a bundle wire.
        self.seg.setup_instance(terminals, ctx);
    }

    fn num_extra_nodes(&self) -> usize {
        self.seg.num_aux_branches()
    }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        self.seg.bind_branches(first_idx);
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        self.seg.set_param(&name.to_lowercase(), value)
    }

    fn lambda_routing(&self) -> Vec<(usize, usize)> {
        self.seg.lambda_routing()
    }

    fn eval(&mut self, x: &[f64], flags: EvalFlags, ctx: &SimContext) {
        // Engage the delay line only in transient runs, when the option is on,
        // and when there is a finite group delay. DC/AC and the default
        // (instantaneous) path are unaffected. Passive ⇒ no perturbation.
        let delay_active = flags.transient && ctx.waveguide_delay && self.seg.tau_g_s() > 0.0;
        self.seg.refresh(
            x,
            PerChannel::zero(),
            0.0,
            PerChannel::zero(),
            delay_active,
            ctx,
        );
    }

    fn load_residual(&self, b: &mut [f64]) {
        self.seg.stamp_residual(b);
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        self.seg.stamp(mat);
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.seg.stamp_residual(b);
    }

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.seg.stamp(mat);
    }

    fn commit_timestep(&mut self, x: &[f64]) {
        self.seg.commit(x);
    }
}
