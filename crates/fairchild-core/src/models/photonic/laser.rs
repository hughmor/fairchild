use crate::device::{CorrelatedNoise, Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

/// Relative-intensity-noise generator for a laser emitting `(re_amp, im_amp)`
/// on branch rows `branches`, shared by every laser model in this file.
///
/// `RIN` is defined on POWER: `S_P(f) = RIN·P²` with `RIN = 10^(rin_db_hz/10)`
/// [1/Hz], one-sided, and flat.  The wires carry the field `A = √P`, so
/// `δA = δP/(2√P)` and the injected PSD is `S_A = RIN·P/4` [W/Hz].
///
/// One fluctuation, two wires: `δA` splits onto `re` and `im` by the emission
/// phase and the two arrive perfectly correlated, which is why this returns a
/// `CorrelatedNoise` rather than two independent sources.
///
/// The injection lands on the laser's branch ROWS, whose RHS entries are the
/// enforced field values — perturbing `b[j]` by δ moves the emitted amplitude
/// by δ.  That is the branch-row spelling of a series voltage source, and it is
/// why the PSD is in W/Hz rather than A²/Hz.
///
/// ponytail: flat.  A real diode laser has a relaxation-oscillation peak at a
/// few GHz and a 1/f tail; both want `rin_f_res` / `rin_damping` parameters and
/// a frequency argument on this hook.  Flat is the right first model and is
/// what a datasheet's single "RIN < −155 dB/Hz" number means anyway.
pub(super) fn rin_source(
    rin_db_hz: Option<f64>,
    re_amp: f64,
    im_amp: f64,
    branches: &[Option<usize>],
) -> Vec<CorrelatedNoise> {
    let Some(rin_db) = rin_db_hz else {
        return Vec::new();
    };
    let p = re_amp * re_amp + im_amp * im_amp;
    if p <= 0.0 || branches.len() < 2 {
        return Vec::new();
    }
    let a = p.sqrt();
    vec![CorrelatedNoise {
        psd: 10f64.powf(rin_db / 10.0) * p / 4.0,
        taps: vec![
            (branches[0], None, re_amp / a),
            (branches[1], None, im_amp / a),
        ],
    }]
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
    re_amp: f64,
    im_amp: f64,
    wavelen_m: f64,
    /// Source-stepping homotopy factor on the FIELD amplitude, 0 → 1. The
    /// wavelength tag is never scaled: it is a label, and a ring detuned by a
    /// ramped lambda would be a worse homotopy path than one that is simply dark.
    src_scale: f64,
    /// Relative intensity noise in dB/Hz (e.g. −155).  `None` = noiseless,
    /// which is the default because 0 dB/Hz is not a sane fallback.
    rin_db_hz: Option<f64>,
    wpc: usize, // 3 (unidir) or 5 (bidir)
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
}

impl Default for NativeCwLaser {
    fn default() -> Self {
        Self::new()
    }
}

impl NativeCwLaser {
    pub fn new() -> Self {
        // Defaults: 1 mW, 0° phase, 1550 nm.
        let p = 1e-3_f64;
        Self {
            re_amp: p.sqrt(),
            im_amp: 0.0,
            wavelen_m: 1550e-9,
            src_scale: 1.0,
            rin_db_hz: None,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
        }
    }
}

impl Device for NativeCwLaser {
    fn num_terminals(&self) -> usize {
        self.nodes.len()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wavelen_m = ctx.lambda_center_m;
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        debug_assert_eq!(
            terminals.len(),
            wpc,
            "fc_cw_laser: expected {wpc} terminals (one channel × wpc); got {}",
            terminals.len()
        );
        self.nodes = terminals.to_vec();
        self.branches = vec![None; wpc];
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
            "wavelength_nm" => {
                self.wavelen_m = value * 1e-9;
                true
            }
            "wavelength_m" => {
                self.wavelen_m = value;
                true
            }
            "rin_db_hz" | "rin" => {
                self.rin_db_hz = Some(value);
                true
            }
            "re_amp" => {
                self.re_amp = value;
                true
            }
            "im_amp" => {
                self.im_amp = value;
                true
            }
            _ => false,
        }
    }

    fn set_source_scale(&mut self, scale: f64) {
        self.src_scale = scale;
    }

    fn eval(&mut self, _x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {}

    fn load_residual(&self, b: &mut [f64]) {
        // Inhomogeneous branch equations: V(out_re_fw) = re_amp, ...
        // Wire order: [re_fw, im_fw, re_bw, im_bw, λ] (5-wire bidir) or
        //             [re,    im,    λ]               (3-wire unidir).
        if self.wpc == 5 {
            if let Some(j) = self.branches[0] {
                b[j] += self.src_scale * self.re_amp;
            }
            if let Some(j) = self.branches[1] {
                b[j] += self.src_scale * self.im_amp;
            }
            // bw wires forced to 0 — no contribution from RHS (branch row
            // already enforces V = 0 because rhs is 0).
            if let Some(j) = self.branches[4] {
                b[j] += self.wavelen_m;
            }
        } else {
            if let Some(j) = self.branches[0] {
                b[j] += self.src_scale * self.re_amp;
            }
            if let Some(j) = self.branches[1] {
                b[j] += self.src_scale * self.im_amp;
            }
            if let Some(j) = self.branches[2] {
                b[j] += self.wavelen_m;
            }
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

    fn correlated_noise_sources(&self, _ctx: &SimContext) -> Vec<CorrelatedNoise> {
        rin_source(self.rin_db_hz, self.re_amp, self.im_amp, &self.branches)
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }
}
