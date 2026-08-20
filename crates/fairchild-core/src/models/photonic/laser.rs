use crate::device::{CorrelatedNoise, Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

/// How many of a laser's `wpc` bundle wires it actually drives: the two forward
/// field components and the λ tag, never the backward pair.
///
/// **A laser absorbs backward light; it does not assert that there is none.**
/// Driving `re_bw = 0` looks like the same thing — a perfect absorber — right
/// up until something upstream sends light back.  The neighbouring device also
/// drives that wire (a waveguide's input-side backward wires are its own
/// output), so two branch rows then pin one node to different values, the block
/// goes rank-deficient, and the solve returns a weighted average of the two
/// answers instead of failing.  Measured on a laser → waveguide → `fc_facet`
/// deck: the returned power came out 4x low, with no diagnostic.
///
/// Leaving the wire alone is both correct and simpler.  When nothing else
/// drives it — a laser straight into a photodetector — nothing stamps into that
/// node's row at all, and `stamp_gmin` pins such a row at `V = 0`, so it still
/// reads exactly zero.
fn emitted_wires(wpc: usize) -> usize {
    if wpc == 5 {
        3
    } else {
        wpc
    }
}

/// Bind branch rows for the wires `emitted_wires` counts: `re`, `im`, and λ
/// (the last wire), leaving the backward pair unbound.
fn bind_emitted(branches: &mut [Option<usize>], wpc: usize, first_idx: usize) {
    branches.fill(None);
    branches[0] = Some(first_idx);
    branches[1] = Some(first_idx + 1);
    branches[wpc - 1] = Some(first_idx + 2);
}

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
/// Drives ONE optical channel: `(re, im, λ)` under unidirectional propagation,
/// and the same three under bidirectional — the backward pair is left to
/// whatever sends light back, see `emitted_wires`.
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
    /// Smallest terminal count that would have worked, recorded when
    /// `setup_instance` refuses. `num_terminals()` otherwise reports the
    /// unconfigured 0, which the caller would quote back as the expectation.
    min_terminals: Option<usize>,
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
            min_terminals: None,
            nodes: Vec::new(),
            branches: Vec::new(),
        }
    }
}

impl Device for NativeCwLaser {
    fn num_terminals(&self) -> usize {
        if let Some(min) = self.min_terminals {
            return min;
        }
        self.nodes.len()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wavelen_m = ctx.lambda_center_m;
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        // debug_assert, so a release build sailed past a mis-wired instance and
        // indexed whatever `nodes` happened to hold. Decline instead.
        if terminals.len() != wpc {
            self.min_terminals = Some(wpc);
            return;
        }
        self.min_terminals = None;
        self.nodes = terminals.to_vec();
        self.branches = vec![None; wpc];
    }

    fn num_extra_nodes(&self) -> usize {
        emitted_wires(self.wpc)
    }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        bind_emitted(&mut self.branches, self.wpc, first_idx);
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
        if let Some(j) = self.branches[0] {
            b[j] += self.src_scale * self.re_amp;
        }
        if let Some(j) = self.branches[1] {
            b[j] += self.src_scale * self.im_amp;
        }
        if let Some(j) = self.branches[self.wpc - 1] {
            b[j] += self.wavelen_m;
        }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        // Branch rows: V(out_wire) - target = 0.  Stamp +1 at (J, out) and
        // +1 at (out, J).  RHS handled in load_residual.  The backward wires
        // get no branch at all — see `bind_emitted`.
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

// ────────────────────────────────────────────────────────────────────────
// Voltage-driven laser (direct intensity modulation)
// ────────────────────────────────────────────────────────────────────────

/// A laser whose optical power follows an electrical input, so one SPICE
/// source produces a time-varying optical waveform with no external modulator.
///
/// ```text
///   P(V) = p_floor_w + max(0, slope_w_v · (V − v_th))      V = V(p) − V(n)
///   A    = √P,   re = A·cos φ₀,   im = A·sin φ₀
/// ```
///
/// That is the L–I curve of a diode laser written against voltage: a hard
/// threshold and a straight line above it.  Drive it from a current source
/// instead by working in `I·r_in` — the input resistance is a real parameter,
/// not a numerical crutch, so `slope_w_v · r_in` is the slope efficiency in W/A.
///
/// | Parameter | Default | Meaning |
/// |---|---|---|
/// | `slope_w_v` / `slope_mw_v` | 1e−3 W/V | dP/dV above threshold. |
/// | `v_th` | 0 | Lasing threshold. |
/// | `p_floor_w` | 1e−12 | Below-threshold floor, −90 dBm. |
/// | `r_in` | 1e6 | Input resistance across (p, n). |
/// | `phi_0_deg`, `wavelength_nm`, `rin_db_hz` | as `fc_cw_laser` | |
///
/// **`p_floor_w` is load-bearing, not cosmetic.**  The wires carry `A = √P`,
/// and `dA/dV = slope/(2√P)` diverges as `P → 0`; a laser switched hard off
/// would hand Newton an unbounded Jacobian entry every time the drive crossed
/// threshold.  The floor caps it at `slope/(2·√p_floor)` and doubles as the
/// spontaneous-emission background a real laser has anyway.  At the default it
/// is 90 dB below a 1 mW output, far under any extinction ratio worth quoting.
///
/// The drive derivative is stamped exactly (no frozen coefficient), so adjoint
/// sensitivities reach through the laser — see `Device::frozen_jacobian_columns`
/// for why that distinction matters.
///
/// ponytail: no chirp.  Direct modulation shifts the emission wavelength with
/// carrier density, and λ here is a static tag on a wire.  Modelling it means
/// making the λ wire drive-dependent, which every downstream device would then
/// see move — worth doing only when something in the tree cares about
/// dispersion penalty.
pub struct NativeDrivenLaser {
    slope_w_v: f64,
    v_th: f64,
    p_floor_w: f64,
    /// Width of the threshold knee in volts — see `power_and_slope`. Set to 0
    /// for the old hard corner, which is a worse Jacobian and no more physical.
    v_knee_v: f64,
    r_in: f64,
    phi_0: f64,
    wavelen_m: f64,
    rin_db_hz: Option<f64>,
    src_scale: f64,
    wpc: usize,
    /// Smallest terminal count that would have worked, recorded when
    /// `setup_instance` refuses. `num_terminals()` otherwise reports the
    /// unconfigured 0, which the caller would quote back as the expectation.
    min_terminals: Option<usize>,
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
    /// Linearisation about the previous iterate: amplitude, its derivative
    /// w.r.t. the drive, and the drive itself.
    a_op: f64,
    da_dv: f64,
    v_op: f64,
}

impl Default for NativeDrivenLaser {
    fn default() -> Self {
        Self::new()
    }
}

impl NativeDrivenLaser {
    pub fn new() -> Self {
        Self {
            slope_w_v: 1e-3,
            v_th: 0.0,
            p_floor_w: 1e-12,
            v_knee_v: 1e-3,
            r_in: 1e6,
            phi_0: 0.0,
            wavelen_m: 1550e-9,
            rin_db_hz: None,
            src_scale: 1.0,
            wpc: 3,
            min_terminals: None,
            nodes: Vec::new(),
            branches: Vec::new(),
            a_op: 0.0,
            da_dv: 0.0,
            v_op: 0.0,
        }
    }

    /// Electrical terminals sit after the bundle wires: `[wires…, p, n]`.
    fn drive_nodes(&self) -> (NodeId, NodeId) {
        (self.nodes[self.wpc], self.nodes[self.wpc + 1])
    }

    /// Optical power and `dP/dV` at drive `v`, with a **smooth** threshold:
    ///
    /// ```text
    /// P(V) = p_floor + slope · s · softplus((V − V_th)/s),   s = v_knee
    /// ```
    ///
    /// The knee is physical — a diode laser's L-I curve bends over a few mV as
    /// spontaneous emission gives way to stimulated, it does not corner — but
    /// it is here because a corner does not converge.
    ///
    /// A hard `max(0, ·)` makes `dP/dV` jump from `0` to `slope` at threshold,
    /// so `dA/dV = slope/(2√P)` jumps from `0` to `slope/(2√p_floor)`: at the
    /// default 1 pW floor that is ~750 for a 1.5 mW/V laser, three orders above
    /// anything else in the Jacobian, and it flips on and off as Newton steps
    /// across the threshold. The iterate ping-pongs and never converges. It
    /// only ever appeared to work because a load capacitance adds `C/h` to the
    /// detector diagonal and damps it — remove the capacitor and a link that
    /// modulates through threshold fails at the falling edge with no diagnostic.
    ///
    /// With the softplus, `dP/dV = slope·σ((V−V_th)/s)` is continuous and `P`
    /// at the knee is at least `p_floor + slope·s·ln2`, so `dA/dV` is bounded
    /// by roughly `½·√(slope/(4·s·ln2))` — 0.37 for the same laser, 2000×
    /// smaller. Away from the knee by more than ~20·s the two forms agree to
    /// machine precision, so nothing above threshold moves.
    fn power_and_slope(&self, v: f64) -> (f64, f64) {
        let s = self.v_knee_v;
        let x = v - self.v_th;
        if s <= 0.0 {
            // Explicitly opted out of the smoothing; the hard corner is back.
            let above = self.slope_w_v * x;
            let p = self.p_floor_w + above.max(0.0);
            return (p, if x > 0.0 { self.slope_w_v } else { 0.0 });
        }
        let u = x / s;
        // ln(1+eᵘ) overflows for large u and underflows to 0 for very negative
        // u; both limits are exact, so take them rather than computing them.
        let softplus = if u > 40.0 {
            u
        } else if u < -40.0 {
            0.0
        } else {
            u.exp().ln_1p()
        };
        let sigma = 1.0 / (1.0 + (-u).exp());
        (
            self.p_floor_w + self.slope_w_v * s * softplus,
            self.slope_w_v * sigma,
        )
    }
}

impl Device for NativeDrivenLaser {
    fn num_terminals(&self) -> usize {
        if let Some(min) = self.min_terminals {
            return min;
        }
        self.nodes.len()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wavelen_m = ctx.lambda_center_m;
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        if terminals.len() != wpc + 2 {
            self.min_terminals = Some(wpc + 2);
            return;
        }
        self.min_terminals = None;
        self.nodes = terminals.to_vec();
        self.branches = vec![None; wpc];
    }

    fn num_extra_nodes(&self) -> usize {
        emitted_wires(self.wpc)
    }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        bind_emitted(&mut self.branches, self.wpc, first_idx);
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "slope_w_v" | "slope" => self.slope_w_v = value,
            "slope_mw_v" => self.slope_w_v = value * 1e-3,
            "v_th" | "vth" => self.v_th = value,
            "p_floor_w" => self.p_floor_w = value.max(0.0),
            "v_knee_v" | "v_knee" => self.v_knee_v = value.max(0.0),
            "r_in" => self.r_in = value,
            "phi_0_deg" | "phase_deg" => self.phi_0 = value.to_radians(),
            "wavelength_nm" => self.wavelen_m = value * 1e-9,
            "wavelength_m" => self.wavelen_m = value,
            "rin_db_hz" | "rin" => self.rin_db_hz = Some(value),
            _ => return false,
        }
        true
    }

    fn set_source_scale(&mut self, scale: f64) {
        self.src_scale = scale;
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        let (p_node, n_node) = self.drive_nodes();
        let v = p_node.map_or(0.0, |i| x[i]) - n_node.map_or(0.0, |i| x[i]);
        let (power, dp_dv) = self.power_and_slope(v);
        self.a_op = power.sqrt();
        // dA/dV = (dP/dV) / (2√P).  Both factors are continuous, so this is too
        // — see `power_and_slope` for why that is load-bearing.
        self.da_dv = if self.a_op > 0.0 {
            dp_dv / (2.0 * self.a_op)
        } else {
            0.0
        };
        self.v_op = v;
    }

    fn load_residual(&self, b: &mut [f64]) {
        // Branch row k enforces V(wire_k) − dA/dV·cos·(V_p − V_n) = const, so
        // the RHS carries the linearisation offset A_op − dA/dV·V_op.
        let (re_w, im_w) = (self.phi_0.cos(), self.phi_0.sin());
        let offset = self.src_scale * (self.a_op - self.da_dv * self.v_op);
        if let Some(j) = self.branches[0] {
            b[j] += offset * re_w;
        }
        if let Some(j) = self.branches[1] {
            b[j] += offset * im_w;
        }
        // The backward wires get no branch at all — see `bind_emitted`.
        let lam = self.wpc - 1;
        if let Some(j) = self.branches[lam] {
            b[j] += self.wavelen_m;
        }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        for (i, out_node) in self.nodes.iter().take(self.wpc).enumerate() {
            if let (Some(out), Some(j)) = (*out_node, self.branches[i]) {
                mat.a[j][out] += 1.0;
                mat.a[out][j] += 1.0;
            }
        }
        let (p_node, n_node) = self.drive_nodes();
        let g = self.src_scale * self.da_dv;
        for (i, w) in [(0usize, self.phi_0.cos()), (1, self.phi_0.sin())] {
            let Some(j) = self.branches[i] else { continue };
            if let Some(p) = p_node {
                mat.a[j][p] -= g * w;
            }
            if let Some(n) = n_node {
                mat.a[j][n] += g * w;
            }
        }
        if self.r_in > 0.0 {
            super::stamp_resistor(mat, p_node, n_node, 1.0 / self.r_in);
        }
    }

    fn correlated_noise_sources(&self, _ctx: &SimContext) -> Vec<CorrelatedNoise> {
        rin_source(
            self.rin_db_hz,
            self.a_op * self.phi_0.cos(),
            self.a_op * self.phi_0.sin(),
            &self.branches,
        )
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }
}
