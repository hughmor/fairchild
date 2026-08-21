use super::{stamp_potential_eq, LambdaSelect};
use crate::delay::DelayLine;
use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

// ────────────────────────────────────────────────────────────────────────
// Behavioural 2×2 optical transfer block (`fc_optical_2x2`)
// ────────────────────────────────────────────────────────────────────────

/// Per-channel 2×2 complex optical transfer matrix — a behavioural stand-in
/// for any 2-in/2-out photonic block whose *response* you know but whose
/// internals you don't want to simulate.
///
/// ```text
///   [ thru ]   [ s11  s12 ] [ in1 ]
///   [ drop ] = [ s21  s22 ] [ in2 ]
/// ```
///
/// applied independently to every WDM channel, with its own matrix per
/// channel. That is the whole point: a cascade of N ring modulators sharing a
/// through bus and a drop bus collapses to one instance with N weights, which
/// (a) removes the rings' free parameters from a fit and (b) removes their
/// resonance, so the transient timestep is set by the electronics instead of a
/// sub-round-trip optical constraint.
///
/// # Terminals — `4·wpc·N + N + 1`
///
/// Optical bundle wires first (all channels of each port, in port order), then
/// the electrical control wires, then one shared control return:
///
/// ```text
///   in1[0..N] in2[0..N] thru[0..N] drop[0..N]  wctl_0 … wctl_{N-1}  ctl_ret
/// ```
///
/// with `wpc` = 3 (default) or 5 (bidirectional). N is recovered as
/// `(len − 1) / (4·wpc + 1)`, which is unique — the shape is fixed rather than
/// inferred among several, because e.g. "static, no control" (`4·wpc·N`) and
/// "shared control pair" (`4·wpc·N + 2`) collide with this form at small N.
///
/// Declare it with vector ports and the width follows the netlist:
///
/// ```text
///   .optical_port     in1  4
///   .optical_port     in2  4
///   .optical_port     thru 4
///   .optical_port     drop 4
///   .electrical_port  wctl 4
///   Xwb in1 in2 thru drop wctl 0 fc_optical_2x2 w=0 dw_dv=0.4
/// ```
///
/// A width disagreement between the optical and control buses is rejected by
/// the parser (`expand_bundle_ports`), which is the only layer that still knows
/// each port's declared width; the assert here is the backstop for hand-written
/// flat netlists.
///
/// # Parameters
///
/// Two ways to fill the matrix. **Weight mode** (the default) takes one real
/// bipolar weight per channel and builds the lossless coupler-form matrix
///
/// ```text
///   s11 = s22 = cos θ,   s12 = s21 = −j sin θ,   θ = ½·acos(−w)
/// ```
///
/// chosen so that `P_drop − P_thru = w · P_in` — i.e. `w` *is* the weight a
/// balanced photodetector pair reads, with `w = −1` all-through, `+1` all-drop,
/// `0` a 50/50 split. Passivity is automatic.
///
/// - `w`, `w_<k>` — weight, clamped to [−1, 1]. Unindexed broadcasts to every
///   channel; indexed overrides one.
/// - `dw_dv`, `dw_dv_<k>` — weight per volt on that channel's control wire:
///   `w_k(V) = w_<k> + dw_dv_<k> · (V(wctl_k) − V(ctl_ret))`, clamped. This is
///   what makes the weight vary *during* a transient, which `set_param` can't
///   do. Default 0 (the control wires are then read but ignored).
///
/// **Explicit-matrix mode** takes the four entries directly, for a block whose
/// response isn't a lossless split (an asymmetric coupler, a measured 2×2, an
/// MZM bias point). Setting any of these switches the channel out of weight
/// mode:
///
/// - `s11_mag_<k>`, `s11_deg_<k>`, … through `s22_deg_<k>` — magnitude and
///   phase in degrees. Unindexed forms broadcast.
///
/// Shared across channels:
///
/// - `il_db` — extra insertion loss in dB applied to every entry (power dB;
///   default 0).
/// - `tau_s` — latency in seconds (default 0 = instantaneous). Engages a
///   [`DelayLine`] in transient analysis. Note the cost: resolving a latency
///   needs a timestep of order `tau_s`, so leave it at 0 when you want the
///   speed and set it when you're studying the delay.
///
///   What is delayed is the *field*: the output is `S(t) · in(t − τ)`, matching
///   `OpticalSegment` — the matrix is evaluated at the current time while the
///   fields carry the latency. So a step on a control voltage reaches the output
///   immediately; only a step on the input field is delayed. Modelling the
///   light already in flight under the *old* transfer would need the matrix
///   history too, which no device here does.
/// - `allow_gain` — set non-zero to permit a matrix whose largest singular
///   value exceeds 1. Off by default: a hand-typed matrix with gain inside a
///   feedback loop diverges silently, which is a miserable thing to debug.
///
/// # Not yet supported
///
/// Bidirectional mode (`.options enable_bidirectional=1`, `wpc = 5`) — the
/// backward-travelling fields would need their own (generally transposed)
/// matrix, and reflection entries only become meaningful there. `setup_instance`
/// rejects it outright rather than leaving the backward wires undriven.
pub struct NativeOptical2x2 {
    /// Per-channel matrix, `[s11, s12, s21, s22]` as (re, im) pairs.
    s: Vec<[(f64, f64); 4]>,
    /// Per-channel static weight and control sensitivity (weight mode only).
    w0: Vec<f64>,
    dw_dv: Vec<f64>,
    /// Channels switched to an explicit matrix (never rebuilt from `w`).
    explicit: Vec<bool>,
    /// Pending broadcast/indexed params, applied once N is known.
    pend: Vec<(String, f64)>,
    il_amp: f64,
    tau_s: f64,
    allow_gain: bool,
    /// Passivity guard runs once, on the first eval (see `eval`).
    checked: bool,
    n_channels: usize,
    wpc: usize,
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
    delay: DelayLine,
    delayed: Vec<f64>,
    lam_src: LambdaSelect,
}

/// Branch rows per channel: thru_re, thru_im, drop_re, drop_im, thru_λ, drop_λ.
const BRANCHES_PER_CHANNEL: usize = 6;
/// Delay-line snapshot values per channel: in1_re, in1_im, in2_re, in2_im.
const DELAY_VALS_PER_CHANNEL: usize = 4;

/// One stamped output row: (branch slot, output node, four weighted inputs).
type StampRow = (usize, NodeId, [(NodeId, f64); 4]);

impl Default for NativeOptical2x2 {
    fn default() -> Self {
        Self::new()
    }
}

impl NativeOptical2x2 {
    pub fn new() -> Self {
        Self {
            s: Vec::new(),
            w0: Vec::new(),
            dw_dv: Vec::new(),
            explicit: Vec::new(),
            pend: Vec::new(),
            il_amp: 1.0,
            tau_s: 0.0,
            allow_gain: false,
            checked: false,
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
            delay: DelayLine::new(),
            delayed: Vec::new(),
            lam_src: LambdaSelect::default(),
        }
    }

    /// Lossless split matrix for balanced-PD weight `w`: `P_drop − P_thru = w·P_in`.
    fn matrix_for_weight(w: f64) -> [(f64, f64); 4] {
        let theta = 0.5 * w.clamp(-1.0, 1.0).mul_add(-1.0, 0.0).acos();
        let (t, k) = (theta.cos(), theta.sin());
        // Coupler form [[t, −jk], [−jk, t]] — unitary for any θ.
        [(t, 0.0), (0.0, -k), (0.0, -k), (t, 0.0)]
    }

    /// Split a param name into (base, channel index) — `"w_3"` → `("w", Some(3))`.
    /// Only a trailing all-digit segment counts, so `dw_dv` keeps its name.
    fn split_index(name: &str) -> (&str, Option<usize>) {
        match name.rsplit_once('_') {
            Some((base, idx)) if !idx.is_empty() && idx.bytes().all(|b| b.is_ascii_digit()) => {
                (base, idx.parse().ok())
            }
            _ => (name, None),
        }
    }

    /// Largest singular value of a channel's 2×2, for the passivity guard.
    /// σ_max² is the larger eigenvalue of SᴴS, which for a 2×2 is a closed form.
    fn sigma_max(m: &[(f64, f64); 4]) -> f64 {
        let sq = |c: (f64, f64)| c.0 * c.0 + c.1 * c.1;
        // Frobenius² and |det|² fully determine the singular values of a 2×2.
        let fro2 = sq(m[0]) + sq(m[1]) + sq(m[2]) + sq(m[3]);
        // det = s11·s22 − s12·s21 (complex)
        let det_re = m[0].0 * m[3].0 - m[0].1 * m[3].1 - (m[1].0 * m[2].0 - m[1].1 * m[2].1);
        let det_im = m[0].0 * m[3].1 + m[0].1 * m[3].0 - (m[1].0 * m[2].1 + m[1].1 * m[2].0);
        let det2 = det_re * det_re + det_im * det_im;
        let disc = (fro2 * fro2 - 4.0 * det2).max(0.0).sqrt();
        (0.5 * (fro2 + disc)).max(0.0).sqrt()
    }

    /// Apply one (name, value) to the sized per-channel vectors. `None` index
    /// broadcasts. Returns false for an unrecognised name.
    fn apply_param(&mut self, name: &str, value: f64) -> bool {
        let (base, idx) = Self::split_index(name);
        let n = self.n_channels;
        let targets: Vec<usize> = match idx {
            Some(k) if k < n => vec![k],
            Some(_) => return true, // out-of-range index: silently ignore, as elsewhere
            None => (0..n).collect(),
        };
        // Weight mode.
        if base == "w" {
            for &k in &targets {
                self.w0[k] = value.clamp(-1.0, 1.0);
            }
            return true;
        }
        if base == "dw_dv" {
            for &k in &targets {
                self.dw_dv[k] = value;
            }
            return true;
        }
        // Explicit-matrix mode: s{11,12,21,22}_{mag,deg}.
        let entry = |b: &str| match b {
            "s11" => Some(0),
            "s12" => Some(1),
            "s21" => Some(2),
            "s22" => Some(3),
            _ => None,
        };
        for (suffix, is_mag) in [("_mag", true), ("_deg", false)] {
            if let Some(stripped) = base.strip_suffix(suffix) {
                if let Some(e) = entry(stripped) {
                    for &k in &targets {
                        let (re, im) = self.s[k][e];
                        let (mag, ph) = (re.hypot(im), im.atan2(re));
                        let (mag, ph) = if is_mag {
                            (value, ph)
                        } else {
                            (mag, value.to_radians())
                        };
                        self.s[k][e] = (mag * ph.cos(), mag * ph.sin());
                        self.explicit[k] = true;
                    }
                    return true;
                }
            }
        }
        false
    }

    /// Rebuild every weight-mode channel's matrix from its control voltage.
    fn refresh(&mut self, x: &[f64], delay_active: bool, ctx: &SimContext) {
        self.delay.set_state(delay_active, ctx.time_s);
        if delay_active {
            self.delayed = self
                .delay
                .sample(self.tau_s, self.n_channels * DELAY_VALS_PER_CHANNEL);
        }
        let ctl_base = 4 * self.wpc * self.n_channels;
        let ret = self.nodes[ctl_base + self.n_channels];
        let v_ret = ret.map_or(0.0, |i| x[i]);
        for k in 0..self.n_channels {
            if self.explicit[k] {
                continue;
            }
            let v = self.nodes[ctl_base + k].map_or(0.0, |i| x[i]) - v_ret;
            let w = (self.w0[k] + self.dw_dv[k] * v).clamp(-1.0, 1.0);
            self.s[k] = Self::matrix_for_weight(w);
        }
    }

    /// Node ids for channel `k`: (in1_re, in1_im, in1_λ, in2_re, in2_im,
    /// thru_re, thru_im, thru_λ, drop_re, drop_im, drop_λ).
    #[allow(clippy::type_complexity)]
    /// The second input port's λ wire — `channel_nodes` skips it because the
    /// transfer matrix never needs it, but λ selection does.
    fn channel_lambda_in2(&self, k: usize) -> NodeId {
        let (wpc, n) = (self.wpc, self.n_channels);
        self.nodes[wpc * n + wpc * k + wpc - 1]
    }

    fn channel_nodes(
        &self,
        k: usize,
    ) -> (
        NodeId,
        NodeId,
        NodeId,
        NodeId,
        NodeId,
        NodeId,
        NodeId,
        NodeId,
        NodeId,
        NodeId,
        NodeId,
    ) {
        let (wpc, n) = (self.wpc, self.n_channels);
        let lam = wpc - 1;
        let (p1, p2, p3, p4) = (0, wpc * n, 2 * wpc * n, 3 * wpc * n);
        (
            self.nodes[p1 + wpc * k],
            self.nodes[p1 + wpc * k + 1],
            self.nodes[p1 + wpc * k + lam],
            self.nodes[p2 + wpc * k],
            self.nodes[p2 + wpc * k + 1],
            self.nodes[p3 + wpc * k],
            self.nodes[p3 + wpc * k + 1],
            self.nodes[p3 + wpc * k + lam],
            self.nodes[p4 + wpc * k],
            self.nodes[p4 + wpc * k + 1],
            self.nodes[p4 + wpc * k + lam],
        )
    }
}

impl Device for NativeOptical2x2 {
    fn num_terminals(&self) -> usize {
        self.nodes.len()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        assert!(
            wpc == 3,
            "fc_optical_2x2: bidirectional propagation is not supported yet \
             (wpc={wpc}); the backward-travelling fields would need their own \
             transfer matrix. Drop `.options enable_bidirectional=1` or use \
             fc_dcoupler for that path."
        );
        self.wpc = wpc;
        // 4 optical bundle ports (wpc·N each) + N control wires + 1 return.
        let stride = 4 * wpc + 1;
        let len = terminals.len();
        assert!(
            len > 1 && (len - 1).is_multiple_of(stride),
            "fc_optical_2x2: terminal count must be 4·wpc·N + N + 1 = {stride}·N + 1 \
             (wpc={wpc}) — four optical bundle ports, then N control wires, then one \
             control return; got {len}. A mismatch here usually means the optical and \
             control buses were declared with different channel counts."
        );
        let n = (len - 1) / stride;
        self.n_channels = n;
        self.nodes = terminals.to_vec();
        self.branches = vec![None; BRANCHES_PER_CHANNEL * n];
        self.lam_src.resize(n);
        // Size the per-channel state, then replay params captured before N was
        // known (parameters are set on the model card, ahead of instancing).
        self.w0 = vec![0.0; n];
        self.dw_dv = vec![0.0; n];
        self.explicit = vec![false; n];
        self.s = vec![Self::matrix_for_weight(0.0); n];
        let pend = std::mem::take(&mut self.pend);
        for (name, value) in &pend {
            self.apply_param(name, *value);
        }
        self.pend = pend;
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
        let lower = name.to_lowercase();
        match lower.as_str() {
            "il_db" | "alpha_db" => {
                // power dB → field amplitude factor
                self.il_amp = 10f64.powf(-value / 20.0);
                return true;
            }
            "tau_s" | "tau" | "latency_s" => {
                self.tau_s = value.max(0.0);
                return true;
            }
            "allow_gain" => {
                self.allow_gain = value != 0.0;
                return true;
            }
            _ => {}
        }
        // Per-channel params: N isn't known until setup_instance, so record
        // them and replay there. Recognise the name now so a typo still errors.
        let (base, _) = Self::split_index(&lower);
        let known = matches!(base, "w" | "dw_dv")
            || matches!(
                base,
                "s11_mag"
                    | "s11_deg"
                    | "s12_mag"
                    | "s12_deg"
                    | "s21_mag"
                    | "s21_deg"
                    | "s22_mag"
                    | "s22_deg"
            );
        if !known {
            return false;
        }
        self.pend.push((lower.clone(), value));
        if self.n_channels > 0 {
            self.apply_param(&lower, value);
        }
        true
    }

    fn eval(&mut self, x: &[f64], flags: EvalFlags, ctx: &SimContext) {
        // Latch which input port carries the λ tag; see LambdaSelect.
        for k in 0..self.n_channels {
            let (_, _, in1_l, ..) = self.channel_nodes(k);
            let in2_l = self.channel_lambda_in2(k);
            self.lam_src.observe(k, in1_l, in2_l, x);
        }
        // Passivity guard, once, on the first eval — the registry applies
        // instance params *after* setup_instance, so this is the earliest point
        // at which the matrix is final. Weight mode is unitary by construction
        // and its clamp keeps it so under any control voltage, hence explicit
        // matrices only.
        if !self.checked {
            self.checked = true;
            if !self.allow_gain {
                for k in 0..self.n_channels {
                    if self.explicit[k] {
                        let sigma = Self::sigma_max(&self.s[k]) * self.il_amp;
                        assert!(
                            sigma <= 1.0 + 1e-9,
                            "fc_optical_2x2: channel {k} matrix has gain (largest \
                             singular value {sigma:.6} > 1). Set allow_gain=1 if that \
                             is deliberate — otherwise it diverges silently in a \
                             feedback path."
                        );
                    }
                }
            }
        }
        let delay_active = flags.transient && self.tau_s > 0.0;
        self.refresh(x, delay_active, ctx);
    }

    fn load_residual(&self, b: &mut [f64]) {
        if !self.delay.is_active() {
            return;
        }
        // Delay mode: the outputs are driven by the history-reconstructed
        // inputs on the RHS instead of by live couplings (see load_jacobian).
        for k in 0..self.n_channels {
            let m = &self.s[k];
            let d = &self.delayed[DELAY_VALS_PER_CHANNEL * k..];
            let (a1_re, a1_im, a2_re, a2_im) = (d[0], d[1], d[2], d[3]);
            let out = |e1: usize, e2: usize| {
                (
                    m[e1].0 * a1_re - m[e1].1 * a1_im + m[e2].0 * a2_re - m[e2].1 * a2_im,
                    m[e1].1 * a1_re + m[e1].0 * a1_im + m[e2].1 * a2_re + m[e2].0 * a2_im,
                )
            };
            let (thru_re, thru_im) = out(0, 1);
            let (drop_re, drop_im) = out(2, 3);
            let base = BRANCHES_PER_CHANNEL * k;
            for (slot, v) in [
                (base, thru_re),
                (base + 1, thru_im),
                (base + 2, drop_re),
                (base + 3, drop_im),
            ] {
                if let Some(j) = self.branches[slot] {
                    b[j] += self.il_amp * v;
                }
            }
        }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let delay_active = self.delay.is_active();
        for k in 0..self.n_channels {
            let (in1_re, in1_im, in1_l, in2_re, in2_im, t_re, t_im, t_l, d_re, d_im, d_l) =
                self.channel_nodes(k);
            let m = &self.s[k];
            let g = self.il_amp;
            let base = BRANCHES_PER_CHANNEL * k;
            // out_re = Re(s_a)·in1_re − Im(s_a)·in1_im + Re(s_b)·in2_re − Im(s_b)·in2_im
            // out_im = Im(s_a)·in1_re + Re(s_a)·in1_im + Im(s_b)·in2_re + Re(s_b)·in2_im
            // stamp_potential_eq writes V(out) − Σk·V(in) = 0, so the couplings
            // carry a minus sign. In delay mode they drop out entirely and the
            // RHS history source drives the output (load_residual).
            let rows: [StampRow; 4] = [
                (
                    base,
                    t_re,
                    [
                        (in1_re, -g * m[0].0),
                        (in1_im, g * m[0].1),
                        (in2_re, -g * m[1].0),
                        (in2_im, g * m[1].1),
                    ],
                ),
                (
                    base + 1,
                    t_im,
                    [
                        (in1_re, -g * m[0].1),
                        (in1_im, -g * m[0].0),
                        (in2_re, -g * m[1].1),
                        (in2_im, -g * m[1].0),
                    ],
                ),
                (
                    base + 2,
                    d_re,
                    [
                        (in1_re, -g * m[2].0),
                        (in1_im, g * m[2].1),
                        (in2_re, -g * m[3].0),
                        (in2_im, g * m[3].1),
                    ],
                ),
                (
                    base + 3,
                    d_im,
                    [
                        (in1_re, -g * m[2].1),
                        (in1_im, -g * m[2].0),
                        (in2_re, -g * m[3].1),
                        (in2_im, -g * m[3].0),
                    ],
                ),
            ];
            for (slot, out, ins) in rows {
                if delay_active {
                    stamp_potential_eq(mat, &self.branches, slot, out, &[]);
                } else {
                    stamp_potential_eq(mat, &self.branches, slot, out, &ins);
                }
            }
            // λ labels pass through from whichever input is lit — in1 when
            // both are (a wavelength tag is not delayed).  Latched, so a
            // closed loop keeps one driver and drop_λ never binds to itself.
            // Same rule as fc_dcoupler; see LambdaSelect.
            let in2_l = self.channel_lambda_in2(k);
            let src_l = self.lam_src.pick(k, in1_l, in2_l);
            stamp_potential_eq(mat, &self.branches, base + 4, t_l, &[(src_l, -1.0)]);
            stamp_potential_eq(mat, &self.branches, base + 5, d_l, &[(src_l, -1.0)]);
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }

    /// `refresh` rebuilds `self.s[k]` from the control voltage and the stamp
    /// then treats it as a constant, so nothing carries `∂S/∂v_ctrl` — the
    /// `dw_dv_<k>` coupling is invisible to the adjoint unless declared.  Only
    /// the channels that actually take a control voltage: an `explicit[k]`
    /// channel has a matrix pinned by parameters and no voltage dependence.
    /// Both inputs reach both outputs — same rule as `fc_dcoupler`, and for the
    /// same reason: a ring fed only through its add port must still carry a
    /// label. Structural, so nothing latches.
    fn lambda_routing(&self) -> Vec<(usize, usize)> {
        let (wpc, n) = (self.wpc, self.n_channels);
        let lam = wpc - 1;
        let (p2, p3, p4) = (wpc * n, 2 * wpc * n, 3 * wpc * n);
        (0..n)
            .flat_map(|k| {
                let (i1, i2) = (wpc * k + lam, p2 + wpc * k + lam);
                let (t, d) = (p3 + wpc * k + lam, p4 + wpc * k + lam);
                [(i1, t), (i1, d), (i2, t), (i2, d)]
            })
            .collect()
    }

    fn frozen_jacobian_columns(&self) -> Vec<usize> {
        let ctl_base = 4 * self.wpc * self.n_channels;
        let mut cols: Vec<usize> = (0..self.n_channels)
            .filter(|&k| !self.explicit[k] && self.dw_dv[k] != 0.0)
            .filter_map(|k| self.nodes[ctl_base + k])
            .collect();
        if !cols.is_empty() {
            cols.extend(self.nodes[ctl_base + self.n_channels]);
        }
        cols
    }

    fn commit_timestep(&mut self, x: &[f64]) {
        if !self.delay.is_active() {
            return;
        }
        let read = |nid: NodeId| nid.map_or(0.0, |i| x[i]);
        let mut snap = vec![0.0; self.n_channels * DELAY_VALS_PER_CHANNEL];
        for k in 0..self.n_channels {
            let (in1_re, in1_im, _, in2_re, in2_im, ..) = self.channel_nodes(k);
            let s = DELAY_VALS_PER_CHANNEL * k;
            snap[s] = read(in1_re);
            snap[s + 1] = read(in1_im);
            snap[s + 2] = read(in2_re);
            snap[s + 3] = read(in2_im);
        }
        self.delay.record(snap, self.tau_s);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// θ = ½·acos(−w) must give P_drop − P_thru = w for any w, with the pair
    /// summing to 1 — that identity is the whole reason `w` is a usable knob.
    #[test]
    fn weight_matrix_is_a_lossless_bipolar_split() {
        for w in [-1.0, -0.6, 0.0, 0.25, 1.0] {
            let m = NativeOptical2x2::matrix_for_weight(w);
            let p_thru = m[0].0 * m[0].0 + m[0].1 * m[0].1;
            let p_drop = m[2].0 * m[2].0 + m[2].1 * m[2].1;
            assert!(
                (p_thru + p_drop - 1.0).abs() < 1e-12,
                "w={w}: power {p_thru}+{p_drop}"
            );
            assert!(
                (p_drop - p_thru - w).abs() < 1e-12,
                "w={w}: got {}",
                p_drop - p_thru
            );
            assert!(
                (NativeOptical2x2::sigma_max(&m) - 1.0).abs() < 1e-12,
                "w={w} should be unitary"
            );
        }
    }

    #[test]
    fn sigma_max_detects_gain() {
        // A diagonal matrix with a 2× entry has σ_max = 2.
        let m = [(2.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.5, 0.0)];
        assert!((NativeOptical2x2::sigma_max(&m) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn split_index_only_takes_trailing_digits() {
        assert_eq!(NativeOptical2x2::split_index("w_3"), ("w", Some(3)));
        assert_eq!(NativeOptical2x2::split_index("w"), ("w", None));
        assert_eq!(NativeOptical2x2::split_index("dw_dv"), ("dw_dv", None));
        assert_eq!(NativeOptical2x2::split_index("dw_dv_2"), ("dw_dv", Some(2)));
        assert_eq!(NativeOptical2x2::split_index("s12_mag"), ("s12_mag", None));
    }
}
