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
//! `docs/photonic-models.md` documents the devices this builds, and the user
//! guide's "Writing custom devices" section walks through adding one.

use super::{dB_per_cm_to_neper_per_m, n_eff_at_lambda, stamp_potential_eq, C0};
use crate::delay::DelayLine;
use crate::device::{NodeId, SimContext};
use crate::mna::MnaMatrix;

/// A quantity a drive reports for the optical bundle: usually one value that
/// applies to every channel, sometimes one per channel.
///
/// `OpticalPerturbation` used to carry plain `f64`s, which made a whole class of
/// physics inexpressible — anything depending on a channel's own wavelength or
/// field. The concrete casualty was multi-wavelength TPA: cross-TPA between
/// distinct frequencies is twice self-TPA, so the loss on channel `j` is
/// `β/A_eff·(I_j + 2·Σ_{k≠j} I_k)`, which is not one number.
///
/// `Uniform` exists so the drives whose effect genuinely is shared pay no
/// allocation for the generality. This is rebuilt on every Newton iteration,
/// once more per control node per finite-difference probe.
#[derive(Clone, Debug, PartialEq)]
pub enum PerChannel {
    /// One value, every channel.
    Uniform(f64),
    /// One value per channel, in channel order.
    Each(Vec<f64>),
}

impl PerChannel {
    /// Channel `k`'s value. Out of range reads 0 rather than panicking: a drive
    /// returning a short vector has a bug, but a wrong number beats a crash
    /// mid-solve, and the segment's own channel count is the authority.
    pub fn at(&self, k: usize) -> f64 {
        match self {
            PerChannel::Uniform(v) => *v,
            PerChannel::Each(v) => v.get(k).copied().unwrap_or(0.0),
        }
    }

    /// Whether every channel shares one value.
    pub fn is_uniform(&self) -> bool {
        matches!(self, PerChannel::Uniform(_))
    }

    fn combine(&self, other: &Self, f: impl Fn(f64, f64) -> f64) -> Self {
        match (self, other) {
            (PerChannel::Uniform(a), PerChannel::Uniform(b)) => PerChannel::Uniform(f(*a, *b)),
            (PerChannel::Each(a), PerChannel::Uniform(b)) => {
                PerChannel::Each(a.iter().map(|x| f(*x, *b)).collect())
            }
            (PerChannel::Uniform(a), PerChannel::Each(b)) => {
                PerChannel::Each(b.iter().map(|y| f(*a, *y)).collect())
            }
            (PerChannel::Each(a), PerChannel::Each(b)) => {
                PerChannel::Each(a.iter().zip(b).map(|(x, y)| f(*x, *y)).collect())
            }
        }
    }

    /// Elementwise sum — how a composite drive adds its two halves.
    pub fn add(&self, other: &Self) -> Self {
        self.combine(other, |a, b| a + b)
    }

    /// `(self − other) / h`: the finite difference the segment's Newton
    /// cross-terms are built from.
    pub fn diff_by(&self, other: &Self, h: f64) -> Self {
        self.combine(other, |a, b| (a - b) / h)
    }

    /// Zero on every channel.
    pub fn zero() -> Self {
        PerChannel::Uniform(0.0)
    }
}

impl Default for PerChannel {
    fn default() -> Self {
        PerChannel::Uniform(0.0)
    }
}

impl From<f64> for PerChannel {
    fn from(v: f64) -> Self {
        PerChannel::Uniform(v)
    }
}

impl From<Vec<f64>> for PerChannel {
    fn from(v: Vec<f64>) -> Self {
        PerChannel::Each(v)
    }
}

/// Input couplings for a `stamp_potential_eq` row: `(node, coefficient)` pairs.
/// Empty in delay mode (the output is driven by a history source on the RHS).
type Couplings<'a> = &'a [(NodeId, f64)];

/// The optical core: bundle propagation + perturbation + delay + stamping,
/// shared by the waveguide and every active phase-shifter / modulator class.
///
/// `Clone` so a composite device can cut `N` identical slices from one
/// configured template rather than repeating its parameter handling `N` times
/// (`fc_tw_ps`). A clone copies the delay history too, which is what you want
/// at construction — it is empty — and is why cloning a *running* segment is
/// not something any caller does.
#[derive(Clone)]
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
    /// Each channel's resolved wavelength (m), handed down by
    /// [`Device::set_resolved_lambda`](crate::device::Device::set_resolved_lambda)
    /// at build time. Zero means "never told", which reads as `wl_ref_m` — the
    /// same answer an undriven λ wire used to bootstrap to.
    lambda_m: Vec<f64>,
    n_channels: usize,
    /// Smallest optical wire count that would have been usable, recorded when
    /// `setup_instance` is handed one it cannot use. Without it
    /// `num_optical_wires()` would report the unconfigured segment's 0 and the
    /// caller's error would quote a nonsense expectation.
    min_wires: Option<usize>,
    wpc: usize, // wires per channel: 3 (unidir) or 5 (bidir)
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
    c_cached: Vec<f64>,
    s_cached: Vec<f64>,
    // ── Newton cross-terms for voltage-dependent coefficients ────────────────
    // c and s depend on the drive's control voltages, so `out = c·in + s·in`
    // is BILINEAR in the unknowns. Stamping c/s frozen (their value at the
    // previous iterate) turns Newton into successive substitution for this
    // coupling: it converges only while the loop's iteration map contracts, and
    // silently fails once a feedback path pushes the spectral radius past 1 —
    // exactly what an unstable fixed point in a photonic recurrent network is.
    // Stamping ∂(out)/∂V restores a true Newton, which does not care about
    // stability. See `set_control_sensitivity`.
    ctrl_nodes: Vec<NodeId>,
    /// Per (channel, control-node): ∂c/∂V and ∂s/∂V at the current iterate.
    dc_dv: Vec<f64>,
    ds_dv: Vec<f64>,
    /// The control-node voltages at the current iterate (for the RHS offset).
    ctrl_v0: Vec<f64>,
    /// The (re, im) source values that c and s multiply in the CURRENT mode —
    /// the live inputs normally, the delay-line history in delay mode.
    src_re: Vec<f64>,
    src_im: Vec<f64>,
    // ── Group-delay line (engaged only when the owning device asks) ──────────
    delay: DelayLine,
    delayed: Vec<f64>,
    /// Whether this eval stamps the delayed form: couplings dropped, output
    /// driven by the reconstructed history (or, in `.ac`, by `ac_stamps`).
    /// False when the delay is engaged but there is no history yet — see
    /// `refresh_with_sens`.
    delayed_stamp: bool,
    /// `SimOptions::waveguide_delay`, cached at `setup_model`. The owning
    /// device decides per-eval whether the delay is engaged; this is the
    /// run-level setting, and it is what the step controller can ask about
    /// before any eval has happened.
    delay_option: bool,
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
            lambda_m: Vec::new(),
            n_channels: 0,
            min_wires: None,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
            c_cached: Vec::new(),
            s_cached: Vec::new(),
            ctrl_nodes: Vec::new(),
            dc_dv: Vec::new(),
            ds_dv: Vec::new(),
            ctrl_v0: Vec::new(),
            src_re: Vec::new(),
            src_im: Vec::new(),
            delay: DelayLine::new(),
            delayed: Vec::new(),
            delayed_stamp: false,
            delay_option: false,
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

    /// Auxiliary branch rows per channel: one per *driven field* wire.
    ///
    /// One fewer than the wire count, because the λ wire is not driven — a
    /// wavelength is resolved before the solve, so there is no equation to
    /// stamp and no row to stamp it into. Getting this wrong by one leaves an
    /// empty branch row, which `gmin` makes non-singular and therefore silent.
    fn bpc(wpc: usize) -> usize {
        wpc - 1
    }

    pub fn refresh_tau(&mut self) {
        self.tau_g_s = self.length_m * self.n_g / C0;
    }

    /// Pick up the band-centre defaults from the context (call from
    /// `Device::setup_model`). Defaults the dispersion reference to the band
    /// centre too, unless the user already overrode `wl_ref_m`.
    pub fn setup_model(&mut self, ctx: &SimContext) {
        if (self.wl_ref_m - 1.55e-6).abs() < 1e-12 {
            self.wl_ref_m = ctx.lambda_center_m;
        }
        self.wpc = ctx.wires_per_channel();
        self.delay_option = ctx.waveguide_delay;
        self.refresh_tau();
    }

    /// The timestep this segment needs, or `None` when its delay is not
    /// engaged.
    ///
    /// Answered from the option and the geometry rather than from
    /// `delay.is_active()`, because the step controller asks before the first
    /// `eval` — and a bound that only appears after the first step is a bound
    /// that missed the first step. See `Device::requested_max_timestep`.
    pub fn requested_max_timestep(&self) -> Option<f64> {
        (self.delay_option && self.tau_g_s > 0.0).then_some(self.tau_g_s / 2.0)
    }

    /// Bind the segment to its `2·wpc·N` optical bundle wires (in block then out
    /// block). Allocates the per-channel branch slots and the c/s caches.
    pub fn setup_instance(&mut self, optical_terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 2 * wpc; // in + out

        // A mis-wired instance used to `assert!` here, which crosses the pyo3
        // boundary as a PanicException naming no element, and aborts the CLI
        // with a backtrace instead of a diagnosis. Leave the segment
        // unconfigured instead: the count is checked against `num_terminals()`
        // in `build_devices_with_footprints`, which can name the element, the
        // model and both counts. `ActiveOpticalDevice` guards before reaching
        // here; passive segment devices (the waveguide) rely on this.
        //
        // The common way to land here is an optical net that was never given
        // its own `.optical_port`, so it stays one scalar wire instead of
        // expanding to `wpc`.
        if optical_terminals.is_empty() || !optical_terminals.len().is_multiple_of(stride) {
            self.min_wires = Some(stride);
            self.n_channels = 0;
            self.nodes.clear();
            self.branches.clear();
            self.lambda_m.clear();
            return;
        }
        self.min_wires = None;
        let n = optical_terminals.len() / stride;
        self.n_channels = n;
        self.nodes = optical_terminals.to_vec();
        self.branches = vec![None; Self::bpc(wpc) * n];
        self.lambda_m = vec![0.0; n];
        self.c_cached = vec![1.0; n];
        self.s_cached = vec![0.0; n];
        let nc = self.ctrl_nodes.len();
        self.dc_dv = vec![0.0; n * nc];
        self.ds_dv = vec![0.0; n * nc];
        self.ctrl_v0 = vec![0.0; nc];
        self.src_re = vec![0.0; n];
        self.src_im = vec![0.0; n];
    }

    /// Declare the drive's electrical terminals, so the segment can stamp the
    /// Newton cross-terms ∂(optical output)/∂(control voltage).
    pub fn set_control_nodes(&mut self, nodes: &[NodeId]) {
        self.ctrl_nodes = nodes.to_vec();
        let n = self.n_channels.max(1);
        self.dc_dv = vec![0.0; n * nodes.len()];
        self.ds_dv = vec![0.0; n * nodes.len()];
        self.ctrl_v0 = vec![0.0; nodes.len()];
    }

    /// `(from, to)` terminal pairs carrying a wavelength label through this
    /// segment: channel `k`'s input λ wire to its output λ wire.
    ///
    /// Terminal numbering is the segment's own, which for every device built on
    /// one is also the device's — optical bundle wires come first, electrical
    /// last. λ sits at `wpc - 1` within a channel, and the output block starts
    /// at `wpc · n_channels`.
    pub fn lambda_routing(&self) -> Vec<(usize, usize)> {
        let lam = self.wpc - 1;
        let out_base = self.wpc * self.n_channels;
        (0..self.n_channels)
            .map(|k| (self.wpc * k + lam, out_base + self.wpc * k + lam))
            .collect()
    }

    /// Number of optical bundle wires this segment occupies (`2·wpc·N`).
    ///
    /// A refused setup leaves the segment with zero channels, so report the
    /// smallest count that would have worked rather than a misleading 0.
    pub fn num_optical_wires(&self) -> usize {
        if let Some(min) = self.min_wires {
            return min;
        }
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

    /// Take the per-channel λ resolved for this instance's input port.
    ///
    /// Terminal layout: λ sits at `wpc - 1` within a channel, and the input
    /// block comes first. A segment's output λ equals its input λ by
    /// construction (that is what `lambda_routing` declares), so reading the
    /// input side is not a choice between two answers.
    pub fn set_resolved_lambda(&mut self, per_terminal: &[f64]) {
        let lam = self.wpc - 1;
        for k in 0..self.n_channels {
            if let Some(&v) = per_terminal.get(self.wpc * k + lam) {
                self.lambda_m[k] = v;
            }
        }
    }

    /// Channel `k`'s wavelength (m), or `wl_ref_m` if resolution never reached
    /// this instance — the same value an undriven λ wire bootstrapped to.
    fn lambda_of(&self, k: usize) -> f64 {
        match self.lambda_m.get(k) {
            Some(&v) if v > 0.0 => v,
            _ => self.wl_ref_m,
        }
    }

    /// Per-channel input optical intensity |A_in|² (W), the input to any
    /// photoconductive / detection back-action.
    ///
    /// Counts BOTH propagation directions. Under `enable_bidirectional` the
    /// backward field enters at the far end of the segment, so its input wires
    /// live in the `out` block (`out_base + wpc·k + 2/3`) — see `stamp`, where
    /// the backward equation reads `out_*_bw` to drive `in_*_bw`. Absorption
    /// does not care which way a photon travels, so a device that heats or
    /// generates carriers from absorbed power must see both; reading the
    /// forward wires only made backward light heat nothing at all.
    pub fn channel_intensities(&self, x: &[f64]) -> Vec<f64> {
        let read = |nid: NodeId| nid.map_or(0.0, |i| x[i]);
        let mag2 = |re: NodeId, im: NodeId| {
            let (re, im) = (read(re), read(im));
            re * re + im * im
        };
        let out_base = self.wpc * self.n_channels;
        (0..self.n_channels)
            .map(|k| {
                let fw = mag2(self.nodes[self.wpc * k], self.nodes[self.wpc * k + 1]);
                let bw = if self.wpc == 5 {
                    mag2(
                        self.nodes[out_base + self.wpc * k + 2],
                        self.nodes[out_base + self.wpc * k + 3],
                    )
                } else {
                    0.0
                };
                fw + bw
            })
            .collect()
    }

    /// Recompute the cached transmission coefficients for every channel from
    /// the segment's own propagation plus the supplied per-segment perturbation:
    /// `Δn_eff` (added to the effective index — a λ-*dependent* phase
    /// `2π·Δn_eff·L/λ`, e.g. free-carrier dispersion), `Δφ` (a λ-*independent*
    /// direct phase added to every channel, e.g. a calibrated thermo-optic
    /// `φ_th = π·P/P_π`), and `Δα` (added to the loss, Neper/m). A passive
    /// segment passes all three as 0.
    ///
    /// `delay_active` engages the group-delay line (the owning device decides
    /// when — e.g. `flags.transient && ctx.waveguide_delay && τ>0`).
    pub fn refresh(
        &mut self,
        x: &[f64],
        dn_eff: PerChannel,
        dphi: f64,
        dalpha_neper_m: PerChannel,
        delay_active: bool,
        ctx: &SimContext,
    ) {
        self.refresh_with_sens(x, dn_eff, dphi, dalpha_neper_m, delay_active, ctx, &[]);
    }

    /// As [`refresh`](Self::refresh), plus the drive's `(d dn_eff/dV, d dphi/dV,
    /// d dalpha/dV)` per control node. Pass an empty slice to skip the Newton
    /// cross-terms (correct for a passive segment, or a drive with no
    /// voltage dependence).
    #[allow(clippy::too_many_arguments)]
    pub fn refresh_with_sens(
        &mut self,
        x: &[f64],
        dn_eff: PerChannel,
        dphi: f64,
        dalpha_neper_m: PerChannel,
        delay_active: bool,
        ctx: &SimContext,
        dpert_dv: &[(PerChannel, f64, PerChannel)],
    ) {
        self.delay.set_state(delay_active, ctx.time_s);
        // Engaged is not the same as "reconstruct from history". Before the
        // first accepted step there is nothing to reconstruct from, and a
        // delayed source read out of an empty buffer is zero — the segment
        // would extinguish its own output for one step and a resonator would
        // ring on the artefact. With no past, the steady-state transfer is the
        // honest answer, and one step later there is a past.
        //
        // `.ac` reaches here under transient flags too, has no history by
        // construction, and does want the delayed form: its delay is the exact
        // `exp(−jΩτ_g)` coupling in `ac_stamps`, which replaces the couplings
        // dropped below. `ctx.discretisation` tells the two apart, because only
        // a time-domain step sets one.
        let time_domain = ctx.discretisation.is_some();
        self.delayed_stamp = delay_active && (!time_domain || !self.delay.is_empty());
        if self.delayed_stamp {
            let width = self.n_channels * self.vals_per_channel();
            self.delayed = self.delay.sample(self.tau_g_s, width);
        }

        let two_pi = 2.0 * std::f64::consts::PI;
        // t_amp is per channel now: Δα may differ by channel (multi-λ TPA), so
        // the amplitude factor cannot be hoisted out of the loop.
        let phi_ref = if self.pin_at_ref {
            two_pi * self.n_eff * self.length_m / self.wl_ref_m
        } else {
            0.0
        };
        let nc = self.ctrl_nodes.len();
        let want_sens = nc > 0 && dpert_dv.len() == nc;
        if want_sens {
            for (i, n) in self.ctrl_nodes.iter().enumerate() {
                self.ctrl_v0[i] = n.map_or(0.0, |j| x[j]);
            }
        }
        let per = self.vals_per_channel();
        for k in 0..self.n_channels {
            let lambda = self.lambda_of(k);
            let n_eff_lam = n_eff_at_lambda(self.n_eff, self.n_g, self.wl_ref_m, lambda);
            // Δn_eff folds into the index (φ_eo = 2π·Δn_eff·L/λ); Δφ adds a
            // wavelength-independent rotation (heater, Pockels-with-fixed-gap…).
            let t_amp = (-(self.alpha_neper_m + dalpha_neper_m.at(k)) * self.length_m / 2.0).exp();
            let phi = two_pi * (n_eff_lam + dn_eff.at(k)) * self.length_m / lambda - phi_ref + dphi;
            let (c, s) = (t_amp * phi.cos(), t_amp * phi.sin());
            self.c_cached[k] = c;
            self.s_cached[k] = s;
            if !want_sens {
                continue;
            }
            // The (re, im) pair that c and s multiply in the current mode.
            // The field the coefficient's *derivative* multiplies, which is a
            // different question from where the output comes from.
            //
            // In the time domain that is the delayed field: the light being
            // modulated now entered a group delay ago. In the frequency domain
            // there is no history to read — `.ac` linearises about an operating
            // point — and reading the empty buffer would put a zero here and
            // silently delete the modulation, leaving a segment that propagates
            // light and does not respond to its drive. The operating-point
            // field at the input is the right small-signal answer, and it is
            // frequency-independent by construction.
            let (sr, si) = if self.delayed_stamp && time_domain {
                (self.delayed[per * k], self.delayed[per * k + 1])
            } else {
                let read = |nid: NodeId| nid.map_or(0.0, |i| x[i]);
                (
                    read(self.nodes[self.wpc * k]),
                    read(self.nodes[self.wpc * k + 1]),
                )
            };
            self.src_re[k] = sr;
            self.src_im[k] = si;
            // dc/dV = (dt/dV)·cos φ − t·sin φ·(dφ/dV) = (dt/dV)/t·c − s·(dφ/dV)
            // ds/dV = (dt/dV)/t·s + c·(dφ/dV)
            let kphi = two_pi * self.length_m / lambda;
            for (i, (dn_dv, ddphi_dv, dalpha_dv)) in dpert_dv.iter().enumerate() {
                let (dn_dv, dalpha_dv) = (dn_dv.at(k), dalpha_dv.at(k));
                let dphi_dv = kphi * dn_dv + ddphi_dv;
                // t_amp = exp(−(α+Δα)·L/2) ⇒ dt/dV = −t·(L/2)·dΔα/dV, so
                // (dt/dV)/t is just −(L/2)·dΔα/dV and needs no division by t.
                let dlnt_dv = -0.5 * self.length_m * dalpha_dv;
                self.dc_dv[k * nc + i] = dlnt_dv * c - s * dphi_dv;
                self.ds_dv[k * nc + i] = dlnt_dv * s + c * dphi_dv;
            }
        }
    }

    /// Whether Newton cross-terms are armed for this segment.
    fn has_sens(&self) -> bool {
        !self.ctrl_nodes.is_empty() && self.dc_dv.len() == self.n_channels * self.ctrl_nodes.len()
    }

    /// Control columns whose `∂(out)/∂V` this segment did *not* stamp, because
    /// its drive model supplied no derivatives and the coefficient stayed
    /// frozen.  Empty when the cross-terms are armed, which is the usual case.
    /// See `Device::frozen_jacobian_columns`.
    pub fn unstamped_control_columns(&self) -> Vec<usize> {
        if self.has_sens() {
            Vec::new()
        } else {
            self.ctrl_nodes.iter().flatten().copied().collect()
        }
    }

    /// ∂(out_re)/∂V and ∂(out_im)/∂V for channel `k`, control node `i`.
    fn dout_dv(&self, k: usize, i: usize) -> (f64, f64) {
        let nc = self.ctrl_nodes.len();
        let (dc, ds) = (self.dc_dv[k * nc + i], self.ds_dv[k * nc + i]);
        let (sr, si) = (self.src_re[k], self.src_im[k]);
        (dc * sr + ds * si, -ds * sr + dc * si)
    }

    /// Stamp the per-channel optical branch equations. In delay mode the output
    /// ports are driven by a history-reconstructed source (see
    /// [`stamp_residual`](Self::stamp_residual)), so their couplings to the live
    /// input nodes are dropped. The delay state is the one set by the most
    /// recent [`refresh`](Self::refresh).
    pub fn stamp(&self, mat: &mut MnaMatrix) {
        let delay_active = self.delayed_stamp;
        let n = self.n_channels;
        let wpc = self.wpc;
        let bpc = Self::bpc(wpc);
        let out_base = wpc * n;
        for k in 0..n {
            let c = self.c_cached[k];
            let s = self.s_cached[k];
            let in_re_fw = self.nodes[wpc * k];
            let in_im_fw = self.nodes[wpc * k + 1];
            let out_re_fw = self.nodes[out_base + wpc * k];
            let out_im_fw = self.nodes[out_base + wpc * k + 1];
            let (re_ins, im_ins): (Couplings, Couplings) = if delay_active {
                (&[], &[])
            } else {
                (
                    &[(in_re_fw, -c), (in_im_fw, -s)],
                    &[(in_re_fw, s), (in_im_fw, -c)],
                )
            };
            stamp_potential_eq(mat, &self.branches, bpc * k, out_re_fw, re_ins);
            stamp_potential_eq(mat, &self.branches, bpc * k + 1, out_im_fw, im_ins);
            // Newton cross-terms: c and s are functions of the control voltage,
            // so the output equation is bilinear. Without these the coupling is
            // stamped frozen and Newton degenerates into successive
            // substitution, which cannot reach a fixed point whose iteration
            // map expands (an unstable operating point in a feedback loop).
            if self.has_sens() {
                for i in 0..self.ctrl_nodes.len() {
                    let Some(cn) = self.ctrl_nodes[i] else {
                        continue;
                    };
                    let (g_re, g_im) = self.dout_dv(k, i);
                    if let Some(j) = self.branches[bpc * k] {
                        mat.a[j][cn] -= g_re;
                    }
                    if let Some(j) = self.branches[bpc * k + 1] {
                        mat.a[j][cn] -= g_im;
                    }
                }
            }
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
                stamp_potential_eq(mat, &self.branches, bpc * k + 2, in_re_bw, bw_re_ins);
                stamp_potential_eq(mat, &self.branches, bpc * k + 3, in_im_bw, bw_im_ins);
            }
        }
    }

    /// Stamp the delay-line history source onto the RHS. No-op when the delay
    /// line is inactive (the instantaneous path is fully homogeneous).
    pub fn stamp_residual(&self, b: &mut [f64]) {
        // The Newton cross-terms carry an inhomogeneous offset −G·V0 (the
        // linearisation is about the previous iterate), independent of the
        // delay line.
        let bpc = Self::bpc(self.wpc);
        if self.has_sens() {
            let nc = self.ctrl_nodes.len();
            for k in 0..self.n_channels {
                for i in 0..nc {
                    if self.ctrl_nodes[i].is_none() {
                        continue;
                    }
                    let (g_re, g_im) = self.dout_dv(k, i);
                    let v0 = self.ctrl_v0[i];
                    if let Some(j) = self.branches[bpc * k] {
                        b[j] -= g_re * v0;
                    }
                    if let Some(j) = self.branches[bpc * k + 1] {
                        b[j] -= g_im * v0;
                    }
                }
            }
        }
        if !self.delayed_stamp {
            return;
        }
        let n = self.n_channels;
        let per = self.vals_per_channel();
        for k in 0..n {
            let c = self.c_cached[k];
            let s = self.s_cached[k];
            let dly_fw_re = self.delayed[per * k];
            let dly_fw_im = self.delayed[per * k + 1];
            if let Some(j) = self.branches[bpc * k] {
                b[j] += c * dly_fw_re + s * dly_fw_im;
            }
            if let Some(j) = self.branches[bpc * k + 1] {
                b[j] += -s * dly_fw_re + c * dly_fw_im;
            }
            if self.wpc == 5 {
                let dly_bw_re = self.delayed[per * k + 2];
                let dly_bw_im = self.delayed[per * k + 3];
                if let Some(j) = self.branches[bpc * k + 2] {
                    b[j] += c * dly_bw_re + s * dly_bw_im;
                }
                if let Some(j) = self.branches[bpc * k + 3] {
                    b[j] += -s * dly_bw_re + c * dly_bw_im;
                }
            }
        }
    }

    /// The envelope delay as a per-frequency complex coupling, `exp(−jΩτ_g)`
    /// times the couplings [`stamp`](Self::stamp) drops in delay mode.
    ///
    /// The optical unknowns are envelope quadratures, so in a small-signal
    /// sweep each of `in_re` and `in_im` is itself a phasor at the *modulation*
    /// frequency `Ω`. The propagation phase and loss are already in `c` and `s`
    /// — they are the carrier's business and do not move with `Ω`. What the
    /// delay adds is `out(Ω) = (c − js)·exp(−jΩτ_g)·in(Ω)`, and that factor is
    /// the whole content of this stamp.
    ///
    /// Empty unless the delay is engaged, in which case the instantaneous
    /// couplings are in `G` already and there is nothing to add.
    pub fn ac_stamps(&self, omega: f64) -> Vec<crate::device::AcStamp> {
        use crate::device::AcStamp;
        if !self.delayed_stamp {
            return Vec::new();
        }
        let (qr, qi) = ((omega * self.tau_g_s).cos(), -(omega * self.tau_g_s).sin());
        let (n, wpc) = (self.n_channels, self.wpc);
        let bpc = Self::bpc(wpc);
        let out_base = wpc * n;
        let mut out = Vec::with_capacity(4 * n);
        let push = |row: Option<usize>, col: NodeId, k: f64, sink: &mut Vec<AcStamp>| {
            if let (Some(r), Some(c)) = (row, col) {
                sink.push(AcStamp {
                    row: r,
                    col: c,
                    re: k * qr,
                    im: k * qi,
                });
            }
        };
        for k in 0..n {
            let (c, s) = (self.c_cached[k], self.s_cached[k]);
            let (in_re, in_im) = (self.nodes[wpc * k], self.nodes[wpc * k + 1]);
            let (r_re, r_im) = (self.branches[bpc * k], self.branches[bpc * k + 1]);
            push(r_re, in_re, -c, &mut out);
            push(r_re, in_im, -s, &mut out);
            push(r_im, in_re, s, &mut out);
            push(r_im, in_im, -c, &mut out);
            if wpc == 5 {
                // Backward direction: the input is the *output* block's wire,
                // matching the couplings `stamp` drops for it.
                let (bw_re, bw_im) = (
                    self.nodes[out_base + wpc * k + 2],
                    self.nodes[out_base + wpc * k + 3],
                );
                let (rb_re, rb_im) = (self.branches[bpc * k + 2], self.branches[bpc * k + 3]);
                push(rb_re, bw_re, -c, &mut out);
                push(rb_re, bw_im, -s, &mut out);
                push(rb_im, bw_re, s, &mut out);
                push(rb_im, bw_im, -c, &mut out);
            }
        }
        out
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
        seg.setup_instance(
            &[Some(0), Some(1), Some(2), Some(3), Some(4), Some(5)],
            &ctx,
        );
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
        seg.refresh(
            &x,
            PerChannel::zero(),
            0.0,
            PerChannel::zero(),
            false,
            &SimContext::default(),
        );
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
            seg.refresh(
                &x,
                PerChannel::zero(),
                0.0,
                PerChannel::zero(),
                true,
                &ctx_delay(t),
            );
            assert!(seg.delay_active());
            seg.commit(&x);
        }
        // Query at t = 3 ⇒ delayed source is the input at t − τ = 1.
        let xq = vec![3.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        seg.refresh(
            &xq,
            PerChannel::zero(),
            0.0,
            PerChannel::zero(),
            true,
            &ctx_delay(3.0),
        );
        assert!(
            (seg.delayed()[0] - 1.0).abs() < 1e-12,
            "got {}",
            seg.delayed()[0]
        );
        // Linear interpolation between samples: t − τ = 1.5 ⇒ 1.5.
        seg.refresh(
            &xq,
            PerChannel::zero(),
            0.0,
            PerChannel::zero(),
            true,
            &ctx_delay(3.5),
        );
        assert!(
            (seg.delayed()[0] - 1.5).abs() < 1e-12,
            "got {}",
            seg.delayed()[0]
        );
    }

    #[test]
    fn delay_preserves_energy_and_stamps_residual() {
        let mut seg = lossless_segment(1.0);
        let x0 = vec![0.6, 0.8, 0.0, 0.0, 0.0, 0.0]; // |A| = 1
        seg.refresh(
            &x0,
            PerChannel::zero(),
            0.0,
            PerChannel::zero(),
            true,
            &ctx_delay(0.0),
        );
        seg.commit(&x0);
        let xq = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        seg.refresh(
            &xq,
            PerChannel::zero(),
            0.0,
            PerChannel::zero(),
            true,
            &ctx_delay(1.0),
        );
        let mut b = vec![0.0; 9];
        seg.stamp_residual(&mut b);
        let out_mag2 = b[6] * b[6] + b[7] * b[7];
        assert!(
            (out_mag2 - 1.0).abs() < 1e-9,
            "lossless delay must conserve |A|²=1, got {out_mag2}"
        );
    }
}
