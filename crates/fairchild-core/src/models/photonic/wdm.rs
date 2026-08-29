use super::spectrum::{AwgSpectrum, ChannelGrid};
use super::{stamp_potential_eq, C0};
use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

/// Optional per-channel spectral response shared by `fc_mux` and `fc_demux`.
///
/// Both devices default to a lossless identity route, so an existing netlist is
/// unaffected. Setting any parameter turns on a **diagonal** filter: channel `k`
/// is multiplied by its own passband, evaluated at whatever wavelength that
/// channel actually carries.
///
/// Diagonal is the whole story for a mux, and deliberately incomplete for a
/// demux:
///
/// - A **mux** combines N inputs onto one fibre. Each input lands in its own
///   channel slot, so there is nowhere for cross-channel leakage to *go* — the
///   only imperfections are per-channel insertion loss and the penalty a
///   detuned laser pays on the skirt. Both are modelled here.
/// - A **demux** physically does leak channel `k` into output port `j ≠ k`, and
///   that is **not** modelled here, because it cannot be: an output port of
///   `fc_demux` is a single channel, and representing the leak would mean
///   summing two different carriers into one complex envelope. For a demux with
///   crosstalk, use [`NativeAwgr`](super::NativeAwgr) with `N−1` input ports
///   left dark — physically the same device, and its output ports are N-channel
///   buses with somewhere for the leakage to live.
///
/// Parameters (all optional): `il_db`, `lambda0_nm`, `df_ghz` (alias
/// `spacing_ghz`), `fwhm_ghz` (alias `bw_ghz`), `shape_p`,
/// `dlambda_dt_pm_per_k`, `t_nom_k`. With only `il_db` set the loss is flat
/// across the band; adding `fwhm_ghz` gives each channel a passband.
#[derive(Clone)]
struct ChannelFilter {
    lambda0_m: f64,
    df_hz: f64,
    fwhm_hz: Option<f64>,
    shape_p: f64,
    il_db: f64,
    dlambda_dt_m_per_k: f64,
    t_nom_k: f64,
    /// Stays false until a parameter is set, which keeps the default identity
    /// route bit-for-bit rather than merely numerically close.
    active: bool,
}

impl ChannelFilter {
    fn new() -> Self {
        Self {
            lambda0_m: 1.55e-6,
            df_hz: 100e9,
            fwhm_hz: None,
            shape_p: 1.0,
            il_db: 0.0,
            dlambda_dt_m_per_k: 0.0,
            t_nom_k: 300.15,
            active: false,
        }
    }

    fn set(&mut self, name: &str, value: f64) -> bool {
        let hit = match name {
            "il_db" => {
                self.il_db = value;
                true
            }
            "lambda0_nm" => {
                self.lambda0_m = value * 1e-9;
                true
            }
            "df_ghz" | "spacing_ghz" => {
                self.df_hz = value * 1e9;
                true
            }
            // 0 means "no passband", matching `fc_awgr`. Taken literally it is a
            // zero-width Gaussian, which darkens every channel but the one
            // sitting exactly on the grid anchor.
            "fwhm_ghz" | "bw_ghz" => {
                self.fwhm_hz = (value > 0.0).then_some(value * 1e9);
                true
            }
            "shape_p" => {
                self.shape_p = value.max(0.1);
                true
            }
            "dlambda_dt_pm_per_k" => {
                self.dlambda_dt_m_per_k = value * 1e-12;
                true
            }
            "t_nom_k" => {
                self.t_nom_k = value;
                true
            }
            _ => false,
        };
        self.active |= hit;
        hit
    }

    /// Field transmission for channel `k` of an `n`-channel bundle at `lambda`.
    /// Exactly 1.0 while no parameter has been set.
    fn amp(&self, lambda: f64, k: usize, n: usize, t_k: f64) -> f64 {
        if !self.active {
            return 1.0;
        }
        let il = 10f64.powf(-self.il_db / 20.0);
        let Some(fwhm_hz) = self.fwhm_hz else {
            return il; // flat loss, no passband
        };
        let spec = AwgSpectrum {
            grid: ChannelGrid {
                f0_hz: if self.lambda0_m > 0.0 {
                    C0 / self.lambda0_m
                } else {
                    0.0
                },
                df_hz: self.df_hz,
                // A mux/demux is a one-shot filter bank, not a cyclic router,
                // so it has no FSR to fold into. Zero means aperiodic; the
                // out-of-band guard then uses the grid span instead.
                fsr_hz: 0.0,
                fwhm_hz,
                shape_p: self.shape_p,
                dlambda_dt_m_per_k: self.dlambda_dt_m_per_k,
                t_nom_k: self.t_nom_k,
            },
            xt_adj: 0.0,
            xt_bg: 0.0,
            il_peak: 1.0,
            il_tilt: 1.0,
        };
        il * spec.power(lambda, k, n, t_k).sqrt()
    }
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
// `expand_bundle_ports` in fairchild-parser.

/// Stamp one channel's route between the bus block and the per-channel block,
/// shared by `fc_mux` (`bus_out = true`) and `fc_demux` (`bus_out = false`).
///
/// The forward field and the λ tag travel from the input side to the output
/// side; the backward pair travels the other way, so the **input** side owns
/// those two wires. Both devices used to stamp all five in the forward
/// direction, which put a mux's driver on the bus's backward wires — the wires
/// the next device down the bus already drives, because a device's `in` port
/// drives its own backward pair. Two drivers on one node leaves the block
/// rank-deficient, and the solve reports nothing: it returns a `gmin`-weighted
/// average of the two answers. Meanwhile the backward light never reached the
/// channel ports at all, so a reflection anywhere past a mux went missing.
/// `newton::check_exclusive_potential_drivers` now refuses that shape by name.
fn stamp_route(
    mat: &mut MnaMatrix,
    branches: &[Option<usize>],
    nodes: &[NodeId],
    n: usize,
    wpc: usize,
    k: usize,
    amp: f64,
    bus_out: bool,
) {
    // The λ wire is the last of the `wpc`, and it is not routed here at all: a
    // wavelength is resolved before the solve, so it has no equation and no
    // branch row. Only the field wires do.
    let bpc = wpc - 1;
    for w in 0..bpc {
        let bus_w = nodes[wpc * k + w];
        let ch_w = nodes[wpc * (n + k) + w];
        // Wires 2 and 3 of a 5-wire bundle are the backward pair. Under
        // unidirectional propagation there is no such pair and every wire runs
        // forward.
        let backward = wpc == 5 && (w == 2 || w == 3);
        let (dst, src) = if bus_out != backward {
            (bus_w, ch_w)
        } else {
            (ch_w, bus_w)
        };
        stamp_potential_eq(mat, branches, bpc * k + w, dst, &[(src, -amp)]);
    }
}

/// Identity-routing combiner: N single-channel optical bundles → 1 N-channel
/// bundle.  Pin 1 (and the first bundle wire block) is the bus output.
/// λ routing for a bundle bridge: slot `k` of the input block becomes slot `k`
/// of the output block.
///
/// Both blocks are `wpc·n` wires wide — a mux's bus carries `n` channels and its
/// `n` single-channel ports carry one each — so the two are the same shape and
/// the mapping is slot for slot. Which block is the input is exactly what
/// `LAMBDA_BASE` already records.
fn bridge_lambda_routing(lambda_base: usize, wpc: usize, n: usize) -> Vec<(usize, usize)> {
    if n == 0 {
        return Vec::new();
    }
    let lam = wpc - 1;
    let in_base = lambda_base * wpc * n;
    let out_base = (1 - lambda_base) * wpc * n;
    (0..n)
        .map(|k| (in_base + wpc * k + lam, out_base + wpc * k + lam))
        .collect()
}

pub struct NativeMux {
    n_channels: usize,
    wpc: usize,
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
    filter: ChannelFilter,
    /// Field transmission per channel, refreshed each eval from `lambda_m`.
    /// All ones until a `ChannelFilter` parameter is set.
    amp: Vec<f64>,
    /// Each channel's resolved λ (m), read from the *input* side's tags.
    /// Zero means "never told", which reads as the filter's own grid centre.
    lambda_m: Vec<f64>,
}

impl Default for NativeMux {
    fn default() -> Self {
        Self::new()
    }
}

impl NativeMux {
    /// Which wire block carries the *input* λ tags, in units of `wpc·n`. A mux
    /// reads from the per-channel block (second), a demux from the bus (first).
    const LAMBDA_BASE: usize = 1;

    pub fn new() -> Self {
        Self {
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
            filter: ChannelFilter::new(),
            amp: Vec::new(),
            lambda_m: Vec::new(),
        }
    }
}

impl Device for NativeMux {
    fn num_terminals(&self) -> usize {
        self.nodes.len()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
        self.filter.lambda0_m = ctx.lambda_center_m;
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 2 * wpc; // bus channel (wpc) + per-channel wires (wpc)
        assert!(
            !terminals.is_empty() && terminals.len().is_multiple_of(stride),
            "fc_mux: terminal count must be a positive multiple of {stride} \
             (wpc={wpc}: bus wires + per-channel wires); got {}",
            terminals.len()
        );
        let n = terminals.len() / stride;
        self.n_channels = n;
        self.nodes = terminals.to_vec();
        self.branches = vec![None; (wpc - 1) * n];
        self.amp = vec![1.0; n];
        self.lambda_m = vec![0.0; n];
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
        self.filter.set(&name.to_lowercase(), value)
    }

    /// A bridge is not slot-for-slot within one bundle, so `OpticalSegment`'s
    /// declaration does not cover it: a mux moves port `k`'s label onto bus slot
    /// `k`, and a demux the other way.
    fn lambda_routing(&self) -> Vec<(usize, usize)> {
        bridge_lambda_routing(Self::LAMBDA_BASE, self.wpc, self.n_channels)
    }

    /// A bridge reads its λ tags from the input side, which `LAMBDA_BASE`
    /// already names: slot `k` of that block.
    fn set_resolved_lambda(&mut self, per_terminal: &[f64]) {
        let (n, wpc) = (self.n_channels, self.wpc);
        let base = Self::LAMBDA_BASE * wpc * n;
        for k in 0..n {
            if let Some(&v) = per_terminal.get(base + wpc * k + wpc - 1) {
                self.lambda_m[k] = v;
            }
        }
    }

    fn eval(&mut self, _x: &[f64], _flags: EvalFlags, ctx: &SimContext) {
        if !self.filter.active {
            return;
        }
        let n = self.n_channels;
        for k in 0..n {
            let lambda = match self.lambda_m.get(k) {
                Some(&v) if v > 0.0 => v,
                _ => self.filter.lambda0_m,
            };
            self.amp[k] = self.filter.amp(lambda, k, n, ctx.temperature);
        }
    }

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        for k in 0..self.n_channels {
            stamp_route(
                mat,
                &self.branches,
                &self.nodes,
                self.n_channels,
                self.wpc,
                k,
                self.amp[k],
                true,
            );
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }
}

/// Identity-routing splitter: 1 N-channel optical bundle → N single-channel
/// bundles.  Pin 1 (and the first bundle wire block) is the bus input.
pub struct NativeDemux {
    n_channels: usize,
    wpc: usize,
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
    filter: ChannelFilter,
    /// Field transmission per channel, refreshed each eval from `lambda_m`.
    /// All ones until a `ChannelFilter` parameter is set.
    amp: Vec<f64>,
    /// Each channel's resolved λ (m), read from the *input* side's tags.
    /// Zero means "never told", which reads as the filter's own grid centre.
    lambda_m: Vec<f64>,
}

impl Default for NativeDemux {
    fn default() -> Self {
        Self::new()
    }
}

impl NativeDemux {
    /// See [`NativeMux::LAMBDA_BASE`]. A demux is fed by the bus block.
    const LAMBDA_BASE: usize = 0;

    pub fn new() -> Self {
        Self {
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
            filter: ChannelFilter::new(),
            amp: Vec::new(),
            lambda_m: Vec::new(),
        }
    }
}

impl Device for NativeDemux {
    fn num_terminals(&self) -> usize {
        self.nodes.len()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
        self.filter.lambda0_m = ctx.lambda_center_m;
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 2 * wpc;
        assert!(
            !terminals.is_empty() && terminals.len().is_multiple_of(stride),
            "fc_demux: terminal count must be a positive multiple of {stride} \
             (wpc={wpc}: bus wires + per-channel wires); got {}",
            terminals.len()
        );
        let n = terminals.len() / stride;
        self.n_channels = n;
        self.nodes = terminals.to_vec();
        self.branches = vec![None; (wpc - 1) * n];
        self.amp = vec![1.0; n];
        self.lambda_m = vec![0.0; n];
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
        self.filter.set(&name.to_lowercase(), value)
    }

    /// A bridge is not slot-for-slot within one bundle, so `OpticalSegment`'s
    /// declaration does not cover it: a mux moves port `k`'s label onto bus slot
    /// `k`, and a demux the other way.
    fn lambda_routing(&self) -> Vec<(usize, usize)> {
        bridge_lambda_routing(Self::LAMBDA_BASE, self.wpc, self.n_channels)
    }

    /// A bridge reads its λ tags from the input side, which `LAMBDA_BASE`
    /// already names: slot `k` of that block.
    fn set_resolved_lambda(&mut self, per_terminal: &[f64]) {
        let (n, wpc) = (self.n_channels, self.wpc);
        let base = Self::LAMBDA_BASE * wpc * n;
        for k in 0..n {
            if let Some(&v) = per_terminal.get(base + wpc * k + wpc - 1) {
                self.lambda_m[k] = v;
            }
        }
    }

    fn eval(&mut self, _x: &[f64], _flags: EvalFlags, ctx: &SimContext) {
        if !self.filter.active {
            return;
        }
        let n = self.n_channels;
        for k in 0..n {
            let lambda = match self.lambda_m.get(k) {
                Some(&v) if v > 0.0 => v,
                _ => self.filter.lambda0_m,
            };
            self.amp[k] = self.filter.amp(lambda, k, n, ctx.temperature);
        }
    }

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        for k in 0..self.n_channels {
            stamp_route(
                mat,
                &self.branches,
                &self.nodes,
                self.n_channels,
                self.wpc,
                k,
                self.amp[k],
                false,
            );
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }
}
