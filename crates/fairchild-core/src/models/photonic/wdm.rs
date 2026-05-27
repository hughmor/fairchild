use super::stamp_potential_eq;
use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

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
    wpc: usize,
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
}

impl Default for NativeMux {
    fn default() -> Self {
        Self::new()
    }
}

impl NativeMux {
    pub fn new() -> Self {
        Self {
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
        }
    }
}

impl Device for NativeMux {
    fn num_terminals(&self) -> usize {
        self.nodes.len()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
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
        self.branches = vec![None; wpc * n];
    }

    fn num_extra_nodes(&self) -> usize {
        self.branches.len()
    }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..self.branches.len() {
            self.branches[i] = Some(first_idx + i);
        }
    }

    fn eval(&mut self, _x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {}

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n = self.n_channels;
        let wpc = self.wpc;
        for k in 0..n {
            for w in 0..wpc {
                let bus_w = self.nodes[wpc * k + w];
                let ch_w = self.nodes[wpc * (n + k) + w];
                // Identity-route every wire (fw, bw, λ) — bus reads from channel.
                stamp_potential_eq(mat, &self.branches, wpc * k + w, bus_w, &[(ch_w, -1.0)]);
            }
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
}

impl Default for NativeDemux {
    fn default() -> Self {
        Self::new()
    }
}

impl NativeDemux {
    pub fn new() -> Self {
        Self {
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
        }
    }
}

impl Device for NativeDemux {
    fn num_terminals(&self) -> usize {
        self.nodes.len()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
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
        self.branches = vec![None; wpc * n];
    }

    fn num_extra_nodes(&self) -> usize {
        self.branches.len()
    }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..self.branches.len() {
            self.branches[i] = Some(first_idx + i);
        }
    }

    fn eval(&mut self, _x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {}

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n = self.n_channels;
        let wpc = self.wpc;
        for k in 0..n {
            for w in 0..wpc {
                let bus_w = self.nodes[wpc * k + w];
                let ch_w = self.nodes[wpc * (n + k) + w];
                // Channels drive FROM bus.
                stamp_potential_eq(mat, &self.branches, wpc * k + w, ch_w, &[(bus_w, -1.0)]);
            }
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }
}
