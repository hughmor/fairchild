use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;
use crate::warn_user;

// ────────────────────────────────────────────────────────────────────────
// Native 3-port circulator (bidir-only)
// ────────────────────────────────────────────────────────────────────────

/// 3-port circulator.  Routes light cyclically: light entering port 1
/// exits port 2; entering port 2 exits port 3; entering port 3 exits
/// port 1.  Requires `enable_bidirectional=1` because the routing
/// fundamentally needs each port to support both an incoming wave (fw,
/// inward to the circulator) and an outgoing wave (bw, outward from the
/// circulator).  Errors at setup_instance if bidir is off.
///
/// Wire convention (consistent across the circulator): at every port,
/// `re_fw`/`im_fw` represent the wave flowing INWARD (toward the device)
/// and `re_bw`/`im_bw` represent the wave flowing OUTWARD.  Internal
/// routing:
///   port_p.bw = port_((p+2) mod 3).fw   — for re and im, every channel
/// (light entering port (p-1) leaves at port p, mod 3).
///
/// λ is tied across all three ports: `port_1.λ = port_0.λ`, `port_2.λ =
/// port_0.λ`.  This works whether the laser drives port 0, 1, or 2 —
/// SPICE branch equations resolve the cycle consistently.
///
/// Bundle-aware: 3·wpc·N terminals for N WDM channels.  Per channel
/// branch count: 6 re/im routing + 2 λ ties = 8.
pub struct NativeCirculator {
    n_channels: usize,
    wpc: usize,
    /// Smallest terminal count that would have been usable, recorded when
    /// `setup_instance` is handed one it cannot use. `num_terminals()` reads
    /// `nodes.len()`, which a refused setup leaves at 0 — a misleading number
    /// to quote back at the user. Declining rather than panicking is what lets
    /// `build_devices_with_footprints` name the element and both counts.
    min_terminals: Option<usize>,
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
}

impl Default for NativeCirculator {
    fn default() -> Self {
        Self::new()
    }
}

impl NativeCirculator {
    pub fn new() -> Self {
        Self {
            n_channels: 0,
            wpc: 5,
            nodes: Vec::new(),
            min_terminals: None,
            branches: Vec::new(),
        }
    }
}

impl Device for NativeCirculator {
    fn num_terminals(&self) -> usize {
        if let Some(min) = self.min_terminals {
            return min;
        }
        self.nodes.len()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        self.wpc = 5;
        let stride = 3 * 5; // 3 ports × 5 wires per channel
        if ctx.wires_per_channel() != 5 {
            // Unidirectional: the netlist supplied 3·3·N wires and this device
            // needs 3·5·N, so the caller's terminal-count error fires with true
            // numbers. It cannot know *why*, so say it here — the count is a
            // symptom of the missing option, not an independent mistake.
            warn_user!(
                "fc_circulator requires bidirectional propagation; set \
                 `.options enable_bidirectional=1` (or via CLI / Python). Without it \
                 each optical port carries {} wires per channel instead of 5, which \
                 is the terminal-count error that follows.",
                ctx.wires_per_channel()
            );
            self.min_terminals = Some(stride);
            return;
        }
        if terminals.is_empty() || !terminals.len().is_multiple_of(stride) {
            self.min_terminals = Some(stride);
            return;
        }
        self.min_terminals = None;
        let n = terminals.len() / stride;
        self.n_channels = n;
        self.nodes = terminals.to_vec();
        self.branches = vec![None; 8 * n];
    }

    fn num_extra_nodes(&self) -> usize {
        self.branches.len()
    }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..self.branches.len() {
            self.branches[i] = Some(first_idx + i);
        }
    }

    fn set_real_param(&mut self, _name: &str, _value: f64) -> bool {
        false
    }

    fn eval(&mut self, _x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {}

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n = self.n_channels;
        let wpc = 5;
        // Per channel: stride = 3 ports × 5 wires = 15.
        for k in 0..n {
            let base = 15 * k;
            // Wires per port: [re_fw, im_fw, re_bw, im_bw, λ].
            let port_wires = |p: usize| -> (NodeId, NodeId, NodeId, NodeId, NodeId) {
                let pb = base + wpc * p;
                (
                    self.nodes[pb],     // re_fw
                    self.nodes[pb + 1], // im_fw
                    self.nodes[pb + 2], // re_bw
                    self.nodes[pb + 3], // im_bw
                    self.nodes[pb + 4], // λ
                )
            };
            let (p0_re_fw, p0_im_fw, p0_re_bw, p0_im_bw, p0_lam) = port_wires(0);
            let (p1_re_fw, p1_im_fw, p1_re_bw, p1_im_bw, p1_lam) = port_wires(1);
            let (p2_re_fw, p2_im_fw, p2_re_bw, p2_im_bw, p2_lam) = port_wires(2);
            let b = 8 * k;
            // port_p.bw = port_((p+2) mod 3).fw
            // port_0.bw = port_2.fw
            stamp_potential_eq(mat, &self.branches, b, p0_re_bw, &[(p2_re_fw, -1.0)]);
            stamp_potential_eq(mat, &self.branches, b + 1, p0_im_bw, &[(p2_im_fw, -1.0)]);
            // port_1.bw = port_0.fw
            stamp_potential_eq(mat, &self.branches, b + 2, p1_re_bw, &[(p0_re_fw, -1.0)]);
            stamp_potential_eq(mat, &self.branches, b + 3, p1_im_bw, &[(p0_im_fw, -1.0)]);
            // port_2.bw = port_1.fw
            stamp_potential_eq(mat, &self.branches, b + 4, p2_re_bw, &[(p1_re_fw, -1.0)]);
            stamp_potential_eq(mat, &self.branches, b + 5, p2_im_bw, &[(p1_im_fw, -1.0)]);
            // λ ties: port_1.λ = port_0.λ, port_2.λ = port_0.λ.
            stamp_potential_eq(mat, &self.branches, b + 6, p1_lam, &[(p0_lam, -1.0)]);
            stamp_potential_eq(mat, &self.branches, b + 7, p2_lam, &[(p0_lam, -1.0)]);
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }
}

use super::stamp_potential_eq;
