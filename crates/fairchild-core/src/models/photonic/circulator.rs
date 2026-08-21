use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;
use crate::warn_user;

// ────────────────────────────────────────────────────────────────────────
// Native 3-port circulator (bidir-only)
// ────────────────────────────────────────────────────────────────────────

/// 3-port circulator.  Routes light cyclically: light entering port 1
/// exits port 2; entering port 2 exits port 3; entering port 3 exits
/// port 1.  Requires `enable_bidirectional=1` because two of the three
/// routes are carried on backward wires that a unidirectional bundle does
/// not have.  Errors at setup_instance if bidir is off.
///
/// # Which wire is which
///
/// The convention is the **along-chain** one every other device uses, not a
/// port-relative one: `fw` means "away from port 1, toward port 3", the same
/// direction at all three ports.  So each port plays one of the two ordinary
/// port roles, and the circulator drops into a chain like anything else:
///
/// | Port | Role | Drives | Reads |
/// |---|---|---|---|
/// | 1 | like a waveguide's `in` | `re_bw`, `im_bw` | `re_fw`, `im_fw`, λ |
/// | 2 | like a waveguide's `out` | `re_fw`, `im_fw`, λ | `re_bw`, `im_bw` |
/// | 3 | like a waveguide's `out` | `re_fw`, `im_fw`, λ | `re_bw`, `im_bw` |
///
/// The three routes then read straight off that table:
///
/// ```text
///   port_2.fw = port_1.fw      light in at 1 leaves at 2
///   port_3.fw = port_2.bw      light in at 2 leaves at 3
///   port_1.bw = port_3.bw      light in at 3 leaves at 1
/// ```
///
/// This used to be port-relative — `fw` meant "into me" at all three ports —
/// which made every port behave like an `in` port.  Wiring port 2 or 3 onward
/// into anything then put two drivers on that bundle's backward wires and none
/// on its forward ones: the block went rank-deficient (silently averaged) and
/// the routed light never left the circulator.  The old convention was
/// documented rather than fixed, on the grounds that a user would drive the
/// wires by name; a circulator exists to be put where light comes back, so it
/// has to compose.  `newton::check_exclusive_potential_drivers` refuses the
/// collision now, and this convention no longer causes one.
///
/// λ is tied from port 1 onto ports 2 and 3, matching their `out` role.
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

    /// The two `out`-role ports take port 1's tag, matching the λ ties stamped
    /// below. Bidirectional-only, so the per-channel stride is 3 ports × 5
    /// wires and λ sits last in each port.
    fn lambda_routing(&self) -> Vec<(usize, usize)> {
        let n = self.n_channels;
        (0..n)
            .flat_map(|k| {
                let base = 15 * k;
                [(base + 4, base + 5 + 4), (base + 4, base + 10 + 4)]
            })
            .collect()
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
            // In at port 1, out at port 2 — both on the forward wires, because
            // that is the direction along the chain from 1 to 2.
            stamp_potential_eq(mat, &self.branches, b, p1_re_fw, &[(p0_re_fw, -1.0)]);
            stamp_potential_eq(mat, &self.branches, b + 1, p1_im_fw, &[(p0_im_fw, -1.0)]);
            // In at port 2, out at port 3. Port 2 is an `out`-role port, so
            // light arriving there comes back on its backward wires.
            stamp_potential_eq(mat, &self.branches, b + 2, p2_re_fw, &[(p1_re_bw, -1.0)]);
            stamp_potential_eq(mat, &self.branches, b + 3, p2_im_fw, &[(p1_im_bw, -1.0)]);
            // In at port 3, out at port 1 — the whole return path, arriving on
            // port 3's backward wires and leaving on port 1's.
            stamp_potential_eq(mat, &self.branches, b + 4, p0_re_bw, &[(p2_re_bw, -1.0)]);
            stamp_potential_eq(mat, &self.branches, b + 5, p0_im_bw, &[(p2_im_bw, -1.0)]);
            // λ ties: the two `out`-role ports take port 1's tag.
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
