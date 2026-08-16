use crate::device::{Device, EvalFlags, NodeId, SimContext, Q_ELECTRON};
use crate::mna::MnaMatrix;

// ────────────────────────────────────────────────────────────────────────
// Native photodetector (PIN — instantaneous responsivity)
// ────────────────────────────────────────────────────────────────────────

/// PIN photodetector with linear responsivity and a shunt resistance.
///
/// Physics:  I_ph = R · (re² + im²) + I_dark
/// Flows cathode → anode (reverse-biased junction convention).  A shunt
/// resistance models the junction impedance.
///
/// The current is a nonlinear function of the optical inputs, so this
/// device contributes both a residual (I_ph) and the linearised Jacobian
/// terms ∂I_ph/∂V(in_re), ∂I_ph/∂V(in_im) and 1/R_shunt for the V/R
/// shunt.  No internal nodes are needed — the photocurrent stamps
/// directly between the electrical terminals.
/// Variable-arity to support WDM: ONE physical photodetector with one shared
/// dark current and shunt, summing photocurrents across N optical channels.
/// Terminal layout (3·N + 2 terminals for N channels):
///
///   [in.0.re, in.0.im, in.0.λ,  in.1.re, in.1.im, in.1.λ,  ...,  in.{N-1}.λ,
///    anode, cathode]
///
/// The parser does NOT replicate this device (see fairchild-parser's
/// `expand_bundle_ports` exception list) — instead the device flattens all
/// channels into one terminal block.  Photocurrent
/// `I_ph = responsivity · Σ_k (re_k² + im_k²) + i_dark` is computed in one
/// place; the shunt `1/r_shunt` is stamped once between anode and cathode.
pub struct NativePhotodetector {
    responsivity: f64,
    i_dark: f64,
    r_shunt: f64,
    r_series: f64,
    n_channels: usize,
    wpc: usize, // 3 (unidir) or 5 (bidir)
    /// Smallest terminal count that would have been usable, recorded when
    /// `setup_instance` is handed one it cannot use. `num_terminals()` reads
    /// `nodes.len()`, which a refused setup leaves at 0 — a misleading number
    /// to quote back at the user. Declining rather than panicking is what lets
    /// `build_devices_with_footprints` name the element and both counts.
    min_terminals: Option<usize>,
    nodes: Vec<NodeId>,
    v_int_idx: Option<usize>,
    has_internal: bool,
    i_ph: f64,
    // Linearisation coefficients per (channel, direction).  In unidir mode
    // only g[k].0/.1 are used; in bidir mode the .2/.3 entries cover the bw
    // wires.  A real PIN absorbs every photon — both fw and bw light heat
    // the same junction and produce one summed photocurrent.
    g_re_fw: Vec<f64>,
    g_im_fw: Vec<f64>,
    g_re_bw: Vec<f64>,
    g_im_bw: Vec<f64>,
    v_re_fw_op: Vec<f64>,
    v_im_fw_op: Vec<f64>,
    v_re_bw_op: Vec<f64>,
    v_im_bw_op: Vec<f64>,
    v_j_op: f64,
}

impl Default for NativePhotodetector {
    fn default() -> Self {
        Self::new()
    }
}

impl NativePhotodetector {
    pub fn new() -> Self {
        Self {
            responsivity: 1.0,
            i_dark: 1e-9,
            r_shunt: 1e6,
            r_series: 0.0,
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            min_terminals: None,
            v_int_idx: None,
            has_internal: false,
            i_ph: 0.0,
            g_re_fw: Vec::new(),
            g_im_fw: Vec::new(),
            g_re_bw: Vec::new(),
            g_im_bw: Vec::new(),
            v_re_fw_op: Vec::new(),
            v_im_fw_op: Vec::new(),
            v_re_bw_op: Vec::new(),
            v_im_bw_op: Vec::new(),
            v_j_op: 0.0,
        }
    }
}

impl Device for NativePhotodetector {
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
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        // Layout: wpc·N (bundle inputs) + 2 (anode, cathode).
        if terminals.len() < wpc + 2 || !(terminals.len() - 2).is_multiple_of(wpc) {
            self.min_terminals = Some(wpc + 2);
            return;
        }
        self.min_terminals = None;
        let n = (terminals.len() - 2) / wpc;
        self.n_channels = n;
        self.nodes = terminals.to_vec();
        self.g_re_fw = vec![0.0; n];
        self.g_im_fw = vec![0.0; n];
        self.g_re_bw = vec![0.0; n];
        self.g_im_bw = vec![0.0; n];
        self.v_re_fw_op = vec![0.0; n];
        self.v_im_fw_op = vec![0.0; n];
        self.v_re_bw_op = vec![0.0; n];
        self.v_im_bw_op = vec![0.0; n];
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "responsivity" => {
                self.responsivity = value;
                true
            }
            "i_dark" | "i_dark_a" => {
                self.i_dark = value;
                true
            }
            "r_shunt" => {
                self.r_shunt = value;
                true
            }
            "r_series" | "r_s" => {
                self.r_series = value.max(0.0);
                true
            }
            // `c_par` (and a bias-dependent C_j(V)) lands with the L2 PD
            // model — it needs device-internal reactive state which the
            // current Device trait doesn't expose.  Accept the keyword to
            // keep schematics forward-compatible; document as no-op.
            "c_par" | "c_j0" | "c_par_f" => true,
            _ => false,
        }
    }

    fn num_extra_nodes(&self) -> usize {
        // Allocate one internal MNA row (the "junction node" between r_series
        // and the parallel R_sh || I_ph stack) when r_series is non-zero.
        if self.r_series > 0.0 {
            1
        } else {
            0
        }
    }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        if self.r_series > 0.0 {
            self.v_int_idx = Some(first_idx);
            self.has_internal = true;
        } else {
            self.v_int_idx = None;
            self.has_internal = false;
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let elec_base = wpc * n;
        let v_j_node = match self.v_int_idx {
            Some(i) => x[i],
            None => self.nodes[elec_base].map_or(0.0, |i| x[i]),
        };
        let v_c = self.nodes[elec_base + 1].map_or(0.0, |i| x[i]);
        let mut p_total = 0.0;
        for k in 0..n {
            let v_re_fw = self.nodes[wpc * k].map_or(0.0, |i| x[i]);
            let v_im_fw = self.nodes[wpc * k + 1].map_or(0.0, |i| x[i]);
            p_total += v_re_fw * v_re_fw + v_im_fw * v_im_fw;
            self.g_re_fw[k] = 2.0 * self.responsivity * v_re_fw;
            self.g_im_fw[k] = 2.0 * self.responsivity * v_im_fw;
            self.v_re_fw_op[k] = v_re_fw;
            self.v_im_fw_op[k] = v_im_fw;
            if wpc == 5 {
                let v_re_bw = self.nodes[wpc * k + 2].map_or(0.0, |i| x[i]);
                let v_im_bw = self.nodes[wpc * k + 3].map_or(0.0, |i| x[i]);
                p_total += v_re_bw * v_re_bw + v_im_bw * v_im_bw;
                self.g_re_bw[k] = 2.0 * self.responsivity * v_re_bw;
                self.g_im_bw[k] = 2.0 * self.responsivity * v_im_bw;
                self.v_re_bw_op[k] = v_re_bw;
                self.v_im_bw_op[k] = v_im_bw;
            }
        }
        self.i_ph = self.responsivity * p_total + self.i_dark;
        self.v_j_op = v_j_node - v_c;
    }

    fn load_residual(&self, b: &mut [f64]) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let mut nonlin_remainder = self.i_ph;
        for k in 0..n {
            nonlin_remainder -=
                self.g_re_fw[k] * self.v_re_fw_op[k] + self.g_im_fw[k] * self.v_im_fw_op[k];
            if wpc == 5 {
                nonlin_remainder -=
                    self.g_re_bw[k] * self.v_re_bw_op[k] + self.g_im_bw[k] * self.v_im_bw_op[k];
            }
        }
        // The shunt is *linear*, so it belongs in the Jacobian and nowhere else.
        // Subtracting `v_j_op/r_shunt` here as well cancelled the stamped `g_sh`
        // at the solution: the load saw `v = I_ph·R` with `r_shunt` inert, while
        // `∂f/∂x` still carried the conductance.  Newton did not care — the two
        // errors agree at the fixed point — but the adjoint did, because a
        // gradient through the detector came out wrong by `R_load/r_shunt`.
        let i_eq = -nonlin_remainder;
        let elec_base = wpc * n;
        let v_j_node = self.v_int_idx.or(self.nodes[elec_base]);
        if let Some(j) = v_j_node {
            b[j] -= i_eq;
        }
        if let Some(c) = self.nodes[elec_base + 1] {
            b[c] += i_eq;
        }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let elec_base = wpc * n;
        let a_idx = self.nodes[elec_base];
        let c_idx = self.nodes[elec_base + 1];
        let g_sh = 1.0 / self.r_shunt;
        let j_idx = self.v_int_idx.or(a_idx);
        if let Some(j) = self.v_int_idx {
            let g_s = 1.0 / self.r_series;
            if let Some(a) = a_idx {
                mat.a[a][a] += g_s;
                mat.a[a][j] -= g_s;
                mat.a[j][a] -= g_s;
            }
            mat.a[j][j] += g_s;
        }
        if let Some(j) = j_idx {
            mat.a[j][j] += g_sh;
            if let Some(c) = c_idx {
                mat.a[j][c] -= g_sh;
            }
        }
        if let Some(c) = c_idx {
            mat.a[c][c] += g_sh;
            if let Some(j) = j_idx {
                mat.a[c][j] -= g_sh;
            }
        }
        for k in 0..n {
            let r_fw = self.nodes[wpc * k];
            let i_fw = self.nodes[wpc * k + 1];
            if let Some(j) = j_idx {
                if let Some(r) = r_fw {
                    mat.a[j][r] -= self.g_re_fw[k];
                }
                if let Some(i) = i_fw {
                    mat.a[j][i] -= self.g_im_fw[k];
                }
            }
            if let Some(c) = c_idx {
                if let Some(r) = r_fw {
                    mat.a[c][r] += self.g_re_fw[k];
                }
                if let Some(i) = i_fw {
                    mat.a[c][i] += self.g_im_fw[k];
                }
            }
            if wpc == 5 {
                let r_bw = self.nodes[wpc * k + 2];
                let i_bw = self.nodes[wpc * k + 3];
                if let Some(j) = j_idx {
                    if let Some(r) = r_bw {
                        mat.a[j][r] -= self.g_re_bw[k];
                    }
                    if let Some(i) = i_bw {
                        mat.a[j][i] -= self.g_im_bw[k];
                    }
                }
                if let Some(c) = c_idx {
                    if let Some(r) = r_bw {
                        mat.a[c][r] += self.g_re_bw[k];
                    }
                    if let Some(i) = i_bw {
                        mat.a[c][i] += self.g_im_bw[k];
                    }
                }
            }
        }
    }

    /// Shot noise on the detected current: `S_i = 2q·I` (one-sided), between
    /// the same terminal pair the photocurrent flows through.
    ///
    /// `i_ph` already carries `responsivity·ΣP + i_dark`, and both terms shot,
    /// so one source covers the pair — a dark-current-limited receiver and a
    /// signal-limited one come out of the same expression.  Nothing here is
    /// frequency-dependent; the receiver's own bandwidth shapes it through the
    /// transfer impedance.
    ///
    /// ponytail: no excess-noise factor.  An APD needs `F(M)·M²`, which is a
    /// second parameter and a second model — add it with the APD, not here.
    fn noise_sources(&self, _ctx: &SimContext, _freq: f64) -> Vec<(NodeId, NodeId, f64)> {
        let s_i = 2.0 * Q_ELECTRON * self.i_ph.abs();
        if s_i == 0.0 {
            return Vec::new();
        }
        let elec_base = self.wpc * self.n_channels;
        let j_idx = self.v_int_idx.or(self.nodes[elec_base]);
        vec![(j_idx, self.nodes[elec_base + 1], s_i)]
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }
}
