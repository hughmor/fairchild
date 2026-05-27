use super::stamp_potential_eq;
use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

// ────────────────────────────────────────────────────────────────────────
// Native testbench MZM (idealised Mach-Zehnder modulator)
// ────────────────────────────────────────────────────────────────────────

/// Testbench MZM — represents an idealised lab-bench Mach-Zehnder modulator
/// (free-space or fibre-pigtailed), NOT a real on-chip MZI.  Use this to
/// model the modulator instrument feeding a chip, or as a reference against
/// which to compare a chip-level MZI built from `fc_splitter` +
/// `fc_pn_th_ps` + `fc_splitter` (the "real" MZI lives in the schematic).
///
/// Pins (4 + 6·N for N optical channels):
///   1  in     — optical bundle input
///   2  out    — optical bundle output
///   3  sig    — electrical drive (sig − gnd is the modulation voltage)
///   4  gnd    — modulation return
///
/// Intensity transmission as a function of V_sig:
///   T(V) = α · [ (1 − 1/E_r) · (1 + cos(π V_sig / V_π)) / 2  +  1/E_r ]
/// with α = 10^(−`alpha_dB`/10) (intensity) and E_r the extinction ratio.
/// T(V) ranges from α (V=0) to α/E_r (V = V_π).  Amplitude transmission
/// `t_amp = √T`.  Wavelength passes through unchanged.
///
/// `f_c` (first-order EO cutoff) is accepted but not yet wired into the
/// signal path — it requires device-internal reactive state which lands
/// with the L2 framework.  At this commit `V_sig` reaches the modulator
/// instantaneously regardless of f_c.
pub struct NativeMzm {
    v_pi: f64,
    alpha_int: f64,
    e_r: f64,
    f_c: f64,
    n_channels: usize,
    wpc: usize,
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
    t_amp_cached: f64,
}

impl Default for NativeMzm {
    fn default() -> Self {
        Self::new()
    }
}

impl NativeMzm {
    pub fn new() -> Self {
        Self {
            v_pi: 3.0,
            alpha_int: 1.0,
            e_r: 1.0e3,
            f_c: 1.0e10,
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
            t_amp_cached: 1.0,
        }
    }
}

impl Device for NativeMzm {
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
            terminals.len() >= stride + 2 && (terminals.len() - 2).is_multiple_of(stride),
            "fc_mzm: terminal count must be {stride}·N + 2 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 2) / stride;
        self.n_channels = n;
        self.nodes = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 };
        self.branches = vec![None; bpc * n];
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
            "v_pi" | "vpi" => {
                self.v_pi = value.max(1e-30);
                true
            }
            "alpha" => {
                self.alpha_int = value.clamp(0.0, 1.0);
                true
            }
            "alpha_db" | "il_db" => {
                self.alpha_int = 10f64.powf(-value / 10.0);
                true
            }
            "e_r" | "er" => {
                self.e_r = value.max(1.0);
                true
            }
            "e_r_db" | "er_db" => {
                self.e_r = 10f64.powf(value / 10.0).max(1.0);
                true
            }
            "f_c" | "fc" => {
                self.f_c = value.max(0.0);
                true
            }
            _ => false,
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let elec_base = 2 * wpc * n;
        let v_sig = self.nodes[elec_base].map_or(0.0, |i| x[i]);
        let v_gnd = self.nodes[elec_base + 1].map_or(0.0, |i| x[i]);
        let v_mod = v_sig - v_gnd;
        let cos_term = (std::f64::consts::PI * v_mod / self.v_pi).cos();
        let inv_er = 1.0 / self.e_r;
        let t_int = self.alpha_int * ((1.0 - inv_er) * 0.5 * (1.0 + cos_term) + inv_er);
        self.t_amp_cached = t_int.max(0.0).sqrt();
    }

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let bpc = if wpc == 5 { 5 } else { 3 };
        let lam = wpc - 1;
        let out_base = wpc * n;
        let t = self.t_amp_cached;
        for k in 0..n {
            let in_re_fw = self.nodes[wpc * k];
            let in_im_fw = self.nodes[wpc * k + 1];
            let in_l = self.nodes[wpc * k + lam];
            let out_re_fw = self.nodes[out_base + wpc * k];
            let out_im_fw = self.nodes[out_base + wpc * k + 1];
            let out_l = self.nodes[out_base + wpc * k + lam];
            stamp_potential_eq(mat, &self.branches, bpc * k, out_re_fw, &[(in_re_fw, -t)]);
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k + 1,
                out_im_fw,
                &[(in_im_fw, -t)],
            );
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k + (bpc - 1),
                out_l,
                &[(in_l, -1.0)],
            );
            if wpc == 5 {
                // MZM is reciprocal: same t_amp applies to backward path.
                let in_re_bw = self.nodes[wpc * k + 2];
                let in_im_bw = self.nodes[wpc * k + 3];
                let out_re_bw = self.nodes[out_base + wpc * k + 2];
                let out_im_bw = self.nodes[out_base + wpc * k + 3];
                stamp_potential_eq(
                    mat,
                    &self.branches,
                    bpc * k + 2,
                    in_re_bw,
                    &[(out_re_bw, -t)],
                );
                stamp_potential_eq(
                    mat,
                    &self.branches,
                    bpc * k + 3,
                    in_im_bw,
                    &[(out_im_bw, -t)],
                );
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
