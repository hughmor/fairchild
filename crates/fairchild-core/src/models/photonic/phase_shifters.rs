use super::{
    dB_per_cm_to_neper_per_m, n_eff_at_lambda, stamp_pn_optical, stamp_pn_ths_jacobian,
    stamp_potential_eq, stamp_resistor,
};
use crate::device::{Device, EvalFlags, NodeId, ReactiveBranchSpec, ReactiveKind, SimContext};
use crate::mna::MnaMatrix;

// ────────────────────────────────────────────────────────────────────────
// Native thermal phase shifter
// ────────────────────────────────────────────────────────────────────────

/// Thermal phase shifter (heater).
///
/// Electrical side: resistive heater with conductance `1/R_heater` between
/// anode and cathode.  Joule power `P = V²/R` is converted to an optical
/// phase shift `φ = π · P / P_pi`, where `P_pi` is the heater power
/// required for a π phase shift.
///
/// Optical side: 3-wire bundle in → 3-wire bundle out, applies `exp(-jφ)`.
/// Wavelength passes through unchanged.
///
/// 8 terminals: [in_re, in_im, in_λ, out_re, out_im, out_λ, anode, cathode]
/// 3 internal branch rows for the three direct-potential outputs.
/// Variable-arity: the parser does NOT replicate this device per channel.
/// One instance handles all N optical channels with one shared heater
/// resistor.  Terminal layout: [in.0.re,...,in.{N-1}.λ, out.0.re,...,
/// out.{N-1}.λ, heat_p, heat_n] = 6N + 2.
pub struct NativeThermalPhaseShifter {
    r_heater: f64,
    p_pi: f64,
    n_channels: usize,
    wpc: usize,
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
    c_cached: f64,
    s_cached: f64,
}

impl NativeThermalPhaseShifter {
    pub fn new() -> Self {
        Self {
            r_heater: 1000.0,
            p_pi: 10e-3,
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
            c_cached: 1.0,
            s_cached: 0.0,
        }
    }
}

impl Device for NativeThermalPhaseShifter {
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
            terminals.len() >= stride + 2 && (terminals.len() - 2) % stride == 0,
            "fc_thermal_ps: terminal count must be {stride}·N + 2 (wpc={wpc}); got {}",
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
            "r_heater" | "r" => {
                self.r_heater = value;
                true
            }
            "p_pi" | "p_pi_w" => {
                self.p_pi = value;
                true
            }
            _ => false,
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let elec_base = 2 * wpc * n;
        let v_a = self.nodes[elec_base].map_or(0.0, |i| x[i]);
        let v_c = self.nodes[elec_base + 1].map_or(0.0, |i| x[i]);
        let v = v_a - v_c;
        let p = v * v / self.r_heater;
        let phi = std::f64::consts::PI * p / self.p_pi;
        self.c_cached = phi.cos();
        self.s_cached = phi.sin();
    }

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let bpc = if wpc == 5 { 5 } else { 3 };
        let lam = wpc - 1;
        let out_base = wpc * n;
        let elec_base = 2 * wpc * n;
        // Electrical: ONE shared heater resistor.
        let g = 1.0 / self.r_heater;
        let p = self.nodes[elec_base];
        let m = self.nodes[elec_base + 1];
        if let Some(a) = p {
            mat.a[a][a] += g;
            if let Some(c) = m {
                mat.a[a][c] -= g;
            }
        }
        if let Some(c) = m {
            mat.a[c][c] += g;
            if let Some(a) = p {
                mat.a[c][a] -= g;
            }
        }
        let c_cos = self.c_cached;
        let s_sin = self.s_cached;
        for k in 0..n {
            let in_re_fw = self.nodes[wpc * k];
            let in_im_fw = self.nodes[wpc * k + 1];
            let in_l = self.nodes[wpc * k + lam];
            let out_re_fw = self.nodes[out_base + wpc * k];
            let out_im_fw = self.nodes[out_base + wpc * k + 1];
            let out_l = self.nodes[out_base + wpc * k + lam];
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k,
                out_re_fw,
                &[(in_re_fw, -c_cos), (in_im_fw, -s_sin)],
            );
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k + 1,
                out_im_fw,
                &[(in_re_fw, s_sin), (in_im_fw, -c_cos)],
            );
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k + (bpc - 1),
                out_l,
                &[(in_l, -1.0)],
            );
            if wpc == 5 {
                let in_re_bw = self.nodes[wpc * k + 2];
                let in_im_bw = self.nodes[wpc * k + 3];
                let out_re_bw = self.nodes[out_base + wpc * k + 2];
                let out_im_bw = self.nodes[out_base + wpc * k + 3];
                stamp_potential_eq(
                    mat,
                    &self.branches,
                    bpc * k + 2,
                    in_re_bw,
                    &[(out_re_bw, -c_cos), (out_im_bw, -s_sin)],
                );
                stamp_potential_eq(
                    mat,
                    &self.branches,
                    bpc * k + 3,
                    in_im_bw,
                    &[(out_re_bw, s_sin), (out_im_bw, -c_cos)],
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

// ────────────────────────────────────────────────────────────────────────
// Native thermal phase shifter with thermal time constant (L2 — path B)
// ────────────────────────────────────────────────────────────────────────

/// Thermal phase shifter with a first-order thermal RC.  Same pin layout
/// as `fc_thermal_ps`.  Adds:
///   - `tau_th` — thermal time constant (s).  The optical phase shift
///     tracks the FILTERED heater power rather than the instantaneous
///     Joule dissipation: `dT/dt = (P − T) / tau_th`, with T in
///     normalised "power-equivalent" units so the steady-state phase
///     equals the L1 model's φ = π · P / P_pi.
///
/// Implementation: T is a *state variable* on the MNA matrix — the
/// device allocates one extra row (above its optical branches) and
/// stamps the BE-discretised state equation directly.  This is the
/// "path B" pattern for nonlinear / nonlinear-coupled state, complementary
/// to the linear-companion path A (used by `fc_pn_ps_cap` for C_j(V)).
/// The previous-timestep T is captured via `commit_timestep` after each
/// successful NR convergence.
pub struct NativeThermalPhaseShifterRc {
    r_heater: f64,
    p_pi: f64,
    tau_th: f64,
    n_channels: usize,
    wpc: usize,
    nodes: Vec<NodeId>,
    /// Optical branch rows (re/im/λ per channel).  Same as L1.
    branches: Vec<Option<usize>>,
    /// State row for T(t) (single MNA index allocated alongside branches).
    t_state_idx: Option<usize>,
    /// Previous-timestep value of T, captured by `commit_timestep`.
    t_old: f64,
    /// Cached operating-point quantities (per NR iteration).
    t_op: f64,
    v_h_op: f64,
    c_cached: f64,
    s_cached: f64,
}

impl NativeThermalPhaseShifterRc {
    pub fn new() -> Self {
        Self {
            r_heater: 1000.0,
            p_pi: 10e-3,
            tau_th: 10e-6, // 10 µs typical waveguide heater
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
            t_state_idx: None,
            t_old: 0.0,
            t_op: 0.0,
            v_h_op: 0.0,
            c_cached: 1.0,
            s_cached: 0.0,
        }
    }
}

impl Device for NativeThermalPhaseShifterRc {
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
            terminals.len() >= stride + 2 && (terminals.len() - 2) % stride == 0,
            "fc_thermal_ps_rc: terminal count must be {stride}·N + 2 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 2) / stride;
        self.n_channels = n;
        self.nodes = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 };
        self.branches = vec![None; bpc * n];
    }

    fn num_extra_nodes(&self) -> usize {
        // Optical branches + 1 state row for T(t).
        self.branches.len() + 1
    }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        let n = self.branches.len();
        for i in 0..n {
            self.branches[i] = Some(first_idx + i);
        }
        self.t_state_idx = Some(first_idx + n);
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "r_heater" | "r" => {
                self.r_heater = value;
                true
            }
            "p_pi" | "p_pi_w" => {
                self.p_pi = value;
                true
            }
            "tau_th" | "tau" => {
                self.tau_th = value.max(1e-30);
                true
            }
            _ => false,
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let elec_base = 2 * wpc * n;
        let v_a = self.nodes[elec_base].map_or(0.0, |i| x[i]);
        let v_c = self.nodes[elec_base + 1].map_or(0.0, |i| x[i]);
        self.v_h_op = v_a - v_c;
        // Read T from the state row.
        self.t_op = self.t_state_idx.map_or(0.0, |i| x[i]);
        // Phase from T (filtered power).
        let phi = std::f64::consts::PI * self.t_op / self.p_pi;
        self.c_cached = phi.cos();
        self.s_cached = phi.sin();
    }

    fn load_residual(&self, b: &mut [f64]) {
        // DC: T = P (steady state).  Linearised: T − P_lin(V_h) = 0
        // where P_lin = 2·V_h_op·V_h/R − V_h_op²/R, so
        //   T − 2·V_h_op·V_h/R + V_h_op²/R = 0
        //   ⇒ rearranged with V_h_op² constant on RHS:  b[t_idx] = −V_h_op²/R.
        if let Some(t_idx) = self.t_state_idx {
            let p_op = self.v_h_op * self.v_h_op / self.r_heater;
            b[t_idx] -= p_op;
        }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let bpc = if wpc == 5 { 5 } else { 3 };
        let lam = wpc - 1;
        let out_base = wpc * n;
        let elec_base = 2 * wpc * n;
        // Electrical: shared heater conductance.
        let g = 1.0 / self.r_heater;
        let p = self.nodes[elec_base];
        let m = self.nodes[elec_base + 1];
        if let Some(a) = p {
            mat.a[a][a] += g;
            if let Some(c) = m {
                mat.a[a][c] -= g;
            }
        }
        if let Some(c) = m {
            mat.a[c][c] += g;
            if let Some(a) = p {
                mat.a[c][a] -= g;
            }
        }
        // State row for T (DC: T = P, i.e., T - P_linearised = 0).
        // Stamp: row = +1·T - (2·V_h_op/R)·V_hp + (2·V_h_op/R)·V_hn = +V_h_op²/R
        if let Some(t_idx) = self.t_state_idx {
            mat.a[t_idx][t_idx] += 1.0;
            let two_vop_over_r = 2.0 * self.v_h_op / self.r_heater;
            if let Some(hp) = p {
                mat.a[t_idx][hp] -= two_vop_over_r;
            }
            if let Some(hn) = m {
                mat.a[t_idx][hn] += two_vop_over_r;
            }
        }
        // Optical branches: identical structure to fc_thermal_ps but using
        // c_cached/s_cached derived from T (state).
        let c_cos = self.c_cached;
        let s_sin = self.s_cached;
        for k in 0..n {
            let in_re_fw = self.nodes[wpc * k];
            let in_im_fw = self.nodes[wpc * k + 1];
            let in_l = self.nodes[wpc * k + lam];
            let out_re_fw = self.nodes[out_base + wpc * k];
            let out_im_fw = self.nodes[out_base + wpc * k + 1];
            let out_l = self.nodes[out_base + wpc * k + lam];
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k,
                out_re_fw,
                &[(in_re_fw, -c_cos), (in_im_fw, -s_sin)],
            );
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k + 1,
                out_im_fw,
                &[(in_re_fw, s_sin), (in_im_fw, -c_cos)],
            );
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k + (bpc - 1),
                out_l,
                &[(in_l, -1.0)],
            );
            if wpc == 5 {
                let in_re_bw = self.nodes[wpc * k + 2];
                let in_im_bw = self.nodes[wpc * k + 3];
                let out_re_bw = self.nodes[out_base + wpc * k + 2];
                let out_im_bw = self.nodes[out_base + wpc * k + 3];
                stamp_potential_eq(
                    mat,
                    &self.branches,
                    bpc * k + 2,
                    in_re_bw,
                    &[(out_re_bw, -c_cos), (out_im_bw, -s_sin)],
                );
                stamp_potential_eq(
                    mat,
                    &self.branches,
                    bpc * k + 3,
                    in_im_bw,
                    &[(out_re_bw, s_sin), (out_im_bw, -c_cos)],
                );
            }
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], alpha: f64) {
        // BE-discretised state-row RHS:  b[t_idx] = T_old·α − V_h_op²/(R·τ).
        // (The −V_h_op² term is the linearisation remainder of P_lin =
        //  2·V_h_op·V_h/R − V_h_op²/R; the 2·V_h_op·V_h part is in the
        //  Jacobian.)  Optical branch rows are homogeneous (no residual).
        if let Some(t_idx) = self.t_state_idx {
            let inv_tau = 1.0 / self.tau_th;
            let p_op = self.v_h_op * self.v_h_op / self.r_heater;
            b[t_idx] += self.t_old * alpha - p_op * inv_tau;
        }
    }

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, alpha: f64) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let bpc = if wpc == 5 { 5 } else { 3 };
        let lam = wpc - 1;
        let out_base = wpc * n;
        let elec_base = 2 * wpc * n;
        // Electrical: shared heater conductance.
        let g = 1.0 / self.r_heater;
        let p = self.nodes[elec_base];
        let m = self.nodes[elec_base + 1];
        if let Some(a) = p {
            mat.a[a][a] += g;
            if let Some(c) = m {
                mat.a[a][c] -= g;
            }
        }
        if let Some(c) = m {
            mat.a[c][c] += g;
            if let Some(a) = p {
                mat.a[c][a] -= g;
            }
        }
        // BE state-row Jacobian: T_new·(α + 1/τ) − 2·V_h_op·V_h/(R·τ) = …
        if let Some(t_idx) = self.t_state_idx {
            let inv_tau = 1.0 / self.tau_th;
            mat.a[t_idx][t_idx] += alpha + inv_tau;
            let two_vop_over_r = 2.0 * self.v_h_op / self.r_heater;
            if let Some(hp) = p {
                mat.a[t_idx][hp] -= two_vop_over_r * inv_tau;
            }
            if let Some(hn) = m {
                mat.a[t_idx][hn] += two_vop_over_r * inv_tau;
            }
        }
        // Optical branches (same as DC, c/s from cached T).
        let c_cos = self.c_cached;
        let s_sin = self.s_cached;
        for k in 0..n {
            let in_re_fw = self.nodes[wpc * k];
            let in_im_fw = self.nodes[wpc * k + 1];
            let in_l = self.nodes[wpc * k + lam];
            let out_re_fw = self.nodes[out_base + wpc * k];
            let out_im_fw = self.nodes[out_base + wpc * k + 1];
            let out_l = self.nodes[out_base + wpc * k + lam];
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k,
                out_re_fw,
                &[(in_re_fw, -c_cos), (in_im_fw, -s_sin)],
            );
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k + 1,
                out_im_fw,
                &[(in_re_fw, s_sin), (in_im_fw, -c_cos)],
            );
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k + (bpc - 1),
                out_l,
                &[(in_l, -1.0)],
            );
            if wpc == 5 {
                let in_re_bw = self.nodes[wpc * k + 2];
                let in_im_bw = self.nodes[wpc * k + 3];
                let out_re_bw = self.nodes[out_base + wpc * k + 2];
                let out_im_bw = self.nodes[out_base + wpc * k + 3];
                stamp_potential_eq(
                    mat,
                    &self.branches,
                    bpc * k + 2,
                    in_re_bw,
                    &[(out_re_bw, -c_cos), (out_im_bw, -s_sin)],
                );
                stamp_potential_eq(
                    mat,
                    &self.branches,
                    bpc * k + 3,
                    in_im_bw,
                    &[(out_re_bw, s_sin), (out_im_bw, -c_cos)],
                );
            }
        }
    }

    fn commit_timestep(&mut self, x: &[f64]) {
        // Snapshot T for use as T_old on the next timestep's BE stamp.
        if let Some(idx) = self.t_state_idx {
            self.t_old = x[idx];
        }
    }
}

// ────────────────────────────────────────────────────────────────────────
// Native PN-junction phase shifter
// ────────────────────────────────────────────────────────────────────────

/// PN-junction phase shifter (carrier-depletion or carrier-injection).
///
/// Electrical side: linearised PN junction.  Modelled as a parallel
/// combination of a small ohmic resistance (1 / G_pn) and a linear
/// capacitance (used in transient via the standard Norton C stamp; in DC
/// this contributes nothing).  Voltage dependence of the cap is ignored at
/// this first-pass level.
///
/// Optical side: phase shift `φ = 2π · L · Δn_eff / λ`, where the
/// effective-index change is linearised as `Δn_eff = δn_dV · V_pn`.  This
/// reproduces the small-signal behaviour of either depletion or carrier
/// injection modulators when calibrated to measurements (parameter
/// `dn_dv`).  Wavelength passes through unchanged.
/// Variable-arity to support WDM: the parser does NOT replicate this device
/// per channel (see `BUNDLE_AWARE_MODELS` in fairchild-parser).  One instance
/// handles all N optical channels with one shared PN junction.  Terminal
/// layout for N channels (6·N + 2 total terminals):
///
///   [in.0.re, in.0.im, in.0.λ,  in.1.re, in.1.im, in.1.λ,  ...,  in.{N-1}.λ,
///    out.0.re, out.0.im, out.0.λ,  ...,  out.{N-1}.λ,
///    anode, cathode]
///
/// The wavelength wires read independently per channel (a WDM laser can
/// drive each channel at a different λ), but the electrical conductance,
/// the EO Δn_eff, and the loss factor are all shared — the single physical
/// device sees one V_pn across one junction regardless of how many
/// wavelengths pass through.
pub struct NativePnPhaseShifter {
    length_m: f64,
    n_eff: f64,
    n_g: f64,
    wl_ref_m: f64,
    dn_dv: f64,
    g_pn: f64,
    alpha_neper_m: f64,
    /// When true (default), subtract the absolute propagation phase at
    /// `wl_ref_m` so the device is "transparent" at λ = λ_ref.  Convenient
    /// for testbench rings where the user wants the ring on-resonance at the
    /// laser wavelength by construction.  Set to false for multi-ring designs
    /// (rings of different L) where you want each ring's natural absolute
    /// resonance position — otherwise all rings cluster at λ_ref regardless
    /// of length.  Set via `pin_at_ref=0|1` SPICE parameter.
    pin_at_ref: bool,
    n_channels: usize,
    wpc: usize, // 3 or 5
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
    c_cached: Vec<f64>,
    s_cached: Vec<f64>,
}

impl NativePnPhaseShifter {
    pub fn new() -> Self {
        // Defaults: SOI rib waveguide PN modulator section, R = 8 µm bent.
        //  n_eff / n_g from `scripts/waveguide_simulations/cband_sweep.csv`
        //   (rib_R8 column at 1550 nm).
        //  alpha = 20 dB/cm is dominated by free-carrier absorption from the
        //   typical 5e17 cm⁻³ slab doping; replace with whatever the
        //   `pn_modulator/` sims report for your specific doping profile.
        //  V_pi_L default → dn_dv = wl_ref/(2·V_pi_L).  Set V_pi_L = 0.015
        //   (V·m) so V_pi = 0.015 / L_um·1e-6.  At the typical 1-mm PN
        //   length, V_pi ≈ 15 V (reverse-bias depletion-mode).
        //  pin_at_ref = false: use physical absolute propagation phase so
        //   ring resonances depend on L.  Pin to a ref-wavelength using
        //   `pin_at_ref=1` for testbench rings designed on-resonance.
        Self {
            length_m: 1e-3,
            n_eff: 2.7654,
            n_g: 4.02,
            wl_ref_m: 1.55e-6,
            dn_dv: 1.55e-6 / (2.0 * 0.015),
            g_pn: 1e-3,
            alpha_neper_m: dB_per_cm_to_neper_per_m(20.0),
            pin_at_ref: false,
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
            c_cached: Vec::new(),
            s_cached: Vec::new(),
        }
    }
}

impl Device for NativePnPhaseShifter {
    fn num_terminals(&self) -> usize {
        self.nodes.len()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wl_ref_m = ctx.lambda_center_m;
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 2 * wpc; // in + out bundle
                              // Layout: wpc·N (in) + wpc·N (out) + 2 (anode, cathode).
        assert!(
            terminals.len() >= stride + 2 && (terminals.len() - 2) % stride == 0,
            "fc_pn_ps: terminal count must be {stride}·N + 2 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 2) / stride;
        self.n_channels = n;
        self.nodes = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 }; // branches per channel
        self.branches = vec![None; bpc * n];
        self.c_cached = vec![1.0; n];
        self.s_cached = vec![0.0; n];
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
            "l_um" => {
                self.length_m = value * 1e-6;
                true
            }
            "l_m" | "length" => {
                self.length_m = value;
                true
            }
            "n_g" => {
                self.n_g = value;
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
            "dn_dv" => {
                self.dn_dv = value;
                true
            }
            "g_pn" => {
                self.g_pn = value;
                true
            }
            "v_pi_l" => {
                // Vπ·L (V·m): solve for dn_dv such that the EO phase shift
                // is π at V = Vπ.  2π·L·dn_dv·Vπ/λ_ref = π →
                // dn_dv = λ_ref / (2·L·Vπ).
                if value > 0.0 {
                    self.dn_dv = self.wl_ref_m / (2.0 * value);
                }
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

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let elec_base = 2 * wpc * n;
        let v_a = self.nodes[elec_base].map_or(0.0, |i| x[i]);
        let v_c = self.nodes[elec_base + 1].map_or(0.0, |i| x[i]);
        let v_pn = v_a - v_c;
        let two_pi = 2.0 * std::f64::consts::PI;
        let t_amp = (-self.alpha_neper_m * self.length_m / 2.0).exp();
        let lam = wpc - 1;
        // Reference absolute propagation phase at λ_ref.  When `pin_at_ref`
        // is on we subtract it so the device is "transparent" at λ = λ_ref;
        // otherwise we use the full absolute phase so each ring's natural
        // resonance position depends on L (correct physics for multi-ring
        // designs).
        let phi_ref = if self.pin_at_ref {
            two_pi * self.n_eff * self.length_m / self.wl_ref_m
        } else {
            0.0
        };
        for k in 0..n {
            let lambda = match self.nodes[wpc * k + lam] {
                Some(i) => {
                    let v = x[i];
                    if v.abs() > 1e-9 {
                        v
                    } else {
                        self.wl_ref_m
                    }
                }
                None => self.wl_ref_m,
            };
            let n_eff_lam = n_eff_at_lambda(self.n_eff, self.n_g, self.wl_ref_m, lambda);
            let phi_abs = two_pi * n_eff_lam * self.length_m / lambda;
            let phi_prop = phi_abs - phi_ref;
            let phi_eo = two_pi * self.length_m * self.dn_dv * v_pn / lambda;
            let phi = phi_prop + phi_eo;
            self.c_cached[k] = t_amp * phi.cos();
            self.s_cached[k] = t_amp * phi.sin();
        }
    }

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let bpc = if wpc == 5 { 5 } else { 3 };
        let lam = wpc - 1;
        let out_base = wpc * n;
        let elec_base = 2 * wpc * n;
        // Electrical: ONE shared PN-junction conductance.
        let g = self.g_pn;
        let anode = self.nodes[elec_base];
        let cath = self.nodes[elec_base + 1];
        if let Some(a) = anode {
            mat.a[a][a] += g;
            if let Some(c) = cath {
                mat.a[a][c] -= g;
            }
        }
        if let Some(c) = cath {
            mat.a[c][c] += g;
            if let Some(a) = anode {
                mat.a[c][a] -= g;
            }
        }
        // Optical: per-channel branch equations.
        for k in 0..n {
            let in_re_fw = self.nodes[wpc * k];
            let in_im_fw = self.nodes[wpc * k + 1];
            let in_l = self.nodes[wpc * k + lam];
            let out_re_fw = self.nodes[out_base + wpc * k];
            let out_im_fw = self.nodes[out_base + wpc * k + 1];
            let out_l = self.nodes[out_base + wpc * k + lam];
            let c = self.c_cached[k];
            let s = self.s_cached[k];
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k,
                out_re_fw,
                &[(in_re_fw, -c), (in_im_fw, -s)],
            );
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k + 1,
                out_im_fw,
                &[(in_re_fw, s), (in_im_fw, -c)],
            );
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k + (bpc - 1),
                out_l,
                &[(in_l, -1.0)],
            );
            if wpc == 5 {
                // Backward path mirrors fw: bw entering at out exits at in
                // with same propagation + EO phase shift (one physical
                // junction; same Δn applies to either direction).
                let in_re_bw = self.nodes[wpc * k + 2];
                let in_im_bw = self.nodes[wpc * k + 3];
                let out_re_bw = self.nodes[out_base + wpc * k + 2];
                let out_im_bw = self.nodes[out_base + wpc * k + 3];
                stamp_potential_eq(
                    mat,
                    &self.branches,
                    bpc * k + 2,
                    in_re_bw,
                    &[(out_re_bw, -c), (out_im_bw, -s)],
                );
                stamp_potential_eq(
                    mat,
                    &self.branches,
                    bpc * k + 3,
                    in_im_bw,
                    &[(out_re_bw, s), (out_im_bw, -c)],
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

// ────────────────────────────────────────────────────────────────────────
// Native PN-junction phase shifter — depletion (L2 with bias-dependent C_j)
// ────────────────────────────────────────────────────────────────────────

/// Depletion-mode PN-junction phase shifter with bias-dependent junction
/// capacitance.  Same optical/electrical port layout and stamping as
/// `fc_pn_ps`, plus:
///   - `C_j(V) = C_j0 / (1 − V_pn/V_bi)^m_j` for V_pn ≤ V_bi/2 (depletion),
///     linearly continued above the singularity to keep NR stable when
///     the user wanders into forward bias.  The integrator owns the
///     companion-model state for this capacitance via
///     `Device::reactive_branches`.
///   - `da/dV` — linear loss-vs-bias coefficient (free-carrier absorption
///     in moderate forward bias).
///
/// Tier convention: this is a separate device class, not a `level=` switch
/// on `fc_pn_ps`.  Forward-injection physics (higher dn/dV, large da/dV,
/// carrier-injection time constants) belongs in a future `fc_pn_ps_inj`.
pub struct NativePnPhaseShifterCap {
    length_m: f64,
    n_eff: f64,
    n_g: f64,
    wl_ref_m: f64,
    /// See `NativePnPhaseShifter::pin_at_ref`.  Default true.
    pin_at_ref: bool,
    dn_dv: f64,
    g_pn: f64,
    alpha_neper_m: f64,
    // Bias-dependent C_j parameters.
    c_j0: f64, // F at V_pn = 0
    v_bi: f64, // V — built-in voltage (knee)
    m_j: f64,  // grading coefficient
    // Linear da/dV loss-vs-bias (Np/m per V).  Default 0.
    da_dv: f64,
    n_channels: usize,
    wpc: usize,
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
    c_cached: Vec<f64>,
    s_cached: Vec<f64>,
    // Cached per-NR-iteration values used to re-feed the integrator's
    // companion model AND the per-channel optical loss factor.
    c_j_cached: f64,
    alpha_eff_neper_m: f64,
}

impl NativePnPhaseShifterCap {
    pub fn new() -> Self {
        // Same baseline defaults as `fc_pn_ps` (bent rib SOI PN modulator).
        // Adds: depletion-mode C_j(V) (`c_j0`, `v_bi`, `m_j`) and a
        // linear reverse-bias loss-vs-bias coefficient (`da_dv`).
        Self {
            length_m: 1e-3,
            n_eff: 2.7654,
            n_g: 4.02,
            wl_ref_m: 1.55e-6,
            pin_at_ref: false,
            dn_dv: 1.55e-6 / (2.0 * 0.015),
            g_pn: 1e-3,
            alpha_neper_m: dB_per_cm_to_neper_per_m(20.0),
            c_j0: 20e-15,
            v_bi: 0.7,
            m_j: 0.5,
            da_dv: 0.0,
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
            c_cached: Vec::new(),
            s_cached: Vec::new(),
            c_j_cached: 20e-15,
            alpha_eff_neper_m: 0.0,
        }
    }
}

impl Device for NativePnPhaseShifterCap {
    fn num_terminals(&self) -> usize {
        self.nodes.len()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wl_ref_m = ctx.lambda_center_m;
        self.wpc = ctx.wires_per_channel();
        self.alpha_eff_neper_m = self.alpha_neper_m;
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 2 * wpc;
        assert!(
            terminals.len() >= stride + 2 && (terminals.len() - 2) % stride == 0,
            "fc_pn_ps_cap: terminal count must be {stride}·N + 2 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 2) / stride;
        self.n_channels = n;
        self.nodes = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 };
        self.branches = vec![None; bpc * n];
        self.c_cached = vec![1.0; n];
        self.s_cached = vec![0.0; n];
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
            "l_um" => {
                self.length_m = value * 1e-6;
                true
            }
            "l_m" | "length" => {
                self.length_m = value;
                true
            }
            "n_g" => {
                self.n_g = value;
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
            "dn_dv" => {
                self.dn_dv = value;
                true
            }
            "g_pn" => {
                self.g_pn = value;
                true
            }
            "v_pi_l" => {
                if value > 0.0 {
                    self.dn_dv = self.wl_ref_m / (2.0 * value);
                }
                true
            }
            "alpha_db_cm" => {
                self.alpha_neper_m = dB_per_cm_to_neper_per_m(value);
                true
            }
            "c_j0" => {
                self.c_j0 = value.max(0.0);
                true
            }
            "v_bi" => {
                self.v_bi = value.max(1e-3);
                true
            }
            "m_j" => {
                self.m_j = value.clamp(0.0, 0.99);
                true
            }
            "da_dv" => {
                self.da_dv = value;
                true
            }
            "pin_at_ref" => {
                self.pin_at_ref = value != 0.0;
                true
            }
            _ => false,
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let elec_base = 2 * wpc * n;
        let v_a = self.nodes[elec_base].map_or(0.0, |i| x[i]);
        let v_c = self.nodes[elec_base + 1].map_or(0.0, |i| x[i]);
        let v_pn = v_a - v_c;

        // C_j(V_pn) with linear continuation above the depletion singularity.
        // Knee chosen at V_bi/2 — matches SPICE diode convention.
        let v_knee = 0.5 * self.v_bi;
        self.c_j_cached = if v_pn < v_knee {
            self.c_j0 / (1.0 - v_pn / self.v_bi).powf(self.m_j)
        } else {
            // Linear extrapolation: c(V_knee) + (dc/dV at knee) · (V_pn − V_knee).
            let c_knee = self.c_j0 / (1.0 - v_knee / self.v_bi).powf(self.m_j);
            let dc_dv = c_knee * self.m_j / (self.v_bi - v_knee);
            c_knee + dc_dv * (v_pn - v_knee)
        };
        // Bias-dependent loss: α(V) = α_0 + (da/dV) · max(0, −V_pn) — only
        // adds extra absorption in reverse bias (free-carrier-like).
        let v_rev = (-v_pn).max(0.0);
        self.alpha_eff_neper_m = self.alpha_neper_m + self.da_dv * v_rev;
        let t_amp = (-self.alpha_eff_neper_m * self.length_m / 2.0).exp();

        let two_pi = 2.0 * std::f64::consts::PI;
        let lam = wpc - 1;
        let phi_ref = if self.pin_at_ref {
            two_pi * self.n_eff * self.length_m / self.wl_ref_m
        } else {
            0.0
        };
        for k in 0..n {
            let lambda = match self.nodes[wpc * k + lam] {
                Some(i) => {
                    let v = x[i];
                    if v.abs() > 1e-9 {
                        v
                    } else {
                        self.wl_ref_m
                    }
                }
                None => self.wl_ref_m,
            };
            let n_eff_lam = n_eff_at_lambda(self.n_eff, self.n_g, self.wl_ref_m, lambda);
            let phi_abs = two_pi * n_eff_lam * self.length_m / lambda;
            let phi_prop = phi_abs - phi_ref;
            let phi_eo = two_pi * self.length_m * self.dn_dv * v_pn / lambda;
            let phi = phi_prop + phi_eo;
            self.c_cached[k] = t_amp * phi.cos();
            self.s_cached[k] = t_amp * phi.sin();
        }
    }

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let bpc = if wpc == 5 { 5 } else { 3 };
        let lam = wpc - 1;
        let out_base = wpc * n;
        let elec_base = 2 * wpc * n;
        let g = self.g_pn;
        let anode = self.nodes[elec_base];
        let cath = self.nodes[elec_base + 1];
        if let Some(a) = anode {
            mat.a[a][a] += g;
            if let Some(c) = cath {
                mat.a[a][c] -= g;
            }
        }
        if let Some(c) = cath {
            mat.a[c][c] += g;
            if let Some(a) = anode {
                mat.a[c][a] -= g;
            }
        }
        for k in 0..n {
            let in_re_fw = self.nodes[wpc * k];
            let in_im_fw = self.nodes[wpc * k + 1];
            let in_l = self.nodes[wpc * k + lam];
            let out_re_fw = self.nodes[out_base + wpc * k];
            let out_im_fw = self.nodes[out_base + wpc * k + 1];
            let out_l = self.nodes[out_base + wpc * k + lam];
            let c = self.c_cached[k];
            let s = self.s_cached[k];
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k,
                out_re_fw,
                &[(in_re_fw, -c), (in_im_fw, -s)],
            );
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k + 1,
                out_im_fw,
                &[(in_re_fw, s), (in_im_fw, -c)],
            );
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k + (bpc - 1),
                out_l,
                &[(in_l, -1.0)],
            );
            if wpc == 5 {
                let in_re_bw = self.nodes[wpc * k + 2];
                let in_im_bw = self.nodes[wpc * k + 3];
                let out_re_bw = self.nodes[out_base + wpc * k + 2];
                let out_im_bw = self.nodes[out_base + wpc * k + 3];
                stamp_potential_eq(
                    mat,
                    &self.branches,
                    bpc * k + 2,
                    in_re_bw,
                    &[(out_re_bw, -c), (out_im_bw, -s)],
                );
                stamp_potential_eq(
                    mat,
                    &self.branches,
                    bpc * k + 3,
                    in_im_bw,
                    &[(out_re_bw, s), (out_im_bw, -c)],
                );
            }
        }
    }

    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> {
        // ONE shared depletion capacitance between anode and cathode (the
        // single physical junction).  Bias-dependent value re-queried per
        // NR iteration; the integrator owns the companion-model state.
        let elec_base = 2 * self.wpc * self.n_channels;
        let anode = self.nodes.get(elec_base).copied().flatten();
        let cath = self.nodes.get(elec_base + 1).copied().flatten();
        vec![ReactiveBranchSpec {
            kind: ReactiveKind::Capacitor,
            pos: anode,
            neg: cath,
            value: self.c_j_cached,
        }]
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }
}

// ────────────────────────────────────────────────────────────────────────
// Native combined PN + thermal phase shifter
// ────────────────────────────────────────────────────────────────────────

/// Waveguide segment with both a PN junction (electro-optic) AND a thermal
/// heater driving it.  Δn contributions sum.  Two electrical interfaces
/// (anode/cathode for the PN, heat_p/heat_n for the heater) are independent
/// — driving either alone produces only that physics's phase shift; driving
/// both produces the sum.
///
/// Variable-arity bundle-aware.  Terminal layout for N channels
/// (6·N + 4 total):
///   [in.0.re..in.{N-1}.λ,  out.0.re..out.{N-1}.λ,
///    anode, cathode, heat_p, heat_n]
///
/// At L1 the physics is the linear sum of fc_pn_ps (small-signal Δn_eff =
/// dn/dV · V_pn) and fc_thermal_ps (instantaneous Joule heating →
/// φ_th = π · P / P_pi_th).  L2 will add bias-dependent C_j, tau_th, and
/// distinct reverse/forward EO coefficients; L3 will add carrier dynamics
/// and self-heating from optical absorption.
pub struct NativePnThermalPhaseShifter {
    length_m: f64,
    n_eff: f64,
    n_g: f64,
    wl_ref_m: f64,
    /// See `NativePnPhaseShifter::pin_at_ref`.  Default true.
    pin_at_ref: bool,
    // PN-side params (small-signal electro-optic).
    dn_dv: f64,
    g_pn: f64,
    // Heater-side params (Joule → phase).
    r_heater: f64,
    p_pi_th: f64,
    // Shared loss along the segment.
    alpha_neper_m: f64,
    n_channels: usize,
    wpc: usize,
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
    c_cached: Vec<f64>,
    s_cached: Vec<f64>,
}

impl NativePnThermalPhaseShifter {
    pub fn new() -> Self {
        // Same baseline as `fc_pn_ps`; adds a linear thermal heater pair.
        Self {
            length_m: 1e-3,
            n_eff: 2.7654,
            n_g: 4.02,
            wl_ref_m: 1.55e-6,
            pin_at_ref: false,
            dn_dv: 1.55e-6 / (2.0 * 0.015),
            g_pn: 1e-3,
            r_heater: 1000.0,
            p_pi_th: 10e-3,
            alpha_neper_m: dB_per_cm_to_neper_per_m(20.0),
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
            c_cached: Vec::new(),
            s_cached: Vec::new(),
        }
    }
}

impl Device for NativePnThermalPhaseShifter {
    fn num_terminals(&self) -> usize {
        self.nodes.len()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wl_ref_m = ctx.lambda_center_m;
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 2 * wpc;
        assert!(
            terminals.len() >= stride + 4 && (terminals.len() - 4) % stride == 0,
            "fc_pn_th_ps: terminal count must be {stride}·N + 4 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 4) / stride;
        self.n_channels = n;
        self.nodes = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 };
        self.branches = vec![None; bpc * n];
        self.c_cached = vec![1.0; n];
        self.s_cached = vec![0.0; n];
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
            "l_um" => {
                self.length_m = value * 1e-6;
                true
            }
            "l_m" | "length" => {
                self.length_m = value;
                true
            }
            "n_g" => {
                self.n_g = value;
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
            "dn_dv" => {
                self.dn_dv = value;
                true
            }
            "g_pn" => {
                self.g_pn = value;
                true
            }
            "v_pi_l" => {
                if value > 0.0 {
                    self.dn_dv = self.wl_ref_m / (2.0 * value);
                }
                true
            }
            "r_heater" | "r" => {
                self.r_heater = value;
                true
            }
            "p_pi" | "p_pi_w" | "p_pi_th" => {
                self.p_pi_th = value;
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

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let elec = 2 * wpc * n;
        let v_a = self.nodes[elec].map_or(0.0, |i| x[i]);
        let v_c = self.nodes[elec + 1].map_or(0.0, |i| x[i]);
        let v_hp = self.nodes[elec + 2].map_or(0.0, |i| x[i]);
        let v_hn = self.nodes[elec + 3].map_or(0.0, |i| x[i]);
        let v_pn = v_a - v_c;
        let v_h = v_hp - v_hn;
        let p_heat = v_h * v_h / self.r_heater;
        let phi_th = std::f64::consts::PI * p_heat / self.p_pi_th;
        let two_pi = 2.0 * std::f64::consts::PI;
        let t_amp = (-self.alpha_neper_m * self.length_m / 2.0).exp();
        let lam = wpc - 1;
        let phi_ref = if self.pin_at_ref {
            two_pi * self.n_eff * self.length_m / self.wl_ref_m
        } else {
            0.0
        };
        for k in 0..n {
            let lambda = match self.nodes[wpc * k + lam] {
                Some(i) => {
                    let v = x[i];
                    if v.abs() > 1e-9 {
                        v
                    } else {
                        self.wl_ref_m
                    }
                }
                None => self.wl_ref_m,
            };
            let n_eff_lam = n_eff_at_lambda(self.n_eff, self.n_g, self.wl_ref_m, lambda);
            let phi_abs = two_pi * n_eff_lam * self.length_m / lambda;
            let phi_prop = phi_abs - phi_ref;
            let phi_eo = two_pi * self.length_m * self.dn_dv * v_pn / lambda;
            let phi = phi_prop + phi_eo + phi_th;
            self.c_cached[k] = t_amp * phi.cos();
            self.s_cached[k] = t_amp * phi.sin();
        }
    }

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let bpc = if wpc == 5 { 5 } else { 3 };
        let lam = wpc - 1;
        let out_base = wpc * n;
        let elec = 2 * wpc * n;
        // Electrical (PN side): ONE g_pn between anode and cathode.
        let anode = self.nodes[elec];
        let cath = self.nodes[elec + 1];
        let g_pn = self.g_pn;
        if let Some(a) = anode {
            mat.a[a][a] += g_pn;
            if let Some(c) = cath {
                mat.a[a][c] -= g_pn;
            }
        }
        if let Some(c) = cath {
            mat.a[c][c] += g_pn;
            if let Some(a) = anode {
                mat.a[c][a] -= g_pn;
            }
        }
        // Electrical (heater side): ONE g_heater between heat_p and heat_n.
        let hp = self.nodes[elec + 2];
        let hn = self.nodes[elec + 3];
        let g_h = 1.0 / self.r_heater;
        if let Some(a) = hp {
            mat.a[a][a] += g_h;
            if let Some(c) = hn {
                mat.a[a][c] -= g_h;
            }
        }
        if let Some(c) = hn {
            mat.a[c][c] += g_h;
            if let Some(a) = hp {
                mat.a[c][a] -= g_h;
            }
        }
        // Optical: per-channel branch equations using cached rotation.
        for k in 0..n {
            let in_re_fw = self.nodes[wpc * k];
            let in_im_fw = self.nodes[wpc * k + 1];
            let in_l = self.nodes[wpc * k + lam];
            let out_re_fw = self.nodes[out_base + wpc * k];
            let out_im_fw = self.nodes[out_base + wpc * k + 1];
            let out_l = self.nodes[out_base + wpc * k + lam];
            let c = self.c_cached[k];
            let s = self.s_cached[k];
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k,
                out_re_fw,
                &[(in_re_fw, -c), (in_im_fw, -s)],
            );
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k + 1,
                out_im_fw,
                &[(in_re_fw, s), (in_im_fw, -c)],
            );
            stamp_potential_eq(
                mat,
                &self.branches,
                bpc * k + (bpc - 1),
                out_l,
                &[(in_l, -1.0)],
            );
            if wpc == 5 {
                let in_re_bw = self.nodes[wpc * k + 2];
                let in_im_bw = self.nodes[wpc * k + 3];
                let out_re_bw = self.nodes[out_base + wpc * k + 2];
                let out_im_bw = self.nodes[out_base + wpc * k + 3];
                stamp_potential_eq(
                    mat,
                    &self.branches,
                    bpc * k + 2,
                    in_re_bw,
                    &[(out_re_bw, -c), (out_im_bw, -s)],
                );
                stamp_potential_eq(
                    mat,
                    &self.branches,
                    bpc * k + 3,
                    in_im_bw,
                    &[(out_re_bw, s), (out_im_bw, -c)],
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

// ────────────────────────────────────────────────────────────────────────
// Native PN+thermal phase shifter, depletion-mode (L2a-with-heater)
// Combines fc_pn_ps_cap (C_j(V), reverse-bias FCA via da_dv) with a heater
// pair (heat_p, heat_n) for thermal trim/biasing on top of the EO shift.
// ────────────────────────────────────────────────────────────────────────

/// Lateral PN modulator + thermal trim heater.  Operates in the reverse-bias
/// regime (V_pn ≤ 0); has depletion C_j(V) and a linear da_dv loss-vs-bias
/// like `fc_pn_ps_cap`, plus heater terminals identical to `fc_pn_th_ps`.
///
/// Terminal layout (N optical channels): 2·wpc·N + 4
///   in.0.re … in.{N-1}.λ, out.0.re … out.{N-1}.λ, anode, cathode, heat_p, heat_n
pub struct NativePnThermalPhaseShifterCap {
    length_m: f64,
    n_eff: f64,
    n_g: f64,
    wl_ref_m: f64,
    pin_at_ref: bool,
    // PN side
    dn_dv: f64,
    g_pn: f64,
    alpha_neper_m: f64,
    c_j0: f64,
    v_bi: f64,
    m_j: f64,
    da_dv: f64,
    // Heater side
    r_heater: f64,
    p_pi_th: f64,
    n_channels: usize,
    wpc: usize,
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
    c_cached: Vec<f64>,
    s_cached: Vec<f64>,
    c_j_cached: f64,
    alpha_eff_neper_m: f64,
}

impl NativePnThermalPhaseShifterCap {
    pub fn new() -> Self {
        // Defaults from pn_modulator.py extraction (5e17/5e17 lateral PN,
        // 100 nm offset, 300 K).
        Self {
            length_m: 1e-3,
            n_eff: 2.7654,
            n_g: 4.02,
            wl_ref_m: 1.55e-6,
            pin_at_ref: false,
            dn_dv: 5.024e-5, // depletion-mode linear coeff
            g_pn: 1e-3,
            alpha_neper_m: 29.78,  // 2.59 dB/cm at V=0 (FCA, sim)
            c_j0: 1.375e-16 * 1e3, // F per µm × µm → F-ish per 1 mm L (= 1.375e-13 F/m); see `length_m` default
            v_bi: 0.917,           // V_bi @ N_A=N_D=5e17, 300 K
            m_j: 0.5,
            da_dv: 7.83, // slope: Δα/ΔV ≈ 7.8 Np/m per V
            r_heater: 1000.0,
            p_pi_th: 10e-3,
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
            c_cached: Vec::new(),
            s_cached: Vec::new(),
            c_j_cached: 1.375e-16,
            alpha_eff_neper_m: 29.78,
        }
    }
}

impl Device for NativePnThermalPhaseShifterCap {
    fn num_terminals(&self) -> usize {
        self.nodes.len()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
        self.alpha_eff_neper_m = self.alpha_neper_m;
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        let stride = 2 * wpc;
        assert!(
            terminals.len() >= stride + 4 && (terminals.len() - 4) % stride == 0,
            "fc_pn_th_ps_cap: terminal count must be {stride}·N + 4 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 4) / stride;
        self.n_channels = n;
        self.nodes = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 };
        self.branches = vec![None; bpc * n];
        self.c_cached = vec![1.0; n];
        self.s_cached = vec![0.0; n];
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
            "l_um" => {
                self.length_m = value * 1e-6;
                true
            }
            "l_m" | "length" => {
                self.length_m = value;
                true
            }
            "n_g" => {
                self.n_g = value;
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
            "dn_dv" => {
                self.dn_dv = value;
                true
            }
            "g_pn" => {
                self.g_pn = value;
                true
            }
            "v_pi_l" => {
                if value > 0.0 {
                    self.dn_dv = self.wl_ref_m / (2.0 * value);
                }
                true
            }
            "alpha_db_cm" => {
                self.alpha_neper_m = dB_per_cm_to_neper_per_m(value);
                true
            }
            "c_j0" => {
                self.c_j0 = value.max(0.0);
                true
            }
            "v_bi" => {
                self.v_bi = value.max(1e-3);
                true
            }
            "m_j" => {
                self.m_j = value.clamp(0.0, 0.99);
                true
            }
            "da_dv" => {
                self.da_dv = value;
                true
            }
            "r_heater" | "r" => {
                self.r_heater = value;
                true
            }
            "p_pi" | "p_pi_w" | "p_pi_th" => {
                self.p_pi_th = value;
                true
            }
            "pin_at_ref" => {
                self.pin_at_ref = value != 0.0;
                true
            }
            _ => false,
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let elec = 2 * wpc * n;
        let v_a = self.nodes[elec].map_or(0.0, |i| x[i]);
        let v_c = self.nodes[elec + 1].map_or(0.0, |i| x[i]);
        let v_hp = self.nodes[elec + 2].map_or(0.0, |i| x[i]);
        let v_hn = self.nodes[elec + 3].map_or(0.0, |i| x[i]);
        let v_pn = v_a - v_c;
        let v_h = v_hp - v_hn;
        let p_heat = v_h * v_h / self.r_heater;
        let phi_th = std::f64::consts::PI * p_heat / self.p_pi_th;

        // C_j(V) with linear continuation past the V_bi/2 knee.
        let v_knee = 0.5 * self.v_bi;
        self.c_j_cached = if v_pn < v_knee {
            self.c_j0 / (1.0 - v_pn / self.v_bi).powf(self.m_j)
        } else {
            let c_knee = self.c_j0 / (1.0 - v_knee / self.v_bi).powf(self.m_j);
            let dc_dv = c_knee * self.m_j / (self.v_bi - v_knee);
            c_knee + dc_dv * (v_pn - v_knee)
        };
        // Bias-dependent loss: α(V) = α_0 + da_dv·max(0, -V).
        let v_rev = (-v_pn).max(0.0);
        self.alpha_eff_neper_m = self.alpha_neper_m + self.da_dv * v_rev;
        let t_amp = (-self.alpha_eff_neper_m * self.length_m / 2.0).exp();

        let two_pi = 2.0 * std::f64::consts::PI;
        let lam = wpc - 1;
        let phi_ref = if self.pin_at_ref {
            two_pi * self.n_eff * self.length_m / self.wl_ref_m
        } else {
            0.0
        };
        for k in 0..n {
            let lambda = match self.nodes[wpc * k + lam] {
                Some(i) => {
                    let v = x[i];
                    if v.abs() > 1e-9 {
                        v
                    } else {
                        self.wl_ref_m
                    }
                }
                None => self.wl_ref_m,
            };
            let n_eff_lam = n_eff_at_lambda(self.n_eff, self.n_g, self.wl_ref_m, lambda);
            let phi_abs = two_pi * n_eff_lam * self.length_m / lambda;
            let phi_prop = phi_abs - phi_ref;
            let phi_eo = two_pi * self.length_m * self.dn_dv * v_pn / lambda;
            let phi = phi_prop + phi_eo + phi_th;
            self.c_cached[k] = t_amp * phi.cos();
            self.s_cached[k] = t_amp * phi.sin();
        }
    }

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        stamp_pn_ths_jacobian(
            mat,
            &self.nodes,
            &self.branches,
            self.n_channels,
            self.wpc,
            self.g_pn,
            self.r_heater,
            &self.c_cached,
            &self.s_cached,
        );
    }

    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> {
        let elec_base = 2 * self.wpc * self.n_channels;
        let anode = self.nodes.get(elec_base).copied().flatten();
        let cath = self.nodes.get(elec_base + 1).copied().flatten();
        vec![ReactiveBranchSpec {
            kind: ReactiveKind::Capacitor,
            pos: anode,
            neg: cath,
            value: self.c_j_cached,
        }]
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }
}

// ────────────────────────────────────────────────────────────────────────
// Native PN phase shifter, forward-bias injection (L2b)
// ────────────────────────────────────────────────────────────────────────

/// Forward-bias carrier-injection PN phase shifter.  Suitable for V_pn ∈
/// [0, ~0.8 V] only — does not model depletion physics.  Use this when the
/// design relies on injected free-carrier index change (typical VOA / slow
/// thermal-trim modulator with carrier dynamics on top).
///
/// Physics differences from fc_pn_ps_cap:
///   - Shockley forward I-V:   I = I_s · (exp(V/(n·V_T)) − 1), linearised at
///     the operating point with g_d = (I + I_s)/(n·V_T).
///   - Diffusion capacitance:  C_d = τ_carrier · g_d  (replaces depletion).
///   - Exponential injection Δn(V):  Δn_inj = K_inj · (exp(V/(n·V_T)) − 1),
///     where K_inj is the linearised forward-bias dn_dv coefficient scaled
///     to give the requested fractional dn at V = ~3·V_T.
///   - Exponential injection loss Δα(V):  same exp shape with K_alpha.
///
/// Terminal layout: 2·wpc·N + 2  (same as `fc_pn_ps`).
pub struct NativePnPhaseShifterInj {
    length_m: f64,
    n_eff: f64,
    n_g: f64,
    wl_ref_m: f64,
    pin_at_ref: bool,
    // Forward-bias diode
    i_sat: f64,
    n_diode: f64,
    tau_carrier: f64,
    // EO + FCA injection coefficients
    dn_dv_inj: f64,     // K_inj = Δn_eff(V→Vt) / 1 (slope-equivalent)
    da_dv_inj: f64,     // Np/m per (exp(V/Vt)-1)
    alpha_neper_m: f64, // background propagation loss
    n_channels: usize,
    wpc: usize,
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
    c_cached: Vec<f64>,
    s_cached: Vec<f64>,
    g_d_cached: f64,
    i_eq_cached: f64,
    c_d_cached: f64,
    v_pn_op: f64,
}

impl NativePnPhaseShifterInj {
    pub fn new() -> Self {
        Self {
            length_m: 1e-3,
            n_eff: 2.7654,
            n_g: 4.02,
            wl_ref_m: 1.55e-6,
            pin_at_ref: false,
            i_sat: 1e-12,
            n_diode: 1.05,
            tau_carrier: 10e-9,
            dn_dv_inj: 1.311e-4, // forward small-signal coeff (sim)
            da_dv_inj: 150.0,    // exp prefactor for FCA injection (Np/m)
            alpha_neper_m: dB_per_cm_to_neper_per_m(1.0),
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
            c_cached: Vec::new(),
            s_cached: Vec::new(),
            g_d_cached: 1e-9,
            i_eq_cached: 0.0,
            c_d_cached: 0.0,
            v_pn_op: 0.0,
        }
    }
}

impl Device for NativePnPhaseShifterInj {
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
            terminals.len() >= stride + 2 && (terminals.len() - 2) % stride == 0,
            "fc_pn_ps_inj: terminal count must be {stride}·N + 2 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 2) / stride;
        self.n_channels = n;
        self.nodes = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 };
        self.branches = vec![None; bpc * n];
        self.c_cached = vec![1.0; n];
        self.s_cached = vec![0.0; n];
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
            "l_um" => {
                self.length_m = value * 1e-6;
                true
            }
            "l_m" | "length" => {
                self.length_m = value;
                true
            }
            "n_g" => {
                self.n_g = value;
                true
            }
            "n_eff" => {
                self.n_eff = value;
                true
            }
            "wl_ref_nm" => {
                self.wl_ref_m = value * 1e-9;
                true
            }
            "wl_ref_m" => {
                self.wl_ref_m = value;
                true
            }
            "i_sat" | "is" => {
                self.i_sat = value.max(0.0);
                true
            }
            "n_diode" | "n" => {
                self.n_diode = value.max(0.5);
                true
            }
            "tau_carrier" | "tau" => {
                self.tau_carrier = value.max(0.0);
                true
            }
            "dn_dv_inj" | "dn_dv" => {
                self.dn_dv_inj = value;
                true
            }
            "da_dv_inj" | "da_dv" => {
                self.da_dv_inj = value;
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

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, ctx: &SimContext) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let elec = 2 * wpc * n;
        let v_a = self.nodes[elec].map_or(0.0, |i| x[i]);
        let v_c = self.nodes[elec + 1].map_or(0.0, |i| x[i]);
        let v_pn = v_a - v_c;
        self.v_pn_op = v_pn;

        // Shockley I-V (clamp exp argument for NR stability)
        let vt = ctx.vt() * self.n_diode;
        let arg = (v_pn / vt).min(40.0).max(-40.0);
        let e = arg.exp();
        let i_diode = self.i_sat * (e - 1.0);
        self.g_d_cached = self.i_sat * e / vt;
        // Norton equivalent: I_eq = I(V_op) - g · V_op
        self.i_eq_cached = i_diode - self.g_d_cached * v_pn;
        // Diffusion capacitance (forward bias only — small under reverse)
        self.c_d_cached = self.tau_carrier * self.g_d_cached;

        // Optical: injection-driven Δn and Δα, only forward bias contributes.
        let inj = (e - 1.0).max(0.0);
        let alpha_eff = self.alpha_neper_m + self.da_dv_inj * inj;
        let t_amp = (-alpha_eff * self.length_m / 2.0).exp();

        let two_pi = 2.0 * std::f64::consts::PI;
        let lam = wpc - 1;
        let phi_ref = if self.pin_at_ref {
            two_pi * self.n_eff * self.length_m / self.wl_ref_m
        } else {
            0.0
        };
        for k in 0..n {
            let lambda = match self.nodes[wpc * k + lam] {
                Some(i) => {
                    let v = x[i];
                    if v.abs() > 1e-9 {
                        v
                    } else {
                        self.wl_ref_m
                    }
                }
                None => self.wl_ref_m,
            };
            let n_eff_lam = n_eff_at_lambda(self.n_eff, self.n_g, self.wl_ref_m, lambda);
            let phi_abs = two_pi * n_eff_lam * self.length_m / lambda;
            let phi_prop = phi_abs - phi_ref;
            // For injection, we add a Soref-Bennett-shaped Δn (negative for
            // more carriers → less n).  Linear coefficient `dn_dv_inj` is the
            // slope at V=0; full exponential form scales as (e-1)/1V.
            let phi_eo = -two_pi * self.length_m * self.dn_dv_inj * inj / lambda;
            let phi = phi_prop + phi_eo;
            self.c_cached[k] = t_amp * phi.cos();
            self.s_cached[k] = t_amp * phi.sin();
        }
    }

    fn load_residual(&self, b: &mut [f64]) {
        let elec = 2 * self.wpc * self.n_channels;
        if let Some(a) = self.nodes[elec] {
            b[a] -= self.i_eq_cached;
        }
        if let Some(c) = self.nodes[elec + 1] {
            b[c] += self.i_eq_cached;
        }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        stamp_pn_optical(
            mat,
            &self.nodes,
            &self.branches,
            self.n_channels,
            self.wpc,
            &self.c_cached,
            &self.s_cached,
        );
        // Electrical: shared diode small-signal conductance
        let elec = 2 * self.wpc * self.n_channels;
        stamp_resistor(mat, self.nodes[elec], self.nodes[elec + 1], self.g_d_cached);
    }

    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> {
        let elec = 2 * self.wpc * self.n_channels;
        let a = self.nodes.get(elec).copied().flatten();
        let c = self.nodes.get(elec + 1).copied().flatten();
        vec![ReactiveBranchSpec {
            kind: ReactiveKind::Capacitor,
            pos: a,
            neg: c,
            value: self.c_d_cached,
        }]
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }
}

// ────────────────────────────────────────────────────────────────────────
// Native PN+thermal phase shifter, forward injection (L2b-with-heater)
// ────────────────────────────────────────────────────────────────────────

pub struct NativePnThermalPhaseShifterInj {
    inj: NativePnPhaseShifterInj,
    r_heater: f64,
    p_pi_th: f64,
}

impl NativePnThermalPhaseShifterInj {
    pub fn new() -> Self {
        Self {
            inj: NativePnPhaseShifterInj::new(),
            r_heater: 1000.0,
            p_pi_th: 10e-3,
        }
    }
}

impl Device for NativePnThermalPhaseShifterInj {
    fn num_terminals(&self) -> usize {
        self.inj.nodes.len()
    }
    fn setup_model(&mut self, ctx: &SimContext) {
        self.inj.setup_model(ctx);
    }
    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        // Same as `_inj` but expects 4 extra electrical pins.
        let wpc = ctx.wires_per_channel();
        self.inj.wpc = wpc;
        let stride = 2 * wpc;
        assert!(
            terminals.len() >= stride + 4 && (terminals.len() - 4) % stride == 0,
            "fc_pn_th_ps_inj: terminal count must be {stride}·N + 4 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 4) / stride;
        self.inj.n_channels = n;
        self.inj.nodes = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 };
        self.inj.branches = vec![None; bpc * n];
        self.inj.c_cached = vec![1.0; n];
        self.inj.s_cached = vec![0.0; n];
    }
    fn num_extra_nodes(&self) -> usize {
        self.inj.branches.len()
    }
    fn bind_extra_nodes(&mut self, idx: usize) {
        self.inj.bind_extra_nodes(idx);
    }
    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "r_heater" | "r" => {
                self.r_heater = value;
                true
            }
            "p_pi" | "p_pi_w" | "p_pi_th" => {
                self.p_pi_th = value;
                true
            }
            _ => self.inj.set_real_param(name, value),
        }
    }
    fn eval(&mut self, x: &[f64], flags: EvalFlags, ctx: &SimContext) {
        self.inj.eval(x, flags, ctx);
        // Add thermal phase from heater (heat_p, heat_n live at the END of nodes).
        let n = self.inj.n_channels;
        let wpc = self.inj.wpc;
        let elec = 2 * wpc * n;
        let v_hp = self.inj.nodes[elec + 2].map_or(0.0, |i| x[i]);
        let v_hn = self.inj.nodes[elec + 3].map_or(0.0, |i| x[i]);
        let v_h = v_hp - v_hn;
        let phi_th = std::f64::consts::PI * (v_h * v_h / self.r_heater) / self.p_pi_th;
        // Rotate the cached c, s by phi_th (no re-eval of EO part needed):
        let cth = phi_th.cos();
        let sth = phi_th.sin();
        for k in 0..n {
            let c = self.inj.c_cached[k];
            let s = self.inj.s_cached[k];
            self.inj.c_cached[k] = c * cth - s * sth;
            self.inj.s_cached[k] = c * sth + s * cth;
        }
    }
    fn load_residual(&self, b: &mut [f64]) {
        self.inj.load_residual(b);
    }
    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        self.inj.load_jacobian(mat);
        // Heater resistor between heat_p and heat_n.
        let elec = 2 * self.inj.wpc * self.inj.n_channels;
        stamp_resistor(
            mat,
            self.inj.nodes[elec + 2],
            self.inj.nodes[elec + 3],
            1.0 / self.r_heater,
        );
    }
    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> {
        self.inj.reactive_branches()
    }
    fn load_residual_tran(&self, b: &mut [f64], a: f64) {
        self.inj.load_residual_tran(b, a);
    }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, a: f64) {
        self.inj.load_jacobian_tran(mat, a);
        let elec = 2 * self.inj.wpc * self.inj.n_channels;
        stamp_resistor(
            mat,
            self.inj.nodes[elec + 2],
            self.inj.nodes[elec + 3],
            1.0 / self.r_heater,
        );
    }
}

// ────────────────────────────────────────────────────────────────────────
// Native PN phase shifter, combined depletion + injection + TPA + static
// self-heating (L3).
// ────────────────────────────────────────────────────────────────────────

/// "Full" PN modulator — captures both bias regimes, two-photon absorption,
/// and static thermal self-heating from absorbed optical power.  Intended
/// for steady-state and slow transient analysis of high-performance
/// modulators; carrier and thermal dynamics belong in `_carrier` (L4).
///
/// Physics included:
///   - Reverse-bias depletion Δn_rev(V) = dn_dv_rev · V  (V ≤ 0)
///   - Forward-bias injection Δn_fwd(V) = -dn_dv_inj · (exp(V/V_T)-1)  (V ≥ 0)
///   - Depletion C_j(V) for V ≤ 0; diffusion C_d for V > 0
///   - Reverse FCA loss α_rev(V) = da_dv_rev · max(0,-V)
///   - Forward FCA loss α_fwd(V) = da_dv_inj · (exp(V/V_T)-1)
///   - TPA loss α_TPA = β_TPA · (|A|²/A_eff)
///   - Static self-heating: ΔT_ss = R_th · α_total · L · |A|²;  Δn_th = dn_dT · ΔT
pub struct NativePnPhaseShifterFull {
    length_m: f64,
    n_eff: f64,
    n_g: f64,
    wl_ref_m: f64,
    pin_at_ref: bool,
    // Reverse
    dn_dv_rev: f64,
    da_dv_rev: f64,
    c_j0: f64,
    v_bi: f64,
    m_j: f64,
    // Forward
    i_sat: f64,
    n_diode: f64,
    tau_carrier: f64,
    dn_dv_inj: f64,
    da_dv_inj: f64,
    // Common
    alpha_neper_m: f64,
    // TPA + thermal
    beta_tpa_m_per_w: f64,
    a_eff_m2: f64,
    r_th_k_per_w: f64,
    dn_dt: f64,
    n_channels: usize,
    wpc: usize,
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
    c_cached: Vec<f64>,
    s_cached: Vec<f64>,
    g_pn_cached: f64,
    i_eq_cached: f64,
    c_eff_cached: f64,
}

impl NativePnPhaseShifterFull {
    pub fn new() -> Self {
        Self {
            length_m: 1e-3,
            n_eff: 2.7654,
            n_g: 4.02,
            wl_ref_m: 1.55e-6,
            pin_at_ref: false,
            dn_dv_rev: 5.024e-5,
            da_dv_rev: 7.83,
            c_j0: 1.375e-13, // F at default L
            v_bi: 0.917,
            m_j: 0.5,
            i_sat: 1e-12,
            n_diode: 1.05,
            tau_carrier: 10e-9,
            dn_dv_inj: 1.311e-4,
            da_dv_inj: 150.0,
            alpha_neper_m: dB_per_cm_to_neper_per_m(1.0),
            beta_tpa_m_per_w: 7.9e-12,
            a_eff_m2: 1.257e-13,
            r_th_k_per_w: 0.0, // user must set
            dn_dt: 1.86e-4,    // crystalline Si
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
            c_cached: Vec::new(),
            s_cached: Vec::new(),
            g_pn_cached: 1e-9,
            i_eq_cached: 0.0,
            c_eff_cached: 1.375e-13,
        }
    }
}

impl Device for NativePnPhaseShifterFull {
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
            terminals.len() >= stride + 2 && (terminals.len() - 2) % stride == 0,
            "fc_pn_ps_full: terminal count must be {stride}·N + 2 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 2) / stride;
        self.n_channels = n;
        self.nodes = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 };
        self.branches = vec![None; bpc * n];
        self.c_cached = vec![1.0; n];
        self.s_cached = vec![0.0; n];
    }
    fn num_extra_nodes(&self) -> usize {
        self.branches.len()
    }
    fn bind_extra_nodes(&mut self, idx: usize) {
        for i in 0..self.branches.len() {
            self.branches[i] = Some(idx + i);
        }
    }
    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "l_um" => {
                self.length_m = value * 1e-6;
                true
            }
            "l_m" | "length" => {
                self.length_m = value;
                true
            }
            "n_g" => {
                self.n_g = value;
                true
            }
            "n_eff" => {
                self.n_eff = value;
                true
            }
            "wl_ref_nm" => {
                self.wl_ref_m = value * 1e-9;
                true
            }
            "wl_ref_m" => {
                self.wl_ref_m = value;
                true
            }
            "dn_dv_rev" | "dn_dv" => {
                self.dn_dv_rev = value;
                true
            }
            "da_dv_rev" | "da_dv" => {
                self.da_dv_rev = value;
                true
            }
            "c_j0" => {
                self.c_j0 = value.max(0.0);
                true
            }
            "v_bi" => {
                self.v_bi = value.max(1e-3);
                true
            }
            "m_j" => {
                self.m_j = value.clamp(0.0, 0.99);
                true
            }
            "i_sat" | "is" => {
                self.i_sat = value.max(0.0);
                true
            }
            "n_diode" | "n" => {
                self.n_diode = value.max(0.5);
                true
            }
            "tau_carrier" | "tau" => {
                self.tau_carrier = value.max(0.0);
                true
            }
            "dn_dv_inj" => {
                self.dn_dv_inj = value;
                true
            }
            "da_dv_inj" => {
                self.da_dv_inj = value;
                true
            }
            "alpha_db_cm" => {
                self.alpha_neper_m = dB_per_cm_to_neper_per_m(value);
                true
            }
            "beta_tpa" | "beta_tpa_m_per_w" => {
                self.beta_tpa_m_per_w = value;
                true
            }
            "a_eff_m2" => {
                self.a_eff_m2 = value.max(1e-20);
                true
            }
            "a_eff_um2" => {
                self.a_eff_m2 = value.max(1e-8) * 1e-12;
                true
            }
            "r_th" | "r_th_k_per_w" => {
                self.r_th_k_per_w = value.max(0.0);
                true
            }
            "dn_dt" => {
                self.dn_dt = value;
                true
            }
            "pin_at_ref" => {
                self.pin_at_ref = value != 0.0;
                true
            }
            _ => false,
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, ctx: &SimContext) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let elec = 2 * wpc * n;
        let v_a = self.nodes[elec].map_or(0.0, |i| x[i]);
        let v_c = self.nodes[elec + 1].map_or(0.0, |i| x[i]);
        let v_pn = v_a - v_c;

        // Electrical: Shockley I-V across the whole regime.
        let vt = ctx.vt() * self.n_diode;
        let e = (v_pn / vt).min(40.0).max(-40.0).exp();
        let i_diode = self.i_sat * (e - 1.0);
        let g_d = self.i_sat * e / vt;
        self.g_pn_cached = g_d.max(1e-15);
        self.i_eq_cached = i_diode - g_d * v_pn;

        // Capacitance: piecewise C_j (reverse) vs C_d (forward).
        let c_j_v = {
            let v_knee = 0.5 * self.v_bi;
            if v_pn < v_knee {
                self.c_j0 / (1.0 - v_pn / self.v_bi).powf(self.m_j)
            } else {
                let c_knee = self.c_j0 / (1.0 - v_knee / self.v_bi).powf(self.m_j);
                let dc_dv = c_knee * self.m_j / (self.v_bi - v_knee);
                c_knee + dc_dv * (v_pn - v_knee)
            }
        };
        let c_d_v = self.tau_carrier * g_d;
        self.c_eff_cached = c_j_v + c_d_v;

        // Optical: sum of contributions.  Compute |A|² at the input for TPA
        // and self-heating (use channel-0; aggregating across WDM is left to L4).
        let v_re_0 = self.nodes[0].map_or(0.0, |i| x[i]);
        let v_im_0 = self.nodes[1].map_or(0.0, |i| x[i]);
        let intensity_w = (v_re_0 * v_re_0 + v_im_0 * v_im_0).max(0.0);
        let alpha_tpa = self.beta_tpa_m_per_w * intensity_w / self.a_eff_m2;
        let inj = (e - 1.0).max(0.0);
        let v_rev = (-v_pn).max(0.0);
        let alpha_fca = self.alpha_neper_m + self.da_dv_rev * v_rev + self.da_dv_inj * inj;
        let alpha_total = alpha_fca + alpha_tpa;
        let t_amp = (-alpha_total * self.length_m / 2.0).exp();

        // Self-heating Δn (static): ΔT = R_th · P_abs;  P_abs ≈ α · L · |A|².
        let p_abs = alpha_total * self.length_m * intensity_w;
        let dn_self = self.dn_dt * self.r_th_k_per_w * p_abs;

        let two_pi = 2.0 * std::f64::consts::PI;
        let lam = wpc - 1;
        let phi_ref = if self.pin_at_ref {
            two_pi * self.n_eff * self.length_m / self.wl_ref_m
        } else {
            0.0
        };
        for k in 0..n {
            let lambda = match self.nodes[wpc * k + lam] {
                Some(i) => {
                    let v = x[i];
                    if v.abs() > 1e-9 {
                        v
                    } else {
                        self.wl_ref_m
                    }
                }
                None => self.wl_ref_m,
            };
            let n_eff_lam = n_eff_at_lambda(self.n_eff, self.n_g, self.wl_ref_m, lambda);
            let phi_abs = two_pi * n_eff_lam * self.length_m / lambda;
            let phi_prop = phi_abs - phi_ref;
            let phi_eo_rev = two_pi * self.length_m * self.dn_dv_rev * v_pn / lambda;
            let phi_eo_inj = -two_pi * self.length_m * self.dn_dv_inj * inj / lambda;
            let phi_self = two_pi * self.length_m * dn_self / lambda;
            let phi = phi_prop + phi_eo_rev + phi_eo_inj + phi_self;
            self.c_cached[k] = t_amp * phi.cos();
            self.s_cached[k] = t_amp * phi.sin();
        }
    }

    fn load_residual(&self, b: &mut [f64]) {
        let elec = 2 * self.wpc * self.n_channels;
        if let Some(a) = self.nodes[elec] {
            b[a] -= self.i_eq_cached;
        }
        if let Some(c) = self.nodes[elec + 1] {
            b[c] += self.i_eq_cached;
        }
    }
    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        stamp_pn_optical(
            mat,
            &self.nodes,
            &self.branches,
            self.n_channels,
            self.wpc,
            &self.c_cached,
            &self.s_cached,
        );
        let elec = 2 * self.wpc * self.n_channels;
        stamp_resistor(
            mat,
            self.nodes[elec],
            self.nodes[elec + 1],
            self.g_pn_cached,
        );
    }
    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> {
        let elec = 2 * self.wpc * self.n_channels;
        let a = self.nodes.get(elec).copied().flatten();
        let c = self.nodes.get(elec + 1).copied().flatten();
        vec![ReactiveBranchSpec {
            kind: ReactiveKind::Capacitor,
            pos: a,
            neg: c,
            value: self.c_eff_cached,
        }]
    }
    fn load_residual_tran(&self, b: &mut [f64], _a: f64) {
        self.load_residual(b);
    }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _a: f64) {
        self.load_jacobian(mat);
    }
}

// ────────────────────────────────────────────────────────────────────────
// Native PN+thermal phase shifter, full (L3-with-heater)
// ────────────────────────────────────────────────────────────────────────

pub struct NativePnThermalPhaseShifterFull {
    full: NativePnPhaseShifterFull,
    r_heater: f64,
    p_pi_th: f64,
}

impl NativePnThermalPhaseShifterFull {
    pub fn new() -> Self {
        Self {
            full: NativePnPhaseShifterFull::new(),
            r_heater: 1000.0,
            p_pi_th: 10e-3,
        }
    }
}

impl Device for NativePnThermalPhaseShifterFull {
    fn num_terminals(&self) -> usize {
        self.full.nodes.len()
    }
    fn setup_model(&mut self, ctx: &SimContext) {
        self.full.setup_model(ctx);
    }
    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.full.wpc = wpc;
        let stride = 2 * wpc;
        assert!(
            terminals.len() >= stride + 4 && (terminals.len() - 4) % stride == 0,
            "fc_pn_th_ps_full: terminal count must be {stride}·N + 4 (wpc={wpc}); got {}",
            terminals.len()
        );
        let n = (terminals.len() - 4) / stride;
        self.full.n_channels = n;
        self.full.nodes = terminals.to_vec();
        let bpc = if wpc == 5 { 5 } else { 3 };
        self.full.branches = vec![None; bpc * n];
        self.full.c_cached = vec![1.0; n];
        self.full.s_cached = vec![0.0; n];
    }
    fn num_extra_nodes(&self) -> usize {
        self.full.branches.len()
    }
    fn bind_extra_nodes(&mut self, idx: usize) {
        self.full.bind_extra_nodes(idx);
    }
    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "r_heater" | "r" => {
                self.r_heater = value;
                true
            }
            "p_pi" | "p_pi_w" | "p_pi_th" => {
                self.p_pi_th = value;
                true
            }
            _ => self.full.set_real_param(name, value),
        }
    }
    fn eval(&mut self, x: &[f64], flags: EvalFlags, ctx: &SimContext) {
        self.full.eval(x, flags, ctx);
        let n = self.full.n_channels;
        let wpc = self.full.wpc;
        let elec = 2 * wpc * n;
        let v_hp = self.full.nodes[elec + 2].map_or(0.0, |i| x[i]);
        let v_hn = self.full.nodes[elec + 3].map_or(0.0, |i| x[i]);
        let v_h = v_hp - v_hn;
        let phi_th = std::f64::consts::PI * (v_h * v_h / self.r_heater) / self.p_pi_th;
        let cth = phi_th.cos();
        let sth = phi_th.sin();
        for k in 0..n {
            let c = self.full.c_cached[k];
            let s = self.full.s_cached[k];
            self.full.c_cached[k] = c * cth - s * sth;
            self.full.s_cached[k] = c * sth + s * cth;
        }
    }
    fn load_residual(&self, b: &mut [f64]) {
        self.full.load_residual(b);
    }
    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        self.full.load_jacobian(mat);
        let elec = 2 * self.full.wpc * self.full.n_channels;
        stamp_resistor(
            mat,
            self.full.nodes[elec + 2],
            self.full.nodes[elec + 3],
            1.0 / self.r_heater,
        );
    }
    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> {
        self.full.reactive_branches()
    }
    fn load_residual_tran(&self, b: &mut [f64], a: f64) {
        self.full.load_residual_tran(b, a);
    }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, a: f64) {
        self.full.load_jacobian_tran(mat, a);
        let elec = 2 * self.full.wpc * self.full.n_channels;
        stamp_resistor(
            mat,
            self.full.nodes[elec + 2],
            self.full.nodes[elec + 3],
            1.0 / self.r_heater,
        );
    }
}
