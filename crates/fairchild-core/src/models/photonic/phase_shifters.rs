use super::{dB_per_cm_to_neper_per_m, n_eff_at_lambda, stamp_pn_optical, stamp_resistor};
use crate::device::{Device, EvalFlags, NodeId, ReactiveBranchSpec, ReactiveKind, SimContext};
use crate::mna::MnaMatrix;

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
    // Ohmic series resistance in the PN contact/via stack (Ω).
    // When non-zero, the terminal voltage V_pn ≠ junction voltage V_j;
    // Newton-Raphson in eval() solves V_j implicitly each iteration.
    r_series: f64,
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

impl Default for NativePnPhaseShifterFull {
    fn default() -> Self {
        Self::new()
    }
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
            r_series: 0.0,
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
            terminals.len() >= stride + 2 && (terminals.len() - 2).is_multiple_of(stride),
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
            "r_series" => {
                self.r_series = value.max(0.0);
                true
            }
            _ => false,
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, ctx: &SimContext) {
        let n = self.n_channels;
        let wpc = self.wpc;
        let elec = 2 * wpc * n;

        // ── 1. Electrical: junction voltage ───────────────────────────────
        // Terminal voltages from the MNA solution vector.
        let v_a = self.nodes[elec].map_or(0.0, |i| x[i]);
        let v_c = self.nodes[elec + 1].map_or(0.0, |i| x[i]);
        // v_pn is the voltage applied *at the device terminals* (anode − cathode).
        // With non-zero series resistance R_s this differs from the actual
        // junction voltage v_junc: v_pn = v_junc + R_s · I_d(v_junc).
        let v_pn = v_a - v_c;

        // ── 2. Electrical: Shockley I-V (with series resistance if set) ───
        // Thermal voltage scaled by ideality: V_T = n_diode · k·T/q.
        let vt = ctx.vt() * self.n_diode;

        // Solve for the junction voltage v_junc by Newton-Raphson on
        //   F(v_j) = v_j + R_s · I_sat · (exp(v_j/V_T) − 1) − v_pn = 0
        // When R_s = 0 the solution is trivially v_j = v_pn (one iteration).
        let v_junc = if self.r_series <= 0.0 {
            v_pn
        } else {
            let mut vj = v_pn; // initial guess: ignore drop across R_s
            for _ in 0..50 {
                let arg = (vj / vt).clamp(-40.0, 40.0);
                let e = arg.exp();
                let id = self.i_sat * (e - 1.0);
                // Residual: how far we are from satisfying the KVL equation.
                let f = vj + self.r_series * id - v_pn;
                // Derivative dF/dv_j = 1 + R_s · (dI_d/dv_j) = 1 + R_s · g_d.
                let gd = self.i_sat * e / vt;
                let df = 1.0 + self.r_series * gd;
                let delta = f / df;
                vj -= delta;
                if delta.abs() < 1e-12 {
                    break;
                }
            }
            vj
        };

        // Evaluate diode quantities at the converged junction voltage.
        let arg = (v_junc / vt).clamp(-40.0, 40.0);
        let e = arg.exp();
        // Shockley current through the junction: I_d = I_sat · (exp(v_j/V_T) − 1).
        let i_diode = self.i_sat * (e - 1.0);
        // Small-signal conductance: g_d = dI_d/dv_j = I_sat · exp(v_j/V_T) / V_T.
        let g_d = self.i_sat * e / vt;

        // ── 3. MNA Norton stamp (accounts for series resistance) ──────────
        // With R_s the effective conductance seen from the terminals differs
        // from g_d.  Linearise: I_d ≈ I_d(v_j_op) + g_d · (v_j − v_j_op).
        // Substituting v_j = v_pn − R_s · I_d and solving:
        //   G_eff = g_d / (1 + g_d · R_s)
        // The Norton current is then: I_eq = I_d − G_eff · v_pn.
        let g_eff = g_d / (1.0 + g_d * self.r_series);
        self.g_pn_cached = g_eff.max(1e-15);
        self.i_eq_cached = i_diode - g_eff * v_pn;

        // ── 4. Piecewise junction capacitance C_j(v_j) + diffusion C_d ───
        // Depletion capacitance (abrupt/graded junction model):
        //   C_j(V) = C_j0 / (1 − V/V_bi)^m_j   for V < V_bi/2
        // Linearly continued above the V_bi/2 knee to avoid the singularity
        // when NR wanders into strong forward bias.
        let c_j_v = {
            let v_knee = 0.5 * self.v_bi;
            if v_junc < v_knee {
                self.c_j0 / (1.0 - v_junc / self.v_bi).powf(self.m_j)
            } else {
                // Tangent at V_bi/2: c(V_knee) + (dc/dV)|_knee · (V − V_knee).
                let c_knee = self.c_j0 / (1.0 - v_knee / self.v_bi).powf(self.m_j);
                let dc_dv = c_knee * self.m_j / (self.v_bi - v_knee);
                c_knee + dc_dv * (v_junc - v_knee)
            }
        };
        // Diffusion (transit-time) capacitance, dominant in strong forward bias:
        //   C_d = τ_carrier · g_d
        // Together they form the total effective cap for the reactive branch.
        let c_d_v = self.tau_carrier * g_d;
        self.c_eff_cached = c_j_v + c_d_v;

        // ── 5. Optical intensity (for TPA and self-heating) ───────────────
        // Use channel-0 forward amplitude squared as a proxy for intensity.
        // (WDM aggregation across channels is deferred to higher-level models.)
        let v_re_0 = self.nodes[0].map_or(0.0, |i| x[i]);
        let v_im_0 = self.nodes[1].map_or(0.0, |i| x[i]);
        let intensity_w = (v_re_0 * v_re_0 + v_im_0 * v_im_0).max(0.0);

        // ── 6. Loss: FCA + TPA ────────────────────────────────────────────
        // `inj` = exp(v_j/V_T) − 1 clamped to ≥ 0 (injection carrier density
        // is proportional to forward current; no carrier injection in reverse).
        let inj = (e - 1.0).max(0.0);
        // `v_rev` = |V_pn| under reverse bias, 0 in forward (depletion FCA).
        let v_rev = (-v_junc).max(0.0);
        // Two-photon absorption (TPA): α_TPA = β_TPA · I / A_eff.
        let alpha_tpa = self.beta_tpa_m_per_w * intensity_w / self.a_eff_m2;
        // Total loss coefficient (Np/m):
        //   α_total = α_0 + α_rev(V) + α_inj(V) + α_TPA(I)
        //   α_rev  = da_dv_rev · max(0, −v_j)   (depletion free-carrier absorption)
        //   α_inj  = da_dv_inj · (exp(v_j/V_T)−1) (injection FCA, forward bias)
        let alpha_fca = self.alpha_neper_m + self.da_dv_rev * v_rev + self.da_dv_inj * inj;
        let alpha_total = alpha_fca + alpha_tpa;
        // Amplitude transmission: field ∝ exp(−α·L/2), power ∝ exp(−α·L).
        let t_amp = (-alpha_total * self.length_m / 2.0).exp();

        // ── 7. Self-heating Δn (quasi-static) ────────────────────────────
        // Absorbed optical power: P_abs ≈ α_total · L · I (watts into heat).
        // Static ΔT = R_th · P_abs, yielding thermo-optic index change
        //   Δn_th = (dn/dT) · ΔT = dn_dt · R_th · P_abs
        let p_abs = alpha_total * self.length_m * intensity_w;
        let dn_self = self.dn_dt * self.r_th_k_per_w * p_abs;

        // ── 8. Per-channel optical phase (WDM-aware) ─────────────────────
        let two_pi = 2.0 * std::f64::consts::PI;
        let lam = wpc - 1; // index of the wavelength wire within each bundle
                           // Optional absolute-phase pinning: subtract φ_0 = 2π n_eff L / λ_ref
                           // so the device is "transparent" at λ = λ_ref (useful for testbench
                           // rings designed to be on-resonance at the laser wavelength).
        let phi_ref = if self.pin_at_ref {
            two_pi * self.n_eff * self.length_m / self.wl_ref_m
        } else {
            0.0
        };
        for k in 0..n {
            // Per-channel wavelength (falls back to λ_ref for undriven ports).
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
            // Group-index dispersion correction:
            //   n_eff(λ) ≈ n_eff_ref + (n_eff_ref − n_g) · (λ − λ_ref) / λ_ref
            let n_eff_lam = n_eff_at_lambda(self.n_eff, self.n_g, self.wl_ref_m, lambda);
            // Absolute propagation phase: φ_prop = 2π · n_eff(λ) · L / λ.
            let phi_abs = two_pi * n_eff_lam * self.length_m / lambda;
            let phi_prop = phi_abs - phi_ref;
            // Depletion EO phase: φ_eo_rev = 2π L · (dn_dv_rev · v_j) / λ.
            // dn_dv_rev is negative for Si PN (carrier depletion reduces n).
            let phi_eo_rev = two_pi * self.length_m * self.dn_dv_rev * v_junc / lambda;
            // Injection EO phase (Soref-Bennett):
            //   Δn_inj ∝ −K · (exp(v_j/V_T)−1)  (negative: more carriers → less n)
            let phi_eo_inj = -two_pi * self.length_m * self.dn_dv_inj * inj / lambda;
            // Thermo-optic phase from absorbed light (negligible at −10 dBm ring power).
            let phi_self = two_pi * self.length_m * dn_self / lambda;
            // Total phase: propagation + depletion EO + injection EO + self-heat.
            let phi = phi_prop + phi_eo_rev + phi_eo_inj + phi_self;
            // MNA optical stamp coefficients:
            //   out_re = t_amp · (cos φ · in_re − sin φ · in_im)   → c_cached
            //   out_im = t_amp · (sin φ · in_re + cos φ · in_im)  (uses s_cached)
            // We store t_amp·cos(φ) and t_amp·sin(φ) so load_jacobian can stamp
            // them directly into the potential-equation rows without extra work.
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

impl Default for NativePnThermalPhaseShifterFull {
    fn default() -> Self {
        Self::new()
    }
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
            terminals.len() >= stride + 4 && (terminals.len() - 4).is_multiple_of(stride),
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
