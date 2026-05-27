use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

/// Small floor conductance for numerical stability.
const GMIN: f64 = 1e-12;

/// SPICE Level 1 (Shichman-Hodges) MOSFET, with full capacitance model.
///
/// DC model parameters (`.model` card): VTO, KP, LAMBDA, GAMMA, PHI.
/// Gate capacitance (`.model`): CGSO, CGDO, CGBO (overlap, F/m), COX (channel, F/m²).
/// Junction capacitance (`.model`): CJ (F/m²), CJSW (F/m), PB, MJ, MJSW, FC.
/// Instance parameters (`M` card): W, L, AS, AD, PS, PD.
///
/// Gate capacitances use the Meyer model (region-dependent linear caps).
/// Junction capacitances use the same depletion-cap model as the diode.
pub struct Mosfet1 {
    // ── DC model parameters ──────────────────────────────────────────────────
    vto: f64,      // threshold voltage (V)
    kp: f64,       // process transconductance (A/V²)
    lambda: f64,   // channel-length modulation (1/V)
    gamma: f64,    // body-effect coefficient (V^0.5)
    phi: f64,      // surface potential (V)
    polarity: f64, // +1 for NMOS, −1 for PMOS

    // ── Gate cap model parameters (Meyer model) ──────────────────────────────
    cgso: f64, // gate-source overlap cap per channel width (F/m)
    cgdo: f64, // gate-drain  overlap cap per channel width (F/m)
    cgbo: f64, // gate-bulk   overlap cap per channel length (F/m)
    cox: f64,  // oxide cap density (F/m²); 0 if unspecified

    // ── Junction cap model parameters ───────────────────────────────────────
    cj: f64,   // zero-bias junction cap per unit area (F/m²)
    cjsw: f64, // zero-bias sidewall cap per unit perimeter (F/m)
    pb: f64,   // junction built-in potential (V)
    mj: f64,   // junction grading coefficient
    #[allow(dead_code)]
    mjsw: f64, // sidewall grading coefficient (default = mj)
    fc: f64,   // forward-bias depletion-cap linearisation boundary

    // ── Instance geometry ────────────────────────────────────────────────────
    w: f64, // channel width  (m)
    l: f64, // channel length (m)
    w_over_l: f64,
    as_: f64, // source diffusion area (m²)
    ad: f64,  // drain  diffusion area (m²)
    ps: f64,  // source diffusion perimeter (m)
    pd: f64,  // drain  diffusion perimeter (m)

    // ── Derived (set in setup_instance) ─────────────────────────────────────
    cgs_ov: f64, // CGSO * W
    cgd_ov: f64, // CGDO * W
    cgb_ov: f64, // CGBO * L
    cox_wl: f64, // COX  * W * L
    cbs0: f64,   // CJ * AS + CJSW * PS  (zero-bias bulk-source cap)
    cbd0: f64,   // CJ * AD + CJSW * PD  (zero-bias bulk-drain  cap)

    // ── Terminal bindings ────────────────────────────────────────────────────
    drain: NodeId,
    gate: NodeId,
    source: NodeId,
    bulk: NodeId,

    // ── DC Newton-Raphson state ──────────────────────────────────────────────
    gm: f64,
    gds: f64,
    gmbs: f64,
    jeq: f64,
    vgs_eff_prev: f64,
    vds_eff_prev: f64,

    // ── Gate cap transient state (Meyer model, linear caps) ──────────────────
    cgs_eval: f64, // Cgs at current NR iterate (set by eval when transient)
    cgd_eval: f64,
    cgb_eval: f64,
    q_gs_prev: f64, // Q_gs = Cgs * Vgs at previous accepted timestep
    q_gd_prev: f64,
    q_gb_prev: f64,

    // ── Junction cap transient state (nonlinear depletion cap) ───────────────
    vbs_eval: f64,  // Vbs = Vb−Vs at current iterate (for stamp use)
    vbd_eval: f64,  // Vbd = Vb−Vd
    cbs_eval: f64,  // Cbs(Vbs)
    cbd_eval: f64,  // Cbd(Vbd)
    q_bs_eval: f64, // Q_bs(Vbs) charge integral at current iterate
    q_bd_eval: f64,
    q_bs_prev: f64, // Q_bs at previous accepted timestep
    q_bd_prev: f64,
}

impl Mosfet1 {
    /// Construct from model-card parameters.
    pub fn from_model_params(is_pmos: bool, params: &[(String, f64)]) -> (Self, Vec<String>) {
        let mut vto = if is_pmos { -0.7 } else { 0.7 };
        let mut kp = 2e-5;
        let mut lambda = 0.0;
        let mut gamma = 0.0;
        let mut phi = 0.6;
        // Gate caps
        let mut cgso = 0.0_f64;
        let mut cgdo = 0.0_f64;
        let mut cgbo = 0.0_f64;
        let mut cox = 0.0_f64;
        // Junction caps
        let mut cj = 0.0_f64;
        let mut cjsw = 0.0_f64;
        let mut pb = 0.8_f64;
        let mut mj = 0.5_f64;
        let mut mjsw = 0.33_f64;
        let mut fc = 0.5_f64;
        let mut unknown = Vec::new();
        for (k, v) in params {
            match k.to_lowercase().as_str() {
                "vto" | "vth0" | "vtho" => vto = *v,
                "kp" => kp = *v,
                "lambda" => lambda = *v,
                "gamma" => gamma = *v,
                "phi" => phi = *v,
                "cgso" => cgso = *v,
                "cgdo" => cgdo = *v,
                "cgbo" => cgbo = *v,
                "cox" => cox = *v,
                "tox" => {
                    // COX = ε₀·ε_SiO2 / TOX; ε_r ≈ 3.9
                    const EPS_OX: f64 = 3.9 * 8.854187817e-12;
                    cox = EPS_OX / *v;
                }
                "cj" => cj = *v,
                "cjsw" => cjsw = *v,
                "pb" => pb = *v,
                "mj" => mj = *v,
                "mjsw" => mjsw = *v,
                "fc" => fc = *v,
                _ => unknown.push(k.clone()),
            }
        }
        let dev = Mosfet1 {
            vto,
            kp,
            lambda,
            gamma,
            phi,
            polarity: if is_pmos { -1.0 } else { 1.0 },
            cgso,
            cgdo,
            cgbo,
            cox,
            cj,
            cjsw,
            pb,
            mj,
            mjsw,
            fc,
            w: 1e-4,
            l: 1e-4,
            w_over_l: 1.0,
            as_: 0.0,
            ad: 0.0,
            ps: 0.0,
            pd: 0.0,
            cgs_ov: 0.0,
            cgd_ov: 0.0,
            cgb_ov: 0.0,
            cox_wl: 0.0,
            cbs0: 0.0,
            cbd0: 0.0,
            drain: None,
            gate: None,
            source: None,
            bulk: None,
            gm: GMIN,
            gds: GMIN,
            gmbs: 0.0,
            jeq: 0.0,
            vgs_eff_prev: 0.0,
            vds_eff_prev: 0.0,
            cgs_eval: 0.0,
            cgd_eval: 0.0,
            cgb_eval: 0.0,
            q_gs_prev: 0.0,
            q_gd_prev: 0.0,
            q_gb_prev: 0.0,
            vbs_eval: 0.0,
            vbd_eval: 0.0,
            cbs_eval: 0.0,
            cbd_eval: 0.0,
            q_bs_eval: 0.0,
            q_bd_eval: 0.0,
            q_bs_prev: 0.0,
            q_bd_prev: 0.0,
        };
        (dev, unknown)
    }

    /// Apply instance parameters (W, L, AS, AD, PS, PD).
    pub fn set_instance_params(&mut self, params: &[(String, f64)]) -> Vec<String> {
        let mut w = 1e-4_f64;
        let mut l = 1e-4_f64;
        let mut as_ = 0.0_f64;
        let mut ad = 0.0_f64;
        let mut ps = 0.0_f64;
        let mut pd = 0.0_f64;
        let mut unknown = Vec::new();
        for (k, v) in params {
            match k.to_lowercase().as_str() {
                "w" => w = *v,
                "l" => l = *v,
                "as" => as_ = *v,
                "ad" => ad = *v,
                "ps" => ps = *v,
                "pd" => pd = *v,
                _ => unknown.push(k.clone()),
            }
        }
        self.w = w;
        self.l = l;
        self.w_over_l = w / l;
        self.as_ = as_;
        self.ad = ad;
        self.ps = ps;
        self.pd = pd;
        // Pre-compute derived quantities.
        self.cgs_ov = self.cgso * w;
        self.cgd_ov = self.cgdo * w;
        self.cgb_ov = self.cgbo * l;
        self.cox_wl = self.cox * w * l;
        self.cbs0 = self.cj * as_ + self.cjsw * ps;
        self.cbd0 = self.cj * ad + self.cjsw * pd;
        unknown
    }

    // ── Depletion cap helpers (same model as ShockleyDiode::cj_depl / q_depl) ─

    fn cj_depl(&self, c0: f64, v: f64) -> f64 {
        if c0 == 0.0 {
            return 0.0;
        }
        let fc_pb = self.fc * self.pb;
        if v < fc_pb {
            c0 * (1.0 - v / self.pb).powf(-self.mj)
        } else {
            let k = (1.0 - self.fc).powf(1.0 + self.mj);
            c0 / k * (1.0 - self.fc * (1.0 + self.mj) + self.mj * v / self.pb)
        }
    }

    fn q_depl(&self, c0: f64, v: f64) -> f64 {
        if c0 == 0.0 {
            return 0.0;
        }
        let fc_pb = self.fc * self.pb;
        if v < fc_pb {
            let x = 1.0 - v / self.pb;
            c0 * self.pb / (1.0 - self.mj) * (1.0 - x.powf(1.0 - self.mj))
        } else {
            let x_fc = 1.0 - self.fc;
            let q_fc = c0 * self.pb / (1.0 - self.mj) * (1.0 - x_fc.powf(1.0 - self.mj));
            let k = x_fc.powf(1.0 + self.mj);
            let f2 = 1.0 - self.fc * (1.0 + self.mj);
            let dv = v - fc_pb;
            q_fc + c0 / k * (f2 * dv + self.mj / (2.0 * self.pb) * (v * v - fc_pb * fc_pb))
        }
    }

    // ── Stamp helpers ────────────────────────────────────────────────────────

    fn stamp_g_pair(mat: &mut MnaMatrix, a: NodeId, b: NodeId, g: f64) {
        if let Some(ai) = a {
            mat.a[ai][ai] += g;
            if let Some(bi) = b {
                mat.a[ai][bi] -= g;
            }
        }
        if let Some(bi) = b {
            mat.a[bi][bi] += g;
            if let Some(ai) = a {
                mat.a[bi][ai] -= g;
            }
        }
    }

    fn stamp_hist(b_vec: &mut [f64], a: NodeId, bnode: NodeId, i_hist: f64) {
        if let Some(ai) = a {
            b_vec[ai] += i_hist;
        }
        if let Some(bi) = bnode {
            b_vec[bi] -= i_hist;
        }
    }
}

impl Device for Mosfet1 {
    fn num_terminals(&self) -> usize {
        4
    }

    fn setup_model(&mut self, _ctx: &SimContext) {}

    fn setup_instance(&mut self, terminals: &[NodeId], _ctx: &SimContext) {
        debug_assert_eq!(terminals.len(), 4, "MOSFET expects [D, G, S, B]");
        self.drain = terminals[0];
        self.gate = terminals[1];
        self.source = terminals[2];
        self.bulk = terminals[3];
    }

    fn eval(&mut self, x: &[f64], flags: EvalFlags, ctx: &SimContext) {
        let pol = self.polarity;
        let vd = self.drain.map_or(0.0, |i| x[i]);
        let vg = self.gate.map_or(0.0, |i| x[i]);
        let vs = self.source.map_or(0.0, |i| x[i]);
        let vb = self.bulk.map_or(0.0, |i| x[i]);

        // Polarity-flipped voltages (PMOS sees inverted potential differences).
        let mut vgs_eff = pol * (vg - vs);
        let mut vds_eff = pol * (vd - vs);
        let vbs_eff = pol * (vb - vs);

        // fetlim: limit Vgs steps above vto (channel exponential blow-up
        // doesn't happen at L1 but the limiter helps when MOSFETs are stacked
        // with diodes/BJTs).  Also clamp Vds sign-changes to keep NR out of
        // the triode/saturation ping-pong basin.
        if ctx.jlim_enabled {
            let vto_eff = pol * self.vto;
            let dvg = vgs_eff - self.vgs_eff_prev;
            if vgs_eff > vto_eff && dvg.abs() > 1.0 {
                vgs_eff = self.vgs_eff_prev + dvg.signum() * (1.0 + (dvg.abs() - 1.0).ln_1p());
            }
            if self.vds_eff_prev.abs() > 1e-6 && self.vds_eff_prev * vds_eff < 0.0 {
                vds_eff = 0.1 * self.vds_eff_prev;
            }
        }
        self.vgs_eff_prev = vgs_eff;
        self.vds_eff_prev = vds_eff;

        // Threshold voltage with body effect.
        let phi_m_vbs = (self.phi - vbs_eff).max(1e-10);
        let vto_eff = pol * self.vto;
        let vth = vto_eff + self.gamma * (phi_m_vbs.sqrt() - self.phi.sqrt());

        let (ids_eff, gm_eff, gds_eff, gmbs_eff) = if vgs_eff < vth {
            (0.0, 0.0, 0.0, 0.0)
        } else {
            let vdsat = vgs_eff - vth;
            let beta = self.kp * self.w_over_l;
            let clm = 1.0 + self.lambda * vds_eff;
            let dvth_dvbs = if self.gamma > 0.0 {
                -self.gamma / (2.0 * phi_m_vbs.sqrt())
            } else {
                0.0
            };

            if vds_eff < vdsat {
                let ids = beta * ((vgs_eff - vth) * vds_eff - 0.5 * vds_eff * vds_eff) * clm;
                let gm = beta * vds_eff * clm;
                let gds = beta * (vdsat - vds_eff) * clm
                    + beta * ((vgs_eff - vth) * vds_eff - 0.5 * vds_eff * vds_eff) * self.lambda;
                let gmbs = -gm * dvth_dvbs;
                (ids, gm, gds, gmbs)
            } else {
                let ids = 0.5 * beta * vdsat * vdsat * clm;
                let gm = beta * vdsat * clm;
                let gds = 0.5 * beta * vdsat * vdsat * self.lambda;
                let gmbs = -gm * dvth_dvbs;
                (ids, gm, gds, gmbs)
            }
        };

        let ids_real = pol * ids_eff;
        let vgs = pol * vgs_eff;
        let vds = pol * vds_eff;
        let vbs = vb - vs;
        let gds_total = gds_eff + GMIN;
        self.gm = gm_eff;
        self.gds = gds_total;
        self.gmbs = gmbs_eff;
        self.jeq = ids_real - gm_eff * vgs - gds_total * vds - gmbs_eff * vbs;

        // ── Capacitance evaluation ──────────────────────────────────────────
        if flags.transient {
            // Gate caps: Meyer model (region-dependent linear caps).
            let vth_capped = vth; // already computed above
            let vdsat = (vgs_eff - vth_capped).max(0.0);
            let (cgs_ch, cgd_ch, cgb_ch) = if vgs_eff < vth_capped {
                // Cutoff: all channel charge to gate-bulk.
                (0.0, 0.0, self.cox_wl)
            } else if vds_eff < vdsat {
                // Triode: split equally between Cgs and Cgd.
                (0.5 * self.cox_wl, 0.5 * self.cox_wl, 0.0)
            } else {
                // Saturation: 2/3 to Cgs, none to Cgd.
                (2.0 / 3.0 * self.cox_wl, 0.0, 0.0)
            };
            self.cgs_eval = self.cgs_ov + cgs_ch;
            self.cgd_eval = self.cgd_ov + cgd_ch;
            self.cgb_eval = self.cgb_ov + cgb_ch;

            // Junction caps: depletion model.
            // Use real (non-pol-flipped) junction voltages.
            // V < 0 → reverse-biased (normal operation).
            let vbs_j = vb - vs;
            let vbd_j = vb - vd;
            self.vbs_eval = vbs_j;
            self.vbd_eval = vbd_j;
            self.cbs_eval = self.cj_depl(self.cbs0, vbs_j);
            self.cbd_eval = self.cj_depl(self.cbd0, vbd_j);
            self.q_bs_eval = self.q_depl(self.cbs0, vbs_j);
            self.q_bd_eval = self.q_depl(self.cbd0, vbd_j);
        }
    }

    fn load_residual(&self, b: &mut [f64]) {
        if let Some(d) = self.drain {
            b[d] -= self.jeq;
        }
        if let Some(s) = self.source {
            b[s] += self.jeq;
        }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let (d, g, s, bk) = (self.drain, self.gate, self.source, self.bulk);
        let gm = self.gm;
        let gds = self.gds;
        let gmbs = self.gmbs;
        let gms = gm + gds + gmbs;

        macro_rules! stamp {
            ($ri:expr, $ci:expr, $val:expr) => {
                if let (Some(r), Some(c)) = ($ri, $ci) {
                    mat.a[r][c] += $val;
                }
            };
        }

        stamp!(d, g, gm);
        stamp!(d, d, gds);
        stamp!(d, s, -gms);
        stamp!(d, bk, gmbs);

        stamp!(s, g, -gm);
        stamp!(s, d, -gds);
        stamp!(s, s, gms);
        stamp!(s, bk, -gmbs);
    }

    fn load_residual_tran(&self, b: &mut [f64], alpha: f64) {
        self.load_residual(b);
        let (d, g, s, bk) = (self.drain, self.gate, self.source, self.bulk);

        // ── Gate caps (linear): i_hist = alpha * Q_prev ─────────────────────
        // Gate-Source
        if self.cgs_eval != 0.0 || self.q_gs_prev != 0.0 {
            Self::stamp_hist(b, g, s, alpha * self.q_gs_prev);
        }
        // Gate-Drain
        if self.cgd_eval != 0.0 || self.q_gd_prev != 0.0 {
            Self::stamp_hist(b, g, d, alpha * self.q_gd_prev);
        }
        // Gate-Bulk
        if self.cgb_eval != 0.0 || self.q_gb_prev != 0.0 {
            Self::stamp_hist(b, g, bk, alpha * self.q_gb_prev);
        }

        // ── Junction caps (nonlinear): i_hist = alpha*(C*V + Q_prev − Q_now) ─
        if self.cbs0 != 0.0 {
            let i_hist = alpha * (self.cbs_eval * self.vbs_eval + self.q_bs_prev - self.q_bs_eval);
            Self::stamp_hist(b, bk, s, i_hist);
        }
        if self.cbd0 != 0.0 {
            let i_hist = alpha * (self.cbd_eval * self.vbd_eval + self.q_bd_prev - self.q_bd_eval);
            Self::stamp_hist(b, bk, d, i_hist);
        }
    }

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, alpha: f64) {
        self.load_jacobian(mat);
        let (d, g, s, bk) = (self.drain, self.gate, self.source, self.bulk);

        // Gate caps.
        if self.cgs_eval != 0.0 {
            Self::stamp_g_pair(mat, g, s, alpha * self.cgs_eval);
        }
        if self.cgd_eval != 0.0 {
            Self::stamp_g_pair(mat, g, d, alpha * self.cgd_eval);
        }
        if self.cgb_eval != 0.0 {
            Self::stamp_g_pair(mat, g, bk, alpha * self.cgb_eval);
        }

        // Junction caps.
        if self.cbs_eval != 0.0 {
            Self::stamp_g_pair(mat, bk, s, alpha * self.cbs_eval);
        }
        if self.cbd_eval != 0.0 {
            Self::stamp_g_pair(mat, bk, d, alpha * self.cbd_eval);
        }
    }

    fn commit_timestep(&mut self, x: &[f64]) {
        let vd = self.drain.map_or(0.0, |i| x[i]);
        let vg = self.gate.map_or(0.0, |i| x[i]);
        let vs = self.source.map_or(0.0, |i| x[i]);
        let vb = self.bulk.map_or(0.0, |i| x[i]);

        // Save gate cap charges (linear caps: Q = C * V).
        self.q_gs_prev = self.cgs_eval * (vg - vs);
        self.q_gd_prev = self.cgd_eval * (vg - vd);
        self.q_gb_prev = self.cgb_eval * (vg - vb);

        // Save junction cap charges (nonlinear depletion cap).
        let vbs_j = vb - vs;
        let vbd_j = vb - vd;
        self.q_bs_prev = self.q_depl(self.cbs0, vbs_j);
        self.q_bd_prev = self.q_depl(self.cbd0, vbd_j);

        // Update fetlim reference.
        self.vgs_eff_prev = self.polarity * (vg - vs);
        self.vds_eff_prev = self.polarity * (vd - vs);
    }

    fn noise_sources(&self, ctx: &SimContext) -> Vec<(NodeId, NodeId, f64)> {
        if self.gm.abs() < 1e-18 {
            return Vec::new();
        }
        const K_BOLTZMANN: f64 = 1.380649e-23;
        let s_i = 8.0 / 3.0 * K_BOLTZMANN * ctx.temperature * self.gm.abs();
        vec![(self.drain, self.source, s_i)]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::device::EvalFlags;

    fn ctx() -> SimContext {
        SimContext::default()
    }

    fn nmos(vto: f64, kp: f64, w_over_l: f64) -> Mosfet1 {
        let (mut m, _) =
            Mosfet1::from_model_params(false, &[("vto".into(), vto), ("kp".into(), kp)]);
        m.w_over_l = w_over_l;
        m.setup_model(&ctx());
        m.setup_instance(&[Some(0), Some(1), Some(2), None], &ctx());
        m
    }

    #[test]
    fn nmos_cutoff_no_current() {
        let mut m = nmos(1.0, 100e-6, 10.0);
        let x = [0.0_f64, 0.0, 0.0];
        m.eval(&x, EvalFlags::dc(), &ctx());
        assert!(m.gm.abs() < 1e-15, "gm in cutoff: {}", m.gm);
        assert!((m.gds - GMIN).abs() < 1e-20, "gds in cutoff: {}", m.gds);
        assert!(
            m.jeq.abs() < 1e-20,
            "jeq in cutoff with VDS=0: {:.3e}",
            m.jeq
        );
    }

    #[test]
    fn nmos_saturation_ids() {
        let kp = 100e-6;
        let mut m = nmos(1.0, kp, 10.0);
        m.vgs_eff_prev = 2.0;
        m.vds_eff_prev = 3.0;
        let x = [3.0_f64, 2.0, 0.0];
        m.eval(&x, EvalFlags::dc(), &ctx());

        let ids_expected = 0.5 * kp * 10.0 * 1.0 * 1.0;
        let ids_from_stamp = m.jeq + m.gm * 2.0 + m.gds * 3.0;
        assert!(
            (ids_from_stamp - ids_expected).abs() < 1e-10,
            "IDS sat: {:.4e} expected {:.4e}",
            ids_from_stamp,
            ids_expected
        );
    }

    #[test]
    fn nmos_triode_ids() {
        let kp = 100e-6;
        let mut m = nmos(1.0, kp, 10.0);
        let x = [0.5_f64, 2.0, 0.0];
        m.eval(&x, EvalFlags::dc(), &ctx());

        let ids_expected = kp * 10.0 * (1.0 * 0.5 - 0.5 * 0.25);
        let ids_from_stamp = m.jeq + m.gm * 2.0 + m.gds * 0.5;
        assert!(
            (ids_from_stamp - ids_expected).abs() < 1e-10,
            "IDS triode: {:.4e} expected {:.4e}",
            ids_from_stamp,
            ids_expected
        );
    }

    #[test]
    fn pmos_saturation_ids() {
        let kp = 100e-6;
        let (mut m, _) =
            Mosfet1::from_model_params(true, &[("vto".into(), -1.0), ("kp".into(), kp)]);
        m.w_over_l = 10.0;
        m.setup_model(&ctx());
        m.setup_instance(&[Some(0), Some(1), Some(2), None], &ctx());

        m.vgs_eff_prev = 2.0;
        m.vds_eff_prev = 2.0;
        let x = [0.0_f64, 0.0, 2.0];
        m.eval(&x, EvalFlags::dc(), &ctx());

        let ids_expected = -500e-6_f64;
        let vgs = 0.0 - 2.0;
        let vds = 0.0 - 2.0;
        let ids_from_stamp = m.jeq + m.gm * vgs + m.gds * vds;
        assert!(
            (ids_from_stamp - ids_expected).abs() < 1e-9,
            "PMOS IDS: {:.4e} expected {:.4e}",
            ids_from_stamp,
            ids_expected
        );
    }

    #[test]
    fn jacobian_matches_numerical_derivative() {
        let kp = 100e-6;
        let eps = 1e-6;
        let (mut m, _) = Mosfet1::from_model_params(
            false,
            &[
                ("vto".into(), 1.0),
                ("kp".into(), kp),
                ("lambda".into(), 0.05),
            ],
        );
        m.w_over_l = 10.0;
        m.setup_model(&ctx());
        m.setup_instance(&[Some(0), Some(1), Some(2), None], &ctx());
        let x0 = [3.0_f64, 2.0, 0.0];

        m.vgs_eff_prev = 2.0;
        m.vds_eff_prev = 3.0;
        m.eval(&x0, EvalFlags::dc(), &ctx());
        let gm_analytic = m.gm;
        let gds_analytic = m.gds;

        let ids0 = m.jeq + gm_analytic * (x0[1] - x0[2]) + gds_analytic * (x0[0] - x0[2]);

        let mut xg = x0;
        xg[1] += eps;
        m.eval(&xg, EvalFlags::dc(), &ctx());
        let ids_g = m.jeq + m.gm * (xg[1] - xg[2]) + m.gds * (xg[0] - xg[2]);
        let gm_fd = (ids_g - ids0) / eps;

        let mut xd = x0;
        xd[0] += eps;
        m.eval(&xd, EvalFlags::dc(), &ctx());
        let ids_d = m.jeq + m.gm * (xd[1] - xd[2]) + m.gds * (xd[0] - xd[2]);
        let gds_fd = (ids_d - ids0) / eps;

        m.eval(&x0, EvalFlags::dc(), &ctx());
        assert!(
            (gm_analytic - gm_fd).abs() / gm_analytic.abs() < 1e-4,
            "gm analytic={:.4e} fd={:.4e}",
            gm_analytic,
            gm_fd
        );
        assert!(
            (gds_analytic - gds_fd).abs().max(1e-14) / (gds_analytic.abs() + GMIN) < 0.01,
            "gds analytic={:.4e} fd={:.4e}",
            gds_analytic,
            gds_fd
        );
    }

    #[test]
    fn cgs_stamps_in_tran() {
        // NMOS with Cgs overlap only. At VGS=2, VDS=3 (saturation):
        // Meyer: Cgs = CGSO*W + (2/3)*COX*W*L
        let (mut m, _) = Mosfet1::from_model_params(
            false,
            &[
                ("vto".into(), 1.0),
                ("kp".into(), 100e-6),
                ("cgso".into(), 1e-9), // 1nF/m overlap
                ("cox".into(), 1e-3),  // 1mF/m²
            ],
        );
        let _ = m.set_instance_params(&[
            ("w".into(), 1e-4), // 100 µm
            ("l".into(), 1e-4), // 100 µm
        ]);
        m.setup_model(&ctx());
        m.setup_instance(&[Some(0), Some(1), Some(2), None], &ctx());
        m.vgs_eff_prev = 2.0;
        m.vds_eff_prev = 3.0;
        let x = [3.0_f64, 2.0, 0.0];
        m.eval(&x, EvalFlags::tran(), &ctx());

        // cgs_ov = 1e-9 * 1e-4 = 1e-13 F
        // cox_wl = 1e-3 * 1e-4 * 1e-4 = 1e-11 F; sat: cgs_ch = 2/3 * 1e-11
        let cgs_expected = 1e-9 * 1e-4 + (2.0 / 3.0) * 1e-3 * 1e-4 * 1e-4;
        assert!(
            (m.cgs_eval - cgs_expected).abs() / cgs_expected < 1e-9,
            "Cgs in sat: {:.4e} expected {:.4e}",
            m.cgs_eval,
            cgs_expected
        );
        // In cutoff (VGS=0 < VTO=1V): channel cap goes to Cgb
        let x_off = [0.0_f64, 0.0, 0.0];
        m.eval(&x_off, EvalFlags::tran(), &ctx());
        assert!(
            (m.cgs_eval - m.cgs_ov).abs() < 1e-25,
            "Cgs should be overlap-only in cutoff: got {:.4e}, expected {:.4e}",
            m.cgs_eval,
            m.cgs_ov
        );
        assert!(
            (m.cgb_eval - (m.cgb_ov + m.cox_wl)).abs() < 1e-25,
            "Cgb in cutoff should include channel cap"
        );
    }

    #[test]
    fn junction_caps_nonzero_when_params_set() {
        let (mut m, _) = Mosfet1::from_model_params(
            false,
            &[
                ("vto".into(), 1.0),
                ("kp".into(), 100e-6),
                ("cj".into(), 0.5e-3), // 0.5 mF/m²
                ("cjsw".into(), 1e-9), // 1 nF/m
                ("pb".into(), 0.8),
                ("mj".into(), 0.5),
            ],
        );
        let _ = m.set_instance_params(&[
            ("w".into(), 10e-6),
            ("l".into(), 0.35e-6),
            ("as".into(), 50e-12), // 50 µm²
            ("ad".into(), 50e-12),
            ("ps".into(), 20e-6), // 20 µm perimeter
            ("pd".into(), 20e-6),
        ]);
        m.setup_model(&ctx());
        m.setup_instance(&[Some(0), Some(1), Some(2), None], &ctx());
        // At Vbs = -1V (reverse biased bulk-source junction), cap should be < cbs0
        let _x = [3.0_f64, 2.0, 0.0, -1.0];
        // But wait, bulk is None (gnd) and we have 3 nodes → need to drop to 3-node test
        // Use 3 nodes with bulk grounded.
        let (mut m2, _) = Mosfet1::from_model_params(
            false,
            &[
                ("vto".into(), 1.0),
                ("kp".into(), 100e-6),
                ("cj".into(), 0.5e-3),
                ("cjsw".into(), 1e-9),
                ("pb".into(), 0.8),
                ("mj".into(), 0.5),
            ],
        );
        let _ = m2.set_instance_params(&[
            ("w".into(), 10e-6),
            ("l".into(), 0.35e-6),
            ("as".into(), 50e-12),
            ("ad".into(), 50e-12),
            ("ps".into(), 20e-6),
            ("pd".into(), 20e-6),
        ]);
        m2.setup_model(&ctx());
        m2.setup_instance(&[Some(0), Some(1), Some(2), None], &ctx());
        // VS=0, VB=gnd=0 → Vbs=0: Cbs = cbs0
        let cbs0_expected = 0.5e-3 * 50e-12 + 1e-9 * 20e-6;
        let x3 = [3.0_f64, 2.0, 0.0];
        m2.vgs_eff_prev = 2.0;
        m2.vds_eff_prev = 3.0;
        m2.eval(&x3, EvalFlags::tran(), &ctx());
        assert!(
            (m2.cbs_eval - cbs0_expected).abs() / cbs0_expected < 1e-9,
            "Cbs at Vbs=0 should equal cbs0: {:.4e} vs {:.4e}",
            m2.cbs_eval,
            cbs0_expected
        );
        assert!(m2.cbs_eval > 0.0, "Cbs should be positive");
    }
}
