use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

/// Small floor conductance for numerical stability.
const GMIN: f64 = 1e-12;

/// SPICE Level 1 (Shichman-Hodges) MOSFET.
///
/// Model parameters (from `.model` card): VTO, KP, LAMBDA, GAMMA, PHI.
/// Instance parameters (from `M` card): W, L.
/// Supports NMOS and PMOS via an internal polarity flag (pol = +1 / −1).
///
/// Physics summary for NMOS (VGS_eff = pol·(VG−VS), etc.):
///   Cutoff  (VGS_eff < Vth): IDS = 0
///   Triode  (VDS_eff < VGS_eff − Vth):
///     IDS = KP·(W/L)·[(VGS_eff−Vth)·VDS_eff − 0.5·VDS_eff²]·(1+λ·VDS_eff)
///   Saturation:
///     IDS = 0.5·KP·(W/L)·(VGS_eff−Vth)²·(1+λ·VDS_eff)
///
/// Real current: pol·IDS_eff (negative for PMOS in active region).
pub struct Mosfet1 {
    // Model params
    vto: f64,       // threshold voltage (V) as given in model card (signed for PMOS)
    kp: f64,        // process transconductance (A/V²)
    lambda: f64,    // channel-length modulation (1/V)
    gamma: f64,     // body-effect coefficient (V^0.5)
    phi: f64,       // surface potential (V)
    polarity: f64,  // +1 for NMOS, −1 for PMOS

    // Instance computed
    w_over_l: f64,

    // Terminal bindings
    drain:  NodeId,
    gate:   NodeId,
    source: NodeId,
    bulk:   NodeId,

    // Cached Newton-Raphson equivalent circuit (set by eval)
    gm:   f64,
    gds:  f64,
    gmbs: f64,
    jeq:  f64, // IDS_real − gm·VGS − gds·VDS − gmbs·VBS

    // fetlim state: previous-iter Vgs_eff and Vds_eff for step-limiting.
    vgs_eff_prev: f64,
    vds_eff_prev: f64,
}

impl Mosfet1 {
    /// Construct from model-card parameters. Returns the device and any unrecognised param names.
    /// `is_pmos`: true for a PMOS model card (`PMOS`).
    pub fn from_model_params(is_pmos: bool, params: &[(String, f64)]) -> (Self, Vec<String>) {
        let mut vto    = if is_pmos { -0.7 } else { 0.7 };
        let mut kp     = 2e-5;
        let mut lambda = 0.0;
        let mut gamma  = 0.0;
        let mut phi    = 0.6;
        let mut unknown = Vec::new();
        for (k, v) in params {
            match k.to_lowercase().as_str() {
                "vto" | "vth0" | "vtho" => vto    = *v,
                "kp"                    => kp     = *v,
                "lambda"                => lambda = *v,
                "gamma"                 => gamma  = *v,
                "phi"                   => phi    = *v,
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
            w_over_l: 1.0,
            drain:  None,
            gate:   None,
            source: None,
            bulk:   None,
            gm:   GMIN,
            gds:  GMIN,
            gmbs: 0.0,
            jeq:  0.0,
            vgs_eff_prev: 0.0,
            vds_eff_prev: 0.0,
        };
        (dev, unknown)
    }

    /// Apply instance parameters (W, L). Returns any unrecognised param names.
    pub fn set_instance_params(&mut self, params: &[(String, f64)]) -> Vec<String> {
        let mut w = 1e-4;
        let mut l = 1e-4;
        let mut unknown = Vec::new();
        for (k, v) in params {
            match k.to_lowercase().as_str() {
                "w" => w = *v,
                "l" => l = *v,
                _ => unknown.push(k.clone()),
            }
        }
        self.w_over_l = w / l;
        unknown
    }
}

impl Device for Mosfet1 {
    fn num_terminals(&self) -> usize { 4 }

    fn setup_model(&mut self, _ctx: &SimContext) {
        // Level 1 has no temperature-scaled parameters in this implementation.
    }

    fn setup_instance(&mut self, terminals: &[NodeId], _ctx: &SimContext) {
        debug_assert_eq!(terminals.len(), 4, "MOSFET expects [D, G, S, B]");
        self.drain  = terminals[0];
        self.gate   = terminals[1];
        self.source = terminals[2];
        self.bulk   = terminals[3];
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, ctx: &SimContext) {
        let pol = self.polarity;
        let vd = self.drain .map_or(0.0, |i| x[i]);
        let vg = self.gate  .map_or(0.0, |i| x[i]);
        let vs = self.source.map_or(0.0, |i| x[i]);
        let vb = self.bulk  .map_or(0.0, |i| x[i]);

        // Polarity-flipped voltages (PMOS sees inverted potential differences).
        let mut vgs_eff = pol * (vg - vs);
        let mut vds_eff = pol * (vd - vs);
        let vbs_eff = pol * (vb - vs);

        // fetlim: limit Vgs steps above vto (channel exponential blow-up
        // doesn't happen at L1 but the limiter helps when MOSFETs are stacked
        // with diodes/BJTs).  Also clamp Vds sign-changes to keep NR out of
        // the triode/saturation ping-pong basin.  The state is reset on every
        // accepted timestep via vgs_eff_prev = vgs_eff at convergence.
        if ctx.jlim_enabled {
            let vto_eff = pol * self.vto;
            let dvg = vgs_eff - self.vgs_eff_prev;
            if vgs_eff > vto_eff && dvg.abs() > 1.0 {
                // Log-compress steps above threshold (mirrors pnjlim).
                vgs_eff = self.vgs_eff_prev
                    + dvg.signum() * (1.0 + (dvg.abs() - 1.0).ln_1p());
            }
            // Vds sign-change clamp.
            if self.vds_eff_prev.abs() > 1e-6
                && self.vds_eff_prev * vds_eff < 0.0
            {
                vds_eff = 0.1 * self.vds_eff_prev;
            }
        }
        self.vgs_eff_prev = vgs_eff;
        self.vds_eff_prev = vds_eff;

        // Threshold voltage with body effect; clamp phi−VBS to avoid sqrt(negative).
        let phi_m_vbs = (self.phi - vbs_eff).max(1e-10);
        let vto_eff   = pol * self.vto; // always positive after sign-flip
        let vth = vto_eff + self.gamma * (phi_m_vbs.sqrt() - self.phi.sqrt());

        let (ids_eff, gm_eff, gds_eff, gmbs_eff) = if vgs_eff < vth {
            // Cutoff: all currents and conductances are zero.
            (0.0, 0.0, 0.0, 0.0)
        } else {
            let vdsat  = vgs_eff - vth;
            let beta   = self.kp * self.w_over_l;
            let clm    = 1.0 + self.lambda * vds_eff;
            // gmbs: dIDS_eff / d(vbs_eff) = −gm · (dVth/dvbs_eff)
            // dVth/dvbs_eff = −gamma/(2·sqrt(phi−vbs_eff))
            let dvth_dvbs = if self.gamma > 0.0 {
                -self.gamma / (2.0 * phi_m_vbs.sqrt())
            } else {
                0.0
            };

            if vds_eff < vdsat {
                // Triode (linear) region.
                let ids = beta * ((vgs_eff - vth) * vds_eff - 0.5 * vds_eff * vds_eff) * clm;
                let gm  = beta * vds_eff * clm;
                let gds = beta * (vdsat - vds_eff) * clm
                        + beta * ((vgs_eff - vth) * vds_eff - 0.5 * vds_eff * vds_eff)
                          * self.lambda;
                let gmbs = -gm * dvth_dvbs;
                (ids, gm, gds, gmbs)
            } else {
                // Saturation region.
                let ids = 0.5 * beta * vdsat * vdsat * clm;
                let gm  = beta * vdsat * clm;
                let gds = 0.5 * beta * vdsat * vdsat * self.lambda;
                let gmbs = -gm * dvth_dvbs;
                (ids, gm, gds, gmbs)
            }
        };

        // Real current (pol-corrected): positive IDS flows D→S for NMOS,
        // and negative IDS (S→D) for PMOS in active region.
        let ids_real = pol * ids_eff;

        // Norton current offset (in real-voltage space; conductances are pol^2 = 1 invariant).
        // After fetlim, vgs_eff / vds_eff may differ from the raw terminal differences;
        // the linearization point must use the limited values for the Norton offset
        // to be self-consistent with the cached ids_real.
        let vgs = pol * vgs_eff;
        let vds = pol * vds_eff;
        let vbs = vb - vs;
        let gds_total = gds_eff + GMIN;
        self.gm   = gm_eff;
        self.gds  = gds_total;
        self.gmbs = gmbs_eff;
        self.jeq  = ids_real - gm_eff * vgs - gds_total * vds - gmbs_eff * vbs;
    }

    fn load_residual(&self, b: &mut [f64]) {
        if let Some(d) = self.drain  { b[d] -= self.jeq; }
        if let Some(s) = self.source { b[s] += self.jeq; }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let (d, g, s, bk) = (self.drain, self.gate, self.source, self.bulk);
        let gm   = self.gm;
        let gds  = self.gds;
        let gmbs = self.gmbs;
        let gms  = gm + gds + gmbs; // dIDS/dVS = −(gm + gds + gmbs)

        macro_rules! stamp {
            ($ri:expr, $ci:expr, $val:expr) => {
                if let (Some(r), Some(c)) = ($ri, $ci) {
                    mat.a[r][c] += $val;
                }
            };
        }

        // IDS flows D→S; stamp rows for drain (−IDS) and source (+IDS).
        stamp!(d, g,  gm);
        stamp!(d, d,  gds);
        stamp!(d, s, -gms);
        stamp!(d, bk, gmbs);

        stamp!(s, g,  -gm);
        stamp!(s, d,  -gds);
        stamp!(s, s,   gms);
        stamp!(s, bk, -gmbs);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::device::EvalFlags;

    fn ctx() -> SimContext { SimContext::default() }

    fn nmos(vto: f64, kp: f64, w_over_l: f64) -> Mosfet1 {
        let (mut m, _) = Mosfet1::from_model_params(false, &[
            ("vto".into(), vto),
            ("kp".into(), kp),
        ]);
        m.w_over_l = w_over_l;
        m.setup_model(&ctx());
        m.setup_instance(&[Some(0), Some(1), Some(2), None], &ctx());
        m
    }

    #[test]
    fn nmos_cutoff_no_current() {
        // VGS = 0 < VTO = 1V → cutoff.
        let mut m = nmos(1.0, 100e-6, 10.0);
        let x = [0.0_f64, 0.0, 0.0]; // VD=0, VG=0, VS=0 (VDS=0 so GMIN term vanishes)
        m.eval(&x, EvalFlags::dc(), &ctx());
        // In cutoff with VDS=0: gm=0, gds=GMIN, jeq=0.
        assert!(m.gm.abs() < 1e-15, "gm in cutoff: {}", m.gm);
        assert!((m.gds - GMIN).abs() < 1e-20, "gds in cutoff: {}", m.gds);
        assert!(m.jeq.abs() < 1e-20, "jeq in cutoff with VDS=0: {:.3e}", m.jeq);
    }

    #[test]
    fn nmos_saturation_ids() {
        // VGS=2V, VTO=1V, VDS=3V → saturation (VDS > VGS-VTO=1V).
        // IDS = 0.5 * KP * (W/L) * (VGS-VTO)^2 = 0.5 * 100e-6 * 10 * 1^2 = 500µA.
        let kp = 100e-6;
        let mut m = nmos(1.0, kp, 10.0);
        // Pre-seed fetlim prev to the operating point so the single-eval test
        // doesn't see step-limiting (mirrors the diode pnjlim test pattern).
        m.vgs_eff_prev = 2.0;
        m.vds_eff_prev = 3.0;
        let x = [3.0_f64, 2.0, 0.0]; // VD=3, VG=2, VS=0
        m.eval(&x, EvalFlags::dc(), &ctx());

        let ids_expected = 0.5 * kp * 10.0 * 1.0 * 1.0;
        // jeq = IDS - gm*VGS - gds*VDS; recover IDS from stamp
        let ids_from_stamp = m.jeq + m.gm * 2.0 + m.gds * 3.0;
        assert!(
            (ids_from_stamp - ids_expected).abs() < 1e-10,
            "IDS sat: {:.4e} expected {:.4e}", ids_from_stamp, ids_expected
        );
    }

    #[test]
    fn nmos_triode_ids() {
        // VGS=2V, VTO=1V, VDS=0.5V → triode (VDS < VGS-VTO=1V).
        // IDS = KP*(W/L)*[(VGS-VTO)*VDS - 0.5*VDS²]
        //     = 100e-6*10*[1*0.5 - 0.5*0.25] = 1e-3 * [0.5 - 0.125] = 375µA.
        let kp = 100e-6;
        let mut m = nmos(1.0, kp, 10.0);
        let x = [0.5_f64, 2.0, 0.0]; // VD=0.5, VG=2, VS=0
        m.eval(&x, EvalFlags::dc(), &ctx());

        let ids_expected = kp * 10.0 * (1.0 * 0.5 - 0.5 * 0.25);
        let ids_from_stamp = m.jeq + m.gm * 2.0 + m.gds * 0.5;
        assert!(
            (ids_from_stamp - ids_expected).abs() < 1e-10,
            "IDS triode: {:.4e} expected {:.4e}", ids_from_stamp, ids_expected
        );
    }

    #[test]
    fn pmos_saturation_ids() {
        // PMOS: VGS=-2V (VG=0,VS=2), VTO=-1V, VDS=-3V (VD=0, VS=3? let's use: VD=0,VG=0,VS=2).
        // pol=-1: VGS_eff = -(0-2) = 2V, VTO_eff = -(-1) = 1V → on.
        // VDS_eff = -(0-2) = 2 > VGS_eff-VTO_eff = 1 → saturation.
        // IDS_eff = 0.5*KP*(W/L)*(2-1)^2 = 0.5*100e-6*10*1 = 500µA.
        // Real IDS = -1 * 500µA = -500µA (current flows S→D, from VS=2 to VD=0).
        let kp = 100e-6;
        let (mut m, _) = Mosfet1::from_model_params(true, &[
            ("vto".into(), -1.0),
            ("kp".into(), kp),
        ]);
        m.w_over_l = 10.0;
        m.setup_model(&ctx());
        // D=node0, G=node1, S=node2, B=gnd
        m.setup_instance(&[Some(0), Some(1), Some(2), None], &ctx());

        // Pre-seed fetlim prev to bypass single-eval step limiting.
        m.vgs_eff_prev = 2.0;
        m.vds_eff_prev = 2.0;
        let x = [0.0_f64, 0.0, 2.0]; // VD=0, VG=0, VS=2
        m.eval(&x, EvalFlags::dc(), &ctx());

        let ids_expected = -500e-6_f64; // negative: flows S→D
        let vgs = 0.0 - 2.0; // -2V
        let vds = 0.0 - 2.0; // -2V
        let ids_from_stamp = m.jeq + m.gm * vgs + m.gds * vds;
        assert!(
            (ids_from_stamp - ids_expected).abs() < 1e-9,
            "PMOS IDS: {:.4e} expected {:.4e}", ids_from_stamp, ids_expected
        );
    }

    #[test]
    fn jacobian_matches_numerical_derivative() {
        // Finite-difference check for gm and gds in saturation.
        // Use lambda=0.05 so gds is real (not just GMIN), making the FD comparison meaningful.
        let kp = 100e-6;
        let eps = 1e-6;
        let (mut m, _) = Mosfet1::from_model_params(false, &[
            ("vto".into(), 1.0),
            ("kp".into(), kp),
            ("lambda".into(), 0.05),
        ]);
        m.w_over_l = 10.0;
        m.setup_model(&ctx());
        m.setup_instance(&[Some(0), Some(1), Some(2), None], &ctx());
        let x0 = [3.0_f64, 2.0, 0.0];

        // Pre-seed prev so fetlim doesn't bias the FD comparison.
        m.vgs_eff_prev = 2.0;
        m.vds_eff_prev = 3.0;
        m.eval(&x0, EvalFlags::dc(), &ctx());
        let gm_analytic  = m.gm;
        let gds_analytic = m.gds;

        let ids0 = m.jeq + gm_analytic * (x0[1] - x0[2]) + gds_analytic * (x0[0] - x0[2]);

        // Perturb VG (+eps on node 1).
        let mut xg = x0;
        xg[1] += eps;
        m.eval(&xg, EvalFlags::dc(), &ctx());
        let ids_g = m.jeq + m.gm * (xg[1] - xg[2]) + m.gds * (xg[0] - xg[2]);
        let gm_fd = (ids_g - ids0) / eps;

        // Perturb VD (+eps on node 0).
        let mut xd = x0;
        xd[0] += eps;
        m.eval(&xd, EvalFlags::dc(), &ctx());
        let ids_d = m.jeq + m.gm * (xd[1] - xd[2]) + m.gds * (xd[0] - xd[2]);
        let gds_fd = (ids_d - ids0) / eps;

        // Re-evaluate at x0 for clean comparison.
        m.eval(&x0, EvalFlags::dc(), &ctx());
        assert!(
            (gm_analytic - gm_fd).abs() / gm_analytic.abs() < 1e-4,
            "gm analytic={:.4e} fd={:.4e}", gm_analytic, gm_fd
        );
        assert!(
            (gds_analytic - gds_fd).abs().max(1e-14) / (gds_analytic.abs() + GMIN) < 0.01,
            "gds analytic={:.4e} fd={:.4e}", gds_analytic, gds_fd
        );
    }
}
