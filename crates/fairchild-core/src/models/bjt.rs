//! Gummel-Poon Level 1 bipolar junction transistor (BJT).
//!
//! Implements the Ebers-Moll transport form with the following model parameters:
//! IS, BF, BR, NF, NR, VAF, VAR, TF, TR.  Series resistances (RB, RC, RE)
//! and junction capacitances (CJE, CJC) are accepted and noted for future use;
//! the current implementation treats them as zero.
//!
//! ## Sign convention and polarity
//!
//! SPICE netlist order: `Q<name> C B E [S] model`.  For NPN in forward active:
//! - IC flows INTO the collector terminal (device sinks IC from C).
//! - IB flows INTO the base terminal (device sinks IB from B).
//! - IE flows OUT OF the emitter terminal (device sources IC+IB into E).
//!
//! PNP uses `polarity = -1`: all terminal-voltage differences and currents are
//! sign-flipped before the shared NPN physics equations run.  The MNA Jacobian
//! stamps are polarity-independent (pol² = 1), so the same `load_jacobian`
//! works for both device types.

use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

const GMIN: f64 = 1e-12;

/// Gummel-Poon Level 1 BJT.
pub struct GummelPoonBjt {
    // ── Model parameters ──────────────────────────────────────────────────────
    is:      f64,   // transport saturation current (A)
    bf:      f64,   // forward beta (current gain)
    br:      f64,   // reverse beta
    nf:      f64,   // forward emission coefficient
    nr:      f64,   // reverse emission coefficient
    vaf:     f64,   // forward Early voltage (V); f64::INFINITY = no Early effect
    var:     f64,   // reverse Early voltage (V)
    tf:      f64,   // forward transit time (s) — B-E diffusion charge
    tr:      f64,   // reverse transit time (s) — B-C diffusion charge
    polarity: f64,  // +1 NPN, -1 PNP
    vcrit:   f64,   // pnjlim critical voltage (derived)

    // ── Terminal bindings ─────────────────────────────────────────────────────
    collector: NodeId,
    base:      NodeId,
    emitter:   NodeId,

    // ── Cached per-NR-iteration quantities (set by eval) ─────────────────────
    vbe_eff:  f64,  // effective junction voltage B-E (after pnjlim, polarity-corrected)
    vbc_eff:  f64,  // effective junction voltage B-C
    gf:       f64,  // dIF/dVBE_eff (≈ transconductance in forward active)
    gr:       f64,  // dIR/dVBC_eff
    gpi:      f64,  // dIBF/dVBE_eff = gf/BF
    gmu:      f64,  // dIBR/dVBC_eff = gr/BR
    jeq_c:    f64,  // Norton offset at collector (see load_residual)
    jeq_b:    f64,  // Norton offset at base

    // ── pnjlim history ────────────────────────────────────────────────────────
    vbe_prev: f64,
    vbc_prev: f64,

    // ── Transient charge history ──────────────────────────────────────────────
    qbe_tprev: f64, // QBE = TF*IF at last accepted timestep
    qbc_tprev: f64, // QBC = TR*IR at last accepted timestep
    // Current-iterate charges (set by eval when transient flag is set)
    qbe_now:   f64,
    qbc_now:   f64,
    cbe_eff:   f64, // dQBE/dVBE_eff = TF*gf (capacitive conductance)
    cbc_eff:   f64, // dQBC/dVBC_eff = TR*gr
}

impl GummelPoonBjt {
    /// Build from model-card parameters.  Returns the device and any
    /// unrecognised parameter names.
    pub fn from_model_params(is_pnp: bool, params: &[(String, f64)]) -> (Self, Vec<String>) {
        let mut is  = 1e-16;
        let mut bf  = 100.0;
        let mut br  = 1.0;
        let mut nf  = 1.0;
        let mut nr  = 1.0;
        let mut vaf = f64::INFINITY;
        let mut var = f64::INFINITY;
        let mut tf  = 0.0;
        let mut tr  = 0.0;
        let mut unknown = Vec::new();
        for (k, v) in params {
            match k.to_lowercase().as_str() {
                "is"           => is  = *v,
                "bf" | "hfe"  => bf  = *v,
                "br" | "hrc"  => br  = *v,
                "nf"           => nf  = *v,
                "nr"           => nr  = *v,
                "vaf" | "va"  => vaf = *v,
                "var" | "vb"  => var = *v,
                "tf"           => tf  = *v,
                "tr"           => tr  = *v,
                // Accepted but not yet modelled.
                "rb" | "rc" | "re" | "ikf" | "ikr" | "ise" | "isc"
                | "ne" | "nc" | "cje" | "cjc" | "vje" | "vjc" | "mje" | "mjc"
                | "fc" | "cjs" | "vjs" | "mjs" | "xtb" | "eg" | "xti"
                | "kf" | "af" | "ptf" | "xcjc" | "tnom" => {}
                _ => unknown.push(k.clone()),
            }
        }
        let dev = GummelPoonBjt {
            is, bf, br, nf, nr, vaf, var, tf, tr,
            polarity: if is_pnp { -1.0 } else { 1.0 },
            vcrit: 0.0,
            collector: None, base: None, emitter: None,
            vbe_eff: 0.0, vbc_eff: 0.0,
            gf: GMIN, gr: GMIN, gpi: GMIN / 100.0, gmu: GMIN / 100.0,
            jeq_c: 0.0, jeq_b: 0.0,
            vbe_prev: 0.0, vbc_prev: 0.0,
            qbe_tprev: 0.0, qbc_tprev: 0.0,
            qbe_now: 0.0, qbc_now: 0.0,
            cbe_eff: 0.0, cbc_eff: 0.0,
        };
        (dev, unknown)
    }

    /// pnjlim: logarithmic compression of large junction-voltage steps.
    fn pnjlim(&self, vnew: f64, vold: f64, vt: f64) -> f64 {
        if vnew > self.vcrit && (vnew - vold).abs() > 2.0 * vt {
            if vnew > vold {
                vold + vt * ((vnew - vold) / vt + 1.0).ln()
            } else {
                vold - vt * ((vold - vnew) / vt + 1.0).ln()
            }
        } else {
            vnew
        }
    }
}

impl Device for GummelPoonBjt {
    fn num_terminals(&self) -> usize { 4 } // C B E S (substrate tied to ground by build_devices)

    fn setup_model(&mut self, ctx: &SimContext) {
        let vt = ctx.vt();
        self.vcrit = vt * (vt / (std::f64::consts::SQRT_2 * self.is)).ln();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], _ctx: &SimContext) {
        // Terminals order: [C, B, E, S] — substrate (S) is ignored in this implementation.
        debug_assert!(terminals.len() >= 3, "BJT expects [C, B, E, S]");
        self.collector = terminals[0];
        self.base      = terminals[1];
        self.emitter   = terminals[2];
        // terminals[3] = substrate — tied to ground by caller; not stamped separately.
    }

    fn eval(&mut self, x: &[f64], flags: EvalFlags, ctx: &SimContext) {
        let pol = self.polarity;
        let vc = self.collector.map_or(0.0, |i| x[i]);
        let vb = self.base     .map_or(0.0, |i| x[i]);
        let ve = self.emitter  .map_or(0.0, |i| x[i]);

        let vt = ctx.vt();

        // Polarity-flipped junction voltages (positive in forward active for both NPN and PNP).
        let vbe_raw = pol * (vb - ve);
        let vbc_raw = pol * (vb - vc);

        let vbe_eff = if ctx.jlim_enabled {
            self.pnjlim(vbe_raw, self.vbe_prev, vt)
        } else { vbe_raw };
        let vbc_eff = if ctx.jlim_enabled {
            self.pnjlim(vbc_raw, self.vbc_prev, vt)
        } else { vbc_raw };
        self.vbe_prev = vbe_eff;
        self.vbc_prev = vbc_eff;
        self.vbe_eff  = vbe_eff;
        self.vbc_eff  = vbc_eff;

        // Forward and reverse junction exponentials.
        let nf_vt = self.nf * vt;
        let nr_vt = self.nr * vt;
        let exp_be = (vbe_eff / nf_vt).exp();
        let exp_bc = (vbc_eff / nr_vt).exp();

        let if_val = self.is * (exp_be - 1.0);
        let ir_val = self.is * (exp_bc - 1.0);

        // Conductances (including GMIN floor).
        let gf = self.is * exp_be / nf_vt + GMIN;
        let gr = self.is * exp_bc / nr_vt + GMIN;
        let gpi = gf / self.bf;
        let gmu = gr / self.br;

        // Early voltage modulation: q1 = 1 / (1 − VBC/VAF − VBE/VAR)
        // Applied only when VAF/VAR < ∞.  Clamped to ≥0.1 to prevent divergence.
        let q1 = if self.vaf.is_finite() || self.var.is_finite() {
            let vbc_o_vaf = if self.vaf.is_finite() { vbc_eff / self.vaf } else { 0.0 };
            let vbe_o_var = if self.var.is_finite() { vbe_eff / self.var } else { 0.0 };
            (1.0 - vbc_o_vaf - vbe_o_var).max(0.1)
        } else { 1.0 };

        // Effective collector and base currents (in NPN-equivalent space).
        let ic_eff = (if_val - ir_val) / q1 - ir_val / self.br;
        let ib_eff = if_val / self.bf + ir_val / self.br;

        // Real currents into each physical terminal: pol * eff.
        // b[C] -= pol*ic_eff; b[B] -= pol*ib_eff; b[E] += pol*(ic_eff + ib_eff)
        // The Norton offset absorbs the linear Jacobian contribution at eval-point voltages.
        // (vb_c, vb_b, ve_c are actual node voltages; see load_residual docs.)
        // Jacobian entries (pol²=1 so pol-independent):
        //   dIC_real/dVB = gf/q1 - gce, dIC_real/dVC = gce, dIC_real/dVE = -gf/q1
        //   dIB_real/dVB = gpi/q1+gmu,  dIB_real/dVC = -gmu, dIB_real/dVE = -gpi/q1
        //
        // jeq_C = pol*ic_eff − Jacobian_C · V_current
        // jeq_B = pol*ib_eff − (gpi+gmu)·vb + gmu·vc + gpi·ve
        self.gf  = gf  / q1; // Early-modulated transconductance
        self.gr  = gr;
        self.gpi = gpi / q1;
        self.gmu = gmu;
        let gce = self.gr + self.gmu;

        self.jeq_c = pol * ic_eff
            - (self.gf - gce) * vb
            - gce * vc
            + self.gf * ve;
        self.jeq_b = pol * ib_eff
            - (self.gpi + self.gmu) * vb
            + self.gmu * vc
            + self.gpi * ve;

        if flags.transient {
            self.cbe_eff = self.tf * self.gf;
            self.cbc_eff = self.tr * self.gr;
            self.qbe_now = self.tf * if_val / q1;
            self.qbc_now = self.tr * ir_val;
        } else {
            self.cbe_eff = 0.0;
            self.cbc_eff = 0.0;
            self.qbe_now = 0.0;
            self.qbc_now = 0.0;
        }
    }

    fn load_residual(&self, b: &mut [f64]) {
        if let Some(c) = self.collector { b[c] -= self.jeq_c; }
        if let Some(bk) = self.base    { b[bk] -= self.jeq_b; }
        if let Some(e) = self.emitter  { b[e]  += self.jeq_c + self.jeq_b; }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let (c, bk, e) = (self.collector, self.base, self.emitter);
        let gf  = self.gf;
        let gpi = self.gpi;
        let gmu = self.gmu;
        let gce = self.gr + gmu;  // combined collector-emitter conductance via VBC

        macro_rules! stamp {
            ($ri:expr, $ci:expr, $val:expr) => {
                if let (Some(r), Some(cc)) = ($ri, $ci) {
                    mat.a[r][cc] += $val;
                }
            };
        }

        // Collector row: IC flows C→E; Jacobian of (device draws pol*IC from C).
        stamp!(c,  bk,  gf - gce);
        stamp!(c,  c,   gce);
        stamp!(c,  e,  -gf);

        // Base row: IB flows into B.
        stamp!(bk, bk,  gpi + gmu);
        stamp!(bk, c,  -gmu);
        stamp!(bk, e,  -gpi);

        // Emitter row: IE = IC + IB sources into E (negative sum of above rows).
        stamp!(e,  bk, -(gf - gce) - (gpi + gmu));
        stamp!(e,  c,  -gce - (-gmu));
        stamp!(e,  e,   gf + gpi);
    }

    fn load_residual_tran(&self, b: &mut [f64], alpha: f64) {
        self.load_residual(b);
        // B-E junction charge companion (transit-time diffusion): TF·IF
        if self.cbe_eff != 0.0 {
            let i_be = alpha * (self.cbe_eff * self.vbe_eff + self.qbe_tprev - self.qbe_now);
            // Companion current flows from E→B (into B, out of E); same polarity as junction.
            let pol = self.polarity;
            if let Some(bk) = self.base    { b[bk] += pol * i_be; }
            if let Some(e)  = self.emitter { b[e]  -= pol * i_be; }
        }
        // B-C junction charge companion: TR·IR
        if self.cbc_eff != 0.0 {
            let i_bc = alpha * (self.cbc_eff * self.vbc_eff + self.qbc_tprev - self.qbc_now);
            let pol = self.polarity;
            if let Some(bk) = self.base      { b[bk] += pol * i_bc; }
            if let Some(c)  = self.collector { b[c]  -= pol * i_bc; }
        }
    }

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, alpha: f64) {
        self.load_jacobian(mat);
        let (c, bk, e) = (self.collector, self.base, self.emitter);

        macro_rules! stamp {
            ($ri:expr, $ci:expr, $val:expr) => {
                if let (Some(r), Some(cc)) = ($ri, $ci) {
                    mat.a[r][cc] += $val;
                }
            };
        }

        // B-E capacitive companion: cbe_eff between base and emitter.
        if self.cbe_eff != 0.0 {
            let c_be = alpha * self.cbe_eff;
            stamp!(bk, bk,  c_be);
            stamp!(bk, e,  -c_be);
            stamp!(e,  bk, -c_be);
            stamp!(e,  e,   c_be);
        }
        // B-C capacitive companion: cbc_eff between base and collector.
        if self.cbc_eff != 0.0 {
            let c_bc = alpha * self.cbc_eff;
            stamp!(bk, bk,  c_bc);
            stamp!(bk, c,  -c_bc);
            stamp!(c,  bk, -c_bc);
            stamp!(c,  c,   c_bc);
        }
    }

    fn commit_timestep(&mut self, x: &[f64]) {
        let pol = self.polarity;
        let vc = self.collector.map_or(0.0, |i| x[i]);
        let vb = self.base     .map_or(0.0, |i| x[i]);
        let ve = self.emitter  .map_or(0.0, |i| x[i]);
        let vbe_eff = pol * (vb - ve);
        let vbc_eff = pol * (vb - vc);
        self.vbe_prev = vbe_eff;
        self.vbc_prev = vbc_eff;
        let vt = 0.02585; // approximate; full ctx not available at commit time
        let exp_be = (vbe_eff / (self.nf * vt)).exp();
        let exp_bc = (vbc_eff / (self.nr * vt)).exp();
        let if_val = self.is * (exp_be - 1.0);
        let ir_val = self.is * (exp_bc - 1.0);
        self.qbe_tprev = self.tf * if_val;
        self.qbc_tprev = self.tr * ir_val;
    }

    fn noise_sources(&self, ctx: &SimContext) -> Vec<(NodeId, NodeId, f64)> {
        // Shot noise on B-E and B-C junctions.
        // i_n_be² = 2q|IB|, flows base→emitter.
        // i_n_ce² = 2q|IC| (collector shot noise), flows collector→emitter.
        const Q_E: f64 = 1.602176634e-19;
        let _ = ctx;
        let ic_approx = self.polarity * (self.jeq_c
            + (self.gf - self.gr - self.gmu) * 0.0
            + (self.gr + self.gmu) * 0.0
            - self.gf * 0.0);
        let ib_approx = self.polarity * self.jeq_b;
        let mut sources = Vec::new();
        if ic_approx.abs() > 1e-20 {
            sources.push((self.collector, self.emitter, 2.0 * Q_E * ic_approx.abs()));
        }
        if ib_approx.abs() > 1e-20 {
            sources.push((self.base, self.emitter, 2.0 * Q_E * ib_approx.abs()));
        }
        sources
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::device::EvalFlags;

    fn ctx() -> SimContext { SimContext::default() }

    fn npn(bf: f64) -> GummelPoonBjt {
        let (mut q, _) = GummelPoonBjt::from_model_params(false, &[
            ("is".into(), 1e-15),
            ("bf".into(), bf),
            ("br".into(), 1.0),
        ]);
        q.setup_model(&ctx());
        // C=0, B=1, E=2, S=None
        q.setup_instance(&[Some(0), Some(1), Some(2), None], &ctx());
        q
    }

    #[test]
    fn npn_forward_active_ic() {
        // VBE = 0.7 V, VCE = 5 V → forward active.
        // IC ≈ IS * exp(VBE/VT) = 1e-15 * exp(0.7/0.025865) ≈ 0.567 mA  (BF=100)
        let bf = 100.0;
        let mut q = npn(bf);
        // Pre-seed pnjlim history so the first eval isn't limited.
        q.vbe_prev = 0.7;
        q.vbc_prev = 0.7 - 5.0;
        let x = [5.0_f64, 0.7, 0.0]; // VC=5, VB=0.7, VE=0
        q.eval(&x, EvalFlags::dc(), &ctx());

        // Recover IC from Norton stamp: IC = jeq_C + Jacobian·V
        let vb = 0.7_f64; let vc = 5.0_f64; let ve = 0.0_f64;
        let gce = q.gr + q.gmu;
        let ic = q.jeq_c + (q.gf - gce) * vb + gce * vc - q.gf * ve;
        // Expected: IS*exp(VBE/VT); use the same VT as the model (300.15 K default).
        let vt = ctx().vt();
        let ic_expected = 1e-15 * (0.7_f64 / vt).exp();
        assert!(
            (ic - ic_expected).abs() / ic_expected < 0.01,
            "IC={:.4e} expected≈{:.4e}", ic, ic_expected
        );
    }

    #[test]
    fn npn_ic_ib_ratio_equals_bf() {
        let bf = 150.0;
        let mut q = npn(bf);
        q.vbe_prev = 0.65;
        q.vbc_prev = -4.0;
        let x = [5.0_f64, 0.65, 0.0];
        q.eval(&x, EvalFlags::dc(), &ctx());

        let vb = 0.65_f64; let vc = 5.0_f64; let ve = 0.0_f64;
        let gce = q.gr + q.gmu;
        let ic = q.jeq_c + (q.gf - gce) * vb + gce * vc - q.gf * ve;
        let ib = q.jeq_b + (q.gpi + q.gmu) * vb - q.gmu * vc - q.gpi * ve;

        // In forward active (VBC << 0), IC/IB ≈ BF.
        let beta_measured = ic / ib;
        assert!(
            (beta_measured - bf).abs() / bf < 0.01,
            "β = {:.1} expected {:.1}", beta_measured, bf
        );
    }

    #[test]
    fn pnp_active_currents_flow_correctly() {
        // PNP in forward active: VEB=0.7V (VE=5, VB=4.3), VCB>0 (VC=0, VB=4.3)
        // IC should be negative (flows OUT of collector = INTO collector node from circuit side)
        let (mut q, _) = GummelPoonBjt::from_model_params(true, &[
            ("is".into(), 1e-15),
            ("bf".into(), 100.0),
        ]);
        q.setup_model(&ctx());
        q.setup_instance(&[Some(0), Some(1), Some(2), None], &ctx());
        q.vbe_prev = -0.7; // raw (VB-VE) for PNP
        q.vbc_prev = 4.3;
        // VC=0, VB=4.3, VE=5.0  → VB-VE = -0.7, VB-VC = 4.3
        let x = [0.0_f64, 4.3, 5.0];
        q.eval(&x, EvalFlags::dc(), &ctx());

        let vb = 4.3_f64; let vc = 0.0_f64; let ve = 5.0_f64;
        let gce = q.gr + q.gmu;
        let ic_into_c = q.jeq_c + (q.gf - gce) * vb + gce * vc - q.gf * ve;
        // PNP: current flows from emitter to collector internally.
        // The collector node SOURCE (external) sees current flowing INTO it
        // when collector is the low-voltage terminal — IC_into_C is negative.
        assert!(ic_into_c < 0.0,
            "PNP IC_into_C should be negative, got {:.4e}", ic_into_c);
        // |IC| / |IB| ≈ BF = 100
        let ib_into_b = q.jeq_b + (q.gpi + q.gmu) * vb - q.gmu * vc - q.gpi * ve;
        let beta = ic_into_c.abs() / ib_into_b.abs();
        assert!((beta - 100.0).abs() / 100.0 < 0.02, "PNP β={:.1}", beta);
    }
}
