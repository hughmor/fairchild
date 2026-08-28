use crate::device::{Device, Discretisation, EvalFlags, NodeId, SimContext};
use crate::mna::{Cell, MnaMatrix, Pattern};
use crate::reactive::ChargeHistory;

/// SPICE's default channel width and length, 100 um each.
///
/// Public because [`crate::binning`] has to pick a bin using the geometry the
/// device will actually evaluate at. Two spellings of one default is two chances
/// for the bin and the model to disagree about which device this is.
pub const DEFAULT_W_M: f64 = 1e-4;
/// See [`DEFAULT_W_M`].
pub const DEFAULT_L_M: f64 = 1e-4;

/// Seed floor before the first `eval`. The operating value is
/// [`SimContext::gmin`] — see the note beside `gds_total` in `eval`.
const GMIN_SEED: f64 = 1e-12;

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
    mj: f64,   // junction grading coefficient (bottom of the junction)
    mjsw: f64, // sidewall grading coefficient (SPICE default 0.33)
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
    // Bottom and sidewall halves of each bulk junction, kept apart because they
    // grade differently: the bottom with MJ, the sidewall with MJSW.  They used
    // to be summed here, which is how MJSW came to be parsed, stored, and never
    // read — a card setting MJ=0.5 MJSW=0.33 silently got 0.5 for both.
    cbs_bot: f64, // CJ   * AS
    cbs_sw: f64,  // CJSW * PS
    cbd_bot: f64, // CJ   * AD
    cbd_sw: f64,  // CJSW * PD

    // ── Terminal bindings ────────────────────────────────────────────────────
    drain: NodeId,
    gate: NodeId,
    source: NodeId,
    bulk: NodeId,

    /// The eight cells `load_jacobian` writes, resolved once against a
    /// pattern, with that pattern's id. `None` until `resolve_cells` runs, and
    /// unused unless the id matches the matrix being stamped — see
    /// [`Device::resolve_cells`]. A ring oscillator's supply row carries a
    /// column per stage, so searching it eight times per transistor per Newton
    /// iteration was the single largest line in the assembly profile.
    jac_cells: Option<(u64, [Option<Cell>; 8])>,

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
    q_gs_eval: f64, // Q_gs = Cgs * Vgs at the current NR iterate
    q_gd_eval: f64,
    q_gb_eval: f64,
    gs_hist: ChargeHistory, // history at the previous accepted timestep
    gd_hist: ChargeHistory,
    gb_hist: ChargeHistory,

    // ── Junction cap transient state (nonlinear depletion cap) ───────────────
    vbs_eval: f64,  // Vbs = Vb−Vs at current iterate (for stamp use)
    vbd_eval: f64,  // Vbd = Vb−Vd
    cbs_eval: f64,  // Cbs(Vbs)
    cbd_eval: f64,  // Cbd(Vbd)
    q_bs_eval: f64, // Q_bs(Vbs) charge integral at current iterate
    q_bd_eval: f64,
    bs_hist: ChargeHistory, // history at the previous accepted timestep
    bd_hist: ChargeHistory,

    /// The integrator's discretisation, captured during `eval` because
    /// `load_*_tran` receives only `alpha` — which is Backward Euler and
    /// nothing else. `None` outside the transient loop.
    disc: Option<Discretisation>,
}

impl Mosfet1 {
    /// The eight `(row, col)` pairs `load_jacobian` writes, in the order its
    /// value array uses.
    ///
    /// One definition, read by both the searching stamp and by
    /// `resolve_cells`, so the resolved cells cannot address a different set
    /// of pairs — or the same pairs in a different order — than the values
    /// they are paired with. Two lists here would be two chances to disagree,
    /// and the disagreement would be a plausible number in the wrong cell.
    fn jac_pairs(d: NodeId, g: NodeId, s: NodeId, bk: NodeId) -> [(NodeId, NodeId); 8] {
        [
            (d, g),
            (d, d),
            (d, s),
            (d, bk),
            (s, g),
            (s, d),
            (s, s),
            (s, bk),
        ]
    }

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
                // Accepted and NOT modelled: `crate::unmodelled` owns that list
                // and the diagnostic that reads it.
                k if crate::unmodelled::is_listed(crate::unmodelled::MOSFET, k) => {}
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
            cbs_bot: 0.0,
            cbs_sw: 0.0,
            cbd_bot: 0.0,
            cbd_sw: 0.0,
            drain: None,
            gate: None,
            source: None,
            bulk: None,
            jac_cells: None,
            gm: GMIN_SEED,
            gds: GMIN_SEED,
            gmbs: 0.0,
            jeq: 0.0,
            vgs_eff_prev: 0.0,
            vds_eff_prev: 0.0,
            cgs_eval: 0.0,
            cgd_eval: 0.0,
            cgb_eval: 0.0,
            q_gs_eval: 0.0,
            q_gd_eval: 0.0,
            q_gb_eval: 0.0,
            gs_hist: ChargeHistory::default(),
            gd_hist: ChargeHistory::default(),
            gb_hist: ChargeHistory::default(),
            vbs_eval: 0.0,
            vbd_eval: 0.0,
            cbs_eval: 0.0,
            cbd_eval: 0.0,
            q_bs_eval: 0.0,
            q_bd_eval: 0.0,
            bs_hist: ChargeHistory::default(),
            bd_hist: ChargeHistory::default(),
            disc: None,
        };
        (dev, unknown)
    }

    /// The bulk-source depletion capacitance from the last `eval`.
    ///
    /// Exposed for the sidewall-grading test, which has to see the two graded
    /// halves summed the way the matrix sees them: an internal-consistency check
    /// between MJ and MJSW would pass with either coefficient used for both.
    pub fn cbs_at_last_eval(&self) -> f64 {
        self.cbs_eval
    }

    /// Apply instance parameters (W, L, AS, AD, PS, PD).
    pub fn set_instance_params(&mut self, params: &[(String, f64)]) -> Vec<String> {
        let mut w = DEFAULT_W_M;
        let mut l = DEFAULT_L_M;
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
        self.cbs_bot = self.cj * as_;
        self.cbs_sw = self.cjsw * ps;
        self.cbd_bot = self.cj * ad;
        self.cbd_sw = self.cjsw * pd;
        unknown
    }

    // ── Depletion cap helpers (same model as ShockleyDiode::cj_depl / q_depl) ─

    /// One graded junction term.  `m` is the grading coefficient — `MJ` for the
    /// bottom of the junction, `MJSW` for its sidewall.
    fn cj_depl_m(&self, c0: f64, v: f64, m: f64) -> f64 {
        if c0 == 0.0 {
            return 0.0;
        }
        let fc_pb = self.fc * self.pb;
        if v < fc_pb {
            c0 * (1.0 - v / self.pb).powf(-m)
        } else {
            let k = (1.0 - self.fc).powf(1.0 + m);
            c0 / k * (1.0 - self.fc * (1.0 + m) + m * v / self.pb)
        }
    }

    fn q_depl_m(&self, c0: f64, v: f64, m: f64) -> f64 {
        if c0 == 0.0 {
            return 0.0;
        }
        let fc_pb = self.fc * self.pb;
        if v < fc_pb {
            let x = 1.0 - v / self.pb;
            c0 * self.pb / (1.0 - m) * (1.0 - x.powf(1.0 - m))
        } else {
            let x_fc = 1.0 - self.fc;
            let q_fc = c0 * self.pb / (1.0 - m) * (1.0 - x_fc.powf(1.0 - m));
            let k = x_fc.powf(1.0 + m);
            let f2 = 1.0 - self.fc * (1.0 + m);
            let dv = v - fc_pb;
            q_fc + c0 / k * (f2 * dv + m / (2.0 * self.pb) * (v * v - fc_pb * fc_pb))
        }
    }

    /// A whole bulk junction: bottom graded with `MJ`, sidewall with `MJSW`.
    fn cj_depl(&self, bot: f64, sw: f64, v: f64) -> f64 {
        self.cj_depl_m(bot, v, self.mj) + self.cj_depl_m(sw, v, self.mjsw)
    }

    fn q_depl(&self, bot: f64, sw: f64, v: f64) -> f64 {
        self.q_depl_m(bot, v, self.mj) + self.q_depl_m(sw, v, self.mjsw)
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
        // A conditioning floor on the *channel*, and deliberately
        // Jacobian-only: `jeq` below subtracts `gds_total·vds`, so this term
        // cancels exactly out of the terminal current and the operating point
        // does not depend on it. That is what a conditioning floor should do,
        // and it is the opposite of what `diode.rs`/`bjt.rs` needed — those have
        // a pn junction, and SPICE's `GMIN` goes *across* a junction and carries
        // current. Level 1 has no body diodes here (`IS`/`JS` are on
        // `docs/model_status.md`'s unmodelled list), so there is no junction on
        // this device to put one across, and a reverse-biased drain-bulk carries
        // no `gmin` leakage where ngspice's does.
        //
        // The value comes from the solve so that `.options gmin=1e-5`, raised to
        // get a stubborn circuit through, actually reaches the MOSFETs.
        let gds_total = gds_eff + ctx.gmin;
        self.gm = gm_eff;
        self.gds = gds_total;
        self.gmbs = gmbs_eff;
        self.jeq = ids_real - gm_eff * vgs - gds_total * vds - gmbs_eff * vbs;

        // ── Capacitance evaluation ──────────────────────────────────────────
        self.disc = ctx.discretisation;

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
            // Meyer caps are linear, so the charge at this iterate is just C·V —
            // but the companion needs it explicitly, not folded into `alpha·Q_prev`.
            self.q_gs_eval = self.cgs_eval * (vg - vs);
            self.q_gd_eval = self.cgd_eval * (vg - vd);
            self.q_gb_eval = self.cgb_eval * (vg - vb);

            // Junction caps: depletion model.
            // Use real (non-pol-flipped) junction voltages.
            // V < 0 → reverse-biased (normal operation).
            let vbs_j = vb - vs;
            let vbd_j = vb - vd;
            self.vbs_eval = vbs_j;
            self.vbd_eval = vbd_j;
            self.cbs_eval = self.cj_depl(self.cbs_bot, self.cbs_sw, vbs_j);
            self.cbd_eval = self.cj_depl(self.cbd_bot, self.cbd_sw, vbd_j);
            self.q_bs_eval = self.q_depl(self.cbs_bot, self.cbs_sw, vbs_j);
            self.q_bd_eval = self.q_depl(self.cbd_bot, self.cbd_sw, vbd_j);
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
        // The eight values, in the one order `JAC_PAIRS` fixes. Both arms
        // below consume this slice, so the fast path cannot drift from the
        // slow one by reordering: there is only one order.
        let vals = [gm, gds, -gms, gmbs, -gm, -gds, gms, -gmbs];

        // Resolved cells, but only if they belong to *this* matrix. A
        // patternless matrix reports id 0 and never matches, so the diagnostic
        // passes that stamp into `MnaMatrix::zeros` fall through correctly.
        if let Some((id, cells)) = &self.jac_cells {
            if *id == mat.pattern_id() {
                for (cell, v) in cells.iter().zip(vals) {
                    if let Some(c) = cell {
                        mat.add(*c, v);
                    }
                }
                return;
            }
        }

        for (&(ri, ci), v) in Self::jac_pairs(d, g, s, bk).iter().zip(vals) {
            if let (Some(r), Some(c)) = (ri, ci) {
                mat.a[r][c] += v;
            }
        }
    }

    fn resolve_cells(&mut self, pattern: &Pattern) {
        let pairs = Self::jac_pairs(self.drain, self.gate, self.source, self.bulk);
        let mut cells = [None; 8];
        for (cell, &(r, c)) in cells.iter_mut().zip(pairs.iter()) {
            // `None` here means ground, or a cell outside the pattern. Either
            // way the searching path handles it: ground is skipped, and an
            // undeclared cell has to be *inserted*, which only `IndexMut` does.
            *cell = pattern.cell(r, c);
            debug_assert!(
                cell.is_some() || r.is_none() || c.is_none(),
                "MOSFET stamps ({r:?}, {c:?}), which its footprint did not declare"
            );
        }
        self.jac_cells = Some((pattern.id(), cells));
    }

    fn load_residual_tran(&self, b: &mut [f64], alpha: f64) {
        self.load_residual(b);
        let (d, g, s, bk) = (self.drain, self.gate, self.source, self.bulk);
        let disc = self.disc;

        // Every cap here is a charge branch, so one companion serves all five:
        // the history current cancels whatever the Jacobian stamp contributes
        // to the residual, under whichever method the integrator chose. Under
        // Backward Euler the gate caps reduce to the old `alpha·Q_prev` and the
        // junction caps to `alpha·(C·V + Q_prev − Q_now)`.
        //
        // Gate caps are linear (Q = C·V), so charge and the Jacobian term
        // coincide; the depletion caps' do not, which is why both are passed.
        let mut cap = |hist: &ChargeHistory, a, bn, q_new, cv| {
            let (i_hist, _) = hist.companion(disc, alpha, q_new, cv);
            Self::stamp_hist(b, a, bn, i_hist);
        };
        cap(&self.gs_hist, g, s, self.q_gs_eval, self.q_gs_eval);
        cap(&self.gd_hist, g, d, self.q_gd_eval, self.q_gd_eval);
        cap(&self.gb_hist, g, bk, self.q_gb_eval, self.q_gb_eval);
        if self.cbs_bot != 0.0 || self.cbs_sw != 0.0 {
            let cv = self.cbs_eval * self.vbs_eval;
            cap(&self.bs_hist, bk, s, self.q_bs_eval, cv);
        }
        if self.cbd_bot != 0.0 || self.cbd_sw != 0.0 {
            let cv = self.cbd_eval * self.vbd_eval;
            cap(&self.bd_hist, bk, d, self.q_bd_eval, cv);
        }
    }

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, alpha: f64) {
        self.load_jacobian(mat);
        let (d, g, s, bk) = (self.drain, self.gate, self.source, self.bulk);
        // Same factor the residual's companion used, from the same place.
        let scale = ChargeHistory::scale(self.disc, alpha);

        // Gate caps.
        if self.cgs_eval != 0.0 {
            Self::stamp_g_pair(mat, g, s, scale * self.cgs_eval);
        }
        if self.cgd_eval != 0.0 {
            Self::stamp_g_pair(mat, g, d, scale * self.cgd_eval);
        }
        if self.cgb_eval != 0.0 {
            Self::stamp_g_pair(mat, g, bk, scale * self.cgb_eval);
        }

        // Junction caps.
        if self.cbs_eval != 0.0 {
            Self::stamp_g_pair(mat, bk, s, scale * self.cbs_eval);
        }
        if self.cbd_eval != 0.0 {
            Self::stamp_g_pair(mat, bk, d, scale * self.cbd_eval);
        }
    }

    /// Small-signal caps for `.ac`/`.noise`: the Meyer gate caps (Cgs/Cgd/Cgb)
    /// and depletion junction caps (Cbs/Cbd) at the operating point — the same
    /// five caps `load_jacobian_tran` stamps. Requires a preceding
    /// `eval(EvalFlags::tran())` to populate the `*_eval` caches.
    fn small_signal_reactances(&self) -> Vec<crate::device::ReactiveBranchSpec> {
        use crate::device::{ReactiveBranchSpec, ReactiveKind};
        let cap = |pos, neg, value| ReactiveBranchSpec {
            kind: ReactiveKind::Capacitor,
            pos,
            neg,
            value,
            // These feed `.ac`/`.noise`, which want the small-signal C itself
            // and not the charge branch's `∂q/∂v`, so zero is correct here
            // rather than merely conservative.
            dvalue_dstate: 0.0,
        };
        let mut v = Vec::new();
        if self.cgs_eval != 0.0 {
            v.push(cap(self.gate, self.source, self.cgs_eval));
        }
        if self.cgd_eval != 0.0 {
            v.push(cap(self.gate, self.drain, self.cgd_eval));
        }
        if self.cgb_eval != 0.0 {
            v.push(cap(self.gate, self.bulk, self.cgb_eval));
        }
        if self.cbs_eval != 0.0 {
            v.push(cap(self.bulk, self.source, self.cbs_eval));
        }
        if self.cbd_eval != 0.0 {
            v.push(cap(self.bulk, self.drain, self.cbd_eval));
        }
        v
    }

    fn commit_timestep(&mut self, x: &[f64]) {
        let vd = self.drain.map_or(0.0, |i| x[i]);
        let vg = self.gate.map_or(0.0, |i| x[i]);
        let vs = self.source.map_or(0.0, |i| x[i]);
        let vb = self.bulk.map_or(0.0, |i| x[i]);

        // Roll each cap's history. Charges are recomputed from the converged
        // solution rather than reused from the last `eval`, which is one NR
        // iterate behind it.
        let disc = self.disc;
        // Gate caps (linear: Q = C·V).
        self.gs_hist.advance(disc, self.cgs_eval * (vg - vs));
        self.gd_hist.advance(disc, self.cgd_eval * (vg - vd));
        self.gb_hist.advance(disc, self.cgb_eval * (vg - vb));
        // Junction caps (nonlinear depletion cap).
        let vbs_j = vb - vs;
        let vbd_j = vb - vd;
        self.bs_hist
            .advance(disc, self.q_depl(self.cbs_bot, self.cbs_sw, vbs_j));
        self.bd_hist
            .advance(disc, self.q_depl(self.cbd_bot, self.cbd_sw, vbd_j));

        // Update fetlim reference.
        self.vgs_eff_prev = self.polarity * (vg - vs);
        self.vds_eff_prev = self.polarity * (vd - vs);
    }

    fn noise_sources(&self, ctx: &SimContext, _freq: f64) -> Vec<(NodeId, NodeId, f64)> {
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
        assert!(
            (m.gds - ctx().gmin).abs() < 1e-20,
            "gds in cutoff: {}",
            m.gds
        );
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
            (gds_analytic - gds_fd).abs().max(1e-14) / (gds_analytic.abs() + GMIN_SEED) < 0.01,
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
