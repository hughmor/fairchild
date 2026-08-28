use crate::device::{Device, Discretisation, EvalFlags, NodeId, SimContext};
use crate::mna::{Cell, MnaMatrix, Pattern};
use crate::reactive::ChargeHistory;

/// `IS`'s default bulk junction saturation current, A — SPICE's value.
const IS_BULK_DEFAULT: f64 = 1e-14;

/// `UO`'s default carrier mobility, cm²/V·s — SPICE's value.
const UO_DEFAULT: f64 = 600.0;

/// `KP` when the card gives neither `KP` nor an oxide capacitance to derive it
/// from. SPICE's fallback.
const KP_FALLBACK: f64 = 2e-5;

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
    vto: f64,    // threshold voltage (V)
    kp: f64,     // process transconductance (A/V²)
    lambda: f64, // channel-length modulation (1/V)
    gamma: f64,  // body-effect coefficient (V^0.5)
    phi: f64,    // surface potential (V)
    /// The external drain and source terminals. The intrinsic channel is
    /// stamped between `drain`/`source`, which are the *internal* nodes and alias
    /// these when the matching series resistance is zero.
    drain_ext: NodeId,
    source_ext: NodeId,
    /// `RD`/`RS` — drain and source ohmic series resistance (Ω). `0` means no
    /// internal node, and the internal node then aliases the external terminal.
    rd: f64,
    rs: f64,
    /// Drain current at the current iterate, signed as the model computes it.
    ///
    /// Kept because flicker noise is driven by `|Id|` and `jeq` is a Norton
    /// *offset* — equal to the current only when every terminal sits at 0 V. The
    /// BJT carried exactly that bug once (see its `noise_sources`), so the
    /// current is stored rather than reconstructed.
    ids_eval: f64,
    /// `IS` — bulk junction saturation current (A), and `JS` — its density
    /// (A/m²). Per SPICE, `JS·area` wins when both `JS` and that junction's area
    /// are given; otherwise `IS` applies to both junctions.
    is_bulk: f64,
    js: f64,
    /// Saturation current of each bulk junction, resolved from `IS`/`JS` and the
    /// areas in `set_instance_params` — the two can differ, because `AS` and `AD`
    /// can.
    isat_bs: f64,
    isat_bd: f64,
    /// Multiplier on both `isat_*` at the operating temperature. Applied where
    /// they are used rather than folded into them, because `set_instance_params`
    /// resolves them from `IS`/`JS` and the areas and may run either side of
    /// `setup_model` — a factor folded in could be applied twice or not at all.
    isat_t_factor: f64,
    /// Bulk junction currents and conductances at the current iterate.
    gbs: f64,
    gbd: f64,
    ibs_eq: f64,
    ibd_eq: f64,
    /// `PB` and the two bulk-capacitance factors at the operating temperature.
    ///
    /// Two factors because `CJ` grades by `MJ` and `CJSW` by `MJSW`, and SPICE
    /// derives each correction from its own coefficient — one factor applied to
    /// both would be wrong for any card where they differ, which is every card
    /// that sets them.
    pb_t: f64,
    cj_t_factor: f64,
    cjsw_t_factor: f64,
    /// `KF`/`AF` — flicker noise coefficient and exponent.
    kf: f64,
    af: f64,
    /// `TNOM` — the temperature this card's parameters were extracted at.
    tnom: f64,
    /// The three temperature-shifted values, derived in `setup_model` from the
    /// nominal ones and `crate::temperature`.
    ///
    /// Held rather than recomputed per eval: `PHI(T)` costs two logs and an
    /// exponential and depends on nothing that moves inside a solve. Derived
    /// idempotently, so a second `setup_model` cannot compound the shift.
    ///
    /// `GAMMA` does not shift — it is a doping/oxide ratio — but it multiplies
    /// `sqrt(PHI(T))`, so the body effect moves with temperature anyway.
    vto_t: f64,
    kp_t: f64,
    phi_t: f64,
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
    jac_cells: Option<(u64, [Option<Cell>; 24])>,

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
    #[allow(clippy::too_many_arguments)]
    fn jac_pairs(
        d: NodeId,
        g: NodeId,
        s: NodeId,
        bk: NodeId,
        d_ext: NodeId,
        s_ext: NodeId,
    ) -> [(NodeId, NodeId); 24] {
        [
            // The intrinsic channel, between the internal nodes.
            (d, g),
            (d, d),
            (d, s),
            (d, bk),
            (s, g),
            (s, d),
            (s, s),
            (s, bk),
            // RD, between the external drain and the internal one.
            (d_ext, d_ext),
            (d_ext, d),
            (d, d_ext),
            (d, d),
            // RS, likewise.
            (s_ext, s_ext),
            (s_ext, s),
            (s, s_ext),
            (s, s),
            // The bulk-source junction, between the bulk and the internal source.
            (bk, bk),
            (bk, s),
            (s, bk),
            (s, s),
            // The bulk-drain junction.
            (bk, bk),
            (bk, d),
            (d, bk),
            (d, d),
        ]
    }

    /// Construct from model-card parameters.
    pub fn from_model_params(is_pmos: bool, params: &[(String, f64)]) -> (Self, Vec<String>) {
        let mut rd = 0.0_f64;
        let mut rs = 0.0_f64;
        let mut kf = 0.0_f64;
        let mut af = 1.0_f64;
        let mut tnom_c = crate::temperature::TNOM_DEFAULT_K - 273.15;
        let mut vto = if is_pmos { -0.7 } else { 0.7 };
        // `None` distinguishes "the card gave a KP" from "the card gave the
        // default", which matters because `UO` only derives `KP` in the second
        // case. SPICE's rule, and the reason it cannot be a bare `2e-5` here.
        let mut kp: Option<f64> = None;
        let mut uo = UO_DEFAULT;
        let mut is_bulk = IS_BULK_DEFAULT;
        let mut js = 0.0_f64;
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
                "kp" => kp = Some(*v),
                // Carrier mobility, cm²/V·s on the card as everywhere in SPICE.
                "uo" | "u0" => uo = *v,
                "is" => is_bulk = *v,
                "js" => js = *v,
                "lambda" => lambda = *v,
                "gamma" => gamma = *v,
                "phi" => phi = *v,
                // Degrees Celsius on the card, like `.temp`.
                "rd" => rd = *v,
                "rs" => rs = *v,
                "kf" => kf = *v,
                "af" => af = *v,
                "tnom" => tnom_c = *v,
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
        // `KP` given wins; otherwise derive it from the mobility and the oxide
        // capacitance, which is what SPICE does. `UO` is cm²/V·s on the card and
        // the expression wants m²/V·s, hence the 1e-4.
        //
        // With no `TOX`/`COX` there is no `COX` to multiply, and SPICE's fallback
        // `KP` applies. That is a real card shape — a deck giving neither `KP` nor
        // an oxide thickness — so it is a default rather than a refusal.
        let kp = kp.unwrap_or(if cox > 0.0 {
            uo * 1e-4 * cox
        } else {
            KP_FALLBACK
        });
        let dev = Mosfet1 {
            vto,
            kp,
            lambda,
            gamma,
            phi,
            polarity: if is_pmos { -1.0 } else { 1.0 },
            drain_ext: None,
            source_ext: None,
            rd,
            rs,
            is_bulk,
            js,
            isat_bs: IS_BULK_DEFAULT,
            isat_bd: IS_BULK_DEFAULT,
            isat_t_factor: 1.0,
            gbs: 0.0,
            gbd: 0.0,
            ibs_eq: 0.0,
            ibd_eq: 0.0,
            pb_t: 0.8,
            cj_t_factor: 1.0,
            cjsw_t_factor: 1.0,
            ids_eval: 0.0,
            kf,
            af,
            tnom: tnom_c + 273.15,
            // Nominal until `setup_model` runs, which is before any eval.
            vto_t: vto,
            kp_t: kp,
            phi_t: phi,
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
        // `JS·area` when both are given, else `IS`. The two junctions resolve
        // independently because `AS` and `AD` can differ.
        self.isat_bs = if self.js > 0.0 && as_ > 0.0 {
            self.js * as_
        } else {
            self.is_bulk
        };
        self.isat_bd = if self.js > 0.0 && ad > 0.0 {
            self.js * ad
        } else {
            self.is_bulk
        };
        self.cbs_bot = self.cj * as_;
        self.cbs_sw = self.cjsw * ps;
        self.cbd_bot = self.cj * ad;
        self.cbd_sw = self.cjsw * pd;
        unknown
    }

    /// One bulk junction's current and slope: the Shockley law plus `gmin`.
    ///
    /// The same law [`crate::models::diode::ShockleyDiode::junction`] uses, and
    /// deliberately the same shape — one law for "what is a pn junction" rather
    /// than a second spelling that can drift from the first.
    ///
    /// # Why not ngspice's reverse branch
    ///
    /// ngspice's MOS1 is flat at exactly `-Isat` from `-3·vt` outward, and inside
    /// `±3·vt` its *total* over the two junctions measures as one junction flat at
    /// `-Isat` plus one plain Shockley — matched to seven digits at every bias
    /// tried, and still present when the bulk-drain junction is held at -5 V. That
    /// asymmetry is a numerical convenience in the reference, not physics.
    ///
    /// Both pure choices sit the same distance from it. Inside the band, at -0.01 V
    /// with `IS = 1e-14`, ngspice reads 1.32e-14 A, Shockley-on-both 6.4e-15 and
    /// flat-on-both 2.0e-14 — a 6.8e-15 A difference either way. Outside the band
    /// they agree: Shockley is within 4e-4 relative at -0.2 V and exact by -0.5 V,
    /// because `exp(v/vt)` underflows toward zero and leaves `-Isat`.
    ///
    /// So this takes the smooth one. Shockley is the junction law, it is C¹ at zero
    /// where the flat branch has a kink, and it needs no second case to explain.
    ///
    /// # No step limiting
    ///
    /// The DC Newton already clamps every node update to `vmax + reltol·|x|`
    /// (`crate::newton`), so the forward exponential cannot be jumped into from a
    /// cold start. A bulk that *converges* forward past ~18 V overflows to a
    /// non-finite the solver reports, and that is a broken deck rather than an
    /// operating point. A limiter here would mean holding a previous junction
    /// voltage — state the outer Newton cannot see, which is the shape that
    /// produced the diode's `RS` error and two others in this tree.
    fn bulk_junction(isat: f64, v: f64, vt: f64, gmin: f64) -> (f64, f64) {
        if isat <= 0.0 {
            return (gmin * v, gmin);
        }
        let e = (v / vt).exp();
        (isat * (e - 1.0) + gmin * v, isat * e / vt + gmin)
    }

    // ── Depletion cap helpers (same model as ShockleyDiode::cj_depl / q_depl) ─

    /// One graded junction term.  `m` is the grading coefficient — `MJ` for the
    /// bottom of the junction, `MJSW` for its sidewall.
    fn cj_depl_m(&self, c0: f64, v: f64, m: f64) -> f64 {
        if c0 == 0.0 {
            return 0.0;
        }
        let fc_pb = self.fc * self.pb_t;
        if v < fc_pb {
            c0 * (1.0 - v / self.pb_t).powf(-m)
        } else {
            let k = (1.0 - self.fc).powf(1.0 + m);
            c0 / k * (1.0 - self.fc * (1.0 + m) + m * v / self.pb_t)
        }
    }

    fn q_depl_m(&self, c0: f64, v: f64, m: f64) -> f64 {
        if c0 == 0.0 {
            return 0.0;
        }
        let fc_pb = self.fc * self.pb_t;
        if v < fc_pb {
            let x = 1.0 - v / self.pb_t;
            c0 * self.pb_t / (1.0 - m) * (1.0 - x.powf(1.0 - m))
        } else {
            let x_fc = 1.0 - self.fc;
            let q_fc = c0 * self.pb_t / (1.0 - m) * (1.0 - x_fc.powf(1.0 - m));
            let k = x_fc.powf(1.0 + m);
            let f2 = 1.0 - self.fc * (1.0 + m);
            let dv = v - fc_pb;
            q_fc + c0 / k * (f2 * dv + m / (2.0 * self.pb_t) * (v * v - fc_pb * fc_pb))
        }
    }

    /// A whole bulk junction: bottom graded with `MJ`, sidewall with `MJSW`.
    fn cj_depl(&self, bot: f64, sw: f64, v: f64) -> f64 {
        // Each area's own factor: `CJ` grades by `MJ` and `CJSW` by `MJSW`, and
        // SPICE derives each correction from its own coefficient.
        self.cj_depl_m(bot * self.cj_t_factor, v, self.mj)
            + self.cj_depl_m(sw * self.cjsw_t_factor, v, self.mjsw)
    }

    fn q_depl(&self, bot: f64, sw: f64, v: f64) -> f64 {
        self.q_depl_m(bot * self.cj_t_factor, v, self.mj)
            + self.q_depl_m(sw * self.cjsw_t_factor, v, self.mjsw)
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

    fn setup_model(&mut self, ctx: &SimContext) {
        // The mobility law and the threshold shift. `.temp` used to reach this
        // model only through `vt`, which Level 1's DC current does not even use —
        // so a 125 C run returned the 27 C drain current to the last bit.
        //
        // Idempotent: derived from the nominal values every time, never from the
        // previous result.
        let t = ctx.temperature;
        self.kp_t = self.kp * crate::temperature::mobility_factor(t, self.tnom);
        self.phi_t = crate::temperature::scaled_phi(self.phi, t, self.tnom);
        self.vto_t = crate::temperature::scaled_vto(
            self.vto,
            self.gamma,
            self.phi,
            self.phi_t,
            t,
            self.tnom,
            self.polarity < 0.0,
        );
        // The bulk junctions. `PB` moves by the same law as `PHI`; the two
        // capacitance corrections each take their own grading coefficient.
        self.pb_t = crate::temperature::scaled_junction_potential(self.pb, t, self.tnom);
        self.cj_t_factor = crate::temperature::junction_cap_factor(self.pb, self.mj, t, self.tnom);
        self.cjsw_t_factor =
            crate::temperature::junction_cap_factor(self.pb, self.mjsw, t, self.tnom);
        // The junction saturation currents take a third law, neither the diode's
        // nor the BJT's - a MOSFET card has no `EG` and no `XTI`, so SPICE puts
        // the temperature-dependent bandgap in the exponent. See
        // `temperature::mos_junction_is_factor`.
        self.isat_t_factor = crate::temperature::mos_junction_is_factor(t, self.tnom);
    }

    fn setup_instance(&mut self, terminals: &[NodeId], _ctx: &SimContext) {
        debug_assert_eq!(terminals.len(), 4, "MOSFET expects [D, G, S, B]");
        self.drain_ext = terminals[0];
        self.gate = terminals[1];
        self.source_ext = terminals[2];
        self.bulk = terminals[3];
        // The intrinsic channel runs between the *internal* nodes. With no series
        // resistance those alias the external terminals, so a card without
        // `RD`/`RS` allocates no extra rows and stamps no extra conductances.
        self.drain = terminals[0];
        self.source = terminals[2];
    }

    /// One internal node per non-zero ohmic series resistance.
    ///
    /// Real rows rather than an analytic elimination, deliberately. `diode.rs`'s
    /// `RS` is the cautionary tale: eliminating a series resistance by iterating
    /// on a junction voltage the outer Newton cannot see read 2.7% low against
    /// ngspice, and the convergence test had no way to notice. A row costs one
    /// unknown; a hidden state costs a silent wrong answer.
    fn num_extra_nodes(&self) -> usize {
        (self.rd > 0.0) as usize + (self.rs > 0.0) as usize
    }

    /// Bind internal drain'/source' nodes. Order is fixed (drain, source) so the
    /// assignment is stable across rebuilds.
    fn bind_extra_nodes(&mut self, first_idx: usize) {
        let mut idx = first_idx;
        if self.rd > 0.0 {
            self.drain = Some(idx);
            idx += 1;
        }
        if self.rs > 0.0 {
            self.source = Some(idx);
        }
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
            let vto_eff = pol * self.vto_t;
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
        let phi_m_vbs = (self.phi_t - vbs_eff).max(1e-10);
        let vto_eff = pol * self.vto_t;
        let vth = vto_eff + self.gamma * (phi_m_vbs.sqrt() - self.phi_t.sqrt());

        let (ids_eff, gm_eff, gds_eff, gmbs_eff) = if vgs_eff < vth {
            (0.0, 0.0, 0.0, 0.0)
        } else {
            let vdsat = vgs_eff - vth;
            let beta = self.kp_t * self.w_over_l;
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
        self.ids_eval = ids_real;

        // The two bulk junctions. Real pn junctions, so `gmin` crosses them —
        // which is what makes this family consistent with the diode and the BJT
        // (see `SimContext::gmin`). Evaluated on every eval, not only in
        // transient: they carry DC current.
        //
        // Polarity-flipped like the channel, so a PMOS's junctions read forward
        // when its bulk is *below* its source.
        let vbs_j = pol * (vb - vs);
        let vbd_j = pol * (vb - vd);
        let (ibs, gbs) =
            Self::bulk_junction(self.isat_bs * self.isat_t_factor, vbs_j, ctx.vt(), ctx.gmin);
        let (ibd, gbd) =
            Self::bulk_junction(self.isat_bd * self.isat_t_factor, vbd_j, ctx.vt(), ctx.gmin);
        self.gbs = gbs;
        self.gbd = gbd;
        // Norton offsets, back in real (unflipped) terms, so `load_residual` can
        // add them without knowing the polarity.
        self.ibs_eq = pol * (ibs - gbs * vbs_j);
        self.ibd_eq = pol * (ibd - gbd * vbd_j);
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
        // The bulk junctions' Norton sources: current out of the bulk and into
        // the internal drain / source.
        for (node, eq) in [(self.source, self.ibs_eq), (self.drain, self.ibd_eq)] {
            if let Some(bk) = self.bulk {
                b[bk] -= eq;
            }
            if let Some(n) = node {
                b[n] += eq;
            }
        }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let (d, g, s, bk) = (self.drain, self.gate, self.source, self.bulk);
        let gm = self.gm;
        let gds = self.gds;
        let gmbs = self.gmbs;
        let gms = gm + gds + gmbs;
        // Series conductances. Zero when the card gives no resistance, and then
        // the internal node aliases the external one so these four cells land on
        // channel cells that already exist — adding zero, which is cheaper than a
        // branch and cannot get the aliasing wrong.
        let g_d = if self.rd > 0.0 { 1.0 / self.rd } else { 0.0 };
        let g_s = if self.rs > 0.0 { 1.0 / self.rs } else { 0.0 };
        // The sixteen values, in the one order `jac_pairs` fixes. Both arms below
        // consume this slice, so the fast path cannot drift from the slow one by
        // reordering: there is only one order.
        let (gbs, gbd) = (self.gbs, self.gbd);
        let vals = [
            gm, gds, -gms, gmbs, -gm, -gds, gms, -gmbs, //
            g_d, -g_d, -g_d, g_d, //
            g_s, -g_s, -g_s, g_s, //
            gbs, -gbs, -gbs, gbs, //
            gbd, -gbd, -gbd, gbd,
        ];

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

        for (&(ri, ci), v) in Self::jac_pairs(d, g, s, bk, self.drain_ext, self.source_ext)
            .iter()
            .zip(vals)
        {
            if let (Some(r), Some(c)) = (ri, ci) {
                mat.a[r][c] += v;
            }
        }
    }

    fn resolve_cells(&mut self, pattern: &Pattern) {
        let pairs = Self::jac_pairs(
            self.drain,
            self.gate,
            self.source,
            self.bulk,
            self.drain_ext,
            self.source_ext,
        );
        let mut cells = [None; 24];
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

    fn noise_sources(&self, ctx: &SimContext, freq: f64) -> Vec<(NodeId, NodeId, f64)> {
        const K_BOLTZMANN: f64 = 1.380649e-23;
        // Channel thermal noise, and flicker if the card asked for it. Both sit
        // drain-to-source and are uncorrelated, so they are one density.
        let thermal = if self.gm.abs() < 1e-18 {
            0.0
        } else {
            8.0 / 3.0 * K_BOLTZMANN * ctx.temperature * self.gm.abs()
        };
        // `KF·|Id|^AF / (f·W·Leff·Cox)`. `Leff` is the drawn `L`: lateral
        // diffusion (`LD`) is on the unmodelled list, so there is nothing to
        // subtract. `validate` has already refused a card with `KF` and no oxide
        // capacitance, so the denominator here is positive whenever `kf > 0`.
        let flicker = if self.kf > 0.0 && freq > 0.0 {
            let norm = self.w * self.l * self.cox;
            self.kf * self.ids_eval.abs().powf(self.af) / (freq * norm)
        } else {
            0.0
        };
        let s_i = thermal + flicker;
        if s_i <= 0.0 {
            return Vec::new();
        }
        vec![(self.drain, self.source, s_i)]
    }

    /// Refuse a card that asks for flicker noise without an oxide capacitance.
    ///
    /// SPICE normalises the flicker density by `W·Leff·Cox`, and `Cox` is zero
    /// unless the card gives `TOX` or `COX`. So `KF` with neither is a division
    /// by zero — which would reach the noise matrix as a non-finite density and
    /// come back as a garbage spectrum rather than as a complaint.
    fn validate(&mut self) -> Result<(), String> {
        // `<=` rather than `!(> 0.0)`: NaN would fall through either way, and a
        // NaN product means the card is broken in a way the message still fits.
        let norm = self.w * self.l * self.cox;
        if self.kf > 0.0 && !norm.is_finite() || (self.kf > 0.0 && norm <= 0.0) {
            return Err(format!(
                "KF={:e} asks for flicker noise, and its density is normalised by \
                 W·L·COX, which is {:e} here. Give the model card a `TOX` (oxide \
                 thickness, m) or a `COX` (capacitance per unit area, F/m²), or \
                 drop `KF`. Dividing by zero would return a non-finite noise \
                 density rather than an error.",
                self.kf, norm
            ));
        }
        Ok(())
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
