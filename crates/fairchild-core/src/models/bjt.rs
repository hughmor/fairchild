//! Gummel-Poon Level 1 bipolar junction transistor (BJT).
//!
//! Implements the Gummel-Poon transport form with the following model
//! parameters: IS, BF, BR, NF, NR, VAF, VAR, IKF, IKR, ISE, NE, ISC, NC, TF,
//! TR, RB, RC, RE, CJE, VJE, MJE, CJC, VJC, MJC, FC — and the instance
//! parameter AREA.  Transit-time diffusion charges (TF·IF/qb, TR·IR) and
//! depletion junction capacitances (CJE, CJC) are both stamped in transient
//! analysis.
//!
//! ## The base charge `qb` is the whole model
//!
//! Everything that makes this Gummel-Poon rather than Ebers-Moll lives in one
//! factor.  Following SPICE3's `BJTload`, because that is the equation set the
//! cards in circulation were extracted against:
//!
//! ```text
//!     q1 = 1 / (1 − VBC/VAF − VBE/VAR)          Early effect
//!     q2 = IF/IKF + IR/IKR                      high-injection knee
//!     qb = q1·(1 + √(1 + 4·q2)) / 2
//!     IC = (IF − IR)/qb − IR/BR − IRC_leak
//! ```
//!
//! Note `q1` is a *reciprocal*: dividing by `qb` multiplies the transport
//! current by `(1 − VBC/VAF − …)`.  Getting that upside down (as this model did
//! until #63) gives a transistor negative output conductance, which is a
//! plausible-looking wrong answer rather than a failure — no golden set `VAF`,
//! so nothing caught it for the model's whole life.
//!
//! `∂qb/∂VBE` and `∂qb/∂VBC` are carried alongside, because the Jacobian's
//! output conductance is entirely made of them: without the `∂qb/∂VBC` term the
//! small-signal `ro` does not match the DC slope of the model's own IC(VCE).
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

use crate::device::{Device, Discretisation, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;
use crate::reactive::ChargeHistory;

/// Seed conductance before the first `eval`, so a device stamped before being
/// evaluated cannot present a singular row. The *operating* `gmin` comes from
/// `SimContext::gmin`.
const GMIN_SEED: f64 = 1e-12;

/// Depletion capacitance at junction voltage `v` (SPICE3 piecewise model).
fn cj_depl(c0: f64, v: f64, vj: f64, mj: f64, fc: f64) -> f64 {
    if c0 == 0.0 {
        return 0.0;
    }
    let fc_vj = fc * vj;
    if v < fc_vj {
        c0 * (1.0 - v / vj).powf(-mj)
    } else {
        let k = (1.0 - fc).powf(1.0 + mj);
        c0 / k * (1.0 - fc * (1.0 + mj) + mj * v / vj)
    }
}

/// Charge integral Q(v) for depletion capacitance (antiderivative of cj_depl).
fn q_depl(c0: f64, v: f64, vj: f64, mj: f64, fc: f64) -> f64 {
    if c0 == 0.0 {
        return 0.0;
    }
    let fc_vj = fc * vj;
    if v < fc_vj {
        let x = 1.0 - v / vj;
        c0 * vj / (1.0 - mj) * (1.0 - x.powf(1.0 - mj))
    } else {
        let x_fc = 1.0 - fc;
        let q_fc = c0 * vj / (1.0 - mj) * (1.0 - x_fc.powf(1.0 - mj));
        let k = x_fc.powf(1.0 + mj);
        let f2 = 1.0 - fc * (1.0 + mj);
        let dv = v - fc_vj;
        q_fc + c0 / k * (f2 * dv + mj / (2.0 * vj) * (v * v - fc_vj * fc_vj))
    }
}

/// The constitutive equations evaluated at one bias point.
///
/// One place computes this and both `eval` and `commit_timestep` read it.  They
/// had drifted: the committed diffusion charge applied no base-charge factor
/// while the evaluated one applied the Early term, so a card with a finite `VAF`
/// integrated a charge that was never evaluated — and it used a different
/// thermal voltage while doing it.
struct Op {
    ic: f64,  // collector current, NPN-equivalent frame
    ib: f64,  // base current
    gf: f64,  // ∂IC/∂VBE
    gce: f64, // −∂IC/∂VBC
    gpi: f64, // ∂IB/∂VBE
    gmu: f64, // ∂IB/∂VBC
    qbe: f64, // B-E diffusion charge TF·IF/qb
    qbc: f64, // B-C diffusion charge TR·IR
    cbe: f64, // ∂QBE/∂VBE
    cbc: f64, // ∂QBC/∂VBC
}

/// Gummel-Poon Level 1 BJT.
pub struct GummelPoonBjt {
    // ── Model parameters ──────────────────────────────────────────────────────
    is: f64,       // transport saturation current (A)
    bf: f64,       // forward beta (current gain)
    br: f64,       // reverse beta
    nf: f64,       // forward emission coefficient
    nr: f64,       // reverse emission coefficient
    vaf: f64,      // forward Early voltage (V); f64::INFINITY = no Early effect
    var: f64,      // reverse Early voltage (V)
    ikf: f64,      // forward high-injection knee current (A); 0 = no roll-off
    ikr: f64,      // reverse high-injection knee current (A)
    ise: f64,      // B-E leakage saturation current (A)
    ne: f64,       // B-E leakage emission coefficient
    isc: f64,      // B-C leakage saturation current (A)
    nc: f64,       // B-C leakage emission coefficient
    tf: f64,       // forward transit time (s) — B-E diffusion charge
    tr: f64,       // reverse transit time (s) — B-C diffusion charge
    rb: f64,       // base ohmic series resistance (Ω); 0 = no internal node
    rc: f64,       // collector ohmic series resistance (Ω)
    re: f64,       // emitter ohmic series resistance (Ω)
    polarity: f64, // +1 NPN, -1 PNP
    vcrit: f64,    // pnjlim critical voltage (derived)
    /// Thermal voltage from the last `setup_model`/`eval`.  `commit_timestep`
    /// gets no `SimContext` and used a hardcoded 0.02585, which is a 0.05 %
    /// error in an exponent — enough to advance a charge the evaluation never
    /// produced.  Cached here so both read the same number.
    vt: f64,
    /// `SimContext::gmin` from the last `eval`, for the same reason `vt` is
    /// stored: `commit_timestep` recomputes `op` from the converged solution and
    /// must do it with the same solve-level constants the eval used.
    gmin: f64,

    // ── Terminal bindings ─────────────────────────────────────────────────────
    // The intrinsic transistor physics operates on the *internal* nodes
    // (collector/base/emitter).  When a series resistance is non-zero, an extra
    // internal node is allocated and the resistor connects it to the matching
    // external terminal (collector_ext/base_ext/emitter_ext).  When the series
    // resistance is zero the internal node simply aliases the external terminal.
    collector: NodeId,
    base: NodeId,
    emitter: NodeId,
    collector_ext: NodeId,
    base_ext: NodeId,
    emitter_ext: NodeId,

    // ── Cached per-NR-iteration quantities (set by eval) ─────────────────────
    vbe_eff: f64, // effective junction voltage B-E (after pnjlim, polarity-corrected)
    vbc_eff: f64, // effective junction voltage B-C
    gf: f64,      // ∂IC/∂VBE_eff (transconductance)
    gce: f64,     // −∂IC/∂VBC_eff — the collector-emitter conductance, and the
    //               only place the Early / high-injection output conductance
    //               enters the matrix
    gpi: f64,     // ∂IB/∂VBE_eff
    gmu: f64,     // ∂IB/∂VBC_eff
    jeq_c: f64,   // Norton offset at collector (see load_residual)
    jeq_b: f64,   // Norton offset at base
    ic_eval: f64, // IC at the last eval, real (polarity-applied) frame
    ib_eval: f64, // IB at the last eval — `.noise` reads these, and used to
    //               read `jeq_b` instead, which is the Norton *offset* and
    //               equals the current only when every node sits at 0 V

    // ── pnjlim history ────────────────────────────────────────────────────────
    vbe_prev: f64,
    vbc_prev: f64,

    // ── Transient charge history ──────────────────────────────────────────────
    qbe_hist: ChargeHistory, // QBE = TF*IF at last accepted timestep
    qbc_hist: ChargeHistory, // QBC = TR*IR at last accepted timestep
    // Current-iterate charges (set by eval when transient flag is set)
    qbe_now: f64,
    qbc_now: f64,
    cbe_eff: f64, // dQBE/dVBE_eff = TF*gf (capacitive conductance)
    cbc_eff: f64, // dQBC/dVBC_eff = TR*gr

    // ── Depletion junction cap model parameters ───────────────────────────────
    cje: f64, // zero-bias B-E depletion capacitance (F)
    vje: f64, // B-E junction potential (V)
    mje: f64, // B-E grading coefficient
    cjc: f64, // zero-bias B-C depletion capacitance (F)
    vjc: f64, // B-C junction potential (V)
    mjc: f64, // B-C grading coefficient
    fc: f64,  // forward-bias cap linearisation coefficient

    /// Instance `AREA`, already folded into the parameters it scales.  Kept only
    /// so the value can be read back; nothing in the physics consults it.
    area: f64,

    // ── Depletion cap transient state ─────────────────────────────────────────
    cje_eval: f64,  // CJE(VBE_eff) at current NR iterate
    cjc_eval: f64,  // CJC(VBC_eff) at current NR iterate
    q_je_eval: f64, // depletion charge at current NR iterate
    q_jc_eval: f64,
    q_je_hist: ChargeHistory, // depletion charge at last committed timestep
    q_jc_hist: ChargeHistory,

    /// The integrator's discretisation, captured during `eval` because
    /// `load_*_tran` receives only `alpha` — which is Backward Euler and
    /// nothing else. `None` outside the transient loop.
    disc: Option<Discretisation>,
}

impl GummelPoonBjt {
    /// Build from model-card parameters.  Returns the device and any
    /// unrecognised parameter names.
    pub fn from_model_params(is_pnp: bool, params: &[(String, f64)]) -> (Self, Vec<String>) {
        let mut is = 1e-16;
        let mut bf = 100.0;
        let mut br = 1.0;
        let mut nf = 1.0;
        let mut nr = 1.0;
        let mut vaf = f64::INFINITY;
        let mut var = f64::INFINITY;
        let mut ikf = 0.0;
        let mut ikr = 0.0;
        let mut ise = 0.0;
        let mut ne = 1.5;
        let mut isc = 0.0;
        let mut nc = 2.0;
        let mut tf = 0.0;
        let mut tr = 0.0;
        let mut rb = 0.0;
        let mut rc = 0.0;
        let mut re = 0.0;
        // Depletion cap params — SPICE3 defaults.
        let mut cje = 0.0;
        let mut vje = 0.75;
        let mut mje = 0.33;
        let mut cjc = 0.0;
        let mut vjc = 0.75;
        let mut mjc = 0.33;
        let mut fc = 0.5;
        let mut unknown = Vec::new();
        for (k, v) in params {
            match k.to_lowercase().as_str() {
                "is" => is = *v,
                "bf" | "hfe" => bf = *v,
                "br" | "hrc" => br = *v,
                "nf" => nf = *v,
                "nr" => nr = *v,
                "vaf" | "va" => vaf = *v,
                "var" | "vb" => var = *v,
                "ikf" | "jbf" => ikf = *v,
                "ikr" | "jbr" => ikr = *v,
                "ise" | "c2" => ise = *v,
                "ne" => ne = *v,
                "isc" | "c4" => isc = *v,
                "nc" => nc = *v,
                "tf" => tf = *v,
                "tr" => tr = *v,
                "rb" => rb = *v,
                "rc" => rc = *v,
                "re" => re = *v,
                "cje" => cje = *v,
                "vje" => vje = *v,
                "mje" => mje = *v,
                "cjc" => cjc = *v,
                "vjc" => vjc = *v,
                "mjc" => mjc = *v,
                "fc" => fc = *v,
                // Accepted and NOT modelled.  The list lives in
                // `crate::unmodelled`, which is what warns about them — keeping
                // a second copy here is how the two drifted apart before.
                k if crate::unmodelled::is_listed(crate::unmodelled::BJT, k) => {}
                _ => unknown.push(k.clone()),
            }
        }
        let dev = GummelPoonBjt {
            is,
            bf,
            br,
            nf,
            nr,
            vaf,
            var,
            ikf,
            ikr,
            ise,
            ne,
            isc,
            nc,
            tf,
            tr,
            rb,
            rc,
            re,
            polarity: if is_pnp { -1.0 } else { 1.0 },
            vcrit: 0.0,
            vt: 0.025864,
            gmin: GMIN_SEED,
            collector: None,
            base: None,
            emitter: None,
            collector_ext: None,
            base_ext: None,
            emitter_ext: None,
            vbe_eff: 0.0,
            vbc_eff: 0.0,
            gf: GMIN_SEED,
            gce: GMIN_SEED,
            gpi: GMIN_SEED / 100.0,
            gmu: GMIN_SEED / 100.0,
            jeq_c: 0.0,
            jeq_b: 0.0,
            ic_eval: 0.0,
            ib_eval: 0.0,
            vbe_prev: 0.0,
            vbc_prev: 0.0,
            qbe_hist: ChargeHistory::default(),
            qbc_hist: ChargeHistory::default(),
            qbe_now: 0.0,
            qbc_now: 0.0,
            cbe_eff: 0.0,
            cbc_eff: 0.0,
            cje,
            vje,
            mje,
            cjc,
            vjc,
            mjc,
            fc,
            area: 1.0,
            cje_eval: 0.0,
            cjc_eval: 0.0,
            q_je_eval: 0.0,
            q_jc_eval: 0.0,
            q_je_hist: ChargeHistory::default(),
            q_jc_hist: ChargeHistory::default(),
            disc: None,
        };
        (dev, unknown)
    }

    /// Apply instance parameters.  `AREA` scales the device: every saturation
    /// current, knee current and junction capacitance is proportional to it, and
    /// every ohmic series resistance inversely proportional — that is what
    /// "N transistors in parallel" means.
    ///
    /// Folded into the fields rather than carried as a multiplier, the way
    /// `Mosfet1::set_instance_params` folds W/L, so nothing downstream has to
    /// remember to apply it.  Called once, before `setup_model`, so `vcrit` sees
    /// the scaled current.
    ///
    /// Returns the names it could not honour; the caller warns about them.
    pub fn set_instance_params(&mut self, params: &[(String, f64)]) -> Vec<String> {
        let mut unknown = Vec::new();
        for (k, v) in params {
            match k.to_lowercase().as_str() {
                "area" => {
                    if *v <= 0.0 {
                        unknown.push(k.clone());
                        continue;
                    }
                    self.area = *v;
                }
                _ => unknown.push(k.clone()),
            }
        }
        let a = self.area;
        if a != 1.0 {
            self.is *= a;
            self.ikf *= a;
            self.ikr *= a;
            self.ise *= a;
            self.isc *= a;
            self.cje *= a;
            self.cjc *= a;
            // A series resistance is a resistance *per device*; N in parallel
            // divide it.  Note this can take a non-zero RB below the threshold
            // `num_extra_nodes` tests — it cannot reach zero from a positive
            // value, so the internal-node count is unchanged.
            self.rb /= a;
            self.rc /= a;
            self.re /= a;
        }
        unknown
    }

    /// The junction currents, their derivatives, and the stored charges at one
    /// bias point — SPICE3 `BJTload`, minus excess phase and the TF bias
    /// modulation (`XTF`/`VTF`/`ITF`), which are not parsed.
    fn op(&self, vbe: f64, vbc: f64, vt: f64, gmin: f64) -> Op {
        // Transport currents.  These are the ones the base charge divides, and
        // the ones the high-injection knee is measured against — not the
        // base currents, which is why IKF is compared with IF and not with IB.
        let nf_vt = self.nf * vt;
        let nr_vt = self.nr * vt;
        let e_be = (vbe / nf_vt).exp();
        let e_bc = (vbc / nr_vt).exp();
        let i_f = self.is * (e_be - 1.0);
        let i_r = self.is * (e_bc - 1.0);
        let gbe = self.is * e_be / nf_vt;
        let gbc = self.is * e_bc / nr_vt;

        // Non-ideal (recombination) leakage: a second, softer exponential in
        // parallel with each junction, and the reason a real beta falls off at
        // LOW current.  It never divides by qb — it is not transport.
        let (i_ben, gben) = if self.ise != 0.0 {
            let ne_vt = self.ne * vt;
            let e = (vbe / ne_vt).exp();
            (self.ise * (e - 1.0), self.ise * e / ne_vt)
        } else {
            (0.0, 0.0)
        };
        let (i_bcn, gbcn) = if self.isc != 0.0 {
            let nc_vt = self.nc * vt;
            let e = (vbc / nc_vt).exp();
            (self.isc * (e - 1.0), self.isc * e / nc_vt)
        } else {
            (0.0, 0.0)
        };

        // ── Base charge ──────────────────────────────────────────────────────
        // `q1` is a reciprocal: dividing the transport current by `qb`
        // MULTIPLIES it by (1 − VBC/VAF − VBE/VAR).  Clamped the way the old
        // code clamped its denominator — a base charge on its way through zero
        // diverges, and 0.1 keeps Newton on the manifold.
        let inv_vaf = if self.vaf.is_finite() {
            1.0 / self.vaf
        } else {
            0.0
        };
        let inv_var = if self.var.is_finite() {
            1.0 / self.var
        } else {
            0.0
        };
        let q1 = 1.0 / (1.0 - inv_vaf * vbc - inv_var * vbe).max(0.1);
        let (qb, dqb_dvbe, dqb_dvbc) = if self.ikf == 0.0 && self.ikr == 0.0 {
            (q1, q1 * q1 * inv_var, q1 * q1 * inv_vaf)
        } else {
            let inv_ikf = if self.ikf > 0.0 { 1.0 / self.ikf } else { 0.0 };
            let inv_ikr = if self.ikr > 0.0 { 1.0 / self.ikr } else { 0.0 };
            let q2 = inv_ikf * i_f + inv_ikr * i_r;
            // √ of a clamped argument: q2 < −1/4 is unphysical (both junctions
            // hard reverse-biased with a knee set), and SPICE clamps rather than
            // producing a NaN that no diagnostic could explain.
            let sqarg = (1.0 + 4.0 * q2).max(0.0).sqrt().max(1.0);
            let qb = q1 * (1.0 + sqarg) / 2.0;
            (
                qb,
                q1 * (qb * inv_var + inv_ikf * gbe / sqarg),
                q1 * (qb * inv_vaf + inv_ikr * gbc / sqarg),
            )
        };

        // ── Terminal currents ────────────────────────────────────────────────
        let it = (i_f - i_r) / qb; // transport, base-charge corrected

        // `gmin·V` on each junction, matching the conductances added to
        // `gpi`/`gmu` below. A conductance that carries no current would only
        // condition the matrix, which is what this used to do.
        let ic = it - i_r / self.br - i_bcn - gmin * vbc;
        let ib = i_f / self.bf + i_ben + i_r / self.br + i_bcn + gmin * vbe + gmin * vbc;

        // ∂IC/∂V: the `it·∂qb/∂V / qb` terms ARE the output conductance.  Drop
        // the ∂qb/∂VBC one and the small-signal ro stops matching the DC slope.
        let gf = (gbe - it * dqb_dvbe) / qb;
        let gce = (gbc + it * dqb_dvbc) / qb + gbc / self.br + gbcn;
        // `gmin` across each junction, at the *terminal* pair.
        //
        // Not folded into `gbe`/`gbc`: those are transport quantities and are
        // divided by `BF`/`BR` on the way out, so a `gmin` added there arrives as
        // `gmin/100` and is not a conductance across anything. It also used to be
        // Jacobian-only — `ic`/`ib` never saw it — so a reverse-biased BJT
        // carried `2·IS = 2e-16` where ngspice carries `gmin·V ≈ 1e-12`.
        //
        // fairchild still reads *half* of ngspice's leakage on a 3-terminal BJT,
        // and that is a different gap: ngspice also puts `gmin` across the
        // collector-substrate junction, whose node defaults to ground. Confirmed
        // by pinning the substrate at the collector potential in ngspice, which
        // removes exactly one `gmin·V`. That junction is not modelled here at all
        // (`docs/model_status.md` — `CJS`/`VJS`/`MJS`/`FCS`), so it is recorded
        // rather than faked.
        let gpi = gbe / self.bf + gben + gmin;
        let gmu = gbc / self.br + gbcn + gmin;

        // Diffusion charge.  Forward carries the base-charge factor (it is the
        // stored minority charge of the transport current); reverse does not,
        // matching SPICE.
        let i_diff = i_f / qb;
        Op {
            ic,
            ib,
            gf,
            gce,
            gpi,
            gmu,
            qbe: self.tf * i_diff,
            qbc: self.tr * i_r,
            cbe: self.tf * (gbe - i_diff * dqb_dvbe) / qb,
            cbc: self.tr * gbc,
        }
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
    fn num_terminals(&self) -> usize {
        4
    } // C B E S (substrate tied to ground by build_devices)

    fn setup_model(&mut self, ctx: &SimContext) {
        let vt = ctx.vt();
        self.vt = vt;
        self.gmin = ctx.gmin;
        self.vcrit = vt * (vt / (std::f64::consts::SQRT_2 * self.is)).ln();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], _ctx: &SimContext) {
        // Terminals order: [C, B, E, S] — substrate (S) is ignored in this implementation.
        debug_assert!(terminals.len() >= 3, "BJT expects [C, B, E, S]");
        self.collector_ext = terminals[0];
        self.base_ext = terminals[1];
        self.emitter_ext = terminals[2];
        // Internal (intrinsic) nodes default to aliasing the external terminals;
        // bind_extra_nodes() re-points them to fresh internal nodes where a series
        // resistance is present.
        self.collector = terminals[0];
        self.base = terminals[1];
        self.emitter = terminals[2];
        // terminals[3] = substrate — tied to ground by caller; not stamped separately.
    }

    /// One internal node per non-zero ohmic series resistance (RB, RC, RE).
    fn num_extra_nodes(&self) -> usize {
        (self.rb > 0.0) as usize + (self.rc > 0.0) as usize + (self.re > 0.0) as usize
    }

    /// Bind internal collector'/base'/emitter' nodes for the series resistances.
    /// Order is fixed (base, collector, emitter) so the assignment is stable.
    fn bind_extra_nodes(&mut self, first_idx: usize) {
        let mut idx = first_idx;
        if self.rb > 0.0 {
            self.base = Some(idx);
            idx += 1;
        }
        if self.rc > 0.0 {
            self.collector = Some(idx);
            idx += 1;
        }
        if self.re > 0.0 {
            self.emitter = Some(idx);
        }
    }

    fn eval(&mut self, x: &[f64], flags: EvalFlags, ctx: &SimContext) {
        let pol = self.polarity;
        let vc = self.collector.map_or(0.0, |i| x[i]);
        let vb = self.base.map_or(0.0, |i| x[i]);
        let ve = self.emitter.map_or(0.0, |i| x[i]);

        let vt = ctx.vt();

        // Polarity-flipped junction voltages (positive in forward active for both NPN and PNP).
        let vbe_raw = pol * (vb - ve);
        let vbc_raw = pol * (vb - vc);

        let vbe_eff = if ctx.jlim_enabled {
            self.pnjlim(vbe_raw, self.vbe_prev, vt)
        } else {
            vbe_raw
        };
        let vbc_eff = if ctx.jlim_enabled {
            self.pnjlim(vbc_raw, self.vbc_prev, vt)
        } else {
            vbc_raw
        };
        self.vbe_prev = vbe_eff;
        self.vbc_prev = vbc_eff;
        self.vbe_eff = vbe_eff;
        self.vbc_eff = vbc_eff;

        self.vt = vt;
        self.gmin = ctx.gmin;
        let op = self.op(vbe_eff, vbc_eff, vt, ctx.gmin);

        // Real currents into each physical terminal: pol * eff.
        // b[C] -= pol*ic; b[B] -= pol*ib; b[E] += pol*(ic + ib)
        // The Norton offset absorbs the linear Jacobian contribution at the
        // eval-point voltages.  Jacobian entries (pol² = 1, so pol-independent):
        //   dIC_real/dVB = gf - gce, dIC_real/dVC = gce, dIC_real/dVE = -gf
        //   dIB_real/dVB = gpi + gmu, dIB_real/dVC = -gmu, dIB_real/dVE = -gpi
        self.gf = op.gf;
        self.gce = op.gce;
        self.gpi = op.gpi;
        self.gmu = op.gmu;
        self.ic_eval = pol * op.ic;
        self.ib_eval = pol * op.ib;

        self.jeq_c = pol * op.ic - (self.gf - self.gce) * vb - self.gce * vc + self.gf * ve;
        self.jeq_b = pol * op.ib - (self.gpi + self.gmu) * vb + self.gmu * vc + self.gpi * ve;

        self.disc = ctx.discretisation;

        if flags.transient {
            self.cbe_eff = op.cbe;
            self.cbc_eff = op.cbc;
            self.qbe_now = op.qbe;
            self.qbc_now = op.qbc;
            self.cje_eval = cj_depl(self.cje, vbe_eff, self.vje, self.mje, self.fc);
            self.cjc_eval = cj_depl(self.cjc, vbc_eff, self.vjc, self.mjc, self.fc);
            self.q_je_eval = q_depl(self.cje, vbe_eff, self.vje, self.mje, self.fc);
            self.q_jc_eval = q_depl(self.cjc, vbc_eff, self.vjc, self.mjc, self.fc);
        } else {
            self.cbe_eff = 0.0;
            self.cbc_eff = 0.0;
            self.qbe_now = 0.0;
            self.qbc_now = 0.0;
            self.cje_eval = 0.0;
            self.cjc_eval = 0.0;
            self.q_je_eval = 0.0;
            self.q_jc_eval = 0.0;
        }
    }

    fn load_residual(&self, b: &mut [f64]) {
        if let Some(c) = self.collector {
            b[c] -= self.jeq_c;
        }
        if let Some(bk) = self.base {
            b[bk] -= self.jeq_b;
        }
        if let Some(e) = self.emitter {
            b[e] += self.jeq_c + self.jeq_b;
        }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let (c, bk, e) = (self.collector, self.base, self.emitter);
        let gf = self.gf;
        let gpi = self.gpi;
        let gmu = self.gmu;
        let gce = self.gce;

        macro_rules! stamp {
            ($ri:expr, $ci:expr, $val:expr) => {
                if let (Some(r), Some(cc)) = ($ri, $ci) {
                    mat.a[r][cc] += $val;
                }
            };
        }

        // Collector row: IC flows C→E; Jacobian of (device draws pol*IC from C).
        stamp!(c, bk, gf - gce);
        stamp!(c, c, gce);
        stamp!(c, e, -gf);

        // Base row: IB flows into B.
        stamp!(bk, bk, gpi + gmu);
        stamp!(bk, c, -gmu);
        stamp!(bk, e, -gpi);

        // Emitter row: IE = IC + IB sources into E (negative sum of above rows).
        stamp!(e, bk, -(gf - gce) - (gpi + gmu));
        stamp!(e, c, -gce - (-gmu));
        stamp!(e, e, gf + gpi);

        // Series ohmic resistances: a conductance 1/R between each external
        // terminal and its internal node.  When R = 0 the internal node aliases
        // the external one and no resistor is stamped.  The `stamp!` macro skips
        // grounded (None) terminals, which correctly yields a conductance-to-
        // ground when an external terminal is ground.
        if self.rb > 0.0 {
            let g = 1.0 / self.rb;
            stamp!(self.base_ext, self.base_ext, g);
            stamp!(self.base_ext, self.base, -g);
            stamp!(self.base, self.base_ext, -g);
            stamp!(self.base, self.base, g);
        }
        if self.rc > 0.0 {
            let g = 1.0 / self.rc;
            stamp!(self.collector_ext, self.collector_ext, g);
            stamp!(self.collector_ext, self.collector, -g);
            stamp!(self.collector, self.collector_ext, -g);
            stamp!(self.collector, self.collector, g);
        }
        if self.re > 0.0 {
            let g = 1.0 / self.re;
            stamp!(self.emitter_ext, self.emitter_ext, g);
            stamp!(self.emitter_ext, self.emitter, -g);
            stamp!(self.emitter, self.emitter_ext, -g);
            stamp!(self.emitter, self.emitter, g);
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], alpha: f64) {
        self.load_residual(b);
        let (disc, pol) = (self.disc, self.polarity);

        // All four charges live in the polarity-flipped ("effective") frame, so
        // the companion is built there and `pol` applies only at stamp time.
        // Companion current flows from the far terminal into the base, the same
        // polarity as the junction it belongs to.
        let mut cap = |hist: &ChargeHistory, far, q_new, cv| {
            let (i_hist, _) = hist.companion(disc, alpha, q_new, cv);
            if let Some(bk) = self.base {
                b[bk] += pol * i_hist;
            }
            if let Some(f) = far {
                b[f] -= pol * i_hist;
            }
        };

        // Transit-time diffusion charge: TF·IF (B-E) and TR·IR (B-C).
        if self.cbe_eff != 0.0 {
            let cv = self.cbe_eff * self.vbe_eff;
            cap(&self.qbe_hist, self.emitter, self.qbe_now, cv);
        }
        if self.cbc_eff != 0.0 {
            let cv = self.cbc_eff * self.vbc_eff;
            cap(&self.qbc_hist, self.collector, self.qbc_now, cv);
        }
        // Depletion charge: CJE (B-E) and CJC (B-C).
        if self.cje_eval != 0.0 {
            let cv = self.cje_eval * self.vbe_eff;
            cap(&self.q_je_hist, self.emitter, self.q_je_eval, cv);
        }
        if self.cjc_eval != 0.0 {
            let cv = self.cjc_eval * self.vbc_eff;
            cap(&self.q_jc_hist, self.collector, self.q_jc_eval, cv);
        }
    }

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, alpha: f64) {
        self.load_jacobian(mat);
        let (c, bk, e) = (self.collector, self.base, self.emitter);
        // Same factor the residual's companions used, from the same place.
        let scale = ChargeHistory::scale(self.disc, alpha);

        macro_rules! stamp {
            ($ri:expr, $ci:expr, $val:expr) => {
                if let (Some(r), Some(cc)) = ($ri, $ci) {
                    mat.a[r][cc] += $val;
                }
            };
        }

        // B-E capacitive companion: cbe_eff between base and emitter.
        if self.cbe_eff != 0.0 {
            let c_be = scale * self.cbe_eff;
            stamp!(bk, bk, c_be);
            stamp!(bk, e, -c_be);
            stamp!(e, bk, -c_be);
            stamp!(e, e, c_be);
        }
        // B-C capacitive companion: cbc_eff between base and collector.
        if self.cbc_eff != 0.0 {
            let c_bc = scale * self.cbc_eff;
            stamp!(bk, bk, c_bc);
            stamp!(bk, c, -c_bc);
            stamp!(c, bk, -c_bc);
            stamp!(c, c, c_bc);
        }
        // B-E depletion cap: CJE
        if self.cje_eval != 0.0 {
            let g_je = scale * self.cje_eval;
            stamp!(bk, bk, g_je);
            stamp!(bk, e, -g_je);
            stamp!(e, bk, -g_je);
            stamp!(e, e, g_je);
        }
        // B-C depletion cap: CJC
        if self.cjc_eval != 0.0 {
            let g_jc = scale * self.cjc_eval;
            stamp!(bk, bk, g_jc);
            stamp!(bk, c, -g_jc);
            stamp!(c, bk, -g_jc);
            stamp!(c, c, g_jc);
        }
    }

    fn commit_timestep(&mut self, x: &[f64]) {
        let pol = self.polarity;
        let vc = self.collector.map_or(0.0, |i| x[i]);
        let vb = self.base.map_or(0.0, |i| x[i]);
        let ve = self.emitter.map_or(0.0, |i| x[i]);
        let vbe_eff = pol * (vb - ve);
        let vbc_eff = pol * (vb - vc);
        self.vbe_prev = vbe_eff;
        self.vbc_prev = vbc_eff;
        // Recomputed analytically from the converged solution rather than reused
        // from the last `eval`, which is one NR iterate behind it — but through
        // the same `op` the eval used, so the charge that gets integrated is the
        // charge the model evaluated.  `vt` comes from the last eval for the
        // same reason.
        let op = self.op(vbe_eff, vbc_eff, self.vt, self.gmin);
        let disc = self.disc;
        self.qbe_hist.advance(disc, op.qbe);
        self.qbc_hist.advance(disc, op.qbc);
        self.q_je_hist
            .advance(disc, q_depl(self.cje, vbe_eff, self.vje, self.mje, self.fc));
        self.q_jc_hist
            .advance(disc, q_depl(self.cjc, vbc_eff, self.vjc, self.mjc, self.fc));
    }

    fn noise_sources(&self, ctx: &SimContext, _freq: f64) -> Vec<(NodeId, NodeId, f64)> {
        // Shot noise on B-E and B-C junctions.
        // i_n_be² = 2q|IB|, flows base→emitter.
        // i_n_ce² = 2q|IC| (collector shot noise), flows collector→emitter.
        const Q_E: f64 = 1.602176634e-19;
        let _ = ctx;
        // The currents from the last eval.  These used to be read off `jeq_c` /
        // `jeq_b`, which are Norton *offsets* — equal to the currents only when
        // every terminal sits at 0 V, so every biased transistor reported a shot
        // noise for a bias point it was not at.
        let mut sources = Vec::new();
        if self.ic_eval.abs() > 1e-20 {
            sources.push((self.collector, self.emitter, 2.0 * Q_E * self.ic_eval.abs()));
        }
        if self.ib_eval.abs() > 1e-20 {
            sources.push((self.base, self.emitter, 2.0 * Q_E * self.ib_eval.abs()));
        }
        sources
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::device::EvalFlags;

    fn ctx() -> SimContext {
        SimContext::default()
    }

    fn npn(bf: f64) -> GummelPoonBjt {
        let (mut q, _) = GummelPoonBjt::from_model_params(
            false,
            &[("is".into(), 1e-15), ("bf".into(), bf), ("br".into(), 1.0)],
        );
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
        let vb = 0.7_f64;
        let vc = 5.0_f64;
        let ve = 0.0_f64;
        let gce = q.gce;
        let ic = q.jeq_c + (q.gf - gce) * vb + gce * vc - q.gf * ve;
        // Expected: IS*exp(VBE/VT); use the same VT as the model (300.15 K default).
        let vt = ctx().vt();
        let ic_expected = 1e-15 * (0.7_f64 / vt).exp();
        assert!(
            (ic - ic_expected).abs() / ic_expected < 0.01,
            "IC={:.4e} expected≈{:.4e}",
            ic,
            ic_expected
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

        let vb = 0.65_f64;
        let vc = 5.0_f64;
        let ve = 0.0_f64;
        let gce = q.gce;
        let ic = q.jeq_c + (q.gf - gce) * vb + gce * vc - q.gf * ve;
        let ib = q.jeq_b + (q.gpi + q.gmu) * vb - q.gmu * vc - q.gpi * ve;

        // In forward active (VBC << 0), IC/IB ≈ BF.
        let beta_measured = ic / ib;
        assert!(
            (beta_measured - bf).abs() / bf < 0.01,
            "β = {:.1} expected {:.1}",
            beta_measured,
            bf
        );
    }

    #[test]
    fn cje_cjc_stamps_in_jacobian_tran() {
        // Verify that CJE and CJC add companion conductances to the Jacobian
        // during transient analysis and are zero for DC.
        let (mut q, _) = GummelPoonBjt::from_model_params(
            false,
            &[
                ("is".into(), 1e-15),
                ("bf".into(), 100.0),
                ("cje".into(), 2e-12),
                ("cjc".into(), 1e-12),
            ],
        );
        q.setup_model(&ctx());
        // C=0, B=1, E=2
        q.setup_instance(&[Some(0), Some(1), Some(2), None], &ctx());
        q.vbe_prev = 0.65;
        q.vbc_prev = -4.35;
        let x = [5.0_f64, 0.65, 0.0];

        // DC eval: depletion caps must be zero in Jacobian.
        q.eval(&x, EvalFlags::dc(), &ctx());
        assert_eq!(q.cje_eval, 0.0, "cje_eval must be zero in DC");
        assert_eq!(q.cjc_eval, 0.0, "cjc_eval must be zero in DC");

        // Transient eval: caps should be non-zero.
        q.eval(&x, EvalFlags::tran(), &ctx());
        assert!(q.cje_eval > 0.0, "cje_eval should be positive in transient");
        assert!(q.cjc_eval > 0.0, "cjc_eval should be positive in transient");

        // The Jacobian diagonal at B (node 1) should include alpha * (cbe_eff + cbc_eff + cje + cjc).
        let mut mat = MnaMatrix::zeros(3);
        let alpha = 1.0 / 1e-9; // 1/h for 1ns step
        q.load_jacobian_tran(&mut mat, alpha);
        // J[B][B] must include contributions from both depletion caps (and transit-time caps).
        // At minimum, alpha*cje_eval + alpha*cjc_eval must be present.
        let min_expected = alpha * (q.cje_eval + q.cjc_eval);
        assert!(
            mat.a[1][1] >= min_expected,
            "J[B][B]={:.3e} < expected depletion contribution {:.3e}",
            mat.a[1][1],
            min_expected
        );
        // Off-diagonal B-E must be negative (conductance from B to E).
        assert!(mat.a[1][2] < 0.0, "J[B][E] should be negative from CJE");
        // Off-diagonal B-C must be negative.
        assert!(mat.a[1][0] < 0.0, "J[B][C] should be negative from CJC");
    }

    #[test]
    fn series_resistances_allocate_internal_nodes_and_stamp() {
        // RB and RC present, RE absent → two internal nodes; the intrinsic
        // collector/base move to fresh internal nodes while the emitter aliases
        // its external terminal.  The external terminals couple to the internal
        // nodes only through the series conductances 1/RB and 1/RC.
        let (mut q, unknown) = GummelPoonBjt::from_model_params(
            false,
            &[
                ("is".into(), 1e-15),
                ("bf".into(), 100.0),
                ("rb".into(), 1000.0),
                ("rc".into(), 50.0),
            ],
        );
        assert!(unknown.is_empty(), "rb/rc must be recognised: {unknown:?}");
        q.setup_model(&ctx());
        // External: C=0, B=1, E=2.
        q.setup_instance(&[Some(0), Some(1), Some(2), None], &ctx());
        assert_eq!(
            q.num_extra_nodes(),
            2,
            "RB and RC each need an internal node"
        );
        // Allocate internal nodes starting at index 3.
        q.bind_extra_nodes(3);
        assert_eq!(q.base, Some(3), "internal base = first extra node");
        assert_eq!(
            q.collector,
            Some(4),
            "internal collector = second extra node"
        );
        assert_eq!(q.emitter, Some(2), "emitter aliases external (RE = 0)");

        q.vbe_prev = 0.7;
        q.vbc_prev = -4.3;
        // x indexed by node: ext C,B,E then int B',C'.
        let x = [5.0_f64, 0.7, 0.0, 0.7, 5.0];
        q.eval(&x, EvalFlags::dc(), &ctx());

        let mut mat = MnaMatrix::zeros(5);
        q.load_jacobian(&mut mat);
        // External base row sees only the RB conductance (intrinsic physics is on
        // the internal base node 3).
        assert!(
            (mat.a[1][1] - 1.0 / 1000.0).abs() < 1e-12,
            "ext-base diagonal should equal 1/RB, got {:.6e}",
            mat.a[1][1]
        );
        assert!(
            (mat.a[1][3] + 1.0 / 1000.0).abs() < 1e-12,
            "ext-base↔int-base coupling should be -1/RB, got {:.6e}",
            mat.a[1][3]
        );
        // External collector row sees only the RC conductance.
        assert!(
            (mat.a[0][0] - 1.0 / 50.0).abs() < 1e-12,
            "ext-collector diagonal should equal 1/RC, got {:.6e}",
            mat.a[0][0]
        );
        assert!(
            (mat.a[0][4] + 1.0 / 50.0).abs() < 1e-12,
            "ext-collector↔int-collector coupling should be -1/RC, got {:.6e}",
            mat.a[0][4]
        );
    }

    #[test]
    fn pnp_active_currents_flow_correctly() {
        // PNP in forward active: VEB=0.7V (VE=5, VB=4.3), VCB>0 (VC=0, VB=4.3)
        // IC should be negative (flows OUT of collector = INTO collector node from circuit side)
        let (mut q, _) =
            GummelPoonBjt::from_model_params(true, &[("is".into(), 1e-15), ("bf".into(), 100.0)]);
        q.setup_model(&ctx());
        q.setup_instance(&[Some(0), Some(1), Some(2), None], &ctx());
        q.vbe_prev = -0.7; // raw (VB-VE) for PNP
        q.vbc_prev = 4.3;
        // VC=0, VB=4.3, VE=5.0  → VB-VE = -0.7, VB-VC = 4.3
        let x = [0.0_f64, 4.3, 5.0];
        q.eval(&x, EvalFlags::dc(), &ctx());

        let vb = 4.3_f64;
        let vc = 0.0_f64;
        let ve = 5.0_f64;
        let gce = q.gce;
        let ic_into_c = q.jeq_c + (q.gf - gce) * vb + gce * vc - q.gf * ve;
        // PNP: current flows from emitter to collector internally.
        // The collector node SOURCE (external) sees current flowing INTO it
        // when collector is the low-voltage terminal — IC_into_C is negative.
        assert!(
            ic_into_c < 0.0,
            "PNP IC_into_C should be negative, got {:.4e}",
            ic_into_c
        );
        // |IC| / |IB| ≈ BF = 100
        let ib_into_b = q.jeq_b + (q.gpi + q.gmu) * vb - q.gmu * vc - q.gpi * ve;
        let beta = ic_into_c.abs() / ib_into_b.abs();
        assert!((beta - 100.0).abs() / 100.0 < 0.02, "PNP β={:.1}", beta);
    }
}
