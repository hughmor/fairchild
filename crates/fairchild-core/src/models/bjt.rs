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

/// The substrate junction's capacitance at `v`, polarity-flipped so forward is
/// positive.
///
/// Not `cj_depl`, and the difference is measured rather than assumed. `cj_depl`
/// takes an `FC` and switches to a straight line at `FC·VJ`. **`FCS` is inert in
/// ngspice** — the capacitance at 0.5 V forward is bit-identical for `FCS` of 0.1,
/// 0.5, 0.9 and absent — and the forward branch is a linearisation about *zero*
/// bias instead:
///
/// ```text
/// v <= 0   CJS·(1 − v/VJS)^−MJS      the depletion law, matched to 5e-8
/// v >  0   CJS·(1 + MJS·v/VJS)       matched to 1.7e-7, and past VJS
/// ```
///
/// The forward form holds out to 2 V, well past `VJS`, where the depletion law is
/// singular. That is what the linearisation is for.
fn cjs_cap(cjs: f64, v: f64, vjs: f64, mjs: f64) -> f64 {
    if cjs <= 0.0 || vjs <= 0.0 {
        return 0.0;
    }
    if v > 0.0 {
        cjs * (1.0 + mjs * v / vjs)
    } else {
        cjs * (1.0 - v / vjs).powf(-mjs)
    }
}

/// The charge whose derivative is [`cjs_cap`], zero at zero bias.
///
/// ```text
/// v <= 0   −CJS·VJS/(1 − MJS)·((1 − v/VJS)^(1 − MJS) − 1)
/// v >  0   CJS·(v + MJS·v²/(2·VJS))
/// ```
///
/// Integrated in closed form rather than accumulated, so a transient's charge
/// cannot drift from the capacitance the same eval reported. `MJS = 1` makes the
/// reverse antiderivative a logarithm, which is a real card value for an abrupt
/// junction, so it is handled rather than divided by zero.
fn q_js_charge(cjs: f64, v: f64, vjs: f64, mjs: f64) -> f64 {
    if cjs <= 0.0 || vjs <= 0.0 {
        return 0.0;
    }
    if v > 0.0 {
        cjs * (v + mjs * v * v / (2.0 * vjs))
    } else if (mjs - 1.0).abs() < 1e-12 {
        -cjs * vjs * (1.0 - v / vjs).ln()
    } else {
        // Spelled the way `q_depl` spells the same antiderivative, so the two can
        // be read against each other.
        cjs * vjs / (1.0 - mjs) * (1.0 - (1.0 - v / vjs).powf(1.0 - mjs))
    }
}

/// The base resistance at one bias point, falling from `rb` towards `rbm`.
///
/// Two laws, and which one applies is decided by `IRB` alone:
///
/// ```text
/// IRB == 0   rbm + (rb − rbm)/qb
/// IRB >  0   rbm + 3·(rb − rbm)·(tan z − z)/(z·tan²z)
///            z = (sqrt(1 + 144/pi²·ib/IRB) − 1) / ((24/pi²)·sqrt(ib/IRB))
/// ```
///
/// `rbm` defaults to `rb`, which makes both laws return `rb` exactly, so a card
/// that does not ask for a variable base resistance does not get one. Measured:
/// ngspice gives bit-identical currents for `RB=10k` and `RB=10k RBM=10k`.
///
/// # The limits are where this needs care
///
/// As `ib → 0`, `z → 0` and the `tan` expression is `0/0`. Its limit is `1/3`, so
/// the whole factor tends to 1 and `rb_eff → rb`. Series-expanded below that
/// threshold rather than evaluated, because the cancellation loses every digit.
///
/// As `ib → ∞`, `z → pi/2` and `tan z → ∞`. Written as
/// `1/(z·tan z) − 1/tan²z` so both terms simply underflow to zero and
/// `rb_eff → rbm`, instead of forming `∞/∞`.
fn base_resistance(rb: f64, rbm: f64, irb: f64, ib: f64, qb: f64) -> f64 {
    if rb <= 0.0 || rbm >= rb {
        return rb;
    }
    if irb <= 0.0 {
        // `qb` is at least 1 for any physical bias, so this only ever lowers it.
        return if qb > 0.0 { rbm + (rb - rbm) / qb } else { rb };
    }
    if ib <= 0.0 {
        return rb;
    }
    let x = ib / irb;
    let sx = x.sqrt();
    let z = ((1.0 + 144.0 / (std::f64::consts::PI * std::f64::consts::PI) * x).sqrt() - 1.0)
        / ((24.0 / (std::f64::consts::PI * std::f64::consts::PI)) * sx);
    if !z.is_finite() || z <= 0.0 {
        return rb;
    }
    // Below this the `tan z − z` cancellation has no digits left; the series limit
    // of the bracket is `1/3 + 2z²/15 + …`, and the leading term is enough here
    // because `z < 1e-4` makes the correction 1e-9 relative.
    const Z_SMALL: f64 = 1e-4;
    let bracket = if z < Z_SMALL {
        1.0 / 3.0
    } else {
        let t = z.tan();
        if !t.is_finite() || t == 0.0 {
            return rbm;
        }
        1.0 / (z * t) - 1.0 / (t * t)
    };
    rbm + 3.0 * (rb - rbm) * bracket
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
    qbe: f64, // B-E diffusion charge TF·(1+xf)·IF/qb
    qbc: f64, // B-C diffusion charge TR·IR
    cbe: f64, // ∂QBE/∂VBE
    /// `dQBE/dVBC_eff` — a **transcapacitance**. The base-emitter diffusion
    /// charge depends on the base-collector voltage twice over: through `VTF`'s
    /// modulation of the transit time, and through the base charge `qb`, which
    /// `VAF` and `IKF` make bias-dependent. Zero for a card that sets none of
    /// those, which is why nothing needed it before.
    cbe_x: f64,
    cbc: f64, // ∂QBC/∂VBC
    /// The base charge factor, which `base_resistance` needs when `IRB` is absent.
    qb: f64,
}

/// Gummel-Poon Level 1 BJT.
pub struct GummelPoonBjt {
    // ── Model parameters ──────────────────────────────────────────────────────
    is: f64,  // transport saturation current (A)
    bf: f64,  // forward beta (current gain)
    br: f64,  // reverse beta
    nf: f64,  // forward emission coefficient
    nr: f64,  // reverse emission coefficient
    vaf: f64, // forward Early voltage (V); f64::INFINITY = no Early effect
    var: f64, // reverse Early voltage (V)
    ikf: f64, // forward high-injection knee current (A); 0 = no roll-off
    ikr: f64, // reverse high-injection knee current (A)
    ise: f64, // B-E leakage saturation current (A)
    ne: f64,  // B-E leakage emission coefficient
    isc: f64, // B-C leakage saturation current (A)
    nc: f64,  // B-C leakage emission coefficient
    tf: f64,  // forward transit time (s) — B-E diffusion charge
    xtf: f64, // TF bias-modulation coefficient
    vtf: f64, // the VBC scale in TF's modulation; 0 disables that term
    itf: f64, // the high-current knee in TF's modulation; 0 disables it
    tr: f64,  // reverse transit time (s) — B-C diffusion charge
    rb: f64,  // base ohmic series resistance (Ω); 0 = no internal node
    /// The floor `rb` falls towards at high base current. Defaults to `rb`, which
    /// disables the variation.
    rbm: f64,
    /// The base current at which `rb` has fallen halfway. `0` selects the
    /// `qb`-driven law instead of the `tan z` one.
    irb: f64,
    /// `rb` at the last `eval`, from [`base_resistance`].
    ///
    /// Not lagged state: it is a function of the iterate's own node voltages, so
    /// once the convergence test says `x` has stopped moving this has too. The
    /// Jacobian stamps `1/rb_eff` and does not differentiate it, which costs
    /// Newton steps and not correctness — the residual defines the answer and
    /// both read the same value.
    rb_eff: f64,
    rc: f64,       // collector ohmic series resistance (Ω)
    re: f64,       // emitter ohmic series resistance (Ω)
    polarity: f64, // +1 NPN, -1 PNP
    vcrit: f64,    // pnjlim critical voltage (derived)
    /// Thermal voltage from the last `setup_model`/`eval`.  `commit_timestep`
    /// gets no `SimContext` and used a hardcoded 0.02585, which is a 0.05 %
    /// error in an exponent — enough to advance a charge the evaluation never
    /// produced.  Cached here so both read the same number.
    vt: f64,
    /// `KF`/`AF` — flicker noise coefficient and exponent.
    kf: f64,
    af: f64,
    /// `TNOM` — the temperature this card's parameters were extracted at.
    tnom: f64,
    /// `EG` — activation energy, eV.
    eg: f64,
    /// `XTI` — the saturation current's temperature exponent.
    xti: f64,
    /// `XTB` — the betas' temperature exponent.
    xtb: f64,
    /// The two junction potentials and zero-bias capacitance factors at the
    /// operating temperature.
    ///
    /// Nominal / 1.0 until `setup_model` runs, which is before any eval. Held
    /// rather than recomputed: each potential costs two logs and an exponential
    /// and depends on nothing that moves inside a solve.
    vje_t: f64,
    vjc_t: f64,
    cje_t_factor: f64,
    cjc_t_factor: f64,
    /// `IS(T)/IS` and `BF(T)/BF`, from `crate::temperature`.
    ///
    /// Factors rather than scaled parameters, so `setup_model` running twice
    /// cannot apply them twice, and so `AREA` still multiplies on top.
    is_t_factor: f64,
    beta_t_factor: f64,
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
    /// The substrate. `Q c b e model` leaves it at ground, `Q c b e s model`
    /// binds it. The junction hangs off the *internal* collector, so a card with
    /// `RC` puts the series resistance between the substrate junction and the
    /// collector pin, which is where SPICE puts it.
    substrate: NodeId,

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
    cbe_eff: f64,   // dQBE/dVBE_eff (capacitive conductance)
    cbe_x_eff: f64, // dQBE/dVBC_eff, the transcapacitance
    cbc_eff: f64,   // dQBC/dVBC_eff = TR*gr

    // ── Depletion junction cap model parameters ───────────────────────────────
    cje: f64, // zero-bias B-E depletion capacitance (F)
    vje: f64, // B-E junction potential (V)
    mje: f64, // B-E grading coefficient
    cjc: f64, // zero-bias B-C depletion capacitance (F)
    vjc: f64, // B-C junction potential (V)
    mjc: f64, // B-C grading coefficient
    fc: f64,  // forward-bias cap linearisation coefficient
    /// The fraction of `CJC` that connects to the **internal** base node, so `RB`
    /// sits in series with it. The rest hangs off the external base pin. Default
    /// 1.0, which is SPICE's and which puts all of it inside `RB`.
    xcjc: f64,
    /// The two halves of the split at this iterate, and the external half's
    /// charge history. The internal half keeps `q_jc_hist`.
    cjcx_eval: f64,
    q_jcx_eval: f64,
    q_jcx_hist: ChargeHistory,
    /// `pol·(V(base_ext) − V(collector))`, the external half's junction voltage.
    vbcx_eff: f64,

    // ── Collector-substrate junction ─────────────────────────────────────────
    iss: f64, // substrate saturation current (A) — default 0, no DC branch
    cjs: f64, // zero-bias collector-substrate capacitance (F) — default 0
    vjs: f64, // substrate junction potential (V) — default 0.75
    mjs: f64, // substrate grading coefficient — default 0, a constant cap
    /// Conductance and Norton offset of the substrate junction at this iterate.
    gsub: f64,
    isub_eq: f64,
    /// Substrate junction voltage at this iterate, polarity-flipped so forward
    /// is positive for both an NPN and a PNP.
    vsub_eff: f64,
    /// Substrate depletion capacitance and charge at this iterate, and the charge
    /// at the last committed timestep.
    cjs_eval: f64,
    q_js_eval: f64,
    q_js_hist: ChargeHistory,

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
        let mut xcjc = 1.0_f64;
        let mut rbm = 0.0_f64;
        let mut irb = 0.0_f64;
        let mut xtf = 0.0_f64;
        let mut vtf = 0.0_f64;
        let mut itf = 0.0_f64;
        // Substrate junction. ngspice's defaults, measured: no capacitance and no
        // DC branch unless the card asks, a 0.75 V potential, and a grading
        // coefficient of zero, which makes `CJS` alone a constant capacitance.
        let mut iss = 0.0_f64;
        let mut cjs = 0.0_f64;
        let mut vjs = 0.75_f64;
        let mut mjs = 0.0_f64;
        let mut kf = 0.0_f64;
        let mut af = 1.0_f64;
        let mut tnom_c = crate::temperature::TNOM_DEFAULT_K - 273.15;
        let mut eg = crate::temperature::EG_DEFAULT;
        let mut xti = crate::temperature::XTI_DEFAULT;
        let mut xtb = crate::temperature::XTB_DEFAULT;
        let mut unknown = Vec::new();
        for (k, v) in params {
            match k.to_lowercase().as_str() {
                "is" => is = *v,
                // Degrees Celsius on the card, like `.temp`.
                "kf" => kf = *v,
                "af" => af = *v,
                "tnom" => tnom_c = *v,
                "eg" => eg = *v,
                "xti" => xti = *v,
                "xtb" => xtb = *v,
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
                "xtf" => xtf = *v,
                "vtf" => vtf = *v,
                "itf" | "jtf" => itf = *v,
                "tr" => tr = *v,
                "rb" => rb = *v,
                "rbm" => rbm = *v,
                "irb" | "jrb" | "irb0" => irb = *v,
                "rc" => rc = *v,
                "re" => re = *v,
                "cje" => cje = *v,
                "vje" => vje = *v,
                "mje" => mje = *v,
                "cjc" => cjc = *v,
                "vjc" => vjc = *v,
                "mjc" => mjc = *v,
                "fc" => fc = *v,
                "xcjc" | "cdis" => xcjc = *v,
                "iss" => iss = *v,
                "cjs" | "ccs" => cjs = *v,
                "vjs" | "pjs" => vjs = *v,
                "mjs" | "ms" => mjs = *v,
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
            xtf,
            vtf,
            itf,
            tr,
            rb,
            // `RBM` defaults to `RB`, and `0` on the card means "not given". A
            // literal `RBM=0` would be a zero-resistance floor, which SPICE also
            // reads as absent, so the two are indistinguishable and match.
            rbm: if rbm > 0.0 { rbm } else { rb },
            irb,
            rb_eff: rb,
            rc,
            re,
            polarity: if is_pnp { -1.0 } else { 1.0 },
            vcrit: 0.0,
            vt: 0.025864,
            vje_t: 0.75,
            vjc_t: 0.75,
            cje_t_factor: 1.0,
            cjc_t_factor: 1.0,
            is_t_factor: 1.0,
            beta_t_factor: 1.0,
            gmin: GMIN_SEED,
            collector: None,
            base: None,
            emitter: None,
            collector_ext: None,
            base_ext: None,
            emitter_ext: None,
            substrate: None,
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
            cbe_x_eff: 0.0,
            cbc_eff: 0.0,
            cje,
            vje,
            mje,
            cjc,
            vjc,
            mjc,
            fc,
            xcjc,
            cjcx_eval: 0.0,
            q_jcx_eval: 0.0,
            q_jcx_hist: ChargeHistory::default(),
            vbcx_eff: 0.0,
            iss,
            cjs,
            vjs,
            mjs,
            gsub: GMIN_SEED,
            isub_eq: 0.0,
            vsub_eff: 0.0,
            cjs_eval: 0.0,
            q_js_eval: 0.0,
            q_js_hist: ChargeHistory::default(),
            kf,
            af,
            tnom: tnom_c + 273.15,
            eg,
            xti,
            xtb,
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
            // `RBM` is a resistance per device like `RB`, and `IRB` is a current,
            // so N in parallel divide the first and multiply the second.
            self.rbm /= a;
            self.irb *= a;
            self.rb_eff = self.rb;
            self.rc /= a;
            self.re /= a;
        }
        unknown
    }

    /// The junction currents, their derivatives, and the stored charges at one
    /// bias point — SPICE3 `BJTload`, minus excess phase and the TF bias
    /// modulation (`XTF`/`VTF`/`ITF`), which are not parsed.
    fn op(&self, vbe: f64, vbc: f64, vt: f64, gmin: f64) -> Op {
        // Temperature-scaled, once, at the top: `op` reads these ten times and a
        // scaling applied at each read is ten chances to miss one. The factors
        // come from `setup_model`; see `crate::temperature`.
        let is_t = self.is * self.is_t_factor;
        let bf_t = self.bf * self.beta_t_factor;
        let br_t = self.br * self.beta_t_factor;
        // Transport currents.  These are the ones the base charge divides, and
        // the ones the high-injection knee is measured against — not the
        // base currents, which is why IKF is compared with IF and not with IB.
        let nf_vt = self.nf * vt;
        let nr_vt = self.nr * vt;
        let e_be = (vbe / nf_vt).exp();
        let e_bc = (vbc / nr_vt).exp();
        let i_f = is_t * (e_be - 1.0);
        let i_r = is_t * (e_bc - 1.0);
        let gbe = is_t * e_be / nf_vt;
        let gbc = is_t * e_bc / nr_vt;

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
        let ic = it - i_r / br_t - i_bcn - gmin * vbc;
        let ib = i_f / bf_t + i_ben + i_r / br_t + i_bcn + gmin * vbe + gmin * vbc;

        // ∂IC/∂V: the `it·∂qb/∂V / qb` terms ARE the output conductance.  Drop
        // the ∂qb/∂VBC one and the small-signal ro stops matching the DC slope.
        let gf = (gbe - it * dqb_dvbe) / qb;
        let gce = (gbc + it * dqb_dvbc) / qb + gbc / br_t + gbcn;
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
        let gpi = gbe / bf_t + gben + gmin;
        let gmu = gbc / br_t + gbcn + gmin;

        // Diffusion charge.  Forward carries the base-charge factor (it is the
        // stored minority charge of the transport current); reverse does not,
        // matching SPICE.
        let i_diff = i_f / qb;

        // TF's bias modulation. `xf` is the *excess* factor: `TF_eff =
        // TF·(1 + xf)`. `XTF = 0` leaves both terms out, which is the default and
        // every card that does not ask.
        //
        // `ITF = 0` and `VTF = 0` each disable their own factor rather than
        // dividing by zero, which is SPICE's convention for both and is why they
        // default to zero rather than to infinity.
        let (xf, tmp) = if self.xtf != 0.0 && i_f > 0.0 {
            let tmp = if self.itf > 0.0 {
                i_f / (i_f + self.itf)
            } else {
                1.0
            };
            let vbc_term = if self.vtf > 0.0 {
                (vbc / (1.44 * self.vtf)).exp()
            } else {
                1.0
            };
            (self.xtf * tmp * tmp * vbc_term, tmp)
        } else {
            (0.0, 1.0)
        };

        // The charge takes `(1 + xf)`; the capacitance takes
        // `(1 + xf·(3 − 2·tmp))`, because `IF·d(xf)/d(vbe) = 2·xf·(1−tmp)·gbe`.
        // Two different factors from one law, which is why the derivative is
        // written out rather than formed as `TF_eff·gbe`.
        let q_factor = 1.0 + xf;
        let c_factor = 1.0 + xf * (3.0 - 2.0 * tmp);
        let cbe = self.tf * (gbe * c_factor - q_factor * i_diff * dqb_dvbe) / qb;
        // The transcapacitance. `VTF`'s term is `xf/(1.44·VTF)`; the `qb` term is
        // there whenever `VAF` or `IKF` is finite and was missing before.
        let dxf_dvbc = if self.vtf > 0.0 {
            xf / (1.44 * self.vtf)
        } else {
            0.0
        };
        let cbe_x = self.tf * i_diff * (dxf_dvbc - q_factor * dqb_dvbc / qb);

        Op {
            ic,
            ib,
            gf,
            gce,
            gpi,
            gmu,
            qbe: self.tf * q_factor * i_diff,
            qbc: self.tr * i_r,
            cbe,
            cbe_x,
            cbc: self.tr * gbc,
            qb,
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
        // Idempotent factors, so a second `setup_model` cannot double-apply and
        // `AREA` — which arrives afterwards — still multiplies on top.
        self.is_t_factor =
            crate::temperature::bjt_is_factor(ctx.temperature, self.tnom, self.eg, self.xti);
        self.beta_t_factor = crate::temperature::beta_factor(ctx.temperature, self.tnom, self.xtb);
        // Both junctions' potentials and zero-bias capacitances. Idempotent, and
        // derived from the nominal values rather than the previous result.
        let (t, tnom) = (ctx.temperature, self.tnom);
        self.vje_t = crate::temperature::scaled_junction_potential(self.vje, t, tnom);
        self.vjc_t = crate::temperature::scaled_junction_potential(self.vjc, t, tnom);
        self.cje_t_factor = crate::temperature::junction_cap_factor(self.vje, self.mje, t, tnom);
        self.cjc_t_factor = crate::temperature::junction_cap_factor(self.vjc, self.mjc, t, tnom);
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
        // The substrate. `Q c b e model` leaves this at ground, `Q c b e s model`
        // binds the named net. It used to be dropped, which cost one `gmin·V` of
        // leakage against ngspice on every reverse-biased BJT.
        self.substrate = terminals.get(3).copied().flatten();
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
        // The base resistance at this iterate. A function of the iterate's own
        // node voltages through `ib` and `qb`, so it is not lagged state and the
        // convergence test sees it stop moving when `x` does.
        self.rb_eff = base_resistance(self.rb, self.rbm, self.irb, op.ib, op.qb);

        self.gf = op.gf;
        self.gce = op.gce;
        self.gpi = op.gpi;
        self.gmu = op.gmu;
        self.ic_eval = pol * op.ic;
        self.ib_eval = pol * op.ib;

        self.jeq_c = pol * op.ic - (self.gf - self.gce) * vb - self.gce * vc + self.gf * ve;
        self.jeq_b = pol * op.ib - (self.gpi + self.gmu) * vb + self.gmu * vc + self.gpi * ve;

        // ── Collector-substrate junction ─────────────────────────────────
        //
        // Real, and unconditional. ngspice puts one `gmin` across it whether or
        // not the card gives a `CJS` or an `ISS` — measured, and it is the whole
        // reason a reverse-biased BJT here used to read `1·gmin·V` against
        // ngspice's `2·gmin·V`.
        //
        // Polarity-flipped like the other two junctions: for an NPN the substrate
        // is p and the collector n, so a substrate *above* the collector is
        // forward. `pol` inverts that for a PNP.
        //
        // Plain Shockley, not the flat reverse branch MOS1 uses. Measured: at
        // −0.05 V with `ISS = 1e-15` ngspice reads 8.553040e-16 against
        // Shockley's 8.553119e-16, where a flat `−ISS` would read 1e-15.
        let vsub = pol * (self.substrate.map_or(0.0, |i| x[i]) - vc);
        self.vsub_eff = vsub;
        let (isub, gsub) = if self.iss > 0.0 {
            let e = (vsub / vt).exp();
            (
                self.iss * (e - 1.0) + ctx.gmin * vsub,
                self.iss * e / vt + ctx.gmin,
            )
        } else {
            (ctx.gmin * vsub, ctx.gmin)
        };
        self.gsub = gsub;
        // Norton offset back in the real frame, so `load_residual` needs no
        // polarity of its own.
        self.isub_eq = pol * (isub - gsub * vsub);

        self.disc = ctx.discretisation;

        if flags.transient {
            self.cbe_eff = op.cbe;
            self.cbe_x_eff = op.cbe_x;
            self.cbc_eff = op.cbc;
            self.qbe_now = op.qbe;
            self.qbc_now = op.qbc;
            self.cje_eval = cj_depl(
                self.cje * self.cje_t_factor,
                vbe_eff,
                self.vje_t,
                self.mje,
                self.fc,
            );
            self.cjc_eval = cj_depl(
                self.cjc * self.cjc_t_factor,
                vbc_eff,
                self.vjc_t,
                self.mjc,
                self.fc,
            );
            self.q_je_eval = q_depl(
                self.cje * self.cje_t_factor,
                vbe_eff,
                self.vje_t,
                self.mje,
                self.fc,
            );
            self.q_jc_eval = q_depl(
                self.cjc * self.cjc_t_factor,
                vbc_eff,
                self.vjc_t,
                self.mjc,
                self.fc,
            );
            // No temperature factor on `CJS`: `TNOM` moves `CJE` and `CJC` here
            // through `cje_t_factor`/`cjc_t_factor`, and the substrate junction
            // is left at its nominal value because nothing has measured the law
            // for it. Recorded in `docs/model_status.md` rather than guessed.
            self.cjs_eval = cjs_cap(self.cjs, vsub, self.vjs, self.mjs);
            self.q_js_eval = q_js_charge(self.cjs, vsub, self.vjs, self.mjs);

            // `XCJC` splits the base-collector depletion capacitance across the
            // base resistance: the internal fraction keeps `cjc_eval` and sees
            // `RB` in series, the rest hangs off the external base pin where it
            // does not. Both halves are evaluated at the *internal* junction
            // voltage, which is the split SPICE makes and which the ngspice
            // agreement at five values of `XCJC` confirms.
            //
            // With `RB = 0` the two base nodes alias, so the split is invisible
            // and a card without a base resistance is unaffected.
            let xcjc = self.xcjc.clamp(0.0, 1.0);
            self.cjcx_eval = (1.0 - xcjc) * self.cjc_eval;
            self.q_jcx_eval = (1.0 - xcjc) * self.q_jc_eval;
            self.cjc_eval *= xcjc;
            self.q_jc_eval *= xcjc;
            self.vbcx_eff = pol * (self.base_ext.map_or(0.0, |i| x[i]) - vc);
        } else {
            self.cbe_eff = 0.0;
            self.cbe_x_eff = 0.0;
            self.cbc_eff = 0.0;
            self.qbe_now = 0.0;
            self.qbc_now = 0.0;
            self.cje_eval = 0.0;
            self.cjc_eval = 0.0;
            self.q_je_eval = 0.0;
            self.q_jc_eval = 0.0;
            self.cjs_eval = 0.0;
            self.q_js_eval = 0.0;
            self.cjcx_eval = 0.0;
            self.q_jcx_eval = 0.0;
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
        // The substrate junction's Norton source: current out of the substrate
        // and into the internal collector.
        if let Some(sb) = self.substrate {
            b[sb] -= self.isub_eq;
        }
        if let Some(c) = self.collector {
            b[c] += self.isub_eq;
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

        // Collector-substrate junction, between the substrate and the *internal*
        // collector. `stamp!` drops the grounded rows, which is what makes a
        // grounded substrate a conductance to ground rather than a missing one.
        let gsub = self.gsub;
        stamp!(self.substrate, self.substrate, gsub);
        stamp!(self.substrate, c, -gsub);
        stamp!(c, self.substrate, -gsub);
        stamp!(c, c, gsub);

        // Series ohmic resistances: a conductance 1/R between each external
        // terminal and its internal node.  When R = 0 the internal node aliases
        // the external one and no resistor is stamped.  The `stamp!` macro skips
        // grounded (None) terminals, which correctly yields a conductance-to-
        // ground when an external terminal is ground.
        if self.rb > 0.0 {
            // `rb_eff`, not `rb`: `RBM`/`IRB` make it fall with base current.
            // Equal to `rb` for a card that gives neither.
            let g = 1.0 / self.rb_eff;
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

        // Transit-time diffusion charge: TF·(1+xf)·IF (B-E) and TR·IR (B-C).
        //
        // `cv` is the term the Jacobian stamp contributes at the linearisation
        // point, which the history current has to cancel. The B-E charge depends
        // on *two* voltages once `VTF`, `VAF` or `IKF` is finite, so `cv` sums
        // both contributions — otherwise the residual and the Jacobian would be
        // linearised about different points and the Newton step would be wrong
        // by the transcapacitance.
        if self.cbe_eff != 0.0 || self.cbe_x_eff != 0.0 {
            let cv = self.cbe_eff * self.vbe_eff + self.cbe_x_eff * self.vbc_eff;
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
        // The external share of CJC, between the base *pin* and the internal
        // collector. Spelled out rather than routed through `cap`, which is
        // hard-wired to the internal base.
        if self.cjcx_eval != 0.0 {
            let cv = self.cjcx_eval * self.vbcx_eff;
            let (i_hist, _) = self.q_jcx_hist.companion(disc, alpha, self.q_jcx_eval, cv);
            if let Some(bx) = self.base_ext {
                b[bx] += pol * i_hist;
            }
            if let Some(c) = self.collector {
                b[c] -= pol * i_hist;
            }
        }
        // Substrate depletion charge: CJS, between the substrate and the internal
        // collector. `cap` above is hard-wired to the base, so this one is spelled
        // out rather than routed through it.
        if self.cjs_eval != 0.0 {
            let cv = self.cjs_eval * self.vsub_eff;
            let (i_hist, _) = self.q_js_hist.companion(disc, alpha, self.q_js_eval, cv);
            if let Some(sb) = self.substrate {
                b[sb] += pol * i_hist;
            }
            if let Some(c) = self.collector {
                b[c] -= pol * i_hist;
            }
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
        // The transcapacitance: the B-E charge's current flows base to emitter,
        // so the *rows* are base and emitter, and it varies with `vbc`, so the
        // *columns* are base and collector. Asymmetric, unlike every other stamp
        // in this device.
        if self.cbe_x_eff != 0.0 {
            let c_x = scale * self.cbe_x_eff;
            stamp!(bk, bk, c_x);
            stamp!(bk, c, -c_x);
            stamp!(e, bk, -c_x);
            stamp!(e, c, c_x);
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
        // The external share of CJC: base pin to internal collector.
        if self.cjcx_eval != 0.0 {
            let g = scale * self.cjcx_eval;
            let bx = self.base_ext;
            stamp!(bx, bx, g);
            stamp!(bx, c, -g);
            stamp!(c, bx, -g);
            stamp!(c, c, g);
        }
        // Collector-substrate depletion cap: CJS
        if self.cjs_eval != 0.0 {
            let g_js = scale * self.cjs_eval;
            let sb = self.substrate;
            stamp!(sb, sb, g_js);
            stamp!(sb, c, -g_js);
            stamp!(c, sb, -g_js);
            stamp!(c, c, g_js);
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
        self.q_je_hist.advance(
            disc,
            q_depl(
                self.cje * self.cje_t_factor,
                vbe_eff,
                self.vje_t,
                self.mje,
                self.fc,
            ),
        );
        self.q_jc_hist.advance(
            disc,
            q_depl(
                self.cjc * self.cjc_t_factor,
                vbc_eff,
                self.vjc_t,
                self.mjc,
                self.fc,
            ),
        );
        let xcjc = self.xcjc.clamp(0.0, 1.0);
        let vbcx_eff = pol * (self.base_ext.map_or(0.0, |i| x[i]) - vc);
        self.q_jcx_hist.advance(
            disc,
            (1.0 - xcjc)
                * q_depl(
                    self.cjc * self.cjc_t_factor,
                    vbcx_eff,
                    self.vjc_t,
                    self.mjc,
                    self.fc,
                ),
        );
        let vsub_eff = pol * (self.substrate.map_or(0.0, |i| x[i]) - vc);
        self.q_js_hist
            .advance(disc, q_js_charge(self.cjs, vsub_eff, self.vjs, self.mjs));
    }

    /// Every reactance this device stamps in transient, for `.ac` and `.noise`.
    ///
    /// This used to be absent. The BJT builds its own companions in
    /// `load_jacobian_tran` and overrode neither this nor `reactive_branches`,
    /// whose default is an empty list — so `.ac` and `.noise` saw a transistor
    /// with no capacitance. Measured before the fix: a 1 kΩ resistor into the base
    /// of a device with `CJE = CJC = CJS = 100p` read `|V(b)| = 1.000000` at 1 kHz,
    /// 1 MHz, 10 MHz and 100 MHz, with the corner at 1.59 MHz. Flat, silently.
    ///
    /// All five reactances, and they must stay all five: the same list
    /// `load_jacobian_tran` stamps, read from the same cached values, so the two
    /// paths cannot disagree about what the device contains.
    fn small_signal_reactances(&self) -> Vec<crate::device::ReactiveBranchSpec> {
        use crate::device::{ReactiveBranchSpec, ReactiveKind};
        let cap = |pos, neg, value| ReactiveBranchSpec {
            kind: ReactiveKind::Capacitor,
            pos,
            neg,
            value,
            // `.ac`/`.noise` want the small-signal C itself, not a charge
            // branch's `∂q/∂v`, so zero is correct rather than conservative.
            dvalue_dstate: 0.0,
        };
        let mut v = Vec::new();
        // The transcapacitance is deliberately absent from this list and is
        // stamped by `load_reactive_jacobian` instead. A `ReactiveBranchSpec` is a
        // two-terminal branch and cannot express a charge whose rows and columns
        // differ.
        for (pos, neg, value) in [
            // Transit-time diffusion capacitance: TF on B-E, TR on B-C.
            (self.base, self.emitter, self.cbe_eff),
            (self.base, self.collector, self.cbc_eff),
            // Depletion capacitance: CJE, CJC, and CJS on collector-substrate.
            (self.base, self.emitter, self.cje_eval),
            (self.base, self.collector, self.cjc_eval),
            (self.base_ext, self.collector, self.cjcx_eval),
            (self.substrate, self.collector, self.cjs_eval),
        ] {
            if value != 0.0 {
                v.push(cap(pos, neg, value));
            }
        }
        v
    }

    /// The one reactance that is not a two-terminal branch: `dQBE/dVBC`.
    ///
    /// `.ac` and `.noise` form their susceptance block as `jw·C`, so what lands
    /// here is the frequency-domain twin of the `scale·cbe_x_eff` that
    /// `load_jacobian_tran` stamps — the same four cells, the same asymmetry.
    fn load_reactive_jacobian(&self, c_mat: &mut [crate::mna::SparseRow]) {
        if self.cbe_x_eff == 0.0 {
            return;
        }
        let c = self.cbe_x_eff;
        for (row, col, val) in [
            (self.base, self.base, c),
            (self.base, self.collector, -c),
            (self.emitter, self.base, -c),
            (self.emitter, self.collector, c),
        ] {
            if let (Some(r), Some(cc)) = (row, col) {
                c_mat[r][cc] += val;
            }
        }
    }

    fn noise_sources(&self, ctx: &SimContext, freq: f64) -> Vec<(NodeId, NodeId, f64)> {
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
            // Flicker rides with the base shot noise: SPICE drives 1/f from the
            // *base* current, not the collector's, and both sit across the same
            // terminal pair. Uncorrelated, so the densities add.
            let ib_mag = self.ib_eval.abs();
            let flicker = if self.kf > 0.0 && freq > 0.0 {
                self.kf * ib_mag.powf(self.af) / freq
            } else {
                0.0
            };
            sources.push((self.base, self.emitter, 2.0 * Q_E * ib_mag + flicker));
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
