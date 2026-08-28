use crate::device::{Device, Discretisation, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;
use crate::reactive::ChargeHistory;

/// Seed conductance before the first `eval`, so a device that is stamped once
/// prior to being evaluated cannot present a singular row. The *operating*
/// `gmin` comes from `SimContext::gmin`; this is only a non-zero starting point.
const GMIN_SEED: f64 = 1e-12;

/// SPICE Shockley diode model.
///
/// DC: Id = Is * (exp(Vd_j / (N·Vt)) − 1), with RS series resistance.
/// Transient: depletion capacitance Cj(V) and transit-time diffusion charge TT·Id.
///
/// All SPICE Level-1 diode parameters are accepted from `.model` cards.
/// Unrecognised parameters (BV, IBV, EG, XTI, …) emit a warning via the registry.
///
/// The instance parameter `AREA` scales the junction: IS and CJO with it, RS
/// against it. It is applied where the parameters are *used* rather than folded
/// into them, because a diode is built before its instance params are applied
/// (`ParamSet::apply` runs after `setup_model`), and the alias path can apply
/// them a second time.
pub struct ShockleyDiode {
    // Model parameters
    is: f64, // saturation current (A)
    n: f64,  // ideality factor
    rs: f64, // series resistance (Ω)
    /// `KF`/`AF` — flicker noise coefficient and exponent. `KF = 0` is off,
    /// which is SPICE's default and means the density is exactly zero rather
    /// than small.
    kf: f64,
    af: f64,
    /// `TNOM` — the temperature the card's parameters were extracted at.
    tnom: f64,
    /// `EG` — activation energy, eV. Silicon's 1.11 by default.
    eg: f64,
    /// `XTI` — the saturation current's temperature exponent.
    xti: f64,
    /// `IS(T)/IS`, from [`crate::temperature::diode_is_factor`].
    ///
    /// A factor rather than a scaled `IS`, so `setup_model` running more than
    /// once cannot apply it twice — and so `AREA`, which arrives afterwards, still
    /// multiplies cleanly on top.
    is_t_factor: f64,
    /// `BV` — reverse breakdown voltage as a positive magnitude. `None` when the
    /// card gives none, which is not the same as `0.0`: a card with `BV=0` breaks
    /// down immediately and a card without one never does.
    bv: Option<f64>,
    /// `IBV` — the reverse current *at* `-BV`. SPICE's default is 1 mA.
    ibv: f64,
    /// `BV` shifted so `I(-BV)` comes out as exactly `IBV` — see
    /// [`Self::adjusted_bv`]. `None` until the first `eval`, and while `bv` is.
    bv_adj: Option<f64>,
    /// The `N·vt` `bv_adj` was derived at, so a `.temp` sweep re-derives it and
    /// a normal solve derives it once. `AREA` is deliberately *not* in the key:
    /// see below.
    ///
    /// # AREA and the knee: a divergence from ngspice, on purpose
    ///
    /// `bv_adj` is derived from the **unit-area** `IS`, so the knee sits at the
    /// same voltage whatever `AREA` says, and the breakdown current scales with
    /// `AREA` like every other current in this model.
    ///
    /// Deriving it from `IS·AREA` instead — which is what ngspice does — makes the
    /// breakdown branch *exactly independent of AREA*: doubling `IS` doubles the
    /// prefactor and lifts `bv_adj` by `vte·ln 2`, and the two cancel to the last
    /// bit. Measured in ngspice-46: `area=1` and `area=2` return the same current
    /// at 4.8, 5.0, 5.1 and 5.3 V, ratio 1.0000, while the forward current
    /// doubles correctly.
    ///
    /// ngspice does not agree with *itself* there: two diodes in parallel give
    /// exactly twice the breakdown current that one with `area=2` gives. This
    /// tree already decided which of those wins — `area_scales_the_diode_exactly`
    /// asserts that "AREA=2 *is* two devices, so it has to agree with two devices
    /// rather than merely being twice something" — and a `area=10` Zener silently
    /// carrying a tenth of its knee current is the failure this codebase exists to
    /// refuse. So: AREA scales breakdown here, `area=N` equals N in parallel, and
    /// the divergence is against ngspice's `area=N` only.
    bv_key: f64,
    cjo: f64,   // zero-bias junction capacitance (F)
    vj: f64,    // junction built-in potential (V)
    mj: f64,    // grading coefficient
    fc: f64,    // forward-bias depletion-cap coefficient
    tt: f64,    // transit time (s)
    area: f64,  // instance AREA multiplier (1.0 = one device)
    vcrit: f64, // pnjlim critical voltage (derived from Is, Vt)

    // Terminal bindings
    anode: NodeId,
    cathode: NodeId,

    // NR state (updated by every eval)
    vd_prev: f64,     // junction Vd from last eval — pnjlim "old" reference
    id_junction: f64, // Id through the Shockley junction (before RS drop)
    gd_junction: f64, // dId/dVd at the junction
    gd_eff: f64,      // effective terminal conductance = gd_j / (1 + gd_j·RS)
    jeq_eff: f64,     // effective Norton source at terminal pair
    vd_j_eval: f64,   // junction voltage at last eval (used by commit_timestep)
    cj_total: f64,    // Cj_depl(Vd_j) + TT·gd_junction (0.0 when DC)
    q_at_vd_j: f64,   // Q_total(Vd_j) at current NR iterate (0.0 when DC)

    // Reactive history (updated by commit_timestep, read by load_*_tran)
    q_hist: ChargeHistory,
    /// The integrator's discretisation, captured during `eval` because
    /// `load_*_tran` receives only `alpha` — which is Backward Euler and
    /// nothing else. `None` outside the transient loop.
    disc: Option<Discretisation>,
}

impl ShockleyDiode {
    /// Construct with explicit model parameters. Call `setup_model` before first `eval`.
    pub fn new(is: f64, n: f64) -> Self {
        ShockleyDiode {
            is,
            n,
            rs: 0.0,
            kf: 0.0,
            af: 1.0,
            tnom: crate::temperature::TNOM_DEFAULT_K,
            eg: crate::temperature::EG_DEFAULT,
            xti: crate::temperature::XTI_DEFAULT,
            is_t_factor: 1.0,
            bv: None,
            ibv: 1e-3,
            bv_adj: None,
            bv_key: f64::NAN,
            cjo: 0.0,
            vj: 1.0,
            mj: 0.5,
            fc: 0.5,
            tt: 0.0,
            area: 1.0,
            vcrit: 0.0,
            anode: None,
            cathode: None,
            vd_prev: 0.0,
            id_junction: 0.0,
            gd_junction: GMIN_SEED,
            gd_eff: GMIN_SEED,
            jeq_eff: 0.0,
            vd_j_eval: 0.0,
            cj_total: 0.0,
            q_at_vd_j: 0.0,
            q_hist: ChargeHistory::default(),
            disc: None,
        }
    }

    /// Build from a model-card key=value list.
    /// Returns the device and any unrecognised param names (warned by the registry).
    pub fn from_params(params: &[(String, f64)]) -> (Self, Vec<String>) {
        let mut is = 1e-14;
        let mut n = 1.0;
        let mut rs = 0.0_f64;
        let mut cjo = 0.0_f64;
        let mut vj = 1.0_f64;
        let mut mj = 0.5_f64;
        let mut fc = 0.5_f64;
        let mut tt = 0.0_f64;
        let mut kf = 0.0_f64;
        let mut af = 1.0_f64;
        let mut tnom_c = crate::temperature::TNOM_DEFAULT_K - 273.15;
        let mut eg = crate::temperature::EG_DEFAULT;
        let mut xti = crate::temperature::XTI_DEFAULT;
        let mut bv: Option<f64> = None;
        let mut ibv = 1e-3_f64;
        let mut unknown = Vec::new();
        for (k, v) in params {
            match k.as_str() {
                "is" => is = *v,
                "n" => n = *v,
                "rs" => rs = *v,
                "cjo" | "cj0" => cjo = *v,
                "vj" => vj = *v,
                "m" | "mj" => mj = *v,
                "fc" => fc = *v,
                "tt" => tt = *v,
                // `BV` is a magnitude in every SPICE dialect, and a card writing
                // it negative means the same diode. `abs` rather than a refusal:
                // `BV=-5` is unambiguous, and erroring on it would reject cards
                // ngspice reads.
                // `TNOM` is degrees Celsius on the card, like `.temp`.
                "kf" => kf = *v,
                "af" => af = *v,
                "tnom" => tnom_c = *v,
                "eg" => eg = *v,
                "xti" => xti = *v,
                "bv" => bv = Some(v.abs()),
                "ibv" => ibv = v.abs(),
                // Accepted and NOT modelled: `crate::unmodelled` owns that list
                // and the diagnostic that reads it.
                k if crate::unmodelled::is_listed(crate::unmodelled::DIODE, k) => {}
                _ => unknown.push(k.clone()),
            }
        }
        let mut d = Self::new(is, n);
        d.rs = rs;
        d.cjo = cjo;
        d.vj = vj;
        d.mj = mj;
        d.fc = fc;
        d.tt = tt;
        d.kf = kf;
        d.af = af;
        d.tnom = tnom_c + 273.15;
        d.eg = eg;
        d.xti = xti;
        d.bv = bv;
        d.ibv = ibv;
        (d, unknown)
    }

    /// Saturation current of the whole instance: `IS(T)·AREA`.
    fn is_eff(&self) -> f64 {
        self.is * self.is_t_factor * self.area
    }

    /// Zero-bias junction capacitance of the whole instance: `CJO·AREA`.
    fn cjo_eff(&self) -> f64 {
        self.cjo * self.area
    }

    /// Series resistance of the whole instance: `RS/AREA` — N junctions in
    /// parallel each carry their own RS.
    fn rs_eff(&self) -> f64 {
        self.rs / self.area
    }

    /// SPICE pnjlim: logarithmically compress large voltage steps.
    /// # Why there is no mirrored version of this for reverse breakdown
    ///
    /// The breakdown exponential is as steep as the forward one, so mirroring
    /// this limiter about `-bv_adj` looks obviously right. It was written, and it
    /// turned out to be a convergence trap that produces a silent wrong answer.
    ///
    /// `vold` is state the outer Newton cannot see. The mirror compresses the walk
    /// into the knee logarithmically while the free node jumps to the supply in a
    /// single step, so the device's stamp keeps saying "barely conducting" and
    /// nothing pulls the node back. In that compressed region the terminal current
    /// is around 1e-11 A — under `abstol` — so the visible unknowns stop moving and
    /// Newton reports success. Measured at `.options vmax=1e6`: a 12 V / 1 kΩ
    /// Zener regulator read `out = 12 V` with the mirror and the correct 5.0501 V
    /// without it.
    ///
    /// This limiter is safe for the reason the mirror was not: while it is active
    /// the current changes by orders of magnitude per iteration, so a stalled walk
    /// cannot pass the convergence test. Reverse breakdown has a flat plateau
    /// under `abstol` instead.
    ///
    /// What bounds a step into breakdown now is the trust region
    /// (`vmax + reltol·|v|`), which covers both exponentials and keeps no state.
    /// `a_zener_regulator_regulates_with_a_loosened_trust_region` is the test that
    /// would have caught the mirror.
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

    /// `BV` shifted so the breakdown branch passes through `(-BV, -IBV)`.
    ///
    /// The card gives a knee voltage *and* a current at the knee, and both have
    /// to hold at once, so the exponential's offset is solved for:
    ///
    /// ```text
    /// IS·(exp((BV − bv_adj)/vte) − 1 + bv_adj/vte) = IBV
    /// ```
    ///
    /// Measured, not assumed: ngspice returns `I(-BV) = IBV` to every printed
    /// digit for any `N`, which is what says an adjustment happens at all. With
    /// `BV=5, IBV=1 mA, IS=10 fA` it gives 4.3449 V, and that predicts ngspice's
    /// current at 4.5 V (−4.02313e−12) exactly.
    ///
    /// # The clamp
    ///
    /// When `IBV` is below the current plain reverse saturation already gives at
    /// `-BV`, there is no offset to solve for — the knee sits under the card's own
    /// leakage floor — and `bv_adj` is `BV` unshifted. ngspice does the same and
    /// lands on `-IS` at `-BV`.
    ///
    /// # Why an iteration here is not the hazard `solve_series` was
    ///
    /// This depends only on the card and the temperature, and it iterates to a
    /// tolerance. It is not state carried between Newton iterations, so nothing
    /// about it is hidden from the outer solve. The map contracts hard: its
    /// derivative is `1/(IBV/IS + 1 − x/vte)`, about 1e-11 for a 1 mA knee on a
    /// 10 fA diode.
    fn adjusted_bv(&self, bv: f64, vte: f64) -> f64 {
        // The *unit-area* IS, deliberately — see the AREA note on `bv_adj`.
        let is = self.is;
        if is <= 0.0 || vte <= 0.0 || self.ibv < is * bv / vte {
            return bv;
        }
        let mut xbv = bv - vte * (1.0 + self.ibv / is).ln();
        for _ in 0..25 {
            let arg = self.ibv / is + 1.0 - xbv / vte;
            if arg <= 0.0 {
                return bv;
            }
            xbv = bv - vte * arg.ln();
            let at_xbv = is * (((bv - xbv) / vte).exp() - 1.0 + xbv / vte);
            if (at_xbv - self.ibv).abs() <= 1e-9 * self.ibv {
                break;
            }
        }
        xbv
    }

    /// `KF·|Id|^AF / f`, the SPICE flicker density. Zero when `KF` is unset.
    ///
    /// `freq <= 0` returns zero rather than an infinity: `.noise` never asks for
    /// DC, and the transient-noise path probes mid-band, but a caller that did
    /// would otherwise get a non-finite matrix instead of an error.
    fn flicker_density(&self, id_mag: f64, freq: f64) -> f64 {
        if self.kf <= 0.0 || freq <= 0.0 {
            return 0.0;
        }
        self.kf * id_mag.powf(self.af) / freq
    }

    /// Derive `bv_adj` when the inputs it depends on have moved. See `bv_key`.
    fn refresh_breakdown(&mut self, nvt: f64) {
        let Some(bv) = self.bv else { return };
        if self.bv_key != nvt {
            self.bv_adj = Some(self.adjusted_bv(bv, nvt));
            self.bv_key = nvt;
        }
    }

    /// Junction current and its slope at `vd_j` — the Shockley law plus `gmin`.
    ///
    /// One place, because `eval` and [`Self::solve_series`] must agree exactly:
    /// the series solve finds the root of a function `eval` then evaluates, and
    /// two spellings of the same law is two chances for them to differ.
    fn junction(&self, vd_j: f64, nvt: f64, gmin: f64) -> (f64, f64) {
        let is = self.is_eff();
        // Reverse breakdown, when the card asked for one: the Shockley
        // exponential mirrored about `-bv_adj`, so current runs away as the
        // junction is pushed past the knee. That runaway is the whole point of a
        // Zener and of an ESD clamp, and without it both block instead.
        if let Some(bv_adj) = self.bv_adj {
            if vd_j <= -bv_adj {
                let e = (-(vd_j + bv_adj) / nvt).exp();
                return (-is * e + gmin * vd_j, is * e / nvt + gmin);
            }
        }
        let exp_term = (vd_j / nvt).exp();
        (
            is * (exp_term - 1.0) + gmin * vd_j,
            is * exp_term / nvt + gmin,
        )
    }

    /// Solve `vd_j + RS·Id(vd_j) = vd_terminal` for the junction voltage.
    ///
    /// # Why this is a solve and not one step
    ///
    /// It used to be `vd_j = vd_terminal − Id·RS` with `Id` from the *previous*
    /// eval — one lagged fixed-point step per outer Newton iteration. `vd_j` is
    /// internal state the outer Newton cannot see, so its convergence test can
    /// be satisfied while the lag is still wide open. That is exactly what
    /// happens when a voltage source pins the diode's terminals: the visible
    /// unknowns stop moving on the first iteration and the lag never closes.
    /// Against ngspice, which gives the junction a real internal node, `RS=10`
    /// read 2.7% low at 0.7 V — and once `gmin` gave the reverse branch a
    /// conductance to converge onto, 100% wrong at 1.0 V.
    ///
    /// # Why scalar Newton is safe here
    ///
    /// `F(v) = v + RS·Id(v) − vd_terminal` has `F' = 1 + RS·gd ≥ 1` and is
    /// increasing, so the root is unique and Newton cannot divide by anything
    /// small. Overshoot into the exponential is the only hazard and `pnjlim` —
    /// the same limiter the outer loop uses — bounds every step. Typically two
    /// or three iterations from the warm start; `ITER_MAX` is a backstop, not a
    /// budget.
    ///
    /// Costs nothing when `RS = 0`: the caller skips it, which is every diode
    /// that does not set the parameter.
    ///
    /// `pnjlim` here is not the `nopnjlim` option's business. That option turns
    /// off junction limiting of the *outer* iterate, which perturbs where the
    /// solve goes; this is a step limiter inside a local root-find and the root
    /// it returns is the same with or without it.
    fn solve_series(&self, vd_terminal: f64, nvt: f64, gmin: f64, rs: f64, vt: f64) -> f64 {
        const ITER_MAX: usize = 100;
        // Warm start from the last junction voltage — across outer iterations
        // this is usually within millivolts, which is why the loop is cheap.
        let mut v = self.vd_prev;
        for _ in 0..ITER_MAX {
            let (id, gd) = self.junction(v, nvt, gmin);
            let f = v + rs * id - vd_terminal;
            let dv = -f / (1.0 + rs * gd);
            let v_new = self.pnjlim(v + dv, v, vt);
            let step = v_new - v;
            v = v_new;
            // Absolute *and* relative: `vd_j` is millivolts near the knee and
            // thousands of volts if a deck reverse-biases hard.
            if step.abs() <= 1e-14 * (1.0 + v.abs()) {
                break;
            }
        }
        v
    }

    /// SPICE depletion capacitance Cj(V).
    ///
    /// Below FC·VJ: Cj = CJO·(1 − V/VJ)^(−MJ).
    /// At and above FC·VJ: linear extrapolation to avoid singularity.
    fn cj_depl(&self, v: f64) -> f64 {
        let cjo = self.cjo_eff();
        if cjo == 0.0 {
            return 0.0;
        }
        let fc_vj = self.fc * self.vj;
        if v < fc_vj {
            cjo * (1.0 - v / self.vj).powf(-self.mj)
        } else {
            let k = (1.0 - self.fc).powf(1.0 + self.mj);
            cjo / k * (1.0 - self.fc * (1.0 + self.mj) + self.mj * v / self.vj)
        }
    }

    /// Charge integral Q(V) = ∫₀ᵛ Cj_depl dV (depletion charge only).
    fn q_depl(&self, v: f64) -> f64 {
        let cjo = self.cjo_eff();
        if cjo == 0.0 {
            return 0.0;
        }
        let fc_vj = self.fc * self.vj;
        if v < fc_vj {
            let x = 1.0 - v / self.vj;
            cjo * self.vj / (1.0 - self.mj) * (1.0 - x.powf(1.0 - self.mj))
        } else {
            // Charge at the FC·VJ boundary
            let x_fc = 1.0 - self.fc;
            let q_fc = cjo * self.vj / (1.0 - self.mj) * (1.0 - x_fc.powf(1.0 - self.mj));
            let k = x_fc.powf(1.0 + self.mj);
            let f2 = 1.0 - self.fc * (1.0 + self.mj);
            let dv = v - fc_vj;
            q_fc + cjo / k * (f2 * dv + self.mj / (2.0 * self.vj) * (v * v - fc_vj * fc_vj))
        }
    }

    /// Total charge Q_total(Vd_j) = Q_depl(Vd_j) + TT·Id(Vd_j).
    fn q_total(&self, vd_j: f64, id: f64) -> f64 {
        self.q_depl(vd_j) + self.tt * id
    }

    /// Stamp a conductance across (anode, cathode) into the Jacobian.
    fn stamp_g(&self, mat: &mut MnaMatrix, g: f64) {
        if let Some(ai) = self.anode {
            mat.a[ai][ai] += g;
            if let Some(ki) = self.cathode {
                mat.a[ai][ki] -= g;
            }
        }
        if let Some(ki) = self.cathode {
            mat.a[ki][ki] += g;
            if let Some(ai) = self.anode {
                mat.a[ki][ai] -= g;
            }
        }
    }
}

impl Device for ShockleyDiode {
    fn num_terminals(&self) -> usize {
        2
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        let vt = ctx.vt();
        // `IS` moves with temperature far more than `vt` does: at 125 C a silicon
        // junction leaks ~9e4 times its 27 C value, and `.temp` used to change
        // only the exponent. Idempotent, so a second `setup_model` is harmless.
        self.is_t_factor = crate::temperature::diode_is_factor(
            ctx.temperature,
            self.tnom,
            self.eg,
            self.xti,
            self.n,
        );
        // Deliberately the unit-area IS: `AREA` arrives after `setup_model`
        // (and can arrive twice, through the alias path). `vcrit` only decides
        // when pnjlim starts compressing steps, and AREA moves it by
        // vt·ln(AREA) — 18 mV at AREA=2, which changes no answer.
        self.vcrit = vt * (vt / (std::f64::consts::SQRT_2 * self.is)).ln();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], _ctx: &SimContext) {
        debug_assert_eq!(
            terminals.len(),
            2,
            "diode expects 2 terminals [anode, cathode]"
        );
        self.anode = terminals[0];
        self.cathode = terminals[1];
    }

    /// `AREA` — the only instance parameter this model honours. Everything else
    /// returns `false`, which is what makes `ParamSet::unconsumed` able to name
    /// it: a diode instance parameter used to reach the netlist and stop there,
    /// with nothing warning.
    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "area" if value > 0.0 => {
                self.area = value;
                true
            }
            _ => false,
        }
    }

    fn eval(&mut self, x: &[f64], flags: EvalFlags, ctx: &SimContext) {
        let v_a = self.anode.map_or(0.0, |i| x[i]);
        let v_k = self.cathode.map_or(0.0, |i| x[i]);
        let vd_terminal = v_a - v_k;

        let vt = ctx.vt();

        let nvt = self.n * vt;
        let rs = self.rs_eff();
        let gmin = ctx.gmin;
        // Before anything reads `junction`, including `solve_series`.
        self.refresh_breakdown(nvt);

        // Junction voltage. With RS this is a solve, not a step — see
        // `solve_series` for the lag it replaces. The series relation limits
        // `vd_j` on its own (it grows logarithmically with the terminal
        // voltage), so `pnjlim` is applied *inside* the solve rather than to its
        // answer: limiting the converged root would put the lag straight back.
        let vd_j = if rs > 0.0 {
            self.solve_series(vd_terminal, nvt, gmin, rs, vt)
        } else if ctx.jlim_enabled {
            self.pnjlim(vd_terminal, self.vd_prev, vt)
        } else {
            vd_terminal
        };
        self.vd_prev = vd_j;
        self.vd_j_eval = vd_j;
        // `gmin` is a real conductance across the junction — it carries current,
        // it is not only a floor under the Jacobian.
        //
        // It used to be added to `gd_junction` alone, and the Norton form below
        // is `jeq = Id − gd·Vd_j`, so at the operating point the terminal current
        // is `gd·Vd_j + jeq = Id` and the `gmin` term cancelled *exactly*. That
        // conditioned the matrix without contributing anything, which is a
        // legitimate technique and is not what SPICE means by `GMIN`: ngspice's
        // reverse-biased diode at −1 V carries `IS + gmin·1 V`, and a deck that
        // raises `.options gmin` sees the leakage follow it. Both are now true
        // here — and the value comes from the solve rather than from a `const`
        // in this file, so `.options gmin=` reaches the junction at all.
        let (id, gd) = self.junction(vd_j, nvt, gmin);
        self.id_junction = id;
        self.gd_junction = gd;

        // Norton equivalent at the terminal pair, accounting for RS.
        // Derivation: linearise Id(Vd_j) and Vd_j = Vd_term - Id·RS simultaneously.
        //   gd_eff = gd_j / (1 + gd_j·RS)
        //   jeq_eff = (Id - gd_j·Vd_j) / (1 + gd_j·RS)
        let denom = 1.0 + self.gd_junction * rs;
        self.gd_eff = self.gd_junction / denom;
        self.jeq_eff = (self.id_junction - self.gd_junction * vd_j) / denom;

        self.disc = ctx.discretisation;

        if flags.transient {
            self.cj_total = self.cj_depl(vd_j) + self.tt * self.gd_junction;
            self.q_at_vd_j = self.q_total(vd_j, self.id_junction);
        } else {
            self.cj_total = 0.0;
            self.q_at_vd_j = 0.0;
        }
    }

    fn load_residual(&self, b: &mut [f64]) {
        if let Some(a) = self.anode {
            b[a] -= self.jeq_eff;
        }
        if let Some(k) = self.cathode {
            b[k] += self.jeq_eff;
        }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        self.stamp_g(mat, self.gd_eff);
    }

    fn load_residual_tran(&self, b: &mut [f64], alpha: f64) {
        self.load_residual(b);
        if self.cj_total == 0.0 {
            return;
        }
        // Companion history current for the nonlinear junction cap, stamped as
        // a current source from cathode → anode (`b[anode] += i_hist`).
        //
        // The charge is `Q_total(Vd_j)`, so the Jacobian's contribution to the
        // residual is `scale·Cj·Vd_j` and the history has to cancel it —
        // `ChargeHistory` does both from one method interpretation. Under
        // Backward Euler this reduces to the old
        // `alpha·(Cj·Vd_j + Q_tprev − Q(Vd_j))` exactly.
        let (i_hist, _) = self.q_hist.companion(
            self.disc,
            alpha,
            self.q_at_vd_j,
            self.cj_total * self.vd_j_eval,
        );
        if let Some(a) = self.anode {
            b[a] += i_hist;
        }
        if let Some(k) = self.cathode {
            b[k] -= i_hist;
        }
    }

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, alpha: f64) {
        self.load_jacobian(mat);
        if self.cj_total == 0.0 {
            return;
        }
        self.stamp_g(mat, ChargeHistory::scale(self.disc, alpha) * self.cj_total);
    }

    /// Small-signal junction capacitance for `.ac`/`.noise`: the depletion +
    /// transit-time charge derivative `Cj_depl(Vd_j) + TT·gd`, between anode and
    /// cathode — identical to the value `load_jacobian_tran` stamps. Requires a
    /// preceding `eval(EvalFlags::tran())` to populate `cj_total`.
    fn small_signal_reactances(&self) -> Vec<crate::device::ReactiveBranchSpec> {
        use crate::device::{ReactiveBranchSpec, ReactiveKind};
        if self.cj_total == 0.0 {
            return Vec::new();
        }
        vec![ReactiveBranchSpec {
            kind: ReactiveKind::Capacitor,
            pos: self.anode,
            neg: self.cathode,
            value: self.cj_total,
            // ponytail: 0 leaves the Jacobian one term short for a *bias-
            // dependent* junction cap — `cj_total` carries both the depletion
            // C_j(V) and the transit-time TT·g_d, so its true `dC/dV` is
            // non-zero. That costs Newton iterations, not accuracy, and makes
            // an adjoint gradient through a diode's charge path wrong the way
            // `fc_pn_ps_cap`'s was (16 %). `jacobian_check_tran` finds it; do
            // the same thing here as `PnCapDrive` does when someone needs it.
            dvalue_dstate: 0.0,
        }]
    }

    fn commit_timestep(&mut self, x: &[f64]) {
        let va = self.anode.map_or(0.0, |i| x[i]);
        let vk = self.cathode.map_or(0.0, |i| x[i]);
        let vd_terminal = va - vk;
        // Use cached id_junction for RS correction; correct when called after eval.
        // On the first call (DC init, before any eval), id_junction=0 → vd_j ≈ vd_terminal.
        let vd_j = vd_terminal - self.id_junction * self.rs_eff();
        // Update pnjlim reference so the next timestep's first NR iter starts unlimted.
        self.vd_prev = vd_j;
        // Recomputed analytically from the converged solution rather than reused
        // from the last `eval`, which is one NR iterate behind it.
        self.q_hist
            .advance(self.disc, self.q_total(vd_j, self.id_junction));
    }

    /// Shot noise: i_n² = 2·q·|Id| (A²/Hz), between anode and cathode.
    ///
    /// Uses the bias-point Id cached by the most recent `eval()`.  Returns
    /// nothing when the junction is essentially off (|Id| < 1e-18 A) so the
    /// matrix doesn't pick up vanishing entries from leakage.
    fn noise_sources(&self, _ctx: &SimContext, freq: f64) -> Vec<(NodeId, NodeId, f64)> {
        let id_mag = self.id_junction.abs();
        if id_mag < 1e-18 {
            return Vec::new();
        }
        const Q: f64 = 1.602176634e-19;
        // Shot noise, plus flicker if the card asked for it. Both sit across the
        // junction, so they are one generator: they are uncorrelated, and adding
        // uncorrelated densities is what a single density means.
        let mut s_i = 2.0 * Q * id_mag;
        s_i += self.flicker_density(id_mag, freq);
        vec![(self.anode, self.cathode, s_i)]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::device::EvalFlags;

    fn ctx() -> SimContext {
        SimContext::default()
    }

    fn make_diode(is: f64, n: f64) -> ShockleyDiode {
        let mut d = ShockleyDiode::new(is, n);
        d.setup_model(&ctx());
        d.setup_instance(&[Some(0), None], &ctx());
        d
    }

    #[test]
    fn shockley_at_zero() {
        let mut d = make_diode(1e-14, 1.0);
        let x = [0.0f64];
        d.eval(&x, EvalFlags::dc(), &ctx());
        let vt = ctx().vt();
        let expected_gd = 1e-14 / vt + ctx().gmin;
        assert!(
            (d.gd_eff - expected_gd).abs() < 1e-18,
            "gd at Vd=0: {:.4e}",
            d.gd_eff
        );
        assert!(d.jeq_eff.abs() < 1e-30, "jeq at Vd=0: {:.4e}", d.jeq_eff);
    }

    #[test]
    fn shockley_forward_bias() {
        let is = 1e-14;
        let mut d = make_diode(is, 1.0);
        let vt = ctx().vt();
        let vd = 0.6;
        d.vd_prev = vd;
        let x = [vd];
        d.eval(&x, EvalFlags::dc(), &ctx());
        // `gmin` is now a real conductance across the junction, so it is in the
        // current as well as the slope. It used to be in the slope only, which
        // cancelled it out of `jeq` exactly — see `eval`.
        let gmin = ctx().gmin;
        let id_expected = is * ((vd / vt).exp() - 1.0) + gmin * vd;
        let gd_expected = is * (vd / vt).exp() / vt + gmin;
        let jeq_expected = id_expected - gd_expected * vd;
        assert!(
            (d.gd_eff - gd_expected).abs() / gd_expected < 1e-9,
            "gd mismatch"
        );
        assert!(
            (d.jeq_eff - jeq_expected).abs() / jeq_expected.abs() < 1e-9,
            "jeq mismatch"
        );
    }

    #[test]
    fn rs_reduces_effective_conductance() {
        let is = 1e-14;
        let mut d = ShockleyDiode::new(is, 1.0);
        d.rs = 1.0; // 1 Ω series resistance
        d.setup_model(&ctx());
        d.setup_instance(&[Some(0), None], &ctx());
        d.vd_prev = 0.6;
        let x = [0.6];
        d.eval(&x, EvalFlags::dc(), &ctx());
        // gd_eff < gd_junction
        assert!(
            d.gd_eff < d.gd_junction,
            "RS should reduce effective conductance"
        );
        let expected_gd_eff = d.gd_junction / (1.0 + d.gd_junction * 1.0);
        assert!(
            (d.gd_eff - expected_gd_eff).abs() / expected_gd_eff < 1e-9,
            "gd_eff formula"
        );
    }

    #[test]
    fn cjo_stamps_capacitive_conductance() {
        let mut d = ShockleyDiode::new(1e-14, 1.0);
        d.cjo = 4e-12; // 4 pF
        d.setup_model(&ctx());
        d.setup_instance(&[Some(0), Some(1)], &ctx());
        d.vd_prev = 0.0;
        let x = [0.0, 0.0];
        d.eval(&x, EvalFlags::tran(), &ctx());
        // At V=0: Cj = CJO (no voltage across the junction)
        let alpha = 1e8; // 1/h where h = 10 ns
        let mut mat = MnaMatrix::zeros(2);
        d.load_jacobian_tran(&mut mat, alpha);
        let g_cap_expected = alpha * 4e-12; // Cj(0) * alpha
                                            // Conductance stamp: a[0][0] includes gd_eff + g_cap; a[0][1] = -(gd_eff + g_cap)
        let g_total = d.gd_eff + g_cap_expected;
        assert!(
            (mat.a[0][0] - g_total).abs() / g_total < 1e-6,
            "Jacobian stamp"
        );
    }

    #[test]
    fn load_residual_stamps_correctly() {
        let is = 1e-14;
        let mut d = ShockleyDiode::new(is, 1.0);
        d.setup_model(&ctx());
        d.setup_instance(&[Some(0), Some(1)], &ctx());
        d.vd_prev = 0.6;
        let x = [0.6, 0.0];
        d.eval(&x, EvalFlags::dc(), &ctx());
        let mut b = [0.0f64; 2];
        d.load_residual(&mut b);
        assert!(
            (b[0] + d.jeq_eff).abs() < 1e-30,
            "b[anode] should be -jeq_eff"
        );
        assert!(
            (b[1] - d.jeq_eff).abs() < 1e-30,
            "b[cathode] should be +jeq_eff"
        );
    }

    #[test]
    fn pnjlim_compresses_large_steps() {
        let d = make_diode(1e-14, 1.0);
        let vt = ctx().vt();
        let limited = d.pnjlim(10.0, 0.0, vt);
        assert!(
            limited < 1.0,
            "pnjlim should limit large step: got {limited}"
        );
        assert!(limited > 0.0, "pnjlim result should be positive");
    }

    #[test]
    fn q_depl_at_zero_is_zero() {
        let mut d = ShockleyDiode::new(1e-14, 1.0);
        d.cjo = 4e-12;
        assert!(d.q_depl(0.0).abs() < 1e-30, "Q(0) should be 0");
    }

    #[test]
    fn q_depl_increases_with_forward_bias() {
        let mut d = ShockleyDiode::new(1e-14, 1.0);
        d.cjo = 4e-12;
        let q0 = d.q_depl(0.0);
        let q1 = d.q_depl(0.3);
        assert!(q1 > q0, "Q should increase with forward bias");
    }

    /// `small_signal_reactances` (consumed by .ac/.noise) reports the depletion
    /// capacitance at the operating point, between anode and cathode, matching
    /// the closed-form Cj = CJO·(1 − Vd/VJ)^(−M). Reverse bias Vd_j = −1 V.
    #[test]
    fn small_signal_cap_matches_depletion_formula_at_reverse_bias() {
        use crate::device::ReactiveKind;
        let (mut d, _) = ShockleyDiode::from_params(&[
            ("is".into(), 1e-14),
            ("cjo".into(), 2e-12),
            ("vj".into(), 0.8),
            ("m".into(), 0.5),
        ]);
        d.setup_model(&ctx());
        d.setup_instance(&[Some(0), Some(1)], &ctx()); // anode=node0, cathode=node1
        d.vd_prev = -1.0; // seed pnjlim so the reverse step is not limited
                          // anode=0 V, cathode=+1 V → Vd_terminal = −1 V (reverse).
        let x = [0.0, 1.0];
        d.eval(&x, EvalFlags::tran(), &ctx());

        let r = d.small_signal_reactances();
        assert_eq!(r.len(), 1, "one junction cap expected");
        assert_eq!(r[0].kind, ReactiveKind::Capacitor);
        assert_eq!(r[0].pos, Some(0));
        assert_eq!(r[0].neg, Some(1));
        // RS=0 and TT=0 here, so Cj = cj_depl(−1) exactly.
        let expected = 2e-12_f64 * (1.0 - (-1.0) / 0.8_f64).powf(-0.5); // 2pF/1.5 = 1.333 pF
        assert!(
            (r[0].value - expected).abs() / expected < 1e-9,
            "Cj={:.6e} expected={:.6e}",
            r[0].value,
            expected
        );

        // Under DC flags the cache is cleared → no small-signal cap reported
        // (so a pure .op stamps no spurious reactance).
        d.eval(&x, EvalFlags::dc(), &ctx());
        assert!(d.small_signal_reactances().is_empty());
    }
}
