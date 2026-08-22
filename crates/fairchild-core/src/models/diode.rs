use crate::device::{Device, Discretisation, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;
use crate::reactive::ChargeHistory;

/// Small conductance added to every p-n junction for numerical conditioning.
const GMIN: f64 = 1e-12;

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
    is: f64,    // saturation current (A)
    n: f64,     // ideality factor
    rs: f64,    // series resistance (Ω)
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
            gd_junction: GMIN,
            gd_eff: GMIN,
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
        (d, unknown)
    }

    /// Saturation current of the whole instance: `IS·AREA`.
    fn is_eff(&self) -> f64 {
        self.is * self.area
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

        // Junction voltage: iterate RS drop using Id from previous NR step.
        let vd_j_raw = vd_terminal - self.id_junction * self.rs_eff();
        let vd_j = if ctx.jlim_enabled {
            self.pnjlim(vd_j_raw, self.vd_prev, vt)
        } else {
            vd_j_raw
        };
        self.vd_prev = vd_j;
        self.vd_j_eval = vd_j;

        let nvt = self.n * vt;
        let exp_term = (vd_j / nvt).exp();
        let is = self.is_eff();
        self.id_junction = is * (exp_term - 1.0);
        self.gd_junction = is * exp_term / nvt + GMIN;

        // Norton equivalent at the terminal pair, accounting for RS.
        // Derivation: linearise Id(Vd_j) and Vd_j = Vd_term - Id·RS simultaneously.
        //   gd_eff = gd_j / (1 + gd_j·RS)
        //   jeq_eff = (Id - gd_j·Vd_j) / (1 + gd_j·RS)
        let denom = 1.0 + self.gd_junction * self.rs_eff();
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
    fn noise_sources(&self, _ctx: &SimContext, _freq: f64) -> Vec<(NodeId, NodeId, f64)> {
        let id_mag = self.id_junction.abs();
        if id_mag < 1e-18 {
            return Vec::new();
        }
        const Q: f64 = 1.602176634e-19;
        let s_i = 2.0 * Q * id_mag;
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
        let expected_gd = 1e-14 / vt + GMIN;
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
        let id_expected = is * ((vd / vt).exp() - 1.0);
        let gd_expected = is * (vd / vt).exp() / vt + GMIN;
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
