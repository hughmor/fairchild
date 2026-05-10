use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

/// Small conductance added to every p-n junction for numerical conditioning.
/// ngspice default is 1e-12 S.
const GMIN: f64 = 1e-12;

/// Shockley ideal-diode model: Id = Is * (exp(Vd / (N·Vt)) − 1).
///
/// Implements `Device` with pnjlim voltage limiting for Newton-Raphson convergence.
pub struct ShockleyDiode {
    // --- model params (set at construction, finalised in setup_model) ---
    is: f64,    // saturation current (A)
    n: f64,     // ideality factor
    vcrit: f64, // pnjlim critical voltage (computed from Is and Vt in setup_model)

    // --- instance bindings (set in setup_instance) ---
    anode: NodeId,
    cathode: NodeId,

    // --- eval state ---
    vd_prev: f64, // operating-point Vd from the last eval (pnjlim "old" value)
    gd: f64,      // Norton conductance cached by eval
    jeq: f64,     // Norton current source: Id(vd_lim) - gd * vd_lim
}

impl ShockleyDiode {
    /// Construct with explicit Is and N. Call `setup_model` before first `eval`.
    pub fn new(is: f64, n: f64) -> Self {
        ShockleyDiode {
            is,
            n,
            vcrit: 0.0,
            anode: None,
            cathode: None,
            vd_prev: 0.0,
            gd: GMIN,
            jeq: 0.0,
        }
    }

    /// Build from a list of model-card key=value pairs.
    /// Unrecognised keys are silently ignored; missing keys keep defaults.
    pub fn from_params(params: &[(String, f64)]) -> Self {
        let mut is = 1e-14;
        let mut n = 1.0;
        for (k, v) in params {
            match k.as_str() {
                "is" => is = *v,
                "n"  => n  = *v,
                _ => {}
            }
        }
        Self::new(is, n)
    }

    /// SPICE pnjlim: logarithmically compress the voltage step when Vd > vcrit.
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

impl Device for ShockleyDiode {
    fn num_terminals(&self) -> usize { 2 }

    fn setup_model(&mut self, ctx: &SimContext) {
        let vt = ctx.vt();
        // vcrit = Vt * ln(Vt / (√2 · Is))
        self.vcrit = vt * (vt / (std::f64::consts::SQRT_2 * self.is)).ln();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], _ctx: &SimContext) {
        debug_assert_eq!(terminals.len(), 2, "diode expects 2 terminals [anode, cathode]");
        self.anode   = terminals[0];
        self.cathode = terminals[1];
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, ctx: &SimContext) {
        let v_a = self.anode.map_or(0.0, |i| x[i]);
        let v_k = self.cathode.map_or(0.0, |i| x[i]);
        let vd_circuit = v_a - v_k;

        let vt = ctx.vt();
        let vd = self.pnjlim(vd_circuit, self.vd_prev, vt);
        self.vd_prev = vd;

        let nvt = self.n * vt;
        let exp_term = (vd / nvt).exp();
        let id = self.is * (exp_term - 1.0);
        self.gd  = self.is * exp_term / nvt + GMIN;
        self.jeq = id - self.gd * vd;
    }

    fn load_residual(&self, b: &mut [f64]) {
        // Norton current: Jeq flows from anode to cathode through companion source.
        if let Some(a) = self.anode   { b[a] -= self.jeq; }
        if let Some(k) = self.cathode { b[k] += self.jeq; }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let a = self.anode;
        let k = self.cathode;
        let g = self.gd;
        // Conductance gd between anode and cathode (same stamp as a resistor).
        if let Some(ai) = a {
            mat.a[ai][ai] += g;
            if let Some(ki) = k { mat.a[ai][ki] -= g; }
        }
        if let Some(ki) = k {
            mat.a[ki][ki] += g;
            if let Some(ai) = a { mat.a[ki][ai] -= g; }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::device::EvalFlags;

    fn ctx() -> SimContext { SimContext::default() }

    #[test]
    fn shockley_at_zero() {
        // At Vd=0, Id=0 and gd=Is/Vt + GMIN.
        let mut d = ShockleyDiode::new(1e-14, 1.0);
        d.setup_model(&ctx());
        d.setup_instance(&[Some(0), None], &ctx());
        let x = [0.0f64];
        d.eval(&x, EvalFlags::dc(), &ctx());
        let vt = ctx().vt();
        let expected_gd = 1e-14 / vt + GMIN;
        assert!((d.gd - expected_gd).abs() < 1e-18, "gd at Vd=0: got {:.4e}", d.gd);
        // Jeq = Id - gd*Vd = 0 - 0 = 0
        assert!(d.jeq.abs() < 1e-30, "jeq at Vd=0: got {:.4e}", d.jeq);
    }

    #[test]
    fn shockley_forward_bias() {
        // At Vd=0.6V: Id ≈ Is·exp(Vd/Vt)
        let is = 1e-14;
        let mut d = ShockleyDiode::new(is, 1.0);
        d.setup_model(&ctx());
        d.setup_instance(&[Some(0), None], &ctx());
        let vt = ctx().vt();
        let vd = 0.6;
        // Drive vd_prev to 0.6 so pnjlim doesn't limit (previous step already there).
        d.vd_prev = vd;
        let x = [vd];
        d.eval(&x, EvalFlags::dc(), &ctx());
        let id_expected = is * ((vd / vt).exp() - 1.0);
        let gd_expected = is * (vd / vt).exp() / vt + GMIN;
        let jeq_expected = id_expected - gd_expected * vd;
        assert!((d.gd - gd_expected).abs() / gd_expected < 1e-9, "gd mismatch");
        assert!((d.jeq - jeq_expected).abs() / jeq_expected.abs() < 1e-9, "jeq mismatch");
    }

    #[test]
    fn load_residual_stamps_into_b() {
        let is = 1e-14;
        let mut d = ShockleyDiode::new(is, 1.0);
        d.setup_model(&ctx());
        // anode=node 0, cathode=node 1
        d.setup_instance(&[Some(0), Some(1)], &ctx());
        d.vd_prev = 0.6;
        let x = [0.6, 0.0];
        d.eval(&x, EvalFlags::dc(), &ctx());
        let mut b = [0.0f64; 2];
        d.load_residual(&mut b);
        // b[anode] -= jeq, b[cathode] += jeq
        assert!((b[0] + d.jeq).abs() < 1e-30, "b[anode] should be -jeq");
        assert!((b[1] - d.jeq).abs() < 1e-30, "b[cathode] should be +jeq");
    }

    #[test]
    fn pnjlim_compresses_large_steps() {
        let mut d = ShockleyDiode::new(1e-14, 1.0);
        d.setup_model(&ctx());
        let vt = ctx().vt();
        // From 0 to 10V: should be compressed to around vcrit + small.
        let limited = d.pnjlim(10.0, 0.0, vt);
        assert!(limited < 1.0, "pnjlim should limit large step: got {limited}");
        assert!(limited > 0.0, "pnjlim result should be positive");
    }
}
