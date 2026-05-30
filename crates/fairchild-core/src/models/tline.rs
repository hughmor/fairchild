//! Lossless transmission line (`T` element), Branin's method of characteristics.
//!
//! A lossless line of characteristic impedance `Z0` and one-way delay `TD` is
//! exactly modelled by two coupled port relations (with port voltages `v1, v2`
//! and the currents `i1, i2` flowing into each port from the external circuit):
//!
//! ```text
//!   v1(t) − Z0·i1(t) = v2(t−TD) + Z0·i2(t−TD)        (port A)
//!   v2(t) − Z0·i2(t) = v1(t−TD) + Z0·i1(t−TD)        (port B)
//! ```
//!
//! Each line is therefore an independent voltage source `E` (the far port's
//! delayed travelling wave) behind a series resistance `Z0`, with the port
//! current as an explicit MNA branch unknown:
//!
//! ```text
//!   V(A+) − V(A−) − Z0·i1 = E1,   E1 = v2(t−TD) + Z0·i2(t−TD)
//!   V(B+) − V(B−) − Z0·i2 = E2,   E2 = v1(t−TD) + Z0·i1(t−TD)
//! ```
//!
//! Making `i1, i2` branch unknowns means the history snapshot `[v1, v2, i1, i2]`
//! is read directly from the solution vector — no back-substitution — and the
//! steady state self-consistently collapses to the DC limit of a lossless line
//! (`v1 = v2`, `i1 = −i2`: an ideal through-connection), because there the
//! delayed terms equal the present ones.
//!
//! The delay is intrinsic to the device and always modelled in transient runs
//! (unlike the photonic waveguide, whose group delay is opt-in). The generic
//! history/interpolation is provided by [`crate::delay::DelayLine`].
//!
//! Operating point: with no history the line presents `Z0` at each port (the
//! `E = 0` seed). For a transient started from rest this is exact; a circuit
//! that relies on a DC bias propagating *through* the line at `.op` will settle
//! to the correct through-connection within a few `TD` of transient but is not
//! captured by a standalone `.op`.

use crate::delay::DelayLine;
use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

pub struct NativeTLine {
    z0: f64,
    td: f64,
    // External terminals: A+, A−, B+, B−.
    a_pos: NodeId,
    a_neg: NodeId,
    b_pos: NodeId,
    b_neg: NodeId,
    // Branch-current unknowns: i1 (port A), i2 (port B).
    br1: Option<usize>,
    br2: Option<usize>,
    delay: DelayLine,
    // History-driven source values for the current step (E1, E2).
    e1: f64,
    e2: f64,
}

impl NativeTLine {
    pub fn new(z0: f64, td: f64) -> Self {
        NativeTLine {
            z0,
            td,
            a_pos: None,
            a_neg: None,
            b_pos: None,
            b_neg: None,
            br1: None,
            br2: None,
            delay: DelayLine::new(),
            e1: 0.0,
            e2: 0.0,
        }
    }

    #[inline]
    fn port_v(&self, x: &[f64], pos: NodeId, neg: NodeId) -> f64 {
        pos.map_or(0.0, |i| x[i]) - neg.map_or(0.0, |i| x[i])
    }
}

impl Device for NativeTLine {
    fn num_terminals(&self) -> usize {
        4
    }

    fn setup_model(&mut self, _ctx: &SimContext) {}

    fn setup_instance(&mut self, terminals: &[NodeId], _ctx: &SimContext) {
        debug_assert!(terminals.len() >= 4, "T-line expects [A+, A-, B+, B-]");
        self.a_pos = terminals[0];
        self.a_neg = terminals[1];
        self.b_pos = terminals[2];
        self.b_neg = terminals[3];
    }

    /// Two branch-current rows (i1, i2).
    fn num_extra_nodes(&self) -> usize {
        2
    }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        self.br1 = Some(first_idx);
        self.br2 = Some(first_idx + 1);
    }

    fn eval(&mut self, _x: &[f64], flags: EvalFlags, ctx: &SimContext) {
        let active = flags.transient && self.td > 0.0;
        self.delay.set_state(active, ctx.time_s);
        if active {
            // Delayed snapshot [v1, v2, i1, i2] at t − TD.
            let d = self.delay.sample(self.td, 4);
            let (v1d, v2d, i1d, i2d) = (d[0], d[1], d[2], d[3]);
            // E1 is the wave arriving at port A from port B one delay ago; E2
            // the reverse.
            self.e1 = v2d + self.z0 * i2d;
            self.e2 = v1d + self.z0 * i1d;
        } else {
            // No history → each port presents Z0 (the operating-point seed).
            self.e1 = 0.0;
            self.e2 = 0.0;
        }
    }

    fn load_residual(&self, b: &mut [f64]) {
        // Branch-equation RHS: V(pos) − V(neg) − Z0·i = E.
        if let Some(j) = self.br1 {
            b[j] += self.e1;
        }
        if let Some(j) = self.br2 {
            b[j] += self.e2;
        }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        // Stamp each port as a voltage source `E` behind series Z0, with the
        // branch current as the unknown (standard MNA source-with-resistance).
        stamp_branch(mat, self.a_pos, self.a_neg, self.br1, self.z0);
        stamp_branch(mat, self.b_pos, self.b_neg, self.br2, self.z0);
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }

    fn commit_timestep(&mut self, x: &[f64]) {
        if !self.delay.is_active() {
            return;
        }
        let v1 = self.port_v(x, self.a_pos, self.a_neg);
        let v2 = self.port_v(x, self.b_pos, self.b_neg);
        let i1 = self.br1.map_or(0.0, |j| x[j]);
        let i2 = self.br2.map_or(0.0, |j| x[j]);
        self.delay.record(vec![v1, v2, i1, i2], self.td);
    }
}

/// Stamp one Branin branch: `V(pos) − V(neg) − Z0·I = E` with branch current
/// `I` (the RHS `E` is added in `load_residual`).  Grounded (`None`) terminals
/// are skipped, which correctly drops the corresponding entries.
fn stamp_branch(mat: &mut MnaMatrix, pos: NodeId, neg: NodeId, br: NodeId, z0: f64) {
    let Some(j) = br else { return };
    // KCL: branch current leaves `pos`, enters `neg`.
    if let Some(p) = pos {
        mat.a[p][j] += 1.0;
        mat.a[j][p] += 1.0;
    }
    if let Some(n) = neg {
        mat.a[n][j] -= 1.0;
        mat.a[j][n] -= 1.0;
    }
    // Series characteristic impedance on the branch diagonal.
    mat.a[j][j] -= z0;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ctx_tran(t: f64) -> SimContext {
        SimContext {
            time_s: t,
            ..Default::default()
        }
    }

    #[test]
    fn dc_seed_presents_z0_at_each_port() {
        // With no history the operating point stamps a Z0 source-behind-resistor
        // at each port: branch row diagonal = −Z0, KCL coupling ±1.
        let mut t = NativeTLine::new(50.0, 1e-9);
        t.setup_instance(&[Some(0), None, Some(1), None], &SimContext::default());
        t.bind_extra_nodes(2); // br1=2, br2=3
        t.eval(&[0.0; 4], EvalFlags::dc(), &SimContext::default());
        let mut mat = MnaMatrix::zeros(4);
        t.load_jacobian(&mut mat);
        assert_eq!(mat.a[2][2], -50.0, "branch-1 diagonal = −Z0");
        assert_eq!(mat.a[3][3], -50.0, "branch-2 diagonal = −Z0");
        assert_eq!(mat.a[0][2], 1.0, "A+ KCL coupling to i1");
        assert_eq!(mat.a[2][0], 1.0, "branch-1 row references V(A+)");
        let mut b = vec![0.0; 4];
        t.load_residual(&mut b);
        assert!(b.iter().all(|&v| v == 0.0), "no source at the DC seed");
    }

    #[test]
    fn delayed_wave_appears_after_td() {
        // Record a step at port A (v1: 0→1) and matched current, then check the
        // history reconstructs E2 (the wave heading to port B) delayed by TD.
        let mut t = NativeTLine::new(50.0, 2.0);
        t.setup_instance(&[Some(0), None, Some(1), None], &SimContext::default());
        t.bind_extra_nodes(2);
        // t=0: v1=0. t=1: v1=1, i1=0.02 (1V into 50Ω matched). Record both.
        t.eval(&[0.0; 4], EvalFlags::tran(), &ctx_tran(0.0));
        t.delay.record(vec![0.0, 0.0, 0.0, 0.0], t.td);
        t.delay.set_state(true, 1.0);
        t.delay.record(vec![1.0, 0.0, 0.02, 0.0], t.td);
        // At t=3 (= 1 + TD), E2 should reflect the port-A wave launched at t=1:
        // E2 = v1(t−TD) + Z0·i1(t−TD) = 1 + 50·0.02 = 2.0.
        t.eval(&[0.0; 4], EvalFlags::tran(), &ctx_tran(3.0));
        assert!(
            (t.e2 - 2.0).abs() < 1e-9,
            "E2 should be the delayed port-A wave (2.0), got {}",
            t.e2
        );
        // Before the wave arrives (t=2, query at t−TD=0): E2 from the zero step.
        t.eval(&[0.0; 4], EvalFlags::tran(), &ctx_tran(2.0));
        assert!(t.e2.abs() < 1e-9, "no wave yet at t<1+TD, got {}", t.e2);
    }
}
