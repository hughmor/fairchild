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
//! Operating point: the same two relations, with the delayed terms equal to the
//! present ones. Adding them gives `i1 + i2 = 0` and subtracting them gives
//! `v1 = v2`, so at DC the line is an ideal through-connection and that is what
//! the branch rows stamp. ngspice agrees: a 1 V source through 50 Ω into a
//! 1 kΩ load across the line draws `1/(50+1000)`, not `1/(50+50)`.
//!
//! The `E = 0` seed the operating point used to stamp made each port a `Z0`
//! resistor and left the far end dead, so any deck biasing a device *through* a
//! line started from the wrong state, and `.tf` read the input resistance of a
//! terminator rather than of a line.

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

    /// The two branch rows at DC, where the delayed terms equal the present
    /// ones: `v1 − v2 = 0` in the first, `i1 + i2 = 0` in the second.
    ///
    /// Both rows are needed. Stamping `v1 = v2` twice leaves the current
    /// undetermined, and `gmin` would make that non-singular rather than an
    /// error.
    fn stamp_dc_rows(&self, mat: &mut MnaMatrix) {
        if let Some(j) = self.br1 {
            for (node, sign) in [
                (self.a_pos, 1.0),
                (self.a_neg, -1.0),
                (self.b_pos, -1.0),
                (self.b_neg, 1.0),
            ] {
                if let Some(i) = node {
                    mat.a[j][i] += sign;
                }
            }
        }
        if let (Some(j), Some(k)) = (self.br2, self.br1) {
            mat.a[j][k] += 1.0;
            mat.a[j][j] += 1.0;
        }
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
            // The DC through-connection is homogeneous: no source, only the
            // `v1 = v2` / `i1 + i2 = 0` rows in `load_jacobian`.
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
        // KCL is the same in both regimes: `i1` leaves port A into the device,
        // `i2` leaves port B.
        stamp_port_kcl(mat, self.a_pos, self.a_neg, self.br1);
        stamp_port_kcl(mat, self.b_pos, self.b_neg, self.br2);
        if self.delay.is_active() {
            // Each port is a voltage source `E` behind series Z0, with the
            // branch current as the unknown (source-with-resistance).
            stamp_branch_row(mat, self.a_pos, self.a_neg, self.br1, self.z0);
            stamp_branch_row(mat, self.b_pos, self.b_neg, self.br2, self.z0);
        } else {
            // DC limit of the same relations: `v1 = v2` and `i1 + i2 = 0`.
            self.stamp_dc_rows(mat);
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }

    /// The delayed coupling `exp(−jωTD)`, which is the frequency-domain twin of
    /// the history source `load_residual` writes.
    ///
    /// `load_jacobian` has already stamped the left-hand sides
    /// `V1 − Z0·I1` and `V2 − Z0·I2`, so what is missing is the right-hand
    /// sides as couplings rather than as a known:
    ///
    /// ```text
    ///   V1 − Z0·I1 − q·(V2 + Z0·I2) = 0
    ///   V2 − Z0·I2 − q·(V1 + Z0·I1) = 0,     q = exp(−jωTD)
    /// ```
    ///
    /// At `ω = 0` this is `q = 1`, and adding and subtracting the two rows
    /// returns `i1 + i2 = 0` and `v1 = v2` — the same through-connection the DC
    /// stamp puts in directly.
    fn ac_stamps(&self, omega: f64) -> Vec<crate::device::AcStamp> {
        use crate::device::AcStamp;
        let (qr, qi) = ((omega * self.td).cos(), -(omega * self.td).sin());
        let mut out = Vec::with_capacity(6);
        // Row `br` couples to the *far* port's voltage and current.
        let mut row = |br: NodeId, pos: NodeId, neg: NodeId, far_br: NodeId| {
            let Some(r) = br else { return };
            if let Some(c) = pos {
                out.push(AcStamp {
                    row: r,
                    col: c,
                    re: -qr,
                    im: -qi,
                });
            }
            if let Some(c) = neg {
                out.push(AcStamp {
                    row: r,
                    col: c,
                    re: qr,
                    im: qi,
                });
            }
            if let Some(c) = far_br {
                out.push(AcStamp {
                    row: r,
                    col: c,
                    re: -qr * self.z0,
                    im: -qi * self.z0,
                });
            }
        };
        row(self.br1, self.b_pos, self.b_neg, self.br2);
        row(self.br2, self.a_pos, self.a_neg, self.br1);
        out
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

/// Couple a port current into the two node KCL rows: the branch current leaves
/// `pos` and enters `neg`.  Grounded (`None`) terminals are skipped, which
/// correctly drops the corresponding entries.
fn stamp_port_kcl(mat: &mut MnaMatrix, pos: NodeId, neg: NodeId, br: NodeId) {
    let Some(j) = br else { return };
    if let Some(p) = pos {
        mat.a[p][j] += 1.0;
    }
    if let Some(n) = neg {
        mat.a[n][j] -= 1.0;
    }
}

/// Stamp one Branin branch row: `V(pos) − V(neg) − Z0·I = E` (the RHS `E` is
/// added in `load_residual`).
fn stamp_branch_row(mat: &mut MnaMatrix, pos: NodeId, neg: NodeId, br: NodeId, z0: f64) {
    let Some(j) = br else { return };
    if let Some(p) = pos {
        mat.a[j][p] += 1.0;
    }
    if let Some(n) = neg {
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
    fn dc_stamps_a_through_connection_not_a_terminator() {
        // At DC the delayed terms equal the present ones, so the two relations
        // collapse to `v1 = v2` and `i1 + i2 = 0`. Z0 must not appear.
        let mut t = NativeTLine::new(50.0, 1e-9);
        t.setup_instance(&[Some(0), None, Some(1), None], &SimContext::default());
        t.bind_extra_nodes(2); // br1=2, br2=3
        t.eval(&[0.0; 4], EvalFlags::dc(), &SimContext::default());
        let mut mat = MnaMatrix::zeros(4);
        t.load_jacobian(&mut mat);
        // Row br1: V(A+) − V(B+) = 0.
        assert_eq!(mat.a[2][0], 1.0, "branch-1 row: +V(A+)");
        assert_eq!(mat.a[2][1], -1.0, "branch-1 row: −V(B+)");
        // Row br2: i1 + i2 = 0.
        assert_eq!(mat.a[3][2], 1.0, "branch-2 row: +i1");
        assert_eq!(mat.a[3][3], 1.0, "branch-2 row: +i2");
        // KCL unchanged.
        assert_eq!(mat.a[0][2], 1.0, "A+ KCL coupling to i1");
        assert_eq!(mat.a[1][3], 1.0, "B+ KCL coupling to i2");
        // No Z0 anywhere: a terminator would put −50 on a branch diagonal.
        for row in 0..4 {
            for col in 0..4 {
                assert!(
                    mat.a[row][col].abs() <= 1.0,
                    "Z0 leaked into the DC stamp at [{row}][{col}]: {}",
                    mat.a[row][col]
                );
            }
        }
        let mut b = vec![0.0; 4];
        t.load_residual(&mut b);
        assert!(b.iter().all(|&v| v == 0.0), "the DC row is homogeneous");
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
