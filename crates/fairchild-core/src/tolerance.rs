//! Per-row convergence tolerances — one interpretation, shared by every solver.
//!
//! # Why this module exists
//!
//! The Newton step test is `|Δx_i| < abs + rel·|x_i|`, and it used to be
//! written out at six sites with `abs = vntol` and `rel = reltol` for **every**
//! row alike. But the unknown vector is not all voltages:
//!
//! | row class | quantity | what `vntol = 1e-6` means there |
//! |---|---|---|
//! | node | volts | 1 µV — correct, this is what it is for |
//! | voltage-source branch | amps | 1 µA — 10⁶× looser than SPICE's `abstol` |
//! | λ wire | **metres** (~1.55e-6) | 1 µm of wavelength; the C-band is 0.035 µm |
//!
//! The λ row is the one where the mismatch is fatal. Newton could stop with λ
//! ~10 pm out, which is a real detuning for a 40 GHz passband and the same
//! order as a PN phase shifter's 13 pm/V tuning — and it did, silently: an 8×8
//! `fc_awgr` read 0 instead of 1.109 while converging happily. `N ≤ 5` hid it.
//!
//! # Why not a unit convention
//!
//! Respelling λ in µm is the obvious fix and it is not enough. It repairs the
//! absolute term (1 µm → 1 pm) but the relative term `reltol·|λ|` is
//! **scale-invariant** — 1e-3 × 1.55 µm is the same 1.55 nm as 1e-3 × 1.55e-6 m.
//! That is still 5× a 40 GHz passband. A unit change cannot fix a relative
//! tolerance, and it would need a `×1e-6` at every consumer of the phase law
//! `φ = 2π·n_eff·L/λ`, which mixes `length_m` with the λ wire.
//!
//! So λ rows get their own absolute tolerance and **no relative term at all**:
//! λ is a label whose absolute precision is what matters.
//!
//! Exempting λ rows from the test entirely was the other candidate and is
//! unsafe — λ sits inside the feedback loop (`LambdaSelect` latching, ring
//! detuning, AWGR routing), so exempting it lets Newton declare convergence
//! while λ is still moving.
//!
//! # What is deliberately left alone
//!
//! Optical field amplitudes (√W, O(0.03)) and device-internal rows (series-R
//! nodes, thermal-RC nodes in kelvin) still use `vntol`. For those the volt
//! tolerance is *tighter* than they need rather than looser, which costs at
//! most an iteration and cannot produce a wrong answer. They are listed here so
//! the omission is a decision and not an oversight.

use fairchild_parser::Netlist;

use crate::mna::CircuitTopology;
use crate::options::SimOptions;

/// The `(absolute, relative)` tolerance pair for each MNA row.
///
/// Build once per solve, after `CircuitTopology::allocate_extra_rows` has
/// settled `topo.size`, then use [`Tolerances::converged`] for the Newton step
/// test and [`Tolerances::bound`] for anything else that needs the same scale
/// (the variable-step LTE norm divides by it).
#[derive(Debug, Clone)]
pub struct Tolerances {
    per_row: Vec<(f64, f64)>,
}

impl Tolerances {
    /// Classify every row of `topo` by the physical quantity it carries.
    ///
    /// Layout, as established by `CircuitTopology`: `[0, n_nodes)` node
    /// voltages, `[n_nodes, vsrc_end)` voltage-source branch currents,
    /// `[vsrc_end, size)` device-internal rows.
    pub fn build(netlist: &Netlist, topo: &CircuitTopology, opts: &SimOptions) -> Self {
        let mut per_row = vec![(opts.vntol, opts.reltol); topo.size];

        // Branch currents are amps, not volts — this is what `abstol` is for.
        let n_nodes = topo.n_nodes();
        for row in per_row.iter_mut().skip(n_nodes).take(topo.vsrc_index.len()) {
            *row = (opts.abstol, opts.reltol);
        }

        // λ wires are metres, and absolute-only (see the module docs). Gated on
        // `optical_nets` membership so an electrical net that happens to be
        // named `foo_wl_0` keeps the voltage tolerance.
        for net in &netlist.optical_nets {
            if !fairchild_parser::is_lambda_wire(net) {
                continue;
            }
            if let Some(&row) = topo.node_index.get(net) {
                per_row[row] = (opts.lambdatol, 0.0);
            }
        }

        Tolerances { per_row }
    }

    /// The convergence bound for row `i` at value `x`: `abs + rel·|x|`.
    pub fn bound(&self, i: usize, x: f64) -> f64 {
        let (abs, rel) = self.per_row[i];
        abs + rel * x.abs()
    }

    /// Did every row's Newton step land inside its own bound?
    ///
    /// The bound is taken at the *new* iterate, matching the test this replaces.
    pub fn converged(&self, x_new: &[f64], x_old: &[f64]) -> bool {
        debug_assert_eq!(
            self.per_row.len(),
            x_new.len(),
            "Tolerances built for a different matrix size — build it after \
             allocate_extra_rows has settled topo.size"
        );
        x_new
            .iter()
            .zip(x_old.iter())
            .enumerate()
            .all(|(i, (n, o))| (n - o).abs() < self.bound(i, *n))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

    fn deck() -> Netlist {
        parse_spice(
            "* one optical port and one voltage source\n\
             .optical_port p\n\
             Xl p fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
             V1 a 0 1\nR1 a 0 1k\n.op\n.end\n",
        )
        .unwrap()
    }

    #[test]
    fn each_row_class_gets_its_own_tolerance() {
        let net = deck();
        let topo = CircuitTopology::build(&net);
        let opts = SimOptions::default();
        let tol = Tolerances::build(&net, &topo, &opts);

        // λ wires: absolute-only, at lambdatol.
        let lam = topo.node_index["p_wl_0"];
        assert_eq!(tol.bound(lam, 1.55e-6), opts.lambdatol);
        assert_eq!(tol.bound(lam, 1e9), opts.lambdatol, "no relative term on λ");

        // Optical field amplitudes keep the voltage tolerance (documented).
        let re = topo.node_index["p_re_0"];
        assert_eq!(tol.bound(re, 0.0), opts.vntol);

        // Plain node: volts.
        let a = topo.node_index["a"];
        assert_eq!(tol.bound(a, 0.0), opts.vntol);

        // Branch current: amps.
        let i_v1 = topo.n_nodes() + topo.vsrc_index["v1"];
        assert_eq!(tol.bound(i_v1, 0.0), opts.abstol);
    }

    /// The bug this exists to prevent: a λ step of 10 pm used to pass.
    #[test]
    fn a_ten_picometre_lambda_step_is_not_converged() {
        let net = deck();
        let topo = CircuitTopology::build(&net);
        let opts = SimOptions::default();
        let tol = Tolerances::build(&net, &topo, &opts);

        let lam = topo.node_index["p_wl_0"];
        let mut x_old = vec![0.0; topo.size];
        let mut x_new = vec![0.0; topo.size];
        x_old[lam] = 1.55e-6;
        x_new[lam] = 1.55e-6 + 10e-12;

        assert!(
            !tol.converged(&x_new, &x_old),
            "10 pm of λ movement must not count as converged"
        );
        // The old uniform test accepted it, which is the whole point.
        let lam_new = x_new[lam];
        assert!((lam_new - x_old[lam]).abs() < opts.vntol + opts.reltol * lam_new.abs());

        // And it does converge once the step is genuinely small.
        x_new[lam] = 1.55e-6 + 1e-14;
        assert!(tol.converged(&x_new, &x_old));
    }
}
