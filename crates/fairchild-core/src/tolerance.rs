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
//! | thermal | kelvin | 1 µK — eight digits on a temperature, nobody's ask |
//!
//! # The row class that is gone
//!
//! A λ wire carried **metres** (~1.55e-6), where `vntol` let Newton stop a
//! micron out — an 8×8 `fc_awgr` read 0 instead of 1.109 while converging
//! happily, and `N ≤ 5` hid it. That earned λ its own absolute-only
//! `lambdatol`, and respelling λ in µm was ruled out because `reltol·|λ|` is
//! scale-invariant and no unit choice can shrink it.
//!
//! λ is not a row any more (see [`crate::lambda`]) — it is a label resolved
//! before the solve — so the class and the option are both gone. This is
//! recorded rather than deleted because the *shape* of the bug recurs: any
//! unknown whose unit is not volts and whose precision matters needs its own
//! row class here, and the failure mode is a converged wrong answer.
//!
//! # The row class that arrived
//!
//! A `thermal` node carries a temperature rise in kelvin and its flow is watts.
//! It takes `temptol` (1 mK). Which rows those are is not declared in the deck:
//! OSDI carries a Verilog-A discipline's units through to the descriptor, so
//! [`crate::device::Device::thermal_nodes`] reads it off the model and
//! `push_device` records the rows on the topology. One statement, in the file
//! that already had to be right.
//!
//! Unlike λ, this class is not about a *wrong* answer from a loose bound —
//! `vntol` on kelvin is tighter than needed. It is about the other half of the
//! same unit-mixing: `vmax` is a trust region in **volts**, and a thermal row
//! allowed to set it clamps a 233 K step to 0.5 and scales every electrical
//! unknown by 1/466. Newton then meets its step test on the clamped deltas
//! rather than on the residual — 10× the iterations and a converged answer that
//! is measurably off. `newton.rs` excludes thermal rows from setting the clamp
//! for that reason; the shape is λ's, and it recurred within one release.
//!
//! # What is deliberately left alone
//!
//! Optical field amplitudes (√W, O(0.03)) and device-internal rows that are not
//! thermal (series-R nodes, OSDI flow branches) still use `vntol`. For those the
//! volt tolerance is *tighter* than they need rather than looser, which costs at
//! most an iteration and cannot produce a wrong answer. They are listed here so
//! the omission is a decision and not an oversight.

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
    pub fn build(topo: &CircuitTopology, opts: &SimOptions) -> Self {
        let mut per_row = vec![(opts.vntol, opts.reltol); topo.size];

        // Branch currents are amps, not volts — this is what `abstol` is for.
        let n_nodes = topo.n_nodes();
        for row in per_row.iter_mut().skip(n_nodes).take(topo.vsrc_index.len()) {
            *row = (opts.abstol, opts.reltol);
        }

        // Thermal rows carry kelvin. `vntol` is a microvolt, which as a bound
        // on a temperature is a demand for eight digits nobody asked for, and
        // the same unit-mixing that made `vntol` meaningless against a 1.55e-6
        // wavelength. Which rows these are comes from the models themselves —
        // see `Device::thermal_nodes` — so a deck states nothing and cannot
        // disagree.
        for &row in &topo.thermal_rows {
            if row < per_row.len() {
                per_row[row] = (opts.temptol, opts.reltol);
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

    #[test]
    fn each_row_class_gets_its_own_tolerance() {
        let net = parse_spice("* two classes\nV1 a 0 1\nR1 a 0 1k\n.op\n").unwrap();
        let topo = CircuitTopology::build(&net);
        let opts = SimOptions::default();
        let tol = Tolerances::build(&topo, &opts);

        // Plain node: volts.
        let a = topo.node_index["a"];
        assert_eq!(tol.bound(a, 0.0), opts.vntol);

        // Branch current: amps.
        let i_v1 = topo.n_nodes() + topo.vsrc_index["v1"];
        assert_eq!(tol.bound(i_v1, 0.0), opts.abstol);

        // And the relative term is on both, at reltol.
        assert_eq!(tol.bound(a, 2.0), opts.vntol + 2.0 * opts.reltol);
        assert_eq!(tol.bound(i_v1, 2.0), opts.abstol + 2.0 * opts.reltol);
    }

    /// A λ wire is no longer a row at all, so there is no λ tolerance to get
    /// wrong. Pinned here because the bug this module was written for was a
    /// *converged* wrong answer, and "the row is gone" is the only remaining
    /// reason it cannot come back.
    #[test]
    fn a_lambda_wire_is_not_a_row_to_give_a_tolerance_to() {
        let net = parse_spice(
            "* one optical port\n\
             .optical_port p\n\
             Xl p fc_cw_laser power_mW=1.0 wavelength_nm=1550\n.op\n",
        )
        .unwrap();
        let reg = crate::device_registry::DeviceRegistry::new();
        let ctx = SimOptions::default().sim_context();
        let topo = CircuitTopology::build_resolved(&net, &ctx, &reg);
        assert!(!topo.node_index.contains_key("p_wl_0"));
        assert!(topo.node_index.contains_key("p_re_0"));
    }
}
