//! Voltage- and current-controlled switches — SPICE `S` and `W`.
//!
//! Both are the same device with a different thing to look at: a resistor whose
//! value is `RON` or `ROFF` depending on whether a control quantity is above or
//! below a threshold, with an optional hysteresis band.
//!
//! ```text
//! S<name> N+ N- NC+ NC- <model> [ON|OFF]   .model <m> SW  (VT= VH= RON= ROFF=)
//! W<name> N+ N- <vsource> <model> [ON|OFF] .model <m> CSW (IT= IH= RON= ROFF=)
//! ```
//!
//! # The switching law, and why it is a hard step
//!
//! Measured against ngspice 46, which is the reference fairchild validates
//! against:
//!
//! ```text
//! ctrl > threshold + hysteresis  → ON
//! ctrl < threshold − hysteresis  → OFF
//! otherwise                      → hold the previous state
//! ```
//!
//! There is no smoothing: at `VT=1, VH=0` ngspice is still OFF at exactly 1.0 V
//! and ON at 1.1 V. That discontinuity is the model, and it is why SPICE
//! switches have a reputation for hurting convergence — a conductance that
//! jumps 12 orders of magnitude between Newton iterations is exactly what
//! Newton is not built for. `VH` is the user's tool against it, and it is worth
//! reaching for before `itl1`.
//!
//! # The hysteresis reference is the *committed* state
//!
//! "The previous state" is deliberately the state at the last accepted
//! timepoint, not the last Newton iterate. Using the iterate would let the
//! switch flip-flop within one NR loop — the same failure the photonic
//! `LambdaSelect` latch exists to prevent. Outside the band the reference does
//! not matter (the control decides), so this only fixes the ambiguous case, and
//! it fixes it deterministically.
//!
//! # Known divergence from ngspice: DC sweeps
//!
//! In ngspice a `.dc` sweep carries switch state from point to point, so
//! sweeping up and sweeping down through the band give different answers —
//! genuine hysteresis. fairchild runs sweep points **in parallel** (rayon), so
//! there is no "previous point" to inherit from and every point starts at the
//! instance's `ON`/`OFF` keyword. Transient runs, where hysteresis actually
//! matters, are unaffected: state is latched at each accepted timestep.
//!
//! With the default `VH = 0` there is no band and no difference at all.

use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

/// What the switch watches.
#[derive(Debug, Clone, Copy)]
pub enum SwitchControl {
    /// Differential voltage across two nodes — the `S` element.
    Voltage { pos: NodeId, neg: NodeId },
    /// Current through a voltage source's branch row — the `W` element.
    ///
    /// Held as an MNA row index rather than a terminal because that is what it
    /// is: `n_nodes + vsrc_index[name]`, resolved when the device is built. The
    /// row is only ever *read*, never stamped into, which is why it does not
    /// appear in [`Device::stamp_pairs`].
    Current { row: NodeId },
}

/// A SPICE `S` / `W` switch.
pub struct Switch {
    ron: f64,
    roff: f64,
    threshold: f64,
    hysteresis: f64,
    control: SwitchControl,

    pos: NodeId,
    neg: NodeId,

    /// State at the last accepted timepoint — the hysteresis reference.
    state: bool,
    /// State at the current Newton iterate; becomes `state` on commit.
    state_eval: bool,
    /// Conductance at the current iterate.
    g: f64,
}

impl Switch {
    /// Build from a `.model` card's parameters.
    ///
    /// `is_current` selects the `CSW` spelling of the threshold parameters
    /// (`IT`/`IH`) over the `SW` one (`VT`/`VH`); everything else is shared.
    /// Returns the unrecognised parameter names alongside the device, matching
    /// how the diode and MOSFET cards report theirs.
    pub fn from_model_params(
        is_current: bool,
        params: &[(String, f64)],
        initial_on: bool,
    ) -> Result<(Self, Vec<String>), String> {
        let mut ron = 1.0;
        let mut roff = 1e12;
        let mut threshold = 0.0;
        let mut hysteresis = 0.0;
        let mut unknown = Vec::new();
        for (k, v) in params {
            match k.to_lowercase().as_str() {
                "ron" => ron = *v,
                "roff" => roff = *v,
                "vt" | "it" | "vthreshold" | "ithreshold" => threshold = *v,
                "vh" | "ih" | "vhysteresis" | "ihysteresis" => hysteresis = *v,
                _ => unknown.push(k.clone()),
            }
        }
        // Rejected rather than clamped: `1/RON` is stamped straight into the
        // Jacobian, so a non-positive value is an infinite or negative
        // conductance — a wrong answer or a singular matrix, either way not
        // something to paper over.
        if ron <= 0.0 {
            return Err(format!("switch RON must be > 0, got {ron}"));
        }
        if roff <= 0.0 {
            return Err(format!("switch ROFF must be > 0, got {roff}"));
        }
        // Negative hysteresis would invert the band and make the switch
        // bistable in the wrong direction; ngspice takes the magnitude.
        let hysteresis = hysteresis.abs();
        let _ = is_current; // only the parameter spelling differs

        Ok((
            Switch {
                ron,
                roff,
                threshold,
                hysteresis,
                // Rebound by `setup_instance`.
                control: SwitchControl::Voltage {
                    pos: None,
                    neg: None,
                },
                pos: None,
                neg: None,
                state: initial_on,
                state_eval: initial_on,
                g: if initial_on { 1.0 / ron } else { 1.0 / roff },
            },
            unknown,
        ))
    }

    /// Point the switch at what it controls on. Called by the builder, which is
    /// the only place that knows a `W`'s controlling branch row.
    pub fn set_control(&mut self, control: SwitchControl) {
        self.control = control;
    }

    fn read(&self, x: &[f64]) -> f64 {
        let at = |n: NodeId| n.map_or(0.0, |i| x[i]);
        match self.control {
            SwitchControl::Voltage { pos, neg } => at(pos) - at(neg),
            SwitchControl::Current { row } => at(row),
        }
    }
}

impl Device for Switch {
    fn num_terminals(&self) -> usize {
        match self.control {
            SwitchControl::Voltage { .. } => 4,
            SwitchControl::Current { .. } => 2,
        }
    }

    fn setup_model(&mut self, _ctx: &SimContext) {}

    fn setup_instance(&mut self, terminals: &[NodeId], _ctx: &SimContext) {
        self.pos = terminals[0];
        self.neg = terminals[1];
        if let SwitchControl::Voltage { .. } = self.control {
            debug_assert_eq!(terminals.len(), 4, "S switch expects [N+, N-, NC+, NC-]");
            self.control = SwitchControl::Voltage {
                pos: terminals[2],
                neg: terminals[3],
            };
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        let ctrl = self.read(x);
        self.state_eval = if ctrl > self.threshold + self.hysteresis {
            true
        } else if ctrl < self.threshold - self.hysteresis {
            false
        } else {
            // Inside the band: hold the committed state, not the last iterate.
            self.state
        };
        self.g = if self.state_eval {
            1.0 / self.ron
        } else {
            1.0 / self.roff
        };
    }

    fn load_residual(&self, _b: &mut [f64]) {
        // A switch is a plain conductance: `i = g·v` with no offset, so the
        // Norton source is zero and the whole device lives in the Jacobian.
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        if let Some(p) = self.pos {
            mat.a[p][p] += self.g;
            if let Some(n) = self.neg {
                mat.a[p][n] -= self.g;
            }
        }
        if let Some(n) = self.neg {
            mat.a[n][n] += self.g;
            if let Some(p) = self.pos {
                mat.a[n][p] -= self.g;
            }
        }
    }

    /// Only the four (N+, N−) cells are ever written — notably *not* the
    /// control nodes, which are read-only here. Declaring that keeps a `W`'s
    /// controlling branch row out of the sparsity pattern.
    fn stamp_pairs(&self) -> Option<Vec<(usize, usize)>> {
        let mut pairs = Vec::with_capacity(4);
        for a in [self.pos, self.neg].into_iter().flatten() {
            for b in [self.pos, self.neg].into_iter().flatten() {
                pairs.push((a, b));
            }
        }
        Some(pairs)
    }

    fn commit_timestep(&mut self, _x: &[f64]) {
        self.state = self.state_eval;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sw(vt: f64, vh: f64, initial_on: bool) -> Switch {
        let params = [
            ("ron".to_string(), 10.0),
            ("roff".to_string(), 1e6),
            ("vt".to_string(), vt),
            ("vh".to_string(), vh),
        ];
        let (mut s, unknown) = Switch::from_model_params(false, &params, initial_on).unwrap();
        assert!(unknown.is_empty());
        s.setup_instance(
            &[Some(0), Some(1), Some(2), Some(3)],
            &SimContext::default(),
        );
        s
    }

    /// The law, pinned against the ngspice 46 measurements in the module docs.
    #[test]
    fn threshold_and_hysteresis_match_ngspice() {
        let ctx = SimContext::default();
        // VH = 0: OFF at exactly VT, ON just above.
        let mut s = sw(1.0, 0.0, false);
        for (v, want_on) in [(0.9, false), (1.0, false), (1.1, true)] {
            s.eval(&[0.0, 0.0, v, 0.0], EvalFlags::dc(), &ctx);
            assert_eq!(s.state_eval, want_on, "VH=0 at Vctrl={v}");
        }

        // VH = 0.5 starting OFF: needs to clear VT+VH before turning on.
        let mut s = sw(1.0, 0.5, false);
        for (v, want_on) in [(1.25, false), (1.5, false), (1.75, true)] {
            s.eval(&[0.0, 0.0, v, 0.0], EvalFlags::dc(), &ctx);
            assert_eq!(s.state_eval, want_on, "rising through the band at {v}");
        }

        // VH = 0.5 starting ON: stays on until it drops below VT−VH.
        let mut s = sw(1.0, 0.5, true);
        for (v, want_on) in [(0.75, true), (0.5, true), (0.25, false)] {
            s.eval(&[0.0, 0.0, v, 0.0], EvalFlags::dc(), &ctx);
            assert_eq!(s.state_eval, want_on, "falling through the band at {v}");
        }
    }

    /// Inside the band the answer must come from the *committed* state, so a
    /// Newton loop cannot make the switch chatter.
    #[test]
    fn the_band_holds_the_committed_state_not_the_last_iterate() {
        let ctx = SimContext::default();
        let mut s = sw(1.0, 0.5, false);

        // Drive it out of the band and back in without committing: the eval
        // state follows the control, but the reference has not moved.
        s.eval(&[0.0, 0.0, 2.0, 0.0], EvalFlags::dc(), &ctx);
        assert!(s.state_eval, "above the band");
        s.eval(&[0.0, 0.0, 1.0, 0.0], EvalFlags::dc(), &ctx);
        assert!(
            !s.state_eval,
            "back inside the band: still the committed OFF"
        );

        // Commit while on, and the band now holds ON.
        s.eval(&[0.0, 0.0, 2.0, 0.0], EvalFlags::dc(), &ctx);
        s.commit_timestep(&[]);
        s.eval(&[0.0, 0.0, 1.0, 0.0], EvalFlags::dc(), &ctx);
        assert!(s.state_eval, "band now holds the committed ON");
    }

    #[test]
    fn a_non_positive_on_resistance_is_rejected() {
        let bad = [("ron".to_string(), 0.0)];
        assert!(Switch::from_model_params(false, &bad, false).is_err());
        let bad = [("roff".to_string(), -1.0)];
        assert!(Switch::from_model_params(false, &bad, false).is_err());
    }

    /// The control nodes are read, never written — a `W`'s branch row must not
    /// end up in the sparsity pattern.
    #[test]
    fn stamp_pairs_covers_only_the_switched_branch() {
        let s = sw(1.0, 0.0, false);
        let pairs = s.stamp_pairs().unwrap();
        assert_eq!(pairs.len(), 4);
        assert!(pairs.iter().all(|&(r, c)| r < 2 && c < 2), "{pairs:?}");
    }
}
