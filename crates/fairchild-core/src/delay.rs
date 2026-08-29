//! Reusable propagation-delay line for time-domain devices.
//!
//! Several devices need to reproduce a signal at one port as a time-delayed
//! copy of a signal at another port: an optical waveguide delays its envelope
//! by the group delay `τ_g = L·n_g/c`; an electrical lossless transmission line
//! delays its travelling waves by `TD`; a long PN phase shifter has the same
//! optical group delay as a plain waveguide. The *bookkeeping* for all of them
//! is identical — record a snapshot of the relevant port quantities at every
//! accepted timestep, then reconstruct them at `t − τ` by interpolation — even
//! though *what* is recorded and *how* the delayed value is stamped differ per
//! device.
//!
//! [`DelayLine`] owns the generic part (the history ring, the interpolation,
//! the trim window, and the active/time state). Each device keeps its own
//! [`DelayLine`] field and supplies the device-specific parts:
//!
//! - the delay `τ` (may be constant or state-dependent),
//! - a *snapshot* function that reads the port quantities to record from the
//!   solution vector (the layout is the device's business),
//! - the stamping that applies the reconstructed delayed snapshot to the
//!   residual / Jacobian (optical re/im rotation, Branin current sources, …).
//!
//! Usage pattern (per device):
//! ```ignore
//! // in eval():
//! let active = flags.transient && delay_option_on && tau > 0.0;
//! self.delay.set_state(active, ctx.time_s);
//! if active { self.delayed = self.delay.sample(tau, width); }
//! // in commit_timestep():
//! if self.delay.is_active() { self.delay.record(self.snapshot(x), tau); }
//! // in load_residual()/load_jacobian(): branch on self.delay.is_active()
//! ```

/// Generic history buffer + linear interpolation for a propagation delay.
///
/// Holds a monotonically-increasing list of timestamps and a parallel list of
/// fixed-width snapshots. The snapshot layout is opaque to `DelayLine` — it is
/// whatever the owning device records; only the width must be consistent.
#[derive(Clone, Debug, Default)]
pub struct DelayLine {
    /// True for the current eval iff the delay model is engaged.
    active: bool,
    /// Absolute time of the current eval (set by [`set_state`](Self::set_state)).
    time_s: f64,
    /// Committed timestamps (monotonic increasing).
    hist_t: Vec<f64>,
    /// Committed snapshots, one per timestamp, all the same width.
    hist_vals: Vec<Vec<f64>>,
}

impl DelayLine {
    pub fn new() -> Self {
        Self::default()
    }

    /// Set whether the delay model is engaged this eval, and record the current
    /// absolute time. Call once at the top of the device's `eval`.
    pub fn set_state(&mut self, active: bool, time_s: f64) {
        self.active = active;
        self.time_s = time_s;
    }

    /// Whether the delay model is engaged for the current step.
    pub fn is_active(&self) -> bool {
        self.active
    }

    /// Reconstruct the recorded snapshot at `time_s − delay_s` by linear
    /// interpolation, clamped to the endpoints.
    ///
    /// Clamping semantics: before the first recorded sample the signal is taken
    /// as the earliest snapshot (the DC/initial value); after the last it holds
    /// the most recent snapshot, so a delay shorter than the current timestep
    /// degrades gracefully to "no delay" rather than extrapolating. Returns a
    /// `width`-long zero vector if no history has been recorded yet.
    pub fn sample(&self, delay_s: f64, width: usize) -> Vec<f64> {
        let tq = self.time_s - delay_s;
        if self.hist_t.is_empty() {
            return vec![0.0; width];
        }
        if tq <= self.hist_t[0] {
            return self.hist_vals[0].clone();
        }
        let last = self.hist_t.len() - 1;
        if tq >= self.hist_t[last] {
            return self.hist_vals[last].clone();
        }
        // Binary search for the bracketing interval [i, i+1].
        let i = match self
            .hist_t
            .binary_search_by(|t| t.partial_cmp(&tq).unwrap_or(std::cmp::Ordering::Less))
        {
            Ok(j) => return self.hist_vals[j].clone(),
            Err(j) => j - 1,
        };
        let (t0, t1) = (self.hist_t[i], self.hist_t[i + 1]);
        let f = if t1 > t0 { (tq - t0) / (t1 - t0) } else { 0.0 };
        let (a, b) = (&self.hist_vals[i], &self.hist_vals[i + 1]);
        let w = a.len().min(b.len());
        (0..w).map(|j| a[j] + f * (b[j] - a[j])).collect()
    }

    /// Record `snapshot` at the current time (set via [`Self::set_state`]) and trim
    /// history older than one full delay window — keeping one sample before the
    /// window so [`sample`](Self::sample) can still bracket `t − delay_s`.
    pub fn record(&mut self, snapshot: Vec<f64>, delay_s: f64) {
        self.hist_t.push(self.time_s);
        self.hist_vals.push(snapshot);
        let cutoff = self.time_s - delay_s;
        let mut drop = 0;
        while drop + 1 < self.hist_t.len() && self.hist_t[drop + 1] <= cutoff {
            drop += 1;
        }
        if drop > 0 {
            self.hist_t.drain(0..drop);
            self.hist_vals.drain(0..drop);
        }
    }

    /// Forget all recorded history (e.g. between independent analysis runs).
    pub fn clear(&mut self) {
        self.hist_t.clear();
        self.hist_vals.clear();
    }

    /// Number of recorded samples (for diagnostics/tests).
    pub fn len(&self) -> usize {
        self.hist_t.len()
    }

    pub fn is_empty(&self) -> bool {
        self.hist_t.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_history_returns_zeros() {
        let mut d = DelayLine::new();
        d.set_state(true, 5.0);
        assert_eq!(d.sample(2.0, 3), vec![0.0, 0.0, 0.0]);
    }

    #[test]
    fn reconstructs_and_interpolates() {
        let mut d = DelayLine::new();
        // Record a ramp val[0] = t for t = 0,1,2,3 (delay window 2).
        for t in 0..=3 {
            d.set_state(true, t as f64);
            d.record(vec![t as f64], 2.0);
        }
        // At t=3, delay 2 → sample at t=1 → exactly 1.0.
        d.set_state(true, 3.0);
        assert!((d.sample(2.0, 1)[0] - 1.0).abs() < 1e-12);
        // At t=3.5, delay 2 → sample at t=1.5 → interpolated 1.5.
        d.set_state(true, 3.5);
        assert!((d.sample(2.0, 1)[0] - 1.5).abs() < 1e-12);
    }

    #[test]
    fn clamps_before_first_and_after_last() {
        let mut d = DelayLine::new();
        d.set_state(true, 1.0);
        d.record(vec![10.0], 5.0);
        d.set_state(true, 2.0);
        d.record(vec![20.0], 5.0);
        // Query well before first sample → earliest value.
        d.set_state(true, 1.0);
        assert_eq!(d.sample(10.0, 1)[0], 10.0);
        // Query at/after last (delay < step) → most recent value.
        d.set_state(true, 2.0);
        assert_eq!(d.sample(0.0, 1)[0], 20.0);
    }

    #[test]
    fn trims_to_delay_window() {
        let mut d = DelayLine::new();
        for t in 0..100 {
            d.set_state(true, t as f64);
            d.record(vec![t as f64], 3.0);
        }
        // With a 3 s window the buffer stays small (a handful of samples), not 100.
        assert!(d.len() <= 6, "history not trimmed: {}", d.len());
        // And it still reconstructs correctly within the window.
        d.set_state(true, 99.0);
        assert!((d.sample(3.0, 1)[0] - 96.0).abs() < 1e-12);
    }
}
