//! `.measure` post-processor.
//!
//! Evaluates SPICE-style measurement directives against a completed
//! `TranResult`.  Each `Measurement` becomes a `(name, value)` scalar that the
//! CLI prints and Python exposes via `SimResult.measurements()`.

use fairchild_parser::{EvalContext, Expr, MeasKind, MeasOp, Measurement};

use crate::tran::TranResult;

/// Result of one `.measure` directive.
pub struct MeasureResult {
    pub name: String,
    pub value: f64,
}

/// Evaluate all `Measurement`s against a transient run.
///
/// Only `.meas tran` directives are evaluated; `.meas ac` and `.meas dc`
/// directives emit a warning and are skipped (AC and DC measurements require
/// the corresponding result type and are not yet implemented).
pub fn evaluate_measurements(
    measurements: &[Measurement],
    result: &TranResult,
) -> Vec<MeasureResult> {
    use fairchild_parser::MeasAnalysis;
    measurements
        .iter()
        .filter(|m| {
            if !matches!(m.analysis, MeasAnalysis::Tran) {
                eprintln!(
                    "warning: .meas '{}' is not a transient measurement and will be skipped \
                     — only `.meas tran` is currently supported",
                    m.name
                );
                return false;
            }
            true
        })
        .map(|m| {
            let v = eval_one(&m.kind, result);
            MeasureResult {
                name: m.name.clone(),
                value: v,
            }
        })
        .collect()
}

fn eval_one(kind: &MeasKind, r: &TranResult) -> f64 {
    match kind {
        MeasKind::FindAt { expr, at } => {
            let ctx = TranSample::at(r, *at);
            expr.eval(&ctx)
        }
        MeasKind::FindWhen { expr, cond, cross } => {
            let mut last_cond = f64::NAN;
            let mut found = 0usize;
            for i in 0..r.time.len() {
                let ctx = TranSample::index(r, i);
                let c = cond.eval(&ctx);
                if !last_cond.is_nan() && cross_zero(last_cond, c) {
                    found += 1;
                    if found >= *cross {
                        // Linear interpolation to the zero-crossing for accuracy.
                        let t0 = r.time[i - 1];
                        let t1 = r.time[i];
                        let alpha = last_cond / (last_cond - c);
                        let t_cross = t0 + alpha * (t1 - t0);
                        let ctx = TranSample::at(r, t_cross);
                        return expr.eval(&ctx);
                    }
                }
                last_cond = c;
            }
            f64::NAN
        }
        MeasKind::DerivAt { expr, at } => {
            let h = (r.time.last().copied().unwrap_or(0.0)
                - r.time.first().copied().unwrap_or(0.0))
                * 1e-6
                + 1e-12;
            let a = TranSample::at(r, at - h).eval_via(expr);
            let b = TranSample::at(r, at + h).eval_via(expr);
            (b - a) / (2.0 * h)
        }
        MeasKind::TrigTarg {
            trig_expr,
            trig_val,
            trig_cross,
            targ_expr,
            targ_val,
            targ_cross,
        } => {
            let t1 = find_cross(r, trig_expr, *trig_val, *trig_cross);
            let t2 = find_cross(r, targ_expr, *targ_val, *targ_cross);
            match (t1, t2) {
                (Some(a), Some(b)) => b - a,
                _ => f64::NAN,
            }
        }
        MeasKind::Aggregate { op, expr, from, to } => {
            // Filter time-slice [from, to].  Default to full range.
            let t0 = from.unwrap_or(r.time.first().copied().unwrap_or(0.0));
            let t1 = to.unwrap_or(r.time.last().copied().unwrap_or(0.0));
            let samples: Vec<(f64, f64)> = r
                .time
                .iter()
                .enumerate()
                .filter(|(_, &t)| t >= t0 && t <= t1)
                .map(|(i, &t)| (t, expr.eval(&TranSample::index(r, i))))
                .collect();
            if samples.is_empty() {
                return f64::NAN;
            }
            match op {
                MeasOp::Max => samples
                    .iter()
                    .map(|(_, v)| *v)
                    .fold(f64::NEG_INFINITY, f64::max),
                MeasOp::Min => samples
                    .iter()
                    .map(|(_, v)| *v)
                    .fold(f64::INFINITY, f64::min),
                MeasOp::Pp => {
                    let mx = samples
                        .iter()
                        .map(|(_, v)| *v)
                        .fold(f64::NEG_INFINITY, f64::max);
                    let mn = samples
                        .iter()
                        .map(|(_, v)| *v)
                        .fold(f64::INFINITY, f64::min);
                    mx - mn
                }
                MeasOp::Avg => {
                    // Time-weighted (trapezoidal) average.
                    if samples.len() < 2 {
                        return samples[0].1;
                    }
                    let span = samples.last().unwrap().0 - samples.first().unwrap().0;
                    if span <= 0.0 {
                        return samples[0].1;
                    }
                    let mut sum = 0.0;
                    for w in samples.windows(2) {
                        let (ta, va) = w[0];
                        let (tb, vb) = w[1];
                        sum += 0.5 * (va + vb) * (tb - ta);
                    }
                    sum / span
                }
                MeasOp::Rms => {
                    if samples.len() < 2 {
                        return samples[0].1.abs();
                    }
                    let span = samples.last().unwrap().0 - samples.first().unwrap().0;
                    if span <= 0.0 {
                        return samples[0].1.abs();
                    }
                    let mut sum = 0.0;
                    for w in samples.windows(2) {
                        let (ta, va) = w[0];
                        let (tb, vb) = w[1];
                        sum += 0.5 * (va * va + vb * vb) * (tb - ta);
                    }
                    (sum / span).sqrt()
                }
                MeasOp::Integ => {
                    // Trapezoidal ∫ expr dt.
                    if samples.len() < 2 {
                        return 0.0;
                    }
                    let mut sum = 0.0;
                    for w in samples.windows(2) {
                        let (ta, va) = w[0];
                        let (tb, vb) = w[1];
                        sum += 0.5 * (va + vb) * (tb - ta);
                    }
                    sum
                }
            }
        }
    }
}

fn cross_zero(a: f64, b: f64) -> bool {
    (a <= 0.0 && b > 0.0) || (a > 0.0 && b <= 0.0)
}

fn find_cross(r: &TranResult, expr: &Expr, val: f64, cross: usize) -> Option<f64> {
    let mut last = f64::NAN;
    let mut count = 0usize;
    for i in 0..r.time.len() {
        let ctx = TranSample::index(r, i);
        let v = expr.eval(&ctx) - val;
        if !last.is_nan() && cross_zero(last, v) {
            count += 1;
            if count >= cross {
                let t0 = r.time[i - 1];
                let t1 = r.time[i];
                let alpha = last / (last - v);
                return Some(t0 + alpha * (t1 - t0));
            }
        }
        last = v;
    }
    None
}

/// EvalContext snapshot over a transient result.
struct TranSample<'a> {
    r: &'a TranResult,
    t: f64,
    /// When set, look up the i'th stored sample directly (faster than
    /// linear-interpolating via voltage_at).
    direct_index: Option<usize>,
}

impl<'a> TranSample<'a> {
    fn at(r: &'a TranResult, t: f64) -> Self {
        TranSample {
            r,
            t,
            direct_index: None,
        }
    }
    fn index(r: &'a TranResult, i: usize) -> Self {
        TranSample {
            r,
            t: r.time[i],
            direct_index: Some(i),
        }
    }
    /// Convenience: evaluate `expr` against this sample.
    fn eval_via(self, expr: &Expr) -> f64 {
        expr.eval(&self)
    }
}

impl<'a> EvalContext for TranSample<'a> {
    fn node_voltage(&self, node: &str) -> f64 {
        if node == "0" || node == "gnd" {
            return 0.0;
        }
        if let Some(i) = self.direct_index {
            return self.r.node_voltages.get(node).map(|s| s[i]).unwrap_or(0.0);
        }
        self.r.voltage_at(node, self.t).unwrap_or(0.0)
    }
    fn branch_current(&self, vsrc: &str) -> f64 {
        if let Some(i) = self.direct_index {
            return self.r.vsrc_currents.get(vsrc).map(|s| s[i]).unwrap_or(0.0);
        }
        self.r.isrc_at(vsrc, self.t).unwrap_or(0.0)
    }
    fn time(&self) -> f64 {
        self.t
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

    #[test]
    fn measure_max_on_pulse() {
        let net = parse_spice(
            "* meas\nV1 in 0 PULSE(0 1 1u 1n 1n 5u 10u)\nR1 in out 1k\nC1 out 0 1u\n\
             .meas tran vmx MAX V(out)\n\
             .tran 100n 8u\n.end\n",
        )
        .unwrap();
        let r = crate::tran_nr_var(&net, 100e-9, 8e-6).unwrap();
        let ms = evaluate_measurements(&net.measurements, &r);
        assert_eq!(ms.len(), 1);
        assert_eq!(ms[0].name, "vmx");
        // Charge is small (τ=1ms ≫ 8µs), so V(out) max stays small but is the
        // peak of the rise.  Just check it's positive and ≤ 1 V.
        assert!(
            ms[0].value > 0.0 && ms[0].value < 1.0,
            "vmx={}",
            ms[0].value
        );
    }

    #[test]
    fn measure_find_at() {
        let net = parse_spice(
            "* meas\nV1 in 0 DC 1\nR1 in out 1k\nC1 out 0 1u\n\
             .meas tran vat FIND V(out) AT=1m\n\
             .tran 1u 5m\n.end\n",
        )
        .unwrap();
        let r = crate::tran_nr_var(&net, 1e-6, 5e-3).unwrap();
        let ms = evaluate_measurements(&net.measurements, &r);
        // V(out)(t=1ms) for RC with τ=1ms starting from DC OP V(out)=1V is constant 1V.
        // Actually with V1=DC 1, the DC OP gives V(out)=1 already so V(out)=1 for all t.
        assert!((ms[0].value - 1.0).abs() < 1e-3, "vat={}", ms[0].value);
    }

    #[test]
    fn measure_trig_targ_delay() {
        // PULSE V1 rises from 0 to 1 starting at t=1µs.  V(in) crosses 0.1V
        // somewhere on the rise; V(out) crosses 0.1V slightly later (RC).
        let net = parse_spice(
            "* delay\nV1 in 0 PULSE(0 1 1u 100n 100n 5u 10u)\n\
             R1 in out 1k\nC1 out 0 1n\n\
             .meas tran tpd TRIG V(in) VAL=0.5 TARG V(out) VAL=0.5\n\
             .tran 10n 3u\n.end\n",
        )
        .unwrap();
        let r = crate::tran_nr_var(&net, 10e-9, 3e-6).unwrap();
        let ms = evaluate_measurements(&net.measurements, &r);
        // Delay should be positive (V(out) lags V(in)).
        assert!(
            ms[0].value > 0.0 && ms[0].value < 1e-6,
            "tpd={} (should be in 0…1µs)",
            ms[0].value
        );
    }
}
