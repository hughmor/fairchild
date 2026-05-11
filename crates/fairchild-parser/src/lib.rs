mod error;
mod spice;

pub use error::ParseError;
pub use spice::parse_spice;

/// A node name. "0" and "gnd" and "GND" all refer to ground.
pub type NodeName = String;

/// A parsed circuit netlist.
#[derive(Debug, Default)]
pub struct Netlist {
    pub title: String,
    pub elements: Vec<Element>,
    pub analyses: Vec<Analysis>,
    pub models: Vec<ModelCard>,
    /// Paths from `.osdi <path>` directives — OSDI shared libraries to load.
    pub osdi_paths: Vec<String>,
}

/// Waveform specification for independent sources.
#[derive(Debug, Clone)]
pub enum Waveform {
    Dc(f64),
    /// PULSE(v0 v1 td tr tf pw per)
    Pulse {
        v0: f64,
        v1: f64,
        td: f64,
        tr: f64,
        tf: f64,
        pw: f64,
        per: f64,
    },
    /// PWL(t0 v0 t1 v1 ...) — piecewise-linear; points must be sorted by time.
    Pwl {
        points: Vec<(f64, f64)>,
    },
}

impl Waveform {
    /// Value used for DC operating-point (t = 0).
    pub fn dc_value(&self) -> f64 {
        match self {
            Waveform::Dc(v) => *v,
            Waveform::Pulse { v0, .. } => *v0,
            Waveform::Pwl { points } => points.first().map(|(_, v)| *v).unwrap_or(0.0),
        }
    }

    /// Next time strictly after `t` at which this waveform has a slope discontinuity.
    ///
    /// Returns `None` for smooth (DC) waveforms or when all breakpoints are in the past.
    pub fn next_breakpoint(&self, t: f64) -> Option<f64> {
        match self {
            Waveform::Dc(_) => None,
            Waveform::Pulse { td, tr, tf, pw, per, .. } => {
                if t < *td {
                    return Some(*td);
                }
                // Offsets from the start of a period where slope changes.
                let offsets: [f64; 4] = [0.0, *tr, tr + pw, tr + pw + tf];
                if *per <= 0.0 {
                    return offsets.iter()
                        .map(|b| td + b)
                        .filter(|&bp| bp > t)
                        .reduce(f64::min);
                }
                let phase = (t - td) % per;
                let base = t - phase; // start of current period
                offsets.iter()
                    .map(|b| base + b)
                    .chain(std::iter::once(base + per))
                    .filter(|&bp| bp > t)
                    .reduce(f64::min)
            }
            Waveform::Pwl { points } => {
                points.iter().map(|(pt, _)| *pt).find(|&pt| pt > t)
            }
        }
    }

    /// Value at time t (seconds).
    pub fn at(&self, t: f64) -> f64 {
        match self {
            Waveform::Dc(v) => *v,
            Waveform::Pulse { v0, v1, td, tr, tf, pw, per } => {
                if t < *td {
                    return *v0;
                }
                // Time within the current period.
                let tp = if *per > 0.0 { (t - td) % per } else { t - td };
                if tp < *tr {
                    v0 + (v1 - v0) * tp / tr
                } else if tp < tr + pw {
                    *v1
                } else if tp < tr + pw + tf {
                    v1 + (v0 - v1) * (tp - tr - pw) / tf
                } else {
                    *v0
                }
            }
            Waveform::Pwl { points } => {
                if points.is_empty() {
                    return 0.0;
                }
                if t <= points[0].0 {
                    return points[0].1;
                }
                if t >= points[points.len() - 1].0 {
                    return points[points.len() - 1].1;
                }
                // Binary search for the segment containing t.
                let idx = points.partition_point(|(pt, _)| *pt <= t);
                let (t0, v0) = points[idx - 1];
                let (t1, v1) = points[idx];
                let frac = (t - t0) / (t1 - t0);
                v0 + (v1 - v0) * frac
            }
        }
    }
}

/// A single circuit element.
#[derive(Debug, Clone)]
pub enum Element {
    Resistor {
        name: String,
        pos: NodeName,
        neg: NodeName,
        resistance: f64,
    },
    Capacitor {
        name: String,
        pos: NodeName,
        neg: NodeName,
        capacitance: f64,
    },
    Inductor {
        name: String,
        pos: NodeName,
        neg: NodeName,
        inductance: f64,
    },
    VoltageSource {
        name: String,
        pos: NodeName,
        neg: NodeName,
        waveform: Waveform,
    },
    CurrentSource {
        name: String,
        pos: NodeName,
        neg: NodeName,
        waveform: Waveform,
    },
    Diode {
        name: String,
        anode: NodeName,
        cathode: NodeName,
        model_name: String,
    },
    Mosfet {
        name: String,
        drain: NodeName,
        gate: NodeName,
        source: NodeName,
        bulk: NodeName,
        model_name: String,
        params: Vec<(String, f64)>,
    },
}

/// A model card parsed from `.model <name> <kind> [param=value ...]`.
#[derive(Debug, Clone)]
pub struct ModelCard {
    pub name: String,
    pub kind: String,
    pub params: Vec<(String, f64)>,
}

/// Frequency spacing for AC sweep.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AcVariation {
    /// Points per decade (logarithmic).
    Dec,
    /// Points per octave (logarithmic).
    Oct,
    /// Total points (linear).
    Lin,
}

/// A requested simulation analysis.
#[derive(Debug, Clone)]
pub enum Analysis {
    Op,
    Tran { step: f64, stop: f64 },
    /// `.ac DEC|OCT|LIN <points> <fstart> <fstop>`
    Ac { variation: AcVariation, points: usize, fstart: f64, fstop: f64 },
}
