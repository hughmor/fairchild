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
}

impl Waveform {
    /// Value used for DC operating-point (t = 0).
    pub fn dc_value(&self) -> f64 {
        match self {
            Waveform::Dc(v) => *v,
            Waveform::Pulse { v0, .. } => *v0,
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
}

/// A model card parsed from `.model <name> <kind> [param=value ...]`.
#[derive(Debug, Clone)]
pub struct ModelCard {
    pub name: String,
    pub kind: String,
    pub params: Vec<(String, f64)>,
}

/// A requested simulation analysis.
#[derive(Debug, Clone)]
pub enum Analysis {
    Op,
    Tran { step: f64, stop: f64 },
}
