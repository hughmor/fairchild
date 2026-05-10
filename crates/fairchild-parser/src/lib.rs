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
    VoltageSource {
        name: String,
        pos: NodeName,
        neg: NodeName,
        dc: f64,
    },
    CurrentSource {
        name: String,
        pos: NodeName,
        neg: NodeName,
        /// Positive convention: current flows from neg to pos (into pos node).
        dc: f64,
    },
}

/// A requested simulation analysis.
#[derive(Debug, Clone)]
pub enum Analysis {
    Op,
}
