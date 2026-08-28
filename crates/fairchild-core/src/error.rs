use thiserror::Error;

#[derive(Debug, Error)]
pub enum SimError {
    #[error("parse error: {0}")]
    Parse(#[from] fairchild_parser::ParseError),
    #[error(
        "singular matrix — the circuit has no unique DC solution. Common causes: \
         a node with no DC path to ground, voltage sources in parallel with \
         different values, or a loop of voltage sources / inductors. Run with \
         verbose to list structurally empty MNA rows by name."
    )]
    SingularMatrix,
    #[error("no .op analysis requested")]
    NoAnalysis,
    #[error("unknown node '{0}'")]
    UnknownNode(String),
    #[error("unknown model '{0}'")]
    UnknownModel(String),
    #[error("Newton-Raphson did not converge after {iters} iterations")]
    NoConvergence { iters: usize },
    #[error("floating node(s) detected: {} have no DC path to ground (R/L/V path). Add a series resistor or .nodeset.", nodes.join(", "))]
    FloatingNodes { nodes: Vec<String> },
    #[error("parameter error: {0}")]
    ParameterError(String),
    #[error(
        "model '{model}' is binned and no bin covers L={l:.4e} W={w:.4e}.\n\
         Bins:\n{windows}\n\
         Picking the nearest bin would be a wrong answer with nothing to read. \
         Either the instance geometry is outside what the PDK characterised, or a \
         bin card is missing from the deck."
    )]
    NoMatchingBin {
        model: String,
        l: f64,
        w: f64,
        windows: String,
    },
    #[error(
        "no AC source: .ac needs at least one source with an `AC <mag> [phase]` spec, \
         e.g. `V1 in 0 DC 0 AC 1`. Without one there is nothing to excite the circuit."
    )]
    NoAcSource,
    #[error(
        "'{first}' and '{second}' both drive node '{node}'. Two devices pinning one \
         potential leave the block rank-deficient, which the solve does not report — it \
         returns a weighted average of the two answers. On an optical bundle this is \
         almost always a port facing the wrong way: a device's `in` port drives the \
         backward wires and reads the forward ones, and its `out` port does the \
         opposite, so two `in` ports (or two `out` ports) wired together collide. \
         Check the port order on both elements."
    )]
    OverdrivenNode {
        node: String,
        first: String,
        second: String,
    },
}
