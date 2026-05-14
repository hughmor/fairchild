use thiserror::Error;

#[derive(Debug, Error)]
pub enum SimError {
    #[error("parse error: {0}")]
    Parse(#[from] fairchild_parser::ParseError),
    #[error("singular matrix — circuit may be floating or have no DC path to ground")]
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
}
