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
}
