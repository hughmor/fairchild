use thiserror::Error;

#[derive(Debug, Error)]
pub enum ParseError {
    #[error("line {line}: {msg}")]
    Syntax { line: usize, msg: String },
    #[error("unknown element type '{letter}' on line {line}")]
    UnknownElement { letter: char, line: usize },
    #[error("expected {expected} fields, got {got} on line {line}")]
    FieldCount { expected: &'static str, got: usize, line: usize },
    #[error("invalid number '{value}' on line {line}: {source}")]
    BadNumber { value: String, line: usize, source: std::num::ParseFloatError },
}

/// Discipline mismatch detected during elaboration.
#[derive(Debug, Error, PartialEq, Eq)]
#[error("discipline mismatch: {element} connects electrical element to optical net '{net}'")]
pub struct DisciplineError {
    pub element: String,
    pub net: String,
}
