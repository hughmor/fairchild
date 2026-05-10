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
