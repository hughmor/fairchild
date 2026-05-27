use super::common::parse_value;
use crate::{ParseError, Waveform};

pub(super) fn parse_waveform(tokens: &[&str], lineno: usize) -> Result<Waveform, ParseError> {
    if tokens.len() < 4 {
        return Err(ParseError::FieldCount {
            expected: "≥4",
            got: tokens.len(),
            line: lineno,
        });
    }
    let rest = tokens[3..].join(" ");
    let rest_lc = rest.to_lowercase();

    if rest_lc.starts_with("pulse") {
        return parse_pulse(&rest_lc, lineno);
    }
    if rest_lc.starts_with("pwl") {
        return parse_pwl(&rest_lc, lineno);
    }
    if rest_lc.starts_with("sin") {
        return parse_sin(&rest_lc, lineno);
    }
    if rest_lc.starts_with("exp") {
        return parse_exp(&rest_lc, lineno);
    }
    if rest_lc.starts_with("sffm") {
        return parse_sffm(&rest_lc, lineno);
    }
    if rest_lc.starts_with("am") {
        return parse_am(&rest_lc, lineno);
    }

    let tok = tokens[3].to_lowercase();
    if tok == "dc" {
        if tokens.len() < 5 {
            return Err(ParseError::FieldCount {
                expected: "≥5 (DC keyword)",
                got: tokens.len(),
                line: lineno,
            });
        }
        Ok(Waveform::Dc(parse_value(tokens[4], lineno)?))
    } else {
        Ok(Waveform::Dc(parse_value(tokens[3], lineno)?))
    }
}

/// Extract the parenthesised parameter list `(p0 p1 ...)` after a waveform
/// keyword, returning the tokens.  The keyword (e.g. `sin`, `exp`) precedes
/// the parens; trailing tokens after `)` are discarded.
pub(super) fn parens_tokens<'a>(
    s: &'a str,
    kind: &str,
    lineno: usize,
) -> Result<Vec<&'a str>, ParseError> {
    let start = s.find('(').ok_or_else(|| ParseError::Syntax {
        line: lineno,
        msg: format!("{}: missing '('", kind.to_uppercase()),
    })?;
    let end = s.rfind(')').ok_or_else(|| ParseError::Syntax {
        line: lineno,
        msg: format!("{}: missing ')'", kind.to_uppercase()),
    })?;
    Ok(s[start + 1..end].split_whitespace().collect())
}

pub(super) fn parse_sin(s: &str, lineno: usize) -> Result<Waveform, ParseError> {
    let parts = parens_tokens(s, "sin", lineno)?;
    let get = |i: usize, default: f64| -> Result<f64, ParseError> {
        parts.get(i).map_or(Ok(default), |t| parse_value(t, lineno))
    };
    if parts.len() < 3 {
        return Err(ParseError::FieldCount {
            expected: "≥3 (SIN vo va freq …)",
            got: parts.len(),
            line: lineno,
        });
    }
    Ok(Waveform::Sin {
        vo: get(0, 0.0)?,
        va: get(1, 0.0)?,
        freq: get(2, 0.0)?,
        td: get(3, 0.0)?,
        theta: get(4, 0.0)?,
        phase: get(5, 0.0)?,
    })
}

pub(super) fn parse_exp(s: &str, lineno: usize) -> Result<Waveform, ParseError> {
    let parts = parens_tokens(s, "exp", lineno)?;
    let get = |i: usize, default: f64| -> Result<f64, ParseError> {
        parts.get(i).map_or(Ok(default), |t| parse_value(t, lineno))
    };
    if parts.len() < 2 {
        return Err(ParseError::FieldCount {
            expected: "≥2 (EXP v1 v2 …)",
            got: parts.len(),
            line: lineno,
        });
    }
    // ngspice defaults: td1=0, tau1=tstep, td2=td1+tstep, tau2=tstep.  We don't
    // know tstep yet at parse time so substitute 0 / a small positive number;
    // the user should pass them explicitly.
    let v1 = get(0, 0.0)?;
    let v2 = get(1, 0.0)?;
    let td1 = get(2, 0.0)?;
    let tau1 = get(3, 1e-9)?;
    let td2 = get(4, td1 + tau1)?;
    let tau2 = get(5, tau1)?;
    Ok(Waveform::Exp {
        v1,
        v2,
        td1,
        tau1,
        td2,
        tau2,
    })
}

pub(super) fn parse_sffm(s: &str, lineno: usize) -> Result<Waveform, ParseError> {
    let parts = parens_tokens(s, "sffm", lineno)?;
    let get = |i: usize, default: f64| -> Result<f64, ParseError> {
        parts.get(i).map_or(Ok(default), |t| parse_value(t, lineno))
    };
    if parts.len() < 5 {
        return Err(ParseError::FieldCount {
            expected: "5 (SFFM vo va fc mdi fs)",
            got: parts.len(),
            line: lineno,
        });
    }
    Ok(Waveform::Sffm {
        vo: get(0, 0.0)?,
        va: get(1, 0.0)?,
        fc: get(2, 0.0)?,
        mdi: get(3, 0.0)?,
        fs: get(4, 0.0)?,
    })
}

pub(super) fn parse_am(s: &str, lineno: usize) -> Result<Waveform, ParseError> {
    let parts = parens_tokens(s, "am", lineno)?;
    let get = |i: usize, default: f64| -> Result<f64, ParseError> {
        parts.get(i).map_or(Ok(default), |t| parse_value(t, lineno))
    };
    if parts.len() < 4 {
        return Err(ParseError::FieldCount {
            expected: "≥4 (AM va vo mf fc [td])",
            got: parts.len(),
            line: lineno,
        });
    }
    Ok(Waveform::Am {
        va: get(0, 0.0)?,
        vo: get(1, 0.0)?,
        mf: get(2, 0.0)?,
        fc: get(3, 0.0)?,
        td: get(4, 0.0)?,
    })
}

pub(super) fn parse_pulse(s: &str, lineno: usize) -> Result<Waveform, ParseError> {
    let start = s.find('(').ok_or_else(|| ParseError::Syntax {
        line: lineno,
        msg: "PULSE: missing '('".into(),
    })?;
    let end = s.rfind(')').ok_or_else(|| ParseError::Syntax {
        line: lineno,
        msg: "PULSE: missing ')'".into(),
    })?;
    let inner = &s[start + 1..end];
    let parts: Vec<&str> = inner.split_whitespace().collect();
    let get = |i: usize, default: f64| -> Result<f64, ParseError> {
        parts.get(i).map_or(Ok(default), |s| parse_value(s, lineno))
    };
    if parts.len() < 2 {
        return Err(ParseError::FieldCount {
            expected: "≥2 (PULSE v0 v1 ...)",
            got: parts.len(),
            line: lineno,
        });
    }
    Ok(Waveform::Pulse {
        v0: get(0, 0.0)?,
        v1: get(1, 0.0)?,
        td: get(2, 0.0)?,
        tr: get(3, 0.0)?,
        tf: get(4, 0.0)?,
        pw: get(5, f64::INFINITY)?,
        per: get(6, f64::INFINITY)?,
    })
}

pub(super) fn parse_pwl(s: &str, lineno: usize) -> Result<Waveform, ParseError> {
    let inner = if let Some(start) = s.find('(') {
        let end = s.rfind(')').ok_or_else(|| ParseError::Syntax {
            line: lineno,
            msg: "PWL: missing ')'".into(),
        })?;
        &s[start + 1..end]
    } else {
        s.strip_prefix("pwl").unwrap_or("").trim()
    };
    let values: Vec<f64> = inner
        .split_whitespace()
        .map(|tok| parse_value(tok, lineno))
        .collect::<Result<_, _>>()?;
    if values.len() < 2 || !values.len().is_multiple_of(2) {
        return Err(ParseError::Syntax {
            line: lineno,
            msg: format!(
                "PWL requires an even number of values (t v pairs), got {}",
                values.len()
            ),
        });
    }
    Ok(Waveform::Pwl {
        points: values.chunks_exact(2).map(|p| (p[0], p[1])).collect(),
    })
}
