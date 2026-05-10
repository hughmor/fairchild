use crate::{Analysis, Element, Netlist, ParseError};

/// Parse a SPICE netlist from a string.
///
/// SPICE format rules implemented here:
/// - Line 1 is always the title (even if it looks like an element).
/// - Lines starting with `*` are comments.
/// - Lines starting with `+` are continuations of the previous line.
/// - Case-insensitive keywords and node names (except node "0"/"gnd").
/// - Element type determined by first character of the element name.
pub fn parse_spice(input: &str) -> Result<Netlist, ParseError> {
    let logical_lines = logical_lines(input);
    let mut netlist = Netlist::default();

    for (lineno, (raw_lineno, line)) in logical_lines.iter().enumerate() {
        if lineno == 0 {
            netlist.title = line.trim().to_string();
            continue;
        }

        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('*') {
            continue;
        }

        let lc = trimmed.to_lowercase();

        if lc.starts_with(".end") && (lc == ".end" || lc.starts_with(".ends")) {
            break;
        } else if lc.starts_with(".op") {
            netlist.analyses.push(Analysis::Op);
        } else if lc.starts_with('.') {
            // Ignore other directives for now (.tran, .param, etc.)
        } else {
            let el = parse_element(trimmed, *raw_lineno)?;
            netlist.elements.push(el);
        }
    }

    Ok(netlist)
}

/// Canonicalise a node name: lowercase; "gnd" → "0".
fn canon_node(s: &str) -> String {
    let s = s.to_lowercase();
    if s == "gnd" { "0".to_string() } else { s }
}

/// Parse an SPICE suffix (k=1e3, meg=1e6, etc.) into a float.
fn parse_value(s: &str, lineno: usize) -> Result<f64, ParseError> {
    let s_lc = s.to_lowercase();
    let (num_part, multiplier) = if let Some(n) = s_lc.strip_suffix("meg") {
        (n, 1e6)
    } else if let Some(n) = s_lc.strip_suffix('k') {
        (n, 1e3)
    } else if let Some(n) = s_lc.strip_suffix('m') {
        (n, 1e-3)
    } else if let Some(n) = s_lc.strip_suffix('u') {
        (n, 1e-6)
    } else if let Some(n) = s_lc.strip_suffix('n') {
        (n, 1e-9)
    } else if let Some(n) = s_lc.strip_suffix('p') {
        (n, 1e-12)
    } else if let Some(n) = s_lc.strip_suffix('f') {
        (n, 1e-15)
    } else if let Some(n) = s_lc.strip_suffix('g') {
        (n, 1e9)
    } else if let Some(n) = s_lc.strip_suffix('t') {
        (n, 1e12)
    } else {
        (s_lc.as_str(), 1.0)
    };

    num_part.parse::<f64>().map(|v| v * multiplier).map_err(|e| ParseError::BadNumber {
        value: s.to_string(),
        line: lineno,
        source: e,
    })
}

fn parse_element(line: &str, lineno: usize) -> Result<Element, ParseError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    let name = tokens[0].to_lowercase();
    let letter = name.chars().next().unwrap();

    match letter {
        'r' => {
            if tokens.len() < 4 {
                return Err(ParseError::FieldCount { expected: "≥4", got: tokens.len(), line: lineno });
            }
            Ok(Element::Resistor {
                name,
                pos: canon_node(tokens[1]),
                neg: canon_node(tokens[2]),
                resistance: parse_value(tokens[3], lineno)?,
            })
        }
        'v' => {
            // V<name> <pos> <neg> [DC] <value>
            let dc = parse_vsrc_value(&tokens, lineno)?;
            Ok(Element::VoltageSource {
                name,
                pos: canon_node(tokens[1]),
                neg: canon_node(tokens[2]),
                dc,
            })
        }
        'i' => {
            let dc = parse_vsrc_value(&tokens, lineno)?;
            Ok(Element::CurrentSource {
                name,
                pos: canon_node(tokens[1]),
                neg: canon_node(tokens[2]),
                dc,
            })
        }
        _ => Err(ParseError::UnknownElement { letter, line: lineno }),
    }
}

/// Extract the DC value from a V/I source token list.
/// Handles: `V1 a b 5`, `V1 a b DC 5`, `V1 a b dc 5.0`
fn parse_vsrc_value(tokens: &[&str], lineno: usize) -> Result<f64, ParseError> {
    if tokens.len() < 4 {
        return Err(ParseError::FieldCount { expected: "≥4", got: tokens.len(), line: lineno });
    }
    let tok = tokens[3].to_lowercase();
    if tok == "dc" {
        if tokens.len() < 5 {
            return Err(ParseError::FieldCount { expected: "≥5 (DC keyword present)", got: tokens.len(), line: lineno });
        }
        parse_value(tokens[4], lineno)
    } else {
        parse_value(tokens[3], lineno)
    }
}

/// Join continuation lines (starting with `+`) and return (original_lineno, joined_line) pairs.
fn logical_lines(input: &str) -> Vec<(usize, String)> {
    let mut result: Vec<(usize, String)> = Vec::new();
    for (i, raw) in input.lines().enumerate() {
        let lineno = i + 1;
        let trimmed = raw.trim_start();
        if trimmed.starts_with('+') {
            if let Some(last) = result.last_mut() {
                last.1.push(' ');
                last.1.push_str(trimmed[1..].trim());
            }
        } else {
            result.push((lineno, raw.to_string()));
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_voltage_divider() {
        let input = "* Voltage divider\nV1 in 0 DC 1.0\nR1 in mid 1k\nR2 mid 0 1k\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(netlist.elements.len(), 3);
        assert_eq!(netlist.analyses.len(), 1);
    }

    #[test]
    fn parse_suffix_k() {
        let v = parse_value("2k", 1).unwrap();
        assert!((v - 2000.0).abs() < 1e-9);
    }

    #[test]
    fn parse_suffix_meg() {
        let v = parse_value("1meg", 1).unwrap();
        assert!((v - 1e6).abs() < 1.0);
    }

    #[test]
    fn gnd_canonical() {
        assert_eq!(canon_node("GND"), "0");
        assert_eq!(canon_node("gnd"), "0");
        assert_eq!(canon_node("0"), "0");
    }
}
