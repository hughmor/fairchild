use crate::{Analysis, Element, ModelCard, Netlist, ParseError, Waveform};

/// Parse a SPICE netlist from a string.
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

        if lc == ".end" || lc.starts_with(".ends") {
            break;
        } else if lc.starts_with(".op") {
            netlist.analyses.push(Analysis::Op);
        } else if lc.starts_with(".tran") {
            netlist.analyses.push(parse_tran(&lc, *raw_lineno)?);
        } else if lc.starts_with(".model") {
            if let Some(card) = parse_model(&lc, *raw_lineno)? {
                netlist.models.push(card);
            }
        } else if lc.starts_with(".osdi") {
            let tokens: Vec<&str> = trimmed.split_whitespace().collect();
            if tokens.len() >= 2 {
                netlist.osdi_paths.push(tokens[1].to_string());
            }
        } else if lc.starts_with('.') {
            // Ignore other directives (.param, .ic, .meas, etc.) for now.
        } else {
            let el = parse_element(trimmed, *raw_lineno)?;
            netlist.elements.push(el);
        }
    }

    Ok(netlist)
}

/// Parse `.tran <step> <stop>` directive.
fn parse_tran(line: &str, lineno: usize) -> Result<Analysis, ParseError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 3 {
        return Err(ParseError::FieldCount { expected: "≥3 (.tran step stop)", got: tokens.len(), line: lineno });
    }
    Ok(Analysis::Tran {
        step: parse_value(tokens[1], lineno)?,
        stop: parse_value(tokens[2], lineno)?,
    })
}

/// Canonicalise a node name: lowercase; "gnd" → "0".
fn canon_node(s: &str) -> String {
    let s = s.to_lowercase();
    if s == "gnd" { "0".to_string() } else { s }
}

/// Parse an SPICE suffix (k, meg, m, u, n, p, f, g, t) into a float.
fn parse_value(s: &str, lineno: usize) -> Result<f64, ParseError> {
    let s_lc = s.to_lowercase();
    // Strip trailing alphabetic suffix that might follow a number (e.g. "1.0v", "5.0a")
    // We support SPICE multiplier suffixes only.
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
        'c' => {
            if tokens.len() < 4 {
                return Err(ParseError::FieldCount { expected: "≥4", got: tokens.len(), line: lineno });
            }
            Ok(Element::Capacitor {
                name,
                pos: canon_node(tokens[1]),
                neg: canon_node(tokens[2]),
                capacitance: parse_value(tokens[3], lineno)?,
            })
        }
        'l' => {
            if tokens.len() < 4 {
                return Err(ParseError::FieldCount { expected: "≥4", got: tokens.len(), line: lineno });
            }
            Ok(Element::Inductor {
                name,
                pos: canon_node(tokens[1]),
                neg: canon_node(tokens[2]),
                inductance: parse_value(tokens[3], lineno)?,
            })
        }
        'v' => {
            let waveform = parse_waveform(&tokens, lineno)?;
            Ok(Element::VoltageSource {
                name,
                pos: canon_node(tokens[1]),
                neg: canon_node(tokens[2]),
                waveform,
            })
        }
        'i' => {
            let waveform = parse_waveform(&tokens, lineno)?;
            Ok(Element::CurrentSource {
                name,
                pos: canon_node(tokens[1]),
                neg: canon_node(tokens[2]),
                waveform,
            })
        }
        'd' => {
            if tokens.len() < 4 {
                return Err(ParseError::FieldCount { expected: "≥4", got: tokens.len(), line: lineno });
            }
            Ok(Element::Diode {
                name,
                anode: canon_node(tokens[1]),
                cathode: canon_node(tokens[2]),
                model_name: tokens[3].to_lowercase(),
            })
        }
        _ => Err(ParseError::UnknownElement { letter, line: lineno }),
    }
}

/// Parse `.model <name> <kind> [params]` — kind is the first letter of device type.
/// Params may be bare `Name=Value` tokens or wrapped in parentheses.
fn parse_model(line: &str, lineno: usize) -> Result<Option<ModelCard>, ParseError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 3 {
        return Ok(None);  // malformed, ignore silently
    }
    let name = tokens[1].to_string();
    let kind = tokens[2].to_lowercase();

    // Join remainder, strip parentheses, then split on whitespace.
    let rest = tokens[3..].join(" ");
    let rest = rest.trim_matches(|c| c == '(' || c == ')').trim().to_string();
    // Handle inner parentheses too: remove all parens.
    let rest = rest.replace('(', " ").replace(')', " ");

    let mut params = Vec::new();
    for tok in rest.split_whitespace() {
        if let Some((k, v)) = tok.split_once('=') {
            if let Ok(val) = parse_value(v, lineno) {
                params.push((k.to_lowercase(), val));
            }
        }
    }

    Ok(Some(ModelCard { name, kind, params }))
}

/// Parse the waveform specification from a V/I source token list.
/// Handles:
///   `V1 a b 5`               → Dc(5)
///   `V1 a b DC 5`            → Dc(5)
///   `V1 a b PULSE(0 1 0 ...)`→ Pulse{...}
fn parse_waveform(tokens: &[&str], lineno: usize) -> Result<Waveform, ParseError> {
    if tokens.len() < 4 {
        return Err(ParseError::FieldCount { expected: "≥4", got: tokens.len(), line: lineno });
    }

    // Rejoin tokens[3..] so PULSE( 0 1 ...) with spaces around parens also works.
    let rest = tokens[3..].join(" ");
    let rest_lc = rest.to_lowercase();

    if rest_lc.starts_with("pulse") {
        return parse_pulse(&rest_lc, lineno);
    }

    // DC value: either "DC value" or just "value"
    let tok = tokens[3].to_lowercase();
    if tok == "dc" {
        if tokens.len() < 5 {
            return Err(ParseError::FieldCount { expected: "≥5 (DC keyword)", got: tokens.len(), line: lineno });
        }
        Ok(Waveform::Dc(parse_value(tokens[4], lineno)?))
    } else {
        Ok(Waveform::Dc(parse_value(tokens[3], lineno)?))
    }
}

/// Parse `PULSE(v0 v1 td tr tf pw per)` into a `Waveform::Pulse`.
/// Accepts both compact `PULSE(0 1 ...)` and spaced forms.
fn parse_pulse(s: &str, lineno: usize) -> Result<Waveform, ParseError> {
    // Extract the parenthesised argument list.
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

    // v0 v1 td tr tf pw per  (7 fields, defaults for omitted trailing ones)
    let get = |i: usize, default: f64| -> Result<f64, ParseError> {
        parts.get(i).map_or(Ok(default), |s| parse_value(s, lineno))
    };

    if parts.len() < 2 {
        return Err(ParseError::FieldCount { expected: "≥2 (PULSE v0 v1 ...)", got: parts.len(), line: lineno });
    }

    Ok(Waveform::Pulse {
        v0:  get(0, 0.0)?,
        v1:  get(1, 0.0)?,
        td:  get(2, 0.0)?,
        tr:  get(3, 0.0)?,
        tf:  get(4, 0.0)?,
        pw:  get(5, f64::INFINITY)?,
        per: get(6, f64::INFINITY)?,
    })
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
    fn parse_rc_tran() {
        let input = "* RC\nV1 in 0 PULSE(0 1 0 1n 1n 10m 20m)\nR1 in out 1k\nC1 out 0 1u\n.tran 1u 5m\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(netlist.elements.len(), 3);
        match &netlist.analyses[0] {
            Analysis::Tran { step, stop } => {
                assert!((step - 1e-6).abs() < 1e-12);
                assert!((stop - 5e-3).abs() < 1e-12);
            }
            _ => panic!("expected Tran analysis"),
        }
    }

    #[test]
    fn parse_pulse_waveform() {
        let input = "* Pulse\nV1 a 0 PULSE(0 1 0 1n 1n 10m 20m)\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        if let Element::VoltageSource { waveform: Waveform::Pulse { v0, v1, tr, .. }, .. } = &netlist.elements[0] {
            assert!((v0 - 0.0).abs() < 1e-12);
            assert!((v1 - 1.0).abs() < 1e-12);
            assert!((tr - 1e-9).abs() < 1e-15);
        } else {
            panic!("expected PULSE VoltageSource");
        }
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

    #[test]
    fn parse_diode_element() {
        let input = "* Diode\nD1 anode cathode myd\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(netlist.elements.len(), 1);
        if let Element::Diode { name, anode, cathode, model_name } = &netlist.elements[0] {
            assert_eq!(name, "d1");
            assert_eq!(anode, "anode");
            assert_eq!(cathode, "cathode");
            assert_eq!(model_name, "myd");
        } else {
            panic!("expected Diode element");
        }
    }

    #[test]
    fn parse_model_card() {
        let input = "* test\n.model myd D (Is=1e-14 N=1)\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(netlist.models.len(), 1);
        let m = &netlist.models[0];
        assert_eq!(m.name, "myd");
        assert_eq!(m.kind, "d");
        let is = m.params.iter().find(|(k, _)| k == "is").map(|(_, v)| *v).unwrap();
        assert!((is - 1e-14).abs() < 1e-20, "is={is}");
        let n = m.params.iter().find(|(k, _)| k == "n").map(|(_, v)| *v).unwrap();
        assert!((n - 1.0).abs() < 1e-12, "n={n}");
    }

    #[test]
    fn parse_model_card_no_parens() {
        let input = "* test\n.model myd D Is=2.52e-9 N=1.752\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        let m = &netlist.models[0];
        assert_eq!(m.kind, "d");
        assert_eq!(m.params.len(), 2);
    }

    #[test]
    fn pulse_waveform_at() {
        let w = Waveform::Pulse { v0: 0.0, v1: 1.0, td: 0.0, tr: 1e-9, tf: 1e-9, pw: 1.0, per: 2.0 };
        assert!((w.at(0.0) - 0.0).abs() < 1e-12);
        assert!((w.at(1e-9) - 1.0).abs() < 1e-6); // fully risen
        assert!((w.at(0.5e-9) - 0.5).abs() < 1e-6); // mid-rise
    }
}
