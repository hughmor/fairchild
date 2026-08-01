use super::common::{canon_node, expand_bus_vectors, parse_value};
use super::waveforms::parse_waveform;
use crate::expr::Expr;
use crate::{BehavioralKind, Element, ModelCard, ParseError};

pub(super) fn parse_element(line: &str, lineno: usize) -> Result<Element, ParseError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    let name = tokens[0].to_lowercase();
    let letter = name.chars().next().unwrap();

    match letter {
        'r' => {
            if tokens.len() < 4 {
                return Err(ParseError::FieldCount {
                    expected: "≥4",
                    got: tokens.len(),
                    line: lineno,
                });
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
                return Err(ParseError::FieldCount {
                    expected: "≥4",
                    got: tokens.len(),
                    line: lineno,
                });
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
                return Err(ParseError::FieldCount {
                    expected: "≥4",
                    got: tokens.len(),
                    line: lineno,
                });
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
                return Err(ParseError::FieldCount {
                    expected: "≥4",
                    got: tokens.len(),
                    line: lineno,
                });
            }
            let mut params = Vec::new();
            for tok in &tokens[4..] {
                if let Some((k, v)) = tok.split_once('=') {
                    if let Ok(val) = parse_value(v, lineno) {
                        params.push((k.to_lowercase(), val));
                    }
                }
            }
            Ok(Element::Diode {
                name,
                anode: canon_node(tokens[1]),
                cathode: canon_node(tokens[2]),
                model_name: tokens[3].to_lowercase(),
                params,
            })
        }
        'b' => {
            // B-element behavioural source: `Bname n+ n- V=<expr>` or `I=<expr>`.
            // The expression may contain spaces, so we re-stitch tokens[3..].
            if tokens.len() < 4 {
                return Err(ParseError::FieldCount {
                    expected: "≥4 (Bname n+ n- V=expr|I=expr)",
                    got: tokens.len(),
                    line: lineno,
                });
            }
            let pos = canon_node(tokens[1]);
            let neg = canon_node(tokens[2]);
            let rest = tokens[3..].join(" ");
            let lc = rest.to_lowercase();
            // Recognise leading `v=`, `v =`, `i=`, `i =`.
            let (kind, expr_str) = if let Some(stripped) = lc.strip_prefix("v=") {
                (
                    BehavioralKind::Voltage,
                    rest[rest.len() - stripped.len()..].to_string(),
                )
            } else if let Some(stripped) = lc.strip_prefix("v =") {
                (
                    BehavioralKind::Voltage,
                    rest[rest.len() - stripped.len()..].to_string(),
                )
            } else if let Some(stripped) = lc.strip_prefix("i=") {
                (
                    BehavioralKind::Current,
                    rest[rest.len() - stripped.len()..].to_string(),
                )
            } else if let Some(stripped) = lc.strip_prefix("i =") {
                (
                    BehavioralKind::Current,
                    rest[rest.len() - stripped.len()..].to_string(),
                )
            } else {
                return Err(ParseError::Syntax {
                    line: lineno,
                    msg: "B-element requires V=<expr> or I=<expr>".into(),
                });
            };
            // Strip optional `{ … }` wrapping.
            let expr_str = expr_str.trim();
            let expr_str = expr_str
                .trim_start_matches('{')
                .trim_end_matches('}')
                .trim();
            let expr = Expr::parse(expr_str).map_err(|e| ParseError::Syntax {
                line: lineno,
                msg: format!("B-element expression: {e}"),
            })?;
            Ok(Element::Behavioral {
                name,
                pos,
                neg,
                kind,
                expr,
            })
        }
        'm' => {
            if tokens.len() < 6 {
                return Err(ParseError::FieldCount {
                    expected: "≥6 (Mname drain gate source bulk model)",
                    got: tokens.len(),
                    line: lineno,
                });
            }
            let mut params = Vec::new();
            for tok in &tokens[6..] {
                if let Some((k, v)) = tok.split_once('=') {
                    if let Ok(val) = parse_value(v, lineno) {
                        params.push((k.to_lowercase(), val));
                    }
                }
            }
            Ok(Element::Mosfet {
                name,
                drain: canon_node(tokens[1]),
                gate: canon_node(tokens[2]),
                source: canon_node(tokens[3]),
                bulk: canon_node(tokens[4]),
                model_name: tokens[5].to_lowercase(),
                params,
            })
        }
        'q' => {
            // Q<name> collector base emitter [substrate] model [param=val ...]
            // Substrate is optional; defaults to "0" (ground) when omitted.
            if tokens.len() < 5 {
                return Err(ParseError::FieldCount {
                    expected: "≥5 (Qname collector base emitter model)",
                    got: tokens.len(),
                    line: lineno,
                });
            }
            // Detect whether tokens[4] is the model name or the substrate node.
            // Strategy: if tokens.len() >= 6 and tokens[4] does NOT contain '='
            // and tokens[5] does NOT look like a value, treat tokens[4] as
            // substrate. Otherwise, treat tokens[4] as the model name.
            let (substrate, model_name, param_start) =
                if tokens.len() >= 6 && !tokens[4].contains('=') && !tokens[5].contains('=') {
                    // Five positional fields: C B E S model
                    (canon_node(tokens[4]), tokens[5].to_lowercase(), 6)
                } else {
                    // Four positional fields: C B E model
                    ("0".to_string(), tokens[4].to_lowercase(), 5)
                };
            let mut params = Vec::new();
            for tok in &tokens[param_start..] {
                if let Some((k, v)) = tok.split_once('=') {
                    if let Ok(val) = parse_value(v, lineno) {
                        params.push((k.to_lowercase(), val));
                    }
                }
            }
            Ok(Element::Bjt {
                name,
                collector: canon_node(tokens[1]),
                base: canon_node(tokens[2]),
                emitter: canon_node(tokens[3]),
                substrate,
                model_name,
                params,
            })
        }
        'k' => {
            // K<name> L1 L2 coupling
            if tokens.len() < 4 {
                return Err(ParseError::FieldCount {
                    expected: "≥4 (Kname L1 L2 coupling)",
                    got: tokens.len(),
                    line: lineno,
                });
            }
            let coupling: f64 = parse_value(tokens[3], lineno)?;
            if !(0.0..=1.0).contains(&coupling) {
                return Err(ParseError::Syntax {
                    line: lineno,
                    msg: format!("coupling must be in [0,1], got {coupling}"),
                });
            }
            Ok(Element::CoupledInductors {
                name,
                l1: tokens[1].to_lowercase(),
                l2: tokens[2].to_lowercase(),
                coupling,
            })
        }
        's' | 'w' => {
            // S<name> N+ N- NC+ NC- MODEL [ON|OFF]
            // W<name> N+ N- VSOURCE  MODEL [ON|OFF]
            let is_current = letter == 'w';
            let (min_tokens, expected) = if is_current {
                (5, "≥5 (Wname n+ n- vsource model [ON|OFF])")
            } else {
                (6, "≥6 (Sname n+ n- nc+ nc- model [ON|OFF])")
            };
            if tokens.len() < min_tokens {
                return Err(ParseError::FieldCount {
                    expected,
                    got: tokens.len(),
                    line: lineno,
                });
            }
            // The trailing ON/OFF keyword is optional and comes after the
            // model name; anything else there is a typo worth naming.
            // W puts the controlling source at 3 and the model at 4; S has
            // four nets, so its model is at 5.
            let model_idx = if is_current { 4 } else { 5 };
            let initial_on = match tokens.get(model_idx + 1) {
                None => false,
                Some(t) if t.eq_ignore_ascii_case("on") => true,
                Some(t) if t.eq_ignore_ascii_case("off") => false,
                Some(t) => {
                    return Err(ParseError::Syntax {
                        line: lineno,
                        msg: format!("switch {name}: expected ON, OFF or nothing after the model name, got '{t}'"),
                    })
                }
            };
            let model_name = tokens[model_idx].to_lowercase();
            if is_current {
                Ok(Element::CurrentSwitch {
                    name,
                    pos: canon_node(tokens[1]),
                    neg: canon_node(tokens[2]),
                    ctrl_vsrc: tokens[3].to_lowercase(),
                    model_name,
                    initial_on,
                })
            } else {
                Ok(Element::VoltageSwitch {
                    name,
                    pos: canon_node(tokens[1]),
                    neg: canon_node(tokens[2]),
                    ctrl_pos: canon_node(tokens[3]),
                    ctrl_neg: canon_node(tokens[4]),
                    model_name,
                    initial_on,
                })
            }
        }
        't' => {
            // T<name> A+ A- B+ B- Z0=<Ω> (TD=<s> | F=<Hz> [NL=<wavelengths>])
            if tokens.len() < 6 {
                return Err(ParseError::FieldCount {
                    expected: "≥6 (Tname A+ A- B+ B- Z0=.. TD=..)",
                    got: tokens.len(),
                    line: lineno,
                });
            }
            let mut z0: Option<f64> = None;
            let mut td: Option<f64> = None;
            let mut freq: Option<f64> = None;
            let mut nl: f64 = 0.25; // default quarter-wave per ngspice
            for tok in &tokens[5..] {
                if let Some((k, v)) = tok.split_once('=') {
                    let val = parse_value(v, lineno)?;
                    match k.to_lowercase().as_str() {
                        "z0" | "zo" => z0 = Some(val),
                        "td" => td = Some(val),
                        "f" => freq = Some(val),
                        "nl" => nl = val,
                        _ => {}
                    }
                }
            }
            let z0 = z0.ok_or_else(|| ParseError::Syntax {
                line: lineno,
                msg: "transmission line requires Z0=<ohms>".to_string(),
            })?;
            // Delay: explicit TD wins; else derive from F (and NL): TD = NL / F.
            let td = match (td, freq) {
                (Some(td), _) => td,
                (None, Some(f)) if f > 0.0 => nl / f,
                _ => {
                    return Err(ParseError::Syntax {
                        line: lineno,
                        msg: "transmission line requires TD=<s> or F=<Hz>".to_string(),
                    })
                }
            };
            Ok(Element::TransmissionLine {
                name,
                a_pos: canon_node(tokens[1]),
                a_neg: canon_node(tokens[2]),
                b_pos: canon_node(tokens[3]),
                b_neg: canon_node(tokens[4]),
                z0,
                td,
            })
        }
        'x' => {
            if tokens.len() < 3 {
                return Err(ParseError::FieldCount {
                    expected: "≥3 (Xname net0 ... model_name)",
                    got: tokens.len(),
                    line: lineno,
                });
            }
            let mut positional: Vec<&str> = Vec::new();
            let mut params: Vec<(String, f64)> = Vec::new();
            for tok in &tokens[1..] {
                if tok.contains('=') {
                    if let Some((k, v)) = tok.split_once('=') {
                        if let Ok(val) = parse_value(v, lineno) {
                            params.push((k.to_lowercase(), val));
                        }
                    }
                } else {
                    positional.push(tok);
                }
            }
            if positional.len() < 2 {
                return Err(ParseError::FieldCount {
                    expected: "≥2 positional (at least one net + model_name)",
                    got: positional.len(),
                    line: lineno,
                });
            }
            let model_name = positional.last().unwrap().to_lowercase();
            let nets: Vec<String> = positional[..positional.len() - 1]
                .iter()
                .flat_map(|s| expand_bus_vectors(s))
                .map(|s| canon_node(&s))
                .collect();
            Ok(Element::XOsdi {
                name,
                nets,
                model_name,
                params,
            })
        }
        _ => Err(ParseError::UnknownElement {
            letter,
            line: lineno,
        }),
    }
}

/// Like `parse_element` but recognises parasitic key=val tokens on R, L, C lines
/// and expands them into equivalent sub-networks with internal `__`-prefixed nodes.
///
/// Supported parasitics:
/// - R: `cpar=<val>` — parallel capacitance
/// - L: `rser=<val>` — series ESR, `cpar=<val>` — parallel winding capacitance
/// - C: `esr=<val>` — series resistance, `esl=<val>` — series inductance,
///   `rpar=<val>` — parallel leakage resistance
pub(super) fn parse_element_expanded(
    line: &str,
    lineno: usize,
) -> Result<Vec<Element>, ParseError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 4 {
        return parse_element(line, lineno).map(|e| vec![e]);
    }

    let name = tokens[0].to_lowercase();
    let letter = name.chars().next().unwrap();

    match letter {
        'r' | 'l' | 'c' => {}
        _ => return parse_element(line, lineno).map(|e| vec![e]),
    }

    let mut rser: Option<f64> = None;
    let mut cpar: Option<f64> = None;
    let mut esr: Option<f64> = None;
    let mut esl: Option<f64> = None;
    let mut rpar: Option<f64> = None;
    for tok in &tokens[4..] {
        if let Some((k, v)) = tok.split_once('=') {
            if let Ok(val) = parse_value(v, lineno) {
                match k.to_lowercase().as_str() {
                    "rser" => rser = Some(val),
                    "cpar" => cpar = Some(val),
                    "esr" => esr = Some(val),
                    "esl" => esl = Some(val),
                    "rpar" => rpar = Some(val),
                    _ => {}
                }
            }
        }
    }

    if rser.is_none() && cpar.is_none() && esr.is_none() && esl.is_none() && rpar.is_none() {
        return parse_element(line, lineno).map(|e| vec![e]);
    }

    let pos: String = canon_node(tokens[1]);
    let neg: String = canon_node(tokens[2]);
    let val: f64 = parse_value(tokens[3], lineno)?;

    let mut elements: Vec<Element> = Vec::new();

    match letter {
        'r' => {
            elements.push(Element::Resistor {
                name: name.clone(),
                pos: pos.clone(),
                neg: neg.clone(),
                resistance: val,
            });
            if let Some(cv) = cpar {
                elements.push(Element::Capacitor {
                    name: format!("__c_{name}_cpar"),
                    pos: pos.clone(),
                    neg: neg.clone(),
                    capacitance: cv,
                });
            }
        }
        'l' => {
            // Series ESR: L → internal node → R → neg.
            let l_neg = if rser.is_some() {
                format!("__{name}_rn")
            } else {
                neg.clone()
            };
            elements.push(Element::Inductor {
                name: name.clone(),
                pos: pos.clone(),
                neg: l_neg.clone(),
                inductance: val,
            });
            if let Some(rv) = rser {
                elements.push(Element::Resistor {
                    name: format!("__r_{name}_rser"),
                    pos: l_neg,
                    neg: neg.clone(),
                    resistance: rv,
                });
            }
            if let Some(cv) = cpar {
                elements.push(Element::Capacitor {
                    name: format!("__c_{name}_cpar"),
                    pos: pos.clone(),
                    neg: neg.clone(),
                    capacitance: cv,
                });
            }
        }
        'c' => {
            // Series chain: pos --[ESL?]--[ESR?]-- c_pos --[C]-- neg
            let mut c_pos = pos.clone();
            if let Some(lv) = esl {
                let int_node = format!("__{name}_esln");
                elements.push(Element::Inductor {
                    name: format!("__l_{name}_esl"),
                    pos: c_pos,
                    neg: int_node.clone(),
                    inductance: lv,
                });
                c_pos = int_node;
            }
            if let Some(rv) = esr {
                let int_node = format!("__{name}_esrn");
                elements.push(Element::Resistor {
                    name: format!("__r_{name}_esr"),
                    pos: c_pos,
                    neg: int_node.clone(),
                    resistance: rv,
                });
                c_pos = int_node;
            }
            elements.push(Element::Capacitor {
                name: name.clone(),
                pos: c_pos,
                neg: neg.clone(),
                capacitance: val,
            });
            if let Some(rv) = rpar {
                elements.push(Element::Resistor {
                    name: format!("__r_{name}_rpar"),
                    pos: pos.clone(),
                    neg: neg.clone(),
                    resistance: rv,
                });
            }
        }
        _ => unreachable!(),
    }

    Ok(elements)
}

pub(super) fn parse_model(line: &str, lineno: usize) -> Result<Option<ModelCard>, ParseError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 3 {
        return Ok(None);
    }
    let name = tokens[1].to_string();
    let kind = tokens[2].to_lowercase();
    let rest = tokens[3..].join(" ");

    let mut params = Vec::new();
    let mut expr_params = Vec::new();
    for tok in split_model_assignments(&rest) {
        let Some((k, v)) = tok.split_once('=') else {
            continue;
        };
        let key = k.to_lowercase();
        // A quoted value (`"…"` or `{…}`) is always an expression; an unquoted
        // value is numeric if it parses as one, otherwise a bare expression.
        if let Some(inner) = strip_expr_quotes(v) {
            expr_params.push((key, inner.to_string()));
        } else if let Ok(val) = parse_value(v, lineno) {
            params.push((key, val));
        } else {
            expr_params.push((key, v.to_string()));
        }
    }
    Ok(Some(ModelCard {
        name,
        kind,
        params,
        expr_params,
    }))
}

/// Split a `.model` parameter list into `key=value` tokens, keeping quoted
/// expression spans (`"…"` / `{…}`) intact and treating SPICE list parentheses
/// (outside quotes) as separators. Whitespace inside a quote is preserved.
fn split_model_assignments(rest: &str) -> Vec<String> {
    let mut toks = Vec::new();
    let mut cur = String::new();
    let mut quote: Option<char> = None; // the opening delimiter
    for c in rest.chars() {
        match quote {
            Some(open) => {
                cur.push(c);
                if (open == '"' && c == '"') || (open == '{' && c == '}') {
                    quote = None;
                }
            }
            None => match c {
                '"' | '{' => {
                    quote = Some(c);
                    cur.push(c);
                }
                '(' | ')' | ' ' | '\t' => {
                    if !cur.is_empty() {
                        toks.push(std::mem::take(&mut cur));
                    }
                }
                _ => cur.push(c),
            },
        }
    }
    if !cur.is_empty() {
        toks.push(cur);
    }
    toks
}

/// If `v` is a quoted/braced expression value, return its inner text; else None.
fn strip_expr_quotes(v: &str) -> Option<&str> {
    let b = v.as_bytes();
    if b.len() >= 2
        && ((b[0] == b'"' && b[b.len() - 1] == b'"') || (b[0] == b'{' && b[b.len() - 1] == b'}'))
    {
        Some(&v[1..v.len() - 1])
    } else {
        None
    }
}
