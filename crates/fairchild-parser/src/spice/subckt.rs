use super::common::{canon_node, parse_value};
use super::directives::is_silent_directive;
use super::element::parse_element_expanded;
use crate::expr::Expr;
use crate::{
    Analysis, BehavioralKind, DcSweepSpec, Element, MeasAnalysis, MeasKind, MeasOp, Measurement,
    ModelCard, ParseError, Waveform,
};
use std::collections::{HashMap, HashSet};

// ─── internal types ──────────────────────────────────────────────────────────

/// Internal representation of a `.subckt ... .ends` block collected in pass 1.
pub(super) struct SubcktDef {
    ports: Vec<String>,               // port names (lowercased), in declaration order
    params: Vec<(String, f64)>,       // default parameter values (header + body .param)
    body_lines: Vec<(usize, String)>, // (original_lineno, raw_line) for pass-2 expansion
}

// ─── pass 1: collect definitions ─────────────────────────────────────────────

/// **Pass 1**: split logical lines into subckt definitions, global `.param`
/// values, and the main-body lines that pass 2 will parse.
///
/// Returns `(subckt_defs, global_params, main_lines)`.  Nested `.subckt`
/// definitions and a stray `.ends` are both hard errors.
pub(super) fn collect_defs(
    lines: &[(usize, String)],
) -> Result<
    (
        HashMap<String, SubcktDef>,
        HashMap<String, f64>,
        Vec<(usize, String)>,
    ),
    ParseError,
> {
    let mut subckt_defs: HashMap<String, SubcktDef> = HashMap::new();
    let mut global_params: HashMap<String, f64> = HashMap::new();
    let mut main_lines: Vec<(usize, String)> = Vec::new();

    let mut in_subckt = false;
    let mut current_name = String::new();
    let mut current_def = SubcktDef {
        ports: vec![],
        params: vec![],
        body_lines: vec![],
    };

    for (lineno, line) in lines {
        let lineno = *lineno;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('*') {
            continue;
        }

        let lc = trimmed.to_lowercase();

        if lc == ".end" {
            // End-of-file marker: stop collecting.
            break;
        } else if lc.starts_with(".ends") {
            if !in_subckt {
                return Err(ParseError::Syntax {
                    line: lineno,
                    msg: ".ends without a matching .subckt".into(),
                });
            }
            subckt_defs.insert(std::mem::take(&mut current_name), current_def);
            current_def = SubcktDef {
                ports: vec![],
                params: vec![],
                body_lines: vec![],
            };
            in_subckt = false;
        } else if lc.starts_with(".subckt") {
            if in_subckt {
                return Err(ParseError::Syntax {
                    line: lineno,
                    msg: "nested .subckt definitions are not supported".into(),
                });
            }
            let (name, ports, params) = parse_subckt_header(trimmed, lineno)?;
            current_name = name;
            current_def = SubcktDef {
                ports,
                params,
                body_lines: vec![],
            };
            in_subckt = true;
        } else if lc.starts_with(".param") {
            let pairs = parse_param_directive(trimmed, lineno)?;
            if in_subckt {
                current_def.params.extend(pairs);
            } else {
                global_params.extend(pairs);
            }
        } else if in_subckt {
            current_def.body_lines.push((lineno, line.clone()));
        } else {
            main_lines.push((lineno, line.clone()));
        }
    }

    if in_subckt {
        return Err(ParseError::Syntax {
            line: 0,
            msg: format!(".subckt '{current_name}' has no matching .ends"),
        });
    }

    Ok((subckt_defs, global_params, main_lines))
}

/// Parse `.subckt <name> <port1> ... [param=default ...]`.
pub(super) fn parse_subckt_header(
    line: &str,
    lineno: usize,
) -> Result<(String, Vec<String>, Vec<(String, f64)>), ParseError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 2 {
        return Err(ParseError::FieldCount {
            expected: "≥2 (.subckt name [ports ...])",
            got: tokens.len(),
            line: lineno,
        });
    }
    let name = tokens[1].to_lowercase();
    let mut ports = Vec::new();
    let mut params = Vec::new();
    for tok in &tokens[2..] {
        if let Some((k, v)) = tok.split_once('=') {
            params.push((k.to_lowercase(), parse_value(v, lineno)?));
        } else {
            ports.push(canon_node(tok));
        }
    }
    Ok((name, ports, params))
}

/// Parse `.param name=value [name2=value2 ...]`.
pub(super) fn parse_param_directive(
    line: &str,
    lineno: usize,
) -> Result<Vec<(String, f64)>, ParseError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    let mut pairs = Vec::new();
    for tok in &tokens[1..] {
        if let Some((k, v)) = tok.split_once('=') {
            pairs.push((k.to_lowercase(), parse_value(v, lineno)?));
        }
    }
    Ok(pairs)
}

// ─── expansion helpers ────────────────────────────────────────────────────────

/// Replace `{param_name}` placeholders using `params`.  Errors on undefined names.
pub(super) fn substitute_params(
    line: &str,
    params: &HashMap<String, f64>,
    lineno: usize,
) -> Result<String, ParseError> {
    if !line.contains('{') {
        return Ok(line.to_string());
    }
    let mut result = String::with_capacity(line.len() + 16);
    let mut chars = line.chars().peekable();
    while let Some(ch) = chars.next() {
        if ch != '{' {
            result.push(ch);
            continue;
        }
        let mut name = String::new();
        let mut closed = false;
        for c in chars.by_ref() {
            if c == '}' {
                closed = true;
                break;
            }
            name.push(c);
        }
        if !closed {
            return Err(ParseError::Syntax {
                line: lineno,
                msg: "unclosed '{' in parameter reference".into(),
            });
        }
        let key = name.to_lowercase();
        match params.get(&key) {
            Some(val) => result.push_str(&format!("{val:e}")),
            None => {
                return Err(ParseError::Syntax {
                    line: lineno,
                    msg: format!("undefined parameter '{name}'"),
                })
            }
        }
    }
    Ok(result)
}

/// Map a single node: port names → call-site nets; ground stays "0"; all
/// others get the `{prefix}.` namespace.
pub(super) fn remap_node(node: &str, port_map: &HashMap<String, String>, prefix: &str) -> String {
    if node == "0" {
        return "0".to_string();
    }
    if let Some(mapped) = port_map.get(node) {
        return mapped.clone();
    }
    format!("{prefix}.{node}")
}

/// Remap every node field (and the element name itself) in a flat element.
pub(super) fn remap_element_nodes(
    el: Element,
    port_map: &HashMap<String, String>,
    prefix: &str,
) -> Element {
    let rn = |n: &str| remap_node(n, port_map, prefix);
    match el {
        Element::Resistor {
            name,
            pos,
            neg,
            resistance,
        } => Element::Resistor {
            name: format!("{prefix}.{name}"),
            pos: rn(&pos),
            neg: rn(&neg),
            resistance,
        },
        Element::Capacitor {
            name,
            pos,
            neg,
            capacitance,
        } => Element::Capacitor {
            name: format!("{prefix}.{name}"),
            pos: rn(&pos),
            neg: rn(&neg),
            capacitance,
        },
        Element::Inductor {
            name,
            pos,
            neg,
            inductance,
        } => Element::Inductor {
            name: format!("{prefix}.{name}"),
            pos: rn(&pos),
            neg: rn(&neg),
            inductance,
        },
        Element::VoltageSource {
            name,
            pos,
            neg,
            waveform,
        } => Element::VoltageSource {
            name: format!("{prefix}.{name}"),
            pos: rn(&pos),
            neg: rn(&neg),
            waveform,
        },
        Element::CurrentSource {
            name,
            pos,
            neg,
            waveform,
        } => Element::CurrentSource {
            name: format!("{prefix}.{name}"),
            pos: rn(&pos),
            neg: rn(&neg),
            waveform,
        },
        Element::Diode {
            name,
            anode,
            cathode,
            model_name,
        } => Element::Diode {
            name: format!("{prefix}.{name}"),
            anode: rn(&anode),
            cathode: rn(&cathode),
            model_name,
        },
        Element::Mosfet {
            name,
            drain,
            gate,
            source,
            bulk,
            model_name,
            params,
        } => Element::Mosfet {
            name: format!("{prefix}.{name}"),
            drain: rn(&drain),
            gate: rn(&gate),
            source: rn(&source),
            bulk: rn(&bulk),
            model_name,
            params,
        },
        Element::Bjt {
            name,
            collector,
            base,
            emitter,
            substrate,
            model_name,
            params,
        } => Element::Bjt {
            name: format!("{prefix}.{name}"),
            collector: rn(&collector),
            base: rn(&base),
            emitter: rn(&emitter),
            substrate: rn(&substrate),
            model_name,
            params,
        },
        Element::XOsdi {
            name,
            nets,
            model_name,
            params,
        } => Element::XOsdi {
            name: format!("{prefix}.{name}"),
            nets: nets.iter().map(|n| rn(n)).collect(),
            model_name,
            params,
        },
        Element::Behavioral {
            name,
            pos,
            neg,
            kind,
            expr,
        } => Element::Behavioral {
            name: format!("{prefix}.{name}"),
            pos: rn(&pos),
            neg: rn(&neg),
            kind,
            expr,
        },
        Element::CoupledInductors {
            name,
            l1,
            l2,
            coupling,
        } => Element::CoupledInductors {
            name: format!("{prefix}.{name}"),
            l1: format!("{prefix}.{l1}"),
            l2: format!("{prefix}.{l2}"),
            coupling,
        },
    }
}

/// Expand one `.subckt` instance into a flat `Vec<Element>`.
///
/// `expanding` is the set of subckt names currently on the call stack (cycle
/// detection).  It is mutated in place and restored before returning.
#[allow(clippy::too_many_arguments)]
pub(super) fn expand_instance(
    def_name: &str,
    inst_name: &str,
    call_nets: &[String],
    call_params: &[(String, f64)],
    def: &SubcktDef,
    subckt_defs: &HashMap<String, SubcktDef>,
    global_params: &HashMap<String, f64>,
    expanding: &mut HashSet<String>,
    call_lineno: usize,
) -> Result<Vec<Element>, ParseError> {
    // Port-count check.
    if call_nets.len() != def.ports.len() {
        return Err(ParseError::SubcktPortCount {
            name: def_name.to_string(),
            expected: def.ports.len(),
            got: call_nets.len(),
            line: call_lineno,
        });
    }

    // Cycle detection.
    if expanding.contains(def_name) {
        return Err(ParseError::SubcktCycle {
            name: def_name.to_string(),
        });
    }
    expanding.insert(def_name.to_string());

    // port_map: def port name → call-site net.
    let port_map: HashMap<String, String> = def
        .ports
        .iter()
        .zip(call_nets.iter())
        .map(|(p, n)| (p.clone(), n.clone()))
        .collect();

    // inst_params: global < def defaults < call overrides.
    let mut inst_params: HashMap<String, f64> = global_params.clone();
    for (k, v) in &def.params {
        inst_params.insert(k.clone(), *v);
    }
    for (k, v) in call_params {
        inst_params.insert(k.clone(), *v);
    }

    let mut result: Vec<Element> = Vec::new();

    for (lineno, body_line) in &def.body_lines {
        let lineno = *lineno;
        let trimmed = body_line.trim();
        if trimmed.is_empty() || trimmed.starts_with('*') {
            continue;
        }
        let lc = trimmed.to_lowercase();

        // Directives consumed by collect_defs — skip.
        if lc == ".end"
            || lc.starts_with(".ends")
            || lc.starts_with(".subckt")
            || lc.starts_with(".param")
        {
            continue;
        }
        // .optical inside subckt body is silently ignored: optical net
        // membership is determined by the caller's top-level declarations.
        if lc.starts_with(".optical") {
            continue;
        }
        if lc.starts_with('.') {
            if !is_silent_directive(&lc) {
                let directive = lc.split_whitespace().next().unwrap_or(&lc).to_string();
                return Err(ParseError::UnsupportedDirective {
                    directive,
                    line: lineno,
                });
            }
            continue;
        }

        let substituted = substitute_params(trimmed, &inst_params, lineno)?;
        for el in parse_element_expanded(&substituted, lineno)? {
            let el = remap_element_nodes(el, &port_map, inst_name);

            // Recurse if this element is a nested subckt instance.
            let is_subckt_inst = if let Element::XOsdi { ref model_name, .. } = el {
                subckt_defs.contains_key(model_name)
            } else {
                false
            };

            if is_subckt_inst {
                if let Element::XOsdi {
                    ref name,
                    ref nets,
                    ref model_name,
                    ref params,
                } = el
                {
                    let nested_def = subckt_defs.get(model_name).unwrap();
                    let nested = expand_instance(
                        model_name,
                        name,
                        nets,
                        params,
                        nested_def,
                        subckt_defs,
                        &inst_params,
                        expanding,
                        lineno,
                    )?;
                    result.extend(nested);
                }
            } else {
                result.push(el);
            }
        }
    }

    expanding.remove(def_name);
    Ok(result)
}
