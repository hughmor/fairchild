use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use crate::{AcVariation, Analysis, DcSweepSpec, Element, ModelCard, Netlist, ParseError, Waveform};

// ─── internal types ──────────────────────────────────────────────────────────

/// Internal representation of a `.subckt ... .ends` block collected in pass 1.
struct SubcktDef {
    ports:      Vec<String>,           // port names (lowercased), in declaration order
    params:     Vec<(String, f64)>,    // default parameter values (header + body .param)
    body_lines: Vec<(usize, String)>,  // (original_lineno, raw_line) for pass-2 expansion
}

// ─── pass 1: collect definitions ─────────────────────────────────────────────

/// **Pass 1**: split logical lines into subckt definitions, global `.param`
/// values, and the main-body lines that pass 2 will parse.
///
/// Returns `(subckt_defs, global_params, main_lines)`.  Nested `.subckt`
/// definitions and a stray `.ends` are both hard errors.
fn collect_defs(
    lines: &[(usize, String)],
) -> Result<(HashMap<String, SubcktDef>, HashMap<String, f64>, Vec<(usize, String)>), ParseError> {
    let mut subckt_defs:  HashMap<String, SubcktDef> = HashMap::new();
    let mut global_params: HashMap<String, f64>       = HashMap::new();
    let mut main_lines:   Vec<(usize, String)>        = Vec::new();

    let mut in_subckt    = false;
    let mut current_name = String::new();
    let mut current_def  = SubcktDef { ports: vec![], params: vec![], body_lines: vec![] };

    for (lineno, line) in lines {
        let lineno  = *lineno;
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
            current_def = SubcktDef { ports: vec![], params: vec![], body_lines: vec![] };
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
            current_def  = SubcktDef { ports, params, body_lines: vec![] };
            in_subckt    = true;
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
fn parse_subckt_header(
    line:   &str,
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
    let mut ports  = Vec::new();
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
fn parse_param_directive(line: &str, lineno: usize) -> Result<Vec<(String, f64)>, ParseError> {
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
fn substitute_params(
    line:   &str,
    params: &HashMap<String, f64>,
    lineno: usize,
) -> Result<String, ParseError> {
    if !line.contains('{') {
        return Ok(line.to_string());
    }
    let mut result = String::with_capacity(line.len() + 16);
    let mut chars  = line.chars().peekable();
    while let Some(ch) = chars.next() {
        if ch != '{' {
            result.push(ch);
            continue;
        }
        let mut name   = String::new();
        let mut closed = false;
        for c in chars.by_ref() {
            if c == '}' { closed = true; break; }
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
            None => return Err(ParseError::Syntax {
                line: lineno,
                msg: format!("undefined parameter '{name}'"),
            }),
        }
    }
    Ok(result)
}

/// Map a single node: port names → call-site nets; ground stays "0"; all
/// others get the `{prefix}.` namespace.
fn remap_node(node: &str, port_map: &HashMap<String, String>, prefix: &str) -> String {
    if node == "0" {
        return "0".to_string();
    }
    if let Some(mapped) = port_map.get(node) {
        return mapped.clone();
    }
    format!("{prefix}.{node}")
}

/// Remap every node field (and the element name itself) in a flat element.
fn remap_element_nodes(
    el:       Element,
    port_map: &HashMap<String, String>,
    prefix:   &str,
) -> Element {
    let rn = |n: &str| remap_node(n, port_map, prefix);
    match el {
        Element::Resistor  { name, pos, neg, resistance } =>
            Element::Resistor  { name: format!("{prefix}.{name}"), pos: rn(&pos), neg: rn(&neg), resistance },
        Element::Capacitor { name, pos, neg, capacitance } =>
            Element::Capacitor { name: format!("{prefix}.{name}"), pos: rn(&pos), neg: rn(&neg), capacitance },
        Element::Inductor  { name, pos, neg, inductance } =>
            Element::Inductor  { name: format!("{prefix}.{name}"), pos: rn(&pos), neg: rn(&neg), inductance },
        Element::VoltageSource { name, pos, neg, waveform } =>
            Element::VoltageSource { name: format!("{prefix}.{name}"), pos: rn(&pos), neg: rn(&neg), waveform },
        Element::CurrentSource { name, pos, neg, waveform } =>
            Element::CurrentSource { name: format!("{prefix}.{name}"), pos: rn(&pos), neg: rn(&neg), waveform },
        Element::Diode { name, anode, cathode, model_name } =>
            Element::Diode { name: format!("{prefix}.{name}"), anode: rn(&anode), cathode: rn(&cathode), model_name },
        Element::Mosfet { name, drain, gate, source, bulk, model_name, params } =>
            Element::Mosfet { name: format!("{prefix}.{name}"), drain: rn(&drain), gate: rn(&gate),
                              source: rn(&source), bulk: rn(&bulk), model_name, params },
        Element::XOsdi { name, nets, model_name, params } =>
            Element::XOsdi { name: format!("{prefix}.{name}"),
                             nets: nets.iter().map(|n| rn(n)).collect(), model_name, params },
    }
}

/// Expand one `.subckt` instance into a flat `Vec<Element>`.
///
/// `expanding` is the set of subckt names currently on the call stack (cycle
/// detection).  It is mutated in place and restored before returning.
#[allow(clippy::too_many_arguments)]
fn expand_instance(
    def_name:     &str,
    inst_name:    &str,
    call_nets:    &[String],
    call_params:  &[(String, f64)],
    def:          &SubcktDef,
    subckt_defs:  &HashMap<String, SubcktDef>,
    global_params: &HashMap<String, f64>,
    expanding:    &mut HashSet<String>,
    call_lineno:  usize,
) -> Result<Vec<Element>, ParseError> {
    // Port-count check.
    if call_nets.len() != def.ports.len() {
        return Err(ParseError::SubcktPortCount {
            name:     def_name.to_string(),
            expected: def.ports.len(),
            got:      call_nets.len(),
            line:     call_lineno,
        });
    }

    // Cycle detection.
    if expanding.contains(def_name) {
        return Err(ParseError::SubcktCycle { name: def_name.to_string() });
    }
    expanding.insert(def_name.to_string());

    // port_map: def port name → call-site net.
    let port_map: HashMap<String, String> = def.ports.iter()
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
        let lineno  = *lineno;
        let trimmed = body_line.trim();
        if trimmed.is_empty() || trimmed.starts_with('*') {
            continue;
        }
        let lc = trimmed.to_lowercase();

        // Directives consumed by collect_defs — skip.
        if lc == ".end" || lc.starts_with(".ends") || lc.starts_with(".subckt") || lc.starts_with(".param") {
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
                return Err(ParseError::UnsupportedDirective { directive, line: lineno });
            }
            continue;
        }

        let substituted = substitute_params(trimmed, &inst_params, lineno)?;
        let el = parse_element(&substituted, lineno)?;
        let el = remap_element_nodes(el, &port_map, inst_name);

        // Recurse if this element is a nested subckt instance.
        let is_subckt_inst = if let Element::XOsdi { ref model_name, .. } = el {
            subckt_defs.contains_key(model_name)
        } else {
            false
        };

        if is_subckt_inst {
            if let Element::XOsdi { ref name, ref nets, ref model_name, ref params } = el {
                let nested_def = subckt_defs.get(model_name).unwrap();
                let nested = expand_instance(
                    model_name, name, nets, params,
                    nested_def, subckt_defs, &inst_params, expanding, lineno,
                )?;
                result.extend(nested);
            }
        } else {
            result.push(el);
        }
    }

    expanding.remove(def_name);
    Ok(result)
}

/// Returns `true` for output-only directives that can be silently ignored.
fn is_silent_directive(lc: &str) -> bool {
    lc.starts_with(".print")
        || lc.starts_with(".plot")
        || lc.starts_with(".probe")
        || lc.starts_with(".measure")
        || lc.starts_with(".meas")
        || lc.starts_with(".backanno")
}

/// Parse `.ic V(n1)=val V(n2)=val ...` or `.nodeset V(n)=val ...`.
///
/// Returns a `Vec<(node, value)>` (node name lowercased, "gnd" canonicalised).
/// Tokens that don't match the `V(<name>)=<value>` shape are silently ignored
/// (they're typically the leading `.ic`/`.nodeset` keyword itself).
fn parse_node_assignments(line: &str) -> Result<Vec<(String, f64)>, ParseError> {
    let raw: String = line.split_whitespace().skip(1).collect::<Vec<_>>().join(" ");
    let mut out = Vec::new();

    // Pre-process: collapse spaces around `=` and inside `V(…)` so we can
    // tokenize on whitespace.
    let cleaned: String = raw.replace(" =", "=").replace("= ", "=");

    for tok in cleaned.split_whitespace() {
        // Accept either V(name)=val or name=val (and case-insensitively v / V).
        let (lhs, rhs) = match tok.split_once('=') {
            Some(p) => p,
            None => continue,
        };
        let lhs_lc = lhs.to_lowercase();
        let name = if let Some(inner) = lhs_lc.strip_prefix("v(").and_then(|s| s.strip_suffix(')')) {
            inner.to_string()
        } else {
            lhs_lc.clone()
        };
        let value: f64 = parse_value(rhs, 0).unwrap_or_else(|_| {
            rhs.parse::<f64>().unwrap_or(0.0)
        });
        out.push((canon_node(&name), value));
    }
    Ok(out)
}

/// Parse `.options key=val key=val ...` into a list of `(key, value)` pairs.
///
/// Bare-flag tokens (no `=`) are stored as `("key", "1")` so `SimOptions::set`
/// can treat them as boolean true.  Quoted values are stripped of surrounding
/// quotes.  Returns an empty list for an empty directive line.
fn parse_options_directive(line: &str) -> Vec<(String, String)> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    let mut pairs = Vec::new();
    for tok in &tokens[1..] {
        if let Some((k, v)) = tok.split_once('=') {
            let v = v.trim_matches('"').trim_matches('\'');
            pairs.push((k.to_lowercase(), v.to_string()));
        } else if !tok.is_empty() {
            pairs.push((tok.to_lowercase(), "1".to_string()));
        }
    }
    pairs
}

// ─── public entry point ───────────────────────────────────────────────────────

/// Parse a SPICE netlist from a string.
///
/// Runs in two passes:
///   1. Collect `.subckt` definitions and `.param` values.
///   2. Parse the main body, expanding subckt instances inline.
pub fn parse_spice(input: &str) -> Result<Netlist, ParseError> {
    let all_lines = logical_lines(input);
    let mut netlist = Netlist::default();

    if all_lines.is_empty() {
        return Ok(netlist);
    }

    // First logical line is the title.
    netlist.title = all_lines[0].1.trim().to_string();

    // Pass 1.
    let (subckt_defs, global_params, main_lines) = collect_defs(&all_lines[1..])?;

    // Pass 2: parse main body.
    let mut expanding: HashSet<String> = HashSet::new();

    for (lineno, line) in &main_lines {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('*') {
            continue;
        }

        let lc = trimmed.to_lowercase();

        if lc == ".end" {
            break;
        } else if lc.starts_with(".optical_bus") {
            // .optical_bus N re_base im_base wl_base
            // Declares an N-channel WDM optical bus, generating 3N net entries:
            //   re_base_0 im_base_0 wl_base_0  re_base_1 im_base_1 wl_base_1  ...
            // Must come before the ".optical" check.
            let tokens: Vec<&str> = trimmed.split_whitespace().collect();
            if tokens.len() >= 5 {
                if let Ok(n) = tokens[1].parse::<usize>() {
                    let (rb, ib, wb) = (tokens[2], tokens[3], tokens[4]);
                    for i in 0..n {
                        netlist.optical_nets.push(canon_node(&format!("{rb}_{i}")));
                        netlist.optical_nets.push(canon_node(&format!("{ib}_{i}")));
                        netlist.optical_nets.push(canon_node(&format!("{wb}_{i}")));
                    }
                }
            }
        } else if lc.starts_with(".optical") {
            // Must come before ".op" check.
            // Supports bus-vector notation: .optical net[0..3] expands to net_0 net_1 net_2 net_3
            let tokens: Vec<&str> = trimmed.split_whitespace().collect();
            for tok in &tokens[1..] {
                for net in expand_bus_vectors(tok) {
                    netlist.optical_nets.push(canon_node(&net));
                }
            }
        } else if lc.starts_with(".options") || lc.starts_with(".option") {
            netlist.options.extend(parse_options_directive(trimmed));
        } else if lc.starts_with(".op") {
            netlist.analyses.push(Analysis::Op);
        } else if lc.starts_with(".tran") {
            netlist.analyses.push(parse_tran(&lc, *lineno)?);
        } else if lc.starts_with(".ac") {
            netlist.analyses.push(parse_ac(&lc, *lineno)?);
        } else if lc.starts_with(".dc") {
            netlist.analyses.push(parse_dc(&lc, *lineno)?);
        } else if lc.starts_with(".model") {
            if let Some(card) = parse_model(&lc, *lineno)? {
                netlist.models.push(card);
            }
        } else if lc.starts_with(".osdi") {
            let tokens: Vec<&str> = trimmed.split_whitespace().collect();
            if tokens.len() >= 2 {
                netlist.osdi_paths.push(tokens[1].to_string());
            }
        } else if lc.starts_with(".ic") {
            netlist.ic.extend(parse_node_assignments(trimmed)?);
        } else if lc.starts_with(".nodeset") {
            netlist.nodeset.extend(parse_node_assignments(trimmed)?);
        } else if lc.starts_with('.') {
            if !is_silent_directive(&lc) {
                let directive = lc.split_whitespace().next().unwrap_or(&lc).to_string();
                return Err(ParseError::UnsupportedDirective { directive, line: *lineno });
            }
        } else {
            // Element or instance line; substitute top-level params first.
            let substituted = substitute_params(trimmed, &global_params, *lineno)?;
            let el = parse_element(&substituted, *lineno)?;

            let is_subckt_inst = if let Element::XOsdi { ref model_name, .. } = el {
                subckt_defs.contains_key(model_name)
            } else {
                false
            };

            if is_subckt_inst {
                if let Element::XOsdi { ref name, ref nets, ref model_name, ref params } = el {
                    let def  = subckt_defs.get(model_name).unwrap();
                    let flat = expand_instance(
                        model_name, name, nets, params,
                        def, &subckt_defs, &global_params, &mut expanding, *lineno,
                    )?;
                    netlist.elements.extend(flat);
                }
            } else {
                netlist.elements.push(el);
            }
        }
    }

    Ok(netlist)
}

/// Recursively expand `.include "file"` lines, substituting them with file
/// content.  `base_dir` is used to resolve relative paths.
fn resolve_includes(
    input:    &str,
    base_dir: Option<&Path>,
    depth:    usize,
) -> Result<String, ParseError> {
    if depth > 16 {
        return Err(ParseError::Syntax {
            line: 0,
            msg: ".include nesting depth > 16 (circular include?)".into(),
        });
    }
    let mut out = String::with_capacity(input.len());
    for (i, raw) in input.lines().enumerate() {
        let lineno = i + 1;
        let lc = raw.trim().to_lowercase();
        if lc.starts_with(".include") {
            let tok: Vec<&str> = raw.trim().splitn(2, char::is_whitespace).collect();
            if tok.len() < 2 {
                return Err(ParseError::Syntax {
                    line: lineno,
                    msg: ".include requires a filename argument".into(),
                });
            }
            let fname = tok[1].trim().trim_matches('"').trim_matches('\'');
            let path: PathBuf = match base_dir {
                Some(dir) => dir.join(fname),
                None      => PathBuf::from(fname),
            };
            let content = std::fs::read_to_string(&path).map_err(|e| ParseError::Syntax {
                line: lineno,
                msg: format!(".include '{}': {e}", path.display()),
            })?;
            let inlined = resolve_includes(&content, path.parent(), depth + 1)?;
            out.push_str(&inlined);
            out.push('\n');
        } else {
            out.push_str(raw);
            out.push('\n');
        }
    }
    Ok(out)
}

/// Parse a SPICE netlist file, resolving `.include` directives relative to
/// the file's parent directory.
pub fn parse_spice_file(path: &Path) -> Result<Netlist, ParseError> {
    let src = std::fs::read_to_string(path).map_err(|e| ParseError::Syntax {
        line: 0,
        msg: format!("cannot read '{}': {e}", path.display()),
    })?;
    let resolved = resolve_includes(&src, path.parent(), 0)?;
    parse_spice(&resolved)
}

// ─── analysis directive parsers ───────────────────────────────────────────────

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

/// Parse `.dc SRC START STOP STEP [SRC2 START2 STOP2 STEP2]`.
fn parse_dc(line: &str, lineno: usize) -> Result<Analysis, ParseError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 5 {
        return Err(ParseError::FieldCount {
            expected: "≥5 (.dc SRC START STOP STEP)",
            got: tokens.len(),
            line: lineno,
        });
    }
    let src   = tokens[1].to_lowercase();
    let start = parse_value(tokens[2], lineno)?;
    let stop  = parse_value(tokens[3], lineno)?;
    let step  = parse_value(tokens[4], lineno)?;
    let nested = if tokens.len() >= 9 {
        Some(DcSweepSpec {
            src:   tokens[5].to_lowercase(),
            start: parse_value(tokens[6], lineno)?,
            stop:  parse_value(tokens[7], lineno)?,
            step:  parse_value(tokens[8], lineno)?,
        })
    } else {
        None
    };
    Ok(Analysis::Dc { src, start, stop, step, nested })
}

fn parse_ac(line: &str, lineno: usize) -> Result<Analysis, ParseError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 5 {
        return Err(ParseError::FieldCount {
            expected: "≥5 (.ac DEC|OCT|LIN points fstart fstop)",
            got: tokens.len(),
            line: lineno,
        });
    }
    let variation = match tokens[1] {
        "dec" => AcVariation::Dec,
        "oct" => AcVariation::Oct,
        "lin" => AcVariation::Lin,
        other => return Err(ParseError::Syntax {
            line: lineno,
            msg: format!("unknown AC variation '{other}'; expected dec, oct, or lin"),
        }),
    };
    let points = tokens[2].parse::<usize>().map_err(|_| ParseError::Syntax {
        line: lineno,
        msg: format!("invalid point count '{}'", tokens[2]),
    })?;
    Ok(Analysis::Ac {
        variation,
        points,
        fstart: parse_value(tokens[3], lineno)?,
        fstop:  parse_value(tokens[4], lineno)?,
    })
}

// ─── element and value parsers ────────────────────────────────────────────────

/// Canonicalise a node name: lowercase; "gnd" → "0".
fn canon_node(s: &str) -> String {
    let s = s.to_lowercase();
    if s == "gnd" { "0".to_string() } else { s }
}

/// Parse an SPICE suffix (k, meg, m, u, n, p, f, g, t) into a float.
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

/// Expand a bus-vector token like `net[M..N]` into individual net names
/// `net_M, net_{M+1}, ..., net_N` (inclusive, underscore-separated).
/// If the token contains no `[M..N]` notation, returns the token unchanged
/// in a single-element vec.
fn expand_bus_vectors(token: &str) -> Vec<String> {
    if let (Some(lb), Some(rb)) = (token.find('['), token.rfind(']')) {
        if lb < rb {
            let base      = &token[..lb];
            let range_str = &token[lb + 1..rb];
            if let Some((lo_s, hi_s)) = range_str.split_once("..") {
                if let (Ok(lo), Ok(hi)) =
                    (lo_s.trim().parse::<usize>(), hi_s.trim().parse::<usize>())
                {
                    return (lo..=hi).map(|i| format!("{}_{}", base, i)).collect();
                }
            }
        }
    }
    vec![token.to_string()]
}

fn parse_element(line: &str, lineno: usize) -> Result<Element, ParseError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    let name   = tokens[0].to_lowercase();
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
                anode:      canon_node(tokens[1]),
                cathode:    canon_node(tokens[2]),
                model_name: tokens[3].to_lowercase(),
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
                drain:      canon_node(tokens[1]),
                gate:       canon_node(tokens[2]),
                source:     canon_node(tokens[3]),
                bulk:       canon_node(tokens[4]),
                model_name: tokens[5].to_lowercase(),
                params,
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
            let mut positional: Vec<&str>      = Vec::new();
            let mut params:     Vec<(String, f64)> = Vec::new();
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
            Ok(Element::XOsdi { name, nets, model_name, params })
        }
        _ => Err(ParseError::UnknownElement { letter, line: lineno }),
    }
}

fn parse_model(line: &str, lineno: usize) -> Result<Option<ModelCard>, ParseError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 3 {
        return Ok(None);
    }
    let name = tokens[1].to_string();
    let kind = tokens[2].to_lowercase();
    let rest = tokens[3..].join(" ");
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

fn parse_waveform(tokens: &[&str], lineno: usize) -> Result<Waveform, ParseError> {
    if tokens.len() < 4 {
        return Err(ParseError::FieldCount { expected: "≥4", got: tokens.len(), line: lineno });
    }
    let rest    = tokens[3..].join(" ");
    let rest_lc = rest.to_lowercase();

    if rest_lc.starts_with("pulse") { return parse_pulse(&rest_lc, lineno); }
    if rest_lc.starts_with("pwl")   { return parse_pwl(&rest_lc, lineno); }
    if rest_lc.starts_with("sin")   { return parse_sin(&rest_lc, lineno); }
    if rest_lc.starts_with("exp")   { return parse_exp(&rest_lc, lineno); }
    if rest_lc.starts_with("sffm")  { return parse_sffm(&rest_lc, lineno); }
    if rest_lc.starts_with("am")    { return parse_am(&rest_lc, lineno); }

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

/// Extract the parenthesised parameter list `(p0 p1 ...)` after a waveform
/// keyword, returning the tokens.  The keyword (e.g. `sin`, `exp`) precedes
/// the parens; trailing tokens after `)` are discarded.
fn parens_tokens<'a>(s: &'a str, kind: &str, lineno: usize) -> Result<Vec<&'a str>, ParseError> {
    let start = s.find('(').ok_or_else(|| ParseError::Syntax {
        line: lineno, msg: format!("{}: missing '('", kind.to_uppercase()) })?;
    let end = s.rfind(')').ok_or_else(|| ParseError::Syntax {
        line: lineno, msg: format!("{}: missing ')'", kind.to_uppercase()) })?;
    Ok(s[start + 1..end].split_whitespace().collect())
}

fn parse_sin(s: &str, lineno: usize) -> Result<Waveform, ParseError> {
    let parts = parens_tokens(s, "sin", lineno)?;
    let get = |i: usize, default: f64| -> Result<f64, ParseError> {
        parts.get(i).map_or(Ok(default), |t| parse_value(t, lineno))
    };
    if parts.len() < 3 {
        return Err(ParseError::FieldCount {
            expected: "≥3 (SIN vo va freq …)", got: parts.len(), line: lineno });
    }
    Ok(Waveform::Sin {
        vo:    get(0, 0.0)?,
        va:    get(1, 0.0)?,
        freq:  get(2, 0.0)?,
        td:    get(3, 0.0)?,
        theta: get(4, 0.0)?,
        phase: get(5, 0.0)?,
    })
}

fn parse_exp(s: &str, lineno: usize) -> Result<Waveform, ParseError> {
    let parts = parens_tokens(s, "exp", lineno)?;
    let get = |i: usize, default: f64| -> Result<f64, ParseError> {
        parts.get(i).map_or(Ok(default), |t| parse_value(t, lineno))
    };
    if parts.len() < 2 {
        return Err(ParseError::FieldCount {
            expected: "≥2 (EXP v1 v2 …)", got: parts.len(), line: lineno });
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
    Ok(Waveform::Exp { v1, v2, td1, tau1, td2, tau2 })
}

fn parse_sffm(s: &str, lineno: usize) -> Result<Waveform, ParseError> {
    let parts = parens_tokens(s, "sffm", lineno)?;
    let get = |i: usize, default: f64| -> Result<f64, ParseError> {
        parts.get(i).map_or(Ok(default), |t| parse_value(t, lineno))
    };
    if parts.len() < 5 {
        return Err(ParseError::FieldCount {
            expected: "5 (SFFM vo va fc mdi fs)", got: parts.len(), line: lineno });
    }
    Ok(Waveform::Sffm {
        vo:  get(0, 0.0)?,
        va:  get(1, 0.0)?,
        fc:  get(2, 0.0)?,
        mdi: get(3, 0.0)?,
        fs:  get(4, 0.0)?,
    })
}

fn parse_am(s: &str, lineno: usize) -> Result<Waveform, ParseError> {
    let parts = parens_tokens(s, "am", lineno)?;
    let get = |i: usize, default: f64| -> Result<f64, ParseError> {
        parts.get(i).map_or(Ok(default), |t| parse_value(t, lineno))
    };
    if parts.len() < 4 {
        return Err(ParseError::FieldCount {
            expected: "≥4 (AM va vo mf fc [td])", got: parts.len(), line: lineno });
    }
    Ok(Waveform::Am {
        va: get(0, 0.0)?,
        vo: get(1, 0.0)?,
        mf: get(2, 0.0)?,
        fc: get(3, 0.0)?,
        td: get(4, 0.0)?,
    })
}

fn parse_pulse(s: &str, lineno: usize) -> Result<Waveform, ParseError> {
    let start = s.find('(').ok_or_else(|| ParseError::Syntax { line: lineno, msg: "PULSE: missing '('".into() })?;
    let end   = s.rfind(')').ok_or_else(|| ParseError::Syntax { line: lineno, msg: "PULSE: missing ')'".into() })?;
    let inner = &s[start + 1..end];
    let parts: Vec<&str> = inner.split_whitespace().collect();
    let get   = |i: usize, default: f64| -> Result<f64, ParseError> {
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

fn parse_pwl(s: &str, lineno: usize) -> Result<Waveform, ParseError> {
    let inner = if let Some(start) = s.find('(') {
        let end = s.rfind(')').ok_or_else(|| ParseError::Syntax { line: lineno, msg: "PWL: missing ')'".into() })?;
        &s[start + 1..end]
    } else {
        s.strip_prefix("pwl").unwrap_or("").trim()
    };
    let values: Vec<f64> = inner
        .split_whitespace()
        .map(|tok| parse_value(tok, lineno))
        .collect::<Result<_, _>>()?;
    if values.len() < 2 || values.len() % 2 != 0 {
        return Err(ParseError::Syntax {
            line: lineno,
            msg: format!("PWL requires an even number of values (t v pairs), got {}", values.len()),
        });
    }
    Ok(Waveform::Pwl { points: values.chunks_exact(2).map(|p| (p[0], p[1])).collect() })
}

// ─── logical-line joiner ──────────────────────────────────────────────────────

/// Join continuation lines and return `(original_lineno, joined_line)` pairs.
fn logical_lines(input: &str) -> Vec<(usize, String)> {
    let mut result: Vec<(usize, String)> = Vec::new();
    for (i, raw) in input.lines().enumerate() {
        let lineno  = i + 1;
        let trimmed = raw.trim_start();

        if trimmed.starts_with('+') {
            if let Some(last) = result.last_mut() {
                last.1.push(' ');
                last.1.push_str(trimmed[1..].trim());
            }
            continue;
        }

        let prev_ends_backslash = result
            .last()
            .map(|(_, s)| s.trim_end().ends_with('\\'))
            .unwrap_or(false);

        if prev_ends_backslash {
            if let Some(last) = result.last_mut() {
                let without_bs = last.1.trim_end().trim_end_matches('\\').trim_end().to_string();
                last.1 = without_bs;
                last.1.push(' ');
                last.1.push_str(trimmed);
            }
        } else {
            result.push((lineno, raw.to_string()));
        }
    }
    result
}

// ─── tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── existing tests (unchanged) ────────────────────────────────────────────

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
        assert!((parse_value("2k", 1).unwrap() - 2000.0).abs() < 1e-9);
    }

    #[test]
    fn parse_suffix_meg() {
        assert!((parse_value("1meg", 1).unwrap() - 1e6).abs() < 1.0);
    }

    #[test]
    fn gnd_canonical() {
        assert_eq!(canon_node("GND"), "0");
        assert_eq!(canon_node("gnd"), "0");
        assert_eq!(canon_node("0"),   "0");
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
        let m  = &netlist.models[0];
        assert_eq!(m.name, "myd");
        assert_eq!(m.kind, "d");
        let is = m.params.iter().find(|(k, _)| k == "is").map(|(_, v)| *v).unwrap();
        assert!((is - 1e-14).abs() < 1e-20, "is={is}");
        let n  = m.params.iter().find(|(k, _)| k == "n").map(|(_, v)| *v).unwrap();
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
        assert!((w.at(1e-9) - 1.0).abs() < 1e-6);
        assert!((w.at(0.5e-9) - 0.5).abs() < 1e-6);
    }

    #[test]
    fn parse_pwl_waveform() {
        let input = "* PWL\nV1 a 0 PWL(0 0 1u 5 2u 5 3u 0)\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        if let Element::VoltageSource { waveform: Waveform::Pwl { points }, .. } = &netlist.elements[0] {
            assert_eq!(points.len(), 4);
            assert!((points[1].0 - 1e-6).abs() < 1e-18);
            assert!((points[1].1 - 5.0).abs() < 1e-12);
        } else {
            panic!("expected PWL VoltageSource");
        }
    }

    #[test]
    fn pwl_waveform_at() {
        let w = Waveform::Pwl { points: vec![(0.0, 0.0), (1e-6, 5.0), (2e-6, 5.0), (3e-6, 0.0)] };
        assert!((w.at(0.0) - 0.0).abs() < 1e-12);
        assert!((w.at(0.5e-6) - 2.5).abs() < 1e-9);
        assert!((w.at(1.5e-6) - 5.0).abs() < 1e-9);
        assert!((w.at(2.5e-6) - 2.5).abs() < 1e-9);
        assert!((w.at(4e-6) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn parse_ac_directive() {
        let input = "* RC\nV1 in 0 DC 1\nR1 in out 1k\nC1 out 0 1u\n.ac dec 20 1 100k\n.end\n";
        let netlist = parse_spice(input).unwrap();
        match &netlist.analyses[0] {
            Analysis::Ac { variation, points, fstart, fstop } => {
                assert_eq!(*variation, crate::AcVariation::Dec);
                assert_eq!(*points, 20);
                assert!((fstart - 1.0).abs() < 1e-12);
                assert!((fstop - 1e5).abs() < 1e-6);
            }
            _ => panic!("expected Ac analysis"),
        }
    }

    #[test]
    fn parse_ac_lin() {
        let input = "* test\nV1 in 0 DC 1\n.ac lin 100 1k 10k\n.end\n";
        let netlist = parse_spice(input).unwrap();
        match &netlist.analyses[0] {
            Analysis::Ac { variation, points, .. } => {
                assert_eq!(*variation, crate::AcVariation::Lin);
                assert_eq!(*points, 100);
            }
            _ => panic!("expected Ac analysis"),
        }
    }

    #[test]
    fn pulse_next_breakpoint_before_td() {
        let w = Waveform::Pulse { v0: 0.0, v1: 1.0, td: 1e-6, tr: 100e-9, tf: 100e-9, pw: 5e-6, per: 10e-6 };
        let bp = w.next_breakpoint(0.0).unwrap();
        assert!((bp - 1e-6).abs() < 1e-18, "expected td=1µs, got {bp}");
    }

    #[test]
    fn pulse_next_breakpoint_at_period_boundary() {
        let td = 0.0_f64; let tr = 100e-9_f64; let pw = 5e-6_f64; let tf = 100e-9_f64; let per = 10e-6_f64;
        let w  = Waveform::Pulse { v0: 0.0, v1: 1.0, td, tr, tf, pw, per };
        let t  = td + per;
        let bp = w.next_breakpoint(t).unwrap();
        assert!((bp - (t + tr)).abs() < 1e-18, "expected t+tr={}, got {bp}", t + tr);
    }

    #[test]
    fn pulse_next_breakpoint_mid_rise() {
        let w  = Waveform::Pulse { v0: 0.0, v1: 1.0, td: 0.0, tr: 100e-9, tf: 100e-9, pw: 5e-6, per: 10e-6 };
        let bp = w.next_breakpoint(50e-9).unwrap();
        assert!((bp - 100e-9).abs() < 1e-18, "expected 100ns, got {bp}");
    }

    #[test]
    fn pulse_next_breakpoint_no_repeat_exhausted() {
        let w = Waveform::Pulse { v0: 0.0, v1: 1.0, td: 0.0, tr: 100e-9, tf: 100e-9, pw: 5e-6, per: 0.0 };
        let after_all = 100e-9 + 5e-6 + 100e-9 + 1e-9;
        assert!(w.next_breakpoint(after_all).is_none());
    }

    #[test]
    fn pwl_next_breakpoint() {
        let w = Waveform::Pwl { points: vec![(0.0, 0.0), (1e-6, 5.0), (2e-6, 5.0), (3e-6, 0.0)] };
        assert!((w.next_breakpoint(-1.0).unwrap() - 0.0).abs() < 1e-18);
        assert!((w.next_breakpoint(0.0).unwrap() - 1e-6).abs() < 1e-18);
        assert!((w.next_breakpoint(1e-6).unwrap() - 2e-6).abs() < 1e-18);
        assert!(w.next_breakpoint(3e-6).is_none());
    }

    #[test]
    fn dc_waveform_no_breakpoints() {
        let w = Waveform::Dc(5.0);
        assert!(w.next_breakpoint(0.0).is_none());
        assert!(w.next_breakpoint(1e6).is_none());
    }

    #[test]
    fn parse_xosdi_element() {
        let input = "* photonic test\n\
                     Xlaser laser_re laser_im cw_laser power_mW=1.0\n\
                     .op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(netlist.elements.len(), 1);
        if let Element::XOsdi { name, nets, model_name, params } = &netlist.elements[0] {
            assert_eq!(name, "xlaser");
            assert_eq!(nets, &["laser_re", "laser_im"]);
            assert_eq!(model_name, "cw_laser");
            assert_eq!(params.len(), 1);
            assert_eq!(params[0].0, "power_mw");
            assert!((params[0].1 - 1.0).abs() < 1e-12);
        } else {
            panic!("expected XOsdi element");
        }
    }

    #[test]
    fn parse_sin_waveform() {
        let input = "* sin\nV1 a 0 SIN(0 1 1k 0 0 0)\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        match &netlist.elements[0] {
            Element::VoltageSource { waveform: Waveform::Sin { vo, va, freq, .. }, .. } => {
                assert!((vo - 0.0).abs() < 1e-12);
                assert!((va - 1.0).abs() < 1e-12);
                assert!((freq - 1e3).abs() < 1e-6);
            }
            _ => panic!("expected SIN"),
        }
    }

    #[test]
    fn parse_exp_waveform() {
        let input = "* exp\nV1 a 0 EXP(0 1 1u 1u 5u 1u)\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        match &netlist.elements[0] {
            Element::VoltageSource { waveform: Waveform::Exp { v1, v2, td1, tau1, td2, tau2 }, .. } => {
                assert_eq!(*v1, 0.0);
                assert!((v2 - 1.0).abs() < 1e-12);
                assert!((td1 - 1e-6).abs() < 1e-12);
                assert!((tau1 - 1e-6).abs() < 1e-12);
                assert!((td2 - 5e-6).abs() < 1e-12);
                assert!((tau2 - 1e-6).abs() < 1e-12);
            }
            _ => panic!("expected EXP"),
        }
    }

    #[test]
    fn parse_sffm_waveform() {
        let input = "* sffm\nV1 a 0 SFFM(0 1 1k 5 100)\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        match &netlist.elements[0] {
            Element::VoltageSource { waveform: Waveform::Sffm { vo, va, fc, mdi, fs }, .. } => {
                assert_eq!(*vo, 0.0);
                assert!((va - 1.0).abs() < 1e-12);
                assert!((fc - 1e3).abs() < 1e-6);
                assert!((mdi - 5.0).abs() < 1e-12);
                assert!((fs - 100.0).abs() < 1e-9);
            }
            _ => panic!("expected SFFM"),
        }
    }

    #[test]
    fn parse_am_waveform() {
        let input = "* am\nV1 a 0 AM(1 0 100 1k 0)\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        match &netlist.elements[0] {
            Element::VoltageSource { waveform: Waveform::Am { va, vo, mf, fc, td }, .. } => {
                assert_eq!(*va, 1.0);
                assert_eq!(*vo, 0.0);
                assert!((mf - 100.0).abs() < 1e-9);
                assert!((fc - 1e3).abs() < 1e-6);
                assert!((td - 0.0).abs() < 1e-12);
            }
            _ => panic!("expected AM"),
        }
    }

    #[test]
    fn sin_at_zero_returns_vo() {
        let w = Waveform::Sin { vo: 0.5, va: 1.0, freq: 1e3, td: 0.0, theta: 0.0, phase: 0.0 };
        // sin(0) = 0 → vo
        assert!((w.at(0.0) - 0.5).abs() < 1e-12);
        // sin(π/2) at t = 0.25 ms with f=1kHz: sin(π/2)=1
        assert!((w.at(0.25e-3) - 1.5).abs() < 1e-6);
    }

    #[test]
    fn sin_pre_delay_is_vo() {
        let w = Waveform::Sin { vo: 0.5, va: 1.0, freq: 1e3, td: 1e-6, theta: 0.0, phase: 0.0 };
        assert!((w.at(0.0) - 0.5).abs() < 1e-12);
        assert!((w.at(0.5e-6) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn exp_pre_td1_is_v1() {
        let w = Waveform::Exp { v1: 0.0, v2: 1.0, td1: 1e-6, tau1: 1e-6, td2: 5e-6, tau2: 1e-6 };
        assert!((w.at(0.0) - 0.0).abs() < 1e-12);
        // At t=td1+tau1 the rise is (1-e^-1) = 0.6321
        assert!((w.at(2e-6) - (1.0 - (-1.0_f64).exp())).abs() < 1e-9);
    }

    #[test]
    fn parse_dc_sweep_single() {
        let input = "* dc\nV1 in 0 DC 0\nR1 in 0 1k\n.dc V1 0 5 0.1\n.end\n";
        let netlist = parse_spice(input).unwrap();
        match &netlist.analyses[0] {
            Analysis::Dc { src, start, stop, step, nested } => {
                assert_eq!(src, "v1");
                assert!((start - 0.0).abs() < 1e-12);
                assert!((stop - 5.0).abs() < 1e-12);
                assert!((step - 0.1).abs() < 1e-12);
                assert!(nested.is_none());
            }
            _ => panic!("expected Dc analysis"),
        }
    }

    #[test]
    fn parse_dc_sweep_nested() {
        let input = "* dc 2d\nV1 in 0 DC 0\nV2 g 0 DC 0\n.dc V1 0 5 0.5 V2 0 2 0.5\n.end\n";
        let netlist = parse_spice(input).unwrap();
        match &netlist.analyses[0] {
            Analysis::Dc { src, nested, .. } => {
                assert_eq!(src, "v1");
                let n = nested.as_ref().expect("nested sweep");
                assert_eq!(n.src, "v2");
                assert!((n.stop - 2.0).abs() < 1e-12);
            }
            _ => panic!("expected Dc analysis"),
        }
    }

    #[test]
    fn parse_options_directive_stores_pairs() {
        let input = "* opts\nV1 in 0 DC 1\nR1 in out 1k\n\
                     .options reltol=1e-5 gmin=1p method=be\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(netlist.options.len(), 3);
        assert_eq!(netlist.options[0], ("reltol".into(), "1e-5".into()));
        assert_eq!(netlist.options[1], ("gmin".into(), "1p".into()));
        assert_eq!(netlist.options[2], ("method".into(), "be".into()));
    }

    #[test]
    fn parse_options_accumulates_across_lines() {
        let input = "* opts\nV1 in 0 DC 1\n\
                     .options reltol=1e-5\n.options vntol=1e-9 itl1=300\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(netlist.options.len(), 3);
    }

    #[test]
    fn parse_options_bare_flag_is_true() {
        let input = "* opts\nV1 in 0 DC 1\n.options uic\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(netlist.options, vec![("uic".into(), "1".into())]);
    }

    #[test]
    fn parse_optical_directive() {
        let input = "* photonic test\n\
                     .optical laser_re laser_im wg_out_re wg_out_im\n\
                     .op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(netlist.optical_nets, vec!["laser_re", "laser_im", "wg_out_re", "wg_out_im"]);
    }

    #[test]
    fn bus_vector_expansion_in_optical() {
        // .optical with bus vector notation
        let input = "* WDM test\n\
                     .optical opt_re[0..2] opt_im[0..2] opt_wl[0..2]\n\
                     .op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(netlist.optical_nets, vec![
            "opt_re_0", "opt_re_1", "opt_re_2",
            "opt_im_0", "opt_im_1", "opt_im_2",
            "opt_wl_0", "opt_wl_1", "opt_wl_2",
        ]);
    }

    #[test]
    fn optical_bus_directive() {
        // .optical_bus N re_base im_base wl_base
        let input = "* WDM test\n\
                     .optical_bus 3 ch_re ch_im ch_wl\n\
                     .op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(netlist.optical_nets, vec![
            "ch_re_0", "ch_im_0", "ch_wl_0",
            "ch_re_1", "ch_im_1", "ch_wl_1",
            "ch_re_2", "ch_im_2", "ch_wl_2",
        ]);
    }

    #[test]
    fn bus_vector_expansion_in_xosdi_nets() {
        // X element with bus vector net arguments
        let input = "* WDM xosdi test\n\
                     .optical ch_re[0..1] ch_im[0..1] ch_wl[0..1]\n\
                     Xmux ch_re[0..1] ch_im[0..1] ch_wl[0..1] out_re out_im out_wl wdm_mux2\n\
                     .op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        if let Element::XOsdi { nets, model_name, .. } = &netlist.elements[0] {
            assert_eq!(model_name, "wdm_mux2");
            assert_eq!(nets, &[
                "ch_re_0", "ch_re_1",
                "ch_im_0", "ch_im_1",
                "ch_wl_0", "ch_wl_1",
                "out_re", "out_im", "out_wl",
            ]);
        } else {
            panic!("expected XOsdi");
        }
    }

    #[test]
    fn discipline_check_clean() {
        use crate::check_disciplines;
        let input = "* clean photonic circuit\n\
                     .optical laser_re laser_im\n\
                     Xlaser laser_re laser_im cw_laser\n\
                     R1 vdd 0 1k\n\
                     .op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert!(check_disciplines(&netlist).is_ok());
    }

    #[test]
    fn discipline_check_mismatch_resistor_on_optical_net() {
        use crate::{check_disciplines, DisciplineError};
        let input = "* BAD: resistor connected to optical net\n\
                     .optical laser_re laser_im\n\
                     R1 laser_re laser_im 50\n\
                     .op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        let err = check_disciplines(&netlist).unwrap_err();
        assert!(matches!(err, DisciplineError { .. }));
        assert_eq!(err.net, "laser_re");
    }

    #[test]
    fn discipline_check_xosdi_mixed_domain_allowed() {
        use crate::check_disciplines;
        let input = "* mixed-domain OK\n\
                     .optical opt_re opt_im\n\
                     Xpd opt_re opt_im ph_a ph_k photodetector\n\
                     R1 ph_a 0 1k\n\
                     .op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert!(check_disciplines(&netlist).is_ok());
    }

    // ── subckt tests ──────────────────────────────────────────────────────────

    #[test]
    fn subckt_basic_expansion() {
        // One instance of a simple resistor divider subckt.
        let input = "\
* Subckt basic
.subckt rdiv in out gnd_node
R1 in out 1k
R2 out gnd_node 1k
.ends rdiv
V1 vdd 0 DC 5
Xdiv1 vdd mid 0 rdiv
.op
.end
";
        let netlist = parse_spice(input).unwrap();
        // V1 + R1 (from Xdiv1) + R2 (from Xdiv1) = 3 elements
        assert_eq!(netlist.elements.len(), 3, "expected 3 flat elements");
        // Check that resistors have the correct namespaced names and nets.
        let names: Vec<&str> = netlist.elements.iter().map(|el| match el {
            Element::Resistor { name, .. } => name.as_str(),
            Element::VoltageSource { name, .. } => name.as_str(),
            _ => "?",
        }).collect();
        assert!(names.contains(&"xdiv1.r1"), "missing xdiv1.r1, got {names:?}");
        assert!(names.contains(&"xdiv1.r2"), "missing xdiv1.r2, got {names:?}");
        // Check node remapping: R1 should connect vdd → mid (port substitution).
        let r1 = netlist.elements.iter().find(|el| matches!(el, Element::Resistor { name, .. } if name == "xdiv1.r1")).unwrap();
        if let Element::Resistor { pos, neg, .. } = r1 {
            assert_eq!(pos, "vdd", "R1 pos should be vdd (port 'in' → vdd)");
            assert_eq!(neg, "mid", "R1 neg should be mid (port 'out' → mid)");
        }
    }

    #[test]
    fn subckt_two_instances() {
        // Two instances of the same subckt produce independent flat elements.
        let input = "\
* Two instances
.subckt inv a b
R1 a b 500
.ends inv
V1 vdd 0 DC 1
Xinv1 vdd n1 inv
Xinv2 vdd n2 inv
.op
.end
";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(netlist.elements.len(), 3); // V1 + 2 × R1
        let res: Vec<_> = netlist.elements.iter().filter(|el| matches!(el, Element::Resistor { .. })).collect();
        assert_eq!(res.len(), 2);
        let n1_r = res.iter().find(|el| {
            if let Element::Resistor { name, .. } = el { name.starts_with("xinv1.") } else { false }
        });
        assert!(n1_r.is_some(), "xinv1.r1 missing");
        if let Some(Element::Resistor { pos, neg, .. }) = n1_r {
            assert_eq!(pos, "vdd");
            assert_eq!(neg, "n1");
        }
    }

    #[test]
    fn subckt_param_default_and_override() {
        let input = "\
* Param default and override
.subckt rvar a b R=1k
R1 a b {R}
.ends rvar
V1 vdd 0 DC 1
Xdef  vdd 0 rvar
Xover vdd 0 rvar R=2k
.op
.end
";
        let netlist = parse_spice(input).unwrap();
        let resistors: Vec<_> = netlist.elements.iter().filter(|el| matches!(el, Element::Resistor { .. })).collect();
        assert_eq!(resistors.len(), 2);
        let r_def  = resistors.iter().find(|el| { if let Element::Resistor { name, .. } = el { name.starts_with("xdef.")  } else { false } }).unwrap();
        let r_over = resistors.iter().find(|el| { if let Element::Resistor { name, .. } = el { name.starts_with("xover.") } else { false } }).unwrap();
        if let Element::Resistor { resistance, .. } = r_def  { assert!((resistance - 1e3).abs() < 1e-9, "default R={resistance}"); }
        if let Element::Resistor { resistance, .. } = r_over { assert!((resistance - 2e3).abs() < 1e-9, "override R={resistance}"); }
    }

    #[test]
    fn subckt_global_param_substitution() {
        let input = "\
* Global .param
.param Rval=4.7k
.subckt rbuf a b
R1 a b {Rval}
.ends rbuf
V1 in 0 DC 1
Xbuf in out rbuf
.op
.end
";
        let netlist = parse_spice(input).unwrap();
        let res = netlist.elements.iter().find(|el| matches!(el, Element::Resistor { .. })).unwrap();
        if let Element::Resistor { resistance, .. } = res {
            assert!((resistance - 4700.0).abs() < 1e-9, "global param R={resistance}");
        }
    }

    #[test]
    fn subckt_nested() {
        // Inner subckt nested inside outer.
        let input = "\
* Nested subckt
.subckt inner a b
R1 a b 100
.ends inner
.subckt outer a b
R2 a b 200
Xin a b inner
.ends outer
V1 vdd 0 DC 1
Xout vdd 0 outer
.op
.end
";
        let netlist = parse_spice(input).unwrap();
        // Flat: V1, xout.r2, xout.xin.r1
        let resistors: Vec<_> = netlist.elements.iter().filter(|el| matches!(el, Element::Resistor { .. })).collect();
        assert_eq!(resistors.len(), 2, "expected 2 flat resistors, got {}", resistors.len());
        let names: Vec<&str> = resistors.iter().map(|el| {
            if let Element::Resistor { name, .. } = el { name.as_str() } else { "" }
        }).collect();
        assert!(names.iter().any(|n| n.starts_with("xout.r2")),   "missing xout.r2: {names:?}");
        assert!(names.iter().any(|n| n.starts_with("xout.xin.")), "missing xout.xin.*: {names:?}");
    }

    #[test]
    fn subckt_cycle_detection() {
        let input = "\
* Cyclic subckt
.subckt cyc a b
Xself a b cyc
.ends cyc
V1 vdd 0 DC 1
Xcyc vdd 0 cyc
.op
.end
";
        let err = parse_spice(input).unwrap_err();
        assert!(matches!(err, ParseError::SubcktCycle { .. }), "expected SubcktCycle, got {err:?}");
    }

    #[test]
    fn subckt_wrong_port_count() {
        let input = "\
* Port count mismatch
.subckt twoport a b
R1 a b 1k
.ends twoport
V1 vdd 0 DC 1
Xbad vdd mid out extra twoport
.op
.end
";
        let err = parse_spice(input).unwrap_err();
        assert!(matches!(err, ParseError::SubcktPortCount { expected: 2, got: 4, .. }), "got {err:?}");
    }

    #[test]
    fn subckt_forward_reference() {
        // Subckt defined AFTER the instance line — collect_defs handles all defs
        // in pass 1 regardless of order.
        let input = "\
* Forward ref
V1 vdd 0 DC 1
Xfwd vdd 0 fwdmod
.op
.subckt fwdmod a b
R1 a b 1k
.ends fwdmod
.end
";
        let netlist = parse_spice(input).unwrap();
        let res = netlist.elements.iter().find(|el| matches!(el, Element::Resistor { .. }));
        assert!(res.is_some(), "forward-referenced subckt not expanded");
    }

    #[test]
    fn unsupported_directive_errors() {
        let cases = [
            "* test\nV1 a 0 DC 1\n.lib mylib.lib\n.op\n.end\n",
            "* test\nV1 a 0 DC 1\n.func myfn(x)=x*x\n.op\n.end\n",
        ];
        for netlist_str in &cases {
            let result = parse_spice(netlist_str);
            assert!(
                matches!(result, Err(ParseError::UnsupportedDirective { .. })),
                "expected UnsupportedDirective for: {netlist_str}, got: {result:?}",
            );
        }
    }

    #[test]
    fn parse_ic_directive() {
        let input = "* ic\nV1 a 0 DC 1\nR1 a out 1k\nC1 out 0 1u\n\
                     .ic V(out)=0.5 V(a)=1.0\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(netlist.ic.len(), 2);
        // Lookup map by name.
        let m: std::collections::HashMap<_, _> = netlist.ic.iter().cloned().collect();
        assert!((m["out"] - 0.5).abs() < 1e-12);
        assert!((m["a"] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn parse_nodeset_directive() {
        let input = "* nodeset\nV1 a 0 DC 1\nR1 a out 1k\n\
                     .nodeset V(out)=0.7\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(netlist.nodeset.len(), 1);
        assert_eq!(netlist.nodeset[0].0, "out");
        assert!((netlist.nodeset[0].1 - 0.7).abs() < 1e-12);
    }

    #[test]
    fn silent_directives_ignored() {
        let input = "\
* test
V1 a 0 DC 1
R1 a 0 1k
.print V(a)
.plot V(a)
.meas tran vmax MAX V(a)
.options reltol=1e-4
.op
.end
";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(netlist.elements.len(), 2);
    }

    #[test]
    fn ends_does_not_terminate_toplevel() {
        // Prior bug: `.ends` was treated the same as `.end` at the top level.
        // With collect_defs, `.ends` at top level is a hard error.
        let input = "* test\nV1 a 0 DC 1\n.ends orphan\n.end\n";
        let err = parse_spice(input).unwrap_err();
        assert!(matches!(err, ParseError::Syntax { .. }), "expected Syntax error for stray .ends, got {err:?}");
    }

    #[test]
    fn subckt_osdi_coexistence() {
        // An X instance whose model_name is NOT in subckt_defs should remain
        // as XOsdi (i.e., fall through to the OSDI device path).
        let input = "\
* OSDI and subckt in same netlist
.subckt mybuf a b
R1 a b 50
.ends mybuf
Xbuf ina inb mybuf
Xlaser l_re l_im cw_laser power_mW=1.0
.op
.end
";
        let netlist = parse_spice(input).unwrap();
        // Xbuf expands to R1; Xlaser stays as XOsdi
        let r = netlist.elements.iter().filter(|el| matches!(el, Element::Resistor { .. })).count();
        let x = netlist.elements.iter().filter(|el| matches!(el, Element::XOsdi { .. })).count();
        assert_eq!(r, 1, "expected 1 resistor from subckt expansion");
        assert_eq!(x, 1, "expected 1 XOsdi remaining");
    }
}
