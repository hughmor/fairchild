use super::common::{canon_node, parse_value};
use super::directives::{expand_and_check, is_silent_directive, parse_func_directive, strip_span};
use super::element::{parse_element_expanded, parse_model};
use crate::expr::FuncTable;
use crate::warn_user;
use crate::{Element, ModelCard, ParseError};
use std::collections::{HashMap, HashSet};

type CollectDefsResult = (
    HashMap<String, SubcktDef>,
    HashMap<String, f64>,
    Vec<(usize, String)>,
    FuncTable,
);
type SubcktHeader = (String, Vec<String>, Vec<(String, f64)>);

/// What one `.subckt` instance flattens into: its elements, plus a private copy
/// of every `.model` card declared in its body (name-mangled per instance).
#[derive(Default)]
pub(super) struct Expansion {
    pub elements: Vec<Element>,
    pub models: Vec<ModelCard>,
}

// ─── internal types ──────────────────────────────────────────────────────────

/// Internal representation of a `.subckt ... .ends` block collected in pass 1.
pub(super) struct SubcktDef {
    ports: Vec<String>,               // port names (lowercased), in declaration order
    params: Vec<(String, f64)>,       // default parameter values (header + body .param)
    body_lines: Vec<(usize, String)>, // (original_lineno, raw_line) for pass-2 expansion
}

impl SubcktDef {
    /// Declared port count — used to pick bundle flatten-vs-replicate semantics.
    pub(super) fn port_count(&self) -> usize {
        self.ports.len()
    }
}

// ─── pass 1: collect definitions ─────────────────────────────────────────────

/// **Pass 1**: split logical lines into subckt definitions, global `.param`
/// values, and the main-body lines that pass 2 will parse.
///
/// Returns `(subckt_defs, global_params, main_lines)`.  Nested `.subckt`
/// definitions and a stray `.ends` are both hard errors.
///
/// `.control … .endc` is consumed here too, and warned about once.  It is
/// imperative script — `run`, `let`, `write`, loops, conditionals — and this is
/// not a shell: a second scripting language inside the simulator would do the job
/// the Python bindings exist for, worse.  Skipping the block costs nothing for
/// the large majority of real blocks, which hold only `run`/`write`/`plot`, all
/// three of which the frontend already does.  What it costs when the block is the
/// only place an analysis was declared is the warning's job to say.
pub(super) fn collect_defs(lines: &[(usize, String)]) -> Result<CollectDefsResult, ParseError> {
    let mut subckt_defs: HashMap<String, SubcktDef> = HashMap::new();
    let mut global_params: HashMap<String, f64> = HashMap::new();
    let mut main_lines: Vec<(usize, String)> = Vec::new();
    // `.func` definitions are collected across the whole deck before pass 2 runs,
    // so a call may precede its definition — more permissive than HSPICE, and one
    // less ordering rule to explain. A definition inside a `.subckt` still lands
    // here: scoping `.func` per subcircuit would need its own name resolution and
    // nothing has asked for it.
    let mut funcs: FuncTable = FuncTable::new();

    let mut in_subckt = false;
    // `.control` state. The verbs are collected only to name them in the
    // warning — nothing reads them, and nothing should start to.
    // `.if` / `.elseif` / `.else` / `.endif`. One frame per open `.if`; a line is
    // kept only while every frame on the stack is emitting.
    let mut cond_stack: Vec<CondFrame> = Vec::new();
    let mut in_control = false;
    let mut control_lineno = 0usize;
    let mut control_verbs: Vec<String> = Vec::new();
    let mut saw_control = false;
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

        // Inside a `.control` block nothing else applies: its lines are shell
        // commands that may look like elements or directives and are neither.
        if in_control {
            if lc.starts_with(".endc") {
                in_control = false;
            } else if let Some(verb) = lc.split_whitespace().next() {
                if !control_verbs.iter().any(|v| v == verb) {
                    control_verbs.push(verb.to_string());
                }
            }
            continue;
        }
        if lc.starts_with(".control") {
            in_control = true;
            saw_control = true;
            control_lineno = lineno;
            continue;
        }
        if lc.starts_with(".endc") {
            return Err(ParseError::Syntax {
                line: lineno,
                msg: ".endc without a matching .control".into(),
            });
        }

        // Conditionals resolve here, before any other arm sees the line: a false
        // branch may hold `.model`, `.subckt`, `.param` or elements, and none of
        // them may be collected. Conditions are evaluated against the parameters
        // and functions defined *above* them, in file order — which is why this
        // lives in the pass that accumulates those rather than in a pass of its
        // own.
        if let Some(kind) = conditional_directive(&lc) {
            if in_subckt {
                // Evaluating it here would use the subcircuit's *default*
                // parameters, once, for every instance — a wrong answer per
                // instance rather than a missing feature. Refuse instead.
                return Err(ParseError::Syntax {
                    line: lineno,
                    msg: format!(
                        "{kind} inside a .subckt is not supported: its condition \
                         would be evaluated once against the subcircuit's default \
                         parameters, not per instance. Select with .if outside the \
                         definition, or give the instances different models"
                    ),
                });
            }
            let parent_live = cond_stack.iter().all(|f| f.live);
            match kind {
                ".if" => {
                    let live =
                        parent_live && eval_condition(&lc, kind, &global_params, &funcs, lineno)?;
                    cond_stack.push(CondFrame { live, taken: live });
                }
                ".elseif" | ".else" | ".endif" => {
                    let Some(frame) = cond_stack.last().copied() else {
                        return Err(ParseError::Syntax {
                            line: lineno,
                            msg: format!("{kind} without a matching .if"),
                        });
                    };
                    if kind == ".endif" {
                        cond_stack.pop();
                    } else {
                        // `parent_live` above included this frame; the enclosing
                        // frames are what gate a sibling branch.
                        let outer_live = cond_stack[..cond_stack.len() - 1].iter().all(|f| f.live);
                        let cond = if kind == ".else" {
                            true
                        } else {
                            // A condition in a branch that cannot run is not
                            // evaluated at all, so a dead branch may reference
                            // parameters that were never defined.
                            outer_live
                                && !frame.taken
                                && eval_condition(&lc, kind, &global_params, &funcs, lineno)?
                        };
                        let live = outer_live && !frame.taken && cond;
                        let last = cond_stack.last_mut().unwrap();
                        last.live = live;
                        last.taken = frame.taken || live;
                    }
                }
                _ => unreachable!("conditional_directive returns only those four"),
            }
            continue;
        }
        if !cond_stack.iter().all(|f| f.live) {
            continue;
        }

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
        } else if lc.starts_with(".func") {
            let (name, def) = parse_func_directive(trimmed, lineno)?;
            funcs.insert(name, def);
        } else if lc.starts_with(".param") {
            // Values may be expressions over the parameters already in scope, so
            // the map in scope goes in: `.param a=2 b={a*3}` resolves left to
            // right, and inside a subcircuit the body's own defaults stack on top
            // of the globals.
            let mut in_scope = global_params.clone();
            if in_subckt {
                in_scope.extend(current_def.params.iter().cloned());
            }
            let pairs = parse_param_directive(trimmed, lineno, &in_scope, &funcs)?;
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
    if !cond_stack.is_empty() {
        return Err(ParseError::Syntax {
            line: 0,
            msg: format!(
                "{} unterminated .if block(s): every .if needs an .endif, and \
                 guessing where one ended would silently drop the rest of the deck",
                cond_stack.len()
            ),
        });
    }
    if in_control {
        // Not recoverable by guessing where the block ended: everything after it
        // would be silently discarded, which is how a deck loses half its
        // circuit and still runs.
        return Err(ParseError::Syntax {
            line: control_lineno,
            msg: ".control block has no matching .endc".into(),
        });
    }
    if saw_control {
        let verbs = if control_verbs.is_empty() {
            "empty".to_string()
        } else {
            control_verbs.join(", ")
        };
        warn_user!(
            ".control block skipped — its commands are not interpreted \
             ({verbs}). fairchild is not an ngspice shell: control flow belongs in \
             Python (fairchild.Circuit) or in CLI flags, and output selection is \
             --probe. An analysis that existed only inside the block will not run: \
             give the deck a .tran/.ac/.dc/.op card, or drive the run from Python \
             (see docs/spice_support.md §4.7)"
        );
    }

    Ok((subckt_defs, global_params, main_lines, funcs))
}

/// Parse `.subckt <name> <port1> ... [param=default ...]`.
pub(super) fn parse_subckt_header(line: &str, lineno: usize) -> Result<SubcktHeader, ParseError> {
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
    in_scope: &HashMap<String, f64>,
    funcs: &FuncTable,
) -> Result<Vec<(String, f64)>, ParseError> {
    let mut pairs: Vec<(String, f64)> = Vec::new();
    // Each value may use the ones before it, on this line or an earlier one.
    let mut scope = in_scope.clone();
    for tok in &split_assignments(line)[1..] {
        let Some((k, v)) = tok.split_once('=') else {
            continue;
        };
        let key = k.trim().to_lowercase();
        let raw = v.trim();
        // A plain number with a SPICE suffix stays the fast path: `1k` is a value,
        // not an expression (`k` would read as an undefined parameter).
        let val = match parse_value(raw, lineno) {
            Ok(n) => n,
            Err(_) => eval_param_expr(
                strip_span(raw),
                &scope,
                funcs,
                lineno,
                &format!(".param '{key}' = '{raw}'"),
            )?,
        };
        scope.insert(key.clone(), val);
        pairs.push((key, val));
    }
    Ok(pairs)
}

/// Split a directive into `key=value` tokens, keeping a braced or quoted value
/// whole and gluing an assignment back together when it was written with spaces
/// around the `=`.
///
/// A value may contain spaces only inside `{…}`, `'…'` or `"…"` — the same rule
/// HSPICE has, and the reason it has it: `.param a = 1 + 2 b = 3` cannot be split
/// unambiguously, so the braces are what make the intent readable.
fn split_assignments(line: &str) -> Vec<String> {
    let mut toks: Vec<String> = Vec::new();
    let mut cur = String::new();
    let mut span: Option<char> = None;
    for c in line.chars() {
        match span {
            Some(open) => {
                cur.push(c);
                let closes = match open {
                    '{' => c == '}',
                    q => c == q,
                };
                if closes {
                    span = None;
                }
            }
            None => match c {
                '{' | '\'' | '"' => {
                    span = Some(c);
                    cur.push(c);
                }
                ' ' | '\t' => {
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
    // Re-glue `name = value` and `name= value` / `name =value`.
    let mut glued: Vec<String> = Vec::new();
    let mut i = 0;
    while i < toks.len() {
        let mut t = toks[i].clone();
        if t == "=" {
            // Orphan `=`: attach to the previous token and pull in the next.
            if let Some(prev) = glued.pop() {
                t = format!("{prev}=");
            }
        }
        while t.ends_with('=') && i + 1 < toks.len() {
            i += 1;
            t.push_str(&toks[i]);
        }
        glued.push(t);
        i += 1;
    }
    glued
}

// ─── expansion helpers ────────────────────────────────────────────────────────

/// Evaluation context that resolves bare variables from a subckt's parameter
/// map, so `{…}` can hold arithmetic and not just a name. `pi` is supplied as a
/// constant because deriving a ring's length from its radius is the single most
/// common PCell expression.
struct ParamCtx<'a>(&'a HashMap<String, f64>);

impl crate::expr::EvalContext for ParamCtx<'_> {
    fn node_voltage(&self, _node: &str) -> f64 {
        f64::NAN // a parameter expression cannot reference the solution
    }
    fn branch_current(&self, _vsrc: &str) -> f64 {
        f64::NAN
    }
    fn time(&self) -> f64 {
        f64::NAN
    }
    fn variable(&self, name: &str) -> f64 {
        match name {
            "pi" => std::f64::consts::PI,
            _ => self.0.get(name).copied().unwrap_or(f64::NAN),
        }
    }
}

/// One open `.if`: whether its current branch emits, and whether any branch has.
#[derive(Clone, Copy)]
struct CondFrame {
    /// This branch is emitting lines.
    live: bool,
    /// Some branch of this `.if` has already been taken, so no later one may be.
    taken: bool,
}

/// Which conditional directive a line is, if any.
///
/// Matched on the whole first word so `.iff` or `.elsewhere` are not conditionals
/// — and `.if` must not swallow `.include`, which is why this is a word match and
/// not a prefix test.
fn conditional_directive(lc: &str) -> Option<&'static str> {
    let word = lc.split_whitespace().next().unwrap_or(lc);
    let word = word.split('(').next().unwrap_or(word);
    match word {
        ".if" => Some(".if"),
        ".elseif" | ".elsif" => Some(".elseif"),
        ".else" => Some(".else"),
        ".endif" => Some(".endif"),
        _ => None,
    }
}

/// Evaluate a `.if` / `.elseif` condition: true when it is a non-zero number.
///
/// The condition is an ordinary parse-time expression, so it may use `.param`
/// values, `.func` calls and the comparison operators the grammar already has.
/// Surrounding parentheses are optional — both spellings are in the wild.
fn eval_condition(
    lc: &str,
    kind: &str,
    params: &HashMap<String, f64>,
    funcs: &FuncTable,
    lineno: usize,
) -> Result<bool, ParseError> {
    let rest = lc[kind.len()..].trim();
    let rest = strip_span(rest).trim();
    if rest.is_empty() {
        return Err(ParseError::Syntax {
            line: lineno,
            msg: format!("{kind} has no condition"),
        });
    }
    let val = eval_param_expr(rest, params, funcs, lineno, &format!("{kind} condition"))?;
    Ok(val != 0.0)
}

/// Evaluate one parse-time expression over `params` and `funcs`.
///
/// The single path for every `{…}`, `'…'`, `.param` value and `.model` value:
/// parse, expand `.func` calls, refuse a call nothing can evaluate, then require
/// a finite result. An undefined name resolves to NaN inside the grammar, so the
/// finiteness check is what turns "you referred to something that does not exist"
/// into an error instead of a plausible number.
fn eval_param_expr(
    src: &str,
    params: &HashMap<String, f64>,
    funcs: &FuncTable,
    lineno: usize,
    what: &str,
) -> Result<f64, ParseError> {
    let parsed = crate::expr::Expr::parse(src).map_err(|e| ParseError::Syntax {
        line: lineno,
        msg: format!("{what}: cannot parse '{src}': {e}"),
    })?;
    let expr = expand_and_check(parsed, funcs, lineno, what)?;
    // Undefined names are caught by name, before evaluation: a comparison would
    // otherwise launder one into a valid answer — `nope == 1` is NaN against 1,
    // which is a perfectly finite `false`, and `.if (nope==1)` would quietly take
    // the other branch.
    let mut vars = Vec::new();
    expr.collect_vars(&mut vars);
    let undefined: Vec<String> = vars
        .into_iter()
        .filter(|v| v != "pi" && !params.contains_key(v))
        .collect();
    if !undefined.is_empty() {
        return Err(ParseError::Syntax {
            line: lineno,
            msg: format!(
                "{what}: undefined parameter(s) {} — define them with .param above \
                 this line, or (for a node voltage or the simulation time) note \
                 that a parse-time expression cannot use the solution",
                undefined
                    .iter()
                    .map(|v| format!("'{v}'"))
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        });
    }
    let val = expr.eval(&ParamCtx(params));
    if !val.is_finite() {
        return Err(ParseError::Syntax {
            line: lineno,
            msg: format!(
                "{what}: '{src}' did not evaluate to a finite number — an \
                 undefined parameter, or a reference to a node voltage or time, \
                 which a parse-time expression cannot use"
            ),
        });
    }
    Ok(val)
}

/// Replace `{…}` placeholders using `params`.
///
/// The contents may be a bare parameter name (`{radius}`) or any arithmetic
/// expression over parameters, numeric literals, `pi`, and the functions
/// `Expr::parse` supports (`{2*pi*radius}`, `{sqrt(p_pi/r_heater)}`). Errors on
/// an undefined name, an unparseable expression, or one that evaluates to a
/// non-finite value — a silently-NaN geometry is far worse than a parse error.
pub(super) fn substitute_params(
    line: &str,
    params: &HashMap<String, f64>,
    funcs: &FuncTable,
    lineno: usize,
) -> Result<String, ParseError> {
    if !line.contains('{') && !line.contains('\'') {
        return Ok(line.to_string());
    }
    let mut result = String::with_capacity(line.len() + 16);
    let mut chars = line.chars().peekable();
    while let Some(ch) = chars.next() {
        let close = match ch {
            '{' => '}',
            // HSPICE spells a parse-time expression `'…'`, and PDK model libraries
            // are written that way throughout. Double quotes stay untouched: on a
            // `.model` line those are a device constitutive map over the device's
            // own bias, which is not a parse-time value.
            '\'' => '\'',
            _ => {
                result.push(ch);
                continue;
            }
        };
        let mut body = String::new();
        let mut closed = false;
        for c in chars.by_ref() {
            if c == close {
                closed = true;
                break;
            }
            body.push(c);
        }
        if !closed {
            return Err(ParseError::Syntax {
                line: lineno,
                msg: format!("unclosed '{ch}' in parameter reference"),
            });
        }
        let key = body.trim().to_lowercase();
        // Fast path: a bare parameter name.
        if let Some(val) = params.get(&key) {
            result.push_str(&format!("{val:e}"));
            continue;
        }
        let val = eval_param_expr(
            &key,
            params,
            funcs,
            lineno,
            &format!("parameter expression '{body}'"),
        )?;
        result.push_str(&format!("{val:e}"));
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
            ac,
        } => Element::VoltageSource {
            name: format!("{prefix}.{name}"),
            pos: rn(&pos),
            neg: rn(&neg),
            waveform,
            ac,
        },
        Element::CurrentSource {
            name,
            pos,
            neg,
            waveform,
            ac,
        } => Element::CurrentSource {
            name: format!("{prefix}.{name}"),
            pos: rn(&pos),
            neg: rn(&neg),
            waveform,
            ac,
        },
        Element::Diode {
            name,
            anode,
            cathode,
            model_name,
            params,
        } => Element::Diode {
            name: format!("{prefix}.{name}"),
            anode: rn(&anode),
            cathode: rn(&cathode),
            model_name,
            params,
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
        Element::VoltageSwitch {
            name,
            pos,
            neg,
            ctrl_pos,
            ctrl_neg,
            model_name,
            initial_on,
        } => Element::VoltageSwitch {
            name: format!("{prefix}.{name}"),
            pos: rn(&pos),
            neg: rn(&neg),
            ctrl_pos: rn(&ctrl_pos),
            ctrl_neg: rn(&ctrl_neg),
            model_name,
            initial_on,
        },
        // `ctrl_vsrc` is an element name, so it is prefixed like the K
        // element's inductor references — not remapped like a net.
        Element::CurrentSwitch {
            name,
            pos,
            neg,
            ctrl_vsrc,
            model_name,
            initial_on,
        } => Element::CurrentSwitch {
            name: format!("{prefix}.{name}"),
            pos: rn(&pos),
            neg: rn(&neg),
            ctrl_vsrc: format!("{prefix}.{ctrl_vsrc}"),
            model_name,
            initial_on,
        },
        Element::TransmissionLine {
            name,
            a_pos,
            a_neg,
            b_pos,
            b_neg,
            z0,
            td,
        } => Element::TransmissionLine {
            name: format!("{prefix}.{name}"),
            a_pos: rn(&a_pos),
            a_neg: rn(&a_neg),
            b_pos: rn(&b_pos),
            b_neg: rn(&b_neg),
            z0,
            td,
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
    funcs: &FuncTable,
    expanding: &mut HashSet<String>,
    call_lineno: usize,
) -> Result<Expansion, ParseError> {
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

    let mut out = Expansion::default();
    // Model cards declared in this body: local name → per-instance mangled name.
    // Each instance gets its OWN card built from its own parameters, which is
    // what makes a `.subckt` usable as a PCell — LEVEL is only ever read from a
    // card, so without this an instance could not carry its own EO model.
    let mut local_models: HashMap<String, String> = HashMap::new();

    // Two passes over the body: cards first, so an element line may reference a
    // model declared below it (SPICE decks are order-independent for `.model`).
    for (lineno, body_line) in &def.body_lines {
        let trimmed = body_line.trim();
        if !trimmed.to_lowercase().starts_with(".model") {
            continue;
        }
        let substituted = substitute_params(trimmed, &inst_params, funcs, *lineno)?;
        if let Some(mut card) = parse_model(&substituted, *lineno)? {
            let local = card.name.to_lowercase();
            let mangled = format!("{inst_name}.{local}");
            local_models.insert(local, mangled.clone());
            card.name = mangled;
            out.models.push(card);
        }
    }

    for (lineno, body_line) in &def.body_lines {
        let lineno = *lineno;
        let trimmed = body_line.trim();
        if trimmed.is_empty() || trimmed.starts_with('*') {
            continue;
        }
        let lc = trimmed.to_lowercase();

        // Directives consumed by collect_defs or the card pass above — skip.
        if lc == ".end"
            || lc.starts_with(".ends")
            || lc.starts_with(".subckt")
            || lc.starts_with(".param")
            || lc.starts_with(".model")
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

        let substituted = substitute_params(trimmed, &inst_params, funcs, lineno)?;
        for el in parse_element_expanded(&substituted, lineno, funcs)? {
            let mut el = remap_element_nodes(el, &port_map, inst_name);

            // Point references at this instance's own copy of a local card.
            // Checked before the subckt lookup so a local model always wins.
            let mut used_local_model = false;
            if let Element::XOsdi {
                ref mut model_name, ..
            } = el
            {
                if let Some(mangled) = local_models.get(model_name.as_str()) {
                    *model_name = mangled.clone();
                    used_local_model = true;
                }
            }

            // Recurse if this element is a nested subckt instance.
            let is_subckt_inst = !used_local_model
                && matches!(&el, Element::XOsdi { model_name, .. }
                    if subckt_defs.contains_key(model_name));

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
                        funcs,
                        expanding,
                        lineno,
                    )?;
                    out.elements.extend(nested.elements);
                    out.models.extend(nested.models);
                }
            } else {
                out.elements.push(el);
            }
        }
    }

    expanding.remove(def_name);
    Ok(out)
}
