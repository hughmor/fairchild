//! Spectre-dialect front end.
//!
//! Foundry model libraries are written in Spectre (`.scs`), and a PDK's own text
//! is the only faithful source for it — a converter that resolves the model
//! authors' conditionals on their behalf produces plausible numbers nobody chose.
//! So fairchild reads the dialect.
//!
//! **This is a front end, not a second parser.** Every Spectre statement is
//! transliterated into the equivalent SPICE statement, and the existing passes do
//! the rest, so exactly one place still decides what a resistor is, what a
//! `.param` means, and how a subcircuit flattens. The two dialects share every
//! concept in play — parameters, functions, conditionals, models, subcircuits,
//! global nets — and differ only in spelling:
//!
//! | Spectre | SPICE |
//! |---|---|
//! | `parameters a=1 b=2*a` | `.param a=1 b={2*a}` |
//! | `R1 (in out) resistor r=1k` | `R1 in out 1k` |
//! | `V1 (in 0) vsource dc=1.8` | `V1 in 0 DC 1.8` |
//! | `X1 (a b) mycell w=1u` | `X1 a b mycell w=1u` |
//! | `dc1 dc`, `tr1 tran stop=1n` | `.op`, `.tran … 1n` |
//! | `vdd!` | a net, plus a `.global vdd!` |
//! | `include "f"` / `include "f" section=s` | `.include "f"` / `.lib "f" s` |
//!
//! Transliteration is **line-aligned**: a statement spanning input lines 10–14
//! emits its SPICE form on output line 10 and blanks on 11–14, so every error
//! message the SPICE passes produce still names the line the user wrote.

use crate::warn_user;
use crate::ParseError;
use std::collections::{HashMap, HashSet};

/// Parse a Spectre-dialect netlist.
///
/// Transliterates to SPICE and hands the result to the SPICE passes, so the
/// resulting [`crate::Netlist`] is indistinguishable from one parsed from a
/// `.sp` deck — which is the point: nothing downstream learns the dialect.
pub fn parse_spectre(input: &str) -> Result<crate::Netlist, ParseError> {
    let spice = to_spice(input)?;
    crate::parse_spice(&spice)
}

/// Does this text look like Spectre rather than SPICE?
///
/// Checked by content, not extension: a `.scs` suffix is a convention and a
/// model library that arrives with another name is still Spectre. `simulator
/// lang=` is the marker the dialect defines for exactly this purpose, and a
/// leading `//` comment is the other unambiguous tell — SPICE has no `//`.
pub fn looks_like_spectre(text: &str) -> bool {
    for line in text.lines().take(200) {
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        let lc = t.to_lowercase();
        if lc.starts_with("simulator lang=spectre") {
            return true;
        }
        if lc.starts_with("simulator lang=spice") {
            return false;
        }
        if t.starts_with("//") {
            return true;
        }
    }
    false
}

/// Which dialect the lexer is currently reading.
#[derive(PartialEq, Clone, Copy)]
enum Lang {
    Spectre,
    Spice,
}

/// One logical statement: its first input line, and its text with continuations
/// joined.
struct Stmt {
    /// 1-based line the statement starts on, for line-aligned output.
    lineno: usize,
    /// How many input lines it consumed, so the output can pad to match.
    span: usize,
    text: String,
    lang: Lang,
}

/// Transliterate Spectre text into SPICE text, line for line.
///
/// The output is fed to the SPICE passes verbatim, so anything this cannot
/// express becomes either a comment (when it carries no meaning for a solve) or
/// an error (when dropping it would change the answer).
pub fn to_spice(input: &str) -> Result<String, ParseError> {
    let stmts = lex(input);
    let mut out: Vec<String> = Vec::new();
    let mut globals: Vec<String> = Vec::new();

    let (hoisted, consumed) = hoist_parameters(&stmts);

    for (i, st) in stmts.iter().enumerate() {
        // Pad to the statement's own line so error messages keep their numbers.
        while out.len() + 1 < st.lineno {
            out.push(String::new());
        }
        let rendered = match st.lang {
            // A `simulator lang=spice` region is SPICE already.
            Lang::Spice => st.text.clone(),
            Lang::Spectre if consumed.contains(&i) => {
                format!("* hoisted to the .subckt header: {}", st.text)
            }
            Lang::Spectre => {
                statement(st, &mut globals, hoisted.get(&i).map_or("", String::as_str))?
            }
        };
        out.push(rendered);
        for _ in 1..st.span {
            out.push(String::new());
        }
    }

    // Spectre marks a global net with a trailing `!` on every use rather than
    // declaring it, so the declaration is synthesised from the uses — appended,
    // not prepended, because a prologue would shift every line number by one and
    // `.global` is order-free on the SPICE side.
    if !globals.is_empty() {
        globals.sort();
        globals.dedup();
        out.push(format!(".global {}", globals.join(" ")));
    }
    Ok(out.join("\n") + "\n")
}

/// Find each `subckt` block's `parameters` statements, so they can be hoisted onto
/// its header line.
///
/// Spectre declares a subcircuit's interface parameters *inside* the body; SPICE
/// declares them on the `.subckt` line. Returns the text to append to each
/// `subckt` statement (by index) and the indices of the `parameters` statements
/// that were consumed, which then render as comments so no line number moves.
///
/// A `parameters` statement anywhere in the block is hoisted, not only the first:
/// Spectre's are overridable wherever they appear, and the SPICE side resolves
/// header parameters in order, so their relative order is what carries the
/// dependencies between them.
fn hoist_parameters(stmts: &[Stmt]) -> (HashMap<usize, String>, HashSet<usize>) {
    let mut hoisted: HashMap<usize, String> = HashMap::new();
    let mut consumed: HashSet<usize> = HashSet::new();
    // Index of the innermost open `subckt`, plus the ones enclosing it.
    let mut open: Vec<usize> = Vec::new();

    for (i, st) in stmts.iter().enumerate() {
        if st.lang == Lang::Spice {
            continue;
        }
        let head = first_word(&st.text).to_lowercase();
        let is_subckt = head == "subckt"
            || (head == "inline" && second_word(&st.text).eq_ignore_ascii_case("subckt"));
        if is_subckt {
            open.push(i);
        } else if head == "ends" || head == "endsubckt" {
            open.pop();
        } else if (head == "parameters" || head == "parameter") && !open.is_empty() {
            let owner = *open.last().unwrap();
            let text = braced_assignments(&st.text);
            if text.is_empty() {
                continue;
            }
            let entry = hoisted.entry(owner).or_default();
            if !entry.is_empty() {
                entry.push(' ');
            }
            entry.push_str(&text);
            consumed.insert(i);
        }
    }
    (hoisted, consumed)
}

// ─── lexer ───────────────────────────────────────────────────────────────────

/// Split Spectre text into logical statements.
///
/// Handles the four lexical features that differ from SPICE: `//` comments to
/// end of line, a trailing `\` continuation, `simulator lang=` switching, and
/// brace-delimited blocks whose contents belong to the statement that opened
/// them. Line-leading `*` comments and `+` continuations are shared with SPICE
/// and handled the same way.
fn lex(input: &str) -> Vec<Stmt> {
    let mut stmts: Vec<Stmt> = Vec::new();
    let mut lang = Lang::Spectre;
    let mut pending: Option<Stmt> = None;
    // Depth of open `{` — while positive, lines join the block's statement.
    let mut depth = 0usize;
    // Set when the previous line ended in `\`, which makes this one a
    // continuation even without a leading `+`.
    let mut joined = false;

    for (i, raw) in input.lines().enumerate() {
        let lineno = i + 1;
        let (body, was_continued) = strip_comments(raw);
        let trimmed = body.trim();

        if trimmed.is_empty() && depth == 0 {
            // Blank or comment-only: flush nothing, but keep line alignment by
            // letting the next statement pad.
            continue;
        }

        // `simulator lang=` switches dialect for everything after it.
        let lc = trimmed.to_lowercase();
        if lc.starts_with("simulator lang=") {
            if let Some(p) = pending.take() {
                stmts.push(p);
            }
            lang = if lc.contains("spice") {
                Lang::Spice
            } else {
                Lang::Spectre
            };
            continue;
        }

        // A control-flow line owns its braces: `if (c) {`, `} else {`, `}`. It must
        // not be joined into a block statement, because the statements inside it
        // are ordinary statements that still need their own lines — and their own
        // line numbers. Only recognised outside a data block, so a `}` closing a
        // `statistics { … }` is not mistaken for an `.endif`.
        if depth == 0 && !joined && control_flow(trimmed).is_some() {
            if let Some(p) = pending.take() {
                stmts.push(p);
            }
            stmts.push(Stmt {
                lineno,
                span: 1,
                text: trimmed.to_string(),
                lang,
            });
            continue;
        }

        let opens = trimmed.matches('{').count();
        let closes = trimmed.matches('}').count();

        // A `+` line continues the previous statement, as in SPICE.
        let is_plus = trimmed.starts_with('+');
        let piece = if is_plus {
            trimmed[1..].trim()
        } else {
            trimmed
        };

        // A `+` line continues the statement before it even when that statement
        // has already been flushed — the same rule the SPICE lexer follows, and
        // the reason it keeps a "last real line" cursor.
        let continues = is_plus || joined || depth > 0;
        let target = match (pending.as_mut(), continues) {
            (Some(_), _) => pending.as_mut(),
            (None, true) => stmts.last_mut(),
            (None, false) => None,
        };
        match target {
            Some(p) if continues || p.text.is_empty() => {
                p.text.push(' ');
                p.text.push_str(piece);
                p.span = lineno - p.lineno + 1;
            }
            _ => {
                if let Some(p) = pending.take() {
                    stmts.push(p);
                }
                pending = Some(Stmt {
                    lineno,
                    span: 1,
                    text: piece.to_string(),
                    lang,
                });
            }
        }

        depth = depth + opens - closes.min(depth + opens);
        joined = was_continued;

        // The statement ends when nothing is holding it open.
        if depth == 0 && !was_continued {
            if let Some(p) = pending.take() {
                stmts.push(p);
            }
        }
    }
    if let Some(p) = pending.take() {
        stmts.push(p);
    }
    stmts
}

/// Which conditional a Spectre line is, if any, and its condition.
///
/// Spectre writes conditionals with braces — `if (cond) {`, `} else if (c) {`,
/// `} else {`, `}` — where SPICE writes directives. The mapping is exact, so the
/// front end translates rather than interprets: the SPICE side already resolves
/// `.if` at parse time, and inside a `.subckt` it resolves it per instance.
fn control_flow(line: &str) -> Option<(&'static str, String)> {
    let t = line.trim();
    let after_brace = t.strip_prefix('}').unwrap_or(t).trim();
    let closes_only = t.starts_with('}') && after_brace.is_empty();
    if closes_only {
        return Some((".endif", String::new()));
    }
    let lc = after_brace.to_lowercase();
    // `else if (…) {` and `elseif (…) {` are the same statement.
    for (prefix, kind) in [
        ("else if", ".elseif"),
        ("elseif", ".elseif"),
        ("if", ".if"),
        ("else", ".else"),
    ] {
        if !lc.starts_with(prefix) {
            continue;
        }
        // A word boundary: `iface` is not `if`, `elsewhere` is not `else`.
        let rest = after_brace[prefix.len()..].trim_start();
        if !after_brace[prefix.len()..].is_empty()
            && !rest.starts_with('(')
            && !rest.starts_with('{')
        {
            continue;
        }
        if kind == ".else" {
            return Some((kind, String::new()));
        }
        // The condition is everything up to the trailing `{`.
        let cond = rest.trim_end().trim_end_matches('{').trim();
        if cond.is_empty() {
            continue;
        }
        return Some((kind, cond.to_string()));
    }
    None
}

/// Strip comments from one raw line; returns the body and whether a trailing
/// `\` asks for the next line to be joined.
///
/// `//` runs to end of line. A `*` in the first column is a comment too — foundry
/// trees use it freely (`***…`, `*+ commented-out continuation`), mixed with `//`
/// in the same file. A `*` anywhere else is multiplication.
fn strip_comments(raw: &str) -> (String, bool) {
    let lead = raw.trim_start();
    if lead.starts_with('*') || lead.starts_with(';') {
        return (String::new(), false);
    }
    let body = match lead.find("//") {
        Some(i) => &lead[..i],
        None => lead,
    };
    let body = body.trim_end();
    match body.strip_suffix('\\') {
        Some(b) => (b.to_string(), true),
        None => (body.to_string(), false),
    }
}

// ─── statement transliteration ───────────────────────────────────────────────

/// Spectre primitives that map onto a SPICE element letter, with the parameter
/// that carries the value.
const PRIMITIVES: &[(&str, char, &str)] = &[
    ("resistor", 'R', "r"),
    ("capacitor", 'C', "c"),
    ("inductor", 'L', "l"),
    ("vsource", 'V', "dc"),
    ("isource", 'I', "dc"),
];

/// Spectre analysis kinds and their SPICE directive.
const ANALYSES: &[&str] = &["dc", "tran", "ac", "noise", "op"];

/// Statement kinds that carry no meaning for a solve. In Spectre these are
/// *named* — `opt1 options …`, `vds_check assert …` — so the keyword is the
/// second word, not the first.
const IGNORED_KINDS: &[&str] = &[
    "options",
    "option",
    "save",
    "statistics",
    "assert",
    "check",
    "checklimit",
    "montecarlo",
    "shell",
    "info",
    "sweep",
    "alter",
    "altergroup",
];

fn statement(st: &Stmt, globals: &mut Vec<String>, hoisted: &str) -> Result<String, ParseError> {
    let text = st.text.trim();
    let head_lc = first_word(text).to_lowercase();

    // Statements that lead with their keyword.
    match head_lc.as_str() {
        "parameters" | "parameter" => return Ok(parameters(text)),
        "include" | "ahdl_include" => return Ok(include(text, &head_lc)),
        "global" => {
            // Spectre also has an explicit `global` statement, alongside `!`.
            for tok in text.split_whitespace().skip(1) {
                globals.push(tok.to_string());
            }
            return Ok(format!("* {text}"));
        }
        "subckt" => return subckt_header(text, hoisted, st.lineno),
        // `inline subckt` is a subcircuit whose body Spectre inlines at each
        // instantiation. Flattening is what the SPICE side does with every
        // subcircuit anyway, so the distinction has no effect here.
        "inline" if second_word(text).eq_ignore_ascii_case("subckt") => {
            let after = text
                .split_once(char::is_whitespace)
                .map(|(_, r)| r)
                .unwrap_or("");
            return subckt_header(after.trim(), hoisted, st.lineno);
        }
        "ends" | "endsubckt" => {
            return Ok(format!(".ends {}", second_word(text))
                .trim_end()
                .to_string());
        }
        "model" => return model(text, st.lineno),
        // A function definition: `real f(real x) { return …; }`.
        "real" | "integer" => return function(text, st.lineno),
        _ => {}
    }

    // Conditionals. The lexer keeps each as its own statement so the statements
    // inside the block keep their own lines — and their own line numbers.
    if let Some((kind, cond)) = control_flow(text) {
        return Ok(if cond.is_empty() {
            kind.to_string()
        } else {
            // The condition keeps the parentheses it was written with — the SPICE
            // side accepts both spellings, and re-wrapping would nest them.
            format!("{kind} {cond}")
        });
    }

    // A keyword-led statement with no meaning for a solve: `statistics { … }` has
    // no instance name in front of it, so its keyword is the first word.
    if IGNORED_KINDS.contains(&head_lc.as_str()) {
        return Ok(skipped(st, text, &head_lc));
    }

    // Otherwise the form is `name kind …`, and the kind decides.
    let kind = second_word(text).to_lowercase();
    if IGNORED_KINDS.contains(&kind.as_str()) {
        return Ok(skipped(st, text, &kind));
    }
    instance_or_analysis(st, text, globals)
}

/// `subckt div (in out)` → `.subckt div in out`, with the body's `parameters`
/// appended.
///
/// Spectre declares a subcircuit's interface parameters in the body; SPICE puts
/// them on the header line. They are the same thing — overridable defaults — and
/// the SPICE side resolves them per instance, so hoisting is a move, not a
/// reinterpretation. `hoisted` is what the pre-pass found in this block.
fn subckt_header(text: &str, hoisted: &str, lineno: usize) -> Result<String, ParseError> {
    let name = second_word(text);
    if name.is_empty() {
        return Err(ParseError::Syntax {
            line: lineno,
            msg: "Spectre 'subckt' with no name".into(),
        });
    }
    // Slice by token, not by searching for the name: `subckt s (a b)` contains
    // the name `s` inside `subckt`, and finding that one cuts the line in half.
    let rest = text
        .splitn(3, char::is_whitespace)
        .nth(2)
        .unwrap_or("")
        .trim();
    let ports = rest.trim_start_matches('(').trim_end_matches(')').trim();
    let mut line = format!(".subckt {name} {ports}");
    if !hoisted.is_empty() {
        line.push(' ');
        line.push_str(hoisted);
    }
    Ok(line.trim_end().to_string())
}

/// `model nch bsim4 type=n version=4.5` → `.model nch bsim4 (type=n version=4.5)`.
fn model(text: &str, lineno: usize) -> Result<String, ParseError> {
    let name = second_word(text);
    let kind = text.split_whitespace().nth(2).unwrap_or("");
    if name.is_empty() || kind.is_empty() {
        return Err(ParseError::Syntax {
            line: lineno,
            msg: format!("cannot read '{text}' as a Spectre model: expected `model <name> <master> [params]`"),
        });
    }
    // A braced `model` body is binning — `model nch bsim4 { 1: lmin=… }` — where
    // geometry picks one card of several. Picking the wrong bin is a wrong answer
    // with nothing to read, and nothing here implements the selection.
    if text.contains('{') && !text.contains("={") {
        return Err(ParseError::Syntax {
            line: lineno,
            msg: format!(
                "model '{name}' is binned (a braced body of numbered sections).                  Bin selection by geometry is not implemented, and guessing a bin                  would be a wrong answer with nothing to read: split the bins into                  separate models, or instantiate the one you want"
            ),
        });
    }
    let params = braced_assignments(text);
    Ok(if params.is_empty() {
        format!(".model {name} {kind}")
    } else {
        format!(".model {name} {kind} ({params})")
    })
}

/// `real f(real a, real b) { return expr; }` → `.func f(a,b)={expr}`.
fn function(text: &str, lineno: usize) -> Result<String, ParseError> {
    let bad = |msg: String| ParseError::Syntax { line: lineno, msg };
    let open = text.find('(').ok_or_else(|| {
        bad(format!(
            "cannot read '{text}' as a Spectre function definition"
        ))
    })?;
    let close = text.find(')').ok_or_else(|| {
        bad(format!(
            "function definition on this line has no closing ')': '{text}'"
        ))
    })?;
    let name = text[..open]
        .split_whitespace()
        .nth(1)
        .unwrap_or("")
        .to_string();
    // Arguments are typed in Spectre and untyped in SPICE: `real x` → `x`.
    let args: Vec<String> = text[open + 1..close]
        .split(',')
        .filter(|a| !a.trim().is_empty())
        .map(|a| a.split_whitespace().last().unwrap_or("").to_string())
        .collect();
    let body_open = text[close..]
        .find('{')
        .map(|i| close + i)
        .ok_or_else(|| bad(format!("function '{name}' has no body")))?;
    let body_close = text
        .rfind('}')
        .ok_or_else(|| bad(format!("function '{name}' has no closing '}}'")))?;
    let body = text[body_open + 1..body_close].trim();
    // One `return` and nothing else. A function with local variables or control
    // flow has no `.func` equivalent, and half-translating it would silently drop
    // the rest of the body.
    let expr = body
        .strip_prefix("return")
        .map(|e| e.trim().trim_end_matches(';').trim())
        .filter(|e| !e.is_empty() && !e.contains(';'))
        .ok_or_else(|| {
            bad(format!(
                "function '{name}' is not a single `return <expression>;`, which is                  all a .func can hold: local variables and control flow inside a                  function are not translated"
            ))
        })?;
    Ok(format!(".func {name}({}) = {{{expr}}}", args.join(", ")))
}

/// The `k=v` assignments of a statement, rendered for a SPICE card.
fn braced_assignments(text: &str) -> String {
    assignments(text)
        .into_iter()
        .filter(|(_, v)| !v.is_empty())
        .map(|(k, v)| format!("{k}={}", braced_value(&v)))
        .collect::<Vec<_>>()
        .join(" ")
}

/// One Spectre value as a SPICE value: a plain number stays as it is, anything
/// else is braced so the SPICE side evaluates it as a parse-time expression
/// rather than reading it as a name.
///
/// One place decides this, because every statement kind needs the same answer and
/// a value that keeps its Spectre quotes reaches the expression parser as a
/// syntax error two layers down from the line that wrote it.
fn braced_value(v: &str) -> String {
    if crate::spice::parse_spice_value(v).is_ok() || v.starts_with('{') {
        v.to_string()
    } else {
        format!("{{{}}}", strip_quotes(v))
    }
}

/// A Spectre expression may be quoted; the braces the SPICE side needs replace
/// the quotes rather than nesting inside them.
fn strip_quotes(v: &str) -> &str {
    v.trim_matches('\'').trim_matches('"')
}

/// `parameters a=1 b=2*a` → `.param a=1 b={2*a}`.
///
/// A value that is not a plain number is wrapped in braces so the SPICE side
/// evaluates it as a parse-time expression — which is what Spectre means by it.
fn parameters(text: &str) -> String {
    format!(".param {}", braced_assignments(text))
        .trim_end()
        .to_string()
}

/// `include "f"` → `.include "f"`; `include "f" section=s` → `.lib "f" s`.
///
/// `ahdl_include` names Verilog-A source, which becomes `.va "f"` — the same
/// request in SPICE spelling, compiled on demand by the consumer. It used to
/// warn and drop the line, which made a foundry PDK a pile of manual compiles
/// before the deck would load at all.
fn include(text: &str, head: &str) -> String {
    let toks: Vec<&str> = text.split_whitespace().collect();
    let path = toks
        .get(1)
        .map(|t| t.trim_matches('"').trim_matches('\''))
        .unwrap_or("");
    if head == "ahdl_include" {
        return format!(".va \"{path}\"");
    }
    let section = assignments(text)
        .into_iter()
        .find(|(k, _)| k == "section")
        .map(|(_, v)| v);
    match section {
        Some(s) => format!(".lib \"{path}\" {s}"),
        None => format!(".include \"{path}\""),
    }
}

/// Statements with no bearing on a solve, kept as comments so the line survives.
fn skipped(st: &Stmt, text: &str, head: &str) -> String {
    // `options` carries `temp=`, which changes the answer, so it is not silent.
    if head == "options" || head == "option" {
        if let Some((_, t)) = assignments(text).into_iter().find(|(k, _)| k == "temp") {
            return format!(".temp {t}");
        }
    }
    warn_user!(
        "line {}: Spectre '{head}' statement is not interpreted and has been skipped",
        st.lineno
    );
    format!("* skipped: {text}")
}

/// Everything else is `name (nodes) kind [params]`, `name nodes kind [params]`,
/// or an analysis statement `name kind [params]`.
fn instance_or_analysis(
    st: &Stmt,
    text: &str,
    globals: &mut Vec<String>,
) -> Result<String, ParseError> {
    let name = first_word(text);
    let rest = text[name.len()..].trim();

    // Nodes may be parenthesised — the modern form — or bare, which foundry
    // wrappers use for the primitive instance inside a subcircuit.
    let (nodes, after) = if let Some(close) = rest.find(')') {
        if rest.starts_with('(') {
            (rest[1..close].to_string(), rest[close + 1..].trim())
        } else {
            (String::new(), rest)
        }
    } else {
        (String::new(), rest)
    };

    let after_toks: Vec<&str> = after.split_whitespace().collect();
    let kind = after_toks.first().copied().unwrap_or("");
    let kind_lc = kind.to_lowercase();

    // An analysis statement has no nodes and a recognised kind.
    if nodes.is_empty() && ANALYSES.contains(&kind_lc.as_str()) {
        return analysis(st, &kind_lc, after);
    }

    let (nodes, kind, params) = if nodes.is_empty() {
        // Bare-node form: the last token that is not `k=v` is the model name.
        let positional: Vec<&str> = after_toks
            .iter()
            .copied()
            .take_while(|t| !t.contains('='))
            .collect();
        if positional.len() < 2 {
            return Err(ParseError::Syntax {
                line: st.lineno,
                msg: format!(
                    "cannot read '{text}' as a Spectre instance: expected \
                     `name (nodes) model [params]` or `name nodes model [params]`"
                ),
            });
        }
        let (model, nodes) = positional.split_last().unwrap();
        let params: Vec<&str> = after_toks[positional.len()..].to_vec();
        (nodes.join(" "), model.to_string(), params.join(" "))
    } else {
        (
            nodes,
            kind.to_string(),
            after[kind.len()..].trim().to_string(),
        )
    };

    for n in nodes.split_whitespace() {
        if n.ends_with('!') {
            globals.push(n.to_string());
        }
    }

    let kind_lc = kind.to_lowercase();
    if let Some((_, letter, value_param)) = PRIMITIVES
        .iter()
        .find(|(spectre, _, _)| *spectre == kind_lc.as_str())
    {
        return Ok(primitive(&name, *letter, value_param, &nodes, &params));
    }

    // Anything else is a subcircuit or model instance: SPICE's X line, which
    // also carries `k=v` parameters unchanged.
    Ok(
        format!("{} {} {} {}", spice_name('X', &name), nodes, kind, params)
            .trim_end()
            .to_string(),
    )
}

/// `R1 (a b) resistor r=1k` → `R1 a b 1k`, with an expression braced.
fn primitive(name: &str, letter: char, value_param: &str, nodes: &str, params: &str) -> String {
    let pairs = assignments(params);
    let value = pairs
        .iter()
        .find(|(k, _)| k == value_param)
        .map(|(_, v)| v.clone())
        .unwrap_or_else(|| "0".to_string());
    let value = braced_value(&value);
    // A source's value keyword is spelled out in SPICE; everything else is
    // positional.
    let head = format!("{} {nodes}", spice_name(letter, name));
    let mut line = match letter {
        'V' | 'I' => format!("{head} DC {value}"),
        _ => format!("{head} {value}"),
    };
    // Extra parameters SPICE understands on the same element keep their form.
    for (k, v) in pairs {
        if k == value_param {
            continue;
        }
        if matches!(k.as_str(), "ac" | "mag" | "phase") {
            continue; // AC spec is handled below, not as a k=v
        }
        line.push_str(&format!(" {k}={}", braced_value(&v)));
    }
    line
}

/// `dc1 dc dev=Vg param=dc start=0 stop=1 step=0.1` → `.dc Vg 0 1 0.1`, and so on.
fn analysis(st: &Stmt, kind: &str, text: &str) -> Result<String, ParseError> {
    let p = assignments(text);
    let get = |k: &str| p.iter().find(|(n, _)| n == k).map(|(_, v)| v.clone());
    match kind {
        "op" => Ok(".op".to_string()),
        "dc" => {
            // A `dc` with no swept device is an operating point.
            match (get("dev"), get("start"), get("stop")) {
                (Some(dev), Some(start), Some(stop)) => {
                    let step = get("step").unwrap_or_else(|| "0.1".into());
                    Ok(format!(".dc {dev} {start} {stop} {step}"))
                }
                _ => Ok(".op".to_string()),
            }
        }
        "tran" => {
            let stop = get("stop").ok_or_else(|| ParseError::Syntax {
                line: st.lineno,
                msg: "Spectre tran analysis needs stop=".into(),
            })?;
            // Spectre's step is advisory (`step=`/`maxstep=`); SPICE needs one.
            let step = get("step")
                .or_else(|| get("maxstep"))
                .unwrap_or_else(|| format!("{{{stop}/1000}}"));
            Ok(format!(".tran {step} {stop}"))
        }
        "ac" => {
            let start = get("start").unwrap_or_else(|| "1".into());
            let stop = get("stop").unwrap_or_else(|| "1e9".into());
            let dec = get("dec").or_else(|| get("log"));
            let lin = get("lin");
            match (dec, lin) {
                (Some(n), _) => Ok(format!(".ac dec {n} {start} {stop}")),
                (None, Some(n)) => Ok(format!(".ac lin {n} {start} {stop}")),
                _ => Ok(format!(".ac dec 10 {start} {stop}")),
            }
        }
        other => Err(ParseError::Syntax {
            line: st.lineno,
            msg: format!("Spectre analysis '{other}' is not supported"),
        }),
    }
}

// ─── small helpers ───────────────────────────────────────────────────────────

fn first_word(s: &str) -> String {
    s.split_whitespace().next().unwrap_or("").to_string()
}

/// The second whitespace token — a Spectre statement's *kind*.
fn second_word(s: &str) -> String {
    s.split_whitespace().nth(1).unwrap_or("").to_string()
}

/// A SPICE element's letter comes from its type, while a Spectre instance carries
/// whatever name the deck gave it. Prepend the letter only when the name does not
/// already start with it, so `R1 (a b) resistor` stays `R1` instead of becoming
/// `RR1` — the name is what a probe and every error message will say, and a deck's
/// own names are the ones its author will look for.
fn spice_name(letter: char, name: &str) -> String {
    let first = name.chars().next().unwrap_or(' ').to_ascii_lowercase();
    if first == letter.to_ascii_lowercase() {
        name.to_string()
    } else {
        format!("{letter}{name}")
    }
}

/// Split `k=v` assignments out of a statement, keeping a braced, quoted or
/// parenthesised value whole and gluing `k = v` written with spaces.
///
/// Shares its rule with the SPICE side: a value may contain spaces only inside
/// a bracket or quote, because `a = 1 + 2 b = 3` has no unambiguous reading.
fn assignments(text: &str) -> Vec<(String, String)> {
    let mut toks: Vec<String> = Vec::new();
    let mut cur = String::new();
    let mut span: Option<char> = None;
    let mut paren = 0usize;
    for c in text.chars() {
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
                '(' => {
                    paren += 1;
                    cur.push(c);
                }
                ')' => {
                    paren = paren.saturating_sub(1);
                    cur.push(c);
                }
                ' ' | '\t' if paren == 0 => {
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

    // Re-glue `name = value`, `name= value`, `name =value`.
    let mut glued: Vec<String> = Vec::new();
    let mut i = 0;
    while i < toks.len() {
        let mut t = toks[i].clone();
        if t == "=" {
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
        .into_iter()
        .filter_map(|t| {
            t.split_once('=')
                .map(|(k, v)| (k.trim().to_lowercase(), v.trim().to_string()))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Analysis, Element};

    /// Every fixture here is hand-written to mimic a construct seen in a foundry
    /// library. **No foundry text appears in this repository** — the PDK is
    /// confidential and this repo is public, so the shapes are reproduced and the
    /// content is invented.
    #[test]
    fn an_rc_deck_parses_like_its_spice_twin() {
        let netlist = parse_spectre(
            "\
// a first-order lowpass
simulator lang=spectre
parameters vdd=1.8 rload=1k
V1 (in 0) vsource dc=vdd
R1 (in out) resistor r=rload
C1 (out 0) capacitor c=1n
tr1 tran stop=1u step=1n
",
        )
        .unwrap();
        assert_eq!(netlist.elements.len(), 3, "{:?}", netlist.elements);
        let r = netlist
            .elements
            .iter()
            .find_map(|e| match e {
                Element::Resistor { resistance, .. } => Some(*resistance),
                _ => None,
            })
            .expect("no resistor");
        assert!((r - 1000.0).abs() < 1e-9, "r={r}");
        assert!(
            matches!(netlist.analyses.first(), Some(Analysis::Tran { .. })),
            "{:?}",
            netlist.analyses
        );
    }

    #[test]
    fn a_source_value_may_be_an_expression_over_parameters() {
        let netlist = parse_spectre(
            "simulator lang=spectre\nparameters vdd=1.8\nV1 (in 0) vsource dc=vdd/2\nR1 (in 0) resistor r=1k\ndc1 dc\n",
        )
        .unwrap();
        let dc = netlist
            .elements
            .iter()
            .find_map(|e| match e {
                Element::VoltageSource { waveform, .. } => Some(format!("{waveform:?}")),
                _ => None,
            })
            .unwrap();
        assert!(dc.contains("0.9"), "expected vdd/2 = 0.9, got {dc}");
    }

    #[test]
    fn line_numbers_survive_transliteration() {
        // The whole point of padding: an error names the line the user wrote, not
        // the line some intermediate text happened to land on. `nope` is
        // undefined, and it sits on input line 5.
        let err = parse_spectre(
            "// header\n\
             simulator lang=spectre\n\
             parameters good=1\n\
             R1 (a 0) resistor r=1k\n\
             R2 (a 0) resistor r=nope*2\n",
        )
        .unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("line 5"), "{msg}");
        assert!(msg.contains("nope"), "{msg}");
    }

    #[test]
    fn continuations_join_a_statement_both_ways() {
        // A trailing `\` (Spectre) and a leading `+` (shared with SPICE) both
        // continue a statement, and a foundry file uses them in the same breath.
        let netlist = parse_spectre(
            "simulator lang=spectre\n\
             parameters a=2 \\\n\
             b=3\n\
             R1 (n1 0) resistor\n\
             + r=1k\n\
             V1 (n1 0) vsource dc=1\n\
             dc1 dc\n",
        )
        .unwrap();
        assert_eq!(netlist.elements.len(), 2, "{:?}", netlist.elements);
        // Both parameters landed, so a later expression can use them.
        let n2 = parse_spectre(
            "simulator lang=spectre\nparameters a=2 \\\nb=3\nR1 (n1 0) resistor r=a*b*1k\nV1 (n1 0) vsource dc=1\ndc1 dc\n",
        )
        .unwrap();
        let r = n2
            .elements
            .iter()
            .find_map(|e| match e {
                Element::Resistor { resistance, .. } => Some(*resistance),
                _ => None,
            })
            .unwrap();
        assert!((r - 6000.0).abs() < 1e-6, "r={r}");
    }

    #[test]
    fn comment_forms_are_both_honoured() {
        // `//` to end of line, and a line-leading `*` — a foundry library mixes
        // them, including `*`-commented continuation lines.
        let netlist = parse_spectre(
            "simulator lang=spectre\n\
             *** a banner comment\n\
             parameters a=1   // trailing comment\n\
             **pro\n\
             R1 (n 0) resistor r=1k\n\
             *+ r=9k\n\
             V1 (n 0) vsource dc=1\n\
             dc1 dc\n",
        )
        .unwrap();
        assert_eq!(netlist.elements.len(), 2);
        let r = netlist
            .elements
            .iter()
            .find_map(|e| match e {
                Element::Resistor { resistance, .. } => Some(*resistance),
                _ => None,
            })
            .unwrap();
        assert!(
            (r - 1000.0).abs() < 1e-9,
            "the commented-out value won: {r}"
        );
    }

    #[test]
    fn a_bang_suffixed_net_becomes_global() {
        // Spectre declares a global net by using `!`, nowhere else. Two instances
        // in different scopes must land on the same node.
        let spice = to_spice(
            "simulator lang=spectre\nR1 (vdd! out) resistor r=1k\nV1 (vdd! 0) vsource dc=1.8\n",
        )
        .unwrap();
        assert!(
            spice
                .lines()
                .any(|l| l.starts_with(".global") && l.contains("vdd!")),
            "no .global synthesised:\n{spice}"
        );
        // And exactly one declaration for two uses.
        assert_eq!(spice.matches(".global").count(), 1, "{spice}");
    }

    #[test]
    fn a_spice_region_passes_through_untouched() {
        // Mixed-language decks are normal: the analysis in Spectre, the output
        // controls in SPICE.
        let netlist = parse_spectre(
            "simulator lang=spectre\n\
             R1 (n 0) resistor r=1k\n\
             V1 (n 0) vsource dc=1\n\
             simulator lang=spice\n\
             .op\n\
             .options reltol=1e-5\n",
        )
        .unwrap();
        assert_eq!(netlist.elements.len(), 2);
        assert!(matches!(netlist.analyses.first(), Some(Analysis::Op)));
        assert!(
            netlist.options.iter().any(|(k, _)| k == "reltol"),
            "{:?}",
            netlist.options
        );
    }

    #[test]
    fn include_and_section_map_to_their_spice_forms() {
        let spice = to_spice(
            "simulator lang=spectre\ninclude \"models.scs\"\ninclude \"corners.scs\" section=tt\n",
        )
        .unwrap();
        assert!(spice.contains(".include \"models.scs\""), "{spice}");
        assert!(spice.contains(".lib \"corners.scs\" tt"), "{spice}");
    }

    /// `ahdl_include` is the whole reason a foundry PDK loads at all: it used
    /// to warn and drop the line, so every model it named became an unknown
    /// model unless the user compiled it by hand first.
    #[test]
    fn ahdl_include_becomes_a_va_directive() {
        let spice = to_spice(
            "simulator lang=spectre\nahdl_include \"bsim.va\"\nahdl_include \"rdiff.va\"\n",
        )
        .unwrap();
        assert!(spice.contains(".va \"bsim.va\""), "{spice}");
        // Order survives: a PDK compiles its sources in the order it lists them.
        assert!(
            spice.find(".va \"bsim.va\"") < spice.find(".va \"rdiff.va\""),
            "{spice}"
        );
    }

    #[test]
    fn options_temp_is_not_dropped() {
        // Temperature silently changing across a translation is a wrong answer,
        // not a cosmetic loss — everything else in `options` is skipped.
        let spice = to_spice("simulator lang=spectre\nopt1 options temp=85 gmin=1e-13\n").unwrap();
        assert!(spice.contains(".temp 85"), "{spice}");
    }

    #[test]
    fn an_unreadable_instance_is_an_error_not_a_guess() {
        let err = parse_spectre("simulator lang=spectre\nfoo bar\n").unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("line 2"), "{msg}");
    }

    /// The shape a foundry wrapper has: an `inline subckt` whose `parameters` sit
    /// in the body, one of them an expression over another, instantiated with an
    /// override. The parameters must reach the `.subckt` header, because that is
    /// the only place SPICE lets a call override them.
    #[test]
    fn a_wrapper_hoists_its_parameters_to_the_header() {
        let net = parse_spectre(
            "simulator lang=spectre\n\
             inline subckt res_wrap (a b)\n\
             parameters w=1u l=2u sheet=100\n\
             parameters rval='sheet*l/w'\n\
             R1 (a b) resistor r=rval\n\
             ends res_wrap\n\
             V1 (in 0) vsource dc=1\n\
             X1 (in 0) res_wrap w=2u\n\
             dc1 dc\n",
        )
        .unwrap();
        let r = sole_resistance(&net);
        // sheet*l/w = 100 * 2u / 2u; the default w=1u would give 200.
        assert!((r - 100.0).abs() < 1e-9, "got {r}, want 100");
    }

    /// A `parameters` line that was hoisted must not also survive where it stood,
    /// and the lines after it must keep their numbers.
    #[test]
    fn a_hoist_leaves_the_line_numbers_alone() {
        let spice = to_spice(
            "simulator lang=spectre\n\
             subckt s (a b)\n\
             parameters w=1u\n\
             R1 (a b) resistor r=w\n\
             ends s\n",
        )
        .unwrap();
        let lines: Vec<&str> = spice.lines().collect();
        assert!(lines[1].starts_with(".subckt s a b w=1u"), "{:?}", lines[1]);
        assert!(lines[2].starts_with('*'), "{:?}", lines[2]);
        assert!(lines[3].starts_with("R1 a b"), "{:?}", lines[3]);
    }

    /// Spectre writes a conditional with braces, SPICE with directives. The block's
    /// statements must keep their own lines, which is why the lexer does not join
    /// a control-flow line into the block it opens.
    #[test]
    fn a_conditional_block_becomes_directives() {
        let spice = to_spice(
            "simulator lang=spectre\n\
             parameters corner=2\n\
             if (corner == 1) {\n\
             R1 (a 0) resistor r=1k\n\
             } else if (corner == 2) {\n\
             R1 (a 0) resistor r=2k\n\
             } else {\n\
             R1 (a 0) resistor r=3k\n\
             }\n",
        )
        .unwrap();
        let lines: Vec<&str> = spice.lines().collect();
        assert_eq!(lines[2], ".if (corner == 1)", "{lines:?}");
        assert_eq!(lines[4], ".elseif (corner == 2)", "{lines:?}");
        assert_eq!(lines[6], ".else", "{lines:?}");
        assert_eq!(lines[8], ".endif", "{lines:?}");
    }

    /// And it has to give the right answer through both layers: the SPICE side
    /// resolves the conditional, so only one branch survives.
    #[test]
    fn a_conditional_selects_the_branch_it_should() {
        let net = parse_spectre(
            "simulator lang=spectre\n\
             parameters corner=2\n\
             V1 (in 0) vsource dc=1\n\
             if (corner == 1) {\n\
             R1 (in 0) resistor r=1k\n\
             } else {\n\
             R1 (in 0) resistor r=2k\n\
             }\n\
             dc1 dc\n",
        )
        .unwrap();
        let r = sole_resistance(&net);
        assert!((r - 2000.0).abs() < 1e-9, "got {r}");
    }

    /// A conditional inside a subcircuit body reaches the SPICE side as a `.if`
    /// there, which resolves it per instance.
    #[test]
    fn a_conditional_inside_a_subckt_selects_per_instance() {
        let net = parse_spectre(
            "simulator lang=spectre\n\
             subckt rsel (a b)\n\
             parameters mode=0\n\
             if (mode == 1) {\n\
             R1 (a b) resistor r=2k\n\
             } else {\n\
             R1 (a b) resistor r=1k\n\
             }\n\
             ends rsel\n\
             V1 (in 0) vsource dc=1\n\
             Xa (in 0) rsel mode=1\n\
             Xb (in 0) rsel\n\
             dc1 dc\n",
        )
        .unwrap();
        let mut got: Vec<f64> = net
            .elements
            .iter()
            .filter_map(|el| match el {
                Element::Resistor { resistance, .. } => Some(*resistance),
                _ => None,
            })
            .collect();
        got.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(got.len(), 2, "{got:?}");
        assert!((got[0] - 1000.0).abs() < 1e-9, "{got:?}");
        assert!((got[1] - 2000.0).abs() < 1e-9, "{got:?}");
    }

    #[test]
    fn a_function_definition_becomes_a_func() {
        let net = parse_spectre(
            "simulator lang=spectre\n\
             real half(real x) { return x/2; }\n\
             parameters r0=4k\n\
             V1 (in 0) vsource dc=1\n\
             R1 (in 0) resistor r='half(r0)'\n\
             dc1 dc\n",
        )
        .unwrap();
        let r = sole_resistance(&net);
        assert!((r - 2000.0).abs() < 1e-9, "got {r}");
    }

    /// A function body that is not a single `return` has no `.func` equivalent,
    /// and translating half of it would drop the rest in silence.
    #[test]
    fn a_multi_statement_function_is_refused() {
        let err =
            to_spice("simulator lang=spectre\nreal f(real x) { real y; y = x*2; return y; }\n")
                .unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("return"), "{msg}");
        assert!(msg.contains("line 2"), "{msg}");
    }

    #[test]
    fn a_model_statement_becomes_a_model_card() {
        let spice =
            to_spice("simulator lang=spectre\nmodel my_d diode is=1e-16 rs='2*1'\n").unwrap();
        assert_eq!(
            spice.lines().nth(1).unwrap(),
            ".model my_d diode (is=1e-16 rs={2*1})"
        );
    }

    /// Binning picks one card of several by geometry. Nothing here implements the
    /// selection, and guessing a bin is a wrong answer with nothing to read.
    #[test]
    fn a_binned_model_is_refused() {
        let err =
            to_spice("simulator lang=spectre\nmodel nch bsim4 {\n1: lmin=1n lmax=2n vth0=0.4\n}\n")
                .unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("binned"), "{msg}");
        assert!(msg.contains("nch"), "{msg}");
    }

    /// A `statistics` block spans lines and holds statements of its own. Skipping
    /// the keyword line only would leak its contents into the netlist.
    #[test]
    fn a_statistics_block_is_skipped_whole() {
        let spice = to_spice(
            "simulator lang=spectre\n\
             statistics {\n\
             process {\n\
             vary rsh dist=gauss std=0.1\n\
             }\n\
             }\n\
             R1 (a 0) resistor r=1k\n",
        )
        .unwrap();
        for (i, line) in spice.lines().enumerate() {
            if i == 6 {
                assert!(line.starts_with("R1 a 0"), "line {i}: {line:?}");
            } else {
                assert!(
                    line.is_empty() || line.starts_with('*'),
                    "line {i} leaked: {line:?}"
                );
            }
        }
    }

    /// The deck's own names are the ones its author will look for in a probe or an
    /// error message, so an instance already named `R1` must not become `RR1`.
    #[test]
    fn an_instance_keeps_the_name_the_deck_gave_it() {
        let net = parse_spectre(
            "simulator lang=spectre\n\
             V1 (in 0) vsource dc=1\n\
             R1 (in 0) resistor r=1k\n\
             rload (in 0) resistor r=2k\n\
             dc1 dc\n",
        )
        .unwrap();
        let mut names: Vec<&str> = net
            .elements
            .iter()
            .filter_map(|el| match el {
                Element::Resistor { name, .. } => Some(name.as_str()),
                _ => None,
            })
            .collect();
        names.sort();
        assert_eq!(names, vec!["r1", "rload"], "{names:?}");
    }

    /// The one resistance in a netlist.
    fn sole_resistance(net: &crate::Netlist) -> f64 {
        let mut found = net.elements.iter().filter_map(|el| match el {
            Element::Resistor { resistance, .. } => Some(*resistance),
            _ => None,
        });
        let r = found.next().expect("no resistor in the netlist");
        assert!(found.next().is_none(), "expected exactly one resistor");
        r
    }

    #[test]
    fn assignments_keeps_a_braced_value_whole() {
        let got = assignments("x r={a + b} c = 1n d=2");
        assert_eq!(
            got,
            vec![
                ("r".into(), "{a + b}".into()),
                ("c".into(), "1n".into()),
                ("d".into(), "2".into())
            ]
        );
    }
}
