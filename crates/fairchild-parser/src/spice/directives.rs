use super::common::{canon_node, parse_value};
use crate::expr::Expr;
use crate::{
    AcVariation, Analysis, DcSweepSpec, MeasAnalysis, MeasKind, MeasOp, Measurement, ParseError,
};

/// `.backanno` (LTspice schematic back-annotation) is the one directive still
/// ignored without a word: it selects nothing, changes nothing, and there is no
/// fairchild mechanism to point a warning at.
pub(super) fn is_silent_directive(lc: &str) -> bool {
    lc.starts_with(".backanno")
}

/// Directives that select *what to report* rather than what the circuit is or
/// what to run.  Returns the directive token, for the warning that names it.
///
/// Output selection belongs to the frontend here — `--probe` on the CLI, numpy
/// indexing in Python — and every signal is available either way, so a deck's
/// version has nowhere to land.  They are ignored, but not in silence: a deck
/// whose `.print tran V(out)` you expect to narrow the output gets every node
/// instead, and that is worth one line of stderr.  `.save` and `.width` used to
/// be hard errors, which is the same class of directive failing a different way
/// for no reason anyone chose.
pub(super) fn select_directive(lc: &str) -> Option<&'static str> {
    const SELECT: &[&str] = &[".print", ".plot", ".probe", ".save", ".width"];
    SELECT.iter().copied().find(|d| lc.starts_with(d))
}

/// Parse a `.measure` / `.meas` directive into a `Measurement`.
///
/// Supported forms (case-insensitive keywords):
///   .meas tran NAME FIND <expr> AT=<t>
///   .meas tran NAME FIND <expr> WHEN <cond> [CROSS=<n>]
///   .meas tran NAME MAX|MIN|AVG|RMS|PP|INTEG <expr> [FROM=<t1>] [TO=<t2>]
///   .meas tran NAME DERIV <expr> AT=<t>
///   .meas tran NAME TRIG <cond1> [VAL=<v1>] [CROSS=<n>] TARG <cond2> [VAL=<v2>] [CROSS=<n>]
///
/// `<expr>` and `<cond>` may be a single tagged reference (e.g. `V(out)`) or
/// any expression accepted by `Expr::parse`.  Comparisons in WHEN/TRIG/TARG
/// either come from the expression itself (`V(out) > 0.5`) or use the
/// `VAL=<v>` keyword which adds an implicit `expr - val` zero-cross.
pub(super) fn parse_measure(line: &str, lineno: usize) -> Result<Measurement, ParseError> {
    let toks: Vec<&str> = line.split_whitespace().collect();
    if toks.len() < 4 {
        return Err(ParseError::FieldCount {
            expected: "≥4 (.meas analysis name op …)",
            got: toks.len(),
            line: lineno,
        });
    }
    let analysis = match toks[1].to_lowercase().as_str() {
        "tran" => MeasAnalysis::Tran,
        "dc" => MeasAnalysis::Dc,
        "ac" => MeasAnalysis::Ac,
        other => {
            return Err(ParseError::Syntax {
                line: lineno,
                msg: format!(".measure analysis '{other}' unsupported (use tran|dc|ac)"),
            })
        }
    };
    let name = toks[2].to_string();
    let op_word = toks[3].to_lowercase();

    // Helper: parse an expression spelled across `parts` (whitespace-joined),
    // also recognising the SPICE shorthand of writing a comparison as
    // `<lhs> <relop> <rhs>` where `<relop>` is a separate token (`>`, `<`,
    // `=`, etc.).
    let parse_expr = |s: &str| -> Result<Expr, ParseError> {
        Expr::parse(s).map_err(|e| ParseError::Syntax {
            line: lineno,
            msg: format!(".meas expr '{s}': {e}"),
        })
    };

    // Helper: split tokens[start..] into (expr_tokens, keyword_pairs).
    // Keyword pairs are recognised by `KEY=VALUE` or `KEY VALUE` for AT/FROM/TO/VAL/CROSS.
    let kw_keys = ["at", "from", "to", "val", "cross"];
    let split_kw = |start: usize,
                    end: usize|
     -> (Vec<&str>, std::collections::HashMap<String, String>) {
        let mut expr_toks: Vec<&str> = Vec::new();
        let mut kws: std::collections::HashMap<String, String> = std::collections::HashMap::new();
        let mut i = start;
        while i < end {
            let t = toks[i];
            if let Some((k, v)) = t.split_once('=') {
                if kw_keys.iter().any(|kk| *kk == k.to_lowercase()) {
                    kws.insert(k.to_lowercase(), v.to_string());
                    i += 1;
                    continue;
                }
            }
            // Two-token form `KEY VALUE`?
            if kw_keys.iter().any(|kk| *kk == t.to_lowercase()) && i + 1 < end {
                kws.insert(t.to_lowercase(), toks[i + 1].to_string());
                i += 2;
                continue;
            }
            expr_toks.push(t);
            i += 1;
        }
        (expr_toks, kws)
    };

    match op_word.as_str() {
        "find" => {
            // Walk forward collecting expr tokens until WHEN or AT=.
            let mut split = toks.len();
            for (i, tok) in toks.iter().enumerate().skip(4) {
                let lc = tok.to_lowercase();
                if lc == "when" || lc == "at" || lc.starts_with("at=") {
                    split = i;
                    break;
                }
            }
            let expr_str = toks[4..split].join(" ");
            let expr = parse_expr(&expr_str)?;
            if split < toks.len() {
                let rest_lc = toks[split].to_lowercase();
                if rest_lc == "when" {
                    let cond_str = toks[split + 1..].join(" ");
                    let cond = parse_expr(&cond_str)?;
                    Ok(Measurement {
                        name,
                        analysis,
                        kind: MeasKind::FindWhen {
                            expr,
                            cond,
                            cross: 1,
                        },
                    })
                } else {
                    let (_e, kws) = split_kw(split, toks.len());
                    let at = kws.get("at").ok_or_else(|| ParseError::Syntax {
                        line: lineno,
                        msg: ".meas FIND requires AT=<time> or WHEN <cond>".into(),
                    })?;
                    let at_v = parse_value(at, lineno)?;
                    Ok(Measurement {
                        name,
                        analysis,
                        kind: MeasKind::FindAt { expr, at: at_v },
                    })
                }
            } else {
                Err(ParseError::Syntax {
                    line: lineno,
                    msg: ".meas FIND requires AT=<time> or WHEN <cond>".into(),
                })
            }
        }

        "deriv" => {
            // .meas tran NAME DERIV <expr> AT=<t>
            // Find AT token.
            let mut split = toks.len();
            for (i, tok) in toks.iter().enumerate().skip(4) {
                let lc = tok.to_lowercase();
                if lc == "at" || lc.starts_with("at=") {
                    split = i;
                    break;
                }
            }
            let expr_str = toks[4..split].join(" ");
            let expr = parse_expr(&expr_str)?;
            let (_e, kws) = split_kw(split, toks.len());
            let at = kws.get("at").ok_or_else(|| ParseError::Syntax {
                line: lineno,
                msg: ".meas DERIV requires AT=<time>".into(),
            })?;
            let at_v = parse_value(at, lineno)?;
            Ok(Measurement {
                name,
                analysis,
                kind: MeasKind::DerivAt { expr, at: at_v },
            })
        }

        "trig" => {
            // .meas tran NAME TRIG <expr1> [VAL=<v1>] [CROSS=<n>] TARG <expr2> [VAL=<v2>] [CROSS=<n>]
            let mut targ_idx = toks.len();
            for (i, tok) in toks.iter().enumerate().skip(4) {
                if tok.to_lowercase() == "targ" {
                    targ_idx = i;
                    break;
                }
            }
            if targ_idx == toks.len() {
                return Err(ParseError::Syntax {
                    line: lineno,
                    msg: ".meas TRIG requires a TARG clause".into(),
                });
            }
            let trig_part = &toks[4..targ_idx];
            let targ_part = &toks[targ_idx + 1..];

            let parse_part = |part: &[&str]| -> Result<(Expr, f64, usize), ParseError> {
                // Collect expr-tokens (until first VAL/CROSS keyword), then kw pairs.
                let mut split = part.len();
                for (i, t) in part.iter().enumerate() {
                    let lc = t.to_lowercase();
                    if lc == "val"
                        || lc.starts_with("val=")
                        || lc == "cross"
                        || lc.starts_with("cross=")
                    {
                        split = i;
                        break;
                    }
                }
                let expr = parse_expr(&part[..split].join(" "))?;
                let mut val = 0.0f64;
                let mut cross = 1usize;
                let mut i = split;
                while i < part.len() {
                    let t = part[i];
                    let (k, v) = if let Some((k, v)) = t.split_once('=') {
                        i += 1;
                        (k.to_lowercase(), v.to_string())
                    } else if i + 1 < part.len() {
                        let pair = (t.to_lowercase(), part[i + 1].to_string());
                        i += 2;
                        pair
                    } else {
                        i += 1;
                        continue;
                    };
                    match k.as_str() {
                        "val" => val = parse_value(&v, lineno)?,
                        "cross" => cross = v.parse().unwrap_or(1),
                        _ => {}
                    }
                }
                Ok((expr, val, cross))
            };

            let (trig_expr, trig_val, trig_cross) = parse_part(trig_part)?;
            let (targ_expr, targ_val, targ_cross) = parse_part(targ_part)?;
            Ok(Measurement {
                name,
                analysis,
                kind: MeasKind::TrigTarg {
                    trig_expr,
                    trig_val,
                    trig_cross,
                    targ_expr,
                    targ_val,
                    targ_cross,
                },
            })
        }

        "max" | "min" | "avg" | "rms" | "pp" | "integ" => {
            let op = match op_word.as_str() {
                "max" => MeasOp::Max,
                "min" => MeasOp::Min,
                "avg" => MeasOp::Avg,
                "rms" => MeasOp::Rms,
                "pp" => MeasOp::Pp,
                "integ" => MeasOp::Integ,
                _ => unreachable!(),
            };
            // FROM/TO are keywords; everything else is expression tokens.
            let mut split = toks.len();
            for (i, tok) in toks.iter().enumerate().skip(4) {
                let lc = tok.to_lowercase();
                if lc == "from" || lc == "to" || lc.starts_with("from=") || lc.starts_with("to=") {
                    split = i;
                    break;
                }
            }
            let expr = parse_expr(&toks[4..split].join(" "))?;
            let (_e, kws) = split_kw(split, toks.len());
            let from = kws
                .get("from")
                .map(|s| parse_value(s, lineno))
                .transpose()?;
            let to = kws.get("to").map(|s| parse_value(s, lineno)).transpose()?;
            Ok(Measurement {
                name,
                analysis,
                kind: MeasKind::Aggregate { op, expr, from, to },
            })
        }

        other => Err(ParseError::Syntax {
            line: lineno,
            msg: format!(
                ".meas op '{other}' unsupported (use FIND/MAX/MIN/AVG/RMS/PP/INTEG/DERIV/TRIG)"
            ),
        }),
    }
}

/// Parse `.ic V(n1)=val V(n2)=val ...` or `.nodeset V(n)=val ...`.
///
/// Returns a `Vec<(node, value)>` (node name lowercased, "gnd" canonicalised).
/// Tokens that don't match the `V(<name>)=<value>` shape are silently ignored
/// (they're typically the leading `.ic`/`.nodeset` keyword itself).
pub(super) fn parse_node_assignments(line: &str) -> Result<Vec<(String, f64)>, ParseError> {
    let raw: String = line
        .split_whitespace()
        .skip(1)
        .collect::<Vec<_>>()
        .join(" ");
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
        let name = if let Some(inner) = lhs_lc.strip_prefix("v(").and_then(|s| s.strip_suffix(')'))
        {
            inner.to_string()
        } else {
            lhs_lc.clone()
        };
        let value: f64 = parse_value(rhs, 0).unwrap_or_else(|_| rhs.parse::<f64>().unwrap_or(0.0));
        out.push((canon_node(&name), value));
    }
    Ok(out)
}

/// Parse `.options key=val key=val ...` into a list of `(key, value)` pairs.
///
/// Bare-flag tokens (no `=`) are stored as `("key", "1")` so `SimOptions::set`
/// can treat them as boolean true.  Quoted values are stripped of surrounding
/// quotes.  Returns an empty list for an empty directive line.
pub(super) fn parse_options_directive(line: &str) -> Vec<(String, String)> {
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
pub(super) fn parse_tran(line: &str, lineno: usize) -> Result<Analysis, ParseError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 3 {
        return Err(ParseError::FieldCount {
            expected: "≥3 (.tran step stop)",
            got: tokens.len(),
            line: lineno,
        });
    }
    let uic = tokens.iter().any(|t| t.eq_ignore_ascii_case("uic"));
    // `.tran step stop [tstart [tmax]] [UIC]`. `UIC` may sit in any trailing
    // slot, so read positionally but skip it rather than parsing it as a number.
    let mut trailing = tokens[3..]
        .iter()
        .filter(|t| !t.eq_ignore_ascii_case("uic"));
    let tstart = match trailing.next() {
        Some(t) => parse_value(t, lineno)?,
        None => 0.0,
    };
    let tmax = match trailing.next() {
        Some(t) => Some(parse_value(t, lineno)?),
        None => None,
    };
    Ok(Analysis::Tran {
        step: parse_value(tokens[1], lineno)?,
        stop: parse_value(tokens[2], lineno)?,
        tstart,
        tmax,
        uic,
    })
}

/// Parse `.dc SRC START STOP STEP [SRC2 START2 STOP2 STEP2]`.
pub(super) fn parse_dc(line: &str, lineno: usize) -> Result<Analysis, ParseError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 5 {
        return Err(ParseError::FieldCount {
            expected: "≥5 (.dc SRC START STOP STEP)",
            got: tokens.len(),
            line: lineno,
        });
    }
    let src = tokens[1].to_lowercase();
    let start = parse_value(tokens[2], lineno)?;
    let stop = parse_value(tokens[3], lineno)?;
    let step = parse_value(tokens[4], lineno)?;
    let nested = if tokens.len() >= 9 {
        Some(DcSweepSpec {
            src: tokens[5].to_lowercase(),
            start: parse_value(tokens[6], lineno)?,
            stop: parse_value(tokens[7], lineno)?,
            step: parse_value(tokens[8], lineno)?,
        })
    } else {
        None
    };
    Ok(Analysis::Dc {
        src,
        start,
        stop,
        step,
        nested,
    })
}

pub(super) fn parse_ac(line: &str, lineno: usize) -> Result<Analysis, ParseError> {
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
        other => {
            return Err(ParseError::Syntax {
                line: lineno,
                msg: format!("unknown AC variation '{other}'; expected dec, oct, or lin"),
            })
        }
    };
    let points = tokens[2].parse::<usize>().map_err(|_| ParseError::Syntax {
        line: lineno,
        msg: format!("invalid point count '{}'", tokens[2]),
    })?;
    Ok(Analysis::Ac {
        variation,
        points,
        fstart: parse_value(tokens[3], lineno)?,
        fstop: parse_value(tokens[4], lineno)?,
    })
}

/// Scan main-body lines for `.options enable_bidirectional=…` and return the
/// resulting flag.  Matches `enable_bidirectional`, `bidirectional`, and
/// `bidirectional_propagation` (mirrors `SimOptions::set`).  Last match wins.
pub(super) fn parse_noise(line: &str, lineno: usize) -> Result<Analysis, ParseError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 7 {
        return Err(ParseError::FieldCount {
            expected: "≥7 (.noise V(node) src DEC|OCT|LIN points fstart fstop)",
            got: tokens.len(),
            line: lineno,
        });
    }
    // tokens[0] = ".noise"; tokens[1] = "v(out[,ref])"; tokens[2] = src;
    // tokens[3] = DEC|OCT|LIN; tokens[4] = pts; tokens[5..7] = fstart, fstop.
    let v_expr = tokens[1];
    let inside = v_expr
        .strip_prefix("v(")
        .and_then(|s| s.strip_suffix(')'))
        .ok_or_else(|| ParseError::Syntax {
            line: lineno,
            msg: format!("expected V(node[,ref]); got '{v_expr}'"),
        })?;
    let (out_pos, out_neg) = if let Some((a, b)) = inside.split_once(',') {
        (canon_node(a.trim()), canon_node(b.trim()))
    } else {
        (canon_node(inside.trim()), "0".to_string())
    };
    let input_src = tokens[2].to_lowercase();
    let variation = match tokens[3] {
        "dec" => AcVariation::Dec,
        "oct" => AcVariation::Oct,
        "lin" => AcVariation::Lin,
        other => {
            return Err(ParseError::Syntax {
                line: lineno,
                msg: format!("unknown variation '{other}'; expected dec, oct, or lin"),
            })
        }
    };
    let points = tokens[4].parse::<usize>().map_err(|_| ParseError::Syntax {
        line: lineno,
        msg: format!("invalid point count '{}'", tokens[4]),
    })?;
    Ok(Analysis::Noise {
        out_pos,
        out_neg,
        input_src,
        variation,
        points,
        fstart: parse_value(tokens[5], lineno)?,
        fstop: parse_value(tokens[6], lineno)?,
    })
}
