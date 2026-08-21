//! One Verilog-A source, any channel count.
//!
//! A Verilog-A module declares its ports at compile time, so an 8-channel
//! optical bundle is 24 declared ports and a 100-channel one is 300. One
//! compiled `.osdi` therefore serves exactly one channel count, and "write once,
//! run at any N" cannot be a property of the artefact.
//!
//! It can be a property of the *source*. N is fixed the moment a deck says
//! `.optical_port bus 8`, so what is needed is elaboration-time polymorphism,
//! not runtime polymorphism — and since fairchild drives the compiler, producing
//! the artefact for the N in front of it is an implementation detail the author
//! never sees.
//!
//! # The dialect
//!
//! Four additions to Verilog-A, kept deliberately small: expansion is textual,
//! so the less of the body this understands, the less it can get wrong.
//!
//! | construct | expands to |
//! |---|---|
//! | `optical_bundle p, q;` | `3·N` (or `5·N`) `inout` ports per bundle, with disciplines |
//! | `N(p)` | the channel count, as a literal |
//! | `E_RE(p,k)` / `E_IM(p,k)` | `p_re_k` / `p_im_k` — the field wires |
//! | `WL(p,k)` | `p_wl_k` — the wavelength *wire*, for passing the tag through |
//! | `LAMBDA(p,k)` | a generated `parameter real wl_k` — the wavelength *value* |
//!
//! `WL` and `LAMBDA` are deliberately different things. Contributing from the λ
//! wire is the one real footgun in optical Verilog-A: propagation phase is
//! thousands of radians, so `∂φ/∂λ = φ/λ` is of order 1e9 per metre, and letting
//! the compiler differentiate that against the λ unknown stops Newton
//! converging at some wavelengths and not others. `LAMBDA` is a parameter, so it
//! has no derivative and cannot do that; `WL` exists only so a model can pass
//! the tag along and stay composable with native devices. The rule the user
//! guide currently states as advice — take the wavelength from a parameter,
//! never off the wire — becomes something the dialect expresses.
//!
//! # Not a Verilog-A parser
//!
//! `scan` looks for `module` and `optical_bundle` and nothing else; a source
//! without `optical_bundle` is not touched at all. `expand` rewrites the port
//! list, unrolls loops whose bound is `N(p)`, and substitutes accessors. Every
//! other line passes through byte for byte.
//!
//! The loop form is restricted on purpose:
//!
//! ```verilog
//! for (k = 0; k < N(bus); k = k + 1) begin … end
//! ```
//!
//! Anything else is refused by name rather than mis-generated, because a wrong
//! expansion is a silently wrong device.

use std::fmt::Write as _;

use crate::error::OsdiError;

/// What a scan found: a module with at least one bundle port.
#[derive(Debug, Clone, PartialEq)]
pub struct BundleModule {
    /// The module's name, as written. Registration uses this — a deck names the
    /// module, never a per-N mangling.
    pub name: String,
    /// Bundle ports, in declaration order.
    pub bundles: Vec<String>,
    /// Ordinary ports, in the order they appear in the module header.
    pub scalars: Vec<String>,
    /// The module header's full port list, in order, so expansion can rebuild
    /// it with each bundle replaced by its wires.
    header: Vec<String>,
}

impl BundleModule {
    /// Terminal count at `n` channels: every bundle contributes `wpc·n`.
    pub fn terminals_at(&self, n: usize, wpc: usize) -> usize {
        self.scalars.len() + self.bundles.len() * wpc * n
    }

    /// Solve `flattened` for a channel count, if one fits exactly.
    ///
    /// How the arity oracle places an instance of a bundle model: the width is
    /// whatever the deck declared, so rather than matching a fixed terminal
    /// count we ask whether *some* N produces this shape.
    pub fn channels_for(&self, flattened: usize, wpc: usize) -> Option<usize> {
        let per_n = self.bundles.len() * wpc;
        if per_n == 0 {
            return None;
        }
        let rest = flattened.checked_sub(self.scalars.len())?;
        if rest % per_n != 0 {
            return None;
        }
        let n = rest / per_n;
        (n >= 1).then_some(n)
    }
}

/// Strip `//` line comments — only for scanning, never for the emitted source.
fn decommented(src: &str) -> String {
    src.lines()
        .map(|l| match l.find("//") {
            Some(i) => &l[..i],
            None => l,
        })
        .collect::<Vec<_>>()
        .join("\n")
}

/// Find the bundle module in `src`, or `None` if this is ordinary Verilog-A.
///
/// Cheap enough to run on every `.va` a deck names: two substring searches
/// before any real work.
pub fn scan(src: &str) -> Result<Option<BundleModule>, OsdiError> {
    if !src.contains("optical_bundle") {
        return Ok(None);
    }
    let flat = decommented(src);

    // `module <name> ( <ports> );`
    let m = flat.find("module ").ok_or_else(|| OsdiError::Dialect {
        detail: "a source using `optical_bundle` has no `module` declaration".into(),
    })?;
    let after = &flat[m + "module ".len()..];
    let open = after.find('(').ok_or_else(|| OsdiError::Dialect {
        detail: "module declaration has no port list".into(),
    })?;
    let close = after.find(')').ok_or_else(|| OsdiError::Dialect {
        detail: "module port list is not closed".into(),
    })?;
    let name = after[..open].trim().to_string();
    let header: Vec<String> = after[open + 1..close]
        .split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();

    // `optical_bundle a, b;` — one or more declarations.
    let mut bundles: Vec<String> = Vec::new();
    for line in flat.lines() {
        let t = line.trim();
        let Some(rest) = t.strip_prefix("optical_bundle") else {
            continue;
        };
        let rest = rest.trim_end().trim_end_matches(';');
        for tok in rest.split(',') {
            let tok = tok.trim();
            if !tok.is_empty() {
                bundles.push(tok.to_string());
            }
        }
    }
    if bundles.is_empty() {
        return Ok(None);
    }
    for b in &bundles {
        if !header.contains(b) {
            return Err(OsdiError::Dialect {
                detail: format!(
                    "`optical_bundle {b}` is not in module `{name}`'s port list — a bundle \
                     must be a port, since its width comes from the deck"
                ),
            });
        }
    }
    let scalars = header
        .iter()
        .filter(|p| !bundles.contains(p))
        .cloned()
        .collect();
    Ok(Some(BundleModule {
        name,
        bundles,
        scalars,
        header,
    }))
}

/// Per-channel wire names for a bundle port, in the order the parser flattens
/// `.optical_port p N` into: `[re, im, (re_bw, im_bw,) λ]` per channel.
///
/// Positional and therefore load-bearing: this order is the ABI between a deck's
/// bundle declaration and a model's port list, and a mismatch is a wrong answer
/// rather than an error. `OpticalSegment::bind` is the other side of it.
fn wires(port: &str, n: usize, wpc: usize) -> Vec<String> {
    let mut out = Vec::with_capacity(n * wpc);
    for k in 0..n {
        out.push(format!("{port}_re_{k}"));
        out.push(format!("{port}_im_{k}"));
        if wpc == 5 {
            out.push(format!("{port}_re_bw_{k}"));
            out.push(format!("{port}_im_bw_{k}"));
        }
        out.push(format!("{port}_wl_{k}"));
    }
    out
}

/// Rewrite `src` for exactly `n` channels, `wpc` wires per channel.
///
/// The result is ordinary Verilog-A: the compiler sees nothing unusual, and the
/// content-hash cache works unchanged once N is part of the key.
pub fn expand(src: &str, m: &BundleModule, n: usize, wpc: usize) -> Result<String, OsdiError> {
    if n == 0 {
        return Err(OsdiError::Dialect {
            detail: format!("module `{}` cannot be built for 0 channels", m.name),
        });
    }
    let mut out = String::with_capacity(src.len() * 2);
    let uses_lambda = src.contains("LAMBDA(");

    let mut lines = src.lines().peekable();
    let mut header_done = false;
    while let Some(line) = lines.next() {
        let t = line.trim();

        // The module header: rebuild the port list with bundles expanded.
        if !header_done && t.starts_with("module ") {
            let mut ports: Vec<String> = Vec::new();
            for p in &m.header {
                if m.bundles.contains(p) {
                    ports.extend(wires(p, n, wpc));
                } else {
                    ports.push(p.clone());
                }
            }
            writeln!(out, "module {}({});", m.name, ports.join(", ")).unwrap();
            // Declare the bundle wires and give them their disciplines.
            for b in &m.bundles {
                let mut field = Vec::new();
                let mut lam = Vec::new();
                for k in 0..n {
                    field.push(format!("{b}_re_{k}"));
                    field.push(format!("{b}_im_{k}"));
                    if wpc == 5 {
                        field.push(format!("{b}_re_bw_{k}"));
                        field.push(format!("{b}_im_bw_{k}"));
                    }
                    lam.push(format!("{b}_wl_{k}"));
                }
                writeln!(out, "    inout {};", field.join(", ")).unwrap();
                writeln!(out, "    inout {};", lam.join(", ")).unwrap();
                writeln!(out, "    optical_field {};", field.join(", ")).unwrap();
                writeln!(out, "    optical_lambda {};", lam.join(", ")).unwrap();
            }
            if uses_lambda {
                // One wavelength parameter per channel. Until λ is resolved
                // before the solve, the deck sets these; after, the elaborator
                // fills them and the author's source does not change.
                writeln!(out, "    parameter real wl_default = 1550e-9 from (0:inf);").unwrap();
                for k in 0..n {
                    writeln!(out, "    parameter real wl_{k} = wl_default from (0:inf);").unwrap();
                }
            }
            header_done = true;
            continue;
        }

        // The bundle declarations themselves are consumed by the header rewrite.
        if t.starts_with("optical_bundle") {
            continue;
        }

        // A loop over a bundle's channels: unroll it.
        if let Some((var, bundle)) = parse_channel_loop(t) {
            if !m.bundles.contains(&bundle) {
                return Err(OsdiError::Dialect {
                    detail: format!(
                        "loop bound `N({bundle})` in module `{}` names something that is not \
                         an optical_bundle port",
                        m.name
                    ),
                });
            }
            let body = take_begin_end(&mut lines, &m.name)?;
            for k in 0..n {
                for bl in &body {
                    out.push_str(&substitute(bl, m, &var, k, n, wpc));
                    out.push('\n');
                }
            }
            continue;
        }

        // Everything else: accessors resolved, otherwise byte for byte. `N(p)`
        // outside a loop is legitimate (a bound, a scale factor), so it is
        // substituted here too; a bare accessor outside a loop has no channel
        // index and is caught by `substitute`.
        out.push_str(&substitute(line, m, "\u{0}", 0, n, wpc));
        out.push('\n');
    }
    Ok(out)
}

/// `for (k = 0; k < N(bus); k = k + 1) begin` → `("k", "bus")`.
fn parse_channel_loop(t: &str) -> Option<(String, String)> {
    let rest = t.strip_prefix("for")?.trim_start().strip_prefix('(')?;
    // The loop head's own `)` — not the one inside `N(bus)`, which is why this
    // counts depth instead of taking the first close paren.
    let mut depth = 1usize;
    let mut end = None;
    for (i, c) in rest.char_indices() {
        match c {
            '(' => depth += 1,
            ')' => {
                depth -= 1;
                if depth == 0 {
                    end = Some(i);
                    break;
                }
            }
            _ => {}
        }
    }
    let head = &rest[..end?];
    let parts: Vec<&str> = head.split(';').collect();
    if parts.len() != 3 {
        return None;
    }
    let var = parts[0].split('=').next()?.trim().to_string();
    let cond = parts[1].trim();
    let n_call = cond.split('<').nth(1)?.trim();
    let bundle = n_call
        .strip_prefix("N(")?
        .trim_end()
        .strip_suffix(')')?
        .trim()
        .to_string();
    Some((var, bundle))
}

/// Collect a `begin … end` body, the `begin` having been on the `for` line.
fn take_begin_end<'a>(
    lines: &mut std::iter::Peekable<impl Iterator<Item = &'a str>>,
    module: &str,
) -> Result<Vec<String>, OsdiError> {
    let mut depth = 1usize;
    let mut body = Vec::new();
    for line in lines.by_ref() {
        let t = line.trim();
        // Count nested blocks so an inner `begin … end` does not close the loop.
        if t == "begin" || t.ends_with(" begin") || t.ends_with(")begin") {
            depth += 1;
        }
        if t == "end" || t == "end;" {
            depth -= 1;
            if depth == 0 {
                return Ok(body);
            }
        }
        body.push(line.to_string());
    }
    Err(OsdiError::Dialect {
        detail: format!(
            "module `{module}`: a channel loop's `begin` is never closed by a matching `end`"
        ),
    })
}

/// Resolve `N(p)` and the per-channel accessors on one line.
fn substitute(line: &str, m: &BundleModule, var: &str, k: usize, n: usize, wpc: usize) -> String {
    let mut s = line.to_string();
    for b in &m.bundles {
        s = s.replace(&format!("N({b})"), &n.to_string());
    }
    // Accessors take (port, index). The index is the loop variable, or a
    // literal, so a model may address one channel without a loop.
    for (call, suffix) in [("E_RE", "re"), ("E_IM", "im"), ("WL", "wl")] {
        s = replace_accessor(&s, call, m, var, k, |port, idx| {
            format!("{port}_{suffix}_{idx}")
        });
    }
    s = replace_accessor(&s, "LAMBDA", m, var, k, |_port, idx| format!("wl_{idx}"));
    if wpc == 5 {
        s = replace_accessor(&s, "E_RE_BW", m, var, k, |port, idx| {
            format!("{port}_re_bw_{idx}")
        });
        s = replace_accessor(&s, "E_IM_BW", m, var, k, |port, idx| {
            format!("{port}_im_bw_{idx}")
        });
    }
    s
}

/// Replace every `CALL(port, index)` using `render`.
fn replace_accessor(
    line: &str,
    call: &str,
    m: &BundleModule,
    var: &str,
    k: usize,
    render: impl Fn(&str, usize) -> String,
) -> String {
    let pat = format!("{call}(");
    let mut out = String::with_capacity(line.len());
    let mut rest = line;
    while let Some(i) = rest.find(&pat) {
        // Only at a token boundary. `OWL(WL(b, k))` contains `WL(` twice — once
        // inside the access function — and matching the inner one swallowed the
        // real accessor, leaving the λ passthrough unexpanded.
        let boundary = i == 0
            || !rest[..i]
                .chars()
                .next_back()
                .is_some_and(|c| c.is_alphanumeric() || c == '_');
        if !boundary {
            let cut = i + 1;
            out.push_str(&rest[..cut]);
            rest = &rest[cut..];
            continue;
        }
        let (before, after) = rest.split_at(i);
        out.push_str(before);
        let inner_start = i + pat.len();
        let Some(close_rel) = rest[inner_start..].find(')') else {
            out.push_str(&rest[i..]);
            return out;
        };
        let args = &rest[inner_start..inner_start + close_rel];
        let mut it = args.split(',');
        let port = it.next().unwrap_or("").trim();
        let idx_tok = it.next().unwrap_or("").trim();
        let idx = if idx_tok == var {
            Some(k)
        } else {
            idx_tok.parse::<usize>().ok()
        };
        match (m.bundles.iter().any(|b| b == port), idx) {
            (true, Some(idx)) => out.push_str(&render(port, idx)),
            // Not ours, or an index we cannot resolve: leave it alone so the
            // compiler reports it against the author's own text.
            _ => out.push_str(&rest[i..inner_start + close_rel + 1]),
        }
        rest = &rest[inner_start + close_rel + 1..];
        let _ = after;
    }
    out.push_str(rest);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const WG: &str = r#"
`include "disciplines.vams"
module wg_wdm(a, b, ctrl);
    optical_bundle a, b;
    electrical ctrl;
    parameter real l_um = 1000.0;
    integer k;
    analog begin
        for (k = 0; k < N(a); k = k + 1) begin
            OF(E_RE(b, k)) <+ OF(E_RE(a, k)) * cos(LAMBDA(a, k));
            OWL(WL(b, k)) <+ OWL(WL(a, k));
        end
    end
endmodule
"#;

    #[test]
    fn a_source_without_the_dialect_is_not_touched() {
        assert_eq!(scan("module plain(a, b); endmodule").unwrap(), None);
    }

    #[test]
    fn scan_finds_the_bundles_and_the_scalars() {
        let m = scan(WG).unwrap().unwrap();
        assert_eq!(m.name, "wg_wdm");
        assert_eq!(m.bundles, vec!["a", "b"]);
        assert_eq!(m.scalars, vec!["ctrl"]);
        // Two bundles × 3 wires × N, plus the one electrical port.
        assert_eq!(m.terminals_at(4, 3), 1 + 2 * 3 * 4);
        assert_eq!(m.channels_for(1 + 2 * 3 * 4, 3), Some(4));
        // A shape no N produces must not be claimed.
        assert_eq!(m.channels_for(1 + 2 * 3 * 4 - 1, 3), None);
    }

    #[test]
    fn a_bundle_that_is_not_a_port_is_refused() {
        let src = "module m(a);\n optical_bundle a, ghost;\nendmodule\n";
        let e = scan(src).unwrap_err();
        assert!(format!("{e}").contains("ghost"), "{e}");
    }

    #[test]
    fn expansion_unrolls_the_loop_and_resolves_every_accessor() {
        let m = scan(WG).unwrap().unwrap();
        let out = expand(WG, &m, 2, 3).unwrap();
        // Ports: both bundles' wires, in channel order, plus the scalar.
        assert!(
            out.contains("module wg_wdm(a_re_0, a_im_0, a_wl_0, a_re_1, a_im_1, a_wl_1, b_re_0")
        );
        assert!(out.contains(", ctrl)"), "the scalar port survives:\n{out}");
        // The loop is gone and each channel got its own statement.
        assert!(!out.contains("for ("), "loop should be unrolled:\n{out}");
        assert!(out.contains("OF(b_re_0) <+ OF(a_re_0) * cos(wl_0);"));
        assert!(out.contains("OF(b_re_1) <+ OF(a_re_1) * cos(wl_1);"));
        // λ: the wire is passed through, the value is a parameter. Confusing
        // the two is the one real footgun in optical Verilog-A.
        assert!(out.contains("OWL(b_wl_1) <+ OWL(a_wl_1);"));
        assert!(out.contains("parameter real wl_1 = wl_default"));
        // Disciplines declared for the generated wires.
        assert!(out.contains("optical_field"));
        assert!(out.contains("optical_lambda"));
    }

    #[test]
    fn one_channel_expands_to_the_scalar_case() {
        let m = scan(WG).unwrap().unwrap();
        let out = expand(WG, &m, 1, 3).unwrap();
        assert!(out.contains("module wg_wdm(a_re_0, a_im_0, a_wl_0, b_re_0, b_im_0, b_wl_0, ctrl)"));
        assert!(out.contains("OF(b_re_0) <+ OF(a_re_0) * cos(wl_0);"));
        assert!(!out.contains("wl_1"), "no second channel at N=1:\n{out}");
    }

    #[test]
    fn bidirectional_expansion_carries_the_backward_pair() {
        let m = scan(WG).unwrap().unwrap();
        let out = expand(WG, &m, 1, 5).unwrap();
        assert!(out.contains("a_re_bw_0"), "backward wires present:\n{out}");
        assert!(out.contains("a_wl_0"));
    }

    #[test]
    fn an_unclosed_channel_loop_is_refused_by_name() {
        let src = "module m(a);\n optical_bundle a;\n analog begin\n \
                   for (k = 0; k < N(a); k = k + 1) begin\n x = 1;\n endmodule\n";
        let m = scan(src).unwrap().unwrap();
        let e = expand(src, &m, 2, 3).unwrap_err();
        assert!(format!("{e}").contains("never closed"), "{e}");
    }

    #[test]
    fn zero_channels_is_refused() {
        let m = scan(WG).unwrap().unwrap();
        assert!(expand(WG, &m, 0, 3).is_err());
    }
}
