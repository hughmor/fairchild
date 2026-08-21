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
//! | `LAMBDA(p,k)` | a generated `parameter real wl_k`, filled by resolution |
//!
//! # There is no λ accessor, and no passthrough to write
//!
//! There used to be two: `WL(p,k)` for the λ *wire* and `LAMBDA(p,k)` for the
//! *value*, with the author expected to pass the tag along by hand
//! (`OWL(WL(b,k)) <+ OWL(WL(a,k));`) and compute from the parameter. The wire
//! was the one real footgun in optical Verilog-A — propagation phase is
//! thousands of radians, so `∂φ/∂λ = φ/λ` is of order 1e9 per metre, and letting
//! the compiler differentiate that against a λ unknown stops Newton converging
//! at some wavelengths and not others.
//!
//! λ is no longer an unknown at all (see `fairchild_core::lambda`), so there is
//! nothing to differentiate and nothing to propagate: the elaborator knows every
//! λ net in the deck and fills `wl_k` from it. The λ ports are still declared —
//! the `X`-line ABI is positional and a bundle is `wpc·N` wires wide either way —
//! they simply carry no equation. A source that still writes `WL` is refused by
//! name rather than silently expanded into a contribution that goes nowhere.
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

    /// Terminal offset of each bundle port, in declaration-list order, for a
    /// module built at `n` channels with `wpc` wires per channel.
    ///
    /// Positional: this walks the *header* order, not the bundle order, because
    /// a scalar port between two bundles shifts everything after it. Getting it
    /// wrong points λ at a field wire, which reads as a wavelength of 0.06 m.
    fn bundle_offsets(&self, n: usize, wpc: usize) -> Vec<usize> {
        let mut offsets = Vec::with_capacity(self.bundles.len());
        let mut at = 0usize;
        for p in &self.header {
            if self.bundles.contains(p) {
                offsets.push(at);
                at += wpc * n;
            } else {
                at += 1;
            }
        }
        offsets
    }

    /// Every terminal of this module that carries a wavelength label: the last
    /// wire of every channel of every bundle port.
    pub fn lambda_terminals(&self, n: usize, wpc: usize) -> Vec<usize> {
        let lam = wpc - 1;
        self.bundle_offsets(n, wpc)
            .into_iter()
            .flat_map(|off| (0..n).map(move |k| off + wpc * k + lam))
            .collect()
    }

    /// How a label moves through this module: within a channel slot, between
    /// every pair of bundle ports, both ways.
    ///
    /// A bundle model has one channel grid — `LAMBDA(p,k)` expands to `wl_k`
    /// whatever `p` is, so slot `k` *is* a wavelength across the whole module —
    /// and that is what makes the routing derivable without asking the author
    /// which port is an input. Resolution takes the label from whichever port a
    /// source actually reached; the cycle the both-ways edges create terminates
    /// because revisiting a net with the value it already has changes nothing.
    pub fn lambda_routing(&self, n: usize, wpc: usize) -> Vec<(usize, usize)> {
        let lam = wpc - 1;
        let offsets = self.bundle_offsets(n, wpc);
        let mut pairs = Vec::new();
        for (i, &a) in offsets.iter().enumerate() {
            for (j, &b) in offsets.iter().enumerate() {
                if i == j {
                    continue;
                }
                for k in 0..n {
                    pairs.push((a + wpc * k + lam, b + wpc * k + lam));
                }
            }
        }
        pairs
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
    // `WL(p,k)` was the λ *wire* accessor, for a passthrough contribution the
    // author had to write. λ is not an unknown any more, so there is no wire to
    // contribute to and the line would expand into a contribution that goes
    // nowhere — a silently dead model rather than a broken one. Refuse by name.
    if decommented(src).contains("WL(") {
        return Err(OsdiError::Dialect {
            detail: format!(
                "module `{}` uses `WL(...)`, which no longer exists: a wavelength is \
                 resolved before the solve rather than solved for, so the λ wire carries \
                 no equation and nothing has to pass it along. Delete the \
                 `OWL(WL(out,k)) <+ OWL(WL(in,k));` line; read the wavelength with \
                 `LAMBDA(port, k)`, which the elaborator fills from the deck's sources.",
                m.name
            ),
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
                // One wavelength parameter per channel, filled from resolution
                // at build time (see `OsdiDevice::set_resolved_lambda`). The
                // default is what a channel no source reaches falls back to,
                // and what a deck sets when it wants to override.
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
    for (call, suffix) in [("E_RE", "re"), ("E_IM", "im")] {
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
        // λ is a parameter the elaborator fills, and there is no λ wire to
        // contribute to: the ports exist for the ABI and carry no equation.
        assert!(out.contains("parameter real wl_1 = wl_default"));
        assert!(
            !out.contains("OWL("),
            "nothing should drive a λ wire:\n{out}"
        );
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

    /// A source written against the old dialect must be refused by name. λ has
    /// no equation any more, so `OWL(WL(b,k)) <+ OWL(WL(a,k));` would expand
    /// into a contribution that goes nowhere — a model that compiles, runs, and
    /// propagates nothing.
    #[test]
    fn the_retired_lambda_wire_accessor_is_refused_by_name() {
        let src = WG.replace(
            "OF(E_RE(b, k)) <+ OF(E_RE(a, k)) * cos(LAMBDA(a, k));",
            "OWL(WL(b, k)) <+ OWL(WL(a, k));",
        );
        let m = scan(&src).unwrap().unwrap();
        let e = expand(&src, &m, 2, 3).unwrap_err();
        let msg = format!("{e}");
        assert!(msg.contains("WL("), "{msg}");
        assert!(msg.contains("LAMBDA"), "the fix must be named: {msg}");
    }

    /// λ terminal indices are positional, and a scalar port between two bundles
    /// shifts every terminal after it. `WG` has its scalar last; this checks the
    /// arithmetic against a header where it is not.
    #[test]
    fn lambda_terminals_and_routing_follow_the_header_order() {
        let src = "module m(a, ctrl, b);\n optical_bundle a, b;\n electrical ctrl;\nendmodule\n";
        let m = scan(src).unwrap().unwrap();
        // N=2, wpc=3: a occupies 0..6, ctrl is 6, b occupies 7..13.
        assert_eq!(m.lambda_terminals(2, 3), vec![2, 5, 9, 12]);
        let routing = m.lambda_routing(2, 3);
        // Both ways, per slot: a↔b for k=0 and k=1.
        assert_eq!(routing.len(), 4);
        assert!(routing.contains(&(2, 9)) && routing.contains(&(9, 2)));
        assert!(routing.contains(&(5, 12)) && routing.contains(&(12, 5)));
    }
}
