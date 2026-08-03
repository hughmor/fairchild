mod bundles;
mod common;
mod directives;
mod element;
mod subckt;
mod waveforms;

pub use bundles::{bundle_arity_for, BundleArity};

// Pull internal helpers into scope so parse_spice + the test module can call
// them unqualified (as they did before the split).
use bundles::{expand_bundle_ports, scan_bidirectional};
use common::{canon_node, expand_bus_vectors, parse_port_decl, parse_value};
use directives::{
    is_silent_directive, parse_ac, parse_dc, parse_measure, parse_node_assignments, parse_noise,
    parse_options_directive, parse_tran,
};
use element::{parse_element_expanded, parse_model};
use subckt::{collect_defs, expand_instance, substitute_params};

use crate::{Analysis, Element, Netlist, ParseError};
use std::collections::HashSet;
use std::path::{Path, PathBuf};

/// Parse a netlist from a string.
///
/// `.include` / `.lib` are resolved relative to the **current working
/// directory**, since a string has no source path. Prefer
/// [`parse_spice_file`] when the deck lives on disk — it resolves includes
/// relative to the deck, which is what a library of `.subckt` files wants.
pub fn parse_spice(input: &str) -> Result<Netlist, ParseError> {
    // Resolving here (rather than only in parse_spice_file) means a
    // programmatically-built deck — Python's `load_str`, a fit script's
    // card + netlist concatenation — can pull in shared PCell files too.
    let resolved = resolve_includes(input, None, 0)?;
    parse_resolved(&resolved)
}

fn parse_resolved(input: &str) -> Result<Netlist, ParseError> {
    let all_lines = logical_lines(input);
    let mut netlist = Netlist::default();

    if all_lines.is_empty() {
        return Ok(netlist);
    }

    // First logical line is the title — by SPICE convention.  But that
    // convention is a foot-gun for programmatically-generated netlists
    // (Python helpers, code-gen): forget the title comment and your first
    // `.optical_port` / `.options` / `R1 …` line gets eaten as the title
    // with no error, producing confusing failures downstream.
    //
    // Be forgiving: only consume the first line as a title when it doesn't
    // look like a directive or an element.  A leading `*` or `;` is the
    // unambiguous "this is a comment / title" marker; anything else that
    // starts with `.` (directive) or an alphabetic character (element
    // prefix) is parsed as a normal body line and the title stays empty.
    let first_trimmed = all_lines[0].1.trim_start();
    let first_is_titlish = first_trimmed.is_empty()
        || first_trimmed.starts_with('*')
        || first_trimmed.starts_with(';');
    let body_start = if first_is_titlish {
        netlist.title = all_lines[0].1.trim().to_string();
        1
    } else {
        // No title: start parsing from line 0.
        0
    };

    // Pass 1.
    let (subckt_defs, global_params, main_lines) = collect_defs(&all_lines[body_start..])?;

    // Pre-scan for `.options enable_bidirectional=…` so we know how to size
    // optical-port wires when we see `.optical_port` directives in pass 2.
    // The flag is sticky — first match wins; later `.options` lines that
    // override the value still take effect through `SimOptions::from_netlist`
    // at the consumer side, but the parser-emitted wire names are fixed at
    // parse time.
    let bidirectional = scan_bidirectional(&main_lines);

    // Pass 2: parse main body.
    let mut expanding: HashSet<String> = HashSet::new();
    let mut current_alter: Option<crate::AlterBlock> = None;

    for (lineno, line) in &main_lines {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('*') {
            continue;
        }

        let lc = trimmed.to_lowercase();

        if lc == ".end" {
            break;
        } else if lc.starts_with(".alter") {
            // Flush any current block, start a new one.  Label defaults to the
            // 1-based ordinal so each block gets a stable name even when the
            // user omits one.
            if let Some(b) = current_alter.take() {
                netlist.alters.push(b);
            }
            let label = lc
                .split_whitespace()
                .nth(1)
                .map(str::to_string)
                .unwrap_or_else(|| format!("alter{}", netlist.alters.len() + 1));
            current_alter = Some(crate::AlterBlock {
                label,
                element_overrides: Vec::new(),
                model_overrides: Vec::new(),
            });
        } else if lc.starts_with(".optical_port") {
            // .optical_port NAME [N]  — declare a bundle PORT that expands to
            // 3·N underlying wires (re/im/λ per channel), or 5·N when
            // bidirectional.  Unlike `.optical` / `.optical_bus` below, this
            // creates a *referenceable* port: subsequent X-element lines may
            // use NAME as a single net token, which `expand_bundle_ports`
            // flattens (bundle-aware devices) or replicates per channel
            // (scalar devices).  All underlying wires are also registered as
            // optical nets so the discipline check works as before.
            let (name, channels) = parse_port_decl(trimmed, ".optical_port", *lineno)?;
            let port = crate::BundlePort {
                name,
                channels,
                kind: crate::BundleKind::Optical { bidirectional },
            };
            for w in port.all_wires() {
                netlist.optical_nets.push(w);
            }
            netlist.bundle_ports.push(port);
        } else if lc.starts_with(".electrical_port") {
            // .electrical_port NAME [N]  — the electrical sibling of
            // `.optical_port`: one plain wire per channel (NAME_0 … NAME_{N-1}),
            // referenceable as a single net token on an X-line.  This is how a
            // bundle-aware device takes one control signal per WDM channel
            // without needing a per-N device variant.
            //
            // Deliberately NOT registered in `optical_nets` — these wires stay
            // in the electrical discipline, so ordinary V/I/R/C may drive them.
            let (name, channels) = parse_port_decl(trimmed, ".electrical_port", *lineno)?;
            netlist.bundle_ports.push(crate::BundlePort {
                name,
                channels,
                kind: crate::BundleKind::Electrical,
            });
        } else if lc.starts_with(".optical_bus") {
            // .optical_bus N re_base im_base wl_base
            // DISCIPLINE ANNOTATION ONLY — tags 3N generated net names as
            // optical:
            //   re_base_0 im_base_0 wl_base_0  re_base_1 im_base_1 wl_base_1 …
            // It creates no referenceable port, so X-lines must spell out the
            // individual wire names.  `.optical_port NAME N` is the modern
            // equivalent and does strictly more; this stays for older netlists.
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
            // DISCIPLINE ANNOTATION ONLY (see `.optical_bus` above) — tags
            // already-named nets as optical.  Must come before ".op" check.
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
        } else if lc.starts_with(".noise") {
            netlist.analyses.push(parse_noise(&lc, *lineno)?);
        } else if lc.starts_with(".temp") {
            // .temp <T1_celsius> [<T2_celsius> …] — one entry per sweep point.
            for tok in lc.split_whitespace().skip(1) {
                let c = parse_value(tok, *lineno)?;
                netlist.temps.push(c + 273.15);
            }
        } else if lc.starts_with(".model") {
            if let Some(card) = parse_model(&lc, *lineno)? {
                if let Some(alter) = current_alter.as_mut() {
                    alter.model_overrides.push(card);
                } else {
                    netlist.models.push(card);
                }
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
        } else if lc.starts_with(".meas") {
            // .measure or .meas — both accepted.  Parse failure here is
            // surfaced as a Syntax error, not silently ignored.
            netlist.measurements.push(parse_measure(trimmed, *lineno)?);
        } else if lc.starts_with('.') {
            if !is_silent_directive(&lc) {
                let directive = lc.split_whitespace().next().unwrap_or(&lc).to_string();
                return Err(ParseError::UnsupportedDirective {
                    directive,
                    line: *lineno,
                });
            }
        } else {
            // Element or instance line; substitute top-level params first.
            let substituted = substitute_params(trimmed, &global_params, *lineno)?;
            for base_el in parse_element_expanded(&substituted, *lineno)? {
                // Bundle-port expansion (B2): any token in an XOsdi nets list
                // that matches a declared `.optical_port` is replaced with its
                // (re,im,λ) underlying wires.  When at least one referenced
                // port has channels > 1, the instance replicates per channel.
                // Returns one XOsdi per channel; for non-XOsdi elements or
                // unreferenced ports, returns a single-element vec.
                // A subckt instance's port count picks its bundle semantics.
                let subckt_ports = match &base_el {
                    Element::XOsdi { model_name, .. } => {
                        subckt_defs.get(model_name).map(|d| d.port_count())
                    }
                    _ => None,
                };
                let expanded_elements =
                    expand_bundle_ports(base_el, &netlist.bundle_ports, subckt_ports, *lineno)?;

                for el in expanded_elements {
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
                            let def = subckt_defs.get(model_name).unwrap();
                            let flat = expand_instance(
                                model_name,
                                name,
                                nets,
                                params,
                                def,
                                &subckt_defs,
                                &global_params,
                                &mut expanding,
                                *lineno,
                            )?;
                            // Per-instance `.model` cards ride along with the
                            // flattened elements that reference them.
                            if let Some(alter) = current_alter.as_mut() {
                                alter.element_overrides.extend(flat.elements);
                                alter.model_overrides.extend(flat.models);
                            } else {
                                netlist.elements.extend(flat.elements);
                                netlist.models.extend(flat.models);
                            }
                        }
                    } else if let Some(alter) = current_alter.as_mut() {
                        alter.element_overrides.push(el);
                    } else {
                        netlist.elements.push(el);
                    }
                }
            }
        }
    }
    // Flush trailing .alter block on EOF (no terminator required).
    if let Some(b) = current_alter.take() {
        netlist.alters.push(b);
    }

    Ok(netlist)
}

/// Recursively expand `.include "file"` and `.lib 'file' section` lines.
///
/// `.include` splices the entire file inline.  `.lib 'file' SECTION` reads the
/// file and splices only the `.lib SECTION` … `.endl [SECTION]` block.  Top-
/// level `.lib SECTION` … `.endl` definition blocks (the 1-arg form, used
/// inside library files) are stripped out: they're only meaningful when
/// referenced by name from a `.lib 'file' SECTION` directive.
///
/// `base_dir` is used to resolve relative paths.
fn resolve_includes(
    input: &str,
    base_dir: Option<&Path>,
    depth: usize,
) -> Result<String, ParseError> {
    if depth > 16 {
        return Err(ParseError::Syntax {
            line: 0,
            msg: ".include nesting depth > 16 (circular include?)".into(),
        });
    }

    let mut out = String::with_capacity(input.len());
    let lines: Vec<&str> = input.lines().collect();
    let mut i = 0usize;
    while i < lines.len() {
        let raw = lines[i];
        let trimmed = raw.trim();
        let lc = trimmed.to_lowercase();
        let lineno = i + 1;

        if lc.starts_with(".include") {
            let tok: Vec<&str> = trimmed.splitn(2, char::is_whitespace).collect();
            if tok.len() < 2 {
                return Err(ParseError::Syntax {
                    line: lineno,
                    msg: ".include requires a filename argument".into(),
                });
            }
            let fname = tok[1].trim().trim_matches('"').trim_matches('\'');
            let path: PathBuf = match base_dir {
                Some(dir) => dir.join(fname),
                None => PathBuf::from(fname),
            };
            let content = std::fs::read_to_string(&path).map_err(|e| ParseError::Syntax {
                line: lineno,
                msg: format!(".include '{}': {e}", path.display()),
            })?;
            let inlined = resolve_includes(&content, path.parent(), depth + 1)?;
            out.push_str(&inlined);
            out.push('\n');
            i += 1;
            continue;
        }

        if lc.starts_with(".lib") {
            let toks: Vec<&str> = trimmed.split_whitespace().collect();
            // Two-arg form: `.lib 'file' section` (include a section).
            // One-arg form: `.lib section` (start of a definition block — skip
            //   until matching `.endl`).
            if toks.len() >= 3 {
                let fname = toks[1].trim_matches('"').trim_matches('\'');
                let section = toks[2];
                let path: PathBuf = match base_dir {
                    Some(dir) => dir.join(fname),
                    None => PathBuf::from(fname),
                };
                let content = std::fs::read_to_string(&path).map_err(|e| ParseError::Syntax {
                    line: lineno,
                    msg: format!(".lib '{}': {e}", path.display()),
                })?;
                let section_text =
                    extract_lib_section(&content, section).ok_or_else(|| ParseError::Syntax {
                        line: lineno,
                        msg: format!(".lib '{}' section '{}' not found", path.display(), section),
                    })?;
                let inlined = resolve_includes(&section_text, path.parent(), depth + 1)?;
                out.push_str(&inlined);
                out.push('\n');
                i += 1;
                continue;
            } else if toks.len() == 2 {
                // Ambiguous in SPICE: `.lib <section>` opens a definition block,
                // but `.lib <file>` is what someone writing an include means. A
                // filename-looking argument gets an error rather than silently
                // swallowing the rest of the deck up to a `.endl` that never comes.
                let arg = toks[1].trim_matches('"').trim_matches('\'');
                if arg.contains('/') || arg.contains('\\') || arg.rsplit('.').count() > 1 {
                    return Err(ParseError::UnsupportedDirective {
                        directive: format!(
                            ".lib {arg} (one argument that looks like a file — use \
                             `.lib '{arg}' <section>` to include a section, or \
                             `.include {arg}` for the whole file)"
                        ),
                        line: lineno,
                    });
                }
                // Definition form: skip lines until matching `.endl`.
                i += 1;
                while i < lines.len() {
                    let l = lines[i].trim().to_lowercase();
                    if l.starts_with(".endl") {
                        i += 1;
                        break;
                    }
                    i += 1;
                }
                continue;
            }
            // `.lib` with no arguments — fall through to be reported as a
            // syntax error by the main parser pass.
        }

        out.push_str(raw);
        out.push('\n');
        i += 1;
    }
    Ok(out)
}

/// Extract the text of `.lib <section_name>` ... `.endl [<section_name>]` from
/// a library file.  Returns `None` if the named section doesn't exist.
///
/// Section matching is case-insensitive.  If `.endl` carries a section name,
/// it must match; if not, the next `.endl` ends the current section.  Nested
/// `.lib` definitions are not supported (rare in practice).
fn extract_lib_section(content: &str, section: &str) -> Option<String> {
    let section_lc = section.to_lowercase();
    let mut out = String::new();
    let mut inside = false;
    for raw in content.lines() {
        let trimmed = raw.trim();
        let lc = trimmed.to_lowercase();

        if !inside {
            if lc.starts_with(".lib") {
                let toks: Vec<&str> = trimmed.split_whitespace().collect();
                if toks.len() == 2 && toks[1].to_lowercase() == section_lc {
                    inside = true;
                }
            }
            continue;
        }

        // inside the target section
        if lc.starts_with(".endl") {
            let toks: Vec<&str> = trimmed.split_whitespace().collect();
            // If the .endl carries a section name, require it to match.
            if toks.len() == 1 || toks[1].to_lowercase() == section_lc {
                return Some(out);
            }
            continue;
        }
        out.push_str(raw);
        out.push('\n');
    }
    if inside {
        Some(out)
    } else {
        None
    }
}

/// Parse a SPICE netlist file, resolving `.include` directives relative to
/// the file's parent directory.
pub fn parse_spice_file(path: &Path) -> Result<Netlist, ParseError> {
    let src = std::fs::read_to_string(path).map_err(|e| ParseError::Syntax {
        line: 0,
        msg: format!("cannot read '{}': {e}", path.display()),
    })?;
    let resolved = resolve_includes(&src, path.parent(), 0)?;
    parse_resolved(&resolved)
}

// ─── analysis directive parsers ───────────────────────────────────────────────

fn logical_lines(input: &str) -> Vec<(usize, String)> {
    let mut result: Vec<(usize, String)> = Vec::new();
    // Index of the last entry a `+` may continue: a real line, not a comment.
    // Without this, interleaving `* …` comment lines inside a continuation block
    // silently reattaches the rest of the block to the comment, so parameters
    // vanish with no error. That has bitten real decks (a `.model` card whose
    // sections were separated by comments); annotating a long parameter list is
    // too natural to punish.
    let mut last_real: Option<usize> = None;
    for (i, raw) in input.lines().enumerate() {
        let lineno = i + 1;
        let trimmed = raw.trim_start();

        if trimmed.starts_with('+') {
            if let Some(idx) = last_real {
                result[idx].1.push(' ');
                result[idx]
                    .1
                    .push_str(trimmed.strip_prefix('+').unwrap_or("").trim());
            }
            continue;
        }
        if trimmed.starts_with('*') || trimmed.starts_with(';') {
            // Keep the entry (the title may be a comment, and line numbers stay
            // meaningful) but don't let it capture later continuations.
            result.push((lineno, raw.to_string()));
            continue;
        }

        let prev_ends_backslash = last_real
            .map(|idx| result[idx].1.trim_end().ends_with('\\'))
            .unwrap_or(false);

        if prev_ends_backslash {
            let idx = last_real.expect("checked above");
            let without_bs = result[idx]
                .1
                .trim_end()
                .trim_end_matches('\\')
                .trim_end()
                .to_string();
            result[idx].1 = without_bs;
            result[idx].1.push(' ');
            result[idx].1.push_str(trimmed);
        } else {
            result.push((lineno, raw.to_string()));
            last_real = Some(result.len() - 1);
        }
    }
    result
}

// ─── tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BehavioralKind, Waveform};

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
            Analysis::Tran { step, stop, uic } => {
                assert!((step - 1e-6).abs() < 1e-12);
                assert!((stop - 5e-3).abs() < 1e-12);
                assert!(!uic);
            }
            _ => panic!("expected Tran analysis"),
        }
    }

    #[test]
    fn parse_pulse_waveform() {
        let input = "* Pulse\nV1 a 0 PULSE(0 1 0 1n 1n 10m 20m)\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        if let Element::VoltageSource {
            waveform: Waveform::Pulse { v0, v1, tr, .. },
            ..
        } = &netlist.elements[0]
        {
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
        assert_eq!(canon_node("0"), "0");
    }

    #[test]
    fn parse_diode_element() {
        let input = "* Diode\nD1 anode cathode myd\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(netlist.elements.len(), 1);
        if let Element::Diode {
            name,
            anode,
            cathode,
            model_name,
            ..
        } = &netlist.elements[0]
        {
            assert_eq!(name, "d1");
            assert_eq!(anode, "anode");
            assert_eq!(cathode, "cathode");
            assert_eq!(model_name, "myd");
        } else {
            panic!("expected Diode element");
        }
    }

    /// `D` used to stop reading at the model name, so trailing `key=value`
    /// pairs were dropped at parse time — silently, and with no way to
    /// parameterise an OSDI model instantiated as a diode.  `M` and `Q`
    /// always kept theirs.
    #[test]
    fn parse_diode_instance_params() {
        let input = "* Diode\nD1 a b myd Is=1e-12 Rs=0.5\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        let Element::Diode { params, .. } = &netlist.elements[0] else {
            panic!("expected Diode element");
        };
        assert_eq!(
            params,
            &[("is".to_string(), 1e-12), ("rs".to_string(), 0.5)]
        );
    }

    #[test]
    fn parse_model_card() {
        let input = "* test\n.model myd D (Is=1e-14 N=1)\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(netlist.models.len(), 1);
        let m = &netlist.models[0];
        assert_eq!(m.name, "myd");
        assert_eq!(m.kind, "d");
        let is = m
            .params
            .iter()
            .find(|(k, _)| k == "is")
            .map(|(_, v)| *v)
            .unwrap();
        assert!((is - 1e-14).abs() < 1e-20, "is={is}");
        let n = m
            .params
            .iter()
            .find(|(k, _)| k == "n")
            .map(|(_, v)| *v)
            .unwrap();
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
        let w = Waveform::Pulse {
            v0: 0.0,
            v1: 1.0,
            td: 0.0,
            tr: 1e-9,
            tf: 1e-9,
            pw: 1.0,
            per: 2.0,
        };
        assert!((w.at(0.0) - 0.0).abs() < 1e-12);
        assert!((w.at(1e-9) - 1.0).abs() < 1e-6);
        assert!((w.at(0.5e-9) - 0.5).abs() < 1e-6);
    }

    #[test]
    fn parse_pwl_waveform() {
        let input = "* PWL\nV1 a 0 PWL(0 0 1u 5 2u 5 3u 0)\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        if let Element::VoltageSource {
            waveform: Waveform::Pwl { points },
            ..
        } = &netlist.elements[0]
        {
            assert_eq!(points.len(), 4);
            assert!((points[1].0 - 1e-6).abs() < 1e-18);
            assert!((points[1].1 - 5.0).abs() < 1e-12);
        } else {
            panic!("expected PWL VoltageSource");
        }
    }

    #[test]
    fn pwl_waveform_at() {
        let w = Waveform::Pwl {
            points: vec![(0.0, 0.0), (1e-6, 5.0), (2e-6, 5.0), (3e-6, 0.0)],
        };
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
            Analysis::Ac {
                variation,
                points,
                fstart,
                fstop,
            } => {
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
            Analysis::Ac {
                variation, points, ..
            } => {
                assert_eq!(*variation, crate::AcVariation::Lin);
                assert_eq!(*points, 100);
            }
            _ => panic!("expected Ac analysis"),
        }
    }

    #[test]
    fn pulse_next_breakpoint_before_td() {
        let w = Waveform::Pulse {
            v0: 0.0,
            v1: 1.0,
            td: 1e-6,
            tr: 100e-9,
            tf: 100e-9,
            pw: 5e-6,
            per: 10e-6,
        };
        let bp = w.next_breakpoint(0.0).unwrap();
        assert!((bp - 1e-6).abs() < 1e-18, "expected td=1µs, got {bp}");
    }

    #[test]
    fn pulse_next_breakpoint_at_period_boundary() {
        let td = 0.0_f64;
        let tr = 100e-9_f64;
        let pw = 5e-6_f64;
        let tf = 100e-9_f64;
        let per = 10e-6_f64;
        let w = Waveform::Pulse {
            v0: 0.0,
            v1: 1.0,
            td,
            tr,
            tf,
            pw,
            per,
        };
        let t = td + per;
        let bp = w.next_breakpoint(t).unwrap();
        assert!(
            (bp - (t + tr)).abs() < 1e-18,
            "expected t+tr={}, got {bp}",
            t + tr
        );
    }

    #[test]
    fn pulse_next_breakpoint_mid_rise() {
        let w = Waveform::Pulse {
            v0: 0.0,
            v1: 1.0,
            td: 0.0,
            tr: 100e-9,
            tf: 100e-9,
            pw: 5e-6,
            per: 10e-6,
        };
        let bp = w.next_breakpoint(50e-9).unwrap();
        assert!((bp - 100e-9).abs() < 1e-18, "expected 100ns, got {bp}");
    }

    #[test]
    fn pulse_next_breakpoint_no_repeat_exhausted() {
        let w = Waveform::Pulse {
            v0: 0.0,
            v1: 1.0,
            td: 0.0,
            tr: 100e-9,
            tf: 100e-9,
            pw: 5e-6,
            per: 0.0,
        };
        let after_all = 100e-9 + 5e-6 + 100e-9 + 1e-9;
        assert!(w.next_breakpoint(after_all).is_none());
    }

    #[test]
    fn pwl_next_breakpoint() {
        let w = Waveform::Pwl {
            points: vec![(0.0, 0.0), (1e-6, 5.0), (2e-6, 5.0), (3e-6, 0.0)],
        };
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
        if let Element::XOsdi {
            name,
            nets,
            model_name,
            params,
        } = &netlist.elements[0]
        {
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
            Element::VoltageSource {
                waveform: Waveform::Sin { vo, va, freq, .. },
                ..
            } => {
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
            Element::VoltageSource {
                waveform:
                    Waveform::Exp {
                        v1,
                        v2,
                        td1,
                        tau1,
                        td2,
                        tau2,
                    },
                ..
            } => {
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
            Element::VoltageSource {
                waveform:
                    Waveform::Sffm {
                        vo,
                        va,
                        fc,
                        mdi,
                        fs,
                    },
                ..
            } => {
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
            Element::VoltageSource {
                waveform: Waveform::Am { va, vo, mf, fc, td },
                ..
            } => {
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
        let w = Waveform::Sin {
            vo: 0.5,
            va: 1.0,
            freq: 1e3,
            td: 0.0,
            theta: 0.0,
            phase: 0.0,
        };
        // sin(0) = 0 → vo
        assert!((w.at(0.0) - 0.5).abs() < 1e-12);
        // sin(π/2) at t = 0.25 ms with f=1kHz: sin(π/2)=1
        assert!((w.at(0.25e-3) - 1.5).abs() < 1e-6);
    }

    #[test]
    fn sin_pre_delay_is_vo() {
        let w = Waveform::Sin {
            vo: 0.5,
            va: 1.0,
            freq: 1e3,
            td: 1e-6,
            theta: 0.0,
            phase: 0.0,
        };
        assert!((w.at(0.0) - 0.5).abs() < 1e-12);
        assert!((w.at(0.5e-6) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn exp_pre_td1_is_v1() {
        let w = Waveform::Exp {
            v1: 0.0,
            v2: 1.0,
            td1: 1e-6,
            tau1: 1e-6,
            td2: 5e-6,
            tau2: 1e-6,
        };
        assert!((w.at(0.0) - 0.0).abs() < 1e-12);
        // At t=td1+tau1 the rise is (1-e^-1) = 0.6321
        assert!((w.at(2e-6) - (1.0 - (-1.0_f64).exp())).abs() < 1e-9);
    }

    #[test]
    fn first_line_directive_is_not_swallowed_as_title() {
        // Programmatically-generated netlists often skip the title comment.
        // A leading `.optical_port` (or any other directive / element) must
        // be parsed as such, not silently consumed as the title.
        let net = parse_spice(
            ".optical_port portA\n\
             .optical_port portB\n\
             Xwg portA portB some_model\n.end\n",
        )
        .unwrap();
        assert_eq!(net.title, "", "no title comment present");
        assert_eq!(net.bundle_ports.len(), 2, "both ports should parse");
        assert_eq!(net.elements.len(), 1, "Xwg should be present");
    }

    #[test]
    fn first_line_comment_is_consumed_as_title() {
        let net = parse_spice(
            "* Real title\n\
             V1 in 0 DC 1\n\
             R1 in 0 1k\n.op\n.end\n",
        )
        .unwrap();
        assert_eq!(net.title, "* Real title");
        // V1 and R1 both present in the body.
        assert_eq!(net.elements.len(), 2);
    }

    #[test]
    fn optical_port_single_channel_expands_3_wires() {
        let net = parse_spice(
            "* port test\n\
             .optical_port portin\n\
             .optical_port portout\n\
             Xwg portin portout some_model\n.end\n",
        )
        .unwrap();
        assert_eq!(net.bundle_ports.len(), 2);
        assert_eq!(net.bundle_ports[0].name, "portin");
        assert_eq!(net.bundle_ports[0].channels, 1);
        // Single XOsdi (max_n = 1); nets list is 6 wires (3 per port).
        assert_eq!(net.elements.len(), 1);
        match &net.elements[0] {
            Element::XOsdi {
                name,
                nets,
                model_name,
                ..
            } => {
                assert_eq!(name, "xwg");
                assert_eq!(model_name, "some_model");
                assert_eq!(
                    nets,
                    &vec![
                        "portin_re_0".to_string(),
                        "portin_im_0".to_string(),
                        "portin_wl_0".to_string(),
                        "portout_re_0".to_string(),
                        "portout_im_0".to_string(),
                        "portout_wl_0".to_string(),
                    ]
                );
            }
            _ => panic!("expected XOsdi"),
        }
        // All 6 underlying wires registered as optical.
        for w in [
            "portin_re_0",
            "portin_im_0",
            "portin_wl_0",
            "portout_re_0",
            "portout_im_0",
            "portout_wl_0",
        ] {
            assert!(net.optical_nets.iter().any(|n| n == w), "missing {w}");
        }
    }

    /// An unrecognised model on a multi-channel bundle is refused, not
    /// replicated.  Replication silently duplicated any electrical port N times
    /// onto the same nodes; nothing in the tree wanted it.
    #[test]
    fn unknown_model_on_a_wdm_bundle_is_refused() {
        let err = parse_spice(
            "* WDM port test\n\
             .optical_port bus_in 4\n\
             .optical_port bus_out 4\n\
             Xwg bus_in bus_out some_model L_um=100\n.end\n",
        )
        .unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("no WDM semantics"), "{msg}");
        assert!(msg.contains("bundle_arity_for"), "{msg}");
    }

    /// The same instance on 1-channel bundles expands in place — the case where
    /// replication and flattening always agreed.
    #[test]
    fn unknown_model_on_single_channel_bundles_expands_in_place() {
        let net = parse_spice(
            "* single channel\n\
             .optical_port bus_in\n\
             .optical_port bus_out\n\
             Xwg bus_in bus_out some_model L_um=100\n.end\n",
        )
        .unwrap();
        assert_eq!(net.elements.len(), 1);
        match &net.elements[0] {
            Element::XOsdi { name, nets, .. } => {
                assert_eq!(name, "xwg");
                assert_eq!(
                    nets,
                    &vec![
                        "bus_in_re_0".to_string(),
                        "bus_in_im_0".to_string(),
                        "bus_in_wl_0".to_string(),
                        "bus_out_re_0".to_string(),
                        "bus_out_im_0".to_string(),
                        "bus_out_wl_0".to_string(),
                    ]
                );
            }
            _ => panic!("expected XOsdi"),
        }
    }

    /// `.electrical_port NAME N` expands to one plain wire per channel, in
    /// place, alongside the optical bundles — this is what lets one
    /// bundle-aware instance take a per-WDM-channel control signal.
    #[test]
    fn electrical_port_expands_one_wire_per_channel() {
        let net = parse_spice(
            "* control bus test\n\
             .optical_port bus 3\n\
             .optical_port out 3\n\
             .electrical_port wctl 3\n\
             Xwb bus out wctl 0 fc_dcoupler kappa_L=0.1\n.end\n",
        )
        .unwrap();
        let ctl = net
            .bundle_ports
            .iter()
            .find(|p| p.name == "wctl")
            .expect("electrical port declared");
        assert_eq!(ctl.channels, 3);
        assert_eq!(ctl.wires_per_channel(), 1);
        assert!(!ctl.is_optical());
        // Electrical wires must NOT be filed as optical nets, or the discipline
        // check would reject an ordinary voltage source driving them.
        for w in ctl.all_wires() {
            assert!(!net.optical_nets.contains(&w), "{w} must stay electrical");
        }
        // fc_dcoupler is bundle-Aware → one instance, bundles flattened in
        // declaration order: 9 optical + 9 optical + 3 control + literal "0".
        assert_eq!(net.elements.len(), 1);
        let Element::XOsdi { nets, .. } = &net.elements[0] else {
            panic!("expected XOsdi");
        };
        assert_eq!(nets.len(), 9 + 9 + 3 + 1);
        assert_eq!(&nets[18..], &["wctl_0", "wctl_1", "wctl_2", "0"]);
    }

    /// The per-instance width check must span kinds: a 4-channel optical bus
    /// with a 2-wire control bus on the same bundle-aware device is meaningless
    /// and must fail at parse time, naming both ports.
    #[test]
    fn electrical_port_width_mismatch_with_optical_errors() {
        let err = parse_spice(
            "* bad widths\n\
             .optical_port bus 4\n\
             .optical_port out 4\n\
             .electrical_port wctl 2\n\
             Xwb bus out wctl 0 fc_dcoupler kappa_L=0.1\n.end\n",
        )
        .unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("same channel count"), "{msg}");
        assert!(msg.contains("bus(optical, 4 ch)"), "{msg}");
        assert!(msg.contains("wctl(electrical, 2 ch)"), "{msg}");
    }

    /// Differing widths are only an error *within one instance*. A netlist may
    /// mix a 4-channel optical bus and an 8-wire control bus that never meet.
    #[test]
    fn differing_bus_widths_are_fine_across_instances() {
        let net = parse_spice(
            "* independent buses\n\
             .optical_port bus 4\n\
             .optical_port out 4\n\
             .electrical_port wide 8\n\
             Xwg bus out fc_waveguide L_um=100\n\
             R1 wide_0 wide_7 1k\n.end\n",
        )
        .unwrap();
        assert_eq!(net.bundle_ports.len(), 3);
        assert_eq!(net.elements.len(), 2);
    }

    #[test]
    fn optical_port_mismatched_channels_error() {
        let res = parse_spice(
            "* bad\n\
             .optical_port a 2\n\
             .optical_port b 4\n\
             Xwg a b model\n.end\n",
        );
        assert!(
            res.is_err(),
            "should error on mismatched channel counts: {:?}",
            res.map(|n| n.elements.len())
        );
    }

    #[test]
    fn parse_alter_blocks_collect_overrides() {
        let input = "* alters\n\
                     V1 in 0 DC 1\nR1 in out 1k\n.op\n\
                     .alter slow\nR1 in out 2k\n\
                     .alter fast\nR1 in out 500\n\
                     .end\n";
        let net = parse_spice(input).unwrap();
        assert_eq!(net.elements.len(), 2, "base has V1 and R1");
        assert_eq!(net.alters.len(), 2);
        assert_eq!(net.alters[0].label, "slow");
        assert_eq!(net.alters[1].label, "fast");
        assert_eq!(net.alters[0].element_overrides.len(), 1);

        // Apply: base R1 = 1k, slow R1 = 2k.
        let mut applied = net.clone();
        applied.apply_alter(&net.alters[0]);
        let r1 = applied
            .elements
            .iter()
            .find_map(|e| match e {
                crate::Element::Resistor {
                    name, resistance, ..
                } if name == "r1" => Some(*resistance),
                _ => None,
            })
            .unwrap();
        assert!(
            (r1 - 2000.0).abs() < 1e-9,
            "expected R1=2k after alter, got {r1}"
        );
    }

    #[test]
    fn parse_dc_sweep_single() {
        let input = "* dc\nV1 in 0 DC 0\nR1 in 0 1k\n.dc V1 0 5 0.1\n.end\n";
        let netlist = parse_spice(input).unwrap();
        match &netlist.analyses[0] {
            Analysis::Dc {
                src,
                start,
                stop,
                step,
                nested,
            } => {
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
        assert_eq!(
            netlist.optical_nets,
            vec!["laser_re", "laser_im", "wg_out_re", "wg_out_im"]
        );
    }

    #[test]
    fn bus_vector_expansion_in_optical() {
        // .optical with bus vector notation
        let input = "* WDM test\n\
                     .optical opt_re[0..2] opt_im[0..2] opt_wl[0..2]\n\
                     .op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(
            netlist.optical_nets,
            vec![
                "opt_re_0", "opt_re_1", "opt_re_2", "opt_im_0", "opt_im_1", "opt_im_2", "opt_wl_0",
                "opt_wl_1", "opt_wl_2",
            ]
        );
    }

    #[test]
    fn optical_bus_directive() {
        // .optical_bus N re_base im_base wl_base
        let input = "* WDM test\n\
                     .optical_bus 3 ch_re ch_im ch_wl\n\
                     .op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        assert_eq!(
            netlist.optical_nets,
            vec![
                "ch_re_0", "ch_im_0", "ch_wl_0", "ch_re_1", "ch_im_1", "ch_wl_1", "ch_re_2",
                "ch_im_2", "ch_wl_2",
            ]
        );
    }

    #[test]
    fn bus_vector_expansion_in_xosdi_nets() {
        // X element with bus vector net arguments
        let input = "* WDM xosdi test\n\
                     .optical ch_re[0..1] ch_im[0..1] ch_wl[0..1]\n\
                     Xmux ch_re[0..1] ch_im[0..1] ch_wl[0..1] out_re out_im out_wl wdm_mux2\n\
                     .op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        if let Element::XOsdi {
            nets, model_name, ..
        } = &netlist.elements[0]
        {
            assert_eq!(model_name, "wdm_mux2");
            assert_eq!(
                nets,
                &[
                    "ch_re_0", "ch_re_1", "ch_im_0", "ch_im_1", "ch_wl_0", "ch_wl_1", "out_re",
                    "out_im", "out_wl",
                ]
            );
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
        let names: Vec<&str> = netlist
            .elements
            .iter()
            .map(|el| match el {
                Element::Resistor { name, .. } => name.as_str(),
                Element::VoltageSource { name, .. } => name.as_str(),
                _ => "?",
            })
            .collect();
        assert!(
            names.contains(&"xdiv1.r1"),
            "missing xdiv1.r1, got {names:?}"
        );
        assert!(
            names.contains(&"xdiv1.r2"),
            "missing xdiv1.r2, got {names:?}"
        );
        // Check node remapping: R1 should connect vdd → mid (port substitution).
        let r1 = netlist
            .elements
            .iter()
            .find(|el| matches!(el, Element::Resistor { name, .. } if name == "xdiv1.r1"))
            .unwrap();
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
        let res: Vec<_> = netlist
            .elements
            .iter()
            .filter(|el| matches!(el, Element::Resistor { .. }))
            .collect();
        assert_eq!(res.len(), 2);
        let n1_r = res.iter().find(|el| {
            if let Element::Resistor { name, .. } = el {
                name.starts_with("xinv1.")
            } else {
                false
            }
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
        let resistors: Vec<_> = netlist
            .elements
            .iter()
            .filter(|el| matches!(el, Element::Resistor { .. }))
            .collect();
        assert_eq!(resistors.len(), 2);
        let r_def = resistors
            .iter()
            .find(|el| {
                if let Element::Resistor { name, .. } = el {
                    name.starts_with("xdef.")
                } else {
                    false
                }
            })
            .unwrap();
        let r_over = resistors
            .iter()
            .find(|el| {
                if let Element::Resistor { name, .. } = el {
                    name.starts_with("xover.")
                } else {
                    false
                }
            })
            .unwrap();
        if let Element::Resistor { resistance, .. } = r_def {
            assert!((resistance - 1e3).abs() < 1e-9, "default R={resistance}");
        }
        if let Element::Resistor { resistance, .. } = r_over {
            assert!((resistance - 2e3).abs() < 1e-9, "override R={resistance}");
        }
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
        let res = netlist
            .elements
            .iter()
            .find(|el| matches!(el, Element::Resistor { .. }))
            .unwrap();
        if let Element::Resistor { resistance, .. } = res {
            assert!(
                (resistance - 4700.0).abs() < 1e-9,
                "global param R={resistance}"
            );
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
        let resistors: Vec<_> = netlist
            .elements
            .iter()
            .filter(|el| matches!(el, Element::Resistor { .. }))
            .collect();
        assert_eq!(
            resistors.len(),
            2,
            "expected 2 flat resistors, got {}",
            resistors.len()
        );
        let names: Vec<&str> = resistors
            .iter()
            .map(|el| {
                if let Element::Resistor { name, .. } = el {
                    name.as_str()
                } else {
                    ""
                }
            })
            .collect();
        assert!(
            names.iter().any(|n| n.starts_with("xout.r2")),
            "missing xout.r2: {names:?}"
        );
        assert!(
            names.iter().any(|n| n.starts_with("xout.xin.")),
            "missing xout.xin.*: {names:?}"
        );
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
        assert!(
            matches!(err, ParseError::SubcktCycle { .. }),
            "expected SubcktCycle, got {err:?}"
        );
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
        assert!(
            matches!(
                err,
                ParseError::SubcktPortCount {
                    expected: 2,
                    got: 4,
                    ..
                }
            ),
            "got {err:?}"
        );
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
        let res = netlist
            .elements
            .iter()
            .find(|el| matches!(el, Element::Resistor { .. }));
        assert!(res.is_some(), "forward-referenced subckt not expanded");
    }

    /// A `*` comment interleaved in a continuation block must not swallow the
    /// rest of the block. Annotating a long parameter list is natural, and the
    /// old behaviour silently dropped everything after the comment.
    #[test]
    fn comment_inside_continuation_block_is_transparent() {
        let net = parse_spice(
            "* deck\n\
             .model m fc_pn_ps LEVEL=4\n\
             * geometry\n\
             + l_m=1e-4 n_g=4.2\n\
             * electro-optic\n\
             + dn_dv=-3.6e-5\n\
             V1 a 0 DC 1\n\
             .op\n.end\n",
        )
        .unwrap();
        let card = &net.models[0];
        for (k, want) in [
            ("level", 4.0),
            ("l_m", 1e-4),
            ("n_g", 4.2),
            ("dn_dv", -3.6e-5),
        ] {
            let got = card
                .params
                .iter()
                .find(|(pk, _)| pk == k)
                .unwrap_or_else(|| panic!("missing {k} in {:?}", card.params))
                .1;
            assert!(
                (got - want).abs() <= 1e-12 * want.abs().max(1.0),
                "{k}: {got}"
            );
        }
    }

    /// PCell arithmetic: `{…}` may hold an expression over the instance's
    /// parameters, so a ring's length can derive from its radius instead of
    /// being pre-computed by the caller.
    #[test]
    fn subckt_param_expressions_evaluate() {
        let net = parse_spice(
            "* pcell arithmetic\n\
             .subckt ring a b radius=8e-6 n=2\n\
             R1 a b {2*pi*radius*n}\n\
             .ends\n\
             X1 p q ring radius=1e-5 n=3\n\
             .op\n.end\n",
        )
        .unwrap();
        let Element::Resistor { resistance, .. } = &net.elements[0] else {
            panic!("expected resistor, got {:?}", net.elements[0]);
        };
        let want = 2.0 * std::f64::consts::PI * 1e-5 * 3.0;
        assert!(
            (resistance - want).abs() < 1e-18,
            "got {resistance}, want {want}"
        );
    }

    #[test]
    fn subckt_param_expression_undefined_errors() {
        let err = parse_spice(
            "* bad param\n\
             .subckt s a b r=1\n\
             R1 a b {2*nope}\n\
             .ends\n\
             X1 p q s\n\
             .op\n.end\n",
        )
        .unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not finite"), "{msg}");
    }

    /// A `.model` card inside a `.subckt` becomes a PRIVATE, per-instance card:
    /// its name is mangled with the instance prefix and references from inside
    /// that instance are retargeted. Two instances therefore carry independent
    /// model parameters — the thing that makes a `.subckt` a real PCell, since
    /// LEVEL is only ever read from a card.
    #[test]
    fn subckt_model_cards_are_per_instance() {
        let net = parse_spice(
            "* per-instance cards\n\
             .subckt ps a b vpn gnd alpha=10 lvl=4\n\
             .model local_ps fc_pn_ps LEVEL={lvl} alpha_db_cm={alpha}\n\
             Xp a b vpn gnd local_ps\n\
             .ends\n\
             X1 i1 o1 v1 0 ps alpha=7\n\
             X2 i2 o2 v2 0 ps alpha=21\n\
             .op\n.end\n",
        )
        .unwrap();
        // Two independent cards, one per instance.
        assert_eq!(net.models.len(), 2, "{:?}", net.models);
        let mut by_name: Vec<(String, f64)> = net
            .models
            .iter()
            .map(|m| {
                let a = m
                    .params
                    .iter()
                    .find(|(k, _)| k == "alpha_db_cm")
                    .map(|(_, v)| *v)
                    .unwrap();
                (m.name.clone(), a)
            })
            .collect();
        by_name.sort_by(|a, b| a.0.cmp(&b.0));
        assert_eq!(by_name[0].0, "x1.local_ps");
        assert_eq!(by_name[1].0, "x2.local_ps");
        assert!((by_name[0].1 - 7.0).abs() < 1e-12, "{:?}", by_name);
        assert!((by_name[1].1 - 21.0).abs() < 1e-12, "{:?}", by_name);
        // LEVEL came through the parameter, so the card dispatches to LEVEL=4.
        for m in &net.models {
            let lvl = m.params.iter().find(|(k, _)| k == "level").unwrap().1;
            assert!((lvl - 4.0).abs() < 1e-12);
        }
        // And each instance's device points at its OWN card.
        let refs: Vec<&str> = net
            .elements
            .iter()
            .filter_map(|e| match e {
                Element::XOsdi { model_name, .. } => Some(model_name.as_str()),
                _ => None,
            })
            .collect();
        assert!(refs.contains(&"x1.local_ps"), "{refs:?}");
        assert!(refs.contains(&"x2.local_ps"), "{refs:?}");
    }

    /// `.include` now resolves for string input too, so a deck assembled in
    /// Python can pull in a shared PCell file.
    #[test]
    fn parse_spice_resolves_includes_from_string() {
        // Process-scoped, like `fc_lib_test_*` below: the fixed path made two
        // concurrent `cargo test` runs delete each other's include file.
        let dir = std::env::temp_dir().join(format!("fc_include_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let inc = dir.join("pcell_r.sp");
        std::fs::write(&inc, ".subckt twice a b r=1\nR1 a b {2*r}\n.ends\n").unwrap();
        let src = format!(
            "* including deck\n.include {}\nX1 p q twice r=50\n.op\n.end\n",
            inc.display()
        );
        let net = parse_spice(&src).unwrap();
        let Element::Resistor { resistance, .. } = &net.elements[0] else {
            panic!("expected the included subckt's resistor");
        };
        assert!((resistance - 100.0).abs() < 1e-12, "got {resistance}");
        std::fs::remove_file(&inc).ok();
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
    fn extract_lib_section_basic() {
        let lib_text = "\
* mylib
.lib tt
.model nmos NMOS Vto=0.5
.model pmos PMOS Vto=-0.5
.endl tt

.lib ff
.model nmos NMOS Vto=0.3
.endl ff
";
        let tt = extract_lib_section(lib_text, "tt").unwrap();
        assert!(tt.contains("Vto=0.5"));
        assert!(tt.contains("nmos"));
        assert!(!tt.contains("Vto=0.3"));

        let ff = extract_lib_section(lib_text, "ff").unwrap();
        assert!(ff.contains("Vto=0.3"));
        assert!(!ff.contains("Vto=0.5"));

        assert!(extract_lib_section(lib_text, "missing").is_none());
    }

    #[test]
    fn parse_spice_file_with_lib_section() {
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!("fc_lib_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let lib_path = dir.join("models.lib");
        let mut f = std::fs::File::create(&lib_path).unwrap();
        writeln!(f, ".lib tt").unwrap();
        writeln!(f, ".model myd D Is=1e-14 N=1").unwrap();
        writeln!(f, ".endl tt").unwrap();
        drop(f);

        let netlist_path = dir.join("test.sp");
        let mut nf = std::fs::File::create(&netlist_path).unwrap();
        writeln!(nf, "* lib test").unwrap();
        writeln!(nf, ".lib 'models.lib' tt").unwrap();
        writeln!(nf, "V1 a 0 DC 1").unwrap();
        writeln!(nf, "R1 a b 1k").unwrap();
        writeln!(nf, "D1 b 0 myd").unwrap();
        writeln!(nf, ".op").unwrap();
        writeln!(nf, ".end").unwrap();
        drop(nf);

        let netlist = parse_spice_file(&netlist_path).unwrap();
        // The .model from the lib file should have been spliced in.
        assert_eq!(netlist.models.len(), 1);
        assert_eq!(netlist.models[0].name, "myd");

        // Cleanup.
        let _ = std::fs::remove_file(&lib_path);
        let _ = std::fs::remove_file(&netlist_path);
        let _ = std::fs::remove_dir(&dir);
    }

    #[test]
    fn extract_lib_section_endl_without_name() {
        let lib_text = "\
.lib tt
.model nmos NMOS Vto=0.5
.endl
";
        let s = extract_lib_section(lib_text, "tt").unwrap();
        assert!(s.contains("Vto=0.5"));
    }

    #[test]
    fn parse_b_element_current() {
        let input = "* b\nV1 in 0 DC 1\nR1 in out 1k\nB1 out 0 I=V(in)*1m\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        let beh = netlist
            .elements
            .iter()
            .find(|e| matches!(e, Element::Behavioral { .. }))
            .unwrap();
        if let Element::Behavioral {
            name,
            pos,
            neg,
            kind,
            ..
        } = beh
        {
            assert_eq!(name, "b1");
            assert_eq!(pos, "out");
            assert_eq!(neg, "0");
            assert_eq!(*kind, BehavioralKind::Current);
        }
    }

    #[test]
    fn parse_b_element_voltage_with_spaces() {
        let input = "* b\nV1 in 0 DC 1\nR1 in out 1k\nB1 out 0 V = V(in) * 2\n.op\n.end\n";
        let netlist = parse_spice(input).unwrap();
        let beh = netlist
            .elements
            .iter()
            .find(|e| matches!(e, Element::Behavioral { .. }))
            .unwrap();
        if let Element::Behavioral { kind, .. } = beh {
            assert_eq!(*kind, BehavioralKind::Voltage);
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
        assert!(
            matches!(err, ParseError::Syntax { .. }),
            "expected Syntax error for stray .ends, got {err:?}"
        );
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
        let r = netlist
            .elements
            .iter()
            .filter(|el| matches!(el, Element::Resistor { .. }))
            .count();
        let x = netlist
            .elements
            .iter()
            .filter(|el| matches!(el, Element::XOsdi { .. }))
            .count();
        assert_eq!(r, 1, "expected 1 resistor from subckt expansion");
        assert_eq!(x, 1, "expected 1 XOsdi remaining");
    }

    // ── passive parasitic expansion ───────────────────────────────────────────

    #[test]
    fn inductor_rser_expands_to_two_elements() {
        let input = "V1 in 0 DC 1\nL1 in out 1m rser=10\n.op\n.end\n";
        let net = parse_spice(input).unwrap();
        // L1 expands into: Inductor(in → __l1_rn) + Resistor(__l1_rn → out)
        let inductors: Vec<_> = net
            .elements
            .iter()
            .filter(|e| matches!(e, Element::Inductor { .. }))
            .collect();
        let resistors: Vec<_> = net
            .elements
            .iter()
            .filter(|e| matches!(e, Element::Resistor { .. }))
            .collect();
        assert_eq!(inductors.len(), 1, "expected 1 inductor");
        assert_eq!(resistors.len(), 1, "expected 1 resistor (ESR)");
        if let Element::Inductor { pos, neg, .. } = inductors[0] {
            assert_eq!(pos, "in");
            assert_eq!(neg, "__l1_rn");
        }
        if let Element::Resistor {
            pos,
            neg,
            resistance,
            ..
        } = resistors[0]
        {
            assert_eq!(pos, "__l1_rn");
            assert_eq!(neg, "out");
            assert!((resistance - 10.0).abs() < 1e-12);
        }
    }

    #[test]
    fn inductor_rser_cpar_expands_to_three_elements() {
        let input = "V1 in 0 DC 1\nL1 in out 1m rser=5 cpar=1p\n.op\n.end\n";
        let net = parse_spice(input).unwrap();
        let n_l = net
            .elements
            .iter()
            .filter(|e| matches!(e, Element::Inductor { .. }))
            .count();
        let n_r = net
            .elements
            .iter()
            .filter(|e| matches!(e, Element::Resistor { .. }))
            .count();
        let n_c = net
            .elements
            .iter()
            .filter(|e| matches!(e, Element::Capacitor { .. }))
            .count();
        assert_eq!(n_l, 1);
        assert_eq!(n_r, 1);
        assert_eq!(n_c, 1);
        // cpar is across original pos/neg (in, out)
        let cap = net
            .elements
            .iter()
            .find(|e| matches!(e, Element::Capacitor { .. }))
            .unwrap();
        if let Element::Capacitor {
            pos,
            neg,
            capacitance,
            ..
        } = cap
        {
            assert_eq!(pos, "in");
            assert_eq!(neg, "out");
            assert!((capacitance - 1e-12).abs() < 1e-24);
        }
    }

    #[test]
    fn capacitor_esr_esl_expands_to_three_elements() {
        let input = "V1 in 0 DC 1\nC1 in 0 100n esr=0.1 esl=2n\n.op\n.end\n";
        let net = parse_spice(input).unwrap();
        let n_l = net
            .elements
            .iter()
            .filter(|e| matches!(e, Element::Inductor { .. }))
            .count();
        let n_r = net
            .elements
            .iter()
            .filter(|e| matches!(e, Element::Resistor { .. }))
            .count();
        let n_c = net
            .elements
            .iter()
            .filter(|e| matches!(e, Element::Capacitor { .. }))
            .count();
        assert_eq!(n_l, 1, "expected 1 inductor (ESL)");
        assert_eq!(n_r, 1, "expected 1 resistor (ESR)");
        assert_eq!(n_c, 1, "expected 1 capacitor");
        // Chain: in --[ESL]--> __c1_esln --[ESR]--> __c1_esrn --[C]--> 0
        let cap = net
            .elements
            .iter()
            .find(|e| matches!(e, Element::Capacitor { .. }))
            .unwrap();
        if let Element::Capacitor { pos, neg, .. } = cap {
            assert_eq!(pos, "__c1_esrn");
            assert_eq!(neg, "0");
        }
    }

    #[test]
    fn capacitor_rpar_adds_parallel_resistor() {
        let input = "V1 in 0 DC 1\nC1 in 0 100n rpar=1G\n.op\n.end\n";
        let net = parse_spice(input).unwrap();
        let n_r = net
            .elements
            .iter()
            .filter(|e| matches!(e, Element::Resistor { .. }))
            .count();
        assert_eq!(n_r, 1, "expected 1 resistor (rpar)");
        let r = net
            .elements
            .iter()
            .find(|e| matches!(e, Element::Resistor { .. }))
            .unwrap();
        if let Element::Resistor {
            pos,
            neg,
            resistance,
            ..
        } = r
        {
            assert_eq!(pos, "in");
            assert_eq!(neg, "0");
            assert!((resistance - 1e9).abs() < 1.0);
        }
    }

    #[test]
    fn resistor_cpar_adds_parallel_cap() {
        let input = "V1 in 0 DC 1\nR1 in 0 1k cpar=1p\n.op\n.end\n";
        let net = parse_spice(input).unwrap();
        let n_c = net
            .elements
            .iter()
            .filter(|e| matches!(e, Element::Capacitor { .. }))
            .count();
        assert_eq!(n_c, 1);
        let c = net
            .elements
            .iter()
            .find(|e| matches!(e, Element::Capacitor { .. }))
            .unwrap();
        if let Element::Capacitor {
            pos,
            neg,
            capacitance,
            ..
        } = c
        {
            assert_eq!(pos, "in");
            assert_eq!(neg, "0");
            assert!((capacitance - 1e-12).abs() < 1e-24);
        }
    }
}
