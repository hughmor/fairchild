//! Netlist sanity-check preflight.
//!
//! Catches obvious-but-fatal netlist errors before NR runs — the kind that
//! otherwise show up as "singular matrix" or "did not converge", forcing
//! the user to drag a full verbose validator across the whole circuit to
//! find the line they typed wrong.
//!
//! Items checked (all emit warnings to stderr; nothing aborts):
//!
//!   1. R = 0, R < 0, or non-finite resistance → stamp produces 1/0 or NaN.
//!   2. L = 0 or non-finite → companion model breaks at first transient step.
//!   3. C < 0 or non-finite → companion model breaks.
//!   4. V/I source with a non-finite amplitude anywhere in its waveform.
//!   5. Two-terminal element with `pos == neg` (shorted to itself).
//!   6. Two elements with the same refdes.
//!   7. Native `fc_*` X-element with a guaranteed-bad param value
//!      (`r_heater=0`, `p_pi=0`, `V_pi_L=0`, `r_shunt=0`).
//!
//! Silenced by `.options nosanitycheck=1` (or `sanity_check=0`).

use std::collections::HashMap;

use fairchild_parser::{Element, Netlist, Waveform};

/// A single sanity-check finding.  Grouped on `category` when reporting,
/// so that hundreds of "R=0" warnings collapse into one line plus a list
/// of refdes — instead of flooding stderr with identical text.
struct Warning {
    /// Short stable key identifying the kind of issue (e.g. "R=0",
    /// "duplicate-refdes", "fc-zero-param").
    category: &'static str,
    /// The element / refdes the finding applies to.
    refdes: String,
    /// Full human-readable detail line, used when this category has
    /// only a few instances and we print it verbatim.
    detail: String,
}

/// Run all sanity checks against `netlist` and print warnings to stderr,
/// grouped by category so a repeated issue (e.g. 70 zero-Ω resistors)
/// emits a single summary line instead of one warning per offender.
/// Returns the number of distinct findings (categories × refdes), which
/// callers can use in tests.
pub fn check_netlist_sanity(netlist: &Netlist) -> usize {
    let mut warnings: Vec<Warning> = Vec::new();

    check_passive_values(netlist, &mut warnings);
    check_source_amplitudes(netlist, &mut warnings);
    check_self_shorts(netlist, &mut warnings);
    check_duplicate_refdes(netlist, &mut warnings);
    check_fc_xosdi_params(netlist, &mut warnings);

    // Group by category, preserving first-seen order.
    const SHOW_FULL: usize = 3;
    let mut order: Vec<&'static str> = Vec::new();
    let mut by_cat: HashMap<&'static str, Vec<&Warning>> = HashMap::new();
    for w in &warnings {
        if !by_cat.contains_key(w.category) {
            order.push(w.category);
        }
        by_cat.entry(w.category).or_default().push(w);
    }
    for cat in order {
        let group = &by_cat[cat];
        if group.len() <= SHOW_FULL {
            for w in group {
                eprintln!("warning: sanity-check: {}", w.detail);
            }
        } else {
            // Print first SHOW_FULL verbatim, then a coalesced tail.
            for w in &group[..SHOW_FULL] {
                eprintln!("warning: sanity-check: {}", w.detail);
            }
            let extras: Vec<&str> = group[SHOW_FULL..]
                .iter()
                .map(|w| w.refdes.as_str())
                .collect();
            let truncated = if extras.len() > 10 {
                format!(
                    "{}, … (+{} more)",
                    extras[..10].join(", "),
                    extras.len() - 10
                )
            } else {
                extras.join(", ")
            };
            eprintln!(
                "warning: sanity-check: …and {} more with the same issue \
                       ('{cat}'): {truncated}",
                extras.len()
            );
        }
    }
    warnings.len()
}

fn check_passive_values(netlist: &Netlist, out: &mut Vec<Warning>) {
    for el in &netlist.elements {
        match el {
            Element::Resistor {
                name, resistance, ..
            } => {
                if !resistance.is_finite() {
                    out.push(Warning {
                        category: "R-nonfinite",
                        refdes: name.clone(),
                        detail: format!(
                            "{name}: R={resistance} is non-finite — \
                                         stamp will produce NaN/Inf"
                        ),
                    });
                } else if *resistance == 0.0 {
                    out.push(Warning {
                        category: "R=0",
                        refdes: name.clone(),
                        detail: format!(
                            "{name}: R=0 → stamp will produce 1/0 = ∞ \
                                         (use a small positive value like 1m, or \
                                         replace with a wire by merging the nets)"
                        ),
                    });
                } else if *resistance < 0.0 {
                    out.push(Warning {
                        category: "R-negative",
                        refdes: name.clone(),
                        detail: format!(
                            "{name}: R={resistance} is negative \
                                         (unusual; may be intentional for active \
                                         modelling — otherwise check sign)"
                        ),
                    });
                }
            }
            Element::Capacitor {
                name, capacitance, ..
            } => {
                if !capacitance.is_finite() {
                    out.push(Warning {
                        category: "C-nonfinite",
                        refdes: name.clone(),
                        detail: format!("{name}: C={capacitance} is non-finite"),
                    });
                } else if *capacitance < 0.0 {
                    out.push(Warning {
                        category: "C-negative",
                        refdes: name.clone(),
                        detail: format!("{name}: C={capacitance} is negative (unphysical)"),
                    });
                }
                // C=0 is a no-op (open circuit) — silent.
            }
            Element::Inductor {
                name, inductance, ..
            } => {
                if !inductance.is_finite() {
                    out.push(Warning {
                        category: "L-nonfinite",
                        refdes: name.clone(),
                        detail: format!("{name}: L={inductance} is non-finite"),
                    });
                } else if *inductance == 0.0 {
                    out.push(Warning {
                        category: "L=0",
                        refdes: name.clone(),
                        detail: format!(
                            "{name}: L=0 → transient companion uses \
                                         h/L which is undefined (DC OP treats it \
                                         as open — fine for `.op`, but transient \
                                         will fail)"
                        ),
                    });
                } else if *inductance < 0.0 {
                    out.push(Warning {
                        category: "L-negative",
                        refdes: name.clone(),
                        detail: format!("{name}: L={inductance} is negative (unphysical)"),
                    });
                }
            }
            _ => {}
        }
    }
}

fn check_source_amplitudes(netlist: &Netlist, out: &mut Vec<Warning>) {
    for el in &netlist.elements {
        let (name, wf) = match el {
            Element::VoltageSource { name, waveform, .. } => (name, waveform),
            Element::CurrentSource { name, waveform, .. } => (name, waveform),
            _ => continue,
        };
        if let Some(bad) = waveform_non_finite(wf) {
            out.push(Warning {
                category: "source-nonfinite",
                refdes: name.clone(),
                detail: format!("{name}: waveform has non-finite amplitude {bad}"),
            });
        }
    }
}

/// If any numeric amplitude in the waveform is NaN/Inf, return a short
/// description.  Timing parameters (td, tr, etc.) are not checked here —
/// they get caught at evaluation time and don't typically come from
/// foundry data.
fn waveform_non_finite(wf: &Waveform) -> Option<String> {
    let chk = |v: f64, what: &str| -> Option<String> {
        if !v.is_finite() {
            Some(format!("{what}={v}"))
        } else {
            None
        }
    };
    match wf {
        Waveform::Dc(v) => chk(*v, "DC"),
        Waveform::Pulse { v0, v1, .. } => chk(*v0, "v0").or_else(|| chk(*v1, "v1")),
        Waveform::Pwl { points } => {
            for (t, v) in points {
                if !t.is_finite() {
                    return Some(format!("PWL t={t}"));
                }
                if !v.is_finite() {
                    return Some(format!("PWL v={v}"));
                }
            }
            None
        }
        Waveform::Sin { vo, va, .. } => chk(*vo, "vo").or_else(|| chk(*va, "va")),
        Waveform::Exp { v1, v2, .. } => chk(*v1, "v1").or_else(|| chk(*v2, "v2")),
        Waveform::Sffm { vo, va, .. } => chk(*vo, "vo").or_else(|| chk(*va, "va")),
        Waveform::Am { vo, va, .. } => chk(*vo, "vo").or_else(|| chk(*va, "va")),
    }
}

fn check_self_shorts(netlist: &Netlist, out: &mut Vec<Warning>) {
    for el in &netlist.elements {
        let (name, pos, neg, kind) = match el {
            Element::Resistor { name, pos, neg, .. } => (name, pos.as_str(), neg.as_str(), "R"),
            Element::Capacitor { name, pos, neg, .. } => (name, pos.as_str(), neg.as_str(), "C"),
            Element::Inductor { name, pos, neg, .. } => (name, pos.as_str(), neg.as_str(), "L"),
            Element::VoltageSource { name, pos, neg, .. } => {
                (name, pos.as_str(), neg.as_str(), "V")
            }
            Element::CurrentSource { name, pos, neg, .. } => {
                (name, pos.as_str(), neg.as_str(), "I")
            }
            _ => continue,
        };
        if pos.eq_ignore_ascii_case(neg) {
            // V<name> a a DC 0 is the SPICE idiom for a 0-V ammeter, so don't
            // flag a zero-volt V source between matching nets.
            if let Element::VoltageSource {
                waveform: Waveform::Dc(v),
                ..
            } = el
            {
                if *v == 0.0 {
                    continue;
                }
            }
            out.push(Warning {
                category: "self-short",
                refdes: name.clone(),
                detail: format!(
                    "{name} ({kind}): both terminals tied to net \
                                 '{pos}' — element has no effect and wastes \
                                 an MNA row"
                ),
            });
        }
    }
}

fn check_duplicate_refdes(netlist: &Netlist, out: &mut Vec<Warning>) {
    let mut seen: HashMap<String, usize> = HashMap::new();
    for el in &netlist.elements {
        let name = element_name(el);
        let key = name.to_ascii_lowercase();
        *seen.entry(key).or_insert(0) += 1;
    }
    let mut dups: Vec<(&String, &usize)> = seen.iter().filter(|(_, &n)| n > 1).collect();
    dups.sort_by(|a, b| a.0.cmp(b.0));
    for (name, count) in dups {
        out.push(Warning {
            category: "duplicate-refdes",
            refdes: name.clone(),
            detail: format!(
                "refdes '{name}' appears {count} times — \
                             duplicate declarations; later ones may overwrite \
                             earlier ones silently"
            ),
        });
    }
}

fn element_name(el: &Element) -> &str {
    match el {
        Element::Resistor { name, .. } => name,
        Element::Capacitor { name, .. } => name,
        Element::Inductor { name, .. } => name,
        Element::CoupledInductors { name, .. } => name,
        Element::VoltageSource { name, .. } => name,
        Element::CurrentSource { name, .. } => name,
        Element::Diode { name, .. } => name,
        Element::Mosfet { name, .. } => name,
        Element::Bjt { name, .. } => name,
        Element::Behavioral { name, .. } => name,
        Element::XOsdi { name, .. } => name,
    }
}

/// Walk every X-element whose model name is a known fairchild native
/// device and verify that its critical params (the ones that appear in
/// `1/param` denominators) aren't zero. The list of "guaranteed-bad zero
/// params" is per-model so we only complain when it's actually fatal —
/// e.g. `n_g=0` is legal-but-unusual on a waveguide, but `r_heater=0` on
/// a thermal phase-shifter is a guaranteed `1/0` in the Jacobian.
fn check_fc_xosdi_params(netlist: &Netlist, out: &mut Vec<Warning>) {
    for el in &netlist.elements {
        let (name, model, params) = match el {
            Element::XOsdi {
                name,
                model_name,
                params,
                ..
            } => (name, model_name, params),
            _ => continue,
        };
        let model_lc = model.to_lowercase();
        let bad_keys = fatal_zero_params(&model_lc);
        if bad_keys.is_empty() {
            continue;
        }
        for (key, value) in params {
            let key_lc = key.to_lowercase();
            if bad_keys.iter().any(|&k| k == key_lc) {
                if *value == 0.0 {
                    out.push(Warning {
                        category: "fc-zero-param",
                        refdes: name.clone(),
                        detail: format!(
                            "{name} ({model}): {key}={value} → \
                                         stamp will produce 1/0 = ∞ (this param \
                                         appears in a reciprocal in the device's \
                                         eval/jacobian)"
                        ),
                    });
                } else if !value.is_finite() {
                    out.push(Warning {
                        category: "fc-nonfinite-param",
                        refdes: name.clone(),
                        detail: format!("{name} ({model}): {key}={value} is non-finite"),
                    });
                } else if *value < 0.0 {
                    out.push(Warning {
                        category: "fc-negative-param",
                        refdes: name.clone(),
                        detail: format!(
                            "{name} ({model}): {key}={value} is negative \
                                         (unphysical for this param)"
                        ),
                    });
                }
            }
        }
    }
}

/// Per-model list of param names that, when zero, cause a divide-by-zero
/// inside the device's eval or load_jacobian.  Names are lowercase to
/// match `set_real_param`'s canonical form.
fn fatal_zero_params(model_lc: &str) -> &'static [&'static str] {
    match model_lc {
        "fc_thermal_ps" | "fc_thermal_ps_rc" => &["r_heater", "r", "p_pi", "p_pi_w"],
        "fc_pn_th_ps" => &["r_heater", "r", "p_pi", "p_pi_w", "p_pi_th", "v_pi_l"],
        "fc_pn_ps" | "fc_pn_ps_cap" => &["v_pi_l"],
        "fc_photodetector" => &["r_shunt"],
        "fc_mzm" => &["v_pi"],
        _ => &[],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

    fn count_warnings(src: &str) -> usize {
        let net = parse_spice(src).unwrap();
        check_netlist_sanity(&net)
    }

    #[test]
    fn r_zero_warns() {
        let n = count_warnings("* r0\nV1 a 0 DC 1\nR1 a b 0.0\nR2 b 0 1k\n.op\n.end\n");
        assert!(n >= 1, "R=0 should be flagged");
    }

    #[test]
    fn r_negative_warns() {
        let n = count_warnings("* rneg\nV1 a 0 DC 1\nR1 a 0 -1k\n.op\n.end\n");
        assert!(n >= 1);
    }

    #[test]
    fn healthy_circuit_quiet() {
        let n = count_warnings("* clean\nV1 a 0 DC 1\nR1 a out 1k\nR2 out 0 1k\n.op\n.end\n");
        assert_eq!(n, 0, "clean circuit should produce no warnings");
    }

    #[test]
    fn duplicate_refdes_warns() {
        let n = count_warnings("* dup\nV1 a 0 DC 1\nR1 a b 1k\nR1 b 0 2k\n.op\n.end\n");
        assert!(n >= 1);
    }

    #[test]
    fn self_short_warns() {
        let n = count_warnings("* self\nV1 a 0 DC 1\nR1 a a 1k\nR2 a 0 1k\n.op\n.end\n");
        assert!(n >= 1);
    }

    #[test]
    fn zero_amp_v_source_self_short_silent() {
        // V1 a a DC 0 is the canonical ammeter idiom — must not warn.
        let n = count_warnings("* ammeter\nV1 a a DC 0\nR1 a 0 1k\nVdd a 0 DC 1\n.op\n.end\n");
        assert_eq!(n, 0);
    }

    #[test]
    fn fc_thermal_ps_zero_r_heater_warns() {
        let n = count_warnings(
            "* ths\n\
             .optical_port ch0\n.optical_port out0\n\
             Xl ch0 fc_cw_laser power_mW=1\n\
             Xt ch0 out0 hp 0 fc_thermal_ps r_heater=0 p_pi=10m\n\
             V1 hp 0 DC 1\n.op\n.end\n",
        );
        assert!(n >= 1, "r_heater=0 on fc_thermal_ps must warn");
    }
}
