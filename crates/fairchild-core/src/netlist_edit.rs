//! Post-parse netlist edits shared by every binding.
//!
//! Both the Python and C bindings let a caller retarget an element value or
//! swap a source waveform without re-emitting SPICE text.  The parameter-name
//! aliases (`value` vs `resistance`, `dc` vs `v`) are API surface, so they live
//! here once rather than in each binding.

use fairchild_parser::{Element, Netlist, Waveform};

/// Set `param` on the element named `element` to `value`.
///
/// Matching is case-insensitive on both names.  Passives accept `value` or
/// their physical name (`resistance` / `capacitance` / `inductance`); sources
/// accept `value`, `dc`, or `v` / `i` and become a DC waveform; MOSFETs and
/// OSDI instances take any instance parameter, appended if not already present.
///
/// Returns `false` if no element matched — bindings should surface that rather
/// than silently ignoring a typo'd name.
pub fn set_element_param(netlist: &mut Netlist, element: &str, param: &str, value: f64) -> bool {
    let el_name = element.to_lowercase();
    let param_lc = param.to_lowercase();
    let mut hit = false;

    for el in &mut netlist.elements {
        match el {
            Element::Resistor {
                name, resistance, ..
            } if name.to_lowercase() == el_name
                && (param_lc == "resistance" || param_lc == "value") =>
            {
                *resistance = value;
                hit = true;
            }
            Element::Capacitor {
                name, capacitance, ..
            } if name.to_lowercase() == el_name
                && (param_lc == "capacitance" || param_lc == "value") =>
            {
                *capacitance = value;
                hit = true;
            }
            Element::Inductor {
                name, inductance, ..
            } if name.to_lowercase() == el_name
                && (param_lc == "inductance" || param_lc == "value") =>
            {
                *inductance = value;
                hit = true;
            }
            Element::VoltageSource { name, waveform, .. }
                if name.to_lowercase() == el_name
                    && (param_lc == "dc" || param_lc == "value" || param_lc == "v") =>
            {
                *waveform = Waveform::Dc(value);
                hit = true;
            }
            Element::CurrentSource { name, waveform, .. }
                if name.to_lowercase() == el_name
                    && (param_lc == "dc" || param_lc == "value" || param_lc == "i") =>
            {
                *waveform = Waveform::Dc(value);
                hit = true;
            }
            Element::XOsdi { name, params, .. } | Element::Mosfet { name, params, .. }
                if name.to_lowercase() == el_name =>
            {
                match params
                    .iter_mut()
                    .find(|(k, _)| k.to_lowercase() == param_lc)
                {
                    Some(entry) => entry.1 = value,
                    None => params.push((param.to_string(), value)),
                }
                hit = true;
            }
            _ => {}
        }
    }
    hit
}

/// Replace the waveform of voltage or current source `name` with a
/// piecewise-linear table.  `points` must be sorted by time.
///
/// Returns `false` if no source matched.
pub fn set_source_pwl(netlist: &mut Netlist, name: &str, points: Vec<(f64, f64)>) -> bool {
    let name_lc = name.to_lowercase();
    let mut hit = false;
    for el in &mut netlist.elements {
        match el {
            Element::VoltageSource { name, waveform, .. }
            | Element::CurrentSource { name, waveform, .. }
                if name.to_lowercase() == name_lc =>
            {
                *waveform = Waveform::Pwl {
                    points: points.clone(),
                };
                hit = true;
            }
            _ => {}
        }
    }
    hit
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

    fn nl() -> Netlist {
        parse_spice("* t\nV1 in 0 DC 1\nR1 in out 1k\nC1 out 0 1p\n.op\n.end\n").unwrap()
    }

    #[test]
    fn value_alias_and_case_insensitivity() {
        let mut n = nl();
        assert!(set_element_param(&mut n, "r1", "value", 2e3));
        assert!(set_element_param(&mut n, "R1", "resistance", 3e3));
        match &n.elements[1] {
            Element::Resistor { resistance, .. } => assert_eq!(*resistance, 3e3),
            other => panic!("expected R1, got {other:?}"),
        }
    }

    #[test]
    fn wrong_param_for_element_kind_is_not_a_match() {
        let mut n = nl();
        // `resistance` on a capacitor must not silently succeed.
        assert!(!set_element_param(&mut n, "c1", "resistance", 1.0));
        assert!(!set_element_param(&mut n, "nosuch", "value", 1.0));
        assert!(set_element_param(&mut n, "c1", "value", 2e-12));
    }

    #[test]
    fn source_becomes_dc_then_pwl() {
        let mut n = nl();
        assert!(set_element_param(&mut n, "v1", "dc", 3.3));
        assert!(matches!(
            &n.elements[0],
            Element::VoltageSource {
                waveform: Waveform::Dc(v),
                ..
            } if *v == 3.3
        ));
        assert!(set_source_pwl(&mut n, "V1", vec![(0.0, 0.0), (1e-9, 1.0)]));
        assert!(matches!(
            &n.elements[0],
            Element::VoltageSource {
                waveform: Waveform::Pwl { .. },
                ..
            }
        ));
        assert!(
            !set_source_pwl(&mut n, "r1", vec![(0.0, 0.0)]),
            "R1 is not a source"
        );
    }
}
