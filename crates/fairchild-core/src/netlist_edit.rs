//! Post-parse netlist edits shared by every binding.
//!
//! Both the Python and C bindings let a caller retarget an element value or
//! swap a source waveform without re-emitting SPICE text.  The parameter-name
//! aliases (`value` vs `resistance`, `dc` vs `v`) are API surface, so they live
//! here once rather than in each binding.

use fairchild_parser::{Element, Netlist, Waveform};

/// Set `param` on the element named `element` to `value`.
///
/// Matching is case-insensitive on both names.  Passives accept `value`, their
/// physical name (`resistance` / `capacitance` / `inductance`), or the bare
/// element letter (`r` / `c` / `l`); sources accept `value`, `dc`, or `v` / `i`
/// and become a DC waveform; MOSFET, BJT, diode and OSDI instances take any
/// instance parameter, appended if not already present.
///
/// The device-with-params arm is what lets a Verilog-A transistor be swept: a
/// `.model`-card OSDI device is instantiated as an ordinary `M`/`Q`/`D` line,
/// so reaching only `X` elements would miss exactly the PDK idiom.
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
                && matches!(param_lc.as_str(), "resistance" | "value" | "r") =>
            {
                *resistance = value;
                hit = true;
            }
            Element::Capacitor {
                name, capacitance, ..
            } if name.to_lowercase() == el_name
                && matches!(param_lc.as_str(), "capacitance" | "value" | "c") =>
            {
                *capacitance = value;
                hit = true;
            }
            Element::Inductor {
                name, inductance, ..
            } if name.to_lowercase() == el_name
                && matches!(param_lc.as_str(), "inductance" | "value" | "l") =>
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
            // `Diode` and `Bjt` are here as of #26: their `params` lists used to
            // be carried and consumed by nothing, so matching them would have
            // returned `true` for an edit that changed nothing. Now `AREA`
            // reaches the device and anything else is named on stderr when the
            // device is built, which is the honest pair — a sweep over `area`
            // works, and a sweep over a parameter no model honours says so.
            Element::XOsdi { name, params, .. }
            | Element::Mosfet { name, params, .. }
            | Element::Diode { name, params, .. }
            | Element::Bjt { name, params, .. }
                if name.to_lowercase() == el_name =>
            {
                match params
                    .iter_mut()
                    .find(|(k, _)| k.to_lowercase() == param_lc)
                {
                    Some(entry) => entry.1 = value,
                    // Stored lower-cased so a second call with different casing
                    // updates this entry instead of appending a rival one.
                    None => params.push((param_lc.clone(), value)),
                }
                hit = true;
            }
            _ => {}
        }
    }
    hit
}

/// Read back what [`set_element_param`] would overwrite.
///
/// Accepts the same names and aliases, so a caller can perturb a parameter and
/// restore it without tracking the nominal itself — which is what the adjoint
/// sensitivity path needs to size its finite-difference step.
///
/// Returns `None` when the element or parameter is not found.  For a device
/// instance that means the parameter is not on the instance line: a `.model`
/// card default is invisible here, because the card is not the element.
pub fn get_element_param(netlist: &Netlist, element: &str, param: &str) -> Option<f64> {
    let el_name = element.to_lowercase();
    let param_lc = param.to_lowercase();

    for el in &netlist.elements {
        match el {
            Element::Resistor {
                name, resistance, ..
            } if name.to_lowercase() == el_name
                && matches!(param_lc.as_str(), "resistance" | "value" | "r") =>
            {
                return Some(*resistance);
            }
            Element::Capacitor {
                name, capacitance, ..
            } if name.to_lowercase() == el_name
                && matches!(param_lc.as_str(), "capacitance" | "value" | "c") =>
            {
                return Some(*capacitance);
            }
            Element::Inductor {
                name, inductance, ..
            } if name.to_lowercase() == el_name
                && matches!(param_lc.as_str(), "inductance" | "value" | "l") =>
            {
                return Some(*inductance);
            }
            Element::VoltageSource { name, waveform, .. }
                if name.to_lowercase() == el_name
                    && (param_lc == "dc" || param_lc == "value" || param_lc == "v") =>
            {
                return dc_level(waveform);
            }
            Element::CurrentSource { name, waveform, .. }
                if name.to_lowercase() == el_name
                    && (param_lc == "dc" || param_lc == "value" || param_lc == "i") =>
            {
                return dc_level(waveform);
            }
            Element::XOsdi { name, params, .. }
            | Element::Mosfet { name, params, .. }
            | Element::Diode { name, params, .. }
            | Element::Bjt { name, params, .. }
                if name.to_lowercase() == el_name =>
            {
                return params
                    .iter()
                    .find(|(k, _)| k.to_lowercase() == param_lc)
                    .map(|(_, v)| *v);
            }
            _ => {}
        }
    }
    None
}

/// The DC level of a waveform.  `set_element_param` replaces the whole
/// waveform with a `Dc`, so that is the only shape that round-trips; a
/// time-varying source has no single value this getter could honestly return.
fn dc_level(waveform: &Waveform) -> Option<f64> {
    match waveform {
        Waveform::Dc(v) => Some(*v),
        _ => None,
    }
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
        parse_spice("* t\nV1 in 0 DC 1\nR1 in out 1k\nC1 out 0 1p\n.op\n").unwrap()
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

    /// The gap this function exists to close: a Verilog-A transistor arrives on
    /// an `M` line via a `.model` card, so reaching only `X` elements misses the
    /// whole PDK idiom. The CLI's private copy used to do exactly that.
    #[test]
    fn mosfet_and_osdi_instance_params_are_reachable() {
        let mut n = parse_spice(
            "* m-line\n.model nch NMOS (VTO=0.7 KP=100u)\n\
             M1 d g 0 0 nch W=10u L=1u\nVg g 0 1\nVd d 0 2\n.op\n",
        )
        .unwrap();
        assert!(set_element_param(&mut n, "M1", "W", 40e-6));
        // Existing key is updated, not duplicated — and casing does not matter.
        assert!(set_element_param(&mut n, "m1", "w", 20e-6));
        match n
            .elements
            .iter()
            .find(|e| matches!(e, Element::Mosfet { .. }))
        {
            Some(Element::Mosfet { params, .. }) => {
                let w: Vec<_> = params.iter().filter(|(k, _)| k == "w").collect();
                assert_eq!(w.len(), 1, "duplicate w entries: {params:?}");
                assert_eq!(w[0].1, 20e-6);
            }
            other => panic!("expected M1, got {other:?}"),
        }
    }

    /// Diode and BJT instance params reach the device as of #26, so a sweep or
    /// an optimiser can drive `AREA` — and the edit is claimed only because
    /// something now consumes it.
    #[test]
    fn diode_and_bjt_instance_params_are_reachable() {
        let mut n = parse_spice(
            "* d/q\n.model dm D (IS=1e-14)\n.model qm NPN (IS=1e-16)\n\
             D1 a 0 dm\nQ1 c b 0 qm\nV1 a 0 1\nVb b 0 0.7\nVc c 0 2\n.op\n",
        )
        .unwrap();
        assert!(set_element_param(&mut n, "d1", "area", 2.0));
        assert!(set_element_param(&mut n, "q1", "area", 3.0));
        assert_eq!(get_element_param(&n, "d1", "area"), Some(2.0));
        assert_eq!(get_element_param(&n, "q1", "area"), Some(3.0));
    }

    /// The bare element letter, which the CLI accepted and the shared version
    /// did not — keeping it is why the CLI could adopt this without regressing.
    #[test]
    fn bare_element_letter_aliases_work() {
        let mut n = nl();
        assert!(set_element_param(&mut n, "r1", "r", 5e3));
        assert!(set_element_param(&mut n, "c1", "c", 5e-12));
        assert!(!set_element_param(&mut n, "r1", "c", 1.0), "wrong letter");
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
