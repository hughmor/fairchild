/// Integration tests for the SPICE parser.
///
/// These tests verify the public API of `fairchild-parser` in isolation —
/// no simulator, no solver.  They serve as the safety net for parser
/// refactors: if all tests here pass before AND after a file split, the
/// split is behaviour-neutral.
use fairchild_parser::{parse_spice, AcVariation, Analysis, Element, Waveform};

// ── Helpers ───────────────────────────────────────────────────────────────

fn parse_ok(s: &str) -> fairchild_parser::Netlist {
    parse_spice(s).unwrap_or_else(|e| panic!("parse failed: {e}"))
}

// ── Element parsing ────────────────────────────────────────────────────────

#[test]
fn resistor_parses() {
    let nl = parse_ok("R1 a b 1k\n.op\n.end\n");
    assert_eq!(nl.elements.len(), 1);
    if let Element::Resistor {
        name,
        pos,
        neg,
        resistance,
    } = &nl.elements[0]
    {
        assert_eq!(name, "r1");
        assert_eq!(pos, "a");
        assert_eq!(neg, "b");
        assert!((resistance - 1000.0).abs() < 1e-9);
    } else {
        panic!("expected Resistor");
    }
}

#[test]
fn capacitor_parses() {
    let nl = parse_ok("C1 a b 1u\n.op\n.end\n");
    if let Element::Capacitor {
        name, capacitance, ..
    } = &nl.elements[0]
    {
        assert_eq!(name, "c1");
        assert!((capacitance - 1e-6).abs() < 1e-18);
    } else {
        panic!("expected Capacitor");
    }
}

#[test]
fn inductor_parses() {
    let nl = parse_ok("L1 a b 1m\n.op\n.end\n");
    if let Element::Inductor {
        name, inductance, ..
    } = &nl.elements[0]
    {
        assert_eq!(name, "l1");
        assert!((inductance - 1e-3).abs() < 1e-15);
    } else {
        panic!("expected Inductor");
    }
}

#[test]
fn voltage_source_dc_parses() {
    let nl = parse_ok("V1 a 0 DC 5\n.op\n.end\n");
    if let Element::VoltageSource {
        name,
        pos,
        neg,
        waveform: Waveform::Dc(v),
        ..
    } = &nl.elements[0]
    {
        assert_eq!(name, "v1");
        assert_eq!(pos, "a");
        assert_eq!(neg, "0");
        assert!((v - 5.0).abs() < 1e-12);
    } else {
        panic!("expected DC VoltageSource");
    }
}

#[test]
fn current_source_dc_parses() {
    let nl = parse_ok("I1 a 0 DC 1m\n.op\n.end\n");
    if let Element::CurrentSource {
        name,
        waveform: Waveform::Dc(i),
        ..
    } = &nl.elements[0]
    {
        assert_eq!(name, "i1");
        assert!((i - 1e-3).abs() < 1e-15);
    } else {
        panic!("expected DC CurrentSource");
    }
}

#[test]
fn diode_parses() {
    let nl = parse_ok("D1 a b myD\n.model myD D\n.op\n.end\n");
    if let Element::Diode {
        name,
        anode,
        cathode,
        model_name,
        ..
    } = &nl.elements[0]
    {
        assert_eq!(name, "d1");
        assert_eq!(anode, "a");
        assert_eq!(cathode, "b");
        assert_eq!(model_name, "myd");
    } else {
        panic!("expected Diode");
    }
}

#[test]
fn nmos_parses() {
    let nl = parse_ok("M1 d g s b nmos\n.model nmos NMOS\n.op\n.end\n");
    if let Element::Mosfet {
        name,
        drain,
        gate,
        source,
        model_name,
        ..
    } = &nl.elements[0]
    {
        assert_eq!(name, "m1");
        assert_eq!(drain, "d");
        assert_eq!(gate, "g");
        assert_eq!(source, "s");
        assert_eq!(model_name, "nmos");
    } else {
        panic!("expected Mosfet");
    }
}

#[test]
fn npn_bjt_parses() {
    let nl = parse_ok("Q1 c b e npn1\n.model npn1 NPN\n.op\n.end\n");
    if let Element::Bjt {
        name,
        collector,
        base,
        emitter,
        model_name,
        ..
    } = &nl.elements[0]
    {
        assert_eq!(name, "q1");
        assert_eq!(collector, "c");
        assert_eq!(base, "b");
        assert_eq!(emitter, "e");
        assert_eq!(model_name, "npn1");
    } else {
        panic!("expected BJT");
    }
}

#[test]
fn coupled_inductor_parses() {
    let nl = parse_ok("L1 a b 1m\nL2 c d 1m\nK1 L1 L2 0.9\n.op\n.end\n");
    let k = nl
        .elements
        .iter()
        .find(|e| matches!(e, Element::CoupledInductors { .. }));
    assert!(k.is_some(), "CoupledInductors element not found");
    if let Some(Element::CoupledInductors {
        name,
        l1,
        l2,
        coupling,
    }) = k
    {
        assert_eq!(name, "k1");
        assert_eq!(l1, "l1");
        assert_eq!(l2, "l2");
        assert!((coupling - 0.9).abs() < 1e-12);
    }
}

#[test]
fn behavioral_source_parses() {
    let nl = parse_ok("V1 in 0 DC 1\nR1 in out 1k\nB1 a 0 V=V(in)*2\n.op\n.end\n");
    assert!(nl
        .elements
        .iter()
        .any(|e| matches!(e, Element::Behavioral { .. })));
}

#[test]
fn xosdi_element_parses() {
    // X elements without a matching .subckt definition parse to XOsdi
    let nl = parse_ok("X1 in out mymod\n.op\n.end\n");
    assert!(nl
        .elements
        .iter()
        .any(|e| matches!(e, Element::XOsdi { .. })));
}

// ── Value suffix parsing ───────────────────────────────────────────────────

#[test]
fn suffix_k_meg_m_u_n_p_f() {
    for (s, expected) in [
        ("1k", 1e3),
        ("1K", 1e3),
        ("1meg", 1e6),
        ("1MEG", 1e6),
        ("1m", 1e-3),
        ("1M", 1e-3),
        ("1u", 1e-6),
        ("1U", 1e-6),
        ("1n", 1e-9),
        ("1N", 1e-9),
        ("1p", 1e-12),
        ("1P", 1e-12),
        ("1f", 1e-15),
        ("1F", 1e-15),
    ] {
        let nl = parse_ok(&format!("R1 a b {s}\n.op\n.end\n"));
        if let Element::Resistor { resistance, .. } = &nl.elements[0] {
            assert!(
                (resistance - expected).abs() / expected < 1e-9,
                "suffix '{s}': got {resistance}, want {expected}"
            );
        }
    }
}

// ── Analysis directive parsing ─────────────────────────────────────────────

#[test]
fn op_analysis() {
    let nl = parse_ok("R1 a 0 1k\n.op\n.end\n");
    assert!(nl.analyses.iter().any(|a| matches!(a, Analysis::Op)));
}

#[test]
fn tran_analysis() {
    let nl = parse_ok("R1 a 0 1k\n.tran 1n 10u\n.end\n");
    if let Some(Analysis::Tran {
        step, stop, uic, ..
    }) = nl.analyses.first()
    {
        assert!((step - 1e-9).abs() < 1e-21);
        assert!((stop - 10e-6).abs() < 1e-18);
        assert!(!uic);
    } else {
        panic!("expected Tran");
    }
}

#[test]
fn dc_analysis() {
    let nl = parse_ok("V1 a 0 DC 0\n.dc V1 0 5 0.1\n.end\n");
    if let Some(Analysis::Dc {
        src,
        start,
        stop,
        step,
        ..
    }) = nl.analyses.first()
    {
        assert_eq!(src.to_lowercase(), "v1");
        assert!((start - 0.0).abs() < 1e-12);
        assert!((stop - 5.0).abs() < 1e-12);
        assert!((step - 0.1).abs() < 1e-12);
    } else {
        panic!("expected Dc");
    }
}

#[test]
fn ac_analysis_dec() {
    let nl = parse_ok("R1 a 0 1k\n.ac DEC 10 1k 1G\n.end\n");
    if let Some(Analysis::Ac {
        variation,
        points,
        fstart,
        fstop,
    }) = nl.analyses.first()
    {
        assert!(matches!(variation, AcVariation::Dec));
        assert_eq!(*points, 10);
        assert!((fstart - 1e3).abs() < 1.0);
        assert!((fstop - 1e9).abs() < 1e3);
    } else {
        panic!("expected Ac");
    }
}

#[test]
fn noise_analysis() {
    let nl = parse_ok("R1 a 0 1k\n.noise V(a) V1 DEC 10 1 1G\n.end\n");
    assert!(nl
        .analyses
        .iter()
        .any(|a| matches!(a, Analysis::Noise { .. })));
}

// ── Waveform parsing ───────────────────────────────────────────────────────

#[test]
fn pulse_waveform() {
    let nl = parse_ok("V1 a 0 PULSE(0 1 10n 1n 1n 50n 100n)\n.op\n.end\n");
    if let Element::VoltageSource {
        waveform:
            Waveform::Pulse {
                v0,
                v1,
                td,
                tr,
                tf,
                pw,
                per,
            },
        ..
    } = &nl.elements[0]
    {
        assert!((v0 - 0.0).abs() < 1e-12);
        assert!((v1 - 1.0).abs() < 1e-12);
        assert!((td - 10e-9).abs() < 1e-21);
        assert!((tr - 1e-9).abs() < 1e-21);
        assert!((tf - 1e-9).abs() < 1e-21);
        assert!((pw - 50e-9).abs() < 1e-21);
        assert!((per - 100e-9).abs() < 1e-21);
    } else {
        panic!("expected PULSE waveform");
    }
}

#[test]
fn sin_waveform() {
    let nl = parse_ok("V1 a 0 SIN(0 1 1k)\n.op\n.end\n");
    if let Element::VoltageSource {
        waveform: Waveform::Sin { vo, va, freq, .. },
        ..
    } = &nl.elements[0]
    {
        assert!((vo - 0.0).abs() < 1e-12);
        assert!((va - 1.0).abs() < 1e-12);
        assert!((freq - 1e3).abs() < 0.1);
    } else {
        panic!("expected SIN waveform");
    }
}

#[test]
fn exp_waveform() {
    let nl = parse_ok("V1 a 0 EXP(0 1 0 1u 1u 2u)\n.op\n.end\n");
    assert!(matches!(
        &nl.elements[0],
        Element::VoltageSource {
            waveform: Waveform::Exp { .. },
            ..
        }
    ));
}

#[test]
fn pwl_waveform() {
    let nl = parse_ok("V1 a 0 PWL(0 0 1u 1 2u 0)\n.op\n.end\n");
    if let Element::VoltageSource {
        waveform: Waveform::Pwl { points },
        ..
    } = &nl.elements[0]
    {
        assert_eq!(points.len(), 3);
        assert!((points[0].0 - 0.0).abs() < 1e-12);
        assert!((points[1].1 - 1.0).abs() < 1e-12);
    } else {
        panic!("expected PWL waveform");
    }
}

// ── Subcircuit expansion ───────────────────────────────────────────────────

#[test]
fn subckt_expands_correctly() {
    let nl = parse_ok(
        ".subckt inv in out vdd\n\
         M1 out in vdd vdd pmos W=1u L=0.35u\n\
         M2 out in 0   0   nmos W=1u L=0.35u\n\
         .ends\n\
         X1 a b vcc inv\n\
         .op\n.end\n",
    );
    // Two MOSFETs expanded from subcircuit
    let mosfets: Vec<_> = nl
        .elements
        .iter()
        .filter(|e| matches!(e, Element::Mosfet { .. }))
        .collect();
    assert_eq!(
        mosfets.len(),
        2,
        "expected 2 MOSFETs after subckt expansion"
    );
}

#[test]
fn param_substitution_in_subckt() {
    let nl = parse_ok(
        ".subckt rpad in out\n\
         .param RVAL=1k\n\
         R1 in out {RVAL}\n\
         .ends\n\
         X1 a b rpad\n\
         .op\n.end\n",
    );
    let r = nl
        .elements
        .iter()
        .find(|e| matches!(e, Element::Resistor { .. }));
    assert!(r.is_some(), "Resistor from subckt not found");
    if let Some(Element::Resistor { resistance, .. }) = r {
        assert!((resistance - 1000.0).abs() < 1e-9);
    }
}

// ── .param directive ───────────────────────────────────────────────────────

#[test]
fn global_param_resolves_in_element() {
    let nl = parse_ok(
        ".param RVAL=2.2k\n\
         R1 a b {RVAL}\n\
         .op\n.end\n",
    );
    if let Element::Resistor { resistance, .. } = &nl.elements[0] {
        assert!((resistance - 2200.0).abs() < 1e-9);
    } else {
        panic!("expected Resistor");
    }
}

// ── .options ──────────────────────────────────────────────────────────────

#[test]
fn options_parses_key_value_pairs() {
    let nl = parse_ok("R1 a 0 1k\n.options reltol=1e-4 abstol=1e-12\n.op\n.end\n");
    assert!(!nl.options.is_empty(), "options map should be non-empty");
    let reltol = nl.options.iter().find(|(k, _)| k == "reltol");
    assert!(reltol.is_some());
}

// ── .ic / .nodeset ─────────────────────────────────────────────────────────

#[test]
fn ic_directive_parses() {
    let nl = parse_ok("C1 a 0 1u\n.ic V(a)=3.3\n.tran 1n 10n\n.end\n");
    assert!(!nl.ic.is_empty(), ".ic should set initial conditions");
}

// ── .measure ──────────────────────────────────────────────────────────────

#[test]
fn measure_find_parses() {
    let nl = parse_ok("R1 a 0 1k\n.tran 1n 10n\n.meas tran vmax FIND V(a) AT=5n\n.end\n");
    assert!(
        !nl.measurements.is_empty(),
        ".meas should add a measurement"
    );
}

// ── Comment and continuation lines ────────────────────────────────────────

#[test]
fn star_comment_ignored() {
    let nl = parse_ok("* this whole line is a comment\nR1 a b 1k\n.op\n.end\n");
    assert_eq!(nl.elements.len(), 1);
}

#[test]
fn plus_continuation_line_joins() {
    // R1 split across two lines via '+'
    let nl = parse_ok("R1 a b\n+ 100\n.op\n.end\n");
    assert_eq!(nl.elements.len(), 1);
    if let Element::Resistor { resistance, .. } = &nl.elements[0] {
        assert!((resistance - 100.0).abs() < 1e-9);
    } else {
        panic!("expected Resistor");
    }
}

// ── Case insensitivity ────────────────────────────────────────────────────

#[test]
fn node_names_case_insensitive() {
    let nl = parse_ok("R1 IN OUT 1k\nR2 in out 2k\n.op\n.end\n");
    // Both should reference the same nodes
    assert_eq!(nl.elements.len(), 2);
    if let (Element::Resistor { pos: p1, .. }, Element::Resistor { pos: p2, .. }) =
        (&nl.elements[0], &nl.elements[1])
    {
        assert_eq!(p1, p2, "node names should be normalized to lowercase");
    }
}

// ── Parasitic expansion ────────────────────────────────────────────────────

#[test]
fn inductor_rser_expands_to_two_elements() {
    let nl = parse_ok("L1 a b 1m rser=10\n.op\n.end\n");
    // Should expand to L + R → 2 elements
    assert_eq!(nl.elements.len(), 2);
    assert!(nl
        .elements
        .iter()
        .any(|e| matches!(e, Element::Inductor { .. })));
    assert!(nl
        .elements
        .iter()
        .any(|e| matches!(e, Element::Resistor { .. })));
}

#[test]
fn capacitor_esr_expands_to_two_elements() {
    let nl = parse_ok("C1 a b 1u esr=5\n.op\n.end\n");
    assert_eq!(nl.elements.len(), 2);
    assert!(nl
        .elements
        .iter()
        .any(|e| matches!(e, Element::Capacitor { .. })));
    assert!(nl
        .elements
        .iter()
        .any(|e| matches!(e, Element::Resistor { .. })));
}

#[test]
fn resistor_cpar_expands_to_two_elements() {
    let nl = parse_ok("R1 a b 1k cpar=1p\n.op\n.end\n");
    assert_eq!(nl.elements.len(), 2);
    assert!(nl
        .elements
        .iter()
        .any(|e| matches!(e, Element::Resistor { .. })));
    assert!(nl
        .elements
        .iter()
        .any(|e| matches!(e, Element::Capacitor { .. })));
}

// ── .model card ───────────────────────────────────────────────────────────

#[test]
fn model_card_nmos_parses() {
    let nl = parse_ok(".model nch NMOS VTO=0.5 KP=100u LAMBDA=0.01\n.op\n.end\n");
    assert_eq!(nl.models.len(), 1);
    let m = &nl.models[0];
    assert_eq!(m.kind, "nmos");
    let vto = m
        .params
        .iter()
        .find(|(k, _)| k == "vto")
        .map(|(_, v)| *v)
        .unwrap();
    assert!((vto - 0.5).abs() < 1e-12);
}

#[test]
fn model_card_pnp_parses() {
    let nl = parse_ok(".model pnp1 PNP IS=1e-15 BF=100\n.op\n.end\n");
    assert_eq!(nl.models.len(), 1);
    assert_eq!(nl.models[0].kind, "pnp");
}

// ── GND aliases ───────────────────────────────────────────────────────────

#[test]
fn gnd_aliases_normalised() {
    for alias in ["gnd", "GND", "0"] {
        let nl = parse_ok(&format!("R1 a {alias} 1k\n.op\n.end\n"));
        if let Element::Resistor { neg, .. } = &nl.elements[0] {
            assert_eq!(neg, "0", "alias '{alias}' should normalise to '0'");
        }
    }
}

// ── Transmission line (T) ──────────────────────────────────────────────────

#[test]
fn parses_transmission_line_z0_td() {
    let nl = parse_ok("T1 in 0 out 0 Z0=50 TD=2n\n.op\n.end\n");
    if let Element::TransmissionLine {
        name,
        a_pos,
        a_neg,
        b_pos,
        b_neg,
        z0,
        td,
    } = &nl.elements[0]
    {
        assert_eq!(name, "t1");
        assert_eq!(a_pos, "in");
        assert_eq!(a_neg, "0");
        assert_eq!(b_pos, "out");
        assert_eq!(b_neg, "0");
        assert!((z0 - 50.0).abs() < 1e-9);
        assert!((td - 2e-9).abs() < 1e-18);
    } else {
        panic!("expected TransmissionLine, got {:?}", nl.elements[0]);
    }
}

#[test]
fn transmission_line_delay_from_freq_and_nl() {
    // TD = NL / F = 0.25 / 1e9 = 250 ps (default NL is a quarter wavelength).
    let nl = parse_ok("T1 a 0 b 0 Z0=75 F=1g\n.op\n.end\n");
    if let Element::TransmissionLine { z0, td, .. } = &nl.elements[0] {
        assert!((z0 - 75.0).abs() < 1e-9);
        assert!((td - 0.25e-9).abs() < 1e-18, "td={td}");
    } else {
        panic!("expected TransmissionLine");
    }
}
