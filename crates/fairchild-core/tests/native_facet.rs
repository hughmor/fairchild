//! `fc_facet` — the one-port terminator / partial reflector / mirror.
//!
//! The bundles are written out wire by wire rather than through
//! `.optical_port`, and driven by plain `V` sources rather than a laser. That
//! is deliberate: `fc_cw_laser` *drives* its port's backward wires to zero, so
//! putting one at the far end of a reflecting chain gives two devices pinning
//! the same node to different values. Any real deck that wants to watch a
//! reflection needs the same care — a laser is a perfect absorber only as long
//! as nothing else is trying to talk on that wire.
//!
//! Wire order under `enable_bidirectional=1` is `[re_fw, im_fw, re_bw, im_bw, λ]`.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

/// Field launched into the port, √W. 1 V of `re` is 1 W, so this is 1 mW.
const A_IN: f64 = 0.031_622_776_601_683_79;

fn run(body: &str) -> fairchild_core::newton::NrResult {
    let src = format!(".options enable_bidirectional=1\n{body}.op\n.end\n");
    let net = parse_spice(&src).expect("netlist should parse");
    dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP should converge")
}

/// Drive one facet directly and read the reflected field back off the same port.
fn reflect(params: &str, drive_re: f64, drive_im: f64) -> (f64, f64) {
    let r = run(&format!(
        "Xf p_re_fw p_im_fw p_re_bw p_im_bw p_wl fc_facet {params}\n\
         V_re p_re_fw 0 DC {drive_re}\n\
         V_im p_im_fw 0 DC {drive_im}\n\
         V_wl p_wl 0 DC 1.55e-6\n"
    ));
    (
        r.node_voltage("p_re_bw").unwrap(),
        r.node_voltage("p_im_bw").unwrap(),
    )
}

/// A facet with no parameters absorbs everything — which is what a terminator
/// is, and the reason it is the default rather than a value someone has to
/// remember to write.
#[test]
fn the_default_facet_is_a_terminator() {
    let (re, im) = reflect("", A_IN, 0.0);
    assert_eq!((re, im), (0.0, 0.0));
}

/// `reflectance` is a POWER fraction, so the returned field carries `√R`.
/// Getting that wrong is a factor-of-R error that looks plausible at every
/// value except 0 and 1.
#[test]
fn reflectance_is_a_power_fraction() {
    for r_power in [0.04, 0.3, 1.0] {
        let (re, im) = reflect(&format!("reflectance={r_power}"), A_IN, 0.0);
        let p_out = re * re + im * im;
        let expect = r_power * A_IN * A_IN;
        assert!(
            (p_out - expect).abs() < 1e-15,
            "R={r_power}: reflected {p_out:.6e} W, expected {expect:.6e} W"
        );
    }
}

/// Reflection phase lands on the field as `e^(−jφ)`, matching the sign
/// convention `OpticalSegment` uses for propagation — so a facet and a length
/// of waveguide compose without a hidden conjugate between them.
#[test]
fn the_reflection_phase_rotates_the_field() {
    for (phi, want_re, want_im) in [
        (0.0, 1.0, 0.0),
        (180.0, -1.0, 0.0),
        (90.0, 0.0, -1.0),
        (-90.0, 0.0, 1.0),
    ] {
        let (re, im) = reflect(&format!("reflectance=1 phase_deg={phi}"), 1.0, 0.0);
        assert!(
            (re - want_re).abs() < 1e-12 && (im - want_im).abs() < 1e-12,
            "φ={phi}°: got ({re:.6}, {im:.6}), expected ({want_re}, {want_im})"
        );
    }
}

/// Leave one of the three out and it takes the remainder, so a mirror can be
/// written `loss=0 transmittance=0` or `reflectance=1` and mean the same thing.
#[test]
fn the_budget_infers_whichever_term_is_left_out() {
    let (re, im) = reflect("transmittance=0.25 loss=0.15", A_IN, 0.0);
    let p_out = re * re + im * im;
    let expect = 0.6 * A_IN * A_IN;
    assert!(
        (p_out - expect).abs() < 1e-15,
        "{p_out:.6e} vs {expect:.6e}"
    );
}

/// Over-unity is a typo, not a design. A facet that silently renormalised
/// would hide exactly the mistake the three-parameter form exists to catch.
#[test]
#[should_panic(expected = "no power left")]
fn an_over_unity_budget_is_rejected() {
    reflect("reflectance=0.9 transmittance=0.5", A_IN, 0.0);
}

#[test]
#[should_panic(expected = "must be 1")]
fn a_fully_specified_budget_must_sum_to_one() {
    reflect("reflectance=0.5 transmittance=0.2 loss=0.2", A_IN, 0.0);
}

/// Unidirectional propagation has no backward wire, so a reflector has nowhere
/// to put the light. Failing loudly beats being a terminator that the deck
/// believes is a mirror.
#[test]
#[should_panic(expected = "enable_bidirectional")]
fn a_reflector_needs_bidirectional_propagation() {
    let net = parse_spice(
        "Xf p_re p_im p_wl fc_facet reflectance=0.5\n\
         V_re p_re 0 DC 1\nV_im p_im 0 DC 0\nV_wl p_wl 0 DC 1.55e-6\n.op\n.end\n",
    )
    .unwrap();
    let _ = dc_op_nr_with_registry(&net, &DeviceRegistry::new());
}

/// A terminator is still legal without bidirectional propagation — there is
/// simply nothing for it to do, and saying so is better than a special case in
/// every deck that wants to be runnable both ways.
#[test]
fn a_terminator_is_legal_without_bidirectional_propagation() {
    let net = parse_spice(
        "Xf p_re p_im p_wl fc_facet\n\
         V_re p_re 0 DC 1\nV_im p_im 0 DC 0\nV_wl p_wl 0 DC 1.55e-6\n.op\n.end\n",
    )
    .unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    assert_eq!(r.node_voltage("p_re").unwrap(), 1.0);
}

/// The round trip through a lossy waveguide: the light pays the propagation
/// loss twice and the reflectance once.
///
/// `alpha_dB_cm` is a POWER loss, so 1000 µm at 2 dB/cm is 0.2 dB each way.
#[test]
fn a_reflected_wave_pays_the_propagation_loss_both_ways() {
    const R_POWER: f64 = 0.3;
    let r = run(&format!(
        "Xwg s_re_fw s_im_fw s_re_bw s_im_bw s_wl f_re_fw f_im_fw f_re_bw f_im_bw f_wl \
           fc_waveguide L_um=1000 n_eff=2.445 n_g=4.2 alpha_dB_cm=2.0\n\
         Xf  f_re_fw f_im_fw f_re_bw f_im_bw f_wl fc_facet reflectance={R_POWER}\n\
         V_re s_re_fw 0 DC {A_IN}\n\
         V_im s_im_fw 0 DC 0\n\
         V_wl s_wl 0 DC 1.55e-6\n"
    ));
    let re = r.node_voltage("s_re_bw").unwrap();
    let im = r.node_voltage("s_im_bw").unwrap();
    let one_way = 10f64.powf(-2.0 * 0.1 / 10.0); // 0.1 cm at 2 dB/cm, power
    let expect = A_IN * A_IN * one_way * R_POWER * one_way;
    let got = re * re + im * im;
    assert!(
        (got - expect).abs() / expect < 1e-9,
        "returned {got:.6e} W, expected {expect:.6e} W"
    );
}
