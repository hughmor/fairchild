//! A user's own OSDI model, carrying a WDM bundle (#52).
//!
//! Before the registry became the arity authority this was impossible: the
//! parser decided dispatch from a hardcoded list of `fc_*` names, so every
//! OSDI model fell through to `Scalar` and a multi-channel bundle on one was a
//! parse error whose only suggested fix was to patch the parser and rebuild.
//!
//! Now a descriptor's fixed terminal count places the instance by shape — if
//! flattening the referenced bundles lands on exactly `num_terminals`, the
//! model takes the whole bus. `wg_wdm2` has its 12 ports written out by hand,
//! which is what stage 1 supports; generating them per N is #55.
//!
//! The two channels carry different loss on purpose. With identical channels a
//! crossed index would be invisible, so the asymmetry is the test: each output
//! must show its OWN channel's attenuation, against a hand-computed budget
//! rather than against a second run that could share the fault.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_osdi::{load_libraries, VaOptions};
use fairchild_parser::{parse_spice_with_arity, PermissiveArity};

mod common;

/// 10 mW at 3 dB/cm and 40 mW at 12 dB/cm, through 1000 µm.
const P0_MW: f64 = 10.0;
const P1_MW: f64 = 40.0;
const A0_DB_CM: f64 = 3.0;
const A1_DB_CM: f64 = 12.0;
const L_UM: f64 = 1000.0;

fn deck() -> String {
    format!(
        "\
.optical_port c0
.optical_port c1
.optical_port bus 2
.optical_port dout 2
Xl0 c0 fc_cw_laser power_mW={P0_MW} wavelength_nm=1550
Xl1 c1 fc_cw_laser power_mW={P1_MW} wavelength_nm=1551
Xmux bus c0 c1 fc_mux
Xw bus dout wg_wdm2 l_um={L_UM} alpha_dB_cm_0={A0_DB_CM} alpha_dB_cm_1={A1_DB_CM}
+ wl_0_nm=1550 wl_1_nm=1551
.op
"
    )
}

/// Power (W) surviving `alpha_dB_cm` over `L_UM`, by hand: dB is POWER dB.
fn expect_mw(p_mw: f64, alpha_db_cm: f64) -> f64 {
    p_mw * 10f64.powf(-alpha_db_cm * (L_UM * 1e-4) / 10.0)
}

#[test]
fn an_osdi_model_carries_a_wdm_bundle_and_keeps_its_channels_apart() {
    let Some(osdi) = common::compiled("wg_wdm2") else {
        return; // no Verilog-A compiler; skips like every test in this directory
    };

    let build = || {
        let mut reg = DeviceRegistry::new();
        load_libraries(
            &[osdi.to_string_lossy().to_string()],
            &[],
            None,
            &VaOptions::from_env(),
            &mut reg,
        )
        .expect("the artefact just compiled must load");
        reg
    };

    // The load boundary's two passes: harvest, then parse knowing what each
    // name resolves to. The first pass is what the parser could never do alone.
    let src = deck();
    let probe = parse_spice_with_arity(&src, &PermissiveArity).expect("probe pass parses");
    let reg = build();
    let _ = probe;
    let netlist = parse_spice_with_arity(&src, &reg)
        .expect("a 12-terminal model should take a 2-channel bundle");

    let r = dc_op_nr_with_registry(&netlist, &reg).expect("DC OP converges");
    let power = |ch: usize| {
        let re = r.node_voltage(&format!("dout_re_{ch}")).unwrap();
        let im = r.node_voltage(&format!("dout_im_{ch}")).unwrap();
        (re * re + im * im) * 1e3
    };

    let (want0, want1) = (expect_mw(P0_MW, A0_DB_CM), expect_mw(P1_MW, A1_DB_CM));
    let (got0, got1) = (power(0), power(1));

    assert!(
        (got0 - want0).abs() < 1e-6,
        "channel 0: got {got0:.6} mW, hand-computed {want0:.6} mW"
    );
    assert!(
        (got1 - want1).abs() < 1e-6,
        "channel 1: got {got1:.6} mW, hand-computed {want1:.6} mW"
    );
    // If the indices were crossed the numbers would be these instead, and both
    // asserts above would still be comparing plausible-looking powers.
    let crossed0 = expect_mw(P0_MW, A1_DB_CM);
    assert!(
        (got0 - crossed0).abs() > 1.0,
        "channel 0 is showing channel 1's loss ({got0:.4} mW ≈ {crossed0:.4} mW)"
    );
}

/// A width the model cannot serve is still refused. 3 channels flattens to 18
/// terminals and one-each gives 6; the descriptor has 12, so neither shape fits
/// and the deck must not quietly become something else.
#[test]
fn a_bundle_width_the_model_cannot_serve_is_refused() {
    let Some(osdi) = common::compiled("wg_wdm2") else {
        return;
    };
    let mut reg = DeviceRegistry::new();
    load_libraries(
        &[osdi.to_string_lossy().to_string()],
        &[],
        None,
        &VaOptions::from_env(),
        &mut reg,
    )
    .unwrap();

    let src = "\
.optical_port bus 3
.optical_port dout 3
Xw bus dout wg_wdm2
.op
";
    let err =
        parse_spice_with_arity(src, &reg).expect_err("a 12-terminal model cannot serve 3 channels");
    let msg = format!("{err}");
    assert!(
        msg.contains("no WDM semantics") || msg.contains("channel"),
        "the refusal should explain itself, got: {msg}"
    );
}
