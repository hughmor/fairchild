//! DC adjoint sensitivity, checked against closed forms and against the
//! brute-force alternative it exists to replace.
//!
//! Every gradient here is verified two ways where a closed form exists, and
//! against a full re-solve finite difference where one does not.  The re-solve
//! reference is deliberately run at a much tighter `reltol` than the default:
//! differencing a converged solution inherits that solution's convergence
//! error, which is the accuracy argument for the adjoint in the first place.
//! `finite_difference_at_default_tolerance_is_the_noisy_one` pins that claim
//! rather than leaving it as a docstring assertion.

use fairchild_core::adjoint::{dc_sensitivity, Output, ParamRef};
use fairchild_core::{dc_op_nr_with_registry_opts, DeviceRegistry, SimOptions};
use fairchild_parser::{parse_spice, Netlist};

/// Tight enough that a re-solve finite difference is trustworthy as a
/// reference; the default `reltol = 1e-3` is not.
fn tight() -> SimOptions {
    SimOptions {
        reltol: 1e-12,
        vntol: 1e-14,
        ..SimOptions::default()
    }
}

fn registry(net: &Netlist) -> DeviceRegistry {
    let mut r = DeviceRegistry::new();
    r.register_builtin_models(&net.models);
    r
}

/// Evaluate one output by actually re-solving with the parameter perturbed —
/// the O(n_params) brute force the adjoint replaces.
fn resolve_output(
    src: &str,
    element: &str,
    param: &str,
    value: f64,
    out: &Output,
    opts: &SimOptions,
) -> f64 {
    let mut net = parse_spice(src).unwrap();
    assert!(
        fairchild_core::set_element_param(&mut net, element, param, value),
        "{element}.{param} did not match any element"
    );
    let reg = registry(&net);
    let r = dc_op_nr_with_registry_opts(&net, &reg, opts).expect("DC OP");
    match out {
        Output::NodeVoltage(n) => r.node_voltage(n).unwrap(),
        Output::NodeVoltageDiff { pos, neg } => {
            r.node_voltage(pos).unwrap() - r.node_voltage(neg).unwrap()
        }
        Output::BranchCurrent(n) => r.topo.vsrc_current(n, &r.x).unwrap(),
        Output::OpticalPower { net: n, channel } => {
            let re = r.node_voltage(&format!("{n}_re_{channel}")).unwrap();
            let im = r.node_voltage(&format!("{n}_im_{channel}")).unwrap();
            re * re + im * im
        }
        Output::Custom(_) => unreachable!("not used as a re-solve reference"),
    }
}

/// Central difference of the whole simulation, at whatever tolerance `opts` sets.
fn fd_reference(
    src: &str,
    element: &str,
    param: &str,
    nominal: f64,
    rel_step: f64,
    out: &Output,
    opts: &SimOptions,
) -> f64 {
    let h = rel_step * nominal.abs();
    let plus = resolve_output(src, element, param, nominal + h, out, opts);
    let minus = resolve_output(src, element, param, nominal - h, out, opts);
    (plus - minus) / (2.0 * h)
}

fn rel_err(got: f64, want: f64) -> f64 {
    (got - want).abs() / want.abs().max(1e-30)
}

// ---------------------------------------------------------------------------
// Linear circuits — the gradient is a closed form, so there is no excuse
// ---------------------------------------------------------------------------

const DIVIDER: &str = "\
* resistive divider
V1 in 0 DC 1
R1 in out 1k
R2 out 0 2k
.op
";

#[test]
fn a_resistive_divider_matches_the_closed_form() {
    let net = parse_spice(DIVIDER).unwrap();
    let opts = tight();
    let outputs = [Output::NodeVoltage("out".into())];
    let params = [
        ParamRef::new("R1", "value"),
        ParamRef::new("R2", "value"),
        ParamRef::new("V1", "dc"),
    ];

    let s = dc_sensitivity(&net, &registry(&net), &opts, &outputs, &params).unwrap();
    assert_eq!(
        s.unreached(&params).len(),
        0,
        "every parameter must resolve"
    );
    // gmin is stamped on every node diagonal, so the OP is a divider
    // loaded by 1e-12 S — that bounds every closed-form comparison here.
    assert!((s.values[0] - 2.0 / 3.0).abs() < 1e-8);

    // v(out) = V1·R2/(R1+R2)
    let (r1, r2, v1) = (1e3, 2e3, 1.0);
    let denom = (r1 + r2) * (r1 + r2);
    for (i, want) in [
        (0usize, -v1 * r2 / denom),
        (1, v1 * r1 / denom),
        (2, r2 / (r1 + r2)),
    ] {
        assert!(
            rel_err(s.grad[0][i], want) < 1e-7,
            "d v(out)/d {} = {}, closed form {want}",
            params[i].param,
            s.grad[0][i]
        );
    }
}

#[test]
fn a_branch_current_output_differentiates_too() {
    let net = parse_spice(DIVIDER).unwrap();
    let outputs = [Output::BranchCurrent("v1".into())];
    let params = [ParamRef::new("R1", "value")];
    let s = dc_sensitivity(&net, &registry(&net), &tight(), &outputs, &params).unwrap();

    // i(V1) = −V1/(R1+R2)  →  d/dR1 = V1/(R1+R2)²
    assert!((s.values[0] - -1.0 / 3e3).abs() < 1e-10);
    assert!(rel_err(s.grad[0][0], 1.0 / 9e6) < 1e-7);
}

/// One adjoint solve per output, shared across every parameter — so two
/// outputs and three parameters is six gradients from two back-substitutions.
#[test]
fn several_outputs_and_parameters_come_out_of_one_solve() {
    let net = parse_spice(DIVIDER).unwrap();
    let outputs = [
        Output::NodeVoltage("out".into()),
        Output::BranchCurrent("v1".into()),
    ];
    let params = [
        ParamRef::new("R1", "value"),
        ParamRef::new("R2", "value"),
        ParamRef::new("V1", "dc"),
    ];
    let s = dc_sensitivity(&net, &registry(&net), &tight(), &outputs, &params).unwrap();
    assert_eq!(s.grad.len(), 2);
    assert_eq!(s.grad[0].len(), 3);
    assert!(rel_err(s.grad[1][1], 1.0 / 9e6) < 1e-7); // d i(V1)/d R2
}

// ---------------------------------------------------------------------------
// Nonlinear — the Jacobian coupling is the whole content of the test
// ---------------------------------------------------------------------------

const DIODE_CLAMP: &str = "\
* series resistor into a diode
.model dmod D (IS=1e-14 N=1.0)
V1 in 0 DC 2
R1 in mid 1k
D1 mid 0 dmod
.op
";

#[test]
fn a_nonlinear_circuit_matches_a_full_finite_difference() {
    let net = parse_spice(DIODE_CLAMP).unwrap();
    let opts = tight();
    let outputs = [Output::NodeVoltage("mid".into())];
    let params = [ParamRef::new("R1", "value"), ParamRef::new("V1", "dc")];
    let s = dc_sensitivity(&net, &registry(&net), &opts, &outputs, &params).unwrap();

    for (i, (el, pn, nominal)) in [("R1", "value", 1e3), ("V1", "dc", 2.0)].iter().enumerate() {
        let want = fd_reference(DIODE_CLAMP, el, pn, *nominal, 1e-6, &outputs[0], &opts);
        assert!(
            rel_err(s.grad[0][i], want) < 1e-5,
            "d v(mid)/d {el}.{pn}: adjoint {}, re-solve FD {want}",
            s.grad[0][i]
        );
    }

    // Cross-check against the small-signal identity.  Differentiating the node
    // equation (V1 − v)/R1 = I_d(v) with respect to R1 gives
    //     dv/dR1 = −I_d / (1 + R1·g_d),   g_d = dI_d/dv = I_d/(N·V_t).
    let v_mid = s.values[0];
    let i_d = (2.0 - v_mid) / 1e3;
    let g_d = i_d / 0.025851999786445024; // I/(N·Vt) at the default temperature
    let want = -i_d / (1.0 + 1e3 * g_d);
    assert!(
        rel_err(s.grad[0][0], want) < 1e-3,
        "d v(mid)/d R1 = {}, small-signal identity {want}",
        s.grad[0][0]
    );
}

/// The accuracy claim, made falsifiable: differencing a converged solve at the
/// default `reltol = 1e-3` is visibly noisy, and the adjoint is not.
#[test]
fn finite_difference_at_default_tolerance_is_the_noisy_one() {
    let net = parse_spice(DIODE_CLAMP).unwrap();
    let out = Output::NodeVoltage("mid".into());
    let params = [ParamRef::new("R1", "value")];

    // Truth, from a re-solve at a tolerance nobody runs production at.
    let truth = fd_reference(DIODE_CLAMP, "R1", "value", 1e3, 1e-6, &out, &tight());

    let adjoint = dc_sensitivity(
        &net,
        &registry(&net),
        &SimOptions::default(),
        std::slice::from_ref(&out),
        &params,
    )
    .unwrap()
    .grad[0][0];

    let brute = fd_reference(
        DIODE_CLAMP,
        "R1",
        "value",
        1e3,
        1e-6,
        &out,
        &SimOptions::default(),
    );

    let e_adjoint = rel_err(adjoint, truth);
    let e_brute = rel_err(brute, truth);
    eprintln!("adjoint rel err {e_adjoint:.3e}, re-solve FD rel err {e_brute:.3e}");
    assert!(
        e_adjoint < 1e-6,
        "the adjoint runs at default tolerance and should still be sharp: {e_adjoint:e}"
    );
    assert!(
        e_adjoint <= e_brute,
        "the adjoint ({e_adjoint:e}) is supposed to beat the brute force ({e_brute:e})"
    );
}

// ---------------------------------------------------------------------------
// Photonic — the point of the exercise
// ---------------------------------------------------------------------------

/// A modulator biased at quadrature.  Both routes to a parameter are used:
/// `Vsig` is stamped straight from the netlist, `V_pi` lives on the device and
/// is reached through `Device::set_real_param`.
const MZM: &str = "\
.optical_port in0
.optical_port out0
Xl0 in0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xmzm in0 out0 vsig 0 fc_mzm V_pi=3.0 alpha=1.0 e_r=1000
Vsig vsig 0 DC 1.5
.op
";

#[test]
fn mzm_slope_efficiency_matches_the_closed_form() {
    let net = parse_spice(MZM).unwrap();
    let opts = tight();
    let outputs = [Output::OpticalPower {
        net: "out0".into(),
        channel: 0,
    }];
    let params = [
        ParamRef::new("Vsig", "dc"),
        ParamRef::with_nominal("Xmzm", "V_pi", 3.0),
    ];
    let s = dc_sensitivity(&net, &registry(&net), &opts, &outputs, &params).unwrap();
    assert_eq!(
        s.unreached(&params).len(),
        0,
        "V_pi must reach the device through set_real_param"
    );

    // P(V) = P_in·α·[(1 − 1/E_r)·(1 + cos(πV/V_pi))/2 + 1/E_r]
    let (p_in, alpha, er, v_pi, v) = (1e-3, 1.0, 1000.0, 3.0, 1.5);
    let k = p_in * alpha * (1.0 - 1.0 / er);
    let arg = std::f64::consts::PI * v / v_pi;
    let d_dv = -k * std::f64::consts::PI * arg.sin() / (2.0 * v_pi);
    let d_dvpi = k * std::f64::consts::PI * v * arg.sin() / (2.0 * v_pi * v_pi);

    assert!(
        rel_err(
            s.values[0],
            p_in * alpha * ((1.0 - 1.0 / er) * (1.0 + arg.cos()) / 2.0 + 1.0 / er)
        ) < 1e-9
    );
    assert!(
        rel_err(s.grad[0][0], d_dv) < 1e-6,
        "dP/dVsig = {} W/V, closed form {d_dv}",
        s.grad[0][0]
    );
    assert!(
        rel_err(s.grad[0][1], d_dvpi) < 1e-6,
        "dP/dV_pi = {} W/V, closed form {d_dvpi}",
        s.grad[0][1]
    );
}

/// Waveguide length moves the output phase, not the output power, on a
/// single-arm link — so the power gradient is the loss term alone, and it has a
/// closed form: P = P_in·10^(−α·L/1e5), dP/dL = −P·ln(10)·α/1e5 with L in µm
/// and α in dB/cm.
#[test]
fn waveguide_length_sensitivity_is_the_loss_slope() {
    let src = "\
.optical_port in0
.optical_port out0
Xl0 in0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xwg in0 out0 fc_waveguide L_um=250 n_g=4.2 alpha_dB_cm=2.0
.op
";
    let net = parse_spice(src).unwrap();
    let outputs = [Output::OpticalPower {
        net: "out0".into(),
        channel: 0,
    }];
    let params = [ParamRef::with_nominal("Xwg", "L_um", 250.0)];
    let s = dc_sensitivity(&net, &registry(&net), &tight(), &outputs, &params).unwrap();
    assert_eq!(s.unreached(&params).len(), 0);

    let want = -s.values[0] * std::f64::consts::LN_10 * 2.0 / 1e5;
    assert!(
        rel_err(s.grad[0][0], want) < 1e-5,
        "dP/dL_um = {} W/µm, closed form {want} (P = {})",
        s.grad[0][0],
        s.values[0]
    );
}

/// A voltage-controlled weight bank, which is the cleanest closed form for an
/// electro-optic gradient there is: the weight mode is defined so that
/// `P_drop − P_thru = w·P_in` with `w = w0 + dw_dv·V`, hence
/// `P_drop = P_in·(1 + w)/2` and `dP_drop/dV = P_in·dw_dv/2` exactly.
///
/// This is the test that fails if `fc_optical_2x2` stops declaring its control
/// column — the gradient goes to zero and nothing else notices.
#[test]
fn a_voltage_controlled_weight_has_an_exact_gradient() {
    let src = ".optical_port bus\n.optical_port dark\n.optical_port thru\n.optical_port drop\n\
         Xl1 bus fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
         Xwb bus dark thru drop wctl 0 fc_optical_2x2 w=0 dw_dv_0=0.5\n\
         Vw wctl 0 DC 0.4\n.op\n";
    let net = parse_spice(src).unwrap();
    let outputs = [Output::OpticalPower {
        net: "drop".into(),
        channel: 0,
    }];
    let params = [ParamRef::new("Vw", "dc")];
    let s = dc_sensitivity(&net, &registry(&net), &tight(), &outputs, &params).unwrap();
    assert_eq!(s.unreached(&params).len(), 0);

    let (p_in, dw_dv) = (1e-3, 0.5);
    assert!(
        rel_err(s.grad[0][0], p_in * dw_dv / 2.0) < 1e-6,
        "dP_drop/dV = {}, closed form {}",
        s.grad[0][0],
        p_in * dw_dv / 2.0
    );
}

/// The reported `fd_error` is a real error bar, not decoration: it has to be
/// small on parameters that behave, and it is the thing to look at when one
/// does not.
#[test]
fn the_reported_fd_error_is_small_on_well_scaled_parameters() {
    let net = parse_spice(DIVIDER).unwrap();
    let params = [ParamRef::new("R1", "value"), ParamRef::new("V1", "dc")];
    let s = dc_sensitivity(
        &net,
        &registry(&net),
        &tight(),
        &[Output::NodeVoltage("out".into())],
        &params,
    )
    .unwrap();
    // Both enter the stamp linearly, so the two step sizes agree exactly and
    // there is nothing for Richardson to remove.
    for (p, e) in params.iter().zip(s.fd_error.iter()) {
        assert!(*e < 1e-9, "{}.{}: fd_error {e:e}", p.element, p.param);
    }
}

// ---------------------------------------------------------------------------
// Honest failure
// ---------------------------------------------------------------------------

/// A diode model parameter reaches neither path today.  The contract is that it
/// is *reported*, not silently returned as a zero gradient — an optimiser
/// cannot tell a wrong zero from a real stationary point.
#[test]
fn an_unreachable_parameter_is_reported_rather_than_zeroed() {
    let net = parse_spice(DIODE_CLAMP).unwrap();
    let outputs = [Output::NodeVoltage("mid".into())];
    let params = [
        ParamRef::with_nominal("D1", "IS", 1e-14),
        ParamRef::new("R1", "value"),
    ];
    let s = dc_sensitivity(&net, &registry(&net), &tight(), &outputs, &params).unwrap();

    assert!(!s.reached[0], "D1.IS is not wired to anything yet");
    assert!(s.reached[1], "R1 is");
    let unreached = s.unreached(&params);
    assert_eq!(unreached.len(), 1);
    assert_eq!(unreached[0].param, "IS");
}

/// A typo'd element name is the same class of problem and gets the same answer.
#[test]
fn an_unknown_element_is_reported_rather_than_zeroed() {
    let net = parse_spice(DIVIDER).unwrap();
    let outputs = [Output::NodeVoltage("out".into())];
    let params = [ParamRef::new("R99", "value")];
    let s = dc_sensitivity(&net, &registry(&net), &tight(), &outputs, &params).unwrap();
    assert!(!s.reached[0]);
    assert_eq!(s.grad[0][0], 0.0);
}
