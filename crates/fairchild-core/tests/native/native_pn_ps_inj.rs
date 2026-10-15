//! Charge-based reactance on `fc_pn_ps_inj`, the forward-bias carrier-injection
//! drive.
//!
//! The diffusion charge is `q = τ·I`, not `C_d(V)·V`. With `C_d = τ·g_d` and
//! `g_d = (I + I_s)/V_T` the two differ by the factor `V/V_T`, which is about
//! 20 at normal forward bias, and which grows with bias. Integrating `C_d·V`
//! therefore gives a device far slower than its own `tau_carrier`, and its
//! transient stops short of the operating point the same deck reports.

use fairchild_core::{
    dc_op_nr_with_registry, tran_nr_with_registry_opts, DeviceRegistry, SimOptions,
};
use fairchild_parser::parse_spice;

const TAU: f64 = 10e-9;
/// Forward current before the step, and after it.
const I0: f64 = 1e-4;
const I1: f64 = 1e-3;

/// Newton's stopping test here is `|Δx| < abstol + reltol·|x|` and nothing
/// else — there is no residual check — so on a forward diode it stops while the
/// node is still ~0.7 mA out of KCL (#139). At the default `reltol` this deck
/// stalls at 0.28 mA whatever the charge model is, which would make the
/// assertion below a test of the stopping rule rather than of the charge.
/// Tightening until the solve is genuinely converged isolates the one variable.
/// Sabotage holds at this tolerance: `C_d·V` reaches 0.45 mA where `τ·I`
/// reaches 1.000 mA. Drop these back to the defaults when #139 lands.
const TOL: &str = "reltol=1e-9 abstol=1e-16 vntol=1e-9";

fn deck(tail: &str) -> String {
    format!(
        "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xpn ch0 out0 a 0 fc_pn_ps_inj L_um=500 tau_carrier={TAU:e}
{tail}
"
    )
}

/// The anode voltage the *static* solve reaches when `i` is forced into the
/// junction.
///
/// This is the model's own I-V map, which no charge model touches. It converts
/// the closed-form `I(t)` below into a voltage to compare against, so the test
/// never needs `I_s` or `V_T` of its own — and so the quantity under test, the
/// time constant, is not borrowed from the thing being tested.
fn v_at_current(i: f64) -> f64 {
    let net = parse_spice(&deck(&format!("Ic 0 a DC {i:e}\n.options {TOL}\n.op"))).expect("parse");
    dc_op_nr_with_registry(&net, &DeviceRegistry::new())
        .expect("DC OP converges")
        .node_voltage("a")
        .expect("node a")
}

/// ABSOLUTE anchor on the diffusion charge.
///
/// `q = τ·I` turns the node equation into `τ·dI/dt = I_ext − I`, so a step in
/// the forced current relaxes with exactly the carrier lifetime, at any bias
/// and for any step size. That is a closed form, and `tau_carrier` is a deck
/// parameter rather than anything the solver computes.
///
/// Against `q = C_d·V` the relaxation is about 20x slower, so the junction is
/// still near its starting current when this says it should be two thirds of
/// the way to the new one.
#[test]
fn injection_charge_relaxes_at_the_carrier_lifetime() {
    let i_of = |t: f64| I1 - (I1 - I0) * (-t / TAU).exp();
    let v0 = v_at_current(I0);

    let net = parse_spice(&deck(&format!(
        "Ic 0 a DC {I1:e}\n.ic v(a)={v0:.12}\n.options uic=1 {TOL}"
    )))
    .expect("parse");
    let opts = SimOptions::from_netlist(&net);
    let res =
        tran_nr_with_registry_opts(&net, TAU / 500.0, 10.0 * TAU, &DeviceRegistry::new(), &opts)
            .expect("transient");

    for k in [1.0, 2.0, 3.0, 10.0] {
        let t = k * TAU;
        let want = v_at_current(i_of(t));
        let got = res.voltage_at("a", t).expect("node a");
        assert!(
            (got - want).abs() < 2e-3,
            "at t = {k}·tau_carrier the junction should carry {:.4} mA, which this \
             model's own DC solve puts at v = {want:.6} V. The transient gives \
             v = {got:.6} V. Integrating C_d·V instead of τ·I lands far short.",
            i_of(t) * 1e3
        );
    }
}
