//! The transient adjoint has to differentiate the *discrete* system the
//! integrator actually solved — not the ODE it approximates.
//!
//! So the reference here is a full re-solve of the same run at a perturbed
//! parameter, on the same step grid, with the tolerance wound down until the
//! difference quotient is trustworthy.  Agreeing with the continuous-time
//! answer to a few digits would prove much less: it would also be consistent
//! with a gradient that had, say, the wrong integrator's history coefficients,
//! since both converge to the same limit as `h → 0`.
//!
//! The one closed-form check is kept as a sanity rail on top of that, at the
//! looseness the discretisation error actually allows.

use fairchild_core::adjoint_tran::TranAdjoint;
use fairchild_core::{DeviceRegistry, Output, ParamRef, SimOptions};
use fairchild_parser::{parse_spice, Netlist};

fn registry_for(netlist: &Netlist) -> DeviceRegistry {
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&netlist.models);
    reg
}

/// Tight enough that the re-solve reference is limited by the finite-difference
/// step rather than by Newton's stopping rule.
fn opts(method: &str) -> SimOptions {
    let mut o = SimOptions {
        reltol: 1e-13,
        vntol: 1e-15,
        abstol: 1e-18,
        ..SimOptions::default()
    };
    assert!(o.set("method", method));
    o
}

/// Weights selecting the value at the final timepoint.
fn at_end(n: usize) -> Vec<f64> {
    let mut w = vec![0.0; n];
    w[n - 1] = 1.0;
    w
}

struct Run {
    adj: TranAdjoint,
    reg: DeviceRegistry,
}

fn run(src: &str, method: &str, step: f64, stop: f64) -> Run {
    let net = parse_spice(src).unwrap();
    let reg = registry_for(&net);
    let adj = TranAdjoint::run(&net, &reg, &opts(method), step, stop).unwrap();
    Run { adj, reg }
}

/// The objective, re-solved from scratch with `element.param` set to `value`.
fn resolve_at(src: &str, method: &str, step: f64, stop: f64, out: &Output, at: (&str, f64)) -> f64 {
    let mut net = parse_spice(src).unwrap();
    let reg = registry_for(&net);
    let (path, value) = at;
    let (element, param) = path.split_once('.').unwrap();
    assert!(
        fairchild_core::netlist_edit::set_element_param(&mut net, element, param, value),
        "could not set {path}"
    );
    let adj = TranAdjoint::run(&net, &reg, &opts(method), step, stop).unwrap();
    let n = adj.time().len();
    adj.weighted(out, &at_end(n)).unwrap().0
}

/// `dL/dp` from the adjoint, and from re-solving the whole run either side of
/// the nominal value.  The comparison the whole module exists to pass.
fn adjoint_vs_resolve(
    src: &str,
    method: &str,
    step: f64,
    stop: f64,
    out: Output,
    path: &str,
    nominal: f64,
    fd_step: f64,
) -> (f64, f64) {
    let r = run(src, method, step, stop);
    let n = r.adj.time().len();
    let (_, seeds) = r.adj.weighted(&out, &at_end(n)).unwrap();
    let (element, param) = path.split_once('.').unwrap();
    let s = r
        .adj
        .gradient(&r.reg, &seeds, &[ParamRef::new(element, param)])
        .unwrap();
    assert!(s.reached[0], "{path} was not reached");

    let plus = resolve_at(src, method, step, stop, &out, (path, nominal + fd_step));
    let minus = resolve_at(src, method, step, stop, &out, (path, nominal - fd_step));
    (s.grad[0], (plus - minus) / (2.0 * fd_step))
}

fn assert_close(tag: &str, got: f64, want: f64, rtol: f64) {
    let err = (got - want).abs() / want.abs().max(1e-30);
    assert!(
        err <= rtol,
        "{tag}: adjoint {got:e} vs reference {want:e} — {err:e} relative, limit {rtol:e}"
    );
}

// ---------------------------------------------------------------------------
// Linear RC — where every term can be checked by hand
// ---------------------------------------------------------------------------

const RC: &str = "* rc step\n\
                  V1 in 0 PULSE(0 1 0 1p 1p 1 2)\n\
                  R1 in out 1k\n\
                  C1 out 0 1n\n\
                  .tran 1u 5u\n.end\n";

/// Backward Euler: the simplest history coupling there is — `∂G_k/∂x_{k-1}`
/// is one term, and getting its sign or its `α` wrong shows up immediately.
#[test]
fn backward_euler_capacitance_gradient_matches_a_full_resolve() {
    let (adj, fd) = adjoint_vs_resolve(
        RC,
        "be",
        2e-7,
        4e-6,
        Output::NodeVoltage("out".into()),
        "C1.c",
        1e-9,
        1e-13,
    );
    assert_close("dv(out)/dC under BE", adj, fd, 1e-6);
}

/// The resistor is purely resistive, so this isolates `∂F/∂p` from the charge
/// path — and it also exercises the `t = 0` term, because the operating point
/// itself depends on `R`.
#[test]
fn backward_euler_resistance_gradient_matches_a_full_resolve() {
    let (adj, fd) = adjoint_vs_resolve(
        RC,
        "be",
        2e-7,
        4e-6,
        Output::NodeVoltage("out".into()),
        "R1.r",
        1e3,
        1e-4,
    );
    assert_close("dv(out)/dR under BE", adj, fd, 1e-6);
}

/// Trapezoidal is the one that needs the second co-state: `i_k` depends on
/// every earlier step through `−i_{k−1}`, so dropping the `ū` recursion leaves
/// a gradient that is wrong by a factor that grows with the number of steps —
/// which is exactly what this catches, and what BE cannot.
#[test]
fn trapezoidal_needs_the_second_costate_and_gets_it_right() {
    let (adj, fd) = adjoint_vs_resolve(
        RC,
        "tr",
        2e-7,
        4e-6,
        Output::NodeVoltage("out".into()),
        "C1.c",
        1e-9,
        1e-13,
    );
    assert_close("dv(out)/dC under TR", adj, fd, 1e-6);
}

/// A rail on top of the re-solve comparison: the RC step response is
/// `v = V(1 − e^{−t/RC})`, so `∂v/∂C = −V·t·e^{−t/RC}/(R·C²)`.  Loose, because
/// the integrator is solving a difference equation and not that exponential —
/// but tight enough to catch a gradient that is right about the discrete system
/// and wrong about the circuit.
#[test]
fn the_capacitance_gradient_is_the_analytic_one_to_within_discretisation_error() {
    let r = run(RC, "tr", 5e-8, 4e-6);
    let n = r.adj.time().len();
    let out = Output::NodeVoltage("out".into());
    let (_, seeds) = r.adj.weighted(&out, &at_end(n)).unwrap();
    let s = r
        .adj
        .gradient(&r.reg, &seeds, &[ParamRef::new("C1", "c")])
        .unwrap();

    let (rr, c, t) = (1e3, 1e-9, *r.adj.time().last().unwrap());
    let analytic = -t * (-t / (rr * c)).exp() / (rr * c * c);
    assert_close("dv/dC vs the exponential", s.grad[0], analytic, 2e-3);
}

/// A time integral, rather than a value at one instant: every timepoint now
/// seeds the backward pass, so a co-state recursion that only happened to work
/// with a single non-zero seed at the end would fail here.
#[test]
fn an_integral_objective_matches_a_full_resolve() {
    let (step, stop) = (2e-7, 4e-6);
    let r = run(RC, "tr", step, stop);
    let n = r.adj.time().len();
    let out = Output::NodeVoltage("out".into());
    let weights = vec![step; n];
    let (_, seeds) = r.adj.weighted(&out, &weights).unwrap();
    let s = r
        .adj
        .gradient(&r.reg, &seeds, &[ParamRef::new("C1", "c")])
        .unwrap();

    let integral = |c: f64| {
        let mut net = parse_spice(RC).unwrap();
        assert!(fairchild_core::netlist_edit::set_element_param(
            &mut net, "C1", "c", c
        ));
        let a = TranAdjoint::run(&net, &registry_for(&net), &opts("tr"), step, stop).unwrap();
        let n = a.time().len();
        a.weighted(&out, &vec![step; n]).unwrap().0
    };
    let d = 1e-13;
    let fd = (integral(1e-9 + d) - integral(1e-9 - d)) / (2.0 * d);
    assert_close("d(∫v dt)/dC", s.grad[0], fd, 1e-6);
}

// ---------------------------------------------------------------------------
// Nonlinear
// ---------------------------------------------------------------------------

/// A diode clamp with a capacitor: the Jacobian now changes every iterate, and
/// the device carries a junction limiter across evaluations.  If the replay
/// disturbed device state the residual it differences would not be `f(x_k)`.
#[test]
fn a_nonlinear_transient_matches_a_full_resolve() {
    const SRC: &str = "* diode clamp\n\
                       .model dmod D (IS=1e-14 N=1.0)\n\
                       V1 in 0 PULSE(0 2 0 1n 1n 1 2)\n\
                       R1 in mid 1k\n\
                       D1 mid 0 dmod\n\
                       C1 mid 0 10p\n\
                       .tran 1n 20n\n.end\n";
    let (adj, fd) = adjoint_vs_resolve(
        SRC,
        "tr",
        1e-9,
        2e-8,
        Output::NodeVoltage("mid".into()),
        "C1.c",
        1e-11,
        1e-15,
    );
    // Looser than the linear cases, and the reference is the reason: the clamp
    // is curved enough in `C` that the re-solve is squeezed from both sides —
    // truncation at the large steps, the solver's own tolerance at the small
    // ones — and it only holds three or four figures.  Sweeping it puts the
    // plateau within 4e-8 of the adjoint; 5e-5 is the width of the plateau, not
    // of the disagreement.
    assert_close("dv(mid)/dC across a diode", adj, fd, 5e-5);
}

// ---------------------------------------------------------------------------
// Electro-optic — the case the whole feature exists for
// ---------------------------------------------------------------------------

/// A modulated link with a low-passed detector: light out of an MZM, current
/// into an RC load.  The gradient has to cross a frozen optical coefficient
/// *and* a charge history to get from `V_pi` to `v(pout, T)`.
///
/// `fc_mzm` reaches its fixed point by successive substitution on a cached
/// transmission amplitude, so `∂(out_re)/∂v_mod` is never stamped.  Without the
/// frozen-column repair this gradient is not merely inaccurate — it is exactly
/// zero, which reads as "the modulator does not affect the detector".
#[test]
fn an_electro_optic_link_has_a_live_transient_gradient() {
    const SRC: &str = ".optical_port in0\n.optical_port out0\n\
                       Xl0 in0 fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
                       Xmzm in0 out0 vsig 0 fc_mzm V_pi=3.0 alpha=1.0 e_r=1000\n\
                       Xpd out0 pout 0 fc_photodetector responsivity=0.8\n\
                       Rl pout 0 1k\n\
                       Cl pout 0 100f\n\
                       Vsig vsig 0 PULSE(0 1.5 0 1n 1n 1 2)\n\
                       .tran 100p 4n\n.end\n";
    let (step, stop) = (1e-10, 3e-9);
    let r = run(SRC, "tr", step, stop);
    let n = r.adj.time().len();
    let out = Output::NodeVoltage("pout".into());
    let (_, seeds) = r.adj.weighted(&out, &at_end(n)).unwrap();
    let s = r
        .adj
        .gradient(
            &r.reg,
            &seeds,
            &[ParamRef::new("Xmzm", "V_pi"), ParamRef::new("Cl", "c")],
        )
        .unwrap();

    assert!(
        s.reached.iter().all(|ok| *ok),
        "unreached: {:?}",
        s.unreached(&[ParamRef::new("Xmzm", "V_pi"), ParamRef::new("Cl", "c")])
    );
    assert!(
        s.grad[0] != 0.0,
        "dv(pout)/dV_pi came out exactly zero — the frozen optical coefficient \
         is not being re-derived"
    );

    // Reference: re-solve either side.  `V_pi` lives on the instance line, so
    // the netlist editor reaches it the same way the adjoint does.
    let solve = |v_pi: f64| {
        let mut net = parse_spice(SRC).unwrap();
        assert!(fairchild_core::netlist_edit::set_element_param(
            &mut net, "Xmzm", "V_pi", v_pi
        ));
        let a = TranAdjoint::run(&net, &registry_for(&net), &opts("tr"), step, stop).unwrap();
        let n = a.time().len();
        a.weighted(&out, &at_end(n)).unwrap().0
    };
    let d = 1e-6;
    let fd = (solve(3.0 + d) - solve(3.0 - d)) / (2.0 * d);
    assert_close("dv(pout)/dV_pi", s.grad[0], fd, 1e-4);
}

/// Optical power at a net, rather than a node voltage — the objective an
/// optical design actually optimises, and a nonlinear one, so its seed has to
/// be taken at each `x_k` rather than once.
#[test]
fn optical_power_is_a_transient_objective_too() {
    const SRC: &str = ".optical_port in0\n.optical_port out0\n\
                       Xl0 in0 fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
                       Xmzm in0 out0 vsig 0 fc_mzm V_pi=3.0 alpha=1.0 e_r=1000\n\
                       Rd vsig drv 500\n\
                       Cd vsig 0 200f\n\
                       Vsig drv 0 PULSE(0 1.5 0 1n 1n 1 2)\n\
                       .tran 100p 4n\n.end\n";
    let (step, stop) = (1e-10, 3e-9);
    let r = run(SRC, "tr", step, stop);
    let n = r.adj.time().len();
    let out = Output::OpticalPower {
        net: "out0".into(),
        channel: 0,
    };
    let (_, seeds) = r.adj.weighted(&out, &at_end(n)).unwrap();
    let s = r
        .adj
        .gradient(&r.reg, &seeds, &[ParamRef::new("Cd", "c")])
        .unwrap();
    assert!(s.reached[0]);

    let solve = |c: f64| {
        let mut net = parse_spice(SRC).unwrap();
        assert!(fairchild_core::netlist_edit::set_element_param(
            &mut net, "Cd", "c", c
        ));
        let a = TranAdjoint::run(&net, &registry_for(&net), &opts("tr"), step, stop).unwrap();
        let n = a.time().len();
        a.weighted(&out, &at_end(n)).unwrap().0
    };
    // The re-solve reference is the noisy side here: the optical power moves by
    // one part in 1e11 over a 0.05 % change in `Cd`, so the difference quotient
    // only clears the solver's own tolerance at the larger steps.
    let d = 1e-16;
    let fd = (solve(200e-15 + d) - solve(200e-15 - d)) / (2.0 * d);
    assert_close("dP(out0)/dC_drive", s.grad[0], fd, 2e-3);
}

// ---------------------------------------------------------------------------
// Contracts
// ---------------------------------------------------------------------------

/// A parameter that reaches nothing must be reported, not returned as a zero
/// gradient — the two are indistinguishable to an optimiser, and one of them
/// looks like a stationary point.
#[test]
fn an_unreachable_parameter_is_reported_rather_than_zeroed() {
    let r = run(RC, "tr", 5e-7, 2e-6);
    let n = r.adj.time().len();
    let out = Output::NodeVoltage("out".into());
    let (_, seeds) = r.adj.weighted(&out, &at_end(n)).unwrap();
    let params = [
        ParamRef::new("C1", "c"),
        ParamRef::with_nominal("Rnope", "r", 1.0),
    ];
    let s = r.adj.gradient(&r.reg, &seeds, &params).unwrap();
    assert!(s.reached[0], "C1 should be reachable");
    assert!(!s.reached[1], "a missing element cannot be reached");
    assert_eq!(s.grad[1], 0.0);
    assert_eq!(s.unreached(&params).len(), 1);
}

/// Under UIC the initial state is imposed rather than solved, so `t = 0` stops
/// being a constraint: it gets no co-state and contributes no `∂G_0/∂p`.  The
/// gradient still has to be right, which is what makes this a separate path
/// worth walking rather than a special case that happens to fall out.
#[test]
fn uic_drops_the_initial_condition_from_the_constraint_set() {
    const SRC: &str = "* rc with uic\n\
                       V1 in 0 PULSE(0 1 0 1p 1p 1 2)\n\
                       R1 in out 1k\n\
                       C1 out 0 1n\n\
                       .ic v(out)=0.25\n\
                       .tran 1u 5u UIC\n.end\n";
    let (step, stop) = (2e-7, 4e-6);
    let mut o = opts("be");
    o.uic = true;
    let net = parse_spice(SRC).unwrap();
    let reg = registry_for(&net);
    let adj = TranAdjoint::run(&net, &reg, &o, step, stop).unwrap();
    let n = adj.time().len();
    let out = Output::NodeVoltage("out".into());
    let (_, seeds) = adj.weighted(&out, &at_end(n)).unwrap();
    let s = adj
        .gradient(&reg, &seeds, &[ParamRef::new("R1", "r")])
        .unwrap();
    assert!(s.reached[0]);

    let solve = |r: f64| {
        let mut net = parse_spice(SRC).unwrap();
        assert!(fairchild_core::netlist_edit::set_element_param(
            &mut net, "R1", "r", r
        ));
        let a = TranAdjoint::run(&net, &registry_for(&net), &o, step, stop).unwrap();
        let m = a.time().len();
        a.weighted(&out, &at_end(m)).unwrap().0
    };
    let d = 1e-4;
    let fd = (solve(1e3 + d) - solve(1e3 - d)) / (2.0 * d);
    assert_close("dv(out)/dR under UIC", s.grad[0], fd, 1e-6);
}

/// Inductance does not fit this formulation, so it is refused rather than
/// silently mis-differentiated.
#[test]
fn inductance_is_rejected_rather_than_mis_differentiated() {
    const SRC: &str = "* rl\nV1 in 0 PULSE(0 1 0 1n 1n 1 2)\n\
                       R1 in out 50\nL1 out 0 1u\n.tran 1n 20n\n.end\n";
    let net = parse_spice(SRC).unwrap();
    let reg = registry_for(&net);
    let msg = match TranAdjoint::run(&net, &reg, &opts("tr"), 1e-9, 2e-8) {
        Err(e) => e.to_string(),
        Ok(_) => panic!("an inductor was accepted"),
    };
    assert!(msg.contains("inductance"), "unhelpful refusal: {msg}");
    assert!(msg.contains("l1"), "should name the offender: {msg}");
}

/// The capture pass re-stamps and probes between the solve and the history
/// advance.  If any of that leaked into the integrator state, the trajectory
/// would drift away from a plain `.tran` — so hold it to being identical.
#[test]
fn capturing_the_adjoint_state_does_not_perturb_the_run() {
    let net = parse_spice(RC).unwrap();
    let reg = registry_for(&net);
    let o = opts("tr");
    let (step, stop) = (2e-7, 4e-6);
    let adj = TranAdjoint::run(&net, &reg, &o, step, stop).unwrap();
    let batch =
        fairchild_core::tran::tran_nr_with_registry_opts(&net, step, stop, &reg, &o).unwrap();

    let out = &batch.node_voltages["out"];
    for (k, t) in adj.time().iter().enumerate() {
        assert_eq!(*t, batch.time[k], "timepoint {k} diverged");
        let v = adj
            .topology()
            .node_voltage("out", &adj.trajectory()[k])
            .unwrap();
        assert_eq!(v, out[k], "value at timepoint {k} diverged");
    }
}

/// **A device-declared reactive branch — the one parameter class no other test
/// here covered, and the one that was wrong.**
///
/// Every other transient test perturbs a netlist `R`/`C`, which never routes
/// through `TranStepper::set_device_param`, or a device parameter that is
/// purely resistive. Deleting the history re-seed entirely left all of them
/// green, and this gradient came out 16 % wrong with nothing to say so.
///
/// The cause was not the re-seed. `stamp_device_branches` stamped `α·C(v)`
/// where the branch carries `q = C(v)·v`, so the true derivative is
/// `α·(C + v·dC/dv)`. Newton converges either way — the missing term is
/// absorbed by successive substitution — so the forward answer was right and
/// only the Jacobian was short. The adjoint needs the real `∂f/∂x`, and
/// `jacobian_check_tran` now measures it: the term was 22 % of the drive
/// node's diagonal, and `v·dC/dv/h` predicted the gap to six figures.
///
/// The interferometer matters: a *single* phase shifter's output power does not
/// depend on phase at all, so `c_j0` legitimately has zero gradient there and a
/// test built that way passes while proving nothing. That was the first version
/// of this test, and its "reference" was two near-identical numbers over 1e-18.
#[test]
fn a_device_declared_capacitance_gradient_matches_a_full_resolve() {
    const SRC: &str = ".optical_port in0\n.optical_port dk\n.optical_port a1\n\
                       .optical_port a2\n.optical_port b1\n.optical_port b2\n\
                       .optical_port out0\n.optical_port ou\n\
                       Xl0 in0 fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
                       Xc1 in0 dk a1 a2 fc_dcoupler kappa_L=0.7853981634\n\
                       Xps a1 b1 drv 0 fc_pn_ps_cap l_um=3000 v_pi_l=0.012 \
                         c_j0=750f v_bi=0.917 m_j=0.5 alpha_dB_cm=2.0 pin_at_ref=1\n\
                       Xref a2 b2 0 0 fc_pn_ps_cap l_um=3000 v_pi_l=0.012 \
                         c_j0=750f v_bi=0.917 m_j=0.5 alpha_dB_cm=2.0 pin_at_ref=1\n\
                       Xc2 b1 b2 out0 ou fc_dcoupler kappa_L=0.7853981634\n\
                       Rd drv nd 25\n\
                       Vsig nd 0 PULSE(-3 -1 0 40p 40p 1 2)\n\
                       .tran 4p 200p\n.end\n";
    // Sample while the edge is still moving: once the drive has settled, the
    // capacitance has no influence left on the final value and the reference
    // difference quotient is pure noise.
    let (step, stop) = (4e-12, 6e-11);
    let out = Output::OpticalPower {
        net: "out0".into(),
        channel: 0,
    };
    let nominal = 750e-15;

    let r = run(SRC, "be", step, stop);
    let n = r.adj.time().len();
    let (_, seeds) = r.adj.weighted(&out, &at_end(n)).unwrap();
    let s = r
        .adj
        .gradient(&r.reg, &seeds, &[ParamRef::new("Xps", "c_j0")])
        .unwrap();
    assert!(s.reached[0], "Xps.c_j0 was not reached");

    // Reference: set it on the netlist and re-solve from scratch — a path that
    // never touches `set_device_param`, which is what makes this an anchor
    // rather than two copies of one mistake agreeing.
    let solve = |c: f64| {
        let mut net = parse_spice(SRC).unwrap();
        assert!(fairchild_core::netlist_edit::set_element_param(
            &mut net, "Xps", "c_j0", c
        ));
        let a = TranAdjoint::run(&net, &registry_for(&net), &opts("be"), step, stop).unwrap();
        let m = a.time().len();
        a.weighted(&out, &at_end(m)).unwrap().0
    };
    let d = nominal * 1e-3;
    let fd = (solve(nominal + d) - solve(nominal - d)) / (2.0 * d);
    assert_close("dP(out0)/dc_j0", s.grad[0], fd, 1e-3);
}

// ---------------------------------------------------------------------------
// The two Jacobian diagnostics
//
// Both were written to chase real defects — `jacobian_check_tran` measured the
// missing `v·dC/dv` term at 22 % of the drive node's diagonal (see the note on
// the interferometer test above), and `charge_lag` explains why a
// bias-dependent capacitance is stamped one step stale. Both then had no
// caller, so a regression in either would have gone unnoticed until the next
// time someone needed them and found they no longer worked.
//
// A linear RC is the case where each has a known exact answer, which is what
// makes them usable as tests rather than as agreement invariants.
// ---------------------------------------------------------------------------

const LINEAR_RC: &str = "* linear RC\n\
                         V1 in 0 PULSE(0 1 0 1n 1n 10m 20m)\n\
                         R1 in out 1k\n\
                         C1 out 0 1u\n\
                         .tran 100u 2m\n.end\n";

fn linear_rc() -> (fairchild_parser::Netlist, DeviceRegistry, SimOptions) {
    let net = fairchild_parser::parse_spice(LINEAR_RC).unwrap();
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);
    let opts = SimOptions::from_netlist(&net);
    (net, reg, opts)
}

#[test]
fn jacobian_check_reports_no_mismatch_on_a_linear_circuit() {
    let (net, reg, opts) = linear_rc();
    let bad = fairchild_core::adjoint_tran::jacobian_check_tran(
        &net, &reg, &opts, 100e-6, 1e-3, 1e-6, 1e-12,
    )
    .unwrap();
    // Every stamp in an R/C/V circuit is a constant, so the stamped Jacobian
    // and the finite-difference one must agree to the tolerance given. A
    // mismatch here is either a wrong stamp or a broken checker.
    let undeclared: Vec<_> = bad.iter().filter(|m| !m.frozen).collect();
    assert!(
        undeclared.is_empty(),
        "linear circuit should have no undeclared Jacobian mismatch, got {:?}",
        undeclared
            .iter()
            .map(|m| (m.row, m.col, m.stamped, m.numeric))
            .collect::<Vec<_>>()
    );
}

#[test]
fn charge_lag_reports_a_device_charge_branch_and_finds_no_lag_at_a_static_bias() {
    // A *netlist* capacitor declares no device branch, so `charge_lag` returns
    // an empty list for an ordinary RC and any assertion over it is vacuous —
    // which is how the first version of this test passed against a deliberately
    // corrupted diagnostic. It needs a device that declares its own charge.
    let src = "* pn phase shifter with Cj(V)\n\
               .optical_port a\n.optical_port b\n\
               Xl a fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
               Xps a b vm 0 fc_pn_ps_cap L_um=500 V_pi_L=2e-3 g_pn=1e-3\n\
               Vm vm 0 PULSE(0 2 0 1n 1n 50n 100n)\n\
               .tran 1n 20n\n.end\n";
    let net = fairchild_parser::parse_spice(src).unwrap();
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);
    reg.register_native_photonics();
    let opts = SimOptions::from_netlist(&net);

    let lag = fairchild_core::adjoint_tran::charge_lag(&net, &reg, &opts, 1e-9, 1e-8).unwrap();
    assert!(
        !lag.is_empty(),
        "fc_pn_ps_cap declares a charge branch; an empty result means the \
         diagnostic stopped seeing device branches at all"
    );
    // t = 10 ns sits on the pulse plateau, so the bias is not moving between
    // steps and the stale-C effect this measures has nothing to bite on. The
    // two companions must therefore be the same number.
    for (name, stamped, refreshed) in &lag {
        assert!(
            (stamped - refreshed).abs() <= 1e-12 * stamped.abs().max(1.0),
            "static bias should show no lag on '{name}': \
             stamped={stamped:.6e} refreshed={refreshed:.6e}"
        );
    }
}
