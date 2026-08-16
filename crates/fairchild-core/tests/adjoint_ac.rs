//! The AC adjoint differentiates a frequency response.
//!
//! Two kinds of reference, and both are needed. A full re-solve of the sweep at
//! a perturbed parameter checks that the adjoint differentiates *the system
//! that was solved*. A closed form checks that the system is the right one —
//! an adjoint can be perfectly self-consistent about a wrong assembly, and only
//! an external anchor sees that.

use fairchild_core::adjoint_ac::{AcAdjoint, AcOutput};
use fairchild_core::{DeviceRegistry, ParamRef, SimOptions};
use fairchild_parser::{parse_spice, Netlist};

fn registry_for(netlist: &Netlist) -> DeviceRegistry {
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&netlist.models);
    reg
}

/// `gmin = 0` matters here and is not tidiness. The default 1e-12 shunts every
/// node, so the circuit fairchild solves is an RC low-pass **in parallel with a
/// 1 TΩ resistor** — and the textbook `|H|² = 1/(1+(ωRC)²)` is not that
/// circuit. Below the corner the shunt is the whole of `d|H|²/dR`'s
/// disagreement: a constant −2e-12 absolute offset, which reads as 2.5e-5
/// *relative* at 1 kHz where the derivative itself is 8e-8. The adjoint was
/// right and the reference was incomplete. Removing gmin makes the closed form
/// exact and the agreement 1e-8 or better across five decades.
fn tight() -> SimOptions {
    SimOptions {
        reltol: 1e-13,
        vntol: 1e-15,
        abstol: 1e-18,
        gmin: 0.0,
        ..SimOptions::default()
    }
}

/// The photonic decks keep gmin. `tight()` zeroes it so the RC's textbook
/// transfer function is exact, but an optical bundle has wires no element
/// conducts to ground — a dark port, a λ tag — and without a floor those rows
/// are structurally empty. The closed forms below are the simulator's own
/// re-solve, not a textbook, so nothing needs gmin gone.
fn tight_gmin() -> SimOptions {
    SimOptions {
        reltol: 1e-13,
        vntol: 1e-15,
        abstol: 1e-18,
        ..SimOptions::default()
    }
}

fn assert_close(tag: &str, got: f64, want: f64, rtol: f64) {
    let err = (got - want).abs() / want.abs().max(1e-30);
    assert!(
        err <= rtol,
        "{tag}: adjoint {got:e} vs reference {want:e} — {err:e} relative, limit {rtol:e}"
    );
}

// ---------------------------------------------------------------------------
// RC low-pass — every term checkable by hand
// ---------------------------------------------------------------------------

const RC: &str = "* rc low pass\n\
                  V1 in 0 DC 0 AC 1\n\
                  R1 in out 1k\n\
                  C1 out 0 1n\n\
                  .end\n";

fn rc_run(freqs: &[f64]) -> (AcAdjoint, DeviceRegistry) {
    let net = parse_spice(RC).unwrap();
    let reg = registry_for(&net);
    let adj = AcAdjoint::run(&net, &reg, &tight(), freqs, Some("V1")).unwrap();
    (adj, reg)
}

/// `|H|² = 1/(1 + (ωRC)²)`, so `∂|H|²/∂R = −2ω²RC²·|H|⁴`.
///
/// This is the anchor: it does not go through the adjoint, the solver, or the
/// stamp, so it cannot agree with them by sharing a mistake.
#[test]
fn magnitude_gradient_matches_the_rc_closed_form() {
    let (r, c) = (1e3, 1e-9);
    for &f in &[1e3, 1e5, 1.0 / (2.0 * std::f64::consts::PI * r * c), 1e7] {
        let (adj, reg) = rc_run(&[f]);
        let out = AcOutput::MagSquared { node: "out".into() };
        let (_, seeds) = adj.weighted(&out, &[1.0]).unwrap();
        let s = adj
            .gradient(&reg, &seeds, &[ParamRef::new("R1", "r")])
            .unwrap();
        assert!(s.reached[0], "R1.r was not reached at {f} Hz");

        let w = 2.0 * std::f64::consts::PI * f;
        let h2 = 1.0 / (1.0 + (w * r * c).powi(2));
        let want = -2.0 * w * w * r * c * c * h2 * h2;
        assert_close(&format!("d|H|²/dR at {f:e} Hz"), s.grad[0], want, 1e-6);
    }
}

/// The same for `C`: `∂|H|²/∂C = −2ω²R²C·|H|⁴`.
#[test]
fn capacitance_gradient_matches_the_rc_closed_form() {
    let (r, c) = (1e3, 1e-9);
    let f = 1e5;
    let (adj, reg) = rc_run(&[f]);
    let out = AcOutput::MagSquared { node: "out".into() };
    let (_, seeds) = adj.weighted(&out, &[1.0]).unwrap();
    let s = adj
        .gradient(&reg, &seeds, &[ParamRef::new("C1", "c")])
        .unwrap();
    assert!(s.reached[0]);

    let w = 2.0 * std::f64::consts::PI * f;
    let h2 = 1.0 / (1.0 + (w * r * c).powi(2));
    let want = -2.0 * w * w * r * r * c * h2 * h2;
    assert_close("d|H|²/dC", s.grad[0], want, 1e-6);
}

/// A band, not a point: the whole reason to use an adjoint here is that one
/// backward pass serves every frequency at once. A co-state recursion that
/// happened to work for a single non-zero seed would fail this.
#[test]
fn a_banded_objective_matches_a_full_resolve() {
    let freqs: Vec<f64> = (0..24).map(|i| 1e3 * 10f64.powf(i as f64 / 8.0)).collect();
    let (adj, reg) = rc_run(&freqs);
    let out = AcOutput::MagSquared { node: "out".into() };
    // Least-squares against a flat target: L = Σ (|H|² − t)².  dL/d|H|² = 2(|H|²−t),
    // which is what the weights carry.
    let target = 0.5;
    let resp = adj.response(&out).unwrap();
    let weights: Vec<f64> = resp.iter().map(|h| 2.0 * (h - target)).collect();
    let (_, seeds) = adj.weighted(&out, &weights).unwrap();
    let s = adj
        .gradient(&reg, &seeds, &[ParamRef::new("R1", "r")])
        .unwrap();
    assert!(s.reached[0]);

    // Reference: the same least-squares scalar, re-solved either side.
    let loss = |rv: f64| {
        let mut net = parse_spice(RC).unwrap();
        assert!(fairchild_core::netlist_edit::set_element_param(
            &mut net, "R1", "r", rv
        ));
        let a = AcAdjoint::run(&net, &registry_for(&net), &tight(), &freqs, Some("V1")).unwrap();
        a.response(&out)
            .unwrap()
            .iter()
            .map(|h| (h - target).powi(2))
            .sum::<f64>()
    };
    let d = 1e3 * 1e-4;
    let fd = (loss(1e3 + d) - loss(1e3 - d)) / (2.0 * d);
    assert_close("dL/dR over a band", s.grad[0], fd, 1e-4);
}

/// Real and imaginary parts differentiate too, and with the right sign — a
/// magnitude objective would hide a transposed-vs-conjugate-transposed mistake
/// because `|·|²` is blind to it.
#[test]
fn the_imaginary_part_carries_a_signed_gradient() {
    let f = 1e5;
    let (adj, reg) = rc_run(&[f]);
    for (tag, out) in [
        ("Re", AcOutput::Real { node: "out".into() }),
        ("Im", AcOutput::Imag { node: "out".into() }),
    ] {
        let (_, seeds) = adj.weighted(&out, &[1.0]).unwrap();
        let s = adj
            .gradient(&reg, &seeds, &[ParamRef::new("R1", "r")])
            .unwrap();
        assert!(s.reached[0]);

        let value = |rv: f64| {
            let mut net = parse_spice(RC).unwrap();
            assert!(fairchild_core::netlist_edit::set_element_param(
                &mut net, "R1", "r", rv
            ));
            let a = AcAdjoint::run(&net, &registry_for(&net), &tight(), &[f], Some("V1")).unwrap();
            a.response(&out).unwrap()[0]
        };
        let d = 1e-1;
        let fd = (value(1e3 + d) - value(1e3 - d)) / (2.0 * d);
        assert_close(&format!("d{tag}(V_out)/dR"), s.grad[0], fd, 1e-4);
    }
}

/// An unreachable parameter is reported, not returned as a confident zero —
/// the same contract the other two adjoints hold.
#[test]
fn an_unreachable_parameter_is_reported_rather_than_zeroed() {
    let (adj, reg) = rc_run(&[1e5]);
    let out = AcOutput::MagSquared { node: "out".into() };
    let (_, seeds) = adj.weighted(&out, &[1.0]).unwrap();
    let s = adj
        .gradient(&reg, &seeds, &[ParamRef::new("R1", "not_a_parameter")])
        .unwrap();
    assert!(
        !s.reached[0],
        "a nonexistent parameter must not report reached"
    );
    assert_eq!(s.grad[0], 0.0);
}

// ---------------------------------------------------------------------------
// The case this exists for: a photonic filter, differentiated w.r.t. geometry
// ---------------------------------------------------------------------------

/// A phase shifter's bias moving an interferometer's transmission — "put the
/// null here", differentiated.
///
/// This is the shape a photonic design problem takes: `|H(f)|²` against a
/// target, with geometry and bias as the parameters. Both are exercised — an
/// optical *length*, whose sensitivity is a near-total cancellation between two
/// much larger phase terms, and a *bias*, which reaches the answer through the
/// DC operating point the small-signal matrices were linearised about. The
/// second is the one an enumerated `∂A/∂p` would silently miss.
#[test]
fn an_interferometer_differentiates_by_geometry_and_by_bias() {
    const MZI: &str = ".optical_port in0\n.optical_port dk\n.optical_port a1\n.optical_port a2\n\
                       .optical_port b1\n.optical_port b2\n.optical_port out0\n.optical_port ou\n\
                       Xl0 in0 fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
                       Vd drv 0 DC -2 AC 1\n\
                       Xc1 in0 dk a1 a2 fc_dcoupler kappa_L=0.7853981634\n\
                       Xps a1 b1 drv 0 fc_pn_ps_cap l_um=2000 v_pi_l=0.012 c_j0=500f \
                         v_bi=0.917 m_j=0.5 alpha_dB_cm=2.0 pin_at_ref=1\n\
                       Xref a2 b2 fc_waveguide l_um=2000 alpha_dB_cm=2.0\n\
                       Xc2 b1 b2 out0 ou fc_dcoupler kappa_L=0.7853981634\n\
                       Rd drv 0 1meg\n\
                       .end\n";
    let freqs = [1e8, 1e9, 5e9];
    let net = parse_spice(MZI).unwrap();
    let reg = registry_for(&net);
    let adj = AcAdjoint::run(&net, &reg, &tight_gmin(), &freqs, Some("Vd")).unwrap();
    let out = AcOutput::MagSquared {
        node: "out0_re_0".into(),
    };
    let (_, seeds) = adj.weighted(&out, &[1.0, 1.0, 1.0]).unwrap();

    for (element, param, nominal, step, rtol) in [
        ("Xps", "l_um", 2000.0, Some(1e-3), 2e-3),
        ("Xps", "c_j0", 500e-15, None, 2e-3),
    ] {
        let mut pr = ParamRef::new(element, param);
        pr.step = step;
        let s = adj.gradient(&reg, &seeds, &[pr]).unwrap();
        assert!(s.reached[0], "{element}.{param} was not reached");

        let loss = |v: f64| {
            let mut n = parse_spice(MZI).unwrap();
            assert!(fairchild_core::netlist_edit::set_element_param(
                &mut n, element, param, v
            ));
            let a =
                AcAdjoint::run(&n, &registry_for(&n), &tight_gmin(), &freqs, Some("Vd")).unwrap();
            a.response(&out).unwrap().iter().sum::<f64>()
        };
        let d = step.unwrap_or(nominal * 1e-4);
        let fd = (loss(nominal + d) - loss(nominal - d)) / (2.0 * d);
        assert_close(&format!("dΣ|H|²/d{param}"), s.grad[0], fd, rtol);
    }
}
