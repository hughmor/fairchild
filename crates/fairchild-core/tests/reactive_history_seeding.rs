//! A transient must start from the operating point it was handed, including
//! the charge on every bias-dependent capacitance.
//!
//! `ReactiveState::new` seeds each branch's history from `reactive_branches()`,
//! which reports a bias-dependent `C` from whatever the device's last `eval`
//! cached. The transient builds *fresh* device instances that the DC solve
//! never touched, so without an eval at the operating point first, a depletion
//! junction was seeded with its constructor default of `c_j_cached = 0` — i.e.
//! `q_prev = 0` — while the first step's companion used the real `C_j(V_op)`.
//!
//! The branch then had to acquire its entire charge in one timestep. It is a
//! large error and a brief one, which is the combination that hides: the
//! waveform is visibly right a few tens of picoseconds later.

use fairchild_core::{tran_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

/// A reverse-biased PN phase shifter driven through a source resistance. Its
/// `C_j(V)` is the bias-dependent branch; nothing else in the deck is reactive.
fn deck(step_s: &str, stop_s: &str) -> String {
    format!(
        "* PN phase shifter held at reverse bias, then ramped\n\
         .optical_port a\n\
         .optical_port b\n\
         Xl a fc_cw_laser power_mW=1\n\
         Vd nd 0 PWL(0 -3 25p -5 1n -5)\n\
         Rd nd n 25\n\
         Xps a b n 0 fc_pn_ps_cap l_um=3000 v_pi_l=0.012 c_j0=750f \
           v_bi=0.917 m_j=0.5 alpha_dB_cm=2.0\n\
         .tran {step_s} {stop_s}\n.end\n"
    )
}

fn run(step: f64, stop: f64) -> (Vec<f64>, Vec<f64>) {
    let net = parse_spice(&deck(&format!("{step:e}"), &format!("{stop:e}"))).unwrap();
    let r = tran_nr_with_registry(&net, step, stop, &DeviceRegistry::new())
        .expect("transient should run");
    (r.time.clone(), r.node_voltages["n"].clone())
}

/// The junction node must not move more in one timestep than its own RC allows.
///
/// With the history seeded at zero charge this node jumped from −2.93 V to
/// −0.26 V on the first 1 ps step — 2.7 V, against a drive that had moved
/// 0.08 V. The bound below is deliberately generous: anything within an order
/// of the drive's own motion is fine, and the bug overshoots by 30×.
#[test]
fn a_bias_dependent_cap_starts_charged_to_its_operating_point() {
    let (t, v) = run(1e-12, 40e-12);
    assert!(t.len() > 10, "need several steps to judge");

    // The drive ramps 2 V over 25 ps, so 0.08 V per 1 ps step. The junction
    // node lags it; it certainly never leads it.
    let first_jump = (v[1] - v[0]).abs();
    assert!(
        first_jump < 0.5,
        "the junction node moved {first_jump:.3} V on the first step, from {:.4} to {:.4}. \
         The drive moved 0.08 V. A capacitance seeded at zero charge has to acquire \
         all of it in one step, which is what this looks like.",
        v[0],
        v[1]
    );

    // And monotone: the drive ramps one way, so the node has no business
    // reversing while it does.
    for k in 1..(25.min(v.len() - 1)) {
        assert!(
            v[k + 1] <= v[k] + 1e-6,
            "junction node reversed at step {k}: {:.4} -> {:.4}",
            v[k],
            v[k + 1]
        );
    }
}

/// **The invariant that cannot be satisfied by accident: an undriven circuit
/// must stay where the DC solve left it.**
///
/// Hold the drive at a constant −3 V. The operating point is a valid solution
/// of the transient equations at every timepoint, so a correctly seeded run is
/// a flat line at the DC answer, to solver tolerance, forever. Any charge the
/// history got wrong has to flow somewhere on the first step, and there is
/// nothing else in the deck for it to hide behind.
///
/// The obvious alternative — comparing two step sizes — was tried and dropped:
/// the artefact is an RC of ~9 ps and both runs have largely recovered by the
/// time they share a sample instant, so it passed with the bug present. A test
/// that needs its probe points tuned until it fails is measuring the tuning.
#[test]
fn an_undriven_transient_stays_at_the_operating_point() {
    let src = "* PN phase shifter parked at reverse bias, nothing moving\n\
         .optical_port a\n\
         .optical_port b\n\
         Xl a fc_cw_laser power_mW=1\n\
         Vd nd 0 DC -3\n\
         Rd nd n 25\n\
         Xps a b n 0 fc_pn_ps_cap l_um=3000 v_pi_l=0.012 c_j0=750f \
           v_bi=0.917 m_j=0.5 alpha_dB_cm=2.0\n\
         .tran 1p 40p\n.end\n";
    let net = parse_spice(src).unwrap();
    let r = tran_nr_with_registry(&net, 1e-12, 40e-12, &DeviceRegistry::new())
        .expect("transient should run");
    let v = &r.node_voltages["n"];

    let v0 = v[0];
    for (k, &vk) in v.iter().enumerate() {
        assert!(
            (vk - v0).abs() < 1e-6,
            "step {k}: {vk:.9} V against the operating point {v0:.9} V. Nothing in this \
             deck moves, so every departure is charge the history did not have."
        );
    }
}
