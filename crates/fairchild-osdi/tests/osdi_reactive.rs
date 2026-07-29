//! A Verilog-A `ddt` charge must reach `.ac` and `.noise`, not just `.tran`.
//!
//! `OsdiDevice` reports no `small_signal_reactances`, so the frequency-domain
//! analyses used to drop its reactive Jacobian entirely — changing a model's
//! junction capacitance by six orders of magnitude left the AC sweep
//! bit-identical. `Device::load_reactive_jacobian` now hands the same ∂q/∂x
//! entries the transient companion uses to the susceptance block.
//!
//! Deliberately *not* fixed by making `OsdiDevice` report two-terminal
//! branches: a Verilog-A charge is a general matrix, and ∂q_i/∂v_j ≠ ∂q_j/∂v_i
//! (transcapacitance) is exactly what a BSIM-class model is made of.
//!
//! Runs against `osdi-mock` — 1 mS in parallel with 1 nF — so it needs no
//! OpenVAF and runs in CI.

use std::f64::consts::PI;
use std::path::PathBuf;
use std::sync::Arc;

use fairchild_core::tran::IntegratorMode;
use fairchild_core::{
    ac::ac_analysis, tran_nr_with_registry_opts, tran_nr_with_registry_var_opts, DeviceRegistry,
    SimOptions,
};
use fairchild_osdi::OsdiLibrary;
use fairchild_parser::parse_spice;
use osdi_mock::{MOCK_C, MOCK_GD};

fn mock_path() -> PathBuf {
    let mut p = std::env::current_exe().expect("current_exe");
    p.pop();
    if p.ends_with("deps") {
        p.pop();
    }
    let ext = if cfg!(target_os = "macos") {
        "dylib"
    } else {
        "so"
    };
    p.push(format!("libosdi_mock.{ext}"));
    p
}

/// `Vin —[R]— out —[mock]— gnd`, so V(out)/Vin = Y_mock⁻¹ / (R + Y_mock⁻¹)
/// with Y_mock = gd + jωC. Solving: V(out) = 1 / (1 + R·(gd + jωC)).
#[test]
fn reactive_jacobian_reaches_ac() {
    let path = mock_path();
    if !path.exists() {
        eprintln!("osdi-mock not found at {path:?}; run `cargo build -p osdi-mock`.");
        return;
    }

    const R: f64 = 1e3;
    let netlist = parse_spice(
        "* mock RC divider\n\
         Vin in 0 DC 0 AC 1\n\
         R1 in out 1k\n\
         Xm out 0 test_conductance\n\
         .ac dec 2 1k 10meg\n\
         .end\n",
    )
    .unwrap();

    let lib = Arc::new(unsafe { OsdiLibrary::open(&path) }.expect("dlopen failed"));
    let mut registry = DeviceRegistry::new();
    lib.register_into(&mut registry);

    let freqs: Vec<f64> = (0..=8).map(|k| 1e3 * 10f64.powf(k as f64 / 2.0)).collect();
    let r = ac_analysis(&netlist, &freqs, Some("Vin"), &registry).expect("AC failed");
    let out = r.voltages.get("out").expect("node out missing");

    let mut moved = false;
    for (i, &f) in freqs.iter().enumerate() {
        let omega = 2.0 * PI * f;
        // 1 / (1 + R·gd + jωRC)
        let (dr, di) = (1.0 + R * MOCK_GD, omega * R * MOCK_C);
        let denom = dr * dr + di * di;
        let (want_re, want_im) = (dr / denom, -di / denom);
        let (got_re, got_im) = out[i];

        assert!(
            (got_re - want_re).abs() < 1e-9 && (got_im - want_im).abs() < 1e-9,
            "f={f:.3e}: got ({got_re:.9}, {got_im:.9}), want ({want_re:.9}, {want_im:.9})"
        );
        // The capacitance has to actually bend the response somewhere in the
        // sweep, or matching the formula proves nothing.
        if got_im.abs() > 0.05 {
            moved = true;
        }
    }
    assert!(moved, "the sweep never left the resistive limit");
}

/// A Verilog-A `ddt` charge must integrate with the *configured* method, not
/// always Backward Euler.
///
/// `tran.rs` and `tran_step.rs` hand `load_*_tran` an `alpha = 1/h`, which can
/// only express BE — Trapezoidal and BDF-2 need history terms a single scalar
/// has nowhere to put. So a Verilog-A `ddt(C*V)` was a different circuit
/// element from a discrete `C` of the same value under the default method
/// (Trapezoidal), by about 0.6 % on a 0.45 τ step.
///
/// osdi-mock is 1 mS in parallel with 1 nF, which is *exactly* a discrete
/// `Rg = 1k` plus `Cg = 1n`. Comparing the two under every method is the
/// strongest available statement, and needs no OpenVAF.
#[test]
fn ddt_honours_the_integration_method() {
    let path = mock_path();
    if !path.exists() {
        eprintln!("osdi-mock not found at {path:?}; run `cargo build -p osdi-mock`.");
        return;
    }

    // 1k into the device, and the identical 1k into its discrete equivalent.
    // tau = 1k‖1k · 1n = 500 ns, stepped at 100 ns so the method is visible.
    let netlist = parse_spice(
        "* mock reactance vs the discrete equivalent\n\
         Vin in 0 PULSE(0 1 0 1n 1n 10u 20u)\n\
         R1 in a 1k\n\
         Xm a 0 test_conductance\n\
         R2 in b 1k\n\
         Rg b 0 1k\n\
         Cg b 0 1n\n\
         .tran 100n 5u\n\
         .end\n",
    )
    .unwrap();

    let lib = Arc::new(unsafe { OsdiLibrary::open(&path) }.expect("dlopen failed"));
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&netlist.models);
    lib.register_into(&mut registry);

    for mode in [
        IntegratorMode::BackwardEuler,
        IntegratorMode::Trapezoidal,
        IntegratorMode::Gear,
    ] {
        for variable_step in [false, true] {
            let mut opts = SimOptions::from_netlist(&netlist);
            opts.method = mode;
            opts.variable_step = variable_step;
            let r = if variable_step {
                tran_nr_with_registry_var_opts(&netlist, 100e-9, 5e-6, &registry, &opts)
            } else {
                tran_nr_with_registry_opts(&netlist, 100e-9, 5e-6, &registry, &opts)
            }
            .expect("transient failed");

            let va = r.node_voltages.get("a").expect("node a");
            let vb = r.node_voltages.get("b").expect("node b");
            let worst = va
                .iter()
                .zip(vb)
                .map(|(x, y)| (x - y).abs())
                .fold(0.0f64, f64::max);
            assert!(
                worst < 1e-12,
                "{mode:?} variable_step={variable_step}: Verilog-A ddt and a discrete C \
                 differ by {worst:.3e} V — the device is not using the configured method"
            );
            // The comparison is only meaningful if the reactance did something.
            let swing = va.iter().cloned().fold(f64::MIN, f64::max)
                - va.iter().cloned().fold(f64::MAX, f64::min);
            assert!(
                swing > 0.1,
                "{mode:?}: no dynamics to compare ({swing:.3e} V)"
            );
        }
    }
}
