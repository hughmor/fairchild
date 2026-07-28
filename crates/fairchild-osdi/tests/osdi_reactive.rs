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

use fairchild_core::{ac::ac_analysis, DeviceRegistry};
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
