//! Regression: an OSDI device must stay Newton-consistent in a transient.
//!
//! `load_residual_tran` used to hand OSDI the *previous timestep* solution as
//! `prev_solve`.  OpenVAF's `load_spice_rhs_tran` linearises about whatever
//! vector it is given (`J·prev_solve − f`), so the resistive residual and the
//! Jacobian were taken about different points.  A nonlinear model then only
//! converged while its operating point sat still: `.op` and `.dc` were fine,
//! and a constant-source `.tran` was fine, but the first moving source blew up
//! with "Newton-Raphson did not converge".
//!
//! Pre-condition: legacy/va-models/build/diode_shockley.osdi must exist.
//! Build it with:
//!   DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib \
//!   openvaf-r legacy/va-models/electronic/diode_shockley.va \
//!     --output legacy/va-models/build/diode_shockley.osdi

use std::path::PathBuf;
use std::sync::Arc;

use fairchild_core::{tran_nr_with_registry, DeviceRegistry};
use fairchild_osdi::OsdiLibrary;
use fairchild_parser::parse_spice;

const STEP: f64 = 2e-5;
const STOP: f64 = 2e-3;

/// Half-wave rectifier driven by a 1 kHz sine — deliberately swings the
/// junction from deep reverse into forward conduction on every cycle.
fn deck(diode_line: &str, model_card: &str) -> String {
    format!(
        "* sine-driven diode\n\
         Vin in 0 SIN(0 5 1k)\n\
         R1 in a 1k\n\
         {diode_line}\n\
         {model_card}\
         .tran {STEP} {STOP}\n\
         .end\n"
    )
}

#[test]
fn osdi_diode_tracks_native_diode_through_a_sine() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../legacy/va-models/build/diode_shockley.osdi");
    if !path.exists() {
        eprintln!("Skipping: {} not found — compile it first.", path.display());
        return;
    }

    let osdi_netlist = parse_spice(&deck("D1 a 0 diode_shockley", "")).unwrap();
    let lib = Arc::new(unsafe { OsdiLibrary::open(&path) }.expect("dlopen failed"));
    let mut registry = DeviceRegistry::new();
    lib.register_into(&mut registry);
    // This is the call that used to fail outright.
    let osdi = tran_nr_with_registry(&osdi_netlist, STEP, STOP, &registry)
        .expect("OSDI diode transient failed to converge");

    // fairchild's own Shockley diode with the same Is/N is the reference.
    let native_netlist =
        parse_spice(&deck("D1 a 0 dmod", ".model dmod d(is=1e-14 n=1)\n")).unwrap();
    let mut native_registry = DeviceRegistry::new();
    native_registry.register_builtin_models(&native_netlist.models);
    let native = tran_nr_with_registry(&native_netlist, STEP, STOP, &native_registry)
        .expect("native diode transient failed");

    let va = osdi.node_voltages.get("a").expect("node a missing");
    let rs = native.node_voltages.get("a").expect("node a missing");
    assert_eq!(va.len(), rs.len());

    // Both models are the same equation, so they should agree closely.  The
    // OSDI one uses Tnom=300.15 K against the native default of 300.15 K, and
    // the built-in adds pnjlim; a few mV of spread over a 5 V swing is the
    // limit of what that leaves.
    let worst = va
        .iter()
        .zip(rs)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0f64, f64::max);
    assert!(worst < 5e-3, "OSDI vs native diode differ by {worst:.3e} V");

    // And it has to actually rectify, or the comparison above is vacuous.
    let peak = va.iter().cloned().fold(f64::MIN, f64::max);
    assert!(
        (0.6..0.85).contains(&peak),
        "forward peak {peak:.4} V is not a conducting diode"
    );
}
