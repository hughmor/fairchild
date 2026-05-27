/// End-to-end test: load diode_shockley.osdi, register into DeviceRegistry,
/// run dc_op_nr_with_registry, compare V(b) against the analytical answer.
///
/// Circuit: Ib=1mA → node "b" → D1 (anode b, cathode 0)
/// Expected: V(b) = Vt * ln(Ib/Is + 1) with Is=1e-14, N=1, T=300.15 K
///
/// Pre-condition: legacy/va-models/build/diode_shockley.osdi must exist.
/// Build it with:
///   DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib \
///   /path/to/openvaf-r legacy/va-models/diode_shockley.va \
///     --output legacy/va-models/build/diode_shockley.osdi

use std::path::PathBuf;
use std::sync::Arc;

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_osdi::OsdiLibrary;
use fairchild_parser::parse_spice;

fn osdi_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../legacy/va-models/build/diode_shockley.osdi")
}

/// Analytical forward-voltage for the Shockley diode: V = Vt * ln(I/Is + 1).
fn analytical_vf(i_a: f64, is: f64, n: f64, temp_k: f64) -> f64 {
    let vt = 1.380649e-23 * temp_k / 1.602176634e-19;
    n * vt * (i_a / is + 1.0).ln()
}

#[test]
fn osdi_diode_current_source_bias() {
    let path = osdi_path();
    if !path.exists() {
        eprintln!(
            "Skipping: {} not found.\n\
             Run: DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib \\\n\
             openvaf-r legacy/va-models/diode_shockley.va \\\n\
             --output legacy/va-models/build/diode_shockley.osdi",
            path.display()
        );
        return;
    }

    // Netlist: 1 mA current source into node b, OSDI diode from b to GND.
    // The model name matches the OSDI descriptor name "diode_shockley".
    let netlist = parse_spice(
        "* OSDI diode DC bias\n\
         Ib 0 b 1m\n\
         D1 b 0 diode_shockley\n\
         .op\n\
         .end\n",
    )
    .unwrap();

    let lib = Arc::new(unsafe { OsdiLibrary::open(&path) }.expect("dlopen failed"));

    let mut registry = DeviceRegistry::new();
    lib.register_into(&mut registry);

    let result = dc_op_nr_with_registry(&netlist, &registry).expect("DC OP failed");

    let vb = result.node_voltage("b").unwrap();

    // Shockley diode defaults: Is=1e-14, N=1, Tnom=300.15 K.
    let expected = analytical_vf(1e-3, 1e-14, 1.0, 300.15);
    let tol = 2e-3; // 2 mV — OSDI eval may use slightly different Vt/temperature

    assert!(
        (vb - expected).abs() < tol,
        "V(b) = {vb:.6}  expected ≈ {expected:.6}  diff = {:.2e}",
        (vb - expected).abs()
    );

    println!("V(b) = {vb:.6} V  (expected {expected:.6} V)");
}

#[test]
fn osdi_diode_series_resistor() {
    let path = osdi_path();
    if !path.exists() { return; }

    // Vdd=5V, R=10k in series with OSDI Shockley diode.
    // V(b) should be in the forward-bias range ~0.55–0.75 V.
    let netlist = parse_spice(
        "* OSDI series R-D\n\
         Vdd a 0 DC 5\n\
         R1 a b 10k\n\
         D1 b 0 diode_shockley\n\
         .op\n\
         .end\n",
    )
    .unwrap();

    let lib = Arc::new(unsafe { OsdiLibrary::open(&path) }.expect("dlopen"));
    let mut registry = DeviceRegistry::new();
    lib.register_into(&mut registry);

    let result = dc_op_nr_with_registry(&netlist, &registry).expect("DC OP failed");
    let vb = result.node_voltage("b").unwrap();

    assert!(vb > 0.50 && vb < 0.80, "V(b) = {vb:.4} V — out of forward-bias range");

    // KCL check: (Vdd - V(b)) / R ≈ Is * (exp(V(b)/Vt) - 1)
    let vt = 1.380649e-23 * 300.15 / 1.602176634e-19;
    let ir = (5.0 - vb) / 10e3;
    let id = 1e-14 * ((vb / vt).exp() - 1.0);
    assert!((ir - id).abs() / ir < 0.01, "KCL error > 1 %: ir={ir:.4e} id={id:.4e}");

    println!("V(b) = {vb:.6} V  KCL residual = {:.2e}", (ir - id).abs());
}
