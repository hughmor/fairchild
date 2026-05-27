/// Load the real diode_shockley.osdi compiled from Verilog-A and sanity-check its descriptor.
///
/// Pre-condition: run
///   DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib \
///   /path/to/openvaf-r legacy/va-models/diode_shockley.va \
///     --output legacy/va-models/build/diode_shockley.osdi
/// or the test will skip with a helpful message.

use std::ffi::CStr;
use std::path::PathBuf;

use fairchild_osdi::OsdiLibrary;

fn osdi_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../legacy/va-models/build/diode_shockley.osdi")
}

#[test]
fn load_shockley_osdi_descriptor() {
    let path = osdi_path();
    if !path.exists() {
        eprintln!(
            "diode_shockley.osdi not found at {}\n\
             Compile it first:\n\
             DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib \\\n\
             /path/to/openvaf-r legacy/va-models/diode_shockley.va \\\n\
             --output legacy/va-models/build/diode_shockley.osdi",
            path.display()
        );
        return;
    }

    let lib = unsafe { OsdiLibrary::open(&path) }.expect("dlopen failed");

    assert_eq!(lib.version, (0, 4), "expected OSDI v0.4");
    assert_eq!(lib.num_descriptors, 1);

    let descs: Vec<_> = lib.descriptors().collect();
    let d = descs[0];

    let name = unsafe { CStr::from_ptr(d.name) }.to_str().unwrap();
    assert_eq!(name, "diode_shockley");

    // The VA model has 2 terminals (anode, cathode) and no internal nodes.
    assert_eq!(d.num_terminals, 2);
    assert_eq!(d.num_nodes, 2);

    // 3 declared parameters (Is, N, Tnom) + OpenVAF may add implicit temperature params.
    assert!(d.num_params >= 3, "expected at least 3 params, got {}", d.num_params);

    // All required function pointers must be populated by OpenVAF.
    assert!(d.setup_model.is_some(), "setup_model missing");
    assert!(d.setup_instance.is_some(), "setup_instance missing");
    assert!(d.eval.is_some(), "eval missing");
    assert!(d.load_spice_rhs_dc.is_some(), "load_spice_rhs_dc missing");
    assert!(d.write_jacobian_array_resist.is_some(), "write_jacobian_array_resist missing");

    // The Shockley diode has a 2x2 Jacobian (I_d flows between anode and cathode):
    // entries: (A,A), (A,C), (C,A), (C,C)  → 4 resistive entries.
    assert_eq!(d.num_resistive_jacobian_entries, 4);
}
