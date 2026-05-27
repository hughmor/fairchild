/// Validates that OpenVAF-compiled optical discipline models load via OSDI.
/// Confirms no compiler fork is needed — custom natures/disciplines are transparent to OSDI.
///
/// Pre-condition: legacy/va-models/build/test_optical_discipline.osdi must exist.
/// Build it with:
///   DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib \
///   openvaf-r legacy/va-models/test_optical_discipline.va \
///     --output legacy/va-models/build/test_optical_discipline.osdi
use std::ffi::CStr;
use std::path::PathBuf;

use fairchild_osdi::OsdiLibrary;

fn osdi_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../legacy/va-models/build/test_optical_discipline.osdi")
}

#[test]
fn load_optical_discipline_osdi() {
    let path = osdi_path();
    if !path.exists() {
        eprintln!(
            "Skipping: {} not found.\n\
             Run: DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib \\\n\
             openvaf-r legacy/va-models/test_optical_discipline.va \\\n\
             --output legacy/va-models/build/test_optical_discipline.osdi",
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
    assert_eq!(name, "optical_wire");

    // 2 terminals (p, n), no internal nodes.
    assert_eq!(d.num_terminals, 2);
    assert!(d.setup_model.is_some(), "setup_model missing");
    assert!(d.setup_instance.is_some(), "setup_instance missing");
    assert!(d.eval.is_some(), "eval missing");

    println!(
        "optical_wire OSDI loaded OK: {} terminals, {} params",
        d.num_terminals, d.num_params
    );
}
