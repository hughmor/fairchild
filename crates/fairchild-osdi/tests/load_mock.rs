//! Integration test: load osdi-mock cdylib and verify registry walk.
//!
//! Run via `cargo test --workspace` (builds osdi-mock first) or:
//!   cargo build -p osdi-mock && cargo test -p fairchild-osdi

use std::ffi::CStr;
use std::path::PathBuf;

use fairchild_osdi::OsdiLibrary;

/// Locate the osdi-mock shared library in the Cargo target directory.
/// The test binary lives at target/{profile}/deps/<name>-<hash>;
/// the cdylib lives at target/{profile}/libosdi_mock.{dylib|so}.
fn mock_path() -> PathBuf {
    let mut p = std::env::current_exe().expect("current_exe");
    p.pop(); // remove test binary name
    if p.ends_with("deps") {
        p.pop(); // step up from deps/ to profile dir
    }
    let ext = if cfg!(target_os = "macos") { "dylib" } else { "so" };
    p.push(format!("libosdi_mock.{ext}"));
    p
}

#[test]
fn load_mock_and_list_models() {
    let path = mock_path();
    if !path.exists() {
        eprintln!(
            "osdi-mock not found at {path:?}\n\
             Run `cargo build -p osdi-mock` first, or use `cargo test --workspace`."
        );
        return; // skip rather than fail when mock is absent
    }

    let lib = unsafe { OsdiLibrary::open(&path) }.expect("failed to open osdi-mock");

    assert_eq!(lib.version, (0, 4), "OSDI version should be 0.4");
    assert_eq!(lib.num_descriptors, 1, "mock exports one descriptor");
    assert_eq!(
        lib.descriptor_size,
        std::mem::size_of::<fairchild_osdi::ffi::OsdiDescriptor>(),
        "OSDI_DESCRIPTOR_SIZE should match sizeof(OsdiDescriptor)"
    );

    let descs: Vec<_> = lib.descriptors().collect();
    assert_eq!(descs.len(), 1);

    let name = unsafe { CStr::from_ptr(descs[0].name) };
    assert_eq!(name.to_str().unwrap(), "test_diode");

    assert_eq!(descs[0].num_nodes, 2);
    assert_eq!(descs[0].num_terminals, 2);
    assert_eq!(descs[0].num_params, 1);

    // Function pointers are all null in the mock
    assert!(descs[0].eval.is_none());
    assert!(descs[0].setup_model.is_none());
}
