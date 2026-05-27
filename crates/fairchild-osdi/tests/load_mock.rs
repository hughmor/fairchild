//! Integration test: load osdi-mock cdylib and verify registry walk.
//!
//! Run via `cargo test --workspace` (builds osdi-mock first) or:
//!   cargo build -p osdi-mock && cargo test -p fairchild-osdi

use std::ffi::CStr;
use std::path::PathBuf;

use fairchild_osdi::OsdiLibrary;

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

#[test]
fn load_mock_and_list_models() {
    let path = mock_path();
    if !path.exists() {
        eprintln!(
            "osdi-mock not found at {path:?}\n\
             Run `cargo build -p osdi-mock` first, or use `cargo test --workspace`."
        );
        return;
    }

    let lib = unsafe { OsdiLibrary::open(&path) }.expect("failed to open osdi-mock");

    assert_eq!(lib.version, (0, 4));
    assert_eq!(lib.num_descriptors, 1);
    assert_eq!(
        lib.descriptor_size,
        std::mem::size_of::<fairchild_osdi::ffi::OsdiDescriptor>(),
    );

    let descs: Vec<_> = lib.descriptors().collect();
    assert_eq!(descs.len(), 1);

    let name = unsafe { CStr::from_ptr(descs[0].name) };
    assert_eq!(name.to_str().unwrap(), "test_conductance");

    assert_eq!(descs[0].num_nodes, 2);
    assert_eq!(descs[0].num_terminals, 2);
    assert_eq!(descs[0].num_params, 0);

    // All key function pointers are now real implementations.
    assert!(descs[0].setup_model.is_some());
    assert!(descs[0].eval.is_some());
    assert!(descs[0].load_spice_rhs_dc.is_some());
    assert!(descs[0].write_jacobian_array_resist.is_some());
    assert_eq!(descs[0].num_resistive_jacobian_entries, 4);
}
