//! Attribution harness for #42: call a model's `setup_model`/`setup_instance`
//! the way OpenVAF's own reference host does, with fairchild's device layer out
//! of the picture.
//!
//! A crash here is in the model's generated code (or in the ABI contract this
//! mirrors); a crash only through `OsdiDevice` is ours. Run against a compiled
//! module:
//!
//! ```bash
//! OSDI_ATTRIB=/path/to/model.osdi cargo test -p fairchild-osdi --test setup_attribution -- --nocapture
//! ```
//!
//! Skipped when the variable is unset, because the `.osdi` is a platform binary
//! that no checkout can carry.

use std::ffi::c_void;
use std::path::PathBuf;

use fairchild_osdi::ffi::{OsdiInitInfo, OsdiSimParas};
use fairchild_osdi::OsdiLibrary;

/// The reference host allocates model and instance state with `max_align_t`
/// alignment, so this mirrors it with a 16-byte-aligned buffer.
fn state(size: usize) -> Vec<u128> {
    vec![0u128; size.div_ceil(16).max(1)]
}

/// An empty simulator-parameter table is a pointer to a *null entry*.
fn empty_paras(slot: &mut *mut i8, slot_str: &mut *mut i8) -> OsdiSimParas {
    OsdiSimParas {
        names: slot as *mut *mut i8,
        vals: std::ptr::null_mut(),
        names_str: slot_str as *mut *mut i8,
        vals_str: std::ptr::null_mut(),
    }
}

/// A debugging harness, not a test: it calls `setup_model`/`setup_instance`
/// by hand and prints what happens. It asserts nothing, so running it in the
/// suite costs time and buys no signal — `#[ignore]` keeps it available
/// (`cargo test -- --ignored`) without pretending it covers anything.
#[test]
#[ignore = "diagnostic harness: prints, asserts nothing"]
fn setup_runs_without_the_device_layer() {
    let Ok(path) = std::env::var("OSDI_ATTRIB") else {
        eprintln!("OSDI_ATTRIB unset — skipping attribution harness");
        return;
    };
    let lib = unsafe { OsdiLibrary::open(&PathBuf::from(&path)) }.expect("load .osdi");
    let desc = lib.descriptors().next().expect("at least one descriptor");

    let mut model = state(desc.model_size as usize);
    let mut inst = state(desc.instance_size as usize);
    let handle = c"attribution".as_ptr() as *mut c_void;

    if let Some(f) = desc.setup_model {
        let (mut a, mut b) = (std::ptr::null_mut(), std::ptr::null_mut());
        let mut paras = empty_paras(&mut a, &mut b);
        let mut res = OsdiInitInfo {
            flags: 0,
            num_errors: 0,
            errors: std::ptr::null_mut(),
        };
        eprintln!("calling setup_model…");
        unsafe {
            f(
                handle,
                model.as_mut_ptr() as *mut c_void,
                &mut paras,
                &mut res,
            )
        };
        eprintln!("setup_model returned, num_errors={}", res.num_errors);
    }

    if let Some(f) = desc.setup_instance {
        let (mut a, mut b) = (std::ptr::null_mut(), std::ptr::null_mut());
        let mut paras = empty_paras(&mut a, &mut b);
        let mut res = OsdiInitInfo {
            flags: 0,
            num_errors: 0,
            errors: std::ptr::null_mut(),
        };
        eprintln!("calling setup_instance…");
        unsafe {
            f(
                handle,
                inst.as_mut_ptr() as *mut c_void,
                model.as_mut_ptr() as *mut c_void,
                300.15,
                desc.num_terminals,
                &mut paras,
                &mut res,
            )
        };
        eprintln!("setup_instance returned, num_errors={}", res.num_errors);
    }
}
