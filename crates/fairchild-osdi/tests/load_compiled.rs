//! What a compiled OSDI library looks like from the outside.
//!
//! The descriptor walk is the first thing that can go wrong with a real model
//! and the last thing a hand-written fixture could honestly check: an imitation
//! reports whatever its author wrote, so agreeing with it proved only that two
//! files matched. These numbers come from the compiler.

use std::ffi::CStr;

use fairchild_osdi::OsdiLibrary;

mod common;

#[test]
fn a_compiled_library_reports_its_one_model() {
    let Some(path) = common::compiled("rc_shunt") else {
        return;
    };

    let lib = unsafe { OsdiLibrary::open(&path) }.expect("failed to open the compiled library");

    // 0.4 is the only version this runtime loads, and the one `openvaf-r` emits.
    assert_eq!(lib.version, (0, 4));
    assert_eq!(lib.num_descriptors, 1);

    // We read a *prefix* of the descriptor: the compiler's is currently larger
    // than ours (360 bytes against 312) while still declaring OSDI 0.4, so
    // demanding equality would reject a library that works. What must hold is
    // that everything we read is there — and that the walk strides by the
    // library's size, which `two_descriptors_are_both_read` covers. The mock
    // reported our own `size_of`, so neither question could ever have come up.
    assert!(
        lib.descriptor_size >= std::mem::size_of::<fairchild_osdi::ffi::OsdiDescriptor>(),
        "library descriptor is {} bytes, shorter than the {} we read",
        lib.descriptor_size,
        std::mem::size_of::<fairchild_osdi::ffi::OsdiDescriptor>()
    );

    let descs: Vec<_> = lib.descriptors().collect();
    assert_eq!(descs.len(), 1);

    let name = unsafe { CStr::from_ptr(descs[0].name) };
    assert_eq!(name.to_str().unwrap(), "rc_shunt");

    assert_eq!(descs[0].num_nodes, 2);
    assert_eq!(descs[0].num_terminals, 2);
    // `gd` and `c` from the source, plus the `$mfactor` every module gets
    // implicitly — and that one is the *instance* parameter, which is the half
    // of the table a model with no parameters could never have exercised.
    assert_eq!(descs[0].num_params, 3);
    assert_eq!(descs[0].num_instance_params, 1);

    // The entry points the device layer calls through.
    assert!(descs[0].setup_model.is_some());
    assert!(descs[0].setup_instance.is_some());
    assert!(descs[0].eval.is_some());
    assert!(descs[0].load_spice_rhs_dc.is_some());
    assert!(descs[0].write_jacobian_array_resist.is_some());

    // Two nodes wired to each other: a full 2×2 block, resistive and reactive.
    assert_eq!(descs[0].num_resistive_jacobian_entries, 4);
    assert_eq!(descs[0].num_reactive_jacobian_entries, 4);
}

/// Two modules in one library. The second descriptor is only in the right place
/// if the walk strides by the size the *library* reports; striding by our own
/// struct would land short of it and read a name out of the middle of the first.
#[test]
fn two_descriptors_are_both_read() {
    let Some(path) = common::compiled("two_models") else {
        return;
    };
    let lib = unsafe { OsdiLibrary::open(&path) }.expect("open");
    assert_eq!(lib.num_descriptors, 2);

    let names: Vec<String> = lib
        .descriptors()
        .map(|d| {
            unsafe { CStr::from_ptr(d.name) }
                .to_string_lossy()
                .into_owned()
        })
        .collect();
    assert!(
        names.contains(&"first_g".to_string()) && names.contains(&"second_g".to_string()),
        "both models must be readable, got {names:?}"
    );

    // And the fields past the name have to be the second model's own, not
    // whatever a wrong stride would put there.
    let second = lib
        .descriptors()
        .find(|d| unsafe { CStr::from_ptr(d.name) }.to_bytes() == b"second_g")
        .expect("second_g");
    assert_eq!(second.num_terminals, 2);
    assert_eq!(
        second.num_params, 3,
        "second_g declares gd and c, plus the implicit $mfactor"
    );
    assert_eq!(
        second.num_reactive_jacobian_entries, 4,
        "second_g has a ddt"
    );
}
