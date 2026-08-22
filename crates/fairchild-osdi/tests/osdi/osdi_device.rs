//! The `OsdiDevice` adapter, against a compiled model.
//!
//! Wraps a real descriptor and checks the three things the adapter owes the
//! solver: the conductance reaches the matrix, a grounded terminal drops out of
//! it, and a parameter written from the deck actually changes the stamp.

use fairchild_core::device::{Device, EvalFlags, SimContext};
use fairchild_core::mna::MnaMatrix;
use fairchild_osdi::{OsdiDevice, OsdiLibrary};

use crate::common;

/// The model's own defaults: `rc_shunt.va` declares 1 mS ∥ 1 nF.
const GD: f64 = 1e-3;
const TOL: f64 = 1e-12;

fn device(path: &std::path::Path) -> (std::sync::Arc<OsdiLibrary>, OsdiDevice) {
    let lib = std::sync::Arc::new(unsafe { OsdiLibrary::open(path) }.expect("open"));
    let desc = lib.descriptors().next().expect("one descriptor") as *const _;
    // SAFETY: the descriptor lives as long as `lib`, which the caller keeps.
    let dev = unsafe { OsdiDevice::new(desc) };
    (lib, dev)
}

#[test]
fn the_conductance_reaches_the_matrix() {
    let Some(path) = common::compiled("rc_shunt") else {
        return;
    };
    let (_lib, mut dev) = device(&path);
    let ctx = SimContext::default();

    dev.setup_model(&ctx);
    dev.setup_instance(&[Some(0), Some(1)], &ctx);
    dev.eval(&[1.0, 0.0], EvalFlags::dc(), &ctx);

    let mut mat = MnaMatrix::zeros(2);
    dev.load_residual(&mut mat.b);
    assert_eq!(
        mat.b,
        vec![0.0, 0.0],
        "a linear element linearised about its own operating point has no residual"
    );

    dev.load_jacobian(&mut mat);
    assert!((mat.a[0][0] - GD).abs() < TOL, "a[0][0] = {}", mat.a[0][0]);
    assert!((mat.a[0][1] + GD).abs() < TOL, "a[0][1] = {}", mat.a[0][1]);
    assert!((mat.a[1][0] + GD).abs() < TOL, "a[1][0] = {}", mat.a[1][0]);
    assert!((mat.a[1][1] - GD).abs() < TOL, "a[1][1] = {}", mat.a[1][1]);
}

#[test]
fn a_grounded_terminal_leaves_one_row() {
    let Some(path) = common::compiled("rc_shunt") else {
        return;
    };
    let (_lib, mut dev) = device(&path);
    let ctx = SimContext::default();

    dev.setup_model(&ctx);
    dev.setup_instance(&[Some(0), None], &ctx); // second terminal to ground
    dev.eval(&[1.0], EvalFlags::dc(), &ctx);

    let mut mat = MnaMatrix::zeros(1);
    dev.load_residual(&mut mat.b);
    dev.load_jacobian(&mut mat);

    assert!((mat.a[0][0] - GD).abs() < TOL, "a[0][0] = {}", mat.a[0][0]);
    assert_eq!(mat.b, vec![0.0]);
}

/// A parameter the deck sets has to change the stamp — the check the old
/// fixture could not make, because it declared no parameters.
///
/// `gd` is a model parameter and `$mfactor` an instance one, and the two live in
/// different halves of the descriptor's table with different access rules. An
/// instance parameter used to be dropped (the id carried a `PARA_KIND` bit that
/// matched no case, and the access flag that says "look in the instance" was
/// missing), so `W`/`L` on a compiled MOSFET silently ran at their defaults.
#[test]
fn a_parameter_written_from_the_deck_changes_the_stamp() {
    let Some(path) = common::compiled("rc_shunt") else {
        return;
    };
    let ctx = SimContext::default();

    for (name, value, want) in [
        ("gd", 4e-3, 4e-3),          // model parameter
        ("$mfactor", 4.0, 4.0 * GD), // instance parameter, by its Verilog-A name
        ("m", 3.0, 3.0 * GD),        // and by the spelling a SPICE deck uses
    ] {
        let (_lib, mut dev) = device(&path);
        dev.setup_model(&ctx);
        dev.setup_instance(&[Some(0), None], &ctx);
        assert!(
            dev.set_real_param(name, value),
            "'{name}' was not applied at all"
        );
        dev.eval(&[1.0], EvalFlags::dc(), &ctx);

        let mut mat = MnaMatrix::zeros(1);
        dev.load_jacobian(&mut mat);
        assert!(
            (mat.a[0][0] - want).abs() < TOL,
            "{name}={value}: stamped {} S, want {want} S",
            mat.a[0][0]
        );
    }
}
