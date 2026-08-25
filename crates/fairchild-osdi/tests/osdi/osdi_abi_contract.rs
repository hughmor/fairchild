//! Three corners of the OSDI ABI that no fixture in this tree used to reach, and
//! that a real foundry model reaches on its first line.
//!
//! Each of these was a silent wrong answer, and together they are why BSIM4
//! loaded, converged, and conducted nothing (#66). They are separated here
//! because each of the three fails *on its own* — verified by reverting one fix
//! at a time against a BSIM4 deck — and a single BSIM4 deck would only say
//! "wrong" without saying which.
//!
//! Every expectation is a closed form, not another simulator: `g · V` with the
//! series resistance written out by hand. Two subsystems agreeing about a
//! shared misreading of the ABI is exactly what happened before.

use std::path::Path;
use std::sync::Arc;

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_osdi::OsdiLibrary;
use fairchild_parser::parse_spice;

use super::common::compiled;

/// Solve `Ib 0 a 1mA` through `<device_line>`, and return V(a).
///
/// A current source drives the device, so V(a) reads the device's conductance
/// directly: `V = I / g`. `None` when there is no Verilog-A compiler.
fn v_at(osdi: &Path, device_line: &str) -> f64 {
    let deck = format!("* abi\nIb 0 a 1m\n{device_line}.op\n");
    let netlist = parse_spice(&deck).expect("parse");
    let lib = Arc::new(unsafe { OsdiLibrary::open(osdi) }.expect("dlopen"));
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&netlist.models);
    lib.register_into(&mut registry);
    registry.register_loaded_model_cards(&netlist.models);
    let result = dc_op_nr_with_registry(&netlist, &registry).expect("DC OP failed");
    result.node_voltage("a").expect("node a")
}

/// 1 mA through a conductance is 1 mV per mS, so every expectation here is an
/// exact rational. The tolerance covers `gmin`, which shunts every node
/// including a device-internal one: at 1 pS against a 1 mA drive that is a
/// relative 1.25e-9 on the series case, and nothing else here is inexact.
fn close(got: f64, want: f64, what: &str) {
    let rel = (got - want).abs() / want.abs();
    assert!(rel < 1e-8, "{what}: got {got:.12e}, want {want:.12e}");
}

#[test]
fn integer_parameter_reaches_the_model_at_its_own_width() {
    let Some(osdi) = compiled("abi_int_param") else {
        return;
    };
    // sel=1 → g·guard = 1e-3 · 7. 1 mA through it is 1/7 V.
    close(
        v_at(&osdi, "Xr a 0 abi_int_param sel=1 guard=7 g=1m\n"),
        1e-3 / (1e-3 * 7.0),
        "sel=1",
    );
    // A negative integer is the other half of BSIM4's `type`, and it exercises
    // the sign the f64 path lost as well as the width.
    close(
        v_at(&osdi, "Xr a 0 abi_int_param sel=-1 guard=7 g=1m\n"),
        1e-3 / (0.5e-3 * 7.0),
        "sel=-1",
    );
    // `guard` sits next to `sel` in the model struct; an 8-byte write into
    // `sel`'s 4-byte slot lands on it. Moving it must move the answer, or the
    // two checks above could pass over a clobbered neighbour.
    close(
        v_at(&osdi, "Xr a 0 abi_int_param sel=1 guard=2 g=1m\n"),
        1e-3 / (1e-3 * 2.0),
        "guard=2",
    );
}

#[test]
fn a_reactive_only_jacobian_entry_does_not_shift_the_resistive_ones() {
    let Some(osdi) = compiled("abi_jac_packing") else {
        return;
    };
    // `c` carries a capacitance and no conductance, and is numbered first, so
    // the packed resistive array starts one entry ahead of `jacobian_entries`.
    // The DC answer is still just I/g: the capacitor contributes nothing here.
    close(
        v_at(&osdi, "Xr c a 0 abi_jac_packing g=1m cval=1p\n"),
        1.0,
        "g=1m",
    );
    close(
        v_at(&osdi, "Xr c a 0 abi_jac_packing g=4m cval=1p\n"),
        0.25,
        "g=4m",
    );
}

#[test]
fn a_collapsed_internal_node_shares_its_neighbours_row() {
    let Some(osdi) = compiled("abi_collapse") else {
        return;
    };
    // rs=0: `ai` collapses into `a`, and the device is just `g`.
    close(v_at(&osdi, "Xr a 0 abi_collapse g=1m\n"), 1.0, "rs=0");
    // rs>0: no collapse, and `ai` is a real unknown — 1 mA through
    // rs + 1/g in series.
    close(
        v_at(&osdi, "Xr a 0 abi_collapse g=1m rs=250\n"),
        1e-3 * (250.0 + 1000.0),
        "rs=250",
    );
}
