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

#[test]
fn a_mixed_case_module_name_is_reachable_from_a_deck() {
    let Some(osdi) = compiled("AbiMixedCase") else {
        return;
    };
    // Verilog-A preserves the case the author wrote; SPICE does not
    // distinguish. So all three spellings must resolve, and the registry is
    // where the folding has to happen — it is the one place that sees both the
    // registrar's name and the deck's.
    for line in [
        "Xr a 0 AbiMixedCase G=1m\n",
        "Xr a 0 abimixedcase G=1m\n",
        "Xr a 0 ABIMIXEDCASE G=1m\n",
    ] {
        close(v_at(&osdi, line), 1.0, line.trim());
    }
    // …and through a `.model` card, which is the foundry idiom and the path
    // `PSP102VA` and friends arrive on.
    close(
        v_at(&osdi, "Xr a 0 mycard\n.model mycard AbiMixedCase (G=4m)\n"),
        0.25,
        "card",
    );
}

/// `$simparam` reaches the model, and carries this simulator's values.
///
/// # What this catches
///
/// `OsdiSimParas` was passed as an empty (if terminated) list, so every
/// `$simparam("gmin", <default>)` in a foundry model took the default written
/// into its own call. Same fault the native models had with a `const GMIN`: the
/// option exists, the deck sets it, and nothing reads it.
///
/// It is invisible in a real model, where `gmin` is one term among many and worth
/// a picoamp. So the fixture makes each lookup the *only* term on its own branch,
/// with a fallback nothing like the real value, and the probe reads the
/// conductance straight off `V = I/g`. "Table missing" and "table present" then
/// differ by orders of magnitude instead of by a tolerance.
///
/// # How the first version of this test failed
///
/// Both branches were on one node, so the total conductance was
/// `simparam(gmin) + simparam(scale)` and `scale = 1` buried `gmin = 1e-6` under
/// a 1e-6 perturbation. It passed with the table sabotaged out of *both* the eval
/// path and the setup path. Separate node pairs, one probe each.
///
/// Two names rather than one, because a table with a broken terminator or a
/// mismatched `vals` array can still answer the first lookup correctly.
#[test]
fn simparam_carries_the_simulators_own_values() {
    let Some(osdi) = compiled("abi_simparam") else {
        return;
    };
    for gmin in [1e-3, 1e-5] {
        let deck = format!(
            "* simparam\n.options gmin={gmin:e}\n\
             Ib 0 a 1m\nIc 0 c 1m\n\
             Xs a 0 c 0 abi_simparam\n.op\n"
        );
        let netlist = parse_spice(&deck).expect("parse");
        let lib = Arc::new(unsafe { OsdiLibrary::open(&osdi) }.expect("dlopen"));
        let mut registry = DeviceRegistry::new();
        registry.register_builtin_models(&netlist.models);
        lib.register_into(&mut registry);
        registry.register_loaded_model_cards(&netlist.models);
        let r = dc_op_nr_with_registry(&netlist, &registry).expect("DC OP failed");

        // V(a) = 1 mA / simparam("gmin"). The 1.0 S fallback would read 1 mV.
        let va = r.node_voltage("a").expect("node a");
        let want_a = 1e-3 / gmin;
        let rel_a = (va - want_a).abs() / want_a;
        assert!(
            rel_a < 1e-6,
            "gmin={gmin:e}: V(a) is {va:.6e} V and 1 mA into simparam(\"gmin\") is \
             {want_a:.6e} V (rel {rel_a:.2e}). 1e-3 V here is the model's own 1.0 S \
             fallback, which is what an empty sim-para table produces."
        );

        // V(c) = 1 mA / simparam("scale"), and `scale` is 1.0. The 1e-9 fallback
        // would read 1e6 V.
        let vc = r.node_voltage("c").expect("node c");
        let rel_c = (vc - 1e-3).abs() / 1e-3;
        assert!(
            rel_c < 1e-6,
            "V(c) is {vc:.6e} V and 1 mA into simparam(\"scale\") = 1.0 S is 1e-3 V \
             (rel {rel_c:.2e}). 1e6 V is the 1e-9 fallback, which means the table \
             answered the first name and not the second — a terminator or a `vals` \
             length that does not match `names`."
        );
    }
}

/// `$bound_step` limits the timestep, and only when it binds.
///
/// # What this catches
///
/// The request was ignored outright, which left the LTE controller as the only
/// thing watching the step size — and LTE measures the error of a step already
/// taken, so a model saying "the *next* step will miss something" had no way to
/// be heard.
///
/// Reading it is a raw `f64` load at an offset a shared library hands over, and
/// the sentinel for "never calls `$bound_step`" is `u32::MAX` rather than zero.
/// Measured: every other fixture here reports 0xffffffff, this one reports a
/// valid 8-aligned 104 in a 112-byte instance. The first version guarded on zero,
/// dereferenced 0xffffffff, and aborted the process on an unrelated test.
///
/// # Why two runs
///
/// A bound that binds and a bound that does not. Asserting only the first cannot
/// tell "the request is honoured" from "the step is always this small", which is
/// the same shape as testing a feature only in its on state.
#[test]
fn bound_step_limits_the_timestep_when_it_binds() {
    let Some(osdi) = compiled("abi_bound_step") else {
        return;
    };

    let spacings = |dtmax: f64| -> (f64, usize) {
        // An RC whose own time constant is 1 us, run for 20 us with a print step
        // far larger than the bound, on the variable-step path where the model's
        // request is consulted.
        let deck = format!(
            "* bound_step\n.options variable_step=1\n\
             V1 in 0 PULSE(0 1 0 1n 1n 1m 2m)\n\
             R1 in n 1k\nC1 n 0 1n\n\
             Xb n 0 abi_bound_step g=1m dtmax={dtmax:e}\n\
             .tran 5u 20u\n"
        );
        let netlist = parse_spice(&deck).expect("parse");
        let lib = Arc::new(unsafe { OsdiLibrary::open(&osdi) }.expect("dlopen"));
        let mut registry = DeviceRegistry::new();
        registry.register_builtin_models(&netlist.models);
        lib.register_into(&mut registry);
        registry.register_loaded_model_cards(&netlist.models);
        let r = fairchild_core::tran_nr_with_registry_var(&netlist, 5e-6, 20e-6, &registry)
            .expect("transient");
        let t = &r.time;
        let worst = t.windows(2).map(|w| w[1] - w[0]).fold(0.0_f64, f64::max);
        (worst, t.len())
    };

    // Binding: 100 ns is well inside what a 1 us RC over 20 us would otherwise
    // take, so the largest gap has to come down to it.
    let (bound_gap, bound_n) = spacings(1e-7);
    assert!(
        bound_gap <= 1.02e-7,
        "with $bound_step(100n) the largest gap between timepoints is \
         {bound_gap:.4e} s. The model asked for 1e-7 and was not heard."
    );

    // Not binding: 1 ms is larger than the whole run, so it must change nothing —
    // otherwise this test is measuring "steps are small" and not "the request is
    // honoured".
    let (loose_gap, loose_n) = spacings(1e-3);
    assert!(
        loose_gap > 2e-7,
        "with $bound_step(1m) — larger than the 20 us run — the step should be \
         set by LTE alone, but the largest gap is {loose_gap:.4e} s. If this is \
         also ~1e-7 the bound is not what produced the first result."
    );
    assert!(
        bound_n > loose_n,
        "a binding bound must cost timepoints: {bound_n} with it against \
         {loose_n} without"
    );
}
