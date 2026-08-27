//! BSIM4 end to end: does a real foundry-class model conduct, and does it
//! conduct the right amount?
//!
//! `docs/user-guide.md` advertises that a `bsim4.osdi` can be dropped in and the
//! deck does not change shape. That claim had never been run: the OSDI path was
//! only ever exercised against fixtures small enough to miss three separate ABI
//! faults (#66, and `osdi_abi_contract.rs` for the faults themselves).
//!
//! # Why this test skips by default
//!
//! BSIM4.8 as Verilog-A (cogenda's VA-BSIM48, the copy in the OpenVAF-Reloaded
//! tree) is CC BY-NC-SA 4.0. NonCommercial and ShareAlike are both incompatible
//! with vendoring it into an Apache-2.0 repository, so it cannot live in
//! `tests/models/` beside the fixtures. Point `FAIRCHILD_BSIM4_VA` at a copy to
//! run it.
//!
//! The reference numbers are ngspice-46 running **the same Verilog-A source**
//! through its own OSDI loader (`pre_osdi`), not ngspice's built-in BSIM4. That
//! is deliberate: comparing against the C model would confound a fault in our
//! OSDI path with a difference between the C model and its Verilog-A
//! translation. Against the same compiled `.osdi`, any disagreement is ours.
//! (For the record the two agree to about 1% here, which is the translation.)
//!
//! Reproduce a reference row with — note `numdgt`, without which `print`
//! rounds to six significant figures and the comparison below cannot be made
//! tighter than 1e-5:
//! ```text
//! .control
//! pre_osdi <bsim4.osdi>
//! .endc
//! .model bmod bsim4va type=1 w=10u l=1u
//! N1 d g 0 0 bmod
//! Vd d 0 DC 1.0
//! Vg g 0 DC 1.8
//! .control
//! set numdgt=12
//! dc Vg 0.4 1.8 0.4
//! print i(Vd)
//! .endc
//! ```

use std::path::{Path, PathBuf};
use std::sync::Arc;

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_osdi::{OsdiLibrary, VaCompiler, VaOptions};
use fairchild_parser::parse_spice;

/// The model source, or `None` with a message saying how to supply it.
fn bsim4_va() -> Option<PathBuf> {
    let Some(raw) = std::env::var_os("FAIRCHILD_BSIM4_VA") else {
        eprintln!(
            "skipping BSIM4 acceptance: set FAIRCHILD_BSIM4_VA to a bsim4.va \
             (CC BY-NC-SA, so it cannot be vendored here)"
        );
        return None;
    };
    let path = PathBuf::from(raw);
    assert!(
        path.exists(),
        "FAIRCHILD_BSIM4_VA points at {} which does not exist",
        path.display()
    );
    Some(path)
}

/// `I(Vd)` at one bias, through the whole `.va` path.
fn i_vd(va: &Path, osdi: &Path, instance: &str, sources: &str) -> f64 {
    let deck = format!(
        "* bsim4 acceptance\n.va {}\n{instance}{sources}.op\n",
        va.display()
    );
    let netlist = parse_spice(&deck).expect("parse");
    let lib = Arc::new(unsafe { OsdiLibrary::open(osdi) }.expect("dlopen"));
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&netlist.models);
    lib.register_into(&mut registry);
    registry.register_loaded_model_cards(&netlist.models);
    let r = dc_op_nr_with_registry(&netlist, &registry).expect("DC OP failed");
    r.vsrc_current("vd").expect("I(Vd)")
}

fn compile(va: &Path) -> PathBuf {
    let opts = VaOptions::from_env();
    let compiler = VaCompiler::find(&opts).expect("no Verilog-A compiler on PATH");
    fairchild_osdi::compile::compile(&compiler, va, &opts).expect("compiling bsim4.va failed")
}

/// ngspice-46 through `pre_osdi` on the same source, printed at `numdgt=12`.
///
/// The two agree to eleven digits, which is the point: this is one compiled
/// `.osdi` evaluated by two simulators, so the only room for disagreement is
/// how each of them stamps it. A tolerance loose enough to be comfortable would
/// have accepted the 18% the packed-Jacobian fault produced.
///
/// There is no floor left to reserve for. `gmin` used to be one — fairchild put
/// it on every node's diagonal and ngspice puts it across pn junctions, so 1 pS
/// against 119 µA left a relative 8e-9 that no amount of correctness could close.
/// With `gmin` across the junctions here too, the same four points agree to
/// between 9e-14 and 6e-13: femtoamps against milliamps, which is round-off at
/// ngspice's twelve printed digits.
///
/// `1e-11` is ~20x above the worst of those and four orders below the old floor.
/// Tightening it is the point — a tolerance sized for a discrepancy that no
/// longer exists is four orders of unearned room for the next one.
const REL: f64 = 1e-11;

#[test]
fn bsim4_nmos_transfer_matches_ngspice_on_the_same_model() {
    let Some(va) = bsim4_va() else { return };
    let osdi = compile(&va);
    let inst = "X1 d g 0 0 bsim4va type=1 W=10u L=1u\n";

    // Vg from 0.4 to 1.8 V at Vd = 1.0 V: threshold, moderate and strong
    // inversion. The 0 V point is left out — there both simulators sit on their
    // own gmin floor, which is not a property of the model.
    for (vg, want) in [
        (0.4, -1.18858575113e-4),
        (0.8, -8.00063082481e-4),
        (1.2, -1.78410742995e-3),
        (1.6, -2.88802191032e-3),
    ] {
        let got = i_vd(
            &va,
            &osdi,
            inst,
            &format!("Vd d 0 DC 1.0\nVg g 0 DC {vg}\n"),
        );
        let rel = (got - want).abs() / want.abs();
        assert!(
            rel < REL,
            "Vg={vg}: got {got:.6e}, ngspice says {want:.6e} (rel {rel:.2e})"
        );
        // And it conducts at all, which is the thing #66 reported.
        assert!(got.abs() > 1e-9, "Vg={vg}: {got:.3e} is the gmin floor");
    }
}

#[test]
fn bsim4_short_channel_with_body_bias_matches_ngspice() {
    let Some(va) = bsim4_va() else { return };
    let osdi = compile(&va);
    // L = 0.18 µm and Vb = −0.5 V: short-channel effects, and a body terminal
    // that is neither the source nor ground — which only reads correctly if the
    // body network's internal nodes collapse onto the `b` terminal.
    let inst = "X1 d g 0 bb bsim4va type=1 W=2u L=0.18u\n";
    for (vd, want) in [
        (0.3, -6.11098820611e-4),
        (0.9, -8.83404114007e-4),
        (1.5, -9.78270702197e-4),
    ] {
        let got = i_vd(
            &va,
            &osdi,
            inst,
            &format!("Vd d 0 DC {vd}\nVg g 0 DC 1.2\nVb bb 0 DC -0.5\n"),
        );
        let rel = (got - want).abs() / want.abs();
        assert!(
            rel < REL,
            "Vd={vd}: got {got:.6e}, ngspice says {want:.6e} (rel {rel:.2e})"
        );
    }
}

#[test]
fn bsim4_pmos_at_85c_matches_ngspice() {
    let Some(va) = bsim4_va() else { return };
    let osdi = compile(&va);
    // `type=-1` is the negative half of the integer parameter, and `.temp 85`
    // is the temperature `setup_instance` is handed.
    let inst = "X1 d g vdd vdd bsim4va type=-1 W=4u L=0.18u\n";
    for (vg, want) in [
        (0.0, 1.685447646794e-3),
        (0.6, 8.751386361639e-4),
        (1.2, 1.720408751930e-4),
    ] {
        let got = i_vd(
            &va,
            &osdi,
            inst,
            &format!("Vdd vdd 0 DC 1.8\nVd d 0 DC 0.9\nVg g 0 DC {vg}\n.temp 85\n"),
        );
        let rel = (got - want).abs() / want.abs();
        assert!(
            rel < REL,
            "Vg={vg}: got {got:.6e}, ngspice says {want:.6e} (rel {rel:.2e})"
        );
    }
}
