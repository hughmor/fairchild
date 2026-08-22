//! A model that emits a diagnostic must still simulate.
//!
//! On macOS/aarch64 it did not: OpenVAF miscompiles every formatted-output and
//! severity task, so the first thing a model said took the process down —
//! `SIGSEGV`, no output, exit 139 (#42). That ruled out every real compact
//! model, all of which talk about their own parameters, which is why the
//! hand-written fixtures in this directory never showed it.
//!
//! `crate::portability` removes those calls from the source before compiling
//! wherever the toolchain cannot be trusted with them. The *transform* is unit
//! tested against explicit platform strings (so the macOS rule is covered from
//! Linux too); this is the end-to-end half: compile a model that talks, solve
//! with it, and check the answer against Ohm's law.
//!
//! On a platform with no rule this passes without any transformation, which is
//! the right shape — the assertion is about the model working, not about the
//! workaround running.

use std::sync::Arc;

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_osdi::OsdiLibrary;
use fairchild_parser::parse_spice;

mod common;

#[test]
fn a_model_that_emits_a_diagnostic_still_solves() {
    let Some(path) = common::compiled("talkative") else {
        return;
    };

    let lib = Arc::new(unsafe { OsdiLibrary::open(&path) }.expect("dlopen"));
    let mut reg = DeviceRegistry::new();
    lib.register_into(&mut reg);

    // 1 V across gd = 1 mS: the source carries 1 mA. `gd` on the instance line
    // also proves the parameter reached the model — a stripped diagnostic must
    // not disturb parameter application, which is the half of `setup_instance`
    // the removed calls sat in the middle of.
    let net = parse_spice("* talks\nV1 a 0 DC 1\nX1 a 0 talkative gd=1m\n.op\n").expect("parse");
    let res = dc_op_nr_with_registry(&net, &reg).expect("DC OP converges");
    let i = res.vsrc_current("v1").expect("source current");
    assert!(
        (i + 1e-3).abs() < 1e-9,
        "1 V across 1 mS must draw 1 mA, got {i:.6e} — the model ran but did not \
         compute (removing a diagnostic must not touch the physics)"
    );
}
