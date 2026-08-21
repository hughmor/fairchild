//! Compiling a test model with the real Verilog-A compiler.
//!
//! These tests used to run against `osdi-mock`, a hand-written cdylib that
//! imitated an OSDI library. It was retired: everything it imitated, a real
//! compiler emits, and the imitation had to be trusted in exactly the place the
//! tests were meant to check. It also grew its own bug — it cached the
//! `prev_solve` pointer past the call that owned it and read freed memory,
//! failing about one run in ten.
//!
//! So the fixtures are Verilog-A source (`tests/models/*.va`) compiled here.
//! Without a compiler on PATH the tests **skip**, the way the ngspice goldens
//! do; CI installs one and asserts it is there, so the skip cannot go quiet.

// Every test binary here includes this module and uses only the part it needs,
// so anything in it is dead code from somebody's point of view.
#![allow(dead_code)]

use std::path::PathBuf;

use fairchild_osdi::{VaCompiler, VaOptions};

/// Compile `tests/models/<name>.va` and return the `.osdi`, or `None` when no
/// Verilog-A compiler is installed.
///
/// The cache lives under the target directory rather than the user's real cache,
/// so a test run cannot be satisfied by an artefact some other run left behind —
/// and so `cargo clean` clears it.
pub fn compiled(name: &str) -> Option<PathBuf> {
    let src = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/models")
        .join(format!("{name}.va"));
    assert!(src.exists(), "missing test fixture {}", src.display());

    let opts = VaOptions {
        cache_dir: Some(cache_dir()),
        ..VaOptions::from_env()
    };
    let compiler = match VaCompiler::find(&opts) {
        Ok(c) => c,
        // The error already names what was looked for and both ways to point at
        // a compiler, so it is the whole message.
        Err(e) => {
            eprintln!("skipping: {e}");
            return None;
        }
    };
    match fairchild_osdi::compile::compile(&compiler, &src, &opts) {
        Ok(osdi) => Some(osdi),
        Err(e) => panic!("compiling {} failed: {e}", src.display()),
    }
}

/// `target/va-test-cache`, found by walking out of the test executable's
/// directory rather than guessing — `CARGO_TARGET_DIR` may point anywhere.
/// Whether a Verilog-A compiler is installed, without compiling anything.
///
/// `compiled` cannot serve as this probe for a bundle-dialect source: that
/// source is a template and is not valid Verilog-A until it is expanded, so
/// asking `compiled` about it fails for the wrong reason.
pub fn have_compiler() -> bool {
    let opts = VaOptions {
        cache_dir: Some(cache_dir()),
        ..VaOptions::from_env()
    };
    match VaCompiler::find(&opts) {
        Ok(_) => true,
        Err(e) => {
            eprintln!("skipping: {e}");
            false
        }
    }
}

pub fn test_cache_dir() -> PathBuf {
    cache_dir()
}

fn cache_dir() -> PathBuf {
    let mut p = std::env::current_exe().expect("current_exe");
    p.pop(); // the test binary's directory (…/target/debug/deps)
    if p.ends_with("deps") {
        p.pop();
    }
    p.push("va-test-cache");
    p
}
