//! `--quiet` has to silence the whole library, not just the CLI's own messages.
//!
//! Most warnings a user sees come from inside the parser and the solver — a
//! skipped `.control` block, an ignored `.print`, an unrecognised `.options` key,
//! a MOSFET card asking for a level that is not implemented. Those printed
//! regardless of `--quiet` until the switch in `fairchild_parser::warn` existed,
//! so the flag documented as "suppress all warning messages" suppressed four of
//! them and left the rest.
//!
//! Driving the real binary is the only way to check this: the warnings go to
//! stderr from library code, which no in-process unit test can capture.

use std::path::PathBuf;
use std::process::Command;

/// A deck that trips one warning from each layer that can raise one.
const DECK: &str = "\
* every warning at once
V1 in 0 PULSE(0 1 0 1n 1n 1u 2u)
R1 in out 1k
C1 out 0 1n
.model nm NMOS (LEVEL=3 VTO=0.7 KP=100u)
.options reltol=1e-4 banana=7
.print tran V(out)
.print tran V(in)
.control
run
plot v(out)
.endc
.tran 1e-8 2e-6
";

fn deck_path() -> PathBuf {
    // Named per test binary and process so a concurrent `cargo test` run cannot
    // delete another's input mid-read.
    let mut p = std::env::temp_dir();
    p.push(format!("fairchild_quiet_{}.sp", std::process::id()));
    std::fs::write(&p, DECK).expect("write deck");
    p
}

fn run(args: &[&str]) -> (String, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_fairchild"))
        .args(args)
        .output()
        .expect("run fairchild");
    assert!(out.status.success(), "fairchild exited {:?}", out.status);
    (
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

#[test]
fn quiet_silences_library_warnings_without_changing_the_answer() {
    let deck = deck_path();
    let path = deck.to_str().unwrap();

    let (loud_out, loud_err) = run(&["-f", path]);
    let (quiet_out, quiet_err) = run(&["-f", path, "--quiet"]);

    // Each of these comes from a different crate or module, and each used to
    // ignore --quiet.
    for expected in [
        ".control block skipped",              // parser, collect_defs
        ".print is ignored",                   // parser, directive dispatch
        ".options 'banana' is not recognised", // core, SimOptions::from_netlist
        "LEVEL=3",                             // core, device_registry
    ] {
        assert!(
            loud_err.contains(expected),
            "expected {expected:?} on stderr without --quiet, got:\n{loud_err}"
        );
    }

    assert_eq!(
        quiet_err, "",
        "--quiet must suppress every warning, still got:\n{quiet_err}"
    );
    assert_eq!(
        loud_out, quiet_out,
        "--quiet may change stderr and nothing else"
    );

    // Two `.print` lines, one warning: the dedup is what keeps a PDK deck's forty
    // of them from burying the rest.
    assert_eq!(
        loud_err.matches(".print is ignored").count(),
        1,
        "expected one .print warning for two .print lines, got:\n{loud_err}"
    );

    let _ = std::fs::remove_file(&deck);
}
