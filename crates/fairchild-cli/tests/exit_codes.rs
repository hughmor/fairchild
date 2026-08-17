//! `run` returns exit codes; it never terminates the process.
//!
//! This is not a style preference. The same function backs the `fairchild`
//! command inside the Python wheel, where `std::process::exit` would take the
//! interpreter down mid-script with no traceback and no flushed output.
//!
//! These tests are also the guard against that regressing, and they work
//! because of how they fail: a `process::exit` reintroduced anywhere below the
//! `run` boundary kills the *test harness* on the first offending case, so the
//! whole suite reports failure rather than one assertion. Verify by sabotage —
//! put an `exit(1)` back in any error arm and `cargo test -p fairchild-cli`
//! dies instead of printing a diff.

use std::io::Write;

/// Exit code for a usage error, as clap and the shell convention agree.
const USAGE: i32 = 2;
/// Exit code for "ran, and failed".
const FAILURE: i32 = 1;

fn deck(body: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::Builder::new()
        .suffix(".sp")
        .tempfile()
        .expect("tempfile");
    f.write_all(body.as_bytes()).expect("write deck");
    f.flush().expect("flush deck");
    f
}

const RC: &str =
    "* rc\nV1 in 0 PULSE(0 1 0 1n 1n 1u 2u)\nR1 in out 1k\nC1 out 0 1n\n.tran 1e-8 1e-7\n.end\n";

#[test]
fn help_and_version_are_not_errors() {
    assert_eq!(fairchild_cli::run(["fairchild", "--help"]), 0);
    assert_eq!(fairchild_cli::run(["fairchild", "--version"]), 0);
}

#[test]
fn malformed_command_line_is_a_usage_error() {
    assert_eq!(fairchild_cli::run(["fairchild", "--no-such-flag"]), USAGE);
    // --file is required, and omitting it is a usage error rather than a
    // simulation failure.
    assert_eq!(fairchild_cli::run(["fairchild"]), USAGE);
}

#[test]
fn unreadable_netlist_fails_without_exiting() {
    let code = fairchild_cli::run(["fairchild", "--file", "/nonexistent/definitely/not/here.sp"]);
    assert_eq!(code, FAILURE);
}

#[test]
fn unwritable_output_fails_without_exiting() {
    let f = deck(RC);
    let code = fairchild_cli::run([
        "fairchild",
        "--file",
        f.path().to_str().unwrap(),
        "--output",
        "/nonexistent/definitely/not/here.csv",
    ]);
    assert_eq!(code, FAILURE);
}

#[test]
fn a_deck_that_simulates_returns_zero() {
    let f = deck(RC);
    let out = tempfile::Builder::new().suffix(".csv").tempfile().unwrap();
    let code = fairchild_cli::run([
        "fairchild",
        "--file",
        f.path().to_str().unwrap(),
        "--output",
        out.path().to_str().unwrap(),
    ]);
    assert_eq!(code, 0);
    let csv = std::fs::read_to_string(out.path()).unwrap();
    assert!(csv.contains("V(out)"), "no V(out) column in:\n{csv}");
}

#[test]
fn check_and_list_modes_return_zero() {
    let f = deck(RC);
    let path = f.path().to_str().unwrap();
    // Each of these used to be its own `process::exit(0)`.
    assert_eq!(
        fairchild_cli::run(["fairchild", "--file", path, "--check"]),
        0
    );
    assert_eq!(
        fairchild_cli::run(["fairchild", "--file", path, "--list-nodes"]),
        0
    );
    assert_eq!(
        fairchild_cli::run(["fairchild", "--file", path, "--list-models"]),
        0
    );
}

#[test]
fn a_failing_analysis_returns_failure_not_a_dead_process() {
    // A floating node: parses and discipline-checks clean, then fails inside
    // the solver — the arm furthest from the `run` boundary, and the one whose
    // `process::exit` was hardest to reach from a test.
    let f = deck("* floating\nI1 a 0 1\n.op\n.end\n");
    let code = fairchild_cli::run(["fairchild", "--file", f.path().to_str().unwrap()]);
    assert_eq!(code, FAILURE);
}
