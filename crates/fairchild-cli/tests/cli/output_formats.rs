//! `--format binary`, and what `--probe` does to a streamed transient.
//!
//! All of this only exists once a process runs, and none of it is reachable
//! from `fairchild-core`'s own tests:
//!
//! * the **file** path patches a reserved point count by seeking inside a real
//!   `BufWriter<File>`, where the core tests use an in-memory `Cursor`;
//! * the **standard output** path cannot seek, so it assembles the rawfile in
//!   memory instead, and the two must produce the same file;
//! * `--probe` on a `.tran` narrows the rawfile itself, and a deck carrying a
//!   `.measure` takes a different route to the same writer.
//!
//! That last one is here because it was wrong. The `.measure` path wrote the
//! collected result straight out and ignored `--probe` in silence — the
//! failure `.print` was fixed for (#72), reintroduced one level down. Nothing
//! but this file would notice it coming back.

use std::process::Command;

use fairchild_core::nutmeg;

/// A transient with a node voltage, a source current, and a rising edge, so a
/// narrowed rawfile is visibly narrower and a decoded one is visibly not zero.
const DECK: &str = "\
* output formats
V1 in 0 PULSE(0 1 0 1n 1n 1u 2u)
R1 in out 1k
C1 out 0 1p
.tran 100p 5n
";

fn deck_path(tag: &str, body: &str) -> std::path::PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!("fc_outfmt_{tag}_{}.sp", std::process::id()));
    std::fs::write(&p, body).expect("write deck");
    p
}

fn out_path(tag: &str) -> std::path::PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!("fc_outfmt_{tag}_{}.raw", std::process::id()));
    p
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_fairchild"))
        .args(args)
        .output()
        .expect("run fairchild")
}

/// Run the deck to a file and return the bytes written.
fn to_file(tag: &str, body: &str, extra: &[&str]) -> Vec<u8> {
    let deck = deck_path(tag, body);
    let out = out_path(tag);
    let mut argv = vec![
        "-f",
        deck.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "-q",
    ];
    argv.extend_from_slice(extra);
    let r = run(&argv);
    assert!(
        r.status.success(),
        "{tag}: {:?}\n{}",
        r.status,
        String::from_utf8_lossy(&r.stderr)
    );
    let bytes = std::fs::read(&out).expect("output file");
    let _ = std::fs::remove_file(&out);
    let _ = std::fs::remove_file(&deck);
    bytes
}

/// Run the deck to standard output and return the bytes written.
fn to_stdout(tag: &str, body: &str, extra: &[&str]) -> Vec<u8> {
    let deck = deck_path(tag, body);
    let mut argv = vec!["-f", deck.to_str().unwrap(), "-q"];
    argv.extend_from_slice(extra);
    let r = run(&argv);
    assert!(
        r.status.success(),
        "{tag}: {:?}\n{}",
        r.status,
        String::from_utf8_lossy(&r.stderr)
    );
    let _ = std::fs::remove_file(&deck);
    r.stdout
}

/// A binary rawfile from the real binary must decode, and carry the run.
///
/// `--format binary` has no other end-to-end check: the core tests exercise the
/// writer, not the flag that selects it.
#[test]
fn format_binary_writes_a_decodable_rawfile() {
    let bytes = to_file("bin", DECK, &["--format", "binary"]);
    let f = nutmeg::read(&bytes).expect("binary rawfile must decode");
    assert_eq!(f.plotname, "Transient Analysis");
    assert!(f.points.len() > 10, "only {} points", f.points.len());
    let v = f.series("v(out)").expect("v(out) must be present");
    assert!(
        v.iter().any(|x| *x > 0.5),
        "v(out) never rises; the values did not survive the round trip"
    );
    // Binary is smaller than the text it replaces, which is half the reason to
    // have it. Asserting it here keeps the flag from silently writing ASCII.
    let ascii = to_file("bin_ascii", DECK, &["--format", "nutmeg"]);
    assert!(
        bytes.len() < ascii.len(),
        "binary is {} bytes against ASCII's {}",
        bytes.len(),
        ascii.len()
    );
}

/// The seekable and non-seekable paths must agree.
///
/// A rawfile written to a file streams and patches its point count in place; to
/// standard output it is assembled in memory, because there is nowhere to seek
/// back to. Two mechanisms, one file — so the count that ends up in the header
/// has to be right in both, and it is the *file* one that no in-process test
/// reaches.
#[test]
fn a_file_and_standard_output_produce_the_same_rawfile() {
    for fmt in ["nutmeg", "binary"] {
        let a = to_file("seek", DECK, &["--format", fmt]);
        let b = to_stdout("nonseek", DECK, &["--format", fmt]);
        assert_eq!(
            a,
            b,
            "{fmt}: a file and standard output differ\n--- file ---\n{}\n--- stdout ---\n{}",
            String::from_utf8_lossy(&a[..200.min(a.len())]),
            String::from_utf8_lossy(&b[..200.min(b.len())]),
        );
        let f = nutmeg::read(&a).unwrap_or_else(|e| panic!("{fmt}: {e}"));
        let stated: usize = String::from_utf8_lossy(&a[..400.min(a.len())])
            .lines()
            .find_map(|l| l.strip_prefix("No. Points:"))
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or_else(|| panic!("{fmt}: no parsable point count"));
        assert_eq!(
            stated,
            f.points.len(),
            "{fmt}: header says {stated}, body carries {}",
            f.points.len()
        );
    }
}

/// `--probe` narrows a `.tran` rawfile — and still does when the deck has a
/// `.measure`.
///
/// The two go through the same writer by different routes: without a
/// `.measure` the run streams, with one it is collected for the measurement and
/// replayed. The collected route used to write the result straight out, so
/// `--probe` was accepted and ignored. Both routes are checked because only one
/// of them was ever wrong.
#[test]
fn probe_narrows_a_transient_rawfile_on_both_routes() {
    let with_measure = format!("{DECK}.measure tran vmax MAX V(out)\n");
    for (tag, body) in [("plain", DECK.to_string()), ("measured", with_measure)] {
        let all = to_file(&format!("{tag}_all"), &body, &["--format", "binary"]);
        let one = to_file(
            &format!("{tag}_one"),
            &body,
            &["--format", "binary", "--probe", "V(out)"],
        );

        let fa = nutmeg::read(&all).unwrap();
        let fo = nutmeg::read(&one).unwrap();
        let names: Vec<&str> = fo.vars.iter().map(|v| v.name.as_str()).collect();
        assert_eq!(
            names,
            vec!["time", "v(out)"],
            "{tag}: --probe did not narrow the rawfile; it carries {names:?}"
        );
        assert!(
            one.len() < all.len(),
            "{tag}: narrowed file is {} bytes against {} for everything",
            one.len(),
            all.len()
        );
        assert_eq!(
            fo.points.len(),
            fa.points.len(),
            "{tag}: narrowing dropped timepoints as well as columns"
        );
        assert_eq!(
            fo.series("v(out)").unwrap(),
            fa.series("v(out)").unwrap(),
            "{tag}: the kept column changed value"
        );
    }
}

/// A `.measure` still reports, through the replay path.
///
/// The measurement is taken from the collected result and the result is then
/// replayed to the writer; a refactor that streamed unconditionally would lose
/// the measurement rather than fail, so it is asserted rather than assumed.
#[test]
fn a_measure_still_reports_when_the_output_is_a_rawfile() {
    let body = format!("{DECK}.measure tran vmax MAX V(out)\n");
    let deck = deck_path("meas", &body);
    let out = out_path("meas");
    let r = run(&[
        "-f",
        deck.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--format",
        "binary",
    ]);
    assert!(r.status.success(), "{:?}", r.status);
    let err = String::from_utf8_lossy(&r.stderr);
    assert!(err.contains("vmax"), "no measurement reported:\n{err}");
    let _ = std::fs::remove_file(&out);
    let _ = std::fs::remove_file(&deck);
}

/// An unmatched probe on a transient is refused, not dropped.
///
/// Same contract as #72, on a different code path: a transient selects columns
/// in the solver's layout rather than in rendered text, so it reports the miss
/// through the run's own error instead of through the CSV filter. The user has
/// to be told either way.
#[test]
fn an_unmatched_probe_fails_a_transient() {
    let deck = deck_path("miss", DECK);
    for fmt in ["csv", "nutmeg", "binary"] {
        let r = run(&[
            "-f",
            deck.to_str().unwrap(),
            "--format",
            fmt,
            "--probe",
            "V(out),V(total_nonsense)",
        ]);
        assert!(
            !r.status.success(),
            "{fmt}: a probe that matches nothing must fail the run"
        );
        let err = String::from_utf8_lossy(&r.stderr);
        assert!(
            err.contains("V(total_nonsense)"),
            "{fmt}: the unmatched name must be spelled out:\n{err}"
        );
    }
    let _ = std::fs::remove_file(&deck);
}

/// A rawfile from a bounded analysis carries every signal, and says so.
///
/// `--probe` selects before writing only where the output is unbounded in time.
/// The seam is real; a silent seam is what this codebase treats as the worst
/// outcome, so it is a warning. This pins the warning, not the seam — close the
/// seam and this test is the one that should be deleted.
#[test]
fn a_bounded_analysis_warns_that_probe_does_not_narrow_its_rawfile() {
    let deck = deck_path(
        "bounded",
        "* op\nV1 in 0 1\nR1 in mid 1k\nR2 mid 0 1k\n.op\n",
    );
    let r = run(&[
        "-f",
        deck.to_str().unwrap(),
        "--format",
        "binary",
        "--probe",
        "V(mid)",
    ]);
    assert!(r.status.success(), "{:?}", r.status);
    let err = String::from_utf8_lossy(&r.stderr);
    assert!(
        err.contains("--probe does not narrow a rawfile") && err.contains(".op"),
        "no warning that the rawfile is not narrowed:\n{err}"
    );

    // And CSV, which does filter, must not carry the warning.
    let r = run(&[
        "-f",
        deck.to_str().unwrap(),
        "--format",
        "csv",
        "--probe",
        "V(mid)",
    ]);
    assert!(r.status.success());
    let err = String::from_utf8_lossy(&r.stderr);
    assert!(
        !err.contains("does not narrow"),
        "CSV filters, so it must not warn:\n{err}"
    );
    let _ = std::fs::remove_file(&deck);
}
