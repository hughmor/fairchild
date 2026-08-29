//! The dialect is the caller's to state, and the transliteration is inspectable.
//!
//! Spectre decides which language a stream is in by how it is entered, not by
//! looking at the text — that is why `simulator lang=` is a statement rather than
//! a file header. Detection is a convenience for the common case, not the only
//! way in.
//!
//! Driven through the binary because both features are CLI surface, and because
//! `--emit-spice` writes to stdout and exits, which no in-process test sees.

use std::io::Write;
use std::process::Command;

/// A Spectre fragment with neither a `simulator lang=` line nor a `//` comment,
/// which is what a foundry library's included files look like. Detection reads
/// this as SPICE.
const BARE_SPECTRE: &str = "\
model dm diode is=1e-14 n=1
v1 (a 0) vsource dc=0.6
d1 (a 0) dm
";

/// The same circuit, announcing itself.
const MARKED_SPECTRE: &str = "\
// a diode
simulator lang=spectre
model dm diode is=1e-14 n=1
v1 (a 0) vsource dc=0.6
d1 (a 0) dm
";

fn write(name: &str, text: &str) -> std::path::PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!("fc_dialect_{}_{name}", std::process::id()));
    let mut f = std::fs::File::create(&p).expect("create");
    f.write_all(text.as_bytes()).expect("write");
    p
}

fn run(path: &std::path::Path, args: &[&str]) -> (bool, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_fairchild"))
        .arg("-f")
        .arg(path)
        .args(args)
        .output()
        .expect("run fairchild");
    let text = format!(
        "{}{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    (out.status.success(), text)
}

/// `--lang spectre` reads a file the detection would call SPICE.
///
/// The detection defaults to SPICE, which is right — most decks are — but it
/// leaves a Spectre fragment carrying no marker unreadable. That is exactly the
/// shape of an included foundry file.
#[test]
fn lang_spectre_reads_a_file_detection_would_not() {
    let p = write("bare.scs", BARE_SPECTRE);

    let (ok_auto, auto) = run(&p, &[]);
    assert!(
        !ok_auto,
        "detection must still read this as SPICE, or the test proves nothing: {auto}"
    );

    let (_, forced) = run(&p, &["--lang", "spectre"]);
    assert!(
        !forced.contains("Mname drain gate source bulk"),
        "`--lang spectre` must transliterate rather than read Spectre as SPICE: {forced}"
    );
}

/// `--lang spice` refuses to transliterate a file that announces itself as
/// Spectre, so the flag is an override in both directions.
#[test]
fn lang_spice_overrides_a_spectre_marker() {
    let p = write("marked.scs", MARKED_SPECTRE);

    let (_, auto) = run(&p, &["--emit-spice"]);
    assert!(
        auto.contains("d1 a 0 dm"),
        "detection reads this as Spectre and transliterates it: {auto}"
    );

    let (_, forced) = run(&p, &["--emit-spice", "--lang", "spice"]);
    assert!(
        forced.contains("v1 (a 0) vsource"),
        "`--lang spice` must pass the text through untransliterated: {forced}"
    );
}

/// A file read as SPICE that looks like Spectre says so when it fails.
///
/// Without this the error is about a malformed `M` line, which mentions neither
/// the dialect nor the fix — a long way from the cause when the real problem is
/// that a `.scs` fragment arrived without its marker.
#[test]
fn a_misdetected_spectre_file_says_what_to_try() {
    let p = write("hint.scs", BARE_SPECTRE);
    let (ok, text) = run(&p, &[]);
    assert!(!ok, "this must fail, or there is nothing to hint about");
    for needle in ["Spectre", "--lang spectre"] {
        assert!(
            text.contains(needle),
            "the failure must name {needle}: {text}"
        );
    }

    // And the hint must not fire when the caller already said. A hint that
    // appears on a deck the user has explicitly typed a flag for is noise.
    let (_, forced) = run(&p, &["--lang", "spice"]);
    assert!(
        !forced.contains("--lang spectre"),
        "an explicit --lang answers the question, so the hint is noise: {forced}"
    );
}

/// `--emit-spice` prints what a Spectre deck became, with includes resolved.
///
/// The transliteration is an intermediate the user never sees, and when a Spectre
/// deck fails the first question is what it turned into. Answering that used to
/// mean writing a program against the parser.
#[test]
fn emit_spice_prints_the_transliteration() {
    let p = write("emit.scs", MARKED_SPECTRE);
    let (ok, text) = run(&p, &["--emit-spice"]);
    assert!(ok, "the dump must succeed: {text}");
    for needle in [".model dm diode", "v1 a 0 DC 0.6", "d1 a 0 dm"] {
        assert!(
            text.contains(needle),
            "the SPICE form must contain `{needle}`: {text}"
        );
    }
    assert!(
        !text.contains("vsource"),
        "the Spectre spelling must be gone: {text}"
    );
}

/// A SPICE deck passes through the dump unchanged, so the flag is usable on
/// either dialect and says so by not lying.
#[test]
fn emit_spice_on_a_spice_deck_changes_nothing() {
    let deck = "* an RC\nV1 in 0 DC 1\nR1 in out 1k\nC1 out 0 1u\n.op\n";
    let p = write("plain.sp", deck);
    let (ok, text) = run(&p, &["--emit-spice"]);
    assert!(ok, "{text}");
    for line in ["V1 in 0 DC 1", "R1 in out 1k", "C1 out 0 1u"] {
        assert!(text.contains(line), "`{line}` must survive: {text}");
    }
}

/// An unknown dialect is refused by name rather than silently defaulting.
#[test]
fn an_unknown_lang_is_refused() {
    let p = write("x.sp", "* rc\nV1 a 0 DC 1\nR1 a 0 1k\n.op\n");
    let (ok, text) = run(&p, &["--lang", "verilog"]);
    assert!(!ok, "an unknown dialect must not be accepted: {text}");
    assert!(
        text.contains("verilog") || text.contains("auto, spectre or spice"),
        "the refusal must say what was wrong: {text}"
    );
}
