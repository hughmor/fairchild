//! A parameter that changes nothing has to *say* it changed nothing.
//!
//! The classification is unit-tested in `fairchild-core`; this drives the real
//! binary, because the diagnostics leave library code through `warn_user!` and
//! no in-process test can see them (the same reason `quiet.rs` exists). It also
//! pins the negative direction, which matters more here than usual: a warning
//! for a parameter that *is* honoured teaches users to ignore the whole class.

use std::path::PathBuf;
use std::process::Command;

/// One card and one instance of each device, mixing honoured parameters with
/// dropped ones.
const DECK: &str = "\
* dropped-parameter diagnostics
.model dm D (IS=1e-14 N=1 BV=50 IBV=1e-3 ISR=1e-12)
.model qm NPN (IS=1e-16 BF=100 VAF=50 IKF=1e-3 XCJC=0.5 PTF=30)
.model nm NMOS (VTO=0.7 KP=100u CJ=0.5m CJSW=1n MJ=0.5 MJSW=0.33 RD=20 RSH=50)
V1 a 0 DC 1.5
R1 a d 1k
D1 d 0 dm area=2 banana=3
Q1 a a 0 qm area=4
M1 a a 0 0 nm W=2u L=1u fruitbat=7
.op
";

fn run(args: &[&str]) -> String {
    let mut deck = std::env::temp_dir();
    deck.push(format!("fairchild_dropped_{}.sp", std::process::id()));
    std::fs::write(&deck, DECK).expect("write deck");
    let path: PathBuf = deck;
    let mut argv = vec!["-f", path.to_str().unwrap()];
    argv.extend_from_slice(args);
    let out = Command::new(env!("CARGO_BIN_EXE_fairchild"))
        .args(&argv)
        .output()
        .expect("run fairchild");
    assert!(out.status.success(), "fairchild exited {:?}", out.status);
    String::from_utf8_lossy(&out.stderr).into_owned()
}

#[test]
fn dropped_parameters_are_named_and_honoured_ones_are_not() {
    let err = run(&[]);

    // Model-card parameters that are accepted and do nothing: named once, with
    // what the deck loses rather than just the key.
    for (key, phrase) in [
        ("ISR ignored", "recombination"),
        ("PTF ignored", "excess phase"),
        ("RSH ignored", "NRD/NRS"),
    ] {
        assert!(err.contains(key), "no diagnostic for {key}:\n{err}");
        let line = err
            .lines()
            .find(|l| l.contains(key))
            .expect("the line we just found");
        assert!(
            line.contains(phrase),
            "'{key}' has to say what is missing, not just the key: {line}"
        );
    }

    // Instance parameters no model can honour: named per instance, with the
    // element that carried them.
    assert!(
        err.contains("d1") && err.contains("banana"),
        "the diode's stray instance parameter is unnamed:\n{err}"
    );
    assert!(
        err.contains("m1") && err.contains("fruitbat"),
        "the MOSFET's stray instance parameter is unnamed:\n{err}"
    );

    // And the negative direction, which is what stops this becoming noise:
    // nothing warns about a parameter that IS stamped. `BV`/`IBV` are on this
    // list rather than the one above because reverse breakdown is modelled now —
    // implementing a parameter is supposed to move it across.
    for honoured in [
        "BV ignored",
        "IBV ignored",
        "KF ignored",
        "AF ignored",
        "RD ignored",
        "RS ignored",
        "IKF ignored",
        "MJSW ignored",
        "VAF ignored",
        "XCJC ignored",
        "'area'",
    ] {
        assert!(
            !err.contains(honoured),
            "warned about '{honoured}', which this simulator honours:\n{err}"
        );
    }

    // One line per card, not one per instance.
    assert_eq!(
        err.matches("ISR ignored").count(),
        1,
        "a card's diagnostic must not repeat per instance:\n{err}"
    );
}

#[test]
fn quiet_silences_them_like_every_other_warning() {
    assert!(run(&["--quiet"]).is_empty());
}

// ── --probe: an unmatched signal name is refused, not dropped (issue #72) ──
//
// Same class as the parameter diagnostics above: the user asked for something
// by name, and the run must either honour it or say so. A silently-missing CSV
// column fails a long way downstream (a KeyError, a plot with one fewer trace)
// instead of at the typo that caused it.

const PROBE_DECK: &str = "\
* probe refusal
V1 in 0 DC 1
R1 in out 1k
R2 out 0 1k
.op
";

fn run_probe(deck: &str, tag: &str, args: &[&str]) -> std::process::Output {
    let mut path = std::env::temp_dir();
    path.push(format!("fairchild_probe_{tag}_{}.sp", std::process::id()));
    std::fs::write(&path, deck).expect("write deck");
    let mut argv = vec!["-f", path.to_str().unwrap()];
    argv.extend_from_slice(args);
    Command::new(env!("CARGO_BIN_EXE_fairchild"))
        .args(&argv)
        .output()
        .expect("run fairchild")
}

#[test]
fn an_unmatched_probe_is_a_named_error_not_a_missing_column() {
    let out = run_probe(
        PROBE_DECK,
        "unmatched",
        &["--probe", "V(out),V(total_nonsense)"],
    );
    assert!(
        !out.status.success(),
        "a probe that matches nothing must fail the run, got {:?}",
        out.status
    );
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains("V(total_nonsense)"),
        "the unmatched name must be spelled out:\n{err}"
    );
    assert!(
        !err.contains("V(out)'"),
        "the probe that DID match must not be blamed:\n{err}"
    );
}

#[test]
fn matching_probes_still_filter_and_succeed() {
    let out = run_probe(PROBE_DECK, "matched", &["--probe", "V(out)"]);
    assert!(out.status.success(), "{:?}", out.status);
    let csv = String::from_utf8_lossy(&out.stdout);
    let header = csv.lines().next().unwrap_or_default();
    assert_eq!(header, "analysis,V(out)", "filtered header:\n{csv}");
}

/// An AC sweep spells the columns `mag_V(x)` / `phase_deg_V(x)`; probing
/// `V(x)` selects that pair rather than silently matching nothing — which was
/// this same bug wearing its frequency-domain hat.
#[test]
fn an_ac_probe_selects_the_mag_phase_pair() {
    let deck = "\
* probe over ac
V1 in 0 DC 0 AC 1
R1 in out 1k
C1 out 0 1n
.ac dec 2 1 100
";
    let out = run_probe(deck, "ac", &["--probe", "V(out)"]);
    assert!(out.status.success(), "{:?}", out.status);
    let csv = String::from_utf8_lossy(&out.stdout);
    let header = csv.lines().next().unwrap_or_default();
    assert_eq!(
        header, "freq_hz,mag_V(out),phase_deg_V(out)",
        "filtered AC header:\n{csv}"
    );
}
