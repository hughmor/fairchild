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
.model dm D (IS=1e-14 N=1 BV=50)
.model qm NPN (IS=1e-16 BF=100 VAF=50 IKF=1e-3 KF=1e-15)
.model nm NMOS (VTO=0.7 KP=100u CJ=0.5m CJSW=1n MJ=0.5 MJSW=0.33 RD=20)
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
        ("BV ignored", "breakdown"),
        ("KF ignored", "flicker"),
        ("RD ignored", "series resistance"),
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
    // nothing warns about a parameter that IS stamped.
    for honoured in ["IKF ignored", "MJSW ignored", "VAF ignored", "'area'"] {
        assert!(
            !err.contains(honoured),
            "warned about '{honoured}', which this simulator honours:\n{err}"
        );
    }

    // One line per card, not one per instance.
    assert_eq!(
        err.matches("BV ignored").count(),
        1,
        "a card's diagnostic must not repeat per instance:\n{err}"
    );
}

#[test]
fn quiet_silences_them_like_every_other_warning() {
    assert!(run(&["--quiet"]).is_empty());
}
