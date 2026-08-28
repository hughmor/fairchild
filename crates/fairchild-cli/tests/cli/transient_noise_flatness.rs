//! The transient-noise flatness probe has to fire on a native device.
//!
//! `TransientNoise::draw` injects one held i.i.d. sample per generator per step,
//! which is white by construction — so a generator whose density is *not* flat
//! across the band the step resolves cannot be realised, only approximated at
//! mid-band. `warn_if_not_flat` exists to say so.
//!
//! Until `KF`/`AF` were modelled, the only way to reach that warning was an OSDI
//! model calling `flicker_noise()`, so nothing in this tree exercised it and its
//! threshold had never been checked against a real density. A native diode with
//! `KF` now produces exactly the shape it is looking for.
//!
//! Driving the real binary is the only way: the warning leaves library code
//! through `warn_user!` to stderr, which no in-process test can capture — the same
//! reason `quiet.rs` and `dropped_parameters.rs` do it this way.
//!
//! # What the threshold means, and why the arithmetic is written out
//!
//! A diode's generator is shot noise **plus** flicker in one density, and the
//! probe compares the total at `1/(20h)` against `1/(2.2h)` — a 9.1× frequency
//! span — warning above 1.05×. So it fires when flicker is a >5% contributor
//! across the resolved band, not merely when `KF` is set. That is the right
//! condition: when shot dominates, the total genuinely is flat and the ZOH
//! injector is exact.
//!
//! At `h = 1 ns` the band is 5.0e7 … 4.55e8 Hz. With `IS=1e-14` at 0.6 V through
//! 100 Ω the diode carries about 8.5e-5 A, so:
//!
//! ```text
//! shot          = 2q·Id                      = 2.72e-23  A²/Hz  (flat)
//! flicker(f)    = KF·Id/f
//! KF = 1e-12:   1.70e-24 at f_lo, 1.87e-25 at f_hi
//! total ratio   = 2.894e-23 / 2.742e-23      = 1.055     -> fires
//! ```
//!
//! and `KF = 1e-14` puts flicker two orders lower, so the ratio is 1.0006 and it
//! correctly stays silent.

use std::path::PathBuf;
use std::process::Command;

/// A forward-biased diode with a resistor, run with transient noise on.
///
/// `IS`/bias/`R1` are fixed because the arithmetic in the module docs depends on
/// the resulting drain current — change them and the thresholds below move.
fn deck(kf: f64) -> String {
    format!(
        "* transient noise flatness\n\
         .options trannoise=1\n\
         .model dm D (IS=1e-14 N=1 KF={kf:e} AF=1)\n\
         V1 a 0 DC 0.6\nR1 a b 100\nD1 b 0 dm\n\
         .tran 1n 200n\n"
    )
}

/// stderr from a run of the real binary.
fn stderr_for(kf: f64) -> String {
    let mut path = std::env::temp_dir();
    path.push(format!("fc_flatness_{}_{kf:e}.sp", std::process::id()));
    std::fs::write(&path, deck(kf)).expect("write deck");
    let path: PathBuf = path;
    let out = Command::new(env!("CARGO_BIN_EXE_fairchild"))
        .args(["-f", path.to_str().unwrap(), "-o", "/dev/null"])
        .output()
        .expect("run fairchild");
    assert!(out.status.success(), "fairchild exited {:?}", out.status);
    String::from_utf8_lossy(&out.stderr).into_owned()
}

fn fired(err: &str) -> bool {
    err.contains("varies by")
}

/// The probe fires when flicker is a material part of the band, and says what to
/// do about it.
#[test]
fn the_probe_fires_when_flicker_is_material() {
    let err = stderr_for(1e-12);
    assert!(
        fired(&err),
        "a native diode with KF=1e-12 has a generator varying 1.055× across the \
         band a 1 ns step resolves, and transient noise can only realise one \
         density. Nothing warned:\n{err}"
    );
    let line = err
        .lines()
        .find(|l| l.contains("varies by"))
        .expect("the line just matched");
    // A warning that does not say what to do about it is not much better than
    // silence — the same standard the unmodelled-parameter diagnostics are held
    // to.
    for phrase in ["mid-band", "shorten the step", ".noise"] {
        assert!(
            line.contains(phrase),
            "the warning must say what the user can do — missing '{phrase}': {line}"
        );
    }
}

/// And it stays silent when the density really is flat.
///
/// This is the half that stops the probe becoming noise people learn to skip. Two
/// separate reasons to be silent, and both have to hold:
///
/// * `KF = 0` — no flicker term at all;
/// * `KF` small enough that shot noise dominates the resolved band, where the
///   total *is* flat and the injector is exact.
#[test]
fn the_probe_stays_silent_when_the_density_is_flat() {
    for kf in [0.0, 1e-14, 1e-16] {
        let err = stderr_for(kf);
        assert!(
            !fired(&err),
            "KF={kf:e} leaves shot noise dominating the 5.0e7…4.55e8 Hz band, so \
             the total density is flat to better than 1.05× and the probe must not \
             fire. A warning here trains users to ignore the class:\n{err}"
        );
    }
}

/// `--quiet` silences it like every other warning.
///
/// Checked here rather than in `quiet.rs` because that file's deck cannot reach
/// this warning: it needs `trannoise=1` and a device carrying `KF`.
#[test]
fn quiet_silences_the_probe() {
    let mut path = std::env::temp_dir();
    path.push(format!("fc_flatness_quiet_{}.sp", std::process::id()));
    std::fs::write(&path, deck(1e-12)).expect("write deck");
    let out = Command::new(env!("CARGO_BIN_EXE_fairchild"))
        .args(["-f", path.to_str().unwrap(), "-o", "/dev/null", "--quiet"])
        .output()
        .expect("run fairchild");
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.is_empty(),
        "--quiet must silence the flatness probe too:\n{err}"
    );
}
