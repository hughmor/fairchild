//! `.tf` and `.pz` against ngspice.
//!
//! The closed-form tests in `tests/circuit/small_signal_analyses.rs` check that
//! these analyses agree with algebra.  This file checks the thing algebra
//! cannot: that we spell the answers the way the rest of the world does.
//! `.tf`'s two resistances and `.pz`'s pole convention (rad/s, and which sign
//! is stable) are conventions, not theorems — a wrong one is self-consistent
//! and every internal test passes anyway.  ngspice is the outside opinion.
//!
//! Skipped, not failed, when ngspice is absent.

use std::io::Write;
use std::process::Command;

use fairchild_core::pz::pole_zero;
use fairchild_core::tf::transfer_function;
use fairchild_core::{DeviceRegistry, SimOptions};
use fairchild_parser::{parse_spice, Analysis, Netlist, PzDrive, PzWant};

use super::ngspice_golden::find_ngspice;

/// Both simulators solve the same linear system; the gap is the adjoint's
/// finite difference against ngspice's direct solve, plus ngspice printing
/// seven significant figures.
const REL: f64 = 1e-5;

fn registry(net: &Netlist) -> DeviceRegistry {
    let mut r = DeviceRegistry::new();
    r.register_builtin_models(&net.models);
    r
}

fn tight() -> SimOptions {
    SimOptions {
        reltol: 1e-12,
        vntol: 1e-14,
        ..SimOptions::default()
    }
}

/// Run `deck` under ngspice in batch mode and hand back stdout.
///
/// `.tf` prints from a bare batch run, but `.pz` does not — batch mode wants a
/// `.print`/`.plot` it has no way to spell for a pole list, so the deck drives
/// it from a `.control` block instead.  That is ngspice's own quirk and not
/// something to model on our side.
fn ngspice_run(deck: &str) -> Option<String> {
    let bin = find_ngspice()?;
    let mut tmp = tempfile::NamedTempFile::new().ok()?;
    write!(tmp, "{deck}").ok()?;
    let out = Command::new(&bin).arg("-b").arg(tmp.path()).output().ok()?;
    Some(String::from_utf8_lossy(&out.stdout).into_owned())
}

/// `name = 1.234e5` off an ngspice `.tf` report.
fn scalar(out: &str, key: &str) -> Option<f64> {
    out.lines()
        .find_map(|l| l.split_once('=').filter(|(k, _)| k.trim() == key))
        .and_then(|(_, v)| v.trim().parse::<f64>().ok())
}

/// `pole(1) = -5.0e4,9.98e5` — ngspice prints complex as `re,im`.
fn roots(out: &str, kind: &str) -> Vec<(f64, f64)> {
    out.lines()
        .filter_map(|l| {
            let (k, v) = l.split_once('=')?;
            let k = k.trim();
            // A single root prints as `all = re,im` rather than `pole(1) = …`.
            if !(k.starts_with(kind) || k == "all") {
                return None;
            }
            let (re, im) = v.trim().split_once(',')?;
            Some((re.trim().parse().ok()?, im.trim().parse().ok()?))
        })
        .collect()
}

fn assert_close(got: f64, want: f64, what: &str) {
    let tol = REL * want.abs().max(1e-9);
    assert!(
        (got - want).abs() <= tol,
        "{what}: fairchild {got:e}, ngspice {want:e} (tol {tol:e})"
    );
}

/// Run one `.tf` deck both ways and compare all three numbers.
fn compare_tf(deck: &str, out_key: &str, in_key: &str) {
    let Some(ng) = ngspice_run(deck) else {
        eprintln!("ngspice not on PATH — skipping");
        return;
    };
    let (gain, r_out, r_in) = (
        scalar(&ng, "transfer_function").unwrap_or_else(|| panic!("{ng}")),
        scalar(&ng, out_key).unwrap_or_else(|| panic!("{ng}")),
        scalar(&ng, in_key).unwrap_or_else(|| panic!("{ng}")),
    );

    let net = parse_spice(deck).unwrap();
    let Some(Analysis::Tf { out, input_src }) = net
        .analyses
        .iter()
        .find(|a| matches!(a, Analysis::Tf { .. }))
        .cloned()
    else {
        panic!("deck declares no .tf")
    };
    let r = transfer_function(&net, &registry(&net), &tight(), &out, &input_src).unwrap();

    assert_close(r.gain, gain, "transfer_function");
    assert_close(r.r_out, r_out, out_key);
    assert_close(r.r_in, r_in, in_key);
}

#[test]
fn tf_divider_matches_ngspice() {
    compare_tf(
        "* divider\n\
         Vin in 0 DC 1\n\
         R1 in out 1k\n\
         R2 out 0 3k\n\
         .tf v(out) Vin\n\
         .end\n",
        "output_impedance_at_v(out)",
        "vin#input_impedance",
    );
}

/// A nonlinear circuit, so the comparison is of two *linearisations* about an
/// operating point and not of two ways to invert the same constant matrix.
/// This is where a wrong Jacobian would show and the divider would not.
#[test]
fn tf_diode_bias_matches_ngspice() {
    compare_tf(
        "* diode bias\n\
         Vin in 0 DC 2\n\
         R1 in out 1k\n\
         D1 out 0 dmod\n\
         .model dmod D(IS=1e-14 N=1)\n\
         .tf v(out) Vin\n\
         .end\n",
        "output_impedance_at_v(out)",
        "vin#input_impedance",
    );
}

/// A current-source input: `.tf`'s transfer is then a transresistance, and the
/// input resistance is a plain port resistance rather than a reciprocal.  The
/// two are computed by different branches here, so ngspice pins both.
#[test]
fn tf_current_input_matches_ngspice() {
    compare_tf(
        "* shunt\n\
         Iin 0 out DC 1m\n\
         R1 out 0 2k\n\
         R2 out 0 2k\n\
         .tf v(out) Iin\n\
         .end\n",
        "output_impedance_at_v(out)",
        "iin#input_impedance",
    );
}

/// Run one `.pz` deck both ways and compare the root sets.
fn compare_pz(body: &str, want: PzWant, kind: &str) {
    let keyword = match want {
        PzWant::Poles => "pol",
        PzWant::Zeros => "zer",
        PzWant::Both => "pz",
    };
    let deck = format!("{body}.control\npz in 0 out 0 vol {keyword}\nprint all\n.endc\n.end\n");
    let Some(ng) = ngspice_run(&deck) else {
        eprintln!("ngspice not on PATH — skipping");
        return;
    };
    let mut want_roots = roots(&ng, kind);
    assert!(!want_roots.is_empty(), "ngspice reported no {kind}s:\n{ng}");

    let net = parse_spice(body).unwrap();
    let res = pole_zero(
        &net,
        &registry(&net),
        &tight(),
        "in",
        "0",
        "out",
        "0",
        PzDrive::Vol,
        want,
    )
    .unwrap();
    let mut got: Vec<(f64, f64)> = match want {
        PzWant::Zeros => res.zeros.iter().map(|r| (r.re, r.im)).collect(),
        _ => res.poles.iter().map(|r| (r.re, r.im)).collect(),
    };

    assert_eq!(
        got.len(),
        want_roots.len(),
        "{kind} count: fairchild {got:?}, ngspice {want_roots:?}"
    );
    // Neither simulator promises an order, so compare as sets.
    let by_place =
        |a: &(f64, f64), b: &(f64, f64)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal);
    got.sort_by(by_place);
    want_roots.sort_by(by_place);
    for (g, w) in got.iter().zip(&want_roots) {
        assert_close(g.0, w.0, &format!("{kind} real part"));
        // ngspice prints a real root's imaginary part as an exact zero, so an
        // absolute floor is the right comparison there, not a relative one.
        if w.1.abs() < 1.0 {
            assert!(g.1.abs() < 1.0, "{kind} should be real: got {g:?}");
        } else {
            assert_close(g.1, w.1, &format!("{kind} imaginary part"));
        }
    }
}

#[test]
fn pz_rc_pole_matches_ngspice() {
    compare_pz(
        "* rc\n\
         Vin in 0 DC 0 AC 1\n\
         R1 in out 1k\n\
         C1 out 0 1n\n",
        PzWant::Poles,
        "pole",
    );
}

/// The inductor case — the one the branch-current construction exists for.  If
/// `.pz` had cleared the `1/s` into a quadratic pencil instead, this is where
/// the spurious roots would appear and the count would not match.
#[test]
fn pz_series_rlc_matches_ngspice() {
    compare_pz(
        "* rlc\n\
         Vin in 0 DC 0 AC 1\n\
         R1 in a 100\n\
         L1 a out 1m\n\
         C1 out 0 1n\n",
        PzWant::Poles,
        "pole",
    );
}

#[test]
fn pz_zero_matches_ngspice() {
    compare_pz(
        "* zero\n\
         Vin in 0 DC 0 AC 1\n\
         R1 in out 1k\n\
         R2 out mid 250\n\
         C1 mid 0 1n\n",
        PzWant::Zeros,
        "zero",
    );
}

/// Two poles and a zero at once, so the pairing cannot be right by accident on
/// a one-root circuit.
#[test]
fn pz_two_pole_network_matches_ngspice() {
    compare_pz(
        "* twopole\n\
         Vin in 0 DC 0 AC 1\n\
         R1 in mid 1k\n\
         C1 mid 0 1n\n\
         R2 mid out 4.7k\n\
         C2 out 0 2.2n\n",
        PzWant::Poles,
        "pole",
    );
}
