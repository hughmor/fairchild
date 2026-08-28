//! Flicker (1/f) noise against ngspice, for all three native model families.
//!
//! `KF` and `AF` were accepted and not modelled on the diode, the BJT *and* the
//! MOSFET (#77 §7), so neither `.noise` nor transient noise could produce 1/f
//! from a native device: those analyses inject exactly the generators the models
//! report, and none of them reported one.
//!
//! # Why the slope *and* the magnitude
//!
//! A 1/f source has two things to get wrong independently. A test that checks
//! only the magnitude at one frequency passes with any spectral shape, and a test
//! that checks only the slope passes with the wrong coefficient. Each test below
//! does both:
//!
//! * **slope**, from a decade of frequency inside the flicker-dominated band,
//!   which is a property of the density's `1/f` and of nothing else;
//! * **magnitude**, against ngspice at a fixed frequency, which is where `KF`,
//!   `AF`, the driving current and any normalising factor all land.
//!
//! `AF` is swept, because `AF = 1` makes the density linear in the current and
//! hides an exponent applied to the wrong quantity.
//!
//! Requires ngspice on PATH; skipped, not failed, without it.

use std::io::Write;
use std::process::Command;

use fairchild_core::{freq_decade, noise_analysis, options::SimOptions, DeviceRegistry};
use fairchild_parser::parse_spice;

/// fairchild's output-referred noise spectrum over `freqs`.
fn fairchild_psd(body: &str, out: &str, src: &str, freqs: &[f64]) -> Vec<f64> {
    let deck = format!("* flicker\n{body}");
    let net = parse_spice(&deck).expect("parse");
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);
    let opts = SimOptions::from_netlist(&net);
    let r = noise_analysis(&net, freqs, out, "0", src, &reg, &opts)
        .unwrap_or_else(|e| panic!("noise failed on\n{deck}\n{e:?}"));
    r.onoise_psd
}

/// ngspice's output-referred noise spectrum at the same frequencies, as a
/// **power** density.
///
/// ngspice's `onoise_spectrum` is an amplitude density in V/√Hz and fairchild
/// reports V²/Hz, so every value is squared here. Worth stating rather than
/// burying: the first run of this file read the two as the same quantity and
/// reported the diode as seven orders wrong, when squaring made it agree to
/// 3e-6.
///
/// `print` rather than `echo $&vec[i]`: the spectrum is read as a whole so the
/// indices cannot slip against fairchild's frequency list.
fn ngspice_psd(
    body: &str,
    out: &str,
    src: &str,
    decade_pts: usize,
    f0: f64,
    f1: f64,
) -> Option<Vec<f64>> {
    let mut tmp = tempfile::NamedTempFile::new().ok()?;
    write!(
        tmp,
        "* flicker\n{body}.control\n\
         noise V({out}) {src} DEC {decade_pts} {f0} {f1}\n\
         setplot noise1\nprint onoise_spectrum\n.endc\n.end\n"
    )
    .ok()?;
    let out_run = Command::new("ngspice")
        .arg("-b")
        .arg(tmp.path())
        .output()
        .ok()?;
    let text = String::from_utf8_lossy(&out_run.stdout);
    let mut vals = Vec::new();
    for line in text.lines() {
        let cols: Vec<&str> = line.split_whitespace().collect();
        // `print` emits `<index> <freq> <value>`.
        if cols.len() == 3 && cols[0].parse::<usize>().is_ok() {
            if let Ok(v) = cols[2].parse::<f64>() {
                vals.push(v * v);
            }
        }
    }
    (!vals.is_empty()).then_some(vals)
}

/// The exponent of a power law fitted through the first and last points.
fn slope(freqs: &[f64], psd: &[f64]) -> f64 {
    let (f0, f1) = (freqs[0], freqs[freqs.len() - 1]);
    let (p0, p1) = (psd[0], psd[psd.len() - 1]);
    (p1 / p0).ln() / (f1 / f0).ln()
}

fn check(what: &str, body: &str, out: &str, src: &str) {
    // A decade well below any corner, so flicker dominates and the slope is the
    // density's own rather than a mixture with the white floor.
    let freqs = freq_decade(1.0, 10.0, 5);
    let fc = fairchild_psd(body, out, src, &freqs);

    let s = slope(&freqs, &fc);
    assert!(
        (s - (-1.0)).abs() < 0.05,
        "{what}: the spectrum's exponent is {s:.3}, and a 1/f density is -1. \
         An exponent near 0 means the flicker term is missing and only the white \
         floor is there."
    );

    let Some(ng) = ngspice_psd(body, out, src, 5, 1.0, 10.0) else {
        eprintln!("ngspice not available — slope checked, magnitude skipped for '{what}'");
        return;
    };
    assert_eq!(
        ng.len(),
        fc.len(),
        "{what}: ngspice returned {} points and fairchild {}",
        ng.len(),
        fc.len()
    );
    for (i, (a, b)) in fc.iter().zip(&ng).enumerate() {
        let rel = (a - b).abs() / b.abs().max(1e-300);
        assert!(
            rel < 5e-3,
            "{what} at f={:.4} Hz: fairchild {a:.6e}, ngspice {b:.6e} (rel {rel:.2e})",
            freqs[i]
        );
    }
}

// ---------------------------------------------------------------------------

/// A forward-biased diode's 1/f noise: `KF·|Id|^AF / f`.
#[test]
fn diode_flicker_matches_ngspice() {
    for af in [1.0, 1.2] {
        let body = format!(
            ".model dm D (IS=1e-14 N=1 KF=1e-14 AF={af})\n\
             V1 a 0 DC 0.6 AC 1\nR1 a b 100\nD1 b 0 dm\n"
        );
        check(&format!("diode flicker AF={af}"), &body, "b", "v1");
    }
}

/// A BJT's 1/f noise is driven by the **base** current, not the collector's.
///
/// `AF` is swept because at `AF = 1` a density built from the wrong current
/// differs only by a constant, which the magnitude check would catch but which is
/// easy to mistake for a coefficient error. At `AF ≠ 1` the two disagree in shape
/// as well.
///
/// # The base resistor is not decoration
///
/// The first version drove the base with `V1 b 0 AC 1` directly. The flicker
/// generator sits base-to-emitter, so an ideal source across those terminals
/// carries its whole current and the output sees only white noise — the slope
/// came back as exactly 0.000. `RB` gives the generator somewhere to develop a
/// voltage. A noise source shorted by the topology contributes nothing, and that
/// is a property of the deck rather than of the model.
#[test]
fn bjt_flicker_matches_ngspice() {
    for af in [1.0, 1.2] {
        let body = format!(
            ".model qm NPN (IS=1e-16 BF=100 KF=1e-14 AF={af})\n\
             V1 in 0 DC 0.75 AC 1\nRB in b 10k\n\
             VC c 0 DC 2\nRC c cc 1k\nQ1 cc b 0 qm\n"
        );
        check(&format!("BJT flicker AF={af}"), &body, "cc", "v1");
    }
}

/// A MOSFET's 1/f noise, against the closed form — **not** against ngspice.
///
/// # Why ngspice is not the anchor here
///
/// Its MOS1 flicker density is bit-identical at `AF` = 0.5, 1.0, 1.2 and 2.0:
/// 3.706770e-11 V²/Hz at every one of them. Over the same sweep the *diode's*
/// `AF` moves as a clean power law (ratios 1, 0.0092, 0.0014), so this is a
/// property of ngspice's MOS1 and not of the deck or the card syntax. It also
/// scales as `W¹·L⁻³` where the documented SPICE3 form gives `W⁰·L⁻²` — `W⁰`
/// because `Id ∝ W/L` cancels the `W` in the denominator.
///
/// Asserting ngspice's number would mean asserting a density that ignores a
/// parameter the card sets. So the anchor is the closed form:
/// `KF·|Id|^AF / (f·W·Leff·Cox)` with `Id` from Level 1's own saturation
/// expression and `Cox = ε_ox/TOX`, and the structural dependencies are checked
/// separately in `the_mosfet_normalisation_is_read`.
#[test]
fn mosfet_flicker_matches_the_closed_form() {
    const EPS_OX: f64 = 3.9 * 8.854187817e-12;
    let (kp, w, l, tox, vgs, vto, rd): (f64, f64, f64, f64, f64, f64, f64) =
        (100e-6, 10e-6, 1e-6, 20e-9, 1.5, 0.7, 1000.0);
    // Level 1 saturation, with no channel-length modulation (LAMBDA unset).
    let ids = 0.5 * kp * (w / l) * (vgs - vto) * (vgs - vto);
    let cox = EPS_OX / tox;

    for af in [1.0, 1.2] {
        let body = format!(
            ".model nm NMOS (VTO={vto} KP={kp:e} TOX={tox:e} KF=1e-24 AF={af})\n\
             V1 g 0 DC {vgs} AC 1\nVD d 0 DC 3\nRD d dd {rd:e}\n\
             M1 dd g 0 0 nm W={w:e} L={l:e}\n"
        );
        let freqs = freq_decade(1.0, 10.0, 5);
        let psd = fairchild_psd(&body, "dd", "v1", &freqs);

        // The slope is the density's own, independent of every coefficient.
        let sl = slope(&freqs, &psd);
        assert!(
            (sl - (-1.0)).abs() < 0.05,
            "MOSFET AF={af}: exponent {sl:.3}, and a 1/f density is -1"
        );

        // The magnitude at 1 Hz. In saturation the drain node's impedance is set
        // by RD against the device's output conductance, which is 1/LAMBDA=inf
        // here, so the transfer is RD² and the closed form is exact.
        let want = rd * rd * 1e-24 * ids.powf(af) / (1.0 * w * l * cox);
        let rel = (psd[0] - want).abs() / want;
        assert!(
            rel < 2e-2,
            "MOSFET AF={af} at 1 Hz: {:.6e} against the closed form {want:.6e} \
             (rel {rel:.2e}). Id={ids:.6e} A, Cox={cox:.6e} F/m².",
            psd[0]
        );
    }
}

/// `W`, `L` and `TOX` all appear in the MOSFET's denominator, so the density has
/// to move when they do. Without the normalisation the ratio below is 1.
#[test]
fn the_mosfet_normalisation_is_read() {
    let psd = |w: f64, l: f64, tox: f64| {
        let body = format!(
            ".model nm NMOS (VTO=0.7 KP=100u TOX={tox:e} KF=1e-24 AF=1)\n\
             V1 g 0 DC 1.5 AC 1\nVD d 0 DC 3\nRD d dd 1k\n\
             M1 dd g 0 0 nm W={w:e} L={l:e}\n"
        );
        fairchild_psd(&body, "dd", "v1", &[1.0])[0]
    };
    let base = psd(10e-6, 1e-6, 20e-9);
    // Doubling W doubles the drain current too (Id goes as W/L), so the density
    // goes as |Id|^AF / W = W/W = flat in W at AF=1. L is the clean knob: it
    // divides the current *and* the normalisation, so the density falls as 1/L².
    let long = psd(10e-6, 2e-6, 20e-9);
    let ratio = base / long;
    assert!(
        (ratio - 4.0).abs() / 4.0 < 0.02,
        "doubling L at AF=1 divides the density by L in the normalisation and by \
         L in the current, so the ratio is 4: got {ratio:.4}. A ratio of 1 means \
         the normalisation is not read at all; 2 means only one of the two is."
    );
    // Thicker oxide is less capacitance, so a larger density.
    let thick = psd(10e-6, 1e-6, 40e-9);
    assert!(
        thick > base * 1.9,
        "doubling TOX halves COX and so doubles the density: {thick:.4e} against \
         {base:.4e}"
    );
}

/// `KF` unset means **no** flicker source, not a small one. Every existing noise
/// golden is this case.
#[test]
fn no_kf_means_no_flicker() {
    let freqs = freq_decade(1.0, 10.0, 5);
    let body = ".model dm D (IS=1e-14 N=1)\n\
                V1 a 0 DC 0.6 AC 1\nR1 a b 100\nD1 b 0 dm\n";
    let psd = fairchild_psd(body, "b", "v1", &freqs);
    let s = slope(&freqs, &psd);
    assert!(
        s.abs() < 1e-6,
        "with no KF the spectrum must be flat across a decade, exponent {s:.2e}"
    );
}

/// `KF` with no oxide capacitance is a division by zero, and is refused by name
/// rather than returning a non-finite density into the noise matrix.
#[test]
fn mosfet_kf_without_an_oxide_is_refused() {
    let deck = "* no tox\n.model nm NMOS (VTO=0.7 KP=100u KF=1e-24 AF=1)\n\
                V1 g 0 DC 1.5 AC 1\nVD d 0 DC 3\nRD d dd 1k\n\
                M1 dd g 0 0 nm W=10u L=1u\n.op\n";
    let net = parse_spice(deck).expect("parse");
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);
    let msg = match fairchild_core::dc_op_nr_with_registry(&net, &reg) {
        Ok(_) => panic!("KF with no TOX/COX must be refused, not solved"),
        Err(e) => e.to_string(),
    };
    for needle in ["KF", "TOX", "COX"] {
        assert!(
            msg.contains(needle),
            "the refusal must name {needle} so the fix is in the message: {msg}"
        );
    }
}
