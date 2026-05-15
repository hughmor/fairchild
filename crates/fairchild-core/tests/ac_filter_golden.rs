//! AC frequency-response golden tests.
//!
//! Validates the small-signal AC machinery by matching the simulator's
//! magnitude / phase response against textbook closed-form formulas at
//! several frequencies — no ngspice dependency, no external models.

use fairchild_core::{ac_analysis, freq_decade, DeviceRegistry};
use fairchild_parser::parse_spice;

/// RC low-pass at three frequencies (DC, corner, decade above) — magnitude
/// and phase must track the analytic Bode plot.
#[test]
fn rc_lowpass_magnitude_and_phase() {
    let net = parse_spice(
        "* RC low-pass: R=1k, C=1µF → f_c = 159.155 Hz\n\
         V1 in 0 DC 1\n\
         R1 in out 1k\n\
         C1 out 0 1u\n\
         .end\n"
    ).unwrap();
    let registry = DeviceRegistry::new();
    let freqs = vec![1.0, 159.1549430918954, 1591.549430918954, 1e5];
    let result = ac_analysis(&net, &freqs, None, &registry).expect("AC analysis");

    // Analytic: |H(jω)| = 1/√(1 + (ωRC)²),  ∠H = −atan(ωRC).
    let rc = 1e-3_f64;
    for (i, &f) in freqs.iter().enumerate() {
        let omega_rc = 2.0 * std::f64::consts::PI * f * rc;
        let mag_expected   = 1.0 / (1.0 + omega_rc * omega_rc).sqrt();
        let phase_expected = -omega_rc.atan().to_degrees();

        let mag   = result.magnitude("out", i).unwrap();
        let phase = result.phase_deg("out", i).unwrap();

        let mag_rel = ((mag - mag_expected) / mag_expected).abs();
        let phase_diff = (phase - phase_expected).abs();
        assert!(mag_rel < 1e-4,
            "f={f:.3} Hz: |H|={mag:.6} expected={mag_expected:.6} ({mag_rel:.2e} rel)");
        assert!(phase_diff < 1e-3,
            "f={f:.3} Hz: ∠H={phase:.4}° expected={phase_expected:.4}° ({phase_diff:.4} diff)");
    }
}

/// Series RLC bandpass: at resonance f₀ = 1/(2π√LC), V(midpoint) magnitude
/// peaks; phase is 0 at f₀ and ±90° at f₀(1±1/(2Q)).
#[test]
fn rlc_resonance_peak_at_f0() {
    // R=1Ω, L=1mH, C=1µF → f₀ ≈ 5.033 kHz, Q ≈ √(L/C)/R ≈ 31.62.
    let net = parse_spice(
        "* RLC bandpass\n\
         V1 in 0 DC 1\n\
         R1 in m1 1\n\
         L1 m1 m2 1m\n\
         C1 m2 0 1u\n\
         .end\n"
    ).unwrap();
    let f0 = 1.0 / (2.0 * std::f64::consts::PI * (1e-3_f64 * 1e-6).sqrt());
    // Sample near resonance and far off-resonance.
    let freqs = freq_decade(1.0, 1e5, 30);
    let result = ac_analysis(&net, &freqs, None, &DeviceRegistry::new()).unwrap();

    // The capacitor-side node (m2) should peak near f0.
    let mut peak_idx = 0usize;
    let mut peak_mag = 0.0_f64;
    for (i, _) in freqs.iter().enumerate() {
        let m = result.magnitude("m2", i).unwrap();
        if m > peak_mag {
            peak_mag = m;
            peak_idx = i;
        }
    }
    let peak_freq = freqs[peak_idx];
    let f_rel = (peak_freq - f0).abs() / f0;
    assert!(f_rel < 0.1,
        "RLC peak at {peak_freq:.2} Hz, expected ≈ f₀ = {f0:.2} Hz ({:.1}% off)",
        f_rel * 100.0);
    // Q ≈ 31 → at resonance the gain should be ~Q (high).
    assert!(peak_mag > 10.0,
        "expected high Q gain at resonance; got |V(m2)| = {peak_mag:.2}");
}
