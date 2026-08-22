//! Time-domain noise: `.options trannoise=1`.
//!
//! Every test here is an identity between a measured variance and a closed
//! form, because that is the only way to tell a correct noise injector from a
//! plausible one — a wrong amplitude still produces a waveform that looks like
//! noise.
//!
//! The three that matter:
//!   * `kT/C`, which is bandwidth-independent and therefore step-independent;
//!   * `S_V/2h` on an unfiltered node, which pins the ZOH scaling directly;
//!   * agreement with `.noise` on a circuit with all three optical/electrical
//!     source kinds, which is what stops the two analyses drifting apart.
//!
//! **The agreement test cannot stand alone**, and that was checked rather than
//! assumed: deleting laser RIN from `NoiseSources` leaves every test in this
//! file green, because both analyses read the same list and so agree about a
//! circuit that is now missing a generator. What catches it is
//! `tests/optical_noise.rs`, which pins the `.noise` PSDs against closed forms.
//! The two files are one chain — absolute PSDs there, absolute variances and
//! frequency/time agreement here — and neither half is redundant.

use fairchild_core::{noise_analysis, tran_nr_with_registry_opts, DeviceRegistry, SimOptions};
use fairchild_parser::parse_spice;

const KB: f64 = 1.380649e-23;

fn run(src: &str, step: f64, stop: f64) -> fairchild_core::TranResult {
    let net = parse_spice(src).expect("netlist should parse");
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&net.models);
    let opts = SimOptions::from_netlist(&net);
    tran_nr_with_registry_opts(&net, step, stop, &registry, &opts).expect("transient should run")
}

/// Mean and variance of a node, discarding `skip` leading points so the
/// deterministic settling transient does not land in the statistics.
fn stats(r: &fairchild_core::TranResult, node: &str, skip: usize) -> (f64, f64) {
    let v = &r.node_voltages[node][skip..];
    let n = v.len() as f64;
    let mean = v.iter().sum::<f64>() / n;
    let var = v.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / (n - 1.0);
    (mean, var)
}

/// The canonical result: a capacitor charged through any resistor settles at
/// `⟨v²⟩ = kT/C`, independent of R — the resistor sets both the noise and the
/// bandwidth, and they cancel exactly.
///
/// This is the strongest single check available. It pins the ZOH amplitude, the
/// injection sign convention, the RC filter, and the claim that observable
/// noise is step-size independent, all in one number that contains no free
/// parameter. Both resistances must land on the same variance — and on the
/// *identical* variance to every digit, not merely a close one: with `h` scaled
/// to `τ` the two runs are the same dimensionless realisation, and R cancels
/// out of the discrete-time answer as exactly as it does out of the continuous
/// one. Equal-looking numbers here are the physics, not a copy-paste.
#[test]
fn a_capacitor_settles_at_kt_over_c() {
    const C: f64 = 1e-12;
    let expect = KB * SimOptions::default().temp_k / C;
    for r_ohm in [1e3, 4e3] {
        let tau = r_ohm * C;
        let step = tau / 50.0;
        let src = format!(
            "* kT/C\n\
             .options trannoise=1 noiseseed=7\n\
             V1 in 0 DC 0\n\
             R1 in out {r_ohm}\n\
             C1 out 0 {C}\n\
             .end\n"
        );
        // 4000 time constants: the variance of a variance estimate over N
        // independent samples is 2/N, so ~2 % here. 12 % leaves headroom for
        // one unlucky seed without hiding a real factor.
        let r = run(&src, step, 4000.0 * tau);
        let (mean, var) = stats(&r, "out", 500);
        assert!(
            mean.abs() < 3.0 * expect.sqrt(),
            "mean {mean:.3e} is offset"
        );
        let rel = (var - expect).abs() / expect;
        assert!(
            rel < 0.12,
            "R={r_ohm}: var {var:.4e} V², expected kT/C = {expect:.4e} V² (rel {rel:.3})"
        );
    }
}

/// An unfiltered node has no bandwidth of its own, so its variance IS the
/// resolved band: `∫₀^{1/2h} S_V df = S_V/2h`. That makes it the direct test of
/// the `√(S/2h)` scaling — and the one measurement that legitimately depends on
/// the timestep, which is why it is asserted at two of them.
#[test]
fn an_unfiltered_node_carries_the_whole_resolved_band() {
    const R1: f64 = 1e3;
    const R2: f64 = 1e3;
    let s_v = 4.0 * KB * SimOptions::default().temp_k * (R1 * R2 / (R1 + R2));
    for step in [1e-9, 4e-9] {
        let src = format!(
            "* white\n\
             .options trannoise=1 noiseseed=3\n\
             V1 in 0 DC 0\n\
             R1 in out {R1}\n\
             R2 out 0 {R2}\n\
             .end\n"
        );
        let r = run(&src, step, 40_000.0 * step);
        let (_, var) = stats(&r, "out", 10);
        let expect = s_v / (2.0 * step);
        let rel = (var - expect).abs() / expect;
        assert!(
            rel < 0.06,
            "h={step:.0e}: var {var:.4e} V², expected S_V/2h = {expect:.4e} V² (rel {rel:.3})"
        );
    }
}

/// The invariant that keeps the two analyses honest: the measured time-domain
/// variance equals the `.noise` PSD integrated over the resolved band.
///
/// Run three times, each biased so a **different** source carries the answer.
/// One receiver would not do: at 1 kΩ and 1 mW, shot noise is 1 % of the total,
/// so dropping it entirely would sit comfortably inside any tolerance loose
/// enough to pass. Shot beats thermal once the load voltage clears `2·V_T`
/// (≈52 mV), which is what the 100 kΩ cases arrange.
///
/// The receiver is deliberately unfiltered so `S_V(f)` is flat and the integral
/// is `S_V/2h`; a shaped response would turn this into a test of the
/// integration rule instead of a test of the two analyses agreeing.
#[test]
fn transient_noise_agrees_with_the_noise_analysis() {
    const STEP: f64 = 2e-11;
    // (label, optical power mW, RIN, load Ω) — dominant term named in the label.
    let cases = [
        ("thermal-dominated", 1.0, Some(-140.0), 1e3),
        ("shot-dominated", 0.01, None, 1e5),
        ("RIN-dominated", 0.01, Some(-120.0), 1e5),
    ];
    for (label, power_mw, rin, r_load) in cases {
        let rin = rin.map_or(String::new(), |v| format!(" rin_db_hz={v}"));
        let src = format!(
            "* direct-detection receiver\n\
             .options trannoise=1 noiseseed=11\n\
             .optical_port opt\n\
             Xlas opt fc_cw_laser power_mw={power_mw}{rin}\n\
             Xpd  opt det bias fc_photodetector responsivity=0.8 r_shunt=1e9 i_dark_a=0\n\
             Rl   det bias {r_load}\n\
             Vb   bias 0 DC 0\n\
             .end\n"
        );
        let net = parse_spice(&src).unwrap();
        let mut registry = DeviceRegistry::new();
        registry.register_builtin_models(&net.models);
        let opts = SimOptions::from_netlist(&net);
        // Flat with frequency, so any probe frequency gives the same number.
        let psd = noise_analysis(&net, &[1e9], "det", "bias", "vb", &registry, &opts)
            .expect("noise analysis")
            .onoise_psd[0];

        let r = run(&src, STEP, 60_000.0 * STEP);
        let (_, var) = stats(&r, "det", 100);
        let expect = psd / (2.0 * STEP);
        let rel = (var - expect).abs() / expect;
        assert!(
            rel < 0.06,
            "{label}: tran var {var:.4e} V² vs .noise {psd:.4e} V²/Hz over {:.3e} Hz \
             = {expect:.4e} V² (rel {rel:.3})",
            1.0 / (2.0 * STEP)
        );
    }
}

/// Off by default, and off means bit-for-bit deterministic. Every transient
/// golden in the tree depends on this.
#[test]
fn noise_is_off_by_default_and_off_is_exactly_reproducible() {
    let quiet = "* rc\nV1 in 0 PULSE(0 1 0 1n 1n 1u 2u)\nR1 in out 1k\nC1 out 0 1n\n";
    let a = run(quiet, 1e-8, 2e-6);
    let b = run(quiet, 1e-8, 2e-6);
    assert_eq!(a.node_voltages["out"], b.node_voltages["out"]);

    // And turning it on actually changes something — otherwise the assertion
    // above would pass just as well on a no-op implementation.
    let noisy = "* rc\n.options trannoise=1\nV1 in 0 PULSE(0 1 0 1n 1n 1u 2u)\n\
                 R1 in out 1k\nC1 out 0 1n\n";
    let n = run(noisy, 1e-8, 2e-6);
    assert_ne!(a.node_voltages["out"], n.node_voltages["out"]);
}

/// A noisy run is still a reproducible one: same seed, same waveform. That is
/// what makes a failure worth debugging and a sweep worth comparing.
#[test]
fn the_seed_selects_the_realisation() {
    let deck = |seed: u32| {
        format!(
            "* rc\n.options trannoise=1 noiseseed={seed}\n\
             V1 in 0 DC 0\nR1 in out 1k\nC1 out 0 1p\n"
        )
    };
    let a = run(&deck(1), 2e-11, 2e-8);
    let again = run(&deck(1), 2e-11, 2e-8);
    let other = run(&deck(2), 2e-11, 2e-8);
    assert_eq!(a.node_voltages["out"], again.node_voltages["out"]);
    assert_ne!(a.node_voltages["out"], other.node_voltages["out"]);
}

/// `noisescale` multiplies the amplitude, so the power goes as its square —
/// the knob is for pulling a deep-BER eye closure into a runnable simulation,
/// and it has to scale the way the person doing the extrapolation assumes.
#[test]
fn noisescale_is_an_amplitude_not_a_power() {
    let deck = |scale: f64| {
        format!(
            "* kT/C\n.options trannoise=1 noiseseed=5 noisescale={scale}\n\
             V1 in 0 DC 0\nR1 in out 1k\nC1 out 0 1p\n"
        )
    };
    let (_, v1) = stats(&run(&deck(1.0), 2e-11, 4e-6), "out", 500);
    let (_, v3) = stats(&run(&deck(3.0), 2e-11, 4e-6), "out", 500);
    let rel = (v3 / v1 - 9.0).abs() / 9.0;
    assert!(rel < 1e-9, "power ratio {:.4} should be 9", v3 / v1);
}

/// Adaptive step control would chase the noise and bias its spectrum, so the
/// combination is refused rather than quietly approximated.
#[test]
fn variable_step_is_refused_rather_than_approximated() {
    let src = "* rc\n.options trannoise=1 variable_step=1\n\
               V1 in 0 DC 1\nR1 in out 1k\nC1 out 0 1n\n";
    let net = parse_spice(src).unwrap();
    let opts = SimOptions::from_netlist(&net);
    let Err(err) = fairchild_core::tran_nr_with_registry_var_opts(
        &net,
        1e-9,
        1e-6,
        &DeviceRegistry::new(),
        &opts,
    ) else {
        panic!("variable step + trannoise must be an error");
    };
    assert!(
        format!("{err}").contains("variable_step=0"),
        "the error should name the fix: {err}"
    );
}
