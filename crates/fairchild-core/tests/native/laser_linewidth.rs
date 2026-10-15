//! A linewidth is measured by coherence, not by whether a waveform looks noisy.
use fairchild_core::{tran_nr_with_registry_opts, DeviceRegistry, SimOptions, TranResult};
use fairchild_parser::parse_spice;

const POWER: f64 = 1e-3;
const LINEWIDTH: f64 = 1e9;

fn run(params: &str, options: &str, step: f64, count: usize) -> TranResult {
    let deck = format!(
        "* laser phase diffusion\n.optical_port a\n.optical_port b\n\
         Xa a fc_cw_laser power_mW=1 phase_deg=47 {params}\n\
         Xb b fc_cw_laser power_mW=1 phase_deg=-31 {params}\n\
         .options {options}\n"
    );
    let net = parse_spice(&deck).unwrap();
    let opts = SimOptions::from_netlist(&net);
    tran_nr_with_registry_opts(
        &net,
        step,
        count as f64 * step,
        &DeviceRegistry::new(),
        &opts,
    )
    .expect("laser transient")
}

fn increments(r: &TranResult, name: &str) -> Vec<f64> {
    let re = &r.node_voltages[&format!("{name}_re_0")];
    let im = &r.node_voltages[&format!("{name}_im_0")];
    (1..re.len())
        .map(|i| {
            (im[i] * re[i - 1] - re[i] * im[i - 1]).atan2(re[i] * re[i - 1] + im[i] * im[i - 1])
        })
        .collect()
}

#[test]
fn linewidth_has_the_wiener_variance_and_exponential_coherence() {
    for h in [1e-11, 2e-11] {
        let r = run("linewidth_hz=1G", "trannoise=1 noiseseed=23", h, 60_000);
        let da = increments(&r, "a");
        let db = increments(&r, "b");
        let variance = da.iter().map(|v| v * v).sum::<f64>() / da.len() as f64;
        let want = std::f64::consts::TAU * LINEWIDTH * h;
        assert!(
            (variance / want - 1.0).abs() < 0.025,
            "h={h:e}: phase-increment variance {variance}, expected {want}"
        );
        let cross = da.iter().zip(&db).map(|(a, b)| a * b).sum::<f64>() / da.len() as f64;
        assert!(
            (cross / want).abs() < 0.025,
            "independent lasers share phase increments: {cross}"
        );
        let re = &r.node_voltages["a_re_0"];
        let im = &r.node_voltages["a_im_0"];
        for i in 0..re.len() {
            assert!(
                ((re[i] * re[i] + im[i] * im[i]) / POWER - 1.0).abs() < 1e-12,
                "pure phase diffusion must preserve power"
            );
        }
        for tau in [40e-12, 160e-12, 320e-12] {
            let lag = (tau / h).round() as usize;
            let measured = (lag..re.len())
                .map(|i| (re[i] * re[i - lag] + im[i] * im[i - lag]) / POWER)
                .sum::<f64>()
                / (re.len() - lag) as f64;
            let expected = (-std::f64::consts::PI * LINEWIDTH * tau).exp();
            eprintln!("h={h:e}, tau={tau:e}: coherence {measured:.5}, expected {expected:.5}");
            assert!(
                (measured - expected).abs() < 0.04,
                "wrong coherence decay at {tau:e}"
            );
        }
    }
}

#[test]
fn linewidth_obeys_the_noise_switch_seed_and_amplitude_scale() {
    let quiet = run("", "", 1e-11, 512);
    let off = run("linewidth_hz=1G", "trannoise=0", 1e-11, 512);
    assert_eq!(quiet.node_voltages, off.node_voltages);
    let reference = run("rin_db_hz=-120", "trannoise=1 noiseseed=5", 1e-11, 512);
    // Pre-linewidth splitmix64/Box-Muller stream: two RIN draws per step,
    // one per laser. An extra draw even for default linewidth=0 changes it.
    let legacy_re = [
        0.02161485597081109,
        0.02490777505566749,
        0.018341359914379986,
        0.0211564743837041,
    ];
    for (i, want) in legacy_re.iter().enumerate() {
        assert!(
            (reference.node_voltages["a_re_0"][i + 1] - want).abs() < 1e-15,
            "zero/default linewidth changed the legacy RIN random stream at step {}",
            i + 1
        );
    }
    let zero = run(
        "linewidth_hz=0 rin_db_hz=-120",
        "trannoise=1 noiseseed=5",
        1e-11,
        512,
    );
    assert_eq!(
        reference.node_voltages, zero.node_voltages,
        "zero linewidth must preserve existing RIN draws"
    );
    let a = run("linewidth_hz=1G", "trannoise=1 noiseseed=5", 1e-11, 512);
    let b = run("linewidth_hz=1G", "trannoise=1 noiseseed=5", 1e-11, 512);
    let c = run("linewidth_hz=1G", "trannoise=1 noiseseed=6", 1e-11, 512);
    assert_eq!(
        a.node_voltages, b.node_voltages,
        "each run must reset phase and RNG state"
    );
    assert!(
        a.node_voltages != c.node_voltages,
        "different seeds must change the phase path"
    );
    let scaled = run(
        "linewidth_hz=1G",
        "trannoise=1 noiseseed=5 noisescale=2",
        1e-11,
        512,
    );
    for (a, b) in increments(&a, "a").iter().zip(increments(&scaled, "a")) {
        assert!(
            (b - 2.0 * a).abs() < 1e-12,
            "noisescale must scale phase increments"
        );
    }
    let muted = run("linewidth_hz=1G", "trannoise=1 noisescale=0", 1e-11, 512);
    assert_eq!(quiet.node_voltages, muted.node_voltages);
}

/// Absolute one-sided small-signal phase PSD: the imaginary field of a real
/// carrier has S_im = P linewidth/(pi f^2). Rotating the carrier splits this
/// coherently across the two quadratures; a direct detector cancels it.
#[test]
fn the_noise_analysis_has_the_same_phase_diffusion_strength() {
    use fairchild_core::noise_analysis;
    for phi in [0.0f64, 47.0, 131.0] {
        let deck = format!(
            ".optical_port a\nXa a fc_cw_laser power_mW=1 phase_deg={phi} linewidth_hz=1G\n\
             Vref ref 0 DC 0 AC 1\nRref ref 0 1k\n"
        );
        let net = parse_spice(&deck).unwrap();
        let opts = SimOptions::default();
        let registry = DeviceRegistry::new();
        let freqs = [1e9, 2e9, 8e9];
        let r = noise_analysis(&net, &freqs, "a_im_0", "0", "vref", &registry, &opts).unwrap();
        for (i, f) in freqs.iter().enumerate() {
            let want =
                POWER * LINEWIDTH / (std::f64::consts::PI * f * f) * phi.to_radians().cos().powi(2);
            assert!(
                (r.onoise_psd[i] / want - 1.0).abs() < 1e-10,
                "wrong phase PSD at {f:e} Hz, phi={phi}"
            );
        }
        assert!(noise_analysis(&net, &[0.0], "a_im_0", "0", "vref", &registry, &opts).is_err());
        let coherent = parse_spice(&deck.replace("linewidth_hz=1G", "linewidth_hz=0")).unwrap();
        let dc =
            noise_analysis(&coherent, &[0.0], "a_im_0", "0", "vref", &registry, &opts).unwrap();
        assert_eq!(
            dc.onoise_psd[0], 0.0,
            "preserve the coherent source's existing DC noise limit"
        );
        let detector =
            "Xpd a det 0 fc_photodetector responsivity=1 i_dark=0 r_shunt=1G\nRl det 0 1k\n";
        let with = parse_spice(&format!("{deck}{detector}")).unwrap();
        let without = parse_spice(&format!(
            "{}{detector}",
            deck.replace("linewidth_hz=1G", "linewidth_hz=0")
        ))
        .unwrap();
        let a = noise_analysis(&with, &freqs, "det", "0", "vref", &registry, &opts).unwrap();
        let b = noise_analysis(&without, &freqs, "det", "0", "vref", &registry, &opts).unwrap();
        for (a, b) in a.onoise_psd.iter().zip(&b.onoise_psd) {
            assert!(
                (a / b - 1.0).abs() < 1e-8,
                "a direct detector must reject pure phase noise"
            );
        }
    }
}

#[test]
fn invalid_linewidth_is_refused_instead_of_silently_clamped() {
    use fairchild_core::{models::NativeCwLaser, Device};
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        let mut laser = NativeCwLaser::new();
        assert!(laser.set_real_param("linewidth_hz", value));
        assert!(laser.validate().unwrap_err().contains("linewidth_hz"));
    }
}

/// RIN must remain radial while the phase wanders. Leaving its taps at the
/// initial carrier angle converts half of the amplitude noise to phase noise.
#[test]
fn rin_rotates_with_the_diffusing_carrier() {
    let h = 1e-11;
    let r = run(
        "linewidth_hz=1G rin_db_hz=-130",
        "trannoise=1 noiseseed=19",
        h,
        20_000,
    );
    let re = &r.node_voltages["a_re_0"];
    let im = &r.node_voltages["a_im_0"];
    let amp: Vec<f64> = re
        .iter()
        .zip(im)
        .skip(1)
        .map(|(re, im)| re.hypot(*im))
        .collect();
    let mean = amp.iter().sum::<f64>() / amp.len() as f64;
    let variance = amp.iter().map(|a| (a - mean).powi(2)).sum::<f64>() / amp.len() as f64;
    let expected = 1e-13 * POWER / (8.0 * h);
    assert!(
        (variance / expected - 1.0).abs() < 0.04,
        "radial RIN variance {variance:e}, expected {expected:e}"
    );
    assert!(
        (mean / POWER.sqrt() - 1.0).abs() < 0.002,
        "RIN mean shifted"
    );
}
