//! Optical noise in `.noise`: photodetector shot noise and laser RIN.
//!
//! The circuit is the textbook direct-detection receiver — laser straight into
//! a PIN, PIN into a load resistor — for which the output noise PSD is
//!
//!     S_V = (4kT/R_L + 2q·I + RIN·I²) · |Z|²,   I = responsivity · P
//!
//! with `Z = R_L ‖ r_shunt`.  Each of the three terms is checked on its own by
//! differencing runs that turn exactly one of them off.

use fairchild_core::{noise_analysis, DeviceRegistry, SimOptions};
use fairchild_parser::parse_spice;

const Q: f64 = 1.602176634e-19;
const KB: f64 = 1.380649e-23;

const R_LOAD: f64 = 1.0e3;
const R_SHUNT: f64 = 1.0e9;
const RESPONSIVITY: f64 = 0.8;
const POWER_MW: f64 = 1.0;

/// `Z` seen by a current injected at the anode: the load in parallel with the
/// PD's own shunt.  Every term below is scaled by `Z²`, so getting it wrong
/// would move all three together — which is why the tests difference runs
/// rather than trusting one absolute number.
fn z_load() -> f64 {
    1.0 / (1.0 / R_LOAD + 1.0 / R_SHUNT)
}

/// One receiver run.  `power_mw = 0` darkens the laser (kills shot noise but
/// leaves every impedance alone); `rin` of `None` omits the parameter entirely.
fn onoise(power_mw: f64, rin_db_hz: Option<f64>, phase_deg: f64) -> f64 {
    let rin = match rin_db_hz {
        Some(v) => format!(" rin_db_hz={v}"),
        None => String::new(),
    };
    let src = format!(
        "* direct-detection receiver\n\
         .optical_port opt\n\
         Xlas opt fc_cw_laser power_mw={power_mw} phi_0_deg={phase_deg}{rin}\n\
         Xpd  opt det bias fc_photodetector responsivity={RESPONSIVITY} \
               r_shunt={R_SHUNT} i_dark_a=0\n\
         Rl   det bias {R_LOAD}\n\
         Vb   bias 0 DC 0\n\
         .end\n"
    );
    let net = parse_spice(&src).expect("netlist should parse");
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&net.models);
    let opts = SimOptions::default();
    let r = noise_analysis(&net, &[1e9], "det", "bias", "vb", &registry, &opts)
        .expect("noise analysis should run");
    r.onoise_psd[0]
}

fn photocurrent(power_mw: f64) -> f64 {
    RESPONSIVITY * power_mw * 1e-3
}

fn assert_close(got: f64, want: f64, rtol: f64, what: &str) {
    let rel = (got - want).abs() / want.abs();
    assert!(
        rel < rtol,
        "{what}: got {got:.6e} want {want:.6e} rel {rel:.2e}"
    );
}

/// Dark receiver: only the load resistor generates noise, so the whole PSD is
/// `4kT/R_L · Z²`.  This pins `Z` before anything optical is asked of it.
#[test]
fn a_dark_receiver_is_thermal_noise_alone() {
    let z = z_load();
    let expect = 4.0 * KB * SimOptions::default().temp_k / R_LOAD * z * z;
    assert_close(onoise(0.0, None, 0.0), expect, 1e-6, "thermal");
}

/// Illuminating the detector adds exactly `2q·I·Z²` — the shot noise of the
/// detected current, and nothing else.
#[test]
fn illumination_adds_shot_noise_of_the_detected_current() {
    let thermal = onoise(0.0, None, 0.0);
    let lit = onoise(POWER_MW, None, 0.0);
    let z = z_load();
    let expect = 2.0 * Q * photocurrent(POWER_MW) * z * z;
    assert_close(lit - thermal, expect, 1e-6, "shot");
}

/// Shot noise is linear in optical power, unlike RIN.  Doubling the power
/// doubles the shot term; a quadratic term hiding in the shot path would show
/// up here as 4×.
#[test]
fn shot_noise_is_linear_in_optical_power() {
    let thermal = onoise(0.0, None, 0.0);
    let one = onoise(POWER_MW, None, 0.0) - thermal;
    let two = onoise(2.0 * POWER_MW, None, 0.0) - thermal;
    assert_close(two / one, 2.0, 1e-6, "shot power scaling");
}

/// RIN adds `RIN·I²·Z²` on top of shot — quadratic in power, which is what
/// makes it the floor a high-power link cannot buy its way out of.
#[test]
fn rin_adds_a_term_quadratic_in_the_photocurrent() {
    const RIN_DB: f64 = -150.0;
    let rin = 10f64.powf(RIN_DB / 10.0);
    let z = z_load();

    for scale in [1.0, 2.0] {
        let p = scale * POWER_MW;
        let without = onoise(p, None, 0.0);
        let with = onoise(p, Some(RIN_DB), 0.0);
        let i = photocurrent(p);
        assert_close(with - without, rin * i * i * z * z, 1e-6, "rin");
    }
}

/// The RIN generator drives BOTH field wires from one intensity fluctuation,
/// so its contribution cannot depend on the emission phase.
///
/// This is the test that pins the coherent sum: treating the `re` and `im`
/// taps as two independent sources gives the same answer at 0° (the `im` tap
/// is dark) and exactly **half** at 45°, where the power splits evenly.
#[test]
fn rin_is_independent_of_the_emission_phase() {
    const RIN_DB: f64 = -150.0;
    let at = |phi: f64| onoise(POWER_MW, Some(RIN_DB), phi) - onoise(POWER_MW, None, phi);
    let reference = at(0.0);
    for phi in [30.0, 45.0, 90.0, 145.0] {
        assert_close(at(phi), reference, 1e-6, &format!("rin at {phi}°"));
    }
}

/// A laser with no `rin_db_hz` is noiseless — 0 dB/Hz would be a RIN of 1/Hz,
/// which is not a defensible default for an unset parameter.
#[test]
fn rin_is_off_unless_asked_for() {
    let z = z_load();
    let expect = (4.0 * KB * SimOptions::default().temp_k / R_LOAD
        + 2.0 * Q * photocurrent(POWER_MW))
        * z
        * z;
    assert_close(onoise(POWER_MW, None, 0.0), expect, 1e-6, "thermal + shot");
    assert!(onoise(POWER_MW, Some(-150.0), 0.0) > onoise(POWER_MW, None, 0.0));
}
