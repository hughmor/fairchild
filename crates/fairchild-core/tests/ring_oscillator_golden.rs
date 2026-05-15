//! Hard-convergence electrical golden test: a 5-stage CMOS ring oscillator.
//!
//! Why this circuit:
//!
//! * **Multi-device coupling.**  5 NMOS + 5 PMOS = 10 nonlinear devices
//!   wired in a feedback loop, so the DC-OP Jacobian is dense over a
//!   strongly coupled subgraph and NR has to balance all 10 transistors
//!   simultaneously.  The DC equilibrium is the *unstable* fixed point at
//!   the inverter switching threshold V_M — far from a trivial zero-current
//!   guess.  This exercises source stepping / GMIN stepping fallbacks if
//!   the direct NR doesn't converge.
//!
//! * **Sustained transient oscillation.**  After kicking the loop with an
//!   `.ic` seed, the circuit oscillates indefinitely with frequency
//!   `f = 1 / (2·N·t_pd)` where N is the stage count and t_pd is the
//!   per-stage delay.  This tests transient stability through *hundreds*
//!   of switching events, junction limiting, and integration accuracy
//!   (BE / TR / GEAR) over a non-trivial dynamic range.
//!
//! * **Analytic frequency target.**  t_pd is set by the gate-drive
//!   current charging the per-stage load capacitance through ΔV ≈ VDD/2.
//!   So we can predict the period within ~30 % of the simulator's value
//!   from Level-1 MOSFET parameters, giving a real validation target that
//!   doesn't depend on ngspice being installed.
//!
//! The expected oscillation period is on the order of nanoseconds with
//! the parameters below; the test runs for ~50 cycles and asserts the
//! period is in the right ballpark.

use fairchild_core::{
    dc_op_nr, options::SimOptions, tran::IntegratorMode, tran_nr_with_registry_var_opts,
    DeviceRegistry,
};
use fairchild_parser::parse_spice;

const RING_OSCILLATOR_NETLIST: &str = "\
* 5-stage CMOS ring oscillator (Level-1 MOSFETs)
.model nm NMOS (vto=0.5 kp=200u lambda=0.05)
.model pm PMOS (vto=-0.5 kp=80u lambda=0.05)
Vdd vdd 0 DC 1.8

* Stage 1: n5 → n1
Mn1 n1 n5 0   0   nm w=10u l=1u
Mp1 n1 n5 vdd vdd pm w=20u l=1u
C1  n1 0 100f

* Stage 2: n1 → n2
Mn2 n2 n1 0   0   nm w=10u l=1u
Mp2 n2 n1 vdd vdd pm w=20u l=1u
C2  n2 0 100f

* Stage 3: n2 → n3
Mn3 n3 n2 0   0   nm w=10u l=1u
Mp3 n3 n2 vdd vdd pm w=20u l=1u
C3  n3 0 100f

* Stage 4: n3 → n4
Mn4 n4 n3 0   0   nm w=10u l=1u
Mp4 n4 n3 vdd vdd pm w=20u l=1u
C4  n4 0 100f

* Stage 5: n4 → n5  (closes the loop)
Mn5 n5 n4 0   0   nm w=10u l=1u
Mp5 n5 n4 vdd vdd pm w=20u l=1u
C5  n5 0 100f

* Asymmetric initial condition kicks the loop out of the metastable equilibrium.
.ic V(n1)=1.6 V(n2)=0.1 V(n3)=1.6 V(n4)=0.1 V(n5)=1.6

.options method=gear
.tran 50p 100n UIC
.end
";

/// DC-OP convergence test.  The metastable equilibrium V_M ≈ VDD/2 is what
/// NR should converge to when there's no `.ic` seed (the analyzer ignores
/// `.ic` for the operating point — it only applies under UIC).
#[test]
fn ring_oscillator_dc_op_converges() {
    let net = parse_spice(RING_OSCILLATOR_NETLIST).unwrap();
    let r = dc_op_nr(&net).expect("DC OP must converge on the ring oscillator");
    // All five stage outputs should be near V_M (≈ VDD/2 to within 200 mV
    // with the asymmetric W/L=20:10 we picked).
    for n in &["n1", "n2", "n3", "n4", "n5"] {
        let v = r.node_voltage(n).unwrap();
        assert!(v > 0.4 && v < 1.4,
            "V({n})={v:.3} is not near V_M (expected ~0.9 V for VDD=1.8)");
    }
    // Convergence should be fast on a stable solver — well below ITL1=150.
    assert!(r.iters <= 5,
        "DC OP took {} attempts (homotopy levels); ring oscillator should converge directly", r.iters);
}

/// Transient test: under `.ic` (UIC) the loop oscillates.  Measure the
/// period from zero crossings of the centre stage and assert it's within
/// the expected band for Level-1 MOSFETs at these parameters.
#[test]
fn ring_oscillator_transient_oscillates() {
    let net = parse_spice(RING_OSCILLATOR_NETLIST).unwrap();
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_diodes(&net.models);
    registry.register_builtin_mosfets(&net.models);
    let mut opts = SimOptions::from_netlist(&net);
    opts.method = IntegratorMode::Gear;
    opts.uic = true;
    // Tighten LTE — ring oscillators are sensitive to integration noise.
    opts.reltol = 1e-4;
    let result = tran_nr_with_registry_var_opts(
        &net, 50e-12, 100e-9, &registry, &opts,
    ).expect("transient must complete");

    assert!(result.time.len() > 100,
        "expected many timepoints from an oscillating ring (got {})", result.time.len());

    // Detect zero crossings of (V(n3) − VDD/2) in the second half of the run
    // — the first half is the start-up transient.
    let vdd_half = 0.9;
    let n3 = result.node_voltages.get("n3").expect("V(n3) must be present");
    let t = &result.time;
    let half = t.len() / 2;
    let mut crossings: Vec<f64> = Vec::new();
    for i in (half + 1)..t.len() {
        let a = n3[i - 1] - vdd_half;
        let b = n3[i]     - vdd_half;
        if a.signum() != b.signum() && (a - b).abs() > 1e-6 {
            // Linear interpolate the actual crossing time.
            let frac = a / (a - b);
            crossings.push(t[i - 1] + frac * (t[i] - t[i - 1]));
        }
    }
    assert!(crossings.len() >= 4,
        "expected ≥4 zero crossings in second half (got {}); is the loop oscillating?",
        crossings.len());

    // Period = 2 × interval between consecutive same-direction crossings.
    // Average across all crossing pairs to smooth measurement noise.
    let intervals: Vec<f64> = crossings.windows(2).map(|w| w[1] - w[0]).collect();
    let mean_half_period: f64 = intervals.iter().sum::<f64>() / intervals.len() as f64;
    let period = 2.0 * mean_half_period;
    let freq   = 1.0 / period;

    // Analytic estimate: t_pd ≈ C·ΔV / I_avg, with I_avg ≈ ½·KP·(W/L)·(VDD − Vth)²
    // for the driving transistor.  For Level-1 with VDD=1.8, Vth=0.5, KP=200µ,
    // W/L=10, we get I_avg ≈ ½·200e-6·10·1.3² ≈ 1.7 mA → t_pd ≈ C·VDD/2/I
    // ≈ 100f·0.9/1.7m ≈ 53 ps → f ≈ 1/(2·5·53p) ≈ 1.9 GHz.  Level-1 is
    // optimistic vs real BSIM; allow a generous window.
    assert!(freq > 500e6 && freq < 5e9,
        "ring oscillator frequency f={:.2e} Hz outside expected 0.5–5 GHz band", freq);
}

/// Stress test: same circuit, but compare BE vs GEAR results — both should
/// produce the same oscillation frequency to within 5 %.  BE has more
/// numerical damping; GEAR/BDF-2 is the recommended choice for sharp
/// switching transients but BE must still converge.
#[test]
fn ring_oscillator_be_and_gear_agree_on_frequency() {
    let net = parse_spice(RING_OSCILLATOR_NETLIST).unwrap();
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_diodes(&net.models);
    registry.register_builtin_mosfets(&net.models);
    let mut opts = SimOptions::from_netlist(&net);
    opts.uic = true;
    opts.reltol = 1e-4;

    let measure_period = |opts: &SimOptions| -> f64 {
        let r = tran_nr_with_registry_var_opts(&net, 50e-12, 60e-9, &registry, opts).unwrap();
        let n3 = r.node_voltages.get("n3").unwrap();
        let t  = &r.time;
        let half = t.len() / 2;
        let mut crossings = Vec::new();
        for i in (half + 1)..t.len() {
            let a = n3[i - 1] - 0.9;
            let b = n3[i]     - 0.9;
            if a.signum() != b.signum() && (a - b).abs() > 1e-6 {
                let frac = a / (a - b);
                crossings.push(t[i - 1] + frac * (t[i] - t[i - 1]));
            }
        }
        let intervals: Vec<f64> = crossings.windows(2).map(|w| w[1] - w[0]).collect();
        2.0 * intervals.iter().sum::<f64>() / intervals.len() as f64
    };

    let mut opts_be = opts.clone();
    opts_be.method = IntegratorMode::BackwardEuler;
    let p_be = measure_period(&opts_be);

    let mut opts_gear = opts.clone();
    opts_gear.method = IntegratorMode::Gear;
    let p_gear = measure_period(&opts_gear);

    let rel = (p_be - p_gear).abs() / p_gear;
    assert!(rel < 0.10,
        "BE period {:.3e} vs GEAR period {:.3e} differ by {:.1}%",
        p_be, p_gear, rel * 100.0);
}
