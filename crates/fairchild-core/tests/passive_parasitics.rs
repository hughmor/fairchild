/// Integration tests for passive parasitic expansion (Item E).
///
/// Tests that `rser=`, `cpar=`, `esr=`, `esl=`, `rpar=` on R/L/C elements
/// expand to the same circuit as the equivalent explicit sub-network.

use fairchild_core::tran::{tran_nr_with_registry, TranResult};
use fairchild_core::device_registry::DeviceRegistry;
use fairchild_parser::parse_spice;

fn run_be(netlist_str: &str, step: f64, stop: f64) -> TranResult {
    let net = parse_spice(netlist_str).expect("parse failed");
    let reg = DeviceRegistry::new();
    tran_nr_with_registry(&net, step, stop, &reg).expect("sim failed")
}

fn voltage_at(res: &TranResult, node: &str, t: f64) -> f64 {
    res.voltage_at(node, t)
        .unwrap_or_else(|| panic!("node '{node}' not found at t={t:.2e}"))
}

// ──────────────────────────────────────────────────────────────────────────────
// Inductor rser: L1 in out 1m rser=10  ≡  L1 in __l1_rn 1m + R __l1_rn out 10

#[test]
fn inductor_rser_equivalent_to_explicit_series_r() {
    // LC tank with inductor ESR injected via parasitic param.
    let with_param = "
* Inductor ESR via rser=
V1  in  0   PULSE(0 1 0 1n 1n 50u 200u)
L1  in  out 1m  rser=10
C1  out 0   1u
.tran 1u 200u
.end
";
    // Exactly equivalent circuit using explicit elements.
    let explicit = "
* Explicit series R for ESR
V1  in  0   PULSE(0 1 0 1n 1n 50u 200u)
L1  in  __l1_rn  1m
Resr __l1_rn out 10
C1  out 0   1u
.tran 1u 200u
.end
";
    let res_p = run_be(with_param, 1e-6, 200e-6);
    let res_e = run_be(explicit,   1e-6, 200e-6);

    for t in [50e-6, 100e-6, 150e-6, 200e-6] {
        let vp = voltage_at(&res_p, "out", t);
        let ve = voltage_at(&res_e, "out", t);
        assert!(
            (vp - ve).abs() < 1e-6,
            "inductor rser: V(out) mismatch at t={t:.1e}: param={vp:.6e} explicit={ve:.6e}"
        );
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Capacitor esr: C1 mid 0 1u esr=10  ≡  R __c1_esrn 0 10 + C1 mid __c1_esrn 1u

#[test]
fn capacitor_esr_equivalent_to_explicit_series_r() {
    let with_param = "
* Capacitor ESR via esr=
V1  in  0   PULSE(0 1 0 1n 1n 50u 200u)
L1  in  mid 1m
C1  mid 0   1u  esr=10
.tran 1u 200u
.end
";
    let explicit = "
* Explicit series R for ESR (between L and C)
V1  in  0   PULSE(0 1 0 1n 1n 50u 200u)
L1  in  mid 1m
Resr mid __c1_esrn 10
C1  __c1_esrn 0  1u
.tran 1u 200u
.end
";
    let res_p = run_be(with_param, 1e-6, 200e-6);
    let res_e = run_be(explicit,   1e-6, 200e-6);

    for t in [50e-6, 100e-6, 150e-6, 200e-6] {
        let vp = voltage_at(&res_p, "mid", t);
        let ve = voltage_at(&res_e, "mid", t);
        assert!(
            (vp - ve).abs() < 1e-6,
            "capacitor esr: V(mid) mismatch at t={t:.1e}: param={vp:.6e} explicit={ve:.6e}"
        );
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Capacitor rpar: leakage resistance across capacitor

#[test]
fn capacitor_rpar_equivalent_to_explicit_parallel_r() {
    // Cap discharging through rpar=1k, C=1µF → τ=1ms
    let with_param = "
* Cap discharging through rpar=1k
V1  in  0   PULSE(1 0 0 1n 1n 5m 20m)
L1  in  out 1u
C1  out 0   1u  rpar=1k
.tran 10u 5m
.end
";
    let explicit = "
* Same circuit, rpar explicit
V1  in  0   PULSE(1 0 0 1n 1n 5m 20m)
L1  in  out 1u
C1  out 0   1u
Rleak out 0 1k
.tran 10u 5m
.end
";
    let res_p = run_be(with_param, 10e-6, 5e-3);
    let res_e = run_be(explicit,   10e-6, 5e-3);

    for t in [1e-3, 2e-3, 3e-3] {
        let vp = voltage_at(&res_p, "out", t);
        let ve = voltage_at(&res_e, "out", t);
        assert!(
            (vp - ve).abs() < 1e-3,
            "capacitor rpar transient: V(out) at t={t:.1e}: param={vp:.4e} explicit={ve:.4e}"
        );
    }
}
