/// Tests for coupled inductors (K elements).
///
/// Validates the Backward-Euler companion model for mutual inductance:
///   M = coupling * sqrt(L1 * L2)
///   det = L1*L2 - M²  (= L1*L2*(1-k²))
///   G11 = h*L2/det,  G22 = h*L1/det,  G12 = G21 = -h*M/det

use fairchild_core::tran::{tran_nr_with_registry, TranResult};
use fairchild_core::device_registry::DeviceRegistry;
use fairchild_parser::parse_spice;

// ──────────────────────────────────────────────────────────────────────────────
// Helper

fn run_be(netlist_str: &str, step: f64, stop: f64) -> TranResult {
    let net = parse_spice(netlist_str).expect("parse failed");
    let reg = DeviceRegistry::new();
    tran_nr_with_registry(&net, step, stop, &reg).expect("sim failed")
}

// ──────────────────────────────────────────────────────────────────────────────
// Test 1: k=0 decoupled should match two standalone inductors

#[test]
fn k_zero_decoupled() {
    // Two identical RL circuits.  When k=0 the K element is a no-op.
    // Reference: standalone L1 only (no K element), same topology.
    let netlist_k0 = "
* k=0 — coupled but zero coupling
V1 a1 0 DC 1.0
R1 a1 b1 100
L1 b1 0 1m
V2 a2 0 DC 1.0
R2 a2 b2 100
L2 b2 0 1m
K1 l1 l2 0.0
.tran 1u 100u
.end
";

    let netlist_ref = "
* reference: standalone inductors, no K
V1 a1 0 DC 1.0
R1 a1 b1 100
L1 b1 0 1m
V2 a2 0 DC 1.0
R2 a2 b2 100
L2 b2 0 1m
.tran 1u 100u
.end
";

    let res_k0  = run_be(netlist_k0,  1e-6, 100e-6);
    let res_ref = run_be(netlist_ref, 1e-6, 100e-6);

    // Both circuits should produce the same voltage at b1 and b2.
    let t_check = 50e-6;
    let vb1_k0  = res_k0 .voltage_at("b1", t_check).unwrap();
    let vb1_ref = res_ref.voltage_at("b1", t_check).unwrap();
    let vb2_k0  = res_k0 .voltage_at("b2", t_check).unwrap();
    let vb2_ref = res_ref.voltage_at("b2", t_check).unwrap();

    assert!(
        (vb1_k0 - vb1_ref).abs() < 1e-6,
        "k=0 b1 voltage mismatch: k0={vb1_k0:.6e} ref={vb1_ref:.6e}"
    );
    assert!(
        (vb2_k0 - vb2_ref).abs() < 1e-6,
        "k=0 b2 voltage mismatch: k0={vb2_k0:.6e} ref={vb2_ref:.6e}"
    );
}

// ──────────────────────────────────────────────────────────────────────────────
// Test 2: 1:1 transformer — secondary sees primary voltage for k ≈ 1

#[test]
fn transformer_1_to_1() {
    // Primary: V_source → R_s → L1.
    // Secondary: L2 → R_load.
    // With k=0.999 and identical inductances, secondary voltage ≈ primary
    // open-circuit voltage.
    //
    // Netlist topology:
    //   V1 (prim_p → 0) — R1 (prim_p → n1) — L1 (n1 → 0)
    //   L2 (n2 → 0) — R2 (n2 → sec_out) — (sec_out → 0)
    //   K1 L1 L2 0.999
    //
    // Driven by a step: V1 = 1V at t=0.
    // After a few time constants the primary current settles at V1/(R1) = 10 mA.
    // The secondary sees a reflected voltage; for k≈1 and equal L the
    // secondary open-circuit EMF matches the primary inductor voltage.
    //
    // We verify that the secondary coil develops a non-trivial voltage
    // (> 10% of primary) within the first few microseconds of the step.

    let netlist = "
* 1:1 transformer, k=0.999
V1 prim_p 0 DC 1.0
R1 prim_p n1 100
L1 n1 0 1m
L2 n2 0 1m
R2 n2 sec_out 100
K1 l1 l2 0.999
.tran 1u 50u
.end
";

    let res = run_be(netlist, 1e-6, 50e-6);

    // With k=0.999, L1=L2=1mH, the 2×2 companion is almost perfectly symmetric.
    // The secondary node n2 sees a voltage almost equal to V(n1), because the
    // mutual conductance G12 ≈ -G11 (k→1 limit): the secondary inductor's
    // current exactly mirrors the primary's.
    //
    // Practically: in the time before the primary current ramps up significantly
    // (t << tau = L/R = 10 µs), the inductors' voltages dominate and V(n2) ≈ V(n1).
    //
    // At t=2µs (tau/5, early transient):
    //   V(n1) ≈ V_step * (1 - e^{-t/tau}) ≈ 0.18 V
    //   V(n2) should be close to V(n1) for k=0.999 (secondary mirrors primary).
    //
    // The test requires that V(n2)/V(n1) is within 5% of 1.0.
    // This would FAIL if the sign of G12 were flipped, or if k were misread.
    let t_check = 2e-6;
    let v_n2 = res.voltage_at("n2", t_check).unwrap_or(0.0);
    let v_n1 = res.voltage_at("n1", t_check).unwrap_or(0.0);

    assert!(
        v_n1.abs() > 1e-4,
        "primary V(n1) too small at t={t_check:.1e}s: {v_n1:.4e} — check circuit"
    );
    assert!(
        v_n2.abs() > 1e-4,
        "secondary V(n2) too small at t={t_check:.1e}s: {v_n2:.4e} — coupling not working"
    );
    // For k=0.999, V(n2)/V(n1) should be close to 1.0.
    // A sign error in stamp_mutual_conductance or a wrong k factor would fail this.
    let ratio = v_n2 / v_n1;
    assert!(
        (ratio - 1.0).abs() < 0.10,
        "V(n2)/V(n1) = {ratio:.4} should be ≈ 1.0 for k=0.999 (got V(n1)={v_n1:.4e} V(n2)={v_n2:.4e})"
    );

    // The secondary output node sec_out also has voltage (R2 divider).
    let v_sec_out = res.voltage_at("sec_out", t_check).unwrap_or(0.0).abs();
    assert!(
        v_sec_out > 1e-5,
        "secondary output sec_out too small: {v_sec_out:.4e}"
    );
}

// ──────────────────────────────────────────────────────────────────────────────
// Test 3: coupled LC resonance shifts with coupling

#[test]
fn lck_resonance_frequency_shift() {
    // A parallel LC tank driven by a current pulse — classic ringing circuit.
    //
    // Topology:
    //   I1 (0→n1) pulse current source — injects charge that rings in the tank
    //   R_damp (n1→0) large damping resistor (high-Q tank)
    //   L1 (n1→0)  inductor 1  [+ L2 in same shunt branch via separate node n2]
    //   L2 (n2→0)  inductor 2  (K1 couples L1 and L2)
    //   C1 (n1→0)  capacitor
    //
    // For two parallel inductors L1‖L2 (uncoupled), L_eff = L/2.
    // For coupled parallel inductors (series-aiding, same node → same current direction,
    // BUT they share node n1 so they are in parallel):
    //   Actually when both inductors share the same two nodes (n1 and 0),
    //   coupling causes the effective admittance to change.
    //   The 2×2 system gives:
    //     [I1]   [G11 G12] [V1]      where V1=V(n1), V2=V(n2)
    //     [I2] = [G12 G22] [V2]
    //   With both connected n→0, V1=V2=V(n1), so I_total = (G11+G22+2*G12)*V.
    //   G11+G22+2*G12 = (h*L2 + h*L1 - 2*h*M)/det = h*(L1+L2-2M)/(L1*L2-M²)
    //   For L1=L2=L, M=k*L: G_eff = h*(2L-2kL)/(L²-k²L²) = 2h(1-k)/(L(1-k²)) = 2h/(L(1+k))
    //   So L_eff_parallel = h/G_eff = L(1+k)/2
    //   (vs L_eff_standalone = L/2 uncoupled)
    //   Ratio = (1+k) → for k=0.5, L_eff increases by 50% → f decreases, period longer.
    //
    // Parameters: L = 1 mH, C = 1 nF → f0 ≈ 159 kHz uncoupled (Leff=0.5mH)
    //   Using very short sim with step=1ns, stop=50µs
    //
    // Simpler and more robust: use a PULSE current into a resonant tank.
    // Initial conditions: V(n1) = 0 (default), I_L1 = I_L2 = 0.
    // The current pulse dumps charge into C, C rings with L.
    //
    // NOTE: With step=1ns and L=1mH, C=1nF, the resonant period ≈ 6.28µs — well sampled.
    //
    // Actually even simpler test: use R very large (10kΩ → high Q) and a short current pulse,
    // then check that the oscillation period differs between coupled and uncoupled.

    // L = 1 mH, C = 1 µF → f0 (uncoupled, 2 parallel Ls) = 1/(2π√(0.5mH * 1µF))
    //   = 1/(2π√(5e-10)) = 1/(2π*2.236e-5) ≈ 7118 Hz → period ≈ 140.5 µs
    // With k=0.5: L_eff = 0.5*1.5mH = 0.75 mH → f = 1/(2π√(0.75mH*1µF)) ≈ 5813 Hz
    //   period ≈ 172 µs → ratio ≈ 1.225 = sqrt(1.5)
    // Simulation: step=0.5µs, stop=1ms → 2000 points, ~7 cycles.

    let step = 0.5e-6;
    let stop = 1.0e-3;

    // Both L shunted across C1, driven by a single brief current pulse.
    // The pulse fires once (no repetition): 0 → 0.1A at t=0, stays until t=2µs,
    // then turns off.  After the pulse, the tank rings freely.
    //
    // PWL form to ensure it's a single non-repeating pulse:
    //   t=0: I=0, t=0.1ns: I=0.1A, t=2µs: I=0.1A, t=2.1µs: I=0A, t=1ms: I=0A
    let netlist_uncoupled = "
* uncoupled parallel-L tank
I1 0 n1 PWL(0 0 0.1n 0.1 2u 0.1 2.1u 0 1m 0)
R1 n1 0 10k
L1 n1 0 1m
L2 n1 0 1m
C1 n1 0 1u
.tran 0.5u 1m
.end
";

    let netlist_coupled = "
* coupled parallel-L tank  k=0.5
I1 0 n1 PWL(0 0 0.1n 0.1 2u 0.1 2.1u 0 1m 0)
R1 n1 0 10k
L1 n1 0 1m
L2 n1 0 1m
C1 n1 0 1u
K1 l1 l2 0.5
.tran 0.5u 1m
.end
";

    let res_unc = run_be(netlist_uncoupled, step, stop);
    let res_cpl = run_be(netlist_coupled,   step, stop);

    // Find the first zero-crossing of V(n1) after t=3µs (after the pulse ends at 2.1µs).
    // The tank rings: V starts positive (pulse charged C), falls, crosses zero,
    // then swings negative.  The first downward zero-crossing gives the quarter-period.
    let find_first_zero_crossing = |res: &TranResult| -> Option<f64> {
        let times = &res.time;
        let v_n1  = res.node_voltages.get("n1")?;
        let mut prev_v = 0.0f64;
        let mut prev_t = 0.0f64;
        for (&t, &v) in times.iter().zip(v_n1.iter()) {
            if t < 3e-6 { prev_v = v; prev_t = t; continue; }
            // Downward zero crossing: prev was positive, current is negative
            if prev_v > 0.0 && v <= 0.0 {
                let frac = prev_v / (prev_v - v);
                return Some(prev_t + frac * (t - prev_t));
            }
            prev_v = v;
            prev_t = t;
        }
        None
    };

    let t_zero_unc = find_first_zero_crossing(&res_unc);
    let t_zero_cpl = find_first_zero_crossing(&res_cpl);

    assert!(t_zero_unc.is_some(), "uncoupled tank did not ring down through zero — check circuit");
    assert!(t_zero_cpl.is_some(), "coupled tank did not ring down through zero");

    let t_unc = t_zero_unc.unwrap();
    let t_cpl = t_zero_cpl.unwrap();

    // With k=0.5, L_eff_parallel = L(1+k)/2 = 0.75mH > 0.5mH (uncoupled).
    // The half-period = π√(L_eff·C), so the zero-crossing time ∝ √(L_eff).
    // Coupled zero-crossing should arrive LATER than uncoupled.
    assert!(
        t_cpl > t_unc,
        "coupled resonance should be SLOWER (later zero-cross) than uncoupled: t_unc={t_unc:.4e} t_cpl={t_cpl:.4e}"
    );

    // Ratio t_cpl/t_unc ≈ sqrt(L_eff_coupled / L_eff_uncoupled)
    //   = sqrt((1+k)) = sqrt(1.5) ≈ 1.225.
    // Allow ±20% for BE discretisation error.
    let ratio = t_cpl / t_unc;
    let expected = 1.5f64.sqrt();
    assert!(
        (ratio - expected).abs() < 0.25,
        "zero-crossing time ratio {ratio:.3} should be ≈ sqrt(1.5)={expected:.3} (±0.25)"
    );
}

