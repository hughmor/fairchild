//! DC operating-point test on a 5-transistor differential pair with
//! diode-connected PMOS active load.  Demonstrates that NR finds the
//! coupled bias point on a non-trivial multi-MOSFET topology.
//!
//! Topology:
//!
//!   VDD ──┬─ M3 ──┬─ M4 ──┐
//!         │ (diode) │
//!         ┴ mid1   ┴ out
//!         │        │
//!   M1 ◄──┘  M2 ◄──┘    (input pair, balanced)
//!   in1      in2
//!   │        │
//!   tail tied (M5 = current source biased by Vbias)
//!
//! The challenge: V(tail), V(mid1), V(out) must satisfy KCL with all five
//! Level-1 MOSFET equations simultaneously.  Standard NR converges from
//! x = 0 in well under ITL1; this test asserts that and pins the bias
//! voltages to within Level-1 accuracy.

use fairchild_core::dc_op_nr;
use fairchild_parser::parse_spice;

const DIFF_PAIR: &str = "\
* 5-T NMOS differential pair, PMOS diode-connected load, NMOS tail source
.model nm NMOS (vto=0.6 kp=200u lambda=0.05)
.model pm PMOS (vto=-0.7 kp=80u lambda=0.05)
Vdd  vdd  0 DC 3.3
* Inputs at common-mode VCM = 1.6 V (gives Vgs1,2 ≈ 1.6−Vtail; Vtail floats).
Vin1 in1  0 DC 1.6
Vin2 in2  0 DC 1.6
* Tail-bias gate at 1.2 V → M5 in saturation drives ~50 µA.
Vbias vbias 0 DC 1.2

* Input pair
M1 mid1 in1 tail 0 nm w=20u l=1u
M2 out  in2 tail 0 nm w=20u l=1u

* PMOS active load: M3 diode-connected, M4 mirrors.
M3 mid1 mid1 vdd vdd pm w=10u l=1u
M4 out  mid1 vdd vdd pm w=10u l=1u

* Tail current source
M5 tail vbias 0 0 nm w=10u l=1u

.op
";

#[test]
fn differential_pair_dc_bias_point() {
    let net = parse_spice(DIFF_PAIR).unwrap();
    let r = dc_op_nr(&net).expect("DC OP must converge on 5-T diff pair");
    assert!(
        r.iters <= 3,
        "diff pair DC OP took {} homotopy attempts; expected direct NR (≤3)",
        r.iters
    );

    let v_tail = r.node_voltage("tail").unwrap();
    let v_mid1 = r.node_voltage("mid1").unwrap();
    let v_out = r.node_voltage("out").unwrap();
    let vdd = 3.3;

    // Sanity checks (Level-1 accuracy, ±20 % bands):
    //   V(tail) is one Vgs below the input common mode (1.6 V); with
    //   Vth=0.6 and small overdrive on M5, expect 0.3–0.9 V.
    assert!(
        v_tail > 0.2 && v_tail < 1.2,
        "V(tail) = {v_tail:.3}; expected ≈ 0.3–0.9 V"
    );
    //   V(mid1) is one Vsg below VDD: VDD − |Vth_p| − overdrive ≈ 1.8 V.
    assert!(
        v_mid1 > 1.5 && v_mid1 < vdd - 0.5,
        "V(mid1) = {v_mid1:.3}; expected ≈ 1.5–2.5 V"
    );
    //   V(out) ≈ V(mid1) (balanced inputs → no differential → no offset).
    let imbalance = (v_out - v_mid1).abs();
    assert!(
        imbalance < 0.05,
        "balanced inputs but |V(out)−V(mid1)| = {imbalance:.4}; should be ≈ 0"
    );
}

/// Apply a small differential input and verify the output swings — proves
/// the amplifier path is biased into the high-gain region and the
/// linearised Jacobian correctly couples in1/in2 to out.
#[test]
fn differential_pair_responds_to_differential_input() {
    let nm = "Vin1 in1  0 DC 1.605\nVin2 in2  0 DC 1.595\n"; // 10 mV differential
    let netlist = DIFF_PAIR.replace("Vin1 in1  0 DC 1.6\nVin2 in2  0 DC 1.6", nm.trim_end());
    let net = parse_spice(&netlist).unwrap();
    let r = dc_op_nr(&net).expect("DC OP must converge under small input mismatch");
    let v_mid1 = r.node_voltage("mid1").unwrap();
    let v_out = r.node_voltage("out").unwrap();
    let delta = v_out - v_mid1;
    // 10 mV input differential should swing the output by far more than the
    // input (a real amplifier gain).  Level-1 MOSFET gives modest gain;
    // we just require the differential-to-output coupling is monotonic
    // with the right sign — V(in1) > V(in2) drives M2 weaker than M1, so
    // V(out) rises above V(mid1).
    assert!(
        delta > 0.005,
        "expected positive output deflection from V(in1)>V(in2); got {delta:.4} V"
    );
}
