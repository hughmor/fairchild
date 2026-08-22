//! Regression tests for fc_thermal_ps_rc — thermal phase shifter with a
//! first-order RC.  Demonstrates path B of the hybrid reactive-state
//! design: the device declares 1 extra MNA state row for T(t) and stamps
//! its BE-discretised state equation directly in load_jacobian_tran.

use fairchild_core::{
    dc_op_nr_with_registry, tran_nr_with_registry_opts, DeviceRegistry, SimOptions,
};
use fairchild_parser::parse_spice;

/// DC: with no time derivative the state equation reduces to T = P, so the
/// L2 device should match fc_thermal_ps's DC output exactly when driven
/// to the same V_h.
#[test]
fn thermal_ps_rc_dc_matches_l1() {
    let l1 = "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xth ch0 out0 vh 0 fc_thermal_ps r_heater=1k p_pi=10m
Vh vh 0 DC 1.0
.op
";
    let l2 = "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xth ch0 out0 vh 0 fc_thermal_ps_rc r_heater=1k p_pi=10m tau_th=10u
Vh vh 0 DC 1.0
.op
";
    let r1 = dc_op_nr_with_registry(&parse_spice(l1).unwrap(), &DeviceRegistry::new()).unwrap();
    let r2 = dc_op_nr_with_registry(&parse_spice(l2).unwrap(), &DeviceRegistry::new()).unwrap();
    let v1_re = r1.node_voltage("out0_re_0").unwrap();
    let v1_im = r1.node_voltage("out0_im_0").unwrap();
    let v2_re = r2.node_voltage("out0_re_0").unwrap();
    let v2_im = r2.node_voltage("out0_im_0").unwrap();
    assert!(
        (v1_re - v2_re).abs() < 1e-6,
        "L1 vs L2 (DC) out.re mismatch: {v1_re} vs {v2_re}"
    );
    assert!((v1_im - v2_im).abs() < 1e-6);
}

/// Transient: drive V_h with a step.  The L2 device's optical output
/// (read at out0_re_0) should NOT change instantaneously — it lags by
/// roughly tau_th.  Compare to L1 which has no lag.
#[test]
fn thermal_ps_rc_transient_lags_l1_by_tau() {
    let l1 = "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xth ch0 out0 vh 0 fc_thermal_ps r_heater=1k p_pi=10m
Vh vh 0 PULSE(0 1.0 1u 100n 100n 100u 200u)
.tran 1u 30u
.options method=be
";
    let l2 = "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xth ch0 out0 vh 0 fc_thermal_ps_rc r_heater=1k p_pi=10m tau_th=10u
Vh vh 0 PULSE(0 1.0 1u 100n 100n 100u 200u)
.tran 1u 30u
.options method=be
";
    let opts = SimOptions::from_netlist(&parse_spice(l1).unwrap());
    let r1 = tran_nr_with_registry_opts(
        &parse_spice(l1).unwrap(),
        1e-6,
        30e-6,
        &DeviceRegistry::new(),
        &opts,
    )
    .unwrap();
    let r2 = tran_nr_with_registry_opts(
        &parse_spice(l2).unwrap(),
        1e-6,
        30e-6,
        &DeviceRegistry::new(),
        &opts,
    )
    .unwrap();
    // At t ≈ 2 µs (1 µs after the edge, well below tau_th = 10 µs):
    //   L1 has fully responded (its phase = π·P/P_pi instantaneous).
    //   L2 has only partially responded (≈ 1 − exp(−1/10) ≈ 9.5% of final).
    let probe_t = 2.0e-6;
    let probe_idx = r1.time.iter().position(|&t| t >= probe_t).unwrap();
    let phi_l1 = r1.node_voltages.get("out0_re_0").unwrap()[probe_idx];
    let phi_l2 = r2.node_voltages.get("out0_re_0").unwrap()[probe_idx];
    // Magnitude of the SVEA out.re changes as cos(φ) does — the L2 should
    // be CLOSER to the input amplitude (less phase shift) than L1.
    let in_amp = 1e-3_f64.sqrt();
    let dl1 = (phi_l1 - in_amp).abs();
    let dl2 = (phi_l2 - in_amp).abs();
    assert!(
        dl2 < dl1,
        "L2 should be less perturbed than L1 at t = 1 µs after edge \
         (tau_th = 10 µs > t).  l1 dev={dl1:.4} l2 dev={dl2:.4}"
    );
}

/// Sanity: after t » tau_th, the L2 device should converge to the same
/// steady-state output as L1.
#[test]
fn thermal_ps_rc_settles_to_l1_steady_state() {
    let l1 = "\
.optical_port ch0
.optical_port out0
Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xth ch0 out0 vh 0 fc_thermal_ps r_heater=1k p_pi=10m
Vh vh 0 DC 1.0
.tran 5u 100u
.options method=be
";
    let l2 = l1.replace(
        "fc_thermal_ps r_heater",
        "fc_thermal_ps_rc tau_th=10u r_heater",
    );
    let opts = SimOptions::from_netlist(&parse_spice(l1).unwrap());
    let r1 = tran_nr_with_registry_opts(
        &parse_spice(l1).unwrap(),
        5e-6,
        100e-6,
        &DeviceRegistry::new(),
        &opts,
    )
    .unwrap();
    let r2 = tran_nr_with_registry_opts(
        &parse_spice(&l2).unwrap(),
        5e-6,
        100e-6,
        &DeviceRegistry::new(),
        &opts,
    )
    .unwrap();
    let v1 = *r1.node_voltages.get("out0_re_0").unwrap().last().unwrap();
    let v2 = *r2.node_voltages.get("out0_re_0").unwrap().last().unwrap();
    assert!(
        (v1 - v2).abs() < 1e-3,
        "at t = 100 µs (10·tau_th), L2 should match L1 steady state; v1={v1} v2={v2}"
    );
}
