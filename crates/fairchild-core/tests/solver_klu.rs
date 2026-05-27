//! End-to-end test: `.options solver=klu` directive routes through to
//! the SuiteSparse KLU backend and produces the same DC operating point
//! as the dense / faer-sparse paths.
//!
//! Gated on the `klu` cargo feature.

#![cfg(feature = "klu")]

use fairchild_core::device_registry::DeviceRegistry;
use fairchild_core::newton::dc_op_nr_with_registry_opts;
use fairchild_core::options::SimOptions;
use fairchild_core::solver::SolverKind;
use fairchild_parser::parse_spice;

#[test]
fn klu_solves_resistor_divider() {
    // V1 1V → 1k → out → 1k → 0   ⇒  V(out) = 0.5 V
    let net =
        parse_spice("* divider\nV1 in 0 DC 1\nR1 in out 1k\nR2 out 0 1k\n.op\n.end\n").unwrap();
    let reg = DeviceRegistry::new();
    let mut opts = SimOptions::from_netlist(&net);
    opts.solver = SolverKind::Klu;
    let r = dc_op_nr_with_registry_opts(&net, &reg, &opts).unwrap();
    let v_out = r.node_voltage("out").unwrap();
    assert!(
        (v_out - 0.5).abs() < 1e-9,
        "V(out) = {v_out} (expected 0.5)"
    );
}

#[test]
fn klu_options_directive_selects_klu_backend() {
    // Parser routes `.options solver=klu` into `SimOptions.solver`.
    let net = parse_spice(
        "* divider with directive\nV1 in 0 DC 2\nR1 in out 1k\nR2 out 0 1k\n\
         .options solver=klu\n.op\n.end\n",
    )
    .unwrap();
    let opts = SimOptions::from_netlist(&net);
    assert_eq!(opts.solver, SolverKind::Klu);

    let reg = DeviceRegistry::new();
    let r = dc_op_nr_with_registry_opts(&net, &reg, &opts).unwrap();
    let v_out = r.node_voltage("out").unwrap();
    assert!(
        (v_out - 1.0).abs() < 1e-9,
        "V(out) = {v_out} (expected 1.0)"
    );
}

#[test]
fn klu_transient_diode_rc_matches_dense() {
    // Exercises the symbolic/numeric split: one DC OP solve + hundreds
    // of NR iterations across timesteps, all reusing the same KLU
    // symbolic factorisation through `klu_refactor`.  Compare the final
    // V(out) at t=stop against the dense path.
    use fairchild_core::tran::tran_nr_with_registry_var_opts;

    let net = parse_spice(
        "* RC + diode rectifier\n\
         Vin in 0 PULSE(0 1 0 1n 1n 10n 20n)\n\
         R1 in n1 1k\nD1 n1 out myd\nR2 out 0 10k\nC1 out 0 1n\n\
         .model myd D (Is=1e-14 N=1)\n.tran 0.5n 30n\n.end\n",
    )
    .unwrap();
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_diodes(&net.models);

    let mut opts_klu = SimOptions::from_netlist(&net);
    opts_klu.solver = SolverKind::Klu;
    let r_klu = tran_nr_with_registry_var_opts(&net, 0.5e-9, 30e-9, &reg, &opts_klu).unwrap();

    let mut opts_dense = SimOptions::from_netlist(&net);
    opts_dense.solver = SolverKind::Dense;
    let r_dense = tran_nr_with_registry_var_opts(&net, 0.5e-9, 30e-9, &reg, &opts_dense).unwrap();

    let v_klu_last = r_klu
        .node_voltages
        .get("out")
        .unwrap()
        .last()
        .copied()
        .unwrap();
    let v_dense_last = r_dense
        .node_voltages
        .get("out")
        .unwrap()
        .last()
        .copied()
        .unwrap();
    assert!(
        (v_klu_last - v_dense_last).abs() < 1e-6,
        "KLU transient diverged from Dense: V(out)_klu={v_klu_last} V(out)_dense={v_dense_last}"
    );
}

#[test]
fn klu_matches_dense_on_nonlinear_diode_circuit() {
    // Forward-biased diode through 1 kΩ — KLU and Dense must agree.
    let net = parse_spice(
        "* diode\nV1 in 0 DC 1.0\nR1 in n1 1k\nD1 n1 0 myd\n\
         .model myd D (Is=1e-14 N=1)\n.op\n.end\n",
    )
    .unwrap();
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_diodes(&net.models);

    let mut opts_klu = SimOptions::from_netlist(&net);
    opts_klu.solver = SolverKind::Klu;
    let r_klu = dc_op_nr_with_registry_opts(&net, &reg, &opts_klu).unwrap();

    let mut opts_dense = SimOptions::from_netlist(&net);
    opts_dense.solver = SolverKind::Dense;
    let r_dense = dc_op_nr_with_registry_opts(&net, &reg, &opts_dense).unwrap();

    let v_klu = r_klu.node_voltage("n1").unwrap();
    let v_dense = r_dense.node_voltage("n1").unwrap();
    assert!(
        (v_klu - v_dense).abs() < 1e-9,
        "KLU vs Dense: V(n1) klu={v_klu} dense={v_dense}"
    );
}
