/// CMOS inverter integration test using Level 1 MOSFET OSDI models.
///
/// Circuit: VDD=5V, NMOS + PMOS complementary pair, no explicit load.
///   MN: drain=out, gate=in, source=0,  bulk=0,   model=nmos_l1
///   MP: drain=out, gate=in, source=vdd, bulk=vdd, model=pmos_l1
///
/// DC OP checks:
///   Vin=0V   → NMOS off, PMOS on  → Vout ≈ VDD = 5V
///   Vin=5V   → PMOS off, NMOS on  → Vout ≈ 0V
///
/// Pre-condition: legacy/va-models/build/nmos_l1.osdi and legacy/va-models/build/pmos_l1.osdi must exist.
/// Build with:
///   DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib \
///   openvaf-r legacy/va-models/nmos_l1.va --output legacy/va-models/build/nmos_l1.osdi
///   openvaf-r legacy/va-models/pmos_l1.va --output legacy/va-models/build/pmos_l1.osdi

use std::path::PathBuf;
use std::sync::Arc;

use fairchild_core::{dc_op_nr_with_registry, tran_nr_with_registry, tran_nr_with_registry_var, DeviceRegistry};
use fairchild_osdi::OsdiLibrary;
use fairchild_parser::parse_spice;

fn nmos_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../legacy/va-models/build/nmos_l1.osdi")
}

fn pmos_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../legacy/va-models/build/pmos_l1.osdi")
}

fn load_cmos_registry() -> Option<DeviceRegistry> {
    let np = nmos_path();
    let pp = pmos_path();

    if !np.exists() || !pp.exists() {
        eprintln!(
            "Skipping: OSDI models not found.\n\
             Build with:\n  \
             DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib \\\n  \
             openvaf-r legacy/va-models/nmos_l1.va --output legacy/va-models/build/nmos_l1.osdi\n  \
             openvaf-r legacy/va-models/pmos_l1.va --output legacy/va-models/build/pmos_l1.osdi"
        );
        return None;
    }

    let nlib = Arc::new(unsafe { OsdiLibrary::open(&np) }.expect("dlopen nmos_l1.osdi failed"));
    let plib = Arc::new(unsafe { OsdiLibrary::open(&pp) }.expect("dlopen pmos_l1.osdi failed"));

    let mut registry = DeviceRegistry::new();
    nlib.register_into(&mut registry);
    plib.register_into(&mut registry);

    Some(registry)
}

// CMOS inverter netlist template; Vin is substituted as a DC value.
fn cmos_netlist(vin_v: f64) -> String {
    format!(
        "* CMOS inverter Level 1 OSDI\n\
         VDD vdd 0 DC 5\n\
         Vin in 0 DC {vin_v}\n\
         MN out in 0 0 nmos_l1\n\
         MP out in vdd vdd pmos_l1\n\
         .op\n\
         .end\n"
    )
}

/// Vin=0 → NMOS off, PMOS on → Vout ≈ VDD = 5V.
#[test]
fn cmos_inverter_input_low() {
    let Some(registry) = load_cmos_registry() else { return; };

    let netlist = parse_spice(&cmos_netlist(0.0)).unwrap();
    let result = dc_op_nr_with_registry(&netlist, &registry).expect("DC OP failed (Vin=0)");

    let vout = result.node_voltage("out").unwrap();
    assert!(
        vout > 4.95,
        "Vin=0: Vout={vout:.4} V — expected ≈ 5V (PMOS on, NMOS off)"
    );
    println!("Vin=0: Vout = {vout:.6} V");
}

/// Vin=VDD → PMOS off, NMOS on → Vout ≈ 0V.
#[test]
fn cmos_inverter_input_high() {
    let Some(registry) = load_cmos_registry() else { return; };

    let netlist = parse_spice(&cmos_netlist(5.0)).unwrap();
    let result = dc_op_nr_with_registry(&netlist, &registry).expect("DC OP failed (Vin=5V)");

    let vout = result.node_voltage("out").unwrap();
    assert!(
        vout < 0.05,
        "Vin=5V: Vout={vout:.4} V — expected ≈ 0V (NMOS on, PMOS off)"
    );
    println!("Vin=5: Vout = {vout:.6} V");
}

/// Switching transient: PULSE input from 0→5V at t=1ns.
/// Verifies Vout transitions from VDD toward 0 after the edge.
#[test]
fn cmos_inverter_switching_transient() {
    let Some(registry) = load_cmos_registry() else { return; };

    // Load capacitor at output (1 pF) to give a measurable time constant.
    // Without it, the gate-overlap caps (2 fF) give τ ≈ 0.4 ps which is too fast.
    let netlist_str =
        "* CMOS inverter transient\n\
         VDD vdd 0 DC 5\n\
         Vin in 0 PULSE(0 5 1n 10p 10p 10n 20n)\n\
         MN out in 0 0 nmos_l1\n\
         MP out in vdd vdd pmos_l1\n\
         CL out 0 1p\n\
         .tran 10p 3n\n\
         .end\n";

    let netlist = parse_spice(netlist_str).unwrap();
    let result = tran_nr_with_registry(&netlist, 10e-12, 3e-9, &registry)
        .expect("transient simulation failed");

    // Before edge: Vout should be high.
    let v_before = result.voltage_at("out", 0.5e-9).unwrap();
    assert!(
        v_before > 4.5,
        "Before edge: Vout={v_before:.4} V — expected ≈ 5V"
    );

    // After edge + settling: Vout should be low.
    let v_after = result.voltage_at("out", 2.9e-9).unwrap();
    assert!(
        v_after < 0.5,
        "After edge: Vout={v_after:.4} V — expected ≈ 0V"
    );

    println!(
        "Switching transient: V(out) at 0.5ns = {v_before:.4}V, at 2.9ns = {v_after:.4}V"
    );
}

/// Same switching transient via the variable-step solver.
/// This exercises the reactive Jacobian at realistic alpha (1/h, not 1.0).
#[test]
fn cmos_inverter_switching_transient_var_step() {
    let Some(registry) = load_cmos_registry() else { return; };

    let netlist_str =
        "* CMOS inverter transient (variable-step)\n\
         VDD vdd 0 DC 5\n\
         Vin in 0 PULSE(0 5 1n 10p 10p 10n 20n)\n\
         MN out in 0 0 nmos_l1\n\
         MP out in vdd vdd pmos_l1\n\
         CL out 0 1p\n\
         .tran 10p 3n\n\
         .end\n";

    let netlist = parse_spice(netlist_str).unwrap();
    let result = tran_nr_with_registry_var(&netlist, 10e-12, 3e-9, &registry)
        .expect("variable-step transient simulation failed");

    let v_before = result.voltage_at("out", 0.5e-9).unwrap();
    assert!(
        v_before > 4.5,
        "Var-step before edge: Vout={v_before:.4} V — expected ≈ 5V"
    );

    let v_after = result.voltage_at("out", 2.9e-9).unwrap();
    assert!(
        v_after < 0.5,
        "Var-step after edge: Vout={v_after:.4} V — expected ≈ 0V"
    );

    println!(
        "Var-step switching: V(out) at 0.5ns = {v_before:.4}V, at 2.9ns = {v_after:.4}V"
    );
}

/// Verify that write_jacobian_array_react covers ALL reactive entries (including
/// drain-drain and drain-gate which are in the resistive-only prefix but also have
/// reactive contributions). Entries with react_ptr_off != u32::MAX are reactive;
/// write_jacobian_array_react iterates them in entry-array order.
#[test]
fn nmos_reactive_jacobian_covers_all_entries() {
    use fairchild_core::device::{Device, EvalFlags, SimContext};
    use fairchild_osdi::OsdiDevice;

    let np = nmos_path();
    if !np.exists() { return; }
    let nlib = Arc::new(unsafe { OsdiLibrary::open(&np) }.expect("dlopen nmos"));

    let desc = nlib.descriptors().next().expect("no descriptors");
    let n_total = desc.num_jacobian_entries as usize;
    let n_react = desc.num_reactive_jacobian_entries as usize;
    let entries = unsafe {
        std::slice::from_raw_parts(desc.jacobian_entries, n_total)
    };

    let ctx = SimContext::default();
    let mut dev = OsdiDevice::from_library(Arc::clone(&nlib), 0).unwrap();
    dev.setup_model(&ctx);
    dev.setup_instance(&[Some(0), Some(1), None, None], &ctx);

    let x = vec![2.5f64, 2.5f64];
    dev.eval(&x, EvalFlags::tran(), &ctx);

    let mut react_buf = vec![0.0f64; n_react];
    unsafe {
        let f = desc.write_jacobian_array_react
            .expect("write_jacobian_array_react must be present");
        f(dev.inst_ptr_raw(), dev.model_ptr_raw(), react_buf.as_mut_ptr());
    }

    // react_buf values are in entry-array order for entries with react_ptr_off != MAX.
    let reactive_entries: Vec<usize> = entries.iter().enumerate()
        .filter(|(_, e)| e.react_ptr_off != u32::MAX)
        .map(|(j, _)| j)
        .collect();

    assert_eq!(reactive_entries.len(), n_react,
        "reactive entry count from react_ptr_off must match num_reactive_jacobian_entries");

    // Entry[0] = (drain,drain): reactive contribution from Cgd must be non-zero.
    let entry0_react_idx = reactive_entries.iter().position(|&j| j == 0)
        .expect("entry[0] (drain,drain) must be in the reactive set");
    assert!(
        react_buf[entry0_react_idx].abs() > 0.0,
        "write_jacobian_array_react[{}] for (drain,drain) = {} — expected non-zero Cgd contribution",
        entry0_react_idx, react_buf[entry0_react_idx]
    );

    // Gate-gate entry must have the largest (Cgs+Cgd) contribution.
    let entry4_react_idx = reactive_entries.iter().position(|&j| j == 4)
        .expect("entry[4] (gate,gate) must be in the reactive set");
    assert!(
        react_buf[entry4_react_idx] > react_buf[entry0_react_idx],
        "gate-gate reactive ({:.2e}) should be larger than drain-drain ({:.2e})",
        react_buf[entry4_react_idx], react_buf[entry0_react_idx]
    );

    println!(
        "Reactive Jacobian: drain-drain={:.2e}, gate-gate={:.2e}",
        react_buf[entry0_react_idx], react_buf[entry4_react_idx]
    );
}
