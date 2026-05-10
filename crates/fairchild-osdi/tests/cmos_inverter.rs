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
/// Pre-condition: va-models/build/nmos_l1.osdi and va-models/build/pmos_l1.osdi must exist.
/// Build with:
///   DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib \
///   openvaf-r va-models/nmos_l1.va --output va-models/build/nmos_l1.osdi
///   openvaf-r va-models/pmos_l1.va --output va-models/build/pmos_l1.osdi

use std::path::PathBuf;
use std::sync::Arc;

use fairchild_core::{dc_op_nr_with_registry, tran_nr_with_registry, DeviceRegistry};
use fairchild_osdi::OsdiLibrary;
use fairchild_parser::parse_spice;

fn nmos_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../va-models/build/nmos_l1.osdi")
}

fn pmos_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../va-models/build/pmos_l1.osdi")
}

fn load_cmos_registry() -> Option<DeviceRegistry> {
    let np = nmos_path();
    let pp = pmos_path();

    if !np.exists() || !pp.exists() {
        eprintln!(
            "Skipping: OSDI models not found.\n\
             Build with:\n  \
             DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib \\\n  \
             openvaf-r va-models/nmos_l1.va --output va-models/build/nmos_l1.osdi\n  \
             openvaf-r va-models/pmos_l1.va --output va-models/build/pmos_l1.osdi"
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
