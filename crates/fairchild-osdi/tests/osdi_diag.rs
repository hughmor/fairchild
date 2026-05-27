/// Diagnostic: verify OSDI Jacobian by running a modified DC OP with print output.
use std::ffi::CStr;
use std::path::PathBuf;
use std::sync::Arc;

use fairchild_core::device::{Device, EvalFlags, SimContext};
use fairchild_core::mna::{CircuitTopology, stamp_netlist};
use fairchild_core::DeviceRegistry;
use fairchild_osdi::{OsdiDevice, OsdiLibrary};
use fairchild_parser::parse_spice;
use indexmap::IndexMap;

fn osdi_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../legacy/va-models/build/diode_shockley.osdi")
}

#[test]
fn osdi_jacobian_value_at_zero_voltage() {
    let path = osdi_path();
    if !path.exists() { return; }

    let lib = Arc::new(unsafe { OsdiLibrary::open(&path) }.expect("dlopen"));

    // Print descriptor metadata once.
    for (i, desc) in lib.descriptors().enumerate() {
        let name = unsafe { CStr::from_ptr(desc.name) }.to_str().unwrap();
        eprintln!("desc[{i}]: {name}");
        eprintln!("  num_nodes={} num_terminals={} num_resistive_jac={}",
            desc.num_nodes, desc.num_terminals, desc.num_resistive_jacobian_entries);
        eprintln!("  node_mapping_offset={}", desc.node_mapping_offset);
        let entries = unsafe {
            std::slice::from_raw_parts(desc.jacobian_entries, desc.num_jacobian_entries as usize)
        };
        for (j, e) in entries.iter().enumerate() {
            eprintln!("  jac_entry[{j}]: ({}, {})", e.nodes.node_1, e.nodes.node_2);
        }
    }

    // Manually run one NR step for the 1-node circuit Ib=1mA → b → D1 → GND.
    let netlist = parse_spice(
        "* diag\nIb 0 b 1m\nD1 b 0 diode_shockley\n.op\n.end\n"
    ).unwrap();

    let ctx = SimContext::default();
    let mut registry = DeviceRegistry::new();
    lib.register_into(&mut registry);

    // Build topology manually to access node count.
    let topo = CircuitTopology::build(&netlist);
    eprintln!("MNA size: {}  n_nodes={}", topo.size, topo.n_nodes());

    // Build devices via registry.
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();
    let x0 = vec![0.0f64; topo.size];

    let mut mat = stamp_netlist(&topo, &netlist, 0.0, &empty, &empty);
    eprintln!("mat.a before device stamp: {:?}", mat.a);
    eprintln!("mat.b before device stamp: {:?}", mat.b);

    // Create device the same way dc_op_nr_with_registry does.
    let anode_idx = topo.node_index.get("b").copied();
    let mut dev = OsdiDevice::from_library(Arc::clone(&lib), 0).unwrap();
    dev.setup_model(&ctx);
    dev.setup_instance(&[anode_idx, None], &ctx);

    dev.eval(&x0, EvalFlags::dc(), &ctx);
    eprintln!("eval called with x={x0:?}");

    dev.load_residual(&mut mat.b);
    eprintln!("mat.b after load_residual: {:?}", mat.b);

    dev.load_jacobian(&mut mat);
    eprintln!("mat.a after load_jacobian: {:?}", mat.a);

    // mat.a[0][0] must be non-zero (the diode's conductance at Vd=0).
    assert!(
        mat.a[0][0].abs() > 0.0,
        "mat.a[0][0] = {} — Jacobian was not stamped!",
        mat.a[0][0]
    );
    eprintln!("gd(Vd=0) = {:.4e}  (expected ≈ {:.4e})", mat.a[0][0], 1e-14 / (1.380649e-23 * 300.15 / 1.602176634e-19));
}
