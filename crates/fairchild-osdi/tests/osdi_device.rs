//! Integration test: OsdiDevice adapter for the test_conductance mock model.
//!
//! Loads osdi-mock, wraps its descriptor in OsdiDevice, and verifies that:
//!   - setup_model initialises gd = 1e-3 S in model memory
//!   - setup_instance writes node_mapping into instance memory
//!   - eval is a no-op for a linear element
//!   - load_residual contributes nothing to b (Jeq = 0)
//!   - load_jacobian stamps the conductance correctly into MnaMatrix

use std::path::PathBuf;

use fairchild_core::device::{Device, EvalFlags, SimContext};
use fairchild_core::mna::MnaMatrix;
use fairchild_osdi::{OsdiDevice, OsdiLibrary};

fn mock_path() -> PathBuf {
    let mut p = std::env::current_exe().expect("current_exe");
    p.pop();
    if p.ends_with("deps") {
        p.pop();
    }
    let ext = if cfg!(target_os = "macos") { "dylib" } else { "so" };
    p.push(format!("libosdi_mock.{ext}"));
    p
}

// MnaMatrix::new is private; this helper mirrors its layout for tests.
fn make_mat(size: usize) -> MnaMatrix {
    MnaMatrix {
        a: vec![vec![0.0f64; size]; size],
        b: vec![0.0f64; size],
    }
}

#[test]
fn osdi_device_conductance_stamp() {
    let path = mock_path();
    if !path.exists() {
        eprintln!("osdi-mock not found — skipping OsdiDevice test");
        return;
    }

    let lib = unsafe { OsdiLibrary::open(&path) }.expect("open mock");
    let desc_ptr = lib.descriptors().next().expect("one descriptor") as *const _;

    let ctx = SimContext::default();

    // --- setup_model ---
    let mut dev = unsafe { OsdiDevice::new(desc_ptr) };
    dev.setup_model(&ctx);

    // --- setup_instance: anode = MNA node 0, cathode = MNA node 1 ---
    dev.setup_instance(&[Some(0), Some(1)], &ctx);

    // --- eval with x = [1.0, 0.0] (1 V across the conductance) ---
    let x = vec![1.0f64, 0.0f64];
    dev.eval(&x, EvalFlags::dc(), &ctx);

    // --- load_residual: Jeq = 0 for linear element, b unchanged ---
    let mut mat = make_mat(2);
    dev.load_residual(&mut mat.b);
    assert_eq!(mat.b, vec![0.0, 0.0], "load_residual should be no-op for linear conductance");

    // --- load_jacobian: stamps gd = 1e-3 S as standard conductance ---
    dev.load_jacobian(&mut mat);

    let gd = 1e-3_f64;
    let tol = 1e-12;
    assert!((mat.a[0][0] - gd).abs() < tol, "a[0][0] = {}", mat.a[0][0]);
    assert!((mat.a[0][1] + gd).abs() < tol, "a[0][1] = {}", mat.a[0][1]);
    assert!((mat.a[1][0] + gd).abs() < tol, "a[1][0] = {}", mat.a[1][0]);
    assert!((mat.a[1][1] - gd).abs() < tol, "a[1][1] = {}", mat.a[1][1]);
}

#[test]
fn osdi_device_ground_terminal() {
    // Cathode connected to ground (NodeId = None): only the anode row/col is stamped.
    let path = mock_path();
    if !path.exists() { return; }

    let lib = unsafe { OsdiLibrary::open(&path) }.expect("open mock");
    let desc_ptr = lib.descriptors().next().unwrap() as *const _;

    let ctx = SimContext::default();
    let mut dev = unsafe { OsdiDevice::new(desc_ptr) };
    dev.setup_model(&ctx);
    dev.setup_instance(&[Some(0), None], &ctx); // cathode = ground
    dev.eval(&[1.0], EvalFlags::dc(), &ctx);

    let mut mat = make_mat(1);
    dev.load_residual(&mut mat.b);
    dev.load_jacobian(&mut mat);

    // Only a[0][0] should be stamped; off-diagonal entries don't exist.
    let gd = 1e-3_f64;
    assert!((mat.a[0][0] - gd).abs() < 1e-12, "a[0][0] = {}", mat.a[0][0]);
    assert_eq!(mat.b, vec![0.0]);
}
