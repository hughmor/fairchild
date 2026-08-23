//! λ through the C ABI (#71): listed after a transient, fetchable by name in
//! both analyses — the same in-the-listing-and-probeable property the core
//! pins in `tests/bundles/lambda_in_signal_listings.rs`, here because this
//! front end keeps its own signal list and its own by-name series lookup, and
//! either could drift from the core's answer on its own.

use std::ffi::{CStr, CString};

use fairchild_c::{
    fc_load_string, fc_op_node, fc_run_op, fc_run_tran, fc_signal, fc_signal_count, fc_signal_name,
    fc_sim_free, fc_sim_new, FcSim, FC_OK,
};

const DECK: &str = "\
.optical_port a
.optical_port b
Xl a fc_cw_laser power_mW=1.0 wavelength_nm=1310
Xw a b fc_waveguide L_um=0.1 alpha_dB_cm=0
.op
";

const WL: f64 = 1310.0e-9;

fn loaded() -> *mut FcSim {
    let sim = fc_sim_new();
    assert!(!sim.is_null());
    let deck = CString::new(DECK).unwrap();
    assert_eq!(unsafe { fc_load_string(sim, deck.as_ptr()) }, FC_OK);
    sim
}

#[test]
fn lambda_answers_at_op_and_is_a_listed_fetchable_tran_signal() {
    let sim = loaded();

    // .op: the by-name probe the core already answered; the ABI must too.
    assert_eq!(unsafe { fc_run_op(sim) }, FC_OK);
    let name = CString::new("b_wl_0").unwrap();
    let mut v = 0.0f64;
    assert_eq!(unsafe { fc_op_node(sim, name.as_ptr(), &mut v) }, FC_OK);
    assert!((v - WL).abs() < 1e-18, "op V(b_wl_0) = {v:e}");

    // .tran: λ is in the signal list...
    assert_eq!(unsafe { fc_run_tran(sim, 1e-12, 5e-12) }, FC_OK);
    let listed: Vec<String> = (0..unsafe { fc_signal_count(sim) })
        .map(|i| {
            let p = unsafe { fc_signal_name(sim, i) };
            assert!(!p.is_null());
            unsafe { CStr::from_ptr(p) }.to_string_lossy().into_owned()
        })
        .collect();
    assert!(
        listed.iter().any(|s| s == "V(b_wl_0)"),
        "λ missing from the tran signal list: {listed:?}"
    );

    // ...and every listed signal is fetchable — the disagreement between the
    // list and the lookup is the bug this file exists to keep dead.
    let mut t_len = 0usize;
    let mut ptr: *const f64 = std::ptr::null();
    let time = CString::new("time").unwrap();
    assert_eq!(
        unsafe { fc_signal(sim, time.as_ptr(), &mut ptr, &mut t_len) },
        FC_OK
    );
    for sig in &listed {
        let c = CString::new(sig.as_str()).unwrap();
        let mut p: *const f64 = std::ptr::null();
        let mut len = 0usize;
        assert_eq!(
            unsafe { fc_signal(sim, c.as_ptr(), &mut p, &mut len) },
            FC_OK,
            "listed signal {sig} is not fetchable"
        );
        assert_eq!(len, t_len, "{sig}: series length differs from time");
    }

    // The λ series is the wavelength at every point.
    let c = CString::new("V(b_wl_0)").unwrap();
    let mut p: *const f64 = std::ptr::null();
    let mut len = 0usize;
    assert_eq!(
        unsafe { fc_signal(sim, c.as_ptr(), &mut p, &mut len) },
        FC_OK
    );
    let series = unsafe { std::slice::from_raw_parts(p, len) };
    assert!(series.iter().all(|&x| (x - WL).abs() < 1e-18), "{series:?}");

    unsafe { fc_sim_free(sim) };
}
