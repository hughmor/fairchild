//! `$abstime` must read the transient clock, not zero.
//!
//! `OsdiSimInfo.abstime` was hardcoded to 0.0, so any Verilog-A model reading
//! `$abstime` — an arbitrary-waveform source, a drift or ageing term, anything
//! time-parametrised — silently saw t = 0 for the whole run. `SimContext`
//! already carried `time_s` (the waveguide group delay uses it); `eval` just
//! never passed it on.
//!
//! Runs against `osdi-mock`, whose `eval` records the `abstime` it was handed
//! into model memory, so this needs no OpenVAF and runs in CI.

use std::path::PathBuf;

use fairchild_core::device::{Device, EvalFlags, SimContext};
use fairchild_osdi::{OsdiDevice, OsdiLibrary};
use osdi_mock::ABSTIME_OFFSET;

fn mock_path() -> PathBuf {
    let mut p = std::env::current_exe().expect("current_exe");
    p.pop();
    if p.ends_with("deps") {
        p.pop();
    }
    let ext = if cfg!(target_os = "macos") {
        "dylib"
    } else {
        "so"
    };
    p.push(format!("libosdi_mock.{ext}"));
    p
}

#[test]
fn eval_forwards_the_simulation_clock_as_abstime() {
    let path = mock_path();
    if !path.exists() {
        eprintln!("osdi-mock not found at {path:?}; run `cargo build -p osdi-mock`.");
        return;
    }

    let lib = unsafe { OsdiLibrary::open(&path) }.expect("dlopen failed");
    let desc = lib.descriptors().next().unwrap() as *const _;
    let mut dev = unsafe { OsdiDevice::new(desc) };

    let model_ptr = dev.model_ptr_raw();
    let recorded = || unsafe { *((model_ptr as *const u8).add(ABSTIME_OFFSET) as *const f64) };

    let mut ctx = SimContext::default();
    dev.setup_model(&ctx);
    dev.setup_instance(&[None, None], &ctx);

    // DC: the clock is 0 and Verilog-A expects to see 0 there.
    dev.eval(&[0.0, 0.0], EvalFlags::dc(), &ctx);
    assert_eq!(recorded(), 0.0, "DC eval should report t = 0");

    for t in [1e-9, 3.25e-6, 0.5] {
        ctx.time_s = t;
        dev.eval(&[0.0, 0.0], EvalFlags::tran(), &ctx);
        assert_eq!(
            recorded(),
            t,
            "eval did not forward SimContext::time_s ({t}) as OSDI abstime"
        );
    }
}
