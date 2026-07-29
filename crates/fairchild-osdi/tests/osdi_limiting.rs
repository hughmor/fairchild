//! A model that calls `$limit()` must not jump through a null pointer.
//!
//! OpenVAF compiles `$limit(V(a,b), "pnjlim", vt, vcrit)` into a call through
//! the library's exported `OSDI_LIM_TABLE`, whose entries ship with
//! `func_ptr = NULL` because the *simulator* is expected to install its own
//! implementations. fairchild never did, and the call is emitted
//! unconditionally — so loading any model that used `$limit` and evaluating it
//! was an immediate SIGSEGV with no diagnostic. Every foundry compact model
//! uses `$limit`, so that was the whole PDK class.
//!
//! Runs against `osdi-mock`, which exports the same table shape, so this needs
//! no OpenVAF and runs in CI.

use std::path::{Path, PathBuf};

use fairchild_osdi::ffi::{FnPnjlim, OsdiLimFunction};
use fairchild_osdi::OsdiLibrary;

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

/// Read the table out of the library we just loaded. `dlopen` of an
/// already-loaded library returns the same handle and the same table, so this
/// observes exactly what `OsdiLibrary::open` wrote.
unsafe fn lim_table(path: &Path) -> &'static [OsdiLimFunction] {
    let c = std::ffi::CString::new(path.to_str().unwrap()).unwrap();
    let h = libc::dlopen(c.as_ptr(), libc::RTLD_LAZY | libc::RTLD_LOCAL);
    assert!(!h.is_null(), "dlopen failed");
    let tbl = libc::dlsym(h, c"OSDI_LIM_TABLE".as_ptr()) as *const OsdiLimFunction;
    let len = *(libc::dlsym(h, c"OSDI_LIM_TABLE_LEN".as_ptr()) as *const u32);
    assert!(!tbl.is_null(), "mock should export OSDI_LIM_TABLE");
    std::slice::from_raw_parts(tbl, len as usize)
}

#[test]
fn loading_installs_pnjlim() {
    let path = mock_path();
    if !path.exists() {
        eprintln!("osdi-mock not found at {path:?}; run `cargo build -p osdi-mock`.");
        return;
    }

    let _lib = unsafe { OsdiLibrary::open(&path) }.expect("dlopen failed");
    let table = unsafe { lim_table(&path) };
    assert_eq!(table.len(), 1);
    assert!(
        !table[0].func_ptr.is_null(),
        "OSDI_LIM_TABLE[0] (\"pnjlim\") is still null — evaluating a model that \
         calls $limit would jump to address 0"
    );

    // And it has to be pnjlim, not just any non-null pointer.
    let pnjlim: FnPnjlim = unsafe { std::mem::transmute(table[0].func_ptr) };
    let (vt, vcrit) = (0.02585, 0.6145);
    let mut check = false;

    // Untouched: below vcrit, so no limiting and no convergence veto.
    let v = unsafe { pnjlim(false, &mut check, 0.3, 0.29, vt, vcrit) };
    assert_eq!(v, 0.3);
    assert!(
        !check,
        "pnjlim vetoed convergence on a step it did not limit"
    );

    // Limited: above vcrit and a step much larger than 2*vt, so it must
    // compress logarithmically and land between the two.
    let (vold, vnew) = (0.65, 5.0);
    let v = unsafe { pnjlim(false, &mut check, vnew, vold, vt, vcrit) };
    assert!(check, "pnjlim limited the step but did not set *check");
    assert!(
        vold < v && v < vnew,
        "limited voltage {v} is not between {vold} and {vnew}"
    );
    let expected = vold + vt * ((vnew - vold) / vt + 1.0).ln();
    assert!((v - expected).abs() < 1e-12, "{v} != {expected}");

    // init: no previous iterate worth limiting against, so start at vcrit.
    check = false;
    let v = unsafe { pnjlim(true, &mut check, 5.0, 0.0, vt, vcrit) };
    assert_eq!(v, vcrit);
    assert!(check);
}
