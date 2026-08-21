//! A model that calls `$limit()` must not jump through a null pointer.
//!
//! OpenVAF compiles `$limit(V(a,b), "pnjlim", vt, vcrit)` into a call through
//! the library's exported `OSDI_LIM_TABLE`, whose entries ship with
//! `func_ptr = NULL` because the *simulator* is expected to install its own
//! implementations. The call is guarded by the `ENABLE_LIM` eval flag, which
//! fairchild sets whenever a model declares limit state — so once limiting is
//! opted into, a null entry is a jump to address 0 on the first `eval`. Every
//! foundry compact model uses `$limit`, so this covers the whole PDK class.
//!
//! The invariant under test is therefore *no entry is ever left null*, not
//! merely "pnjlim gets installed": OpenVAF does not validate limiter names at
//! all — it forwards whatever string literal the model wrote — so fairchild
//! cannot enumerate them, and anything unrecognised has to degrade to an
//! identity limiter rather than a crash.
//!
//! `lim_probe.va` asks for all three cases in one module, and the compiler puts
//! them in the table: the name and arity we implement, a name nobody
//! implements, and a known name at an arity we do not. The retired mock built
//! that table by hand, which meant the shape being tested was our own guess at
//! it.

use std::path::Path;

use fairchild_osdi::ffi::{FnPnjlim, OsdiLimFunction};
use fairchild_osdi::OsdiLibrary;

mod common;

/// Read the table out of the library we just loaded. `dlopen` of an
/// already-loaded library returns the same handle and the same table, so this
/// observes exactly what `OsdiLibrary::open` wrote.
unsafe fn lim_table(path: &Path) -> &'static [OsdiLimFunction] {
    let c = std::ffi::CString::new(path.to_str().unwrap()).unwrap();
    let h = libc::dlopen(c.as_ptr(), libc::RTLD_LAZY | libc::RTLD_LOCAL);
    assert!(!h.is_null(), "dlopen failed");
    let tbl = libc::dlsym(h, c"OSDI_LIM_TABLE".as_ptr()) as *const OsdiLimFunction;
    let len = *(libc::dlsym(h, c"OSDI_LIM_TABLE_LEN".as_ptr()) as *const u32);
    assert!(
        !tbl.is_null(),
        "a model calling $limit must export OSDI_LIM_TABLE"
    );
    std::slice::from_raw_parts(tbl, len as usize)
}

#[test]
fn no_lim_table_entry_is_left_null() {
    let Some(path) = common::compiled("lim_probe") else {
        return;
    };
    let _lib = unsafe { OsdiLibrary::open(&path) }.expect("dlopen failed");
    let table = unsafe { lim_table(&path) };

    // One entry per `$limit` call in the model, in source order.
    assert_eq!(table.len(), 3, "lim_probe.va makes three $limit calls");
    let names: Vec<(String, u32)> = table
        .iter()
        .map(|e| {
            let n = unsafe { std::ffi::CStr::from_ptr(e.name) }
                .to_string_lossy()
                .into_owned();
            (n, e.num_args)
        })
        .collect();
    assert_eq!(
        names,
        vec![
            ("pnjlim".to_string(), 2),
            ("no_such_limiter".to_string(), 2),
            ("pnjlim".to_string(), 1),
        ],
        "the compiler forwards names and arities verbatim; the installer has to \
         cope with all three"
    );

    // THE invariant: not one entry may be left null, whether or not fairchild
    // recognises it. OpenVAF guards the call with ENABLE_LIM, which fairchild
    // sets whenever a model declares limit state — so a null here is a jump to
    // address 0 on the first eval.
    for (i, entry) in table.iter().enumerate() {
        let name = unsafe { std::ffi::CStr::from_ptr(entry.name) };
        assert!(
            !entry.func_ptr.is_null(),
            "OSDI_LIM_TABLE[{i}] ({name:?}, {} args) left null",
            entry.num_args
        );
    }

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

    // [1] an unimplemented name and [2] a known name at the wrong arity both
    // get the identity fallback: pass the proposed value straight through and
    // leave `check` alone, since nothing was limited. Degrading to "no
    // limiting" costs convergence robustness; it does not change the answer,
    // and it does not kill the process.
    for i in [1usize, 2] {
        let fallback: FnPnjlim = unsafe { std::mem::transmute(table[i].func_ptr) };
        check = false;
        // The value that real pnjlim would have compressed.
        let v = unsafe { fallback(false, &mut check, 5.0, 0.65, vt, vcrit) };
        assert_eq!(v, 5.0, "entry {i} should be the identity limiter");
        assert!(!check, "entry {i} limited nothing but vetoed convergence");
        // Even on the init call, where real pnjlim returns vcrit.
        let v = unsafe { fallback(true, &mut check, 5.0, 0.65, vt, vcrit) };
        assert_eq!(
            v, 5.0,
            "entry {i} should be the identity limiter on init too"
        );
        assert!(!check);
    }
}

/// And the model has to actually run: a `$limit` call means `num_states > 0`,
/// which turns on `ENABLE_LIM` and makes every one of those entries live.
#[test]
fn a_limiting_model_solves() {
    let Some(path) = common::compiled("lim_probe") else {
        return;
    };
    let deck = format!(
        "* a model that limits\n\
         .osdi {}\n\
         V1 in 0 DC 1\n\
         Xl in 0 lim_probe\n\
         .op\n\
         .end\n",
        path.display()
    );
    let netlist = fairchild_parser::parse_spice(&deck).expect("parse");
    let mut registry = fairchild_core::DeviceRegistry::new();
    registry.register_builtin_models(&netlist.models);
    fairchild_osdi::load_libraries(
        &netlist.osdi_paths,
        &netlist.va_sources,
        None,
        &Default::default(),
        &mut registry,
    )
    .expect("load");
    let r = fairchild_core::dc_op_nr_with_registry(&netlist, &registry).expect("DC OP failed");

    // Below vcrit nothing is limited, so all three $limit calls pass V through
    // and the current is 1 mS · 1 V.
    let i = r.vsrc_current("v1").expect("v1");
    assert!(
        (i.abs() - 1e-3).abs() < 1e-9,
        "supply current {i} A, want 1 mA"
    );
}
