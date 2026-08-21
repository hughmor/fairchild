//! `.va` source in, registered device out — the whole path fairchild now owns.
//!
//! The compiler here is a shell stub that answers `--version`, implements
//! `--print-expansion` as `cat`, and "compiles" by copying a prepared `.osdi`
//! to the requested output. Stubbing it is what makes the cache testable: the
//! stub counts its own invocations, so a hit and a miss are distinguishable,
//! which they would not be if a real compile ran each time.
//!
//! The artefact it copies is real, compiled once from `tests/models` by the
//! installed compiler — so `dlopen` and descriptor registration are exercised
//! against the genuine article. Without a compiler these tests skip, like every
//! other test in this directory.
//!
//! What this deliberately does not prove: that OpenVAF compiles any particular
//! model correctly. That is the other tests' job.

use std::path::{Path, PathBuf};

use fairchild_core::DeviceRegistry;
use fairchild_osdi::{load_libraries, VaOptions};

mod common;

fn scratch(tag: &str) -> PathBuf {
    // Per-process: parallel `cargo test` runs must not delete each other's.
    let dir = std::env::temp_dir().join(format!("fc_vae2e_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// A stub `openvaf-r`. `counter` gains a line per real compile, so a test can
/// tell a cache hit from a recompile. `artefact` is the `.osdi` it "produces".
fn stub_compiler(dir: &Path, artefact: &Path, counter: &Path) -> PathBuf {
    let path = dir.join("stub-openvaf-r");
    std::fs::write(
        &path,
        format!(
            "#!/bin/sh\n\
             for a in \"$@\"; do case \"$a\" in --version) echo 'stub-openvaf 1.0'; exit 0;; esac; done\n\
             src=''; out=''; expand=0; prev=''\n\
             for a in \"$@\"; do\n\
             \x20 case \"$prev\" in -o) out=\"$a\";; esac\n\
             \x20 case \"$a\" in --print-expansion) expand=1;; *.va) src=\"$a\";; esac\n\
             \x20 prev=\"$a\"\n\
             done\n\
             if [ \"$expand\" = 1 ]; then cat \"$src\"; exit 0; fi\n\
             echo compiled >> '{}'\n\
             cp '{}' \"$out\"\n",
            counter.display(),
            artefact.display(),
        ),
    )
    .unwrap();
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755)).unwrap();
    }
    settle(&path);
    path
}

/// Wait until the OS will run a file we just wrote.
///
/// Four tests here each write an executable and each spawn one, in parallel
/// threads. On Linux exec returns `ETXTBSY` while any process holds a write
/// descriptor to that inode, and a sibling thread that forks mid-write hands its
/// child exactly such a descriptor. One successful exec settles it for good.
/// `compile.rs`'s own `write_stub` carries the long version of this note.
fn settle(path: &Path) {
    let deadline = std::time::Instant::now() + std::time::Duration::from_secs(10);
    loop {
        match std::process::Command::new(path).arg("--version").output() {
            Ok(_) => return,
            Err(e)
                if e.kind() == std::io::ErrorKind::ExecutableFileBusy
                    && std::time::Instant::now() < deadline =>
            {
                std::thread::sleep(std::time::Duration::from_millis(5));
            }
            Err(e) => panic!("stub '{}' will not run: {e}", path.display()),
        }
    }
}

fn compiles(counter: &Path) -> usize {
    std::fs::read_to_string(counter)
        .map(|s| s.lines().count())
        .unwrap_or(0)
}

/// A `.va` source becomes a registered device with no `.osdi` line anywhere —
/// the step this change exists to remove from the user's hands.
#[test]
fn a_va_source_registers_its_device() {
    let Some(artefact) = common::compiled("rc_shunt") else {
        return;
    };
    let dir = scratch("register");
    let counter = dir.join("compiles");
    let stub = stub_compiler(&dir, &artefact, &counter);
    std::fs::write(dir.join("m.va"), "module m(a,b); end\n").unwrap();

    let opts = VaOptions {
        compiler: Some(stub),
        cache_dir: Some(dir.join("cache")),
        ..Default::default()
    };
    let mut registry = DeviceRegistry::new();
    let loaded = load_libraries(&[], &["m.va".to_string()], Some(&dir), &opts, &mut registry)
        .expect("a .va source should compile and load");

    assert_eq!(loaded.len(), 1);
    assert_eq!(
        loaded[0].extension().and_then(|e| e.to_str()),
        Some("osdi"),
        "the artefact is a real .osdi file, portable to any OSDI simulator"
    );
    assert!(
        registry.get("rc_shunt").is_some(),
        "the compiled model must reach the registry, not merely the disk"
    );
    assert_eq!(compiles(&counter), 1);
}

/// Unchanged source, unchanged artefact: the second run must not recompile.
/// Edited source must. Both directions, because only one of them is a bug that
/// shows up as a wrong answer.
#[test]
fn the_cache_hits_on_unchanged_source_and_misses_on_edited() {
    let Some(artefact) = common::compiled("rc_shunt") else {
        return;
    };
    let dir = scratch("cache");
    let counter = dir.join("compiles");
    let stub = stub_compiler(&dir, &artefact, &counter);
    let src = dir.join("m.va");
    std::fs::write(&src, "module m(a,b); end\n").unwrap();

    let opts = VaOptions {
        compiler: Some(stub.clone()),
        cache_dir: Some(dir.join("cache")),
        ..Default::default()
    };
    let compiler = fairchild_osdi::VaCompiler::find(&opts).unwrap();

    let first = fairchild_osdi::compile::compile(&compiler, &src, &opts).unwrap();
    assert_eq!(compiles(&counter), 1);

    // The expansion runs again — it is the cache key — but the artefact it
    // keys is already there, so nothing is recompiled.
    let again = fairchild_osdi::compile::compile(&compiler, &src, &opts).unwrap();
    assert_eq!(again, first, "same source must reuse the same artefact");
    assert_eq!(compiles(&counter), 1, "unchanged source must not recompile");

    std::fs::write(&src, "module m(a,b); // now different\nend\n").unwrap();
    let edited = fairchild_osdi::compile::compile(&compiler, &src, &opts).unwrap();
    assert_ne!(
        edited, first,
        "edited source must produce a different artefact — a stale .osdi is a \
         silently wrong device"
    );
    assert_eq!(compiles(&counter), 2);
}

/// A deck may mix both routes. `.va` sources load first, then `.osdi`
/// artefacts, and the returned order says which is which.
#[test]
fn va_sources_and_osdi_artefacts_load_in_a_defined_order() {
    let Some(artefact) = common::compiled("rc_shunt") else {
        return;
    };
    let dir = scratch("mixed");
    let counter = dir.join("compiles");
    let stub = stub_compiler(&dir, &artefact, &counter);
    std::fs::write(dir.join("m.va"), "module m(a,b); end\n").unwrap();
    // A pre-built artefact, i.e. the offline route the explicit `.osdi` keeps.
    let prebuilt = dir.join("prebuilt.osdi");
    std::fs::copy(&artefact, &prebuilt).unwrap();

    let opts = VaOptions {
        compiler: Some(stub),
        cache_dir: Some(dir.join("cache")),
        ..Default::default()
    };
    let mut registry = DeviceRegistry::new();
    let loaded = load_libraries(
        &["prebuilt.osdi".to_string()],
        &["m.va".to_string()],
        Some(&dir),
        &opts,
        &mut registry,
    )
    .expect("both routes load");

    assert_eq!(loaded.len(), 2);
    assert!(loaded[1].ends_with("prebuilt.osdi"), "{loaded:?}");
    assert!(registry.get("rc_shunt").is_some());
}

/// `--no-va-compile` refuses a `.va` source, and does so whatever the cache
/// holds: a flag whose effect depends on an unseen directory is not offline
/// reproducibility. It must never load the circuit with the device absent.
#[test]
fn no_va_compile_refuses_and_names_the_source() {
    let dir = scratch("refuse");
    let counter = dir.join("compiles");
    let stub = stub_compiler(&dir, Path::new("/dev/null"), &counter);
    std::fs::write(dir.join("m.va"), "module m(a,b); end\n").unwrap();

    let opts = VaOptions {
        compiler: Some(stub),
        cache_dir: Some(dir.join("cache")),
        no_compile: true,
        ..Default::default()
    };
    let mut registry = DeviceRegistry::new();
    let err = load_libraries(&[], &["m.va".to_string()], Some(&dir), &opts, &mut registry)
        .expect_err("told not to compile, so no device — that is an error");
    let msg = err.to_string();
    assert!(msg.contains("m.va"), "{msg}");
    assert!(msg.contains("no-va-compile"), "{msg}");
    assert_eq!(compiles(&counter), 0, "the compiler must not have run");
}
