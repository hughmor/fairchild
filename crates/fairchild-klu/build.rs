//! Locate the SuiteSparse KLU library on the host and emit the link
//! directive `cargo:rustc-link-lib=klu`.
//!
//! Probing order:
//!   1. `KLU_LIB_DIR` env var (explicit override).
//!   2. `pkg-config --libs klu` (the canonical path on Linux + macOS Homebrew).
//!   3. Homebrew default at `/opt/homebrew/opt/suite-sparse/lib`.
//!   4. Linux fallbacks: `/usr/lib`, `/usr/local/lib`.
//!
//! A clear panic with install instructions is emitted if KLU cannot be
//! located — the user opted into the `klu` cargo feature, so a missing
//! system library is a build error, not a silent fallback.

use std::env;
use std::process::Command;

fn main() {
    println!("cargo:rerun-if-env-changed=KLU_LIB_DIR");
    println!("cargo:rerun-if-env-changed=PKG_CONFIG_PATH");

    // 1. Explicit override via env.
    if let Ok(dir) = env::var("KLU_LIB_DIR") {
        println!("cargo:rustc-link-search=native={dir}");
        println!("cargo:rustc-link-lib=klu");
        // KLU depends on AMD, COLAMD, BTF, and SuiteSparse_config — pkg-config
        // pulls them transitively but for the manual override we link
        // them explicitly so static + dynamic both work.
        println!("cargo:rustc-link-lib=amd");
        println!("cargo:rustc-link-lib=colamd");
        println!("cargo:rustc-link-lib=btf");
        println!("cargo:rustc-link-lib=suitesparseconfig");
        return;
    }

    // 2. pkg-config — the recommended path on macOS (Homebrew installs
    //    klu.pc) and on Linux distributions.
    if let Ok(out) = Command::new("pkg-config").args(["--libs", "klu"]).output() {
        if out.status.success() {
            let s = String::from_utf8_lossy(&out.stdout);
            let mut linked_klu = false;
            for tok in s.split_whitespace() {
                if let Some(path) = tok.strip_prefix("-L") {
                    println!("cargo:rustc-link-search=native={path}");
                } else if let Some(lib) = tok.strip_prefix("-l") {
                    println!("cargo:rustc-link-lib={lib}");
                    if lib == "klu" { linked_klu = true; }
                }
            }
            if linked_klu {
                return;
            }
        }
    }

    // 3. Homebrew default on Apple Silicon (suite-sparse installs here even
    //    when its pkg-config file is missing from PKG_CONFIG_PATH).
    let brew_lib = "/opt/homebrew/opt/suite-sparse/lib";
    if std::path::Path::new(brew_lib).exists() {
        println!("cargo:rustc-link-search=native={brew_lib}");
        println!("cargo:rustc-link-lib=klu");
        return;
    }

    // 4. Linux fallback.
    for dir in ["/usr/lib/x86_64-linux-gnu", "/usr/lib", "/usr/local/lib"] {
        if std::path::Path::new(&format!("{dir}/libklu.so")).exists()
            || std::path::Path::new(&format!("{dir}/libklu.so.0")).exists()
        {
            println!("cargo:rustc-link-search=native={dir}");
            println!("cargo:rustc-link-lib=klu");
            return;
        }
    }

    panic!(
        "Could not locate SuiteSparse KLU on this system.\n\n\
         The `klu` feature requires a system install of SuiteSparse:\n  \
         macOS:  brew install suite-sparse\n  \
         Debian: sudo apt install libsuitesparse-dev\n  \
         Fedora: sudo dnf install suitesparse-devel\n\n\
         If KLU lives in a non-standard location, set KLU_LIB_DIR=<path>.\n"
    );
}
