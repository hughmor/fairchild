//! Optical models through OSDI, including the bundle dialect.
//!
//! One test binary per subject, not per file. Cargo compiles and links one
//! binary per file in `tests/`, so these 6 files were 6 links against the
//! whole crate at ~25 s each — to run tests that finish in seconds.
//!
//! Cargo takes `tests/<dir>/main.rs` as one target named for the directory, and
//! does not compile the files beside it on their own. Name filters read the
//! same way as before: `cargo test --test photonic <name>`.

// `tests/common/` is shared by both groups, so it stays where it is and each
// group reaches it from the crate root rather than keeping a copy.
#[path = "../common/mod.rs"]
mod common;

mod bundle_dialect_e2e;
mod load_optical_osdi;
mod mrm_wdm_example;
mod optical_wdm_bundle;
mod photonic_models;
mod ring_resonator;
