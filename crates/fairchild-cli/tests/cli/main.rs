//! Tests that have to drive the real binary.
//!
//! Everything here checks something that only exists once a process runs: a
//! warning that leaves library code through `warn_user!` to stderr, a `--quiet`
//! switch that has to reach the whole library, a `--probe` name that must be
//! refused rather than dropped. No in-process test can see any of it.
//!
//! One binary, not one per file. Cargo compiles and links each `tests/*.rs` as
//! its own binary against the whole crate, so a file per subject is a link per
//! subject — the reason `crates/*/tests/<subject>/main.rs` is the pattern
//! everywhere else in this tree. These three were two separate binaries before
//! `transient_noise_flatness` would have made a third.

mod dropped_parameters;
mod quiet;
mod transient_noise_flatness;
