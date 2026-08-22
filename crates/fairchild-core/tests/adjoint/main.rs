//! Parameter sensitivity: DC, AC and transient adjoints.
//!
//! One test binary per subject, not per file. Cargo compiles and links one
//! binary per file in `tests/`, so the 4 files in `adjoint/` were 4 links
//! against the whole crate at ~25 s each — to run tests that finish in
//! seconds. `cargo test --workspace` was spending 95 % of its 45 minutes in
//! the linker, and none of it on the tests.
//!
//! Cargo takes `tests/<dir>/main.rs` as one target named for the directory,
//! and does not compile the files beside it on their own — so declaring them
//! as modules here collapses that without touching a single test. Name
//! filters still read the same way: `cargo test --test adjoint <name>`.
mod adjoint_ac;
mod adjoint_dc;
mod adjoint_jacobian;
mod adjoint_tran;
