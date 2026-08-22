//! WDM bundles, wavelength labels, and bidirectional propagation.
//!
//! One test binary per subject, not per file. Cargo compiles and links one
//! binary per file in `tests/`, so the 8 files in `bundles/` were 8 links
//! against the whole crate at ~25 s each — to run tests that finish in
//! seconds. `cargo test --workspace` was spending 95 % of its 45 minutes in
//! the linker, and none of it on the tests.
//!
//! Cargo takes `tests/<dir>/main.rs` as one target named for the directory,
//! and does not compile the files beside it on their own — so declaring them
//! as modules here collapses that without touching a single test. Name
//! filters still read the same way: `cargo test --test bundles <name>`.
mod bidirectional_composition;
mod bidirectional_endtoend;
mod bidirectional_option;
mod bundle_arity_registry;
mod dark_channel_transient;
mod lambda_is_a_label;
mod lambda_is_not_an_unknown;
mod lambda_resolution;
