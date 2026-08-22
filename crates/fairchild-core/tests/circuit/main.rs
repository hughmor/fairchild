//! Electrical elements, the solver, and the frontends.
//!
//! One test binary per subject, not per file. Cargo compiles and links one
//! binary per file in `tests/`, so the 14 files in `circuit/` were 14 links
//! against the whole crate at ~25 s each — to run tests that finish in
//! seconds. `cargo test --workspace` was spending 95 % of its 45 minutes in
//! the linker, and none of it on the tests.
//!
//! Cargo takes `tests/<dir>/main.rs` as one target named for the directory,
//! and does not compile the files beside it on their own — so declaring them
//! as modules here collapses that without touching a single test. Name
//! filters still read the same way: `cargo test --test circuit <name>`.
mod controlled_sources;
mod coupled_inductors;
mod differential_pair_dc;
mod global_nets;
mod inductor_is_a_dc_short;
mod model_parameter_diagnostics;
mod no_false_convergence_on_stalled_line_search;
mod passive_parasitics;
mod pcell_subckt;
mod reactive_history_seeding;
mod singular_is_not_nonconvergence;
mod solver_klu;
mod spectre_wrapper;
mod transient_noise;
