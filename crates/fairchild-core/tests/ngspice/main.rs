//! Golden comparisons against an external simulator.
//!
//! One test binary per subject, not per file. Cargo compiles and links one
//! binary per file in `tests/`, so the 12 files in `ngspice/` were 12 links
//! against the whole crate at ~25 s each — to run tests that finish in
//! seconds. `cargo test --workspace` was spending 95 % of its 45 minutes in
//! the linker, and none of it on the tests.
//!
//! Cargo takes `tests/<dir>/main.rs` as one target named for the directory,
//! and does not compile the files beside it on their own — so declaring them
//! as modules here collapses that without touching a single test. Name
//! filters still read the same way: `cargo test --test ngspice <name>`.
mod ac_filter_golden;
mod ngspice_bjt_golden;
mod ngspice_coupled_inductors_golden;
mod ngspice_diode_breakdown_golden;
mod ngspice_diode_golden;
mod ngspice_diode_tran_golden;
mod ngspice_gmin_golden;
mod ngspice_golden;
mod ngspice_method_golden;
mod ngspice_mosfet_golden;
mod ngspice_noise_golden;
mod ngspice_switch_golden;
mod ngspice_temperature_golden;
mod ngspice_tf_pz_golden;
mod ngspice_tline_golden;
mod ngspice_tran_golden;
mod ring_oscillator_golden;
