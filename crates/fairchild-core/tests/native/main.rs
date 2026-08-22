//! Native photonic devices.
//!
//! One test binary per subject, not per file. Cargo compiles and links one
//! binary per file in `tests/`, so the 20 files in `native/` were 20 links
//! against the whole crate at ~25 s each — to run tests that finish in
//! seconds. `cargo test --workspace` was spending 95 % of its 45 minutes in
//! the linker, and none of it on the tests.
//!
//! Cargo takes `tests/<dir>/main.rs` as one target named for the directory,
//! and does not compile the files beside it on their own — so declaring them
//! as modules here collapses that without touching a single test. Name
//! filters still read the same way: `cargo test --test native <name>`.
mod native_awgr;
mod native_circulator;
mod native_driven_laser;
mod native_facet;
mod native_grating_coupler;
mod native_mrr_add_port;
mod native_mrr_modulator;
mod native_mzm;
mod native_optical_2x2;
mod native_pd_r_series;
mod native_photonic_expr;
mod native_photonic_level;
mod native_pn_ps_cap;
mod native_pn_ps_full;
mod native_pn_th_ps;
mod native_splitter_asymmetric;
mod native_thermal_ps_rc;
mod native_wdm_mrr;
mod native_wdm_mux;
mod optical_noise;
