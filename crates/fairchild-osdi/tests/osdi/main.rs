//! The OSDI ABI: loading, parameters, stamping, limiting, diagnostics.
//!
//! One test binary per subject, not per file. Cargo compiles and links one
//! binary per file in `tests/`, so these 14 files were 14 links against the
//! whole crate at ~25 s each — to run tests that finish in seconds.
//!
//! Cargo takes `tests/<dir>/main.rs` as one target named for the directory, and
//! does not compile the files beside it on their own. Name filters read the
//! same way as before: `cargo test --test osdi <name>`.

// `tests/common/` is shared by both groups, so it stays where it is and each
// group reaches it from the crate root rather than keeping a copy.
#[path = "../common/mod.rs"]
mod common;

mod bsim4_acceptance;
mod cmos_inverter;
mod load_compiled;
mod load_real_osdi;
mod osdi_abi_contract;
mod osdi_abstime;
mod osdi_dc_op;
mod osdi_device;
mod osdi_diag;
mod osdi_limiting;
mod osdi_model_card;
mod osdi_reactive;
mod osdi_tran;
mod setup_attribution;
mod thermal_discipline;
mod va_compile;
