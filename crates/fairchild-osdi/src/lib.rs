//! `fairchild-osdi` — OSDI (OpenVAF-compiled Verilog-A) compatibility layer.
//!
//! # Status (Phase B5+)
//!
//! As of the photonic refactor in Phase B, **OSDI is no longer the
//! recommended path** for new photonic devices.  Native Rust implementations
//! (in `fairchild_core::models::photonic`) are:
//!
//! - faster (no FFI / `dlopen` per simulation),
//! - simpler to author (a Rust struct + `Device` impl is shorter than the
//!   Verilog-A source the equivalent OpenVAF model needs),
//! - decoupled from the OpenVAF 23.5 codegen bug that crashes on a
//!   potential-only optical discipline (worked around in the photonic
//!   `.vams` file but still a maintenance hazard), and
//! - decoupled from IEEE-1735 / Cadence-encrypted PDK models, which OpenVAF
//!   cannot compile at all.
//!
//! This crate remains supported as a **compatibility shim** for two use
//! cases:
//!
//! 1. Loading clear-text third-party Verilog-A models that the user already
//!    has compiled and doesn't want to port to Rust (analog mixed-signal
//!    BSIM/PSP/etc.).
//! 2. Loading the legacy fairchild photonic models that haven't yet been
//!    migrated to the B1 discipline scheme (the MRR/MZI/PN-PS family).
//!
//! For new model authoring, prefer `Device` impls registered in
//! `fairchild_core::DeviceRegistry`.  See `models/photonic.rs` for examples.

pub mod device;
pub mod error;
pub mod ffi;
mod loader;

pub use device::OsdiDevice;
pub use error::OsdiError;
pub use loader::OsdiLibrary;
