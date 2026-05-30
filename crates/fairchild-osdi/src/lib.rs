//! `fairchild-osdi` — OSDI (OpenVAF-compiled Verilog-A) runtime.
//!
//! # Status (2026-05-30: un-deprecated for electrical models)
//!
//! OSDI is the **supported path for electrical device models distributed as
//! Verilog-A** — most importantly the foundry transistor models (BSIM, PSP,
//! HiCUM, …) that the industry ships as compiled `.osdi` shared objects via
//! OpenVAF. fairchild will not hand-write BSIM in Rust; the OSDI loader is how
//! those models are consumed. The loader is exercised in CI by the `osdi-mock`
//! fixture (a tiny conductance model) — keep it: it is the only proof the
//! `dlopen`/FFI/stamp path works end-to-end, and removing it would leave this
//! crate with zero coverage.
//!
//! **Division of labour (the deliberate architecture):**
//!
//! - *Electrical / mixed-signal* (transistors, foundry PDKs): **OSDI**. The
//!   Verilog-A → OpenVAF → `.osdi` flow is the right delegation; it reaches
//!   models we could never maintain by hand.
//! - *Photonic* (waveguides, couplers, phase shifters, detectors, …): **native
//!   Rust** `Device` impls in `fairchild_core::models::photonic`. OSDI/Verilog-A
//!   cannot express fairchild's optical abstractions (complex-envelope bundle
//!   ports, centre wavelength, bidirectional wires), and OpenVAF 23.5 has a
//!   codegen bug on potential-only optical disciplines. So photonics stay
//!   native; the two paths coexist, which is exactly what a real EPDA tool
//!   needs.
//!
//! It also loads clear-text third-party Verilog-A models the user has already
//! compiled (it cannot compile IEEE-1735 / Cadence-encrypted PDKs — neither can
//! OpenVAF). Authoring a *new photonic* model is still a native Rust `Device`
//! impl; see `models/photonic.rs`.

pub mod device;
pub mod error;
pub mod ffi;
mod loader;

pub use device::OsdiDevice;
pub use error::OsdiError;
pub use loader::OsdiLibrary;
