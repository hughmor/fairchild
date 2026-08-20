//! `fairchild-osdi` — OSDI (OpenVAF-compiled Verilog-A) runtime.
//!
//! # Status
//!
//! OSDI is the **supported path for device models distributed as Verilog-A** —
//! most importantly the foundry transistor models (BSIM, PSP, HiCUM, …) that
//! the industry ships as compiled `.osdi` shared objects via OpenVAF.
//! fairchild will not hand-write BSIM in Rust; this loader is how those models
//! are consumed. It is exercised in CI by the `osdi-mock` fixture (1 mS in
//! parallel with 1 nF) — keep it: it is the only proof the `dlopen`/FFI/stamp
//! path works end-to-end, for both the resistive and the reactive Jacobian,
//! and removing it would leave this crate with zero coverage.
//!
//! **Optical models work too.** fairchild carries a complex envelope on three
//! ordinary real MNA unknowns per channel (`re`, `im`, `wl`), so a custom
//! `optical_field` / `optical_lambda` discipline is metadata OSDI passes
//! through untouched — no compiler fork, and exact interoperation with the
//! native devices, which use the same wires and units. (OpenVAF 23.5 does
//! panic on a discipline declaring only a potential; the fix is a one-line
//! unused placeholder flow nature, not a fork. See
//! `examples/verilog_a/models/optical.vams`.)
//!
//! **Division of labour.** What a Verilog-A optical model cannot reach is the
//! rest of the abstraction layer: WDM bundle awareness, bidirectional
//! propagation, `crate::delay::DelayLine` group delay, and
//! `PhotonicActiveModel` composition. It is single-channel and forward-only.
//! So photonics needing those stay as native `Device` impls in
//! `fairchild_core::models::photonic`; anything else can go either way. The
//! two paths coexist, which is what a real EPDA tool needs.
//!
//! **The compile is ours to drive.** A deck naming `.va` source — via `.va`,
//! or Spectre's `ahdl_include` — gets it compiled and cached on the way in, so
//! a PDK full of `ahdl_include` lines is not a pile of manual compiles first.
//! `crate::compile` is that path, and explains why it invokes a compiler
//! rather than linking one. The explicit `.osdi` route stays: it is the
//! offline, no-toolchain, reproducible one, and it belongs in CI.
//!
//! Clear-text third-party Verilog-A loads fine; IEEE-1735 / Cadence-encrypted
//! PDKs do not, and cannot — that is upstream of OpenVAF, not a fairchild
//! limitation.
//!
//! Authoring guide and worked examples: `docs/user-guide.md` §14 and
//! `examples/verilog_a/`.

pub mod compile;
pub mod device;
pub mod error;
pub mod ffi;
mod loader;

pub use compile::{load_libraries, VaCompiler, VaOptions};
pub use device::OsdiDevice;
pub use error::OsdiError;
pub use loader::OsdiLibrary;
