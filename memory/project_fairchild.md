---
name: fairchild project context
description: Electronic-photonic co-simulation suite: Rust core, ngspice validation, SiPh target
type: project
---

Open-source EO co-simulation suite targeting SiPh. Rust core, PyO3 Python bindings, CLI.

**Why:** No open-source time-domain EO co-simulator exists. Cadence Spectre Photonics is the only working implementation; it is proprietary and extremely expensive.

**Plan document:** `/Users/hugh/Local/src/fairchild/PLAN.md`

**Codebase layout (as of 2026-05-10, Day 3):**
- `crates/fairchild-parser/` — SPICE parser: R, V(DC/PULSE), I, C, L; .op, .tran
- `crates/fairchild-core/` — MNA assembler, DC solver (faer LU), BE transient solver
- `crates/fairchild-osdi/` — OSDI v0.4 runtime: #[repr(C)] FFI structs + dlopen loader
- `crates/osdi-mock/` — cdylib test fixture; exports one "test_diode" descriptor
- `crates/fairchild-cli/` — stub binary
- `tests/golden/` — reference netlists (DC and transient) for ngspice comparison
- `crates/fairchild-core/tests/ngspice_golden.rs` — DC integration test harness
- `crates/fairchild-core/tests/ngspice_tran_golden.rs` — transient test harness (uses .meas)
- `crates/fairchild-osdi/tests/load_mock.rs` — OSDI registry-walk integration test

**Current state (Phase 0, Day 3 complete):**
- DC OP solver for R, V, I — 5 ngspice golden tests passing
- Fixed-step Backward Euler transient solver for R, C, L, V(PULSE), I
- RC and RL step responses validated against ngspice .meas at t=τ,2τ,5τ (1% tolerance)
- OSDI v0.4 runtime: dlopen + symbol resolution + registry walk — 1 integration test passing
- Total: 19 tests all green
- Rust 1.95.0, faer 0.24, libc 0.2, ngspice 45.2 at /opt/homebrew/bin/ngspice

**Phase 0 falsifiability gate (Weeks 1-6):** CMOS inverter transient via BSIM4 OSDI model.

**Next priorities:**
1. Build OpenVAF-Reloaded on macOS to produce a real v0.4 .osdi file (a simple VA model)
2. Load and call setup_model/setup_instance/eval on a real .osdi model
3. Integrate OSDI device into MNA stamp_netlist → transient loop (Newton-Raphson)
4. Validate OSDI diode model against ngspice golden (forward-bias I-V and transient)
5. Phase 1 solver hardening: homotopy, variable-step GEAR (after OSDI works end-to-end)

**Key architectural decisions:**
- Verilog-A compiler: OpenVAF-Reloaded (arpadbuermen fork, GitHub mob branch) → OSDI v0.4
- Pre-built .osdi files from VA-Models v2.2 are v0.3 ELF (Linux x86-64, useless on macOS)
- ngspice 45.2 requires OSDI >= v0.4; must build OpenVAF-Reloaded to get v0.4 on macOS
- Solver: port VACASK (C++, FOSDEM 2025, arpadbuermen/VACASK on Codeberg) to Rust
- Differentiable: design adjoint sensitivity from the start
- Photonic signal: complex envelope (SVEA), two nodes per optical port, MIT Sorace-Agaskar 2015

**OSDI v0.4 runtime notes:**
- OsdiDescriptor size = 312 bytes (verified by compile-time assertion in ffi.rs)
- Exported symbols: OSDI_VERSION_MAJOR, OSDI_VERSION_MINOR, OSDI_NUM_DESCRIPTORS,
  OSDI_DESCRIPTOR_SIZE, OSDI_DESCRIPTORS (array), osdi_log (fn ptr filled by simulator)
- Descriptor array strided by OSDI_DESCRIPTOR_SIZE bytes (not sizeof) for forward-compat
- Do NOT use libloading on macOS for data symbols: returns false IncompatibleSize errors
  due to macOS Mach-O n_size check. Raw libc::dlopen/dlsym works correctly.
- OSDI_DESCRIPTORS is in __DATA_CONST,__const on macOS (read-only after init)
- osdi_log is in __DATA,__common (BSS, writable — simulator fills it in at load time)

**MNA sign conventions (validated against ngspice):**
- Current source `I n+ n-` (SPICE): b[n+] -= dc, b[n-] += dc
- Capacitor BE companion (C/h conductance + I_hist): stamp_current_source(neg, pos, I_hist) → b[pos] += I_hist
- Inductor BE companion (h/L conductance + I_hist): stamp_current_source(pos, neg, I_hist) → b[pos] -= I_hist
- Both companion signs caught by tests on Day 2 — test-first approach continues to pay off

**How to apply:** Use this context when implementing any part of the simulator. The plan doc is the authoritative reference.
