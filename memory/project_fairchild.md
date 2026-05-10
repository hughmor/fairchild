---
name: fairchild project context
description: Electronic-photonic co-simulation suite: Rust core, ngspice validation, SiPh target
type: project
---

Open-source EO co-simulation suite targeting SiPh. Rust core, PyO3 Python bindings, CLI.

**Why:** No open-source time-domain EO co-simulator exists. Cadence Spectre Photonics is the only working implementation; it is proprietary and extremely expensive.

**Plan document:** `/Users/hugh/Local/src/fairchild/PLAN.md`

**Codebase layout (as of 2026-05-10, Day 2):**
- `crates/fairchild-parser/` — SPICE parser: R, V(DC/PULSE), I, C, L; .op, .tran
- `crates/fairchild-core/` — MNA assembler, DC solver (faer LU), BE transient solver
- `crates/fairchild-cli/` — stub binary
- `tests/golden/` — reference netlists (DC and transient) for ngspice comparison
- `crates/fairchild-core/tests/ngspice_golden.rs` — DC integration test harness
- `crates/fairchild-core/tests/ngspice_tran_golden.rs` — transient test harness (uses .meas)

**Current state (Phase 0, Day 2 complete):**
- DC OP solver for R, V, I — 5 ngspice golden tests passing
- Fixed-step Backward Euler transient solver for R, C, L, V(PULSE), I
- RC and RL step responses validated against ngspice .meas at t=τ,2τ,5τ (1% tolerance)
- Total: 18 tests all green
- Rust 1.95.0, faer 0.24, ngspice 45.2 at /opt/homebrew/bin/ngspice

**Phase 0 falsifiability gate (Weeks 1-6):** CMOS inverter transient via BSIM4 OSDI model.

**Next priorities:**
1. OSDI v0.4 runtime (dlopen .osdi shared libraries from OpenVAF-Reloaded)
2. Build OpenVAF-Reloaded; compile a simple Verilog-A model to .osdi
3. Integrate OSDI model into MNA + transient loop
4. Phase 1 solver hardening: homotopy, variable-step GEAR (after OSDI works)

**Key architectural decisions:**
- Verilog-A compiler: OpenVAF-Reloaded (arpadbuermen fork on Codeberg) → OSDI v0.4 shared libs
- Solver: port VACASK (C++, FOSDEM 2025, arpadbuermen/VACASK on Codeberg) to Rust
- Differentiable: design adjoint sensitivity from the start
- Photonic signal: complex envelope (SVEA), two nodes per optical port, MIT Sorace-Agaskar 2015 convention

**MNA sign conventions (validated against ngspice):**
- Current source `I n+ n-` (SPICE): b[n+] -= dc, b[n-] += dc
- Capacitor BE companion (C/h conductance + I_hist): stamp_current_source(neg, pos, I_hist) → b[pos] += I_hist
- Inductor BE companion (h/L conductance + I_hist): stamp_current_source(pos, neg, I_hist) → b[pos] -= I_hist
- Both companion signs caught by tests on Day 2 — test-first approach continues to pay off

**How to apply:** Use this context when implementing any part of the simulator. The plan doc is the authoritative reference.
