---
name: fairchild project context
description: Electronic-photonic co-simulation suite: Rust core, ngspice validation, SiPh target
type: project
---

Open-source EO co-simulation suite targeting SiPh. Rust core, PyO3 Python bindings, CLI.

**Why:** No open-source time-domain EO co-simulator exists. Cadence Spectre Photonics is the only working implementation; it is proprietary and extremely expensive.

**Plan document:** `/Users/hugh/Local/src/fairchild/PLAN.md`

**Codebase layout (as of 2026-05-10):**
- `crates/fairchild-parser/` — SPICE netlist parser (R, V, I, .op)
- `crates/fairchild-core/` — MNA assembler + DC solver (faer 0.24 dense LU)
- `crates/fairchild-cli/` — stub binary
- `tests/golden/` — reference netlists for ngspice comparison
- `crates/fairchild-core/tests/ngspice_golden.rs` — integration test harness

**Current state (Phase 0, Day 1 complete):**
- DC operating-point solver for resistive circuits (R, V, I) working
- 5 ngspice golden comparison tests passing (voltage divider, current divider, Wheatstone bridge, ladder, multi-source)
- Rust 1.95.0, faer 0.24, ngspice 45.2 at /opt/homebrew/bin/ngspice

**Phase 0 falsifiability gate (Weeks 1-6):** CMOS inverter transient via BSIM4 OSDI model.

**Next priorities:**
1. GEAR order-1 transient integrator (BDF-1 = Backward Euler)
2. Capacitor and inductor MNA stamps for transient
3. ngspice golden tests for RC/RL/RLC transient responses
4. OSDI runtime for loading OpenVAF-Reloaded compiled models

**Key architectural decisions:**
- Verilog-A compiler: OpenVAF-Reloaded (arpadbuermen fork on Codeberg) → OSDI v0.4 shared libs
- Solver: port VACASK (C++, FOSDEM 2025, arpadbuermen/VACASK on Codeberg) to Rust
- Differentiable: design adjoint sensitivity from the start
- Photonic signal: complex envelope (SVEA), two nodes per optical port, MIT Sorace-Agaskar 2015 convention

**MNA sign convention (validated 2026-05-10 against ngspice):**
- Current source `I n+ n-`: SPICE = current flows from n+ through source to n-
  → b[n+] -= dc, b[n-] += dc
- ngspice tests immediately caught a sign bug — validates test-first approach

**How to apply:** Use this context when implementing any part of the simulator. The plan doc is the authoritative reference.
