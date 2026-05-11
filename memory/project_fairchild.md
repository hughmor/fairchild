---
name: fairchild project context
description: Electronic-photonic co-simulation suite: Rust core, ngspice validation, SiPh target
type: project
---

Open-source EO co-simulation suite targeting SiPh. Rust core, PyO3 Python bindings, CLI.

**Why:** No open-source time-domain EO co-simulator exists. Cadence Spectre Photonics is the only working implementation; it is proprietary and extremely expensive.

**Plan document:** `/Users/hugh/Local/src/fairchild/PLAN.md`

**Current state (Phase 1.5 complete — 2026-05-11):**

Phase 0 complete, Phase 1 mostly complete. Committed work:
- DC NR solver (dc_op_nr): Newton-Raphson + GMIN stepping + source stepping homotopy
- Fixed-step transient: Backward Euler (run_tran) + Trapezoidal Rule (run_tran_tr)
- Nonlinear transient: tran_nr (fixed-step BE+NR), tran_nr_tr (fixed-step TR+NR)
- Variable-step transient: tran_nr_var (BE+LTE predictor-corrector, adaptive h)
- AC small-signal solver: ac_analysis (linearize at DC OP, complex admittance matrix)
- Built-in diode (ShockleyDiode with pnjlim)
- Built-in MOSFET Level 1 (Shichman-Hodges, NMOS+PMOS, body effect via GAMMA/PHI)
- OSDI v0.4 runtime: dlopen, descriptor walk, OsdiDevice Device-trait adapter
- Output: Nutmeg rawfile + CSV for all analysis types (NrResult, TranResult, AcResult)
- DeviceRegistry: OSDI factory + built-in MOSFET card storage
- CLI (fairchild binary): clap-based, -f input, --format csv/nutmeg, --output, AC flags
- README.md, docs/user-guide.md, docs/osdi-reactive-jacobian-findings.md
- examples/ (5 SPICE netlists), examples/compare_ngspice.py, scripts/benchmark.py
- Cleanup: removed MnaSystem legacy wrapper; all compiler warnings fixed

**Open branches:**
- `pulse-breakpoints`: Waveform::next_breakpoint in parser (WIP — not yet wired into tran_nr_with_registry_var)
- `osdi-reactive-jac-investigation`: investigation of OSDI reactive Jacobian bug. Status: root cause identified (drain-drain diagonal missing from reactive stamp). Fix: flat-buffer aliasing via load_jacobian_tran. See docs/osdi-reactive-jacobian-findings.md.

**Remaining Phase 1 items (tracked in PLAN.md Phase 1.5 cleanup section):**
- PULSE breakpoint insertion in variable-step solver (branch: pulse-breakpoints)
- OSDI reactive Jacobian fix (load_jacobian_tran flat-buffer aliasing)
- @cross event detection (zero-crossing step refinement)
- KLU sparse solver via SuiteSparse FFI (optional)

**Next major phase: Phase 2 — Photonic discipline**
Goal: first optical circuit — CW laser → waveguide → coupler → photodetector.

**Key architectural decisions:**
- Verilog-A compiler: OpenVAF-Reloaded (arpadbuermen fork) → OSDI v0.4
- Solver: ported from VACASK (C++, FOSDEM 2025) to Rust
- Photonic signal: complex envelope (SVEA), two nodes per optical port, MIT Sorace-Agaskar 2015
- Differentiable: adjoint sensitivity designed in from the start

**Codebase layout:**
- `crates/fairchild-parser/` — SPICE parser: R, L, C, V(DC/PULSE), I, D, M; .op, .tran, .model, .osdi
- `crates/fairchild-core/` — MNA, Newton-Raphson, transient, AC, device models
  - `src/newton.rs` — dc_op_nr, NrResult (with write_csv/write_nutmeg)
  - `src/tran.rs` — tran_nr, tran_nr_tr, tran_nr_var, TranResult
  - `src/ac.rs` — ac_analysis, AcResult
  - `src/models/diode.rs` — ShockleyDiode
  - `src/models/mosfet1.rs` — Mosfet1 (Level 1 NMOS/PMOS)
  - `src/device_registry.rs` — DeviceRegistry: OSDI factory + mosfet_cards
- `crates/fairchild-osdi/` — OSDI v0.4 runtime: ffi.rs + OsdiDevice
- `crates/fairchild-cli/` — CLI binary (clap)
- `tests/golden/` — reference SPICE netlists for ngspice comparison
- `examples/` — 5 example SPICE netlists
- `scripts/benchmark.py`, `examples/compare_ngspice.py`
- `docs/user-guide.md`, `docs/osdi-reactive-jacobian-findings.md`

**Build commands:**
- `cargo build` / `cargo build --release`
- `cargo test` (all 29 unit tests + integration tests)
- `./target/debug/fairchild -f examples/rc_step.sp`
- `python examples/compare_ngspice.py --release` (requires ngspice + matplotlib)
- `python scripts/benchmark.py` (requires ngspice)

**OSDI reactive Jacobian fix (NOT YET DONE):**
- Root cause: load_spice_rhs_tran embeds G_react*V in b-vector, but write_jacobian_array_react
  misses the (drain,drain) diagonal entry → unmatched stamp → voltage drift in transient.
- Fix: use load_jacobian_tran (aliasing path) with a per-device flat buffer; zero buffer,
  call load_jacobian_tran(inst, model, alpha), scatter-add into mat.a.
- Affects: CMOS inverter switching transient with OSDI MOSFET models.

**How to apply:** Use this context when implementing any part of the simulator. The plan doc is the authoritative reference.
