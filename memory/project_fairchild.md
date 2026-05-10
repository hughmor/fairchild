---
name: fairchild project context
description: Electronic-photonic co-simulation suite: Rust core, ngspice validation, SiPh target
type: project
---

Open-source EO co-simulation suite targeting SiPh. Rust core, PyO3 Python bindings, CLI.

**Why:** No open-source time-domain EO co-simulator exists. Cadence Spectre Photonics is the only working implementation; it is proprietary and extremely expensive.

**Plan document:** `/Users/hugh/Local/src/fairchild/PLAN.md`

**Codebase layout (as of 2026-05-10, Day 5):**
- `crates/fairchild-parser/` — SPICE parser: R, V(DC/PULSE), I, C, L, D; .op, .tran, .model
- `crates/fairchild-core/` — MNA assembler, DC solver (faer LU), BE transient, NR solver
  - `src/device.rs` — Device trait, SimContext, NodeId, EvalFlags
  - `src/models/diode.rs` — ShockleyDiode with pnjlim
  - `src/newton.rs` — dc_op_nr: Newton-Raphson loop with pnjlim + global VMAX damping
- `crates/fairchild-osdi/` — OSDI v0.4 runtime: #[repr(C)] FFI structs + dlopen loader
  - `src/device.rs` — OsdiDevice: Device trait adapter for OSDI descriptors
- `crates/osdi-mock/` — cdylib test fixture; exports "test_conductance" with real fn ptrs
- `crates/fairchild-cli/` — stub binary
- `tests/golden/` — reference netlists (DC, transient, diode) for ngspice comparison
- `crates/fairchild-core/tests/ngspice_golden.rs` — DC integration test harness
- `crates/fairchild-core/tests/ngspice_tran_golden.rs` — transient test harness (uses .meas)
- `crates/fairchild-core/tests/ngspice_diode_golden.rs` — nonlinear DC diode test harness
- `crates/fairchild-osdi/tests/load_mock.rs` — OSDI registry-walk integration test
- `crates/fairchild-osdi/tests/osdi_device.rs` — OsdiDevice MNA stamp integration tests

**Current state (Phase 0, Day 6 complete):**
- DC OP solver for R, V, I — 5 ngspice golden tests passing
- Fixed-step Backward Euler transient solver for R, C, L, V(PULSE), I (run_tran)
- RC and RL step responses validated against ngspice .meas at t=τ,2τ,5τ (1% tolerance)
- OSDI v0.4 runtime: dlopen + symbol resolution + registry walk — 1 integration test passing
- Newton-Raphson DC solver (dc_op_nr) with ShockleyDiode device:
  - pnjlim voltage limiting + global VMAX=0.5V backup damping
  - 2 ngspice golden tests: current-source bias and R-D series (both within 0.1% of ngspice)
- OsdiDevice wrapper: Device trait adapter for OSDI v0.4 descriptors
  - setup_model/setup_instance/eval/load_residual/load_jacobian all implemented
  - Uses copy-based Jacobian path (write_jacobian_array_resist) — no pointer aliasing
  - Handles ground terminals (NodeId = None) correctly
- osdi-mock extended: "test_conductance" model with real function pointers (1 mS conductance)
- tran_nr: nonlinear transient solver (BE companion + NR per step)
  - Starts from dc_op_nr operating point; capacitors seeded from DC voltages
  - Same pnjlim + VMAX damping as dc_op_nr; build_devices() shared helper
  - Validated: pure RC matches run_tran; R-D steady state matches dc_op_nr
  - ngspice golden: half-wave rectifier V(cap) at t=550µs within 2%
- Total: 35 tests all green
- Rust 1.95.0, faer 0.24, libc 0.2, ngspice 45.2 at /opt/homebrew/bin/ngspice

**Phase 0 falsifiability gate (Weeks 1-6):** CMOS inverter transient via BSIM4 OSDI model.

**Next priorities:**
1. Build OpenVAF-Reloaded on macOS (arpadbuermen fork, mob branch) to produce a real v0.4 .osdi
2. Validate OsdiDevice + real .osdi diode model against ngspice golden (forward-bias I-V)
3. Integrate OsdiDevice into dc_op_nr/tran_nr (load OSDI devices from netlist alongside built-ins)
   — requires parser extension: `.osdi path` directive + element-to-model name resolution
4. Jacobian aliasing optimisation: refactor MnaMatrix to contiguous Vec<f64>, then point OSDI
   jacobian_ptr_resist array directly into MnaMatrix (eliminates copy on critical NR path)

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

**Device trait and NR solver notes:**
- Device trait shape matches PLAN.md §2.5: setup_model(&SimContext), setup_instance(&[NodeId]),
  eval(&[f64], EvalFlags, &SimContext), load_residual(&mut [f64]), load_jacobian(&mut MnaMatrix)
- NodeId = Option<usize>; None = connected to ground (excluded from MNA matrix)
- pnjlim formula: vcrit = Vt * ln(Vt / (√2 * Is)); compress step logarithmically above vcrit
- GMIN = 1e-12 S added to every diode conductance stamp (matches ngspice default)
- Global voltage step damping VMAX=0.5V as backup (pnjlim handles the critical exponential range)
- Convergence: |Δx| < VNTOL(1μV) + RELTOL(0.1%) * |x|; MAX_ITER=150

**MNA sign conventions (validated against ngspice):**
- Current source `I n+ n-` (SPICE): b[n+] -= dc, b[n-] += dc
- Capacitor BE companion (C/h conductance + I_hist): stamp_current_source(neg, pos, I_hist) → b[pos] += I_hist
- Inductor BE companion (h/L conductance + I_hist): stamp_current_source(pos, neg, I_hist) → b[pos] -= I_hist
- Diode companion: gd conductance (anode↔cathode) + Jeq = Id - gd*Vd current from anode to cathode
  → b[anode] -= Jeq, b[cathode] += Jeq

**How to apply:** Use this context when implementing any part of the simulator. The plan doc is the authoritative reference.
