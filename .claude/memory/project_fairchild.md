# fairchild — Project Context

## Goal

First open-source **time-domain electro-optic co-simulator**: a Rust SPICE engine that co-simulates electronic and photonic circuits in the same Newton-Raphson loop. Differentiable (adjoint-method gradients) and Python-accessible via PyO3.

**Why it matters**: All open-source photonic simulators (SAX, Simphony, Photontorch) are frequency-domain S-matrix only. Cadence Spectre Photonics is the only time-domain EO co-simulator — it costs ~$100k/seat/yr.

---

## Phase Status

| Phase | Status | Summary |
|-------|--------|---------|
| 0 — Foundation | ✅ done | Cargo workspace, MNA, Newton-Raphson, SPICE parser, OSDI runtime |
| 1 — Solver hardening | ✅ done | Homotopy, variable-step BE+LTE, Trapezoidal Rule, AC, Nutmeg output, Level 1 MOSFET, Shockley diode |
| 1.5 — Consolidation | ✅ done | CLI, docs, examples, ngspice validation, OSDI reactive Jacobian fix, PWL, `.ac` directive |
| 2 — Photonic discipline | 🔄 in progress | 3/4 ring resonator tests pass; sweep assertions need 2-line fix — see `phase_2.md` |
| 3 — Python bindings | 📋 planned | See `phase_3.md` |
| 4 — Differentiable sim | 📋 planned | See `phase_4.md` |
| 5+ — PDK, model zoo | 📋 planned | See PLAN.md Parts 3/5-7 |

---

## Crate Structure

```
crates/
  fairchild-core/     # DAE solver, MNA, Newton-Raphson, transient (BE/TR/var-step), AC
  fairchild-parser/   # SPICE parser: R,L,C,V,I,M,D; DC/PULSE/PWL; .op/.tran/.ac/.model
  fairchild-cli/      # Binary: -f netlist.sp, --format nutmeg, --output, --ac-* flags
  fairchild-osdi/     # OSDI v0.4 runtime: dlopen, OsdiDevice wraps *const OsdiDescriptor
```

---

## Key Architectural Decisions

**Device trait** — all devices (built-in Rust + OSDI-loaded) implement the same trait:
`eval()` → `load_residual()` → `load_jacobian()`, and the transient variants
`load_residual_tran(b, alpha)` / `load_jacobian_tran(mat, alpha)` / `commit_timestep(x)`.

**MNA matrix** — `Vec<Vec<f64>>` (dense; KLU deferred to later). `a[row][col]` = conductance, `b[row]` = current source. Node 0 = ground, eliminated. Voltage source branches appended after node rows.

**OSDI Jacobian path** — use the **copy path** (`write_jacobian_array_resist` / `write_jacobian_array_react`), NOT the aliasing path (`load_jacobian_resist`). The aliasing path crashes (SIGSEGV) with OpenVAF-compiled models. Copy path is sufficient and correct.

**OSDI reactive Jacobian** — `write_jacobian_array_react` writes `n_react` values in traversal order of `jacobian_entries[0..n_total]` where `react_ptr_off != u32::MAX`. Iterate ALL entries, skip those with `react_ptr_off == MAX`, consume `jac_buf` in that order. Use `alpha = 1/h` (not 1.0). Stamp: `mat.a[mr][mc] += alpha * jac_buf[react_idx]`.

**Reactive history (`x_tprev`)** — `OsdiDevice` holds `x_tprev` (solution at last accepted timestep). `commit_timestep(x)` snapshots it. `load_residual_tran` passes `x_tprev[1..]` as `prev_solve` to `load_spice_rhs_tran`. The `[0]` guard element (=0.0) handles OSDI's ground sentinel (UINT32_MAX → -1 index).

**Optical discipline (Phase 2)** — SVEA two-wire convention: `(o_re, o_im)` complex envelope per port, maps to two MNA nodes. Carrier frequency `ω₀` is a `SimContext` parameter, not a nodal variable.

**Norton-equivalent optical models** — VA models MUST use flow contributions (Oamp) only, never potential contributions (Ophase<+). Potential contributions add internal OSDI branch nodes that the runtime doesn't allocate (confirmed: Ophase<+ on a port adds internal nodes, causing `num_nodes > num_terminals`). Use Norton equivalent: `Oamp(p) <+ G*(target - Ophase(p))` with G=1e6.

**XOsdi element** — SPICE `X<name> net0 ... netN model_name [key=val ...]` syntax for arbitrary-port OSDI instances. Parsed in fairchild-parser. Discipline checking via `.optical net1 net2...` directive + `check_disciplines()`. XOsdi devices instantiated in `build_devices()` (newton.rs) by model name lookup.

---

## Build

```bash
cargo build --release
cargo test
./target/release/fairchild -f examples/rc_step.sp
```

OpenVAF-Reloaded: `/Users/hugh/Local/src/OpenVAF-Reloaded/target/release/openvaf-r`
Runtime env: `DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib`
