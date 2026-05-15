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
| 2 — Photonic discipline | ✅ done | Ring resonator sweep validated against CMT (Δλ=0.02 nm); example + plot script in repo |
| 2.5 — Photonic model library | ✅ done | 20 models: PN PS L1/L2, thermo PS L1/L2, PD L2, MRR mod L1/L2/L3, heater MRR L1/L2, add-drop MRR (PN L1/L2/L3, heater L1/L2), MZI (PN L1/L2, thermo L1/L2) |
| 3 — Python bindings | ✅ done | `fairchild-py`: Circuit/SimResult/WaveformSource, op/tran/AC, numpy arrays, maturin |
| 3.5 — Parser improvements | ✅ done | `.subckt`/`.ends` two-pass flattening, `.param` global params, `{param}` substitution, `.include` pre-processing, unsupported-directive errors; branch `feature/subckt-support` |
| 3.6 — Compound photonic subckts | ✅ done | All-pass MRR (PN+thermo) and add-drop MRR (PN) compound subckts verified; MZI subckts created but have known two-DC NR convergence bug (use monolithic .osdi models) |
| 3.7 — WDM Tier 2 + KiCAD integration | ✅ done | Bus vector expansion `net[M..N]`→`net_M..net_N`; `.optical_bus N re im wl`; `scripts/kicad_to_fairchild.py` post-processor; `kicad_integration.md` setup guide (not committed); 2-channel smoketest verified ✓ |
| **A — Foundation rewrite (Tier 0)** | ✅ done 2026-05-15 | All 10 SoTU Tier-0 items shipped on `feature/subckt-support`. SimOptions struct with `.options` parse + CLI/Python kwargs; `.dc` sweep as first-class analysis; SIN/EXP/SFFM/AM source waveforms; `.ic` / `.nodeset` / UIC; `.lib 'file' section` + `.endl`; B-element behavioral sources (full expression grammar in `fairchild-parser::expr`); `.measure` post-processing (FIND/MAX/AVG/RMS/PP/INTEG/DERIV/TRIG); floating-node connectivity check before LU. **Three-surface rule established**: every new feature must land on parser + CLI + Python. See `sotu.md` (uncommitted, repo root) for the State-of-the-Union write-up that motivated this phase. |
| 4 — Differentiable sim | 📋 planned | See `phase_4.md` |
| 5+ — PDK, model zoo | 📋 planned | See PLAN.md Parts 3/5-7 |

---

## Crate Structure

```
crates/
  fairchild-core/     # DAE solver, MNA, Newton-Raphson, transient (BE/TR/var-step), AC
  fairchild-parser/   # SPICE parser: R,L,C,V,I,M,D; DC/PULSE/PWL; .op/.tran/.ac/.model
  fairchild-cli/      # Binary: -f netlist.sp, --format nutmeg, --probe, --param, --check
  fairchild-osdi/     # OSDI v0.4 runtime: dlopen, OsdiDevice wraps *const OsdiDescriptor
  fairchild-py/       # PyO3 Python package: Circuit/SimResult/WaveformSource, maturin build
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

**Optical discipline (Phase 2)** — 3-wire convention: `(o_re, o_im, o_lambda)` per optical port. All three use the single `optical` discipline (`Ophase`/`Oamp`). `o_lambda` stores wavelength in µm. All lambda terminals on the same optical bus share one `wl` net. Declared via `.optical` SPICE directive (no separate `.optical_lambda`). CRITICAL: `Optical_Phase` nature must have `units = "rad"` (different from `Optical_Amplitude`'s `"sqrt(W)"`) — same units triggers OpenVAF flow-node generation and NR divergence.

**Norton-equivalent optical models** — VA models MUST use flow contributions (Oamp) only, never potential contributions (Ophase<+). Potential contributions add internal OSDI branch nodes that the runtime doesn't allocate (confirmed: Ophase<+ on a port adds internal nodes, causing `num_nodes > num_terminals`). Use Norton equivalent: `Oamp(p) <+ G*(target - Ophase(p))` with G=1e6.

**XOsdi element** — SPICE `X<name> net0 ... netN model_name [key=val ...]` syntax for arbitrary-port OSDI instances. Parsed in fairchild-parser. Discipline checking via `.optical net1 net2...` directive + `check_disciplines()` (no separate `.optical_lambda` — lambda wires are declared in `.optical` alongside amplitude wires). XOsdi devices instantiated in `build_devices()` (newton.rs) by model name lookup.

**OSDI parameter matching** — `osdi_param_name_matches()` uses `eq_ignore_ascii_case` (case-insensitive). The SPICE parser lowercases all tokens, but Verilog-A identifiers preserve case (e.g. `Vpi_L`, `L_arm_um`). Without case-insensitive matching, any parameter with an uppercase letter is silently ignored and the model uses its default. Fixed in phase 2.5 session.

**`.include` pre-processing** — `resolve_includes()` in `fairchild-parser/src/spice.rs` handles `.include "path"` before the logical-line tokenizer sees the input. Recursive includes supported (depth limit 16). CLI and Python bindings use `parse_spice_file(path)` which calls `resolve_includes` first; the raw `parse_spice(&str)` remains unchanged (no include support, for programmatic use). Relative paths resolved from the referencing file's parent directory.

**Two-DC compound subckt NR bug** — A compound subckt where BOTH outputs of DC1 feed BOTH inputs of DC2 (i.e., DC1.b1→DC2.a1, DC1.b2→DC2.a2) causes NR divergence. This pattern appears in MZI modulators. Diagnostic confirmed the bug occurs even with kappa=0 (identity passthrough), ruling out coupler math. Root cause: suspected competing cross-Jacobian entries from the two OSDI instances on shared nodes. **Workaround: use monolithic `.osdi` MZI models** (mzi_modulator_pn_l1.osdi, mzi_modulator_thermo_l1.osdi). The compound `.spc` files (mzi_pn_l1.spc, mzi_thermo_l1.spc) are kept in the repo with warning banners. MRR subckts (DC + PS feedback) work correctly — the topology avoids the two-DC pattern.

**WDM bus vector expansion** — `net[M..N]` in X-element nets or `.optical` lists expands to `net_M, net_{M+1}, ..., net_N` (inclusive, `_`-separated). The `.optical_bus N re_base im_base wl_base` directive generates 3N channel-interleaved optical net entries. This is parser-level syntactic sugar only — existing VA models are still single-channel, with no cross-channel physics (XPM/FWM). `scripts/kicad_to_fairchild.py` auto-generates the wrapper netlist from a KiCAD SPICE export.

**KiCAD integration** — `kicad_integration.md` (not committed, lives at repo root) is the setup guide for KiCAD symbol library + circuit creation. It was written referencing the old `Spice_*` property system; KiCAD 7+ uses the `Sim.*` namespace (`Sim.Device`, `Sim.Name`, `Sim.Pins`, `Sim.Params`). Needs verification against what user sees in KiCAD 10 Symbol Editor and possible doc update.

---

## Build

```bash
cargo build --release
cargo test
./target/release/fairchild -f examples/rc_step.sp
```

OpenVAF-Reloaded: `/Users/hugh/Local/src/OpenVAF-Reloaded/target/release/openvaf-r`
Runtime env: `DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib`
