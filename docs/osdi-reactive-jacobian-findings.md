# OSDI Reactive Jacobian — Investigation Findings

**Branch**: `osdi-reactive-jac-investigation`  
**Status**: Blocked — root cause identified, fix approach known, implementation deferred.

---

## Problem

`cmos_inverter_switching_transient` fails: V(out) drifts ~10 mV per timestep above VDD
even in DC steady-state (before any input edge). By t=1 ns, V(out) ≈ 6.0 V, causing
LTE to saturate at ~2.0 and reject every step indefinitely.

---

## Root Cause

OSDI's `load_spice_rhs_tran` outputs the **SPICE Newton-form RHS**, not a raw current
delta. Specifically, for a B-E companion model of a capacitor C between nodes i and j:

```
b[i] += alpha * C * V_i(k)   −   alpha * C * V_j(k)   −   prev_Q/h
```

This embeds `G_react * V_k` (i.e. the reactive companion conductance times the current
iterate) into the b vector. For Newton's method to be consistent, the Jacobian matrix A
must receive the matching `alpha * C` stamp on the diagonal.

### Why the diagonal is missing

`write_jacobian_array_react` writes values for only the **last `num_reactive_jacobian_entries`
slots** in the `jacobian_entries` array. OpenVAF lays out the array as:

```
jacobian_entries[0 .. n_resist):    resistive-only entries  (e.g. drain-drain g_ds)
jacobian_entries[react_start .. n_total):  reactive entries
                                    where react_start = n_total − n_react
```

When `n_resist + n_react > n_total` there is an overlap region. For the NMOS/PMOS Level 1
models (n_resist=6, n_react=7, n_total=9, react_start=2):

| entry | nodes | in resist | in react |
|-------|-------|-----------|----------|
| [0]   | (drain, drain) | ✓ | **✗** |
| [1]   | (drain, gate)  | ✓ | **✗** |
| [2]   | (drain, source)| ✓ | ✓ |
| [3]   | (gate, drain)  | ✓ | ✓ |
| [4]   | (gate, gate)   | ✓ | ✓ |
| [5]   | (gate, source) | ✓ | ✓ |
| [6]   | (source, drain)| ✗ | ✓ |
| [7]   | (source, gate) | ✗ | ✓ |
| [8]   | (source, source)| ✗ | ✓ |

The **(drain, drain)** diagonal is in the resistive-only prefix — `write_jacobian_array_react`
never writes it. Yet `load_spice_rhs_tran` contributes `G_cgd * V_drain` to `b[drain]`.

With alpha = 1/10ps = 1e11, G_cgd = 2e-15 * 1e11 = 0.2 mS:

```
unmatched b contribution per MOSFET = 0.2 mS × 5 V = 1 mA
delta_V per step ≈ 1 mA / 100 mS (A[out][out]) ≈ 10 mV/step  ← observed ✓
```

---

## What Was Tried

1. **Wrong reactive-start offset (original code)**: `skip(n_resist)` = skip 6 → only
   3 reactive entries reached (indices 6,7,8). Fixed to `react_start = n_total − n_react = 2`.

2. **Correct reactive-start offset**: Now 7 entries are stamped, but the (drain,drain)
   diagonal is not among them (it is at entry[0], resistive-only). Drift reduced slightly
   but still ~10 mV/step.

3. **Forcing nodes from LTE**: Voltage-source-constrained nodes correctly excluded from
   LTE denominator (prevents spurious LTE rejections at PULSE edges). Working.

4. **`x_tprev` / `commit_timestep`**: Added per-device prev-timestep snapshot so
   `load_spice_rhs_tran` uses the correct previous solution, not the DC OP. Working.

---

## Fix Required

Use the **aliasing Jacobian path** instead of the copy path:

```rust
pub type FnLoadJacobianAlpha =
    unsafe extern "C" fn(inst: *mut c_void, model: *mut c_void, alpha: f64);
// descriptor.load_jacobian_tran  — stamps resistive + alpha*reactive in one call
```

This function writes directly into pre-registered pointer slots within instance memory.
To use it without a stable global MNA matrix (our `Vec<Vec<f64>>` reallocates each NR
iteration), the plan is:

1. At `setup_instance` time, inspect `jacobian_ptr_resist_offset` (and the equivalent
   reactive offset) in the `OsdiDescriptor` to find where the pointer array lives inside
   instance memory.
2. Allocate a **per-device flat Jacobian buffer** of size `num_jacobian_entries`.
3. Write the buffer's element addresses into the instance's pointer array.
4. Each NR iteration: zero the flat buffer, call `load_jacobian_tran(inst, model, alpha)`,
   then scatter-add the flat buffer into `mat.a` using `jacobian_entries[i].nodes`.

This is the "copy via temporary" approach — allocates once per device at setup, copies
once per NR iteration. Cheaper than the current two-call (resist + react) copy path and
produces the correct complete Jacobian.

### Descriptor fields needed

```rust
// in OsdiDescriptor (ffi.rs):
pub jacobian_ptr_resist_offset: u32,   // offset into instance memory of *mut f64 array
// reactive offset may be at jacobian_ptr_resist_offset + n_resist * 8 or separate field
pub load_jacobian_tran: Option<FnLoadJacobianAlpha>,  // already mapped, is Some ✓
```

Need to verify the reactive pointer offset field name in the OSDI v0.4 spec.

---

## Next Experiments

1. Print `jacobian_ptr_resist_offset` for the NMOS model to confirm its value and locate
   the reactive pointer offset field.
2. Implement the "flat buffer aliasing" path in `OsdiDevice::setup_instance` and
   `load_jacobian_tran`.
3. Verify: with correct Jacobian, V(out) should stay at 5.000 V during pre-pulse phase,
   and `cmos_inverter_switching_transient` should pass.

---

## Current Test Status

- `cmos_inverter_input_low` ✓  
- `cmos_inverter_input_high` ✓  
- `cmos_inverter_switching_transient` ✗ (V(out) drifts above VDD → LTE loop)
