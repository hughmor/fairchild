# fairchild-c — C API (`libfairchild_c`)

Embed fairchild in a C/C++ program, in the spirit of `libngspice`. Two layers
over the same solver:

| Layer        | Use it for                                            | Entry points |
|--------------|-------------------------------------------------------|--------------|
| **Batch**    | run an analysis, read the waveforms                   | `fc_load_*` → `fc_run_tran` → `fc_signal` |
| **Stepping** | mixed-signal co-simulation — the host drives the clock | `fc_stepper_new` → `fc_get_node` / `fc_set_source` / `fc_step` |

The full contract is in [`include/fairchild.h`](include/fairchild.h).

## Build

```sh
cargo build -p fairchild-c --release       # -> target/release/libfairchild_c.{so,dylib,a}
cc -O2 -o prog prog.c -I crates/fairchild-c/include -L target/release -lfairchild_c
```

`cdylib` and `staticlib` are both produced. No Python runtime is linked — pyo3
lives only in `fairchild-py`. The `_c` suffix avoids an artifact collision with
that crate, which also builds as `libfairchild`; rename or symlink freely.

`--features klu` links SuiteSparse KLU; `--no-default-features` drops OSDI
(netlists with `.osdi` directives then fail loudly rather than silently ignoring
them).

## Mixed-signal loop

See [`examples/mixed_signal.c`](examples/mixed_signal.c) for a runnable version.
The shape is:

```c
fc_sim *sim = fc_sim_new();
fc_load_string(sim, netlist);
fc_stepper *st = fc_stepper_new(sim, 1e-10);  /* solves the operating point */
fc_sim_free(sim);                             /* the stepper snapshotted it */

for (long cycle = 0; cycle < n; cycle++) {
    fc_advance_to(st, (cycle + 1) * clock_period, NULL);  /* analog runs   */
    double v; fc_get_node(st, "out", &v);                 /* analog -> digital */
    fc_set_source(st, "VDRIVE", decide(v));               /* digital -> analog */
}
fc_stepper_free(st);
```

`fc_set_source` is a zero-order hold, matching ngspice's `GetVSRCData`
semantics: the value stands until the next write.

## Things worth knowing

* **No global state.** Every handle is independent, so N concurrent simulations
  in one process just work — no `dlmopen` tricks. A single handle is *not*
  thread-safe; use one per thread.
* **The step size is fixed** for a stepper's lifetime. `fc_advance_to` lands on
  the first grid point at or past the target, so digital events quantise to the
  analog grid. Pick the step from the shortest edge you need to resolve. (Off-grid
  event times would need companion models rebuilt for a new `h` — see the
  `ponytail:` note on `TranStepper`.)
* **Errors, not crashes.** NULL handles and NULL strings return `FC_ERR_ARG`;
  panics are caught and become `FC_ERR_PANIC`. On `FC_ERR_SIM` from `fc_step`
  the stepper still holds the last accepted timepoint, so a host can back off a
  drive level and retry.
* **Borrowed pointers.** `fc_signal`, `fc_signal_name`, and `fc_error` hand back
  pointers into the handle's storage. Don't free them; they die with the next
  call on that handle.

## Cost of the boundary

Measured on an M-series mac, release build, per call:

| Call            | Cost      |
|-----------------|-----------|
| `fc_time`       | ~1.2 ns   |
| `fc_get_node`   | ~21 ns    |
| `fc_set_source` | ~36 ns    |

Against a timestep, which costs ~0.6 µs for a trivial RC and ~58 µs for a
201-stage ring oscillator, a read/write pair per step is under 10% on the
smallest circuit and unmeasurable on real ones. Stepping is in fact *cheaper*
per step than the batch path (~5–8%), because it doesn't accumulate every node
at every timepoint into a result table.

If a host ever crosses a wide bus every step, the fix is resolve-once handles
rather than a cleverer name lookup — an exact-match fast path was tried and
measured slower, because parsed netlist names are already lower case.
