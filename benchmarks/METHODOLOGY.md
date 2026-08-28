# Benchmark Methodology

Rules governing how fairchild benchmark numbers are produced. Without disclosed
rules, every comparison plot is open to "you tuned the options".

## Rules

1. **Default options only.** Each simulator runs with its factory defaults — no
   tweaking `reltol`, no swapping integrator, no undocumented flags. A second
   "tuned" row is allowed in the table as a separate entry but must be labelled.

2. **Same netlist, both simulators.** fairchild and ngspice receive identical
   netlists. The only difference is the wrapper (fairchild's `-f` flag vs
   ngspice's batch `.control` block) needed to get parseable output.

3. **Median of 3 runs for timing.** Wall-clock is `time.perf_counter()` from
   process spawn to exit, including binary load and parse. Process startup
   dominates for small circuits; the timing delta is meaningful only at
   ≥20 nodes.

4. **No selective reporting.** If a circuit fails (parse error, non-convergence,
   NaN output), it appears as a red "FAILED" panel in the accuracy figure and
   a missing data point in the scaling plot — not a row that's quietly omitted.
   Failures are a feature of the report: they surface prioritisation pressure.

5. **RMS error metric.** Accuracy is `sqrt(mean((V_fc − V_ng)²))` where
   `V_fc` is linearly interpolated onto the ngspice timepoints within the
   overlapping time window. Both simulators use their native adaptive
   timesteppers, so point-by-point subtraction would be meaningless.
   Max absolute error is also computed but not shown in the main figure to
   avoid visual clutter; it is in the JSON blob.

6. **Commercial simulators (HSPICE, PrimeSim) are opt-in.** They require an
   environment variable pointing to the licence (`HSPICE_BIN`, `PRIMESIM_BIN`).
   CI never tries to run them. Their results are read from a checked-in JSON
   blob (`benchmarks/results/hspice.json`) that the user updates manually.

7. **CI environment.** Nightly runs execute on `ubuntu-latest` (2 vCPU) in
   `.github/workflows/benchmarks.yml`. Local numbers differ; prefer the CI
   numbers for published comparisons.

## Simulator versions

| Simulator | Version in CI | Notes |
|-----------|---------------|-------|
| fairchild | HEAD of `main` | Built `--release` with `--features klu` when SuiteSparse is installed |
| ngspice   | 42+ (from apt) | Default integrator (GEAR); KLU backend enabled |

## Circuit corpus

See `benchmarks/circuits/` for the full list. The short version:

- **`analog_small`** — RC step, RLC, diode rectifier, CMOS inverter, BJT CE
  amp, Schmitt trigger, ring oscillators (3/5/11 stage). Proves the solver
  matches ngspice on circuits every designer recognises.
- **`analog_scaling`** — ring oscillator family (3 → 51 stages). Provides the
  log-log scaling plot. Circuits are synthetic; accuracy is trivial. The point
  is the wall-clock curve.

## Reproducibility

```bash
cargo build --release
python3 benchmarks/run_all.py --output /tmp/results.json
python3 benchmarks/plot.py
```

Plots land in `docs/plots/`. The raw JSON from `run_all.py` can be diffed
between runs to audit regressions.

## The OSDI-vs-native comparison

`benchmarks/osdi_vs_native.py` answers a different question from `run_all.py`:
not "how does fairchild compare to ngspice" but "what does the compiled
Verilog-A path cost against a native Rust model of the same arithmetic". That
number sets the cost of the whole bring-your-own-PDK story, since BSIM, PSP,
HiSIM, HICUM and MEXTRAM all arrive through OSDI.

It follows the rules above, plus three of its own:

5. **Correctness gates the measurement.** Each model pair must agree before a
   single timing is reported, and the script exits non-zero if one does not. A
   ratio between two models that compute different things is not a measurement of
   overhead. All three pairs currently agree exactly — the Verilog-A Shockley
   diode and the native one are bit-identical.

6. **Interleaved A/B, not batched.** This machine drifts up to ~45% run to run,
   which is larger than the effect. The two sides alternate within each rep.

7. **Report absolute overhead, not only the ratio.** A ratio is relative to the
   native side, so a model whose native counterpart is slow looks *better*. The
   `nonlinear_exp` pair reads 6.1× against `resistive`'s 7.4× purely because a
   native diode is slower than a native resistor; per device the two overheads are
   231 and 232 µs, i.e. the same.

### What it found

Three pairs, chosen so the result can be attributed: a bare conductance, the same
model plus one `ddt()`, and a Shockley diode with an `exp`.

| | DC, 2048 devices | transient, 2048 devices | transient overhead |
|---|---|---|---|
| resistive (no `ddt`) | 2.40× | 7.39× | 232 µs/device |
| reactive (one `ddt`) | 1.96× | 15.05× | 594 µs/device |
| nonlinear (one `exp`) | 3.25× | 6.13× | 231 µs/device |

* **The overhead is the ABI call, not the arithmetic.** An `exp` per eval costs
  the same 231 µs/device as a bare multiply. So OSDI amortises well for a complex
  model — BSIM4 does thousands of flops per eval — and badly for a trivial one.
* **A `ddt()` term more than doubles it**, because the reactive residual and
  Jacobian are a second pair of ABI calls every timestep.
* **DC barely pays.** A handful of Newton iterations against a transient's
  hundreds of steps.
* **Unexplained, and a lead rather than a conclusion:** the per-device overhead
  *rises* with device count instead of flattening. From 512 to 2048 the OSDI side
  grows as N^1.5–1.7 where the native side grows as N^0.7–0.9. Something in the
  OSDI path scales with more than the device's own footprint.
