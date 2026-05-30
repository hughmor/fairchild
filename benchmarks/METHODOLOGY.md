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
