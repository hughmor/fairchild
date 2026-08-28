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
  OSDI path scales with more than the device's own footprint. Probably *not* the
  same cause as the internal-node pathology below — those OSDI fixtures are
  two-terminal and allocate no internal nodes.

## What a device's internal node costs (#99, fixed)

Worth recording because the first three measurements were all wrong, and each was
wrong in a way that looked convincing.

The question came out of #77 §2's `RD`/`RS` work: internal unknowns belong in the
matrix rather than being eliminated (three separate silent wrong answers say so),
so what does a row cost? Measured on 256 MOSFETs with a 2 kΩ load each,
`.tran 5n 2u`, the only change being `RD=50 RS=50`:

| | before | after |
|---|---|---|
| `.op` | 1.37× | 1.37× |
| `.tran` | **10.3×** | **1.7×** |

### The eliminations, in order

Each of these looked like the answer and was not:

1. **"It is the step controller."** No — the timepoint count is *identical*, 402
   with and without, at every device count.
2. **"It is solver scaling with row count."** No — a plain linear RC ladder at
   fixed step is linear: 86–105 µs per node from 256 to 4096 nodes.
3. **"It is the O(rows²) clique footprint."** No — nnz tracks rows
   proportionally, 2.9 per row either way.
4. **"It is stiffness."** No — sweeping the series resistance over a 100× range
   of conductance ratio gives 13.52, 14.06, 14.04, 14.06. Flat.
5. **"It is the reactive path."** No — present with no capacitors in the deck.

### What it actually was

Profiling, not reasoning. Two `sample` runs of the same circuit differing only by
`RD=50 RS=50`:

```
no internal nodes   ->  20 `simplicial` frames,  0 dense-LU frames
with internal nodes ->  18 `supernodal` frames, 471 `lu_in_place_recursion`
```

An internal node changed the sparsity structure enough to flip faer's flop-ratio
heuristic from simplicial to **supernodal**, whose dense-block kernel is the wrong
shape of work for a circuit matrix. KLU — which never takes a dense path — showed
1.2–1.6× on the same decks, which is what said the work was fill-in and not the
rows. `SUPERNODAL_THRESHOLD` in `solver.rs` now pins the choice.

### Confirming the fix did not cost anything

Every deck in `benchmarks/circuits/` plus four photonic examples, interleaved A/B,
5 reps: worst ratio **1.010×** (`rc_step`, which is noise — 4.9 → 5.0 ms), and
most are 0.93–0.98×, i.e. *faster*. Output **byte-identical** on all 16.

### And it predates the change that exposed it

The BJT's `RB`/`RC`/`RE` internal nodes are years old and #98 did not touch the
BJT's stamping. Against the pre-#98 binary the BJT ladder is **4.98×**; against
the post-#98 one, **5.05×**. So #98 exposed the pathology for a second device
family rather than introducing it — worth stating because the issue was filed
before that was known.
