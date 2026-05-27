# fairchild — Benchmark Results

This page documents fairchild's accuracy and performance relative to ngspice
on a representative set of analog and mixed-signal circuits.

Results are generated nightly by `.github/workflows/benchmarks.yml` and
can be reproduced locally:

```bash
cargo build --release
python3 benchmarks/run_all.py --output /tmp/results.json
```

---

## Circuits

| Circuit | Description | Analysis | Nodes |
|---------|-------------|----------|-------|
| RC step response | 1 kΩ / 1 µF, 1V step, 5ms | Transient | 2 |
| RLC resonator | Q≈10, f₀≈5 kHz, 1ms | Transient | 3 |
| Diode rectifier | Half-wave, 1 MHz, Shockley model | Transient | 2 |
| CMOS inverter | Level 1 N+P MOS, 100 MHz switching | Transient | 2 |
| BJT CE amplifier | NPN Gummel-Poon Level 1, 200 ns | Transient | 3 |
| Ring oscillator (3-stage) | CMOS Level 1, 10 ns | Transient | 4 |
| Ring oscillator (11-stage) | CMOS Level 1, 25 ns | Transient | 12 |

---

## Accuracy vs ngspice

The accuracy comparison runs fairchild and ngspice on the same netlist and
compares the output voltage at a representative sample time. Relative error
is `|V_fc − V_ng| / |V_ng|`. Tolerances listed are the same as the Rust
integration tests in `crates/fairchild-core/tests/`.

| Circuit | fairchild V | ngspice V | rel err | tol |
|---------|------------|-----------|---------|-----|
| RC step (t=2ms) | — | — | — | 1% |
| RLC resonator (t=0.5ms) | — | — | — | 1% |
| Diode rectifier (t=2µs) | — | — | — | 0.1% |
| CMOS inverter (t=60ns) | — | — | — | 0.2% |
| BJT CE amp (t=100ns) | — | — | — | 0.5% |

*Dashes indicate results not yet generated. Run `benchmarks/run_all.py` to populate.*

---

## Wall-clock performance

Measured on the CI runner (2-core ubuntu-latest) with `cargo build --release`.
Both fairchild and ngspice run in single-threaded mode.

| Circuit | fairchild | ngspice | ratio |
|---------|-----------|---------|-------|
| RC step | — | — | — |
| RLC resonator | — | — | — |
| Diode rectifier | — | — | — |
| CMOS inverter | — | — | — |
| BJT CE amplifier | — | — | — |
| Ring osc 3-stage | — | — | — |
| Ring osc 11-stage | — | — | — |

*Populated by nightly CI run. See Actions → Benchmarks for the latest artifact.*

---

## Known gaps

- **MOSFET Level 1 junction capacitances (Cgs, Cgd, Cbs, Cbd)** are not yet
  stamped. Switching edge timing in the CMOS inverter benchmark is unphysically
  sharp; the period comparison against ngspice is still valid.
- **BJT junction capacitances (CJE, CJC)** are not yet stamped. The CE
  amplifier benchmark compares DC steady-state, not edge timing.
- Results are compared at a single sample time, not over the full waveform.
  Full waveform RMS comparison (`compare_ngspice.py`) is the more complete
  check.

---

## Methodology

- fairchild uses variable-step Gear-2 for transient (`.options method=gear`).
  Fixed-step Backward Euler is used when `method=` is not specified.
- ngspice uses its default adaptive integrator (GEAR by default).
- Timing is wall-clock from process spawn to exit, including binary load time.
  For small circuits, process startup dominates; the timing delta is meaningful
  only for ring oscillator and larger circuits.
- The CI runner has no SuiteSparse installed; fairchild runs with the default
  faer-sparse backend. Add `--features klu` locally to compare the KLU backend.
