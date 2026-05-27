# fairchild — Benchmark Results

Fairchild accuracy and performance vs ngspice on a representative set of analog
circuits. All numbers produced by `benchmarks/plot.py` on an Apple M-series
laptop (single-threaded, fairchild built `--release`). See
[`benchmarks/METHODOLOGY.md`](../benchmarks/METHODOLOGY.md) for the rules.

Reproduce locally:
```bash
cargo build --release
python3 benchmarks/plot.py
python3 benchmarks/run_all.py --output /tmp/results.json
```

---

## Accuracy vs ngspice

![Accuracy overlay](plots/accuracy_analog.png)

Six panels, each showing the fairchild waveform (blue, solid) vs ngspice (red,
dashed). RMS error is computed by interpolating fairchild onto ngspice's
adaptive timepoints within the shared time window.

| Circuit | Waveform RMS error | Point-sample rel. error | Notes |
|---------|--------------------|------------------------|-------|
| RC step response | 0.02 mV | 6×10⁻⁶ | Linear — essentially identical |
| RLC resonator | 16.15 mV | 9.4×10⁻⁴ | Small phase drift from different adaptive steppers |
| Diode rectifier | 0.05 mV | 9.8×10⁻⁵ | Nonlinear — excellent agreement |
| CMOS inverter | 201 mV | 0 (sampled in steady state) | Edge timing differs; see note below |
| BJT CE amplifier | 503 mV | 0 (sampled in steady state) | Edge timing differs; see note below |
| Schmitt trigger | 324 mV | 2.8×10⁻⁴ (in HIGH state) | Hysteresis behavior matches; edge timing differs |

**Note on switching circuits.** The CMOS inverter and BJT CE amp RMS errors are
dominated by edge-timing offset, not waveform shape disagreement. Fairchild's
MOSFET Level 1 and BJT Gummel-Poon models do not yet stamp junction
capacitances (Cgs/Cgd/Cbs/Cbd for MOSFET; CJE/CJC for BJT). This makes
switching edges unphysically fast; the steady-state voltages and oscillation
frequency are correct. Implementing junction caps is the next correctness item;
see `sotu.md` §9b.

---

## Performance vs ngspice

![Scaling wall-clock](plots/scaling_wall_time.png)

CMOS ring oscillator family (3 → 51 stages, Level-1 MOSFET), transient
simulation. Wall-clock is end-to-end (binary load + parse + simulation).
Median of 3 runs.

| Circuit | Fairchild | ngspice | Ratio |
|---------|-----------|---------|-------|
| Ring osc 3-stage  (7 nodes)  |   4 ms | 16 ms | 4.0× |
| Ring osc 5-stage  (11 nodes) |   7 ms | 25 ms | 3.6× |
| Ring osc 11-stage (23 nodes) |   9 ms | 30 ms | 3.3× |
| Ring osc 21-stage (43 nodes) |  40 ms | 103 ms | 2.6× |
| Ring osc 51-stage (103 nodes)| 377 ms | 556 ms | 1.5× |

Fairchild's advantage is largest on small circuits (lower startup overhead)
and converges toward ~1.5× at 50+ stages. The scaling exponent is similar
for both simulators, indicating equivalent algorithmic complexity.

---

## Known gaps

- **MOSFET Cgs/Cgd/Cbs/Cbd** not yet stamped → switching edge timing is wrong above ~1 MHz.
- **BJT CJE/CJC** not yet stamped → same limitation at ~100 MHz.
- Ring oscillator benchmark uses synthetic circuits; no real-circuit scaling data yet.
- ngspice comparison is on macOS developer hardware, not the CI runner.
  CI numbers (ubuntu-latest) will differ; see Actions → Benchmarks for nightly results.

---

## Methodology

See [`benchmarks/METHODOLOGY.md`](../benchmarks/METHODOLOGY.md) for the full
rules. Short version: default options, identical netlists, median of 3 runs,
failed circuits shown as red panels rather than omitted.
