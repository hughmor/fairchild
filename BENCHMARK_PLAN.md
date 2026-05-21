# Fairchild Benchmark & Competitiveness Plan

*Companion to `sotu.md` §7. This file is the concrete actionable plan for
making the GitHub repo's claim of "credible alternative to ngspice / Spectre
Photonics" visibly true.*

---

## 0. The user question this answers

A first-time visitor to the GitHub repo lands on the README, scrolls for
30 seconds, and asks themselves: **"Is this real?"** Today there are
five plots showing fairchild matches ngspice on three-resistor circuits.
That doesn't answer the question — it's too small a sample to distinguish
"works on tutorials" from "works on real designs". The answer the visitor
wants is a single benchmarks page they can scan in 60 seconds and walk
away from with one of: *yes credible*, *yes for photonics specifically*,
*not yet*.

This plan is the bare-minimum content for that page, the rules for
producing it honestly, and the build infrastructure to keep it fresh.

---

## 1. What we benchmark, and against whom

### 1a. The four categories

| Category | What it measures | Primary baseline | Stretch baselines |
|---|---|---|---|
| **Analog accuracy** | RMS / max error on a fixed set of circuits | ngspice | HSPICE, PrimeSim |
| **Analog performance** | Wall-clock + peak RSS, scaling with N | ngspice | HSPICE, PrimeSim |
| **Convergence robustness** | Pass/fail + iteration count on hard nonlinear circuits | ngspice | HSPICE, PrimeSim |
| **Photonic accuracy** | Transmission, eye opening, phase vs analytic + commercial | Analytic (CMT, S-matrix) | Lumerical INTERCONNECT, SAX, Spectre Photonics if available |

### 1b. Why these four

- **Accuracy** is the threshold question. If you can't trust the numbers,
  speed is irrelevant.
- **Performance** is what determines whether a designer reaches for
  fairchild instead of ngspice on the second simulation, not just the
  first.
- **Convergence** is what separates a research solver from a tool.
  Half of ngspice's reputation comes from "it converges on circuits
  Spice2 didn't".
- **Photonic** is the differentiator and the reason for existence.
  Skipping a serious photonic baseline forfeits the entire moat.

### 1c. Baseline simulators we can realistically use

| Simulator | Availability | Notes |
|---|---|---|
| **ngspice** | Open source on every dev machine | Already used by `compare_ngspice.py`. The credibility floor. |
| **Xyce** | Open source (Sandia) | Secondary open peer; tests fairchild against another implementation of the same standard. Costs nothing. |
| **HSPICE** | User has access | Gold standard analog signoff. Probably needs a separate machine; license file constraints. |
| **Synopsys PrimeSim** | User can try to obtain | Modern Synopsys flagship; sparse-LU and multi-rate stories worth comparing on. |
| **Cadence Spectre / Spectre Photonics** | Not available today | The photonic gold standard. Defer; revisit if academic licence appears. |
| **Lumerical INTERCONNECT** | Free Ansys evaluation tier | Frequency-domain S-matrix; the open-photonic baseline. Fair comparison only on steady-state spectra, not on time-domain dynamics. |
| **SAX** | Open source (JAX) | Pure freq-domain peer in Python. Easiest to run. Comparison highlights time-domain differentiation. |
| **Simphony** | Open source | Other freq-domain peer. Optional, redundant with SAX. |
| **Photontorch** | Open source | RNN-based; comparison only meaningful on circuits both can express. Skip in v1. |

The benchmark table in the README needs **at least three columns**
(fairchild, ngspice, HSPICE) for analog and **at least two**
(fairchild numerical, analytic closed-form) for photonic. PrimeSim,
Xyce, INTERCONNECT, SAX go in extended tables that the README links
out to but does not embed.

---

## 2. The circuit corpus

### 2a. Analog — small (the smoke-test floor)

These prove fairchild matches ngspice on textbook circuits. They are
where the project is today.

| File | Type | Nodes | Why |
|---|---|---|---|
| `rc_step.sp` | TR | 3 | First-order LPF settling |
| `rlc_resonator.sp` | TR | 4 | Second-order ringing / damping |
| `rc_bode.sp` | AC | 3 | Magnitude + phase vs analytic |
| `diode_rectifier.sp` | TR | 3 | Nonlinear half-wave |
| `cmos_inverter.sp` | TR | 4 | MOSFET-L1 + Miller |
| `nmos_dc_sweep.sp` | DC | 4 | DC sweep robustness |
| `ring_oscillator.sp` | TR + DC OP | 16 | Metastable DC + oscillation |

All shipping today. Keep them; they're the regression floor.

### 2b. Analog — medium (the credibility tier)

These prove fairchild handles things a designer recognises. None exist
in the repo today.

| Circuit | Source | Why it matters |
|---|---|---|
| **uA741 op-amp** | Public-domain SPICE deck, ships with ngspice tutorials | The most-simulated analog circuit in history. Tests BJT badly — if BJT isn't shipping, this fails on parse. |
| **Folded-cascode CMOS op-amp** | Razavi textbook netlist | All-MOSFET version of above. ~25 nodes. Common designer reference. |
| **Bandgap reference** | Razavi / Hastings | Startup-sensitive nonlinear DC. Tests homotopy. |
| **PLL charge pump** | Behavioral + analog mix | Stresses B-element + integrators. |
| **5th-order Chebyshev LC filter** | Standard filter design | Tests `K`-coupled inductors (currently unsupported — generates a Tier-1 backlog row). |
| **Schmitt trigger** | Razavi | Convergence on hysteretic nonlinearity. |
| **CMOS NAND chain (8 / 16 / 32 stages)** | Synthesised | Performance-scaling-friendly. Tests sparse LU vs dense. |

These ship in `benchmarks/circuits/analog/medium/` as part of the first
benchmark PR. Three of them (uA741, Chebyshev, PLL charge pump) will
fail on parse today — that failure surfacing on a public page is
exactly the prioritisation signal we want.

### 2c. Analog — large (the scaling tier)

The point of these is the log-log scaling plot. Synthesise them.

| Circuit family | N range | Generator |
|---|---|---|
| **Ring oscillator** | 3, 5, 11, 21, 51, 101 stages | Trivial Python loop |
| **CMOS inverter chain** | 8, 16, 32, 64, 128 stages | Python loop |
| **Resistor ladder R-network** | 50, 100, 500, 1 000, 5 000 nodes | Tests sparse LU, no NR |
| **Capacitively-loaded clock distribution H-tree** | 4-level → 7-level | Tests reactive scaling |

These are not "real circuits" — they're scaling probes. The point is
not their accuracy (it's trivial); the point is the wall-clock curve.

### 2d. Photonic — analytic comparisons

| Circuit | Comparison | Why |
|---|---|---|
| **MRR wavelength sweep** | CMT closed form | Already in repo; tighten error bars |
| **MZI bias sweep** | cos² closed form | Test `fc_splitter` + 2× `fc_pn_ps` + `fc_splitter` topology |
| **WDM 4-channel MRR** | Per-channel CMT | Tests bundle-aware shared-electrical correctness |
| **Bidirectional ring (drop port)** | Sorace-Agaskar 2015 reported eq. | First public-result comparison |

### 2e. Photonic — commercial / open comparisons

| Circuit | Baseline | Notes |
|---|---|---|
| **MRR through-port spectrum** | Lumerical INTERCONNECT (free eval) | Spectra should overlay to < 0.5 % on resonance |
| **8-channel WDM transceiver** | INTERCONNECT for steady-state spectrum; SAX for differentiable freq-domain | Pure fairchild can produce a time-domain eye that neither baseline can |
| **25 Gb/s NRZ eye through MRR mod + waveguide + PD** | Analytic eye from the same MRR transfer | Reference open-source eye; published Sorace-Agaskar 2015 reports similar |
| **Sorace-Agaskar 2015 PDH loop replication** | Published optical / electrical waveforms | Highest-credibility photonic comparison available open source |

The Sorace-Agaskar 2015 paper is the *single most leverage-bearing*
comparison we can make: it is the seminal time-domain EO co-simulation
paper, it has published waveforms, and reproducing them in fairchild
is the photonic equivalent of "we boot to a login prompt" — a
necessary prerequisite for being taken seriously by the photonic
research community.

---

## 3. The methodology rules

Honest comparisons need disclosed rules. Without them, every plot is
open to "but you tuned the options". The rules below ship as
`benchmarks/METHODOLOGY.md`.

1. **Each simulator runs with its default options.** No tweaking
   `reltol`, no swapping integrator. The point is the out-of-the-box
   experience. A second optional run with "best tuned options" is
   allowed *as a separate plot row* but must be labelled.
2. **Same machine, same toolchain, no GPU.** Benchmarks run in the CI
   container or on a dedicated bare-metal box; never on the developer
   laptop.
3. **Three runs, report median.** Avoids noise from JIT warmup,
   filesystem caching, etc.
4. **No selective reporting.** If a benchmark fails (parse error,
   non-convergence, NaN), that's a row in the table marked "fail",
   not a row that's missing. Failures are a feature of the report.
5. **Output diff metric is RMS error on the probe channel, plus
   max absolute error.** Both are reported because RMS hides large
   localised errors.
6. **Wall-clock is end-to-end** — including parse time. Designers
   notice startup time; subtracting it would advantage fairchild
   unfairly on small circuits.
7. **HSPICE and PrimeSim runs require an environment variable
   pointing to the licence.** They're skipped (not failed) when
   absent. Their results are read from a checked-in JSON blob that
   the user updates locally; CI never tries to run them.
8. **Every benchmark commit also commits the raw output files** so
   the diffs are reproducible after the fact.

---

## 4. The plots — what they look like

Six plots ship to `docs/benchmarks/`; the README embeds the first
three thumbnails and links to the rest.

### 4a. `accuracy_analog_small.png`

6-panel figure: rc_step, rlc_resonator, diode_rectifier, cmos_inverter,
nmos_dc_sweep, ring_osc_5. Each panel overlays 2–4 curves (fairchild,
ngspice, optional HSPICE/PrimeSim). Below each panel: RMS / max error
in units the designer thinks in (mV / µA / dB).

### 4b. `accuracy_analog_medium.png`

4-panel figure: uA741 unity-gain transient, folded-cascode bias sweep,
bandgap startup, Schmitt trigger hysteresis. Same overlay format.
**Some panels will be blank or show fairchild errored** — that's the
point, prioritisation pressure.

### 4c. `scaling_log_log.png`

Log-log wall-clock vs node count. Lines: fairchild-dense,
fairchild-sparse, ngspice, HSPICE, PrimeSim. Two facets (one for the
ring-oscillator chain, one for the resistor ladder). The reader
should be able to read "at N=1000 fairchild-sparse is 1.4× ngspice"
from this directly.

### 4d. `convergence_bar.png`

Ten "hard convergence" circuits as rows. Bars per simulator showing:
green (converged in default options), yellow (converged with
non-default options), red (failed). Annotated with NR iteration count
when converged. Shows where fairchild's diagnostics + Armijo damping
buy something the others don't.

### 4e. `photonic_mrr_overlay.png`

MRR through-port transmission vs wavelength: CMT analytic + fairchild
+ Lumerical INTERCONNECT overlaid. Inset zoom on resonance.

### 4f. `photonic_wdm_eye.png`

25 Gb/s NRZ eye after MRR modulator + PD: fairchild time-domain
overlaid on the analytic eye. No other open-source tool can produce
the time-domain version, so the "differentiator" framing is justified
in caption.

---

## 5. The three README tables

### 5a. Feature parity

Embedded directly in the README under a "Simulator comparison"
section. Categories: SPICE primitives, directives, analyses,
specialty (photonic / S-param / adjoint).

```
| Feature       | fairchild | ngspice | Xyce | HSPICE | PrimeSim | Spectre |
|---------------|-----------|---------|------|--------|----------|---------|
| R/L/C/V/I/D   | ✓         | ✓       | ✓    | ✓      | ✓        | ✓       |
| MOSFET-L1     | ✓         | ✓       | ✓    | ✓      | ✓        | ✓       |
| BJT GP        | ✗         | ✓       | ✓    | ✓      | ✓        | ✓       |
| K-coupled L   | ✗         | ✓       | ✓    | ✓      | ✓        | ✓       |
| .tran .dc .ac | ✓         | ✓       | ✓    | ✓      | ✓        | ✓       |
| .noise        | ✓         | ✓       | ✓    | ✓      | ✓        | ✓       |
| .measure      | ✓         | ✓       | ✓    | ✓      | ✓        | ✓       |
| BSIM4 OSDI    | ✓ (compat)| ✓       | ✓    | n/a    | ✓        | ✓       |
| Photonic time | ✓         | ✗       | ✗    | ✗      | ✗        | ✓       |
| Adjoint       | (Phase 4) | ✗       | ✗    | ✗      | (limited)| ✗       |
| Python API    | ✓         | ✗       | ✓    | ✗      | limited  | limited |
| Cost          | free      | free    | free | $$$    | $$$      | $$$$    |
```

### 5b. Accuracy summary

One row per benchmark, columns: RMS error vs ngspice, max error,
relative to ngspice's reported uncertainty band.

### 5c. Performance summary

Single-number wall-time on three target circuits per simulator.
The 3-circuit cut: 5-stage ring osc (small), uA741 (medium), 32-stage
inverter chain (large). One row per simulator.

---

## 6. Build infrastructure — the concrete first PR

Minimum landable slice that gets the benchmark page live.

### 6a. New directory tree

```
benchmarks/
├── METHODOLOGY.md               ← rules of §3
├── circuits/
│   ├── analog_small/            ← move from examples/electronic/
│   ├── analog_medium/           ← new
│   ├── analog_scaling/          ← synthesised, generators in scripts/
│   └── photonic/                ← move featured examples + new
├── results/
│   ├── ngspice/                 ← ngspice CSV outputs, committed
│   ├── fairchild/               ← fairchild CSV outputs, committed
│   ├── hspice/                  ← from HSPICE_BIN env, committed by hand
│   └── primesim/                ← same
├── run_all.py                   ← single entry point: fairchild + ngspice always; HSPICE / PrimeSim if env vars set
├── plot.py                      ← generates the 6 plots from results/
└── compare.py                   ← computes the accuracy + perf tables

docs/benchmarks/
├── index.md                     ← embeds plots + tables; linked from README
└── (auto-generated png + svg)
```

### 6b. CI workflow

Single `.github/workflows/benchmarks.yml`:

- Runs on every push to `main` and nightly via `schedule`.
- Steps: cargo build, maturin build, install ngspice from apt,
  `python benchmarks/run_all.py`, `python benchmarks/plot.py`.
- On `main`: pushes generated PNGs to a `gh-pages` branch and
  updates a README badge with the date + commit-hash of the latest
  benchmark run.
- On PR: posts a comment with the accuracy + perf table delta vs
  `main`.

### 6c. Two GitHub badges in the README

- `Benchmarks: 2026-05-17 #abc1234` — links to `docs/benchmarks/`.
- `vs ngspice: 0.3 % RMS, 1.2× wall` — the single number summary,
  refreshed nightly.

---

## 7. Sequencing

1. **Week 1** — Move `compare_ngspice.py` / `benchmark.py` into
   `benchmarks/`. Add `benchmarks/METHODOLOGY.md`. Write `run_all.py`
   and `plot.py` for the small analog circuits already in the repo.
   Wire up CI. Ship the first version of `docs/benchmarks/index.md`
   showing only the small analog circuits + the photonic MRR overlay.
   *This alone moves the README a long way.*
2. **Week 2** — Generate the analog-medium circuits (uA741, folded
   cascode, bandgap, Chebyshev, Schmitt trigger). Several will fail
   on parse. Add the failure rows to the accuracy table. The failures
   inform Tier-1 priority order: the most-frequently-blocking missing
   element rises to top of backlog.
3. **Week 3** — Generate the scaling circuits and ship the log-log
   scaling plot. This is where sparse LU is exercised for the first
   time at the size it was designed for.
4. **Week 4** — Photonic: WDM-MRR overlay vs Lumerical INTERCONNECT
   (free eval); Sorace-Agaskar PDH reproduction.
5. **Week 5 onwards** — User-owned HSPICE / PrimeSim runs, committed
   as JSON blobs. Update the comparison table.

After week 4 the GitHub front door looks credible. The work after
that is incremental, every accuracy and convergence win shows up as
a row in the public table on the next nightly run.

---

## 8. What this plan deliberately does *not* try to do

- **Beat HSPICE on speed.** Not a realistic outcome v1; the goal is
  parity-within-a-small-multiple. Honest reporting of "we are 1.5× HSPICE
  on this circuit family" is more valuable than a misleading win on a
  cherry-picked one.
- **Compete with FDTD on rigour.** fairchild is a circuit simulator;
  Meep / Lumerical FDTD solve Maxwell's equations. Different problems.
  The comparison story is "fairchild is faster than FDTD on EO co-sim
  because it doesn't solve Maxwell's equations", not "fairchild is more
  accurate than FDTD".
- **Cover BSIM4 + foundry PDK transistor accuracy.** Out of scope for
  the v1 benchmark page. That's a separate "PDK compatibility" story
  later.
- **Run benchmarks across operating systems.** Linux x86_64 only in
  v1. macOS / Windows results would be valuable but cost CI minutes
  for marginal trust gain.
- **Add yet-another-frontend.** No web UI, no interactive plot. PNGs
  and tables in GitHub-rendered markdown.

---

## 9. The single most important graph

If only one graph could ship: **`scaling_log_log.png`** (§4c).
Reason: it answers all three of the visitor's 60-second questions in
one image — accuracy is implicit (the runs converged), performance is
explicit (the y-axis), scaling is the *headline* (the x-axis). A
designer who sees fairchild-sparse parallel to ngspice from 100 nodes
upward will read that as "this is a real solver" and reach for it the
next time they have a photonic-adjacent question. A designer who sees
fairchild diverging from ngspice at N=200 will know exactly what to
expect — and that honest signal builds more trust than a polished
demo on three nodes.

That graph, and the work it implies (synthesising the scaling
circuits, getting sparse LU exercised for the first time at scale,
profiling what slows down at N=1000), is the highest-leverage benchmark
work on the list.
