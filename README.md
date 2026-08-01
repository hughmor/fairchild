# fairchild

A SPICE-compatible analog circuit simulator written in Rust, with a native
time-domain electro-optic co-simulation discipline.

**Why it exists.** Every open-source photonic simulator (SAX, Simphony,
Photontorch) is frequency-domain S-matrix only. Cadence Spectre Photonics is
the only time-domain electro-optic co-simulator and costs ~$100k/seat/year.
fairchild aims to be the first credible open-source alternative — built around
a real SPICE engine so the same Newton-Raphson loop handles electronics and
photonics in lockstep.

**Status.** Solver foundations complete (`.options`, `.dc`, `.ic`, `.measure`,
`.lib`, `.noise`, `.alter`, `.temp`, B-element behavioral sources, GEAR/BDF-2,
sparse LU, junction limiting, Armijo-damped Newton, verbose diagnostics).
Photonic discipline rebuilt around 14 native Rust devices, bundle-port
syntax with WDM as the default, and optional bidirectional propagation.
KiCad schematic capture wired up around native devices. Python bindings
cover every analysis the CLI does. See [`docs/user-guide.md`](docs/user-guide.md)
for the full feature reference, [`docs/model_status.md`](docs/model_status.md)
for the per-parameter contract (what is parsed vs actually stamped vs
validated), and [`docs/benchmarks.md`](docs/benchmarks.md) for accuracy and
performance vs ngspice.

---

## What works today

### SPICE compatibility

| Category | Coverage |
|---|---|
| Elements | R, L, C, V, I, D, K (coupled inductors), MOSFET (Level 1), BJT (Gummel-Poon), B (behavioral), S/W (switches), T (lossless line), X (subckt / OSDI) |
| Sources | DC, PULSE, PWL, SIN, EXP, SFFM, AM |
| Analyses | `.op`, `.dc`, `.tran`, `.ac`, `.noise` |
| Solvers | NR with `pnjlim` / `fetlim`; BE, TR, GEAR (BDF-2); dense or sparse LU |
| Directives | `.options`, `.ic`, `.nodeset`, `.measure`, `.lib`/`.endl`, `.include`, `.param`, `.subckt`/`.ends`, `.temp` (sweep), `.alter`, `.model`, `.osdi` |
| Output | CSV (stdout / file), Nutmeg rawfile (ngspice-compatible) |
| Verilog-A | compiled to OSDI v0.4 by OpenVAF-Reloaded; `.osdi` + `.model` cards, electrical **and** optical models ([guide](docs/user-guide.md#14-verilog-a-models-osdi)) |

What's not yet supported: lossy transmission lines (lossless `T` and
`S`/`W` switches are supported), `.disto`, `.pz`, native `.mc` Monte Carlo,
PSF/FSDB binary output.

### Electro-optic co-simulation

Native Rust photonic devices using a slowly-varying-envelope `(re, im, λ)`
representation; bundle-port syntax (`.optical_port NAME [N]`) so a 4-port
device is a 4-port symbol, not a 12-pin one; WDM via the parser's bus vector
expansion.

| Device | Card name |
|---|---|
| CW laser | `fc_cw_laser` |
| Waveguide | `fc_waveguide` |
| 2×2 directional coupler | `fc_dcoupler` |
| Behavioural 2×2 transfer block | `fc_optical_2x2` |
| Y-splitter | `fc_splitter` |
| Grating coupler | `fc_grating_coupler` |
| 3-port circulator (bidir) | `fc_circulator` |
| WDM mux / demux | `fc_mux` / `fc_demux` |
| PN-junction phase shifter | `fc_pn_ps` |
| PN phase shifter + C_j(V) | `fc_pn_ps_cap` |
| Thermo-optic phase shifter | `fc_thermal_ps` |
| Thermo-optic PS + τ_th | `fc_thermal_ps_rc` |
| Combined PN + thermal PS | `fc_pn_th_ps` |
| Idealised testbench MZM | `fc_mzm` |
| Photodetector | `fc_photodetector` |

Every bundle-aware device handles WDM by default: N-channel optical bus
in, N parallel propagation paths inside, one shared electrical interface.
Bidirectional propagation is enabled with `.options enable_bidirectional=1`.
Higher-level structures (micro-ring resonators, MZIs) are composed in the
netlist from these primitives; see `examples/photonic/`.

### Python bindings

`pip install` produces a `fairchild` package exposing every analysis. The
same SimOptions knobs available on the CLI are kwargs to `Circuit.run`.

```python
import fairchild
c = fairchild.Circuit()
c.load("examples/photonic/native_mrr_modulator.sp")
result = c.run("tran", step=5e-9, stop=2e-6, method="gear", reltol=1e-4)
import matplotlib.pyplot as plt
plt.plot(result.time(), result["V(pd_anode)"])
```

`Circuit.sweep(param, values, analysis, ...)` parametrises any element value;
this is the Monte Carlo and corner-sweep path today (a native `.mc` directive
is on the roadmap).

---

## Quickstart

### Build

```bash
cargo build --release
```

For Python bindings:

```bash
cd crates/fairchild-py
maturin develop --release
```

The Rust toolchain is pinned by `rust-toolchain.toml`, so rustup installs the
matching version on the first build. Enable the pre-commit hook once per clone —
it runs the same `cargo fmt` and `cargo clippy -D warnings` that CI fails on:

```bash
git config core.hooksPath .githooks
```

### Run a simulation

```bash
# DC operating point
./target/release/fairchild -f examples/electronic/nmos_dc_sweep.sp

# Transient with GEAR integrator and tightened tolerance
./target/release/fairchild -f examples/electronic/rc_step.sp \
    --opt method=gear --opt reltol=1e-5

# Nutmeg output for ngspice / spyci
./target/release/fairchild -f examples/electronic/rlc_resonator.sp \
    --format nutmeg -o rlc.raw

# Photonic transient — single-channel MRR modulator
./target/release/fairchild -f examples/photonic/native_mrr_modulator.sp \
    --probe "V(pd_anode),V(vmod)"
```

See [`docs/user-guide.md`](docs/user-guide.md) for full directive, element,
and CLI reference.

---

## Example circuits

```
examples/
├── electronic/
│   ├── rc_step.sp                ← RC step response (τ = 1 ms)
│   ├── rlc_resonator.sp          ← RLC series resonator (f₀ ≈ 5 kHz)
│   ├── rc_bode.sp                ← AC magnitude/phase
│   ├── diode_rectifier.sp        ← half-wave rectifier
│   ├── cmos_inverter.sp          ← Level 1 NMOS + PMOS, VDD = 3.3 V
│   ├── nmos_dc_sweep.sp          ← .dc sweep over V_GS
│   ├── ring_oscillator.sp        ← 5-stage CMOS ring; tests DC + GEAR
│   └── bjt_ce_amplifier.sp       ← NPN common-emitter amp (BJT GP L1)
├── photonic/
│   ├── native_mrr_modulator.{sp,py}        ← electro-optic micro-ring
│   ├── native_mrr_wavelength_sweep.py      ← parametric λ sweep
│   ├── native_wdm_mrr_modulator.{sp,py}    ← 2-channel WDM through one ring
│   └── legacy/                             ← archived early photonic examples
└── verilog_a/
    ├── models/*.va                         ← 8 maintained Verilog-A models
    ├── rectifier.sp, cmos_inverter.sp      ← electrical: VA device + native parts
    ├── eam_link.sp, va_link.sp             ← optical: VA in a native link, and all-VA
    └── check.py                            ← runs them all, asserts the physics
```

Photonic examples come with a `README.md` describing the topology of each
circuit. The WDM example shows two lasers detuned ±50 pm sharing one ring;
the same V_pn pulse produces very different transmission profiles per
channel — no MUX/DEMUX device needed because the bundle-port abstraction
makes each wavelength its own independent SVEA channel.

---

## Benchmarks

Head-to-head accuracy and performance vs ngspice. See [`docs/benchmarks.md`](docs/benchmarks.md) for the full table and methodology.

### Accuracy (RMS error vs ngspice)

![Accuracy overlay](docs/plots/accuracy_analog.png)

Each panel carries a residual strip (fairchild − ngspice, on ngspice's
timebase), because two curves drawn on top of each other look identical at
1 mV and at 100 mV alike.

Linear circuits (RC, RLC, diode) match ngspice to sub-1 mV RMS. The switching
circuits' larger RMS is **edge timing, not offset** — which is exactly what the
residual strips show: the error is a spike at each transition and flat between
them. A fixed-step run resolves a 1 ns edge to within one step, so the two
simulators disagree only about *when* the edge lands, and a finer `.tran` step
shrinks it. BJT CJE/CJC and MOSFET Meyer + depletion caps are all stamped now
(see [`docs/model_status.md`](docs/model_status.md)).

### Performance scaling

![Scaling plot](docs/plots/scaling_wall_time.png)

Transient wall-clock on the CMOS ring oscillator family (Level-1 MOSFET,
3–499 stages, n ≈ 7–999 nodes), each fairchild solver backend forced in turn.
At 499 stages:

| backend | wall-clock | vs ngspice |
|---|---|---|
| fairchild — KLU | 2.96 s | **6.5× faster** |
| fairchild — sparse LU (faer) | 6.38 s | 3.0× faster |
| ngspice (default) | 19.2 s | — |

KLU is the fastest backend at every size and pulls away as circuits grow; dense
LU is capped at 250 nodes because it is O(N³). The MNA matrix is stored sparse,
so both sparse backends allocate and factorise O(nnz) rather than O(n²).

Reproduce:
```bash
cargo build --release
python3 benchmarks/plot.py    # generates docs/plots/*.png
python3 benchmarks/run_all.py  # generates point-sample JSON for tables
```

## Validation

Golden tests in `crates/fairchild-core/tests/` cover:

- ngspice comparison on R/L/C, diode, MOSFET-L1 circuits (skipped, not
  failed, when ngspice is absent on PATH).
- Hard convergence: 5-stage CMOS ring oscillator (DC OP at metastable
  point; transient oscillation; BE vs GEAR agreement).
- AC: RC low-pass Bode magnitude/phase vs analytic; RLC resonance peak.
- DC: differential pair bias point + differential response.
- Photonic: MRR resonance vs CMT, MRR wavelength sweep, WDM symmetry +
  asymmetry-under-bias.

```bash
cargo test --release
```

---

## Project structure

```
crates/
  fairchild-core/     DAE solver, MNA, Newton-Raphson, transient (BE/TR/GEAR/var-step),
                      AC, .dc, .noise, .measure, SimOptions, sparse/dense LU,
                      native photonic devices.
  fairchild-parser/   SPICE parser: all directives + B-element expression grammar
                      + bus vector expansion + `.optical_port` bundles.
  fairchild-cli/      Binary `fairchild`: `-f netlist.sp`, `--format`, `--probe`,
                      `--param`, `--opt key=val`, `--check`, `--list-nodes`.
  fairchild-osdi/     OSDI v0.4 runtime: loads Verilog-A compiled by OpenVAF,
                      electrical (BSIM et al.) and optical alike.
  fairchild-py/       PyO3 Python package: Circuit / SimResult / WaveformSource.
examples/             Ready-to-run SPICE netlists per discipline (electronic, photonic).
benchmarks/           Head-to-head comparison circuits + scripts vs ngspice.
                      `run_all.py` → JSON, `plot.py` → docs/plots/*.png.
docs/                 user-guide.md, benchmarks.md, pn_phase_shifter_tiers.md,
                      docs/plots/ (generated accuracy + scaling figures).
scripts/              kicad_to_fairchild.py (KiCad netlist post-processor;
                      native fc_* devices).
legacy/               Superseded pre-Phase-B Verilog-A models. Historical only —
                      see legacy/README.md before reusing any of them.
```

---

## Roadmap

The major work ahead, in rough order:

1. **BJT CJE/CJC junction capacitances** — not yet stamped; required for
   correct BJT switching edge timing. (MOSFET Meyer gate caps and Cbs/Cbd are
   fully implemented; add CGSO/CGDO/CJ/CJSW to model cards to activate them.)
2. **Real-netlist test corpus on CI** — drop a foundry opamp and a published
   EO transceiver into the regression suite. Every failure becomes a Tier-0
   backlog item.
3. **Remaining analog elements** — lossy transmission lines (LTRA; lossless
   `T` and `S`/`W` switches already supported).
4. **Adjoint sensitivity** (the original Phase 4 differentiator).
5. **Tier-2 moats**: envelope-following, S-parameter Touchstone blocks with
   time-domain convolution, harmonic balance / PSS, WDM cross-channel
   nonlinearity.

---

## License

MIT
