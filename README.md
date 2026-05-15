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
sparse LU, junction limiting). Photonic discipline rebuilt around native Rust
devices and bundle-port syntax. Python bindings cover every analysis the CLI
does. See [`sotu.md`](sotu.md) for the project state-of-the-union and
[`PLAN.md`](PLAN.md) for the architectural plan.

---

## What works today

### SPICE compatibility

| Category | Coverage |
|---|---|
| Elements | R, L, C, V, I, D, MOSFET (Level 1), B (behavioral), X (subckt / OSDI) |
| Sources | DC, PULSE, PWL, SIN, EXP, SFFM, AM |
| Analyses | `.op`, `.dc`, `.tran`, `.ac`, `.noise` |
| Solvers | NR with `pnjlim` / `fetlim`; BE, TR, GEAR (BDF-2); dense or sparse LU |
| Directives | `.options`, `.ic`, `.nodeset`, `.measure`, `.lib`/`.endl`, `.include`, `.param`, `.subckt`/`.ends`, `.temp` (sweep), `.alter`, `.model`, `.osdi` |
| Output | CSV (stdout / file), Nutmeg rawfile (ngspice-compatible) |

What's not yet supported: BJT, coupled inductors `K`, switches `S`/`W`,
transmission lines, `.disto`, `.pz`, native `.mc` Monte Carlo, PSF/FSDB binary
output. See [`sotu.md`](sotu.md) §3 for the live status list.

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
| Y-splitter | `fc_splitter` |
| PN-junction phase shifter | `fc_pn_ps` |
| Thermo-optic phase shifter | `fc_thermal_ps` |
| Photodetector | `fc_photodetector` |

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
│   └── ring_oscillator.sp        ← 5-stage CMOS ring; tests DC + GEAR
└── photonic/
    ├── native_mrr_modulator.{sp,py}        ← electro-optic micro-ring
    ├── native_mrr_wavelength_sweep.py      ← parametric λ sweep
    ├── native_wdm_mrr_modulator.{sp,py}    ← 2-channel WDM through one ring
    └── legacy/                              ← pre-Phase-B OSDI examples
```

Photonic examples come with a `README.md` describing the topology of each
circuit. The WDM example shows two lasers detuned ±50 pm sharing one ring;
the same V_pn pulse produces very different transmission profiles per
channel — no MUX/DEMUX device needed because the bundle-port abstraction
makes each wavelength its own independent SVEA channel.

---

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
  fairchild-osdi/     OSDI v0.4 runtime (compatibility shim — see crate docstring
                      for deprecation rationale).
  fairchild-py/       PyO3 Python package: Circuit / SimResult / WaveformSource.
examples/             Ready-to-run SPICE netlists + Python driver scripts.
docs/                 user-guide.md, photonic_models.md (legacy OSDI catalog),
                      generated comparison plots.
scripts/              kicad_to_fairchild.py (KiCad netlist post-processor —
                      due for rewrite around native devices).
va-models/            Legacy Verilog-A photonic models + build scripts.
```

---

## Roadmap

See [`sotu.md`](sotu.md) for the live status list.

The major work ahead, in rough order:

1. **Real-netlist test corpus on CI** — drop a foundry opamp and a published
   EO transceiver into the regression suite. Every failure becomes a Tier-0
   backlog item.
2. **Remaining analog elements** — BJT (Gummel-Poon), `K` coupled inductors,
   switches, transmission lines.
3. **KiCad integration around native `fc_*` devices** — symbol library
   rewrite, post-processor rewrite, action-plugin prototype. The current
   `scripts/kicad_to_fairchild.py` and `kicad_integration.md` are built
   around OSDI and need redoing.
4. **Adjoint sensitivity** (the original Phase 4 differentiator).
5. **Tier-2 moats**: envelope-following, S-parameter Touchstone blocks with
   time-domain convolution, harmonic balance / PSS, WDM cross-channel
   nonlinearity.

---

## License

MIT
