# giona — PN/thermal micro-ring modulator characterisation

Model-fitting work against a specific silicon-photonic neuron chip ("giona"):
banks of cascaded add-drop micro-ring modulators (one ring per wavelength
channel, shared through + drop bus) feeding a balanced-PD neuron. Everything
here exists to produce **one deliverable**:

> **[`giona_pn_th_ps.inc`](giona_pn_th_ps.inc)** — a `.model … fc_pn_th_ps LEVEL=4`
> card whose defaults reproduce the measured device. Concatenate it ahead of a
> netlist (the parser has no `.include`).

This is not an example — for those see [`examples/`](../../examples). Fitting
scripts here assume the chip's topology, PCB network, and dataset layout.

## Fits

| script | dataset | fits | result |
|---|---|---|---|
| `ringfit.py` | May sparse joint IV+spectra (neuron2 capture → neuron7 params) | staged: passive → thermal → EO → injection | `results/giona_neuron7_pn_th_ps_full_fit.json` |
| `fit_may_injection.py` | ↑ same, forward-bias slice | `dn_di`, `da_di` with the diode pinned from the IV alone | `results/giona_dn_di_fit.json` |
| `fit_jul_neuron3.py` | Jul 2026 dense 50×100 HC×JV sweep (3.6 GB) | passive + thermal + linear EO from tracked notch positions | `results/giona_neuron3_pn_th_ps_fit.json` |
| `fit_transient.py` | (machinery + synthetic recovery selftest) | time-domain fitting vs an AWG-drive/scope-PD capture | — |
| `expt_forward_mismatch.py` | synthetic | model-**form** adequacy: linear vs full PN | 4.6× residual floor gap |
| `sweep_mod_bank.py` | — | 8-ring cascade spectrum from `netlists/giona_mod_bank_full.sp` | `results/giona_mod_bank_spectrum.png` |
| `build_frontend.py` | — | generates the chip front-end netlist (source bank → 8-ring bank → 2 mm bus → 1:8 log tree → 8 programmable 2×2 weight blocks) | `netlists/giona_frontend.sp`, `netlists/mrm_wdm8.sp` |

`ringfit.py` is also the shared library: dataset → observables
(`load_sweep`, `extract_data`), netlist assembly, ring wavelength sweep,
the staged fitter, and all plotting. The other scripts import from it.

```bash
maturin develop --release -m crates/fairchild-py/Cargo.toml   # once
.venv/bin/python experiments/giona/ringfit.py --model fc_pn_th_ps_full
.venv/bin/python experiments/giona/fit_may_injection.py
.venv/bin/python experiments/giona/fit_jul_neuron3.py --cache   # rebuild npz first
```

## Layout

- `data/` — raw lightlab `NdSweeper` pickles + extraction caches (`.npz`).
  Gitignored: ~13 GB.
- `results/` — fitted params (`.json`, committed) and plots (`.png`, ignored).
- `netlists/` — KiCad SPICE exports of the chip, hand-written testbenches, and
  the output of `build_frontend.py`. Gitignored (machine-generated, and the chip
  layout isn't ours to publish) — regenerate with
  `.venv/bin/python experiments/giona/build_frontend.py`.
- `kicad/` — the KiCad working project for the chip. Gitignored. Its
  `sym-lib-table` resolves the shared symbol library out of
  `examples/kicad_photonics/`, so there is only one copy of it.

## Known caveats

Carried forward into the model card's header, repeated here because they bound
what the fits mean:

1. **The diode trio is assembly-effective.** `i_sat`/`n_diode`/`r_series` were
   fitted through the PCB network (10 kΩ series, 2 kΩ shunt) and `n_diode`
   sits at its 5.0 bound. `dn_di`/`da_di` are only meaningful *paired with*
   that I(V) — refit both together from a clean on-die IV.
2. **The July captures have a broken junction path.** All three (neuron 3/7/8)
   show a dead-linear ~15.9 kΩ IV matching to 0.2% across devices, symmetric
   through 0 V, with working heaters. Only their passive + thermal + linear-EO
   content is trusted; they cannot speak to forward bias.
3. **κL/α are degenerate** under-vs-over-coupled without an independent loss
   or drop-port measurement.
4. **Absolute `n_eff` must be trimmed per ring** — fab variation exceeds
   anything a shared card can carry.
