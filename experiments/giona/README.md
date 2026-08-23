# giona — PN/thermal micro-ring characterisation and the photonic RNN

Model-fitting and network work against a specific silicon-photonic neuron chip
("giona"): eight add-drop micro-ring **modulators** on a shared WDM bus feeding
balanced photodetectors, and 64 add-drop **weight rings** — eight banks of
eight — that close the recurrence.

This is not an example. For those see [`examples/`](../../examples). Everything
here assumes the chip's topology, its PCB network, and the dataset layout.

## How this directory is organised

One directory per experiment. Inside each, `netlists/` and `results/` hold that
experiment's decks and outputs and nothing else, so a reader can take any one
directory and know that everything in it belongs together.

| directory | what it is |
|---|---|
| [`common/`](common/) | the shared library, the shared driver, and the device card every experiment reads |
| [`pn_modulator_fit/`](pn_modulator_fit/) | fitting the `fc_pn_th_ps` card to the measured modulator |
| [`va_mrm_model/`](va_mrm_model/) | the Verilog-A add-drop modulator, against the native cell and against the capture |
| [`va_weight_ring_model/`](va_weight_ring_model/) | the Verilog-A N-doped-heater weight ring |
| [`mod_bank/`](mod_bank/) | the 8-ring cascade spectrum |
| [`frontend/`](frontend/) | the chip front-end deck generator |
| [`rnn_characterization/`](rnn_characterization/) | activation, rest point, loop gain — the physics-to-recurrence map |
| [`wta/`](wta/) | winner-take-all |
| [`mpc/`](mpc/) | model-predictive control (in progress) |

Two directories are shared rather than owned, because the data and the chip
project are inputs to everything:

- `data/` — raw lightlab `NdSweeper` pickles and extraction caches. Gitignored,
  ~13 GB.
- `kicad/` — the KiCad working project. Gitignored. Its `sym-lib-table` resolves
  the shared symbol library out of `examples/kicad_photonics/`, so there is only
  one copy of it.

`results/*.json` is committed; `results/*.png` and every `netlists/` are not.

## Running anything

```bash
source ~/.zshrc                                                # FAIRCHILD_OPENVAF, openvaf-r
maturin develop --release -m crates/fairchild-py/Cargo.toml    # once, and after any Rust change
.venv/bin/python experiments/giona/<experiment>/<script>.py
```

Scripts import the shared library by inserting `../common` on `sys.path`; run
them from anywhere.

## Known caveats

Carried forward into the model card's header, repeated here because they bound
what every fit downstream means.

1. **The diode trio is assembly-effective.** `i_sat`/`n_diode`/`r_series` were
   fitted through the PCB network (10 kΩ series, 2 kΩ shunt) and `n_diode` sits
   at its 5.0 bound. `dn_di`/`da_di` are only meaningful *paired with* that
   I(V) — refit both together from a clean on-die IV.
2. **The July captures have a broken junction path.** All three (neuron 3/7/8)
   show a dead-linear ~15.9 kΩ IV matching to 0.2 % across devices, symmetric
   through 0 V, with working heaters. Only their passive + thermal + linear-EO
   content is trusted.
3. **κL/α are degenerate** under-vs-over-coupled without an independent loss or
   drop-port measurement.
4. **Absolute `n_eff` must be trimmed per ring** — fab variation exceeds
   anything a shared card can carry.
