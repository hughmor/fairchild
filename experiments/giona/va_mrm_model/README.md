# va_mrm_model — the Verilog-A add-drop modulator

[`examples/verilog_a/models/mrm_addrop.va`](../../../examples/verilog_a/models/mrm_addrop.va)
is the eight modulators. This directory is how it earned its defaults.

| script | asks |
|---|---|
| `compare_va_mrm.py` | does it reproduce `examples/photonic/pcells/mrm.sp`, the native cell it replaces? |
| `va_vs_may_data.py` | does it reproduce the May capture the cell's card was fitted to? |

```bash
.venv/bin/python experiments/giona/va_mrm_model/compare_va_mrm.py
.venv/bin/python experiments/giona/va_mrm_model/va_vs_may_data.py --fit
```

Two results worth not rediscovering:

- In **card-compatible mode** (eleven overrides, listed in the model header) the
  Verilog-A ring is the native cell to 0.00 pm across reverse bias to −4 V,
  forward to +0.9 V and both notch flanks. That mode is also a Rust regression
  test, so the compatibility path is tested rather than claimed.
- On its own defaults it is **not** the cell, deliberately: it counts carriers
  once and runs depletion through Soref-Bennett, so the response is √-shaped
  where the card is linear. The two are 0.1 pm apart at −1 V, where the capture
  lives, and 23.5 pm apart at −4 V, where it does not.

`--fit` re-derives `vol_dep`/`vol_inj`/`r_th` from the capture and **does not
beat** the values derived from the waveguide volume, Soref-Bennett and the
heater slope. That is a statement about the capture, not the fitter: over its
1 V window there is 13.5 pm of shift against a 1 pm floor.
