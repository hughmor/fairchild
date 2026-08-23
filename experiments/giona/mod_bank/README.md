# mod_bank — the 8-ring cascade spectrum

```bash
.venv/bin/python experiments/giona/mod_bank/sweep_mod_bank.py
```

Wavelength sweep across `netlists/giona_mod_bank_full.sp`: eight modulator rings
on one bus, each trimmed to its own channel, read at the through and drop ports.
This is the check that the comb and the ring trims line up before anything tries
to use the bank as a weight or a neuron layer.

→ `results/giona_mod_bank_spectrum.png`
