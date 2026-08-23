# frontend — the chip front-end deck

Generates the front end: WDM source bank → 8-ring modulator bank → 2 mm bus →
1:1 tap → 1:8 log tree → 8 banks of 8 weight rings.

```bash
.venv/bin/python experiments/giona/frontend/build_frontend.py
```

Writes two decks into `netlists/`, identical except for the weight bank and
sharing every net name a script would read:

| deck | weights |
|---|---|
| `giona_frontend.sp` | 64 real `ring_nheater` rings. `VWij` is a heater **voltage** across 2.2 kΩ, not a weight — zero volts is full drop. |
| `giona_frontend_idealW.sp` | `fc_optical_2x2` blocks at `w=0 dw_dv=1`, so `V(Wij)` **is** the weight in [−1,+1]. For getting an experiment working before the weight-to-voltage map is. |

Both rings are bundle-aware, so one instance carries all eight channels and the
deck states the width. The hand-generated 8-channel subckt this directory used
to emit is gone, along with the 121-solve numerical `n_eff` trim it needed — the
dispersion is a straight line, so the comb formula carries it in closed form.

⚠️ `netlists/mrm_wdm8.sp` is a stale leftover of that generator, still
`.include`d by `rnn_characterization/netlists/giona_rnn_perfectW.sp`. That deck
needs porting onto the Verilog-A rings; until then, do not delete the leftover.
