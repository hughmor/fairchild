# rnn_characterization — physics to recurrence

The map from device coefficients to the recurrence the network actually
implements, and the measurements behind every term in it.

**[`rnn_math.md`](rnn_math.md) is the document.** It derives the fixed point,
the bias row-sum rule, the small-signal division at the modulator node, the
junction-current-to-optical response, the loop gain `G`, and the time constants
— with every quoted number traced to the script that produced it.

| script | what |
|---|---|
| `rnn_explore.py` | activation sweeps, rest-point trade-off, in-situ gain (`--space`, `--insitu`, `--radii`) |
| `rnn_plots.py` | the figures: activation, rest point, dynamics across the oscillation threshold |

```bash
.venv/bin/python experiments/giona/rnn_characterization/rnn_explore.py --space
.venv/bin/python experiments/giona/rnn_characterization/rnn_plots.py
```

The driver itself is shared: [`../common/rnn_drive.py`](../common/rnn_drive.py).

Two facts from `rnn_math.md` §9 that bound everything built on top:

- Trimmed, `|G|` is 1.8–3.7 at 30 mW/channel. **Untrimmed it is 0.24–1.0, with
  the opposite sign.**
- `netlists/giona_rnn_perfectW.sp` hardcodes `Iht<k> ... DC 0`, so every run to
  date has been untrimmed. Trimming the heaters is the first thing to change
  before reading anything into those transients.
