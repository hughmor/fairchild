# wta — winner-take-all

```bash
.venv/bin/python experiments/giona/wta/rnn_wta.py
```

Drives the 8-neuron recurrence with a competitive weight matrix and watches one
channel win. Uses the shared driver, [`../common/rnn_drive.py`](../common/rnn_drive.py),
over `../rnn_characterization/netlists/giona_rnn_perfectW.sp`.

WTA needs `Re(G·λ) > 1` — it is an *instability*, so below threshold the
behaviour does not exist rather than degrading. That is what pushed these runs
to >30 mW/channel.

Worth re-testing before spending more optical power: those runs were
**untrimmed**, where `|G| ≈ 0.24–1.0` against 1.8–3.7 trimmed
(`../rnn_characterization/rnn_math.md` §9). The gain shortfall may have been a
trimming problem rather than a power one.
