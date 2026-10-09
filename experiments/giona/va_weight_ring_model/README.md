# va_weight_ring_model — the Verilog-A N-doped-heater weight ring

[`examples/verilog_a/models/ring_nheater.va`](../../../examples/verilog_a/models/ring_nheater.va)
is the 64 weight rings: an add-drop ring with no junction, tuned by doping its
own waveguide n-type and running current through it.

```bash
.venv/bin/python experiments/giona/va_weight_ring_model/nheater_readout.py
```

`nheater_readout.py` drives the one loop worth building the model for. Light
lands in the resistance twice and with opposite signs — absorbed photons make
carriers, which conduct and pull R down; what everything dissipates warms
silicon whose tempco pushes R up. The carriers win, so ramping the heater
current with a laser on the bus finds the resonance as a **dip in resistance**,
with no detector in the circuit.

Driven in current mode on purpose: the self-heating loop amplifies the readout
at constant current and opposes it at constant voltage, and current drive runs
away above 2.73 mA on these defaults.

Matched to a published characterisation of this kind of device — 134 Ω/mW
self-heating, 0.25 nm/mW tuning, Q 5900, 12.1 nm FSR — which between them fix
`kappa_l` and `n_g` and leave the carriers no say: the linewidth is 31:1
coupling-dominated. Everything not so pinned is flagged PLACEHOLDER at its
declaration, because there is no capture of the giona weight rings to fit to.
