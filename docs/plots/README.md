# Generated figures

Committed because the README and the guides reference them, and a reader
should not have to run anything to see what the simulator does. Everything
here is reproducible from the repository — nothing is hand-drawn or touched up.

Set `MPLBACKEND=Agg` unless you want windows.

## Benchmarks

Both come from one script, which runs every circuit in `benchmarks/circuits/`
against ngspice and plots the comparison:

```bash
cargo build --release
python3 benchmarks/plot.py
```

| File | What |
|---|---|
| `accuracy_analog.png` | fairchild vs ngspice waveform overlay, with a residual strip under each panel |
| `scaling_wall_time.png` | transient wall-clock vs circuit size, each solver backend forced in turn |

## Photonic examples

Each is the figure its example script writes beside itself; copy it here after
regenerating. Every one of these scripts also takes `--selftest`, which asserts
the physics rather than plotting it.

```bash
python3 examples/photonic/receiver_noise_budget.py
python3 examples/photonic/native_mrr_wavelength_sweep.py
python3 examples/photonic/native_weight_bank.py
cp examples/photonic/{receiver_noise_budget,native_mrr_wavelength_sweep,native_weight_bank}.png docs/plots/

# The eye example twice, once per receiver front end. `--tia` needs
# examples/verilog_a/build/va_tia.osdi, so build the Verilog-A models first.
python3 examples/photonic/noisy_eye_and_ber.py --tia --png eye_tia.png
python3 examples/photonic/noisy_eye_and_ber.py       --png eye_res.png
cp examples/photonic/eye_tia.png docs/plots/noisy_eye_and_ber.png
cp examples/photonic/eye_res.png docs/plots/noisy_eye_rin_limited.png
```

| File | What |
|---|---|
| `receiver_noise_budget.png` | thermal / shot / RIN crossovers and the SNR ceiling at `1/(RIN·B)` |
| `noisy_eye_and_ber.png` | The headline figure: NRZ and PAM-4 eyes through an MZM built from primitives, read by the Verilog-A TIA. Noise at true amplitude — both rails equal, because the amplifier dominates |
| `noisy_eye_rin_limited.png` | The same link through a load resistor instead. Noise piles onto the `1` rail, which is what RIN-limited looks like and why the rails must be measured separately |
| `native_mrr_wavelength_sweep.png` | micro-ring through-port transmission, resonance shifting under bias |
| `native_weight_bank.png` | a 4-channel WDM weight bank: per-channel weights, passivity, balanced readout |

`examples/*/*.png` is gitignored, which is why these are copies rather than
symlinks — the examples directory stays clean for people running them.
