#!/usr/bin/env python3
"""nheater_readout.py — the weight ring reading its own resonance out.

`examples/verilog_a/models/ring_nheater.va` has no photodiode. It has a ring
whose waveguide is doped n-type so it can be used as a resistor, and doping a
waveguide puts free carriers in the optical mode, and free carriers absorb. What
they absorb becomes heat in the same silicon whose resistance is being measured.

So: put a laser on the bus, sweep the heater voltage, measure the CURRENT. The
heater walks the resonance across the laser line; when it crosses, the
circulating power spikes by the cavity enhancement, the extra absorbed power
warms the ring further, and the resistance rises. The resonance appears in an
electrical measurement with no detector in the circuit.

    .venv/bin/python experiments/giona/nheater_readout.py

Writes results/nheater_readout.png.
"""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault(
    "FAIRCHILD_OPENVAF",
    "/Users/hugh/Local/src/OpenVAF-Reloaded/target/release/openvaf-r",
)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fairchild as fc

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
VA = HERE.parents[1] / "examples" / "verilog_a" / "models" / "ring_nheater.va"

DECK = f""".va {VA}
.optical_port src
.optical_port th
.optical_port ad
.optical_port dr
XL src fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xr src th ad dr hp 0 tr ring_nheater
VH hp 0 DC 0
.op
"""
# 1 uW stands in for dark: the model has no zero-power path, and a microwatt
# absorbs a millionth of what a milliwatt does.
P_DARK_MW = 1e-3
POWERS_MW = (0.5, 1.0, 2.0, 4.0)
V = np.linspace(0.05, 2.4, 130)


def sweep(ckt, p_mw: float):
    ckt.set_param("XL", "power_mW", float(p_mw))
    r_ohm, drop, thru, dT = (np.empty(len(V)) for _ in range(4))
    for i, v in enumerate(V):
        ckt.set_param("VH", "dc", float(v))
        r = ckt.run("op")
        i_h = -float(r["I(vh)"][0])
        r_ohm[i] = v / i_h
        dT[i] = float(r["V(tr)"][0])
        for name, out in (("dr", drop), ("th", thru)):
            out[i] = (float(r[f"V({name}_re_0)"][0]) ** 2
                      + float(r[f"V({name}_im_0)"][0]) ** 2)
    return r_ohm, drop, thru, dT


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    ckt = fc.Circuit()
    ckt.load_str(DECK)

    dark, _, _, dT_dark = sweep(ckt, P_DARK_MW)
    runs = {p: sweep(ckt, p) for p in POWERS_MW}

    print(f"{'P_in mW':>8} {'V_res':>7} {'R_res':>9} {'excess ppm':>11}"
          f" {'dT opt K':>9} {'dI uA':>8}")
    for p, (r_ohm, drop, _, dT) in runs.items():
        k = int(np.argmax(drop))
        exc = (r_ohm - dark) / dark
        j = int(np.argmax(exc))
        i_lit, i_dark = V[j] / r_ohm[j], V[j] / dark[j]
        print(f"{p:8.2f} {V[k]:7.3f} {r_ohm[k]:9.2f} {exc[j] * 1e6:11.1f}"
              f" {dT[j] - dT_dark[j]:9.4f} {(i_lit - i_dark) * 1e6:8.2f}")

    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6))
    cmap = plt.get_cmap("viridis")
    for n, (p, (r_ohm, drop, thru, dT)) in enumerate(runs.items()):
        c = cmap(n / max(len(runs) - 1, 1))
        ax[0].plot(V, drop / (p * 1e-3), color=c, label=f"{p:.1f} mW")
        ax[1].plot(V, r_ohm, color=c, label=f"{p:.1f} mW")
        ax[2].plot(V, (r_ohm - dark) / dark * 1e6, color=c, label=f"{p:.1f} mW")
    ax[0].plot(V, runs[1.0][2] / 1e-3, "k:", lw=1, label="thru, 1 mW")
    ax[1].plot(V, dark, "k--", lw=1.2, label="dark")

    ax[0].set_ylabel("drop / input")
    ax[0].set_title("the optics: the heater tunes the ring\nthrough the laser line")
    ax[1].set_ylabel("heater resistance (ohm)")
    ax[1].set_title("the electrics: R rises with its own\ntemperature, light or no light")
    ax[2].set_ylabel("(R − R_dark) / R_dark  (ppm)")
    ax[2].set_title("the readout: what the light adds,\nand it peaks on resonance")
    for a in ax:
        a.set_xlabel("heater voltage (V)")
        a.legend(fontsize=8)
        a.grid(alpha=0.25)
    fig.suptitle("ring_nheater — the doping that makes the resistor also makes the "
                 "absorption, so the ring reads its own resonance out", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = RESULTS / "nheater_readout.png"
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
