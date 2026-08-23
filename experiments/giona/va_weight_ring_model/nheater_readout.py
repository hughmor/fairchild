#!/usr/bin/env python3
"""nheater_readout.py — the weight ring reading its own resonance out.

`examples/verilog_a/models/ring_nheater.va` has no photodiode. It has a ring
whose waveguide is doped n-type so it can be used as a resistor, and light lands
in that resistance twice, with opposite signs: absorbed photons make carriers,
which conduct and pull R DOWN, and what everything dissipates warms the silicon,
whose positive tempco pushes R up. The carriers win by nearly two orders of
magnitude, so the resonance shows up as a DIP.

Driven in CURRENT mode, deliberately. The self-heating closes a loop around the
readout and its sign depends on the drive: at constant current a falling R means
less power means cooler means R falls further, so the signal is amplified
(1.6x at 1.25 mA); at constant voltage the same fall means more power and the
loop fights it. Current mode also runs away above 2.73 mA on these defaults,
which is presumably why the published sweep stops at 1.25.

    .venv/bin/python experiments/giona/va_weight_ring_model/nheater_readout.py

Writes results/nheater_readout.png.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fairchild as fc

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
VA = HERE.parents[2] / "examples" / "verilog_a" / "models" / "ring_nheater.va"

# The model's own n_eff already puts the cold resonance about a third of a
# nanometre blue of 1550, so the sweep crosses it near the middle. Nothing to
# trim here.
#
# Worth knowing before reading the contrast: tuning is 0.251 nm/mW and 1.25 mA
# reaches 4.3 mW, so the ENTIRE current range is worth about 1.1 nm — roughly
# four linewidths of this 257 pm ring. The sweep still never reaches true
# anti-resonance, so the ratio below is set as much by how far the heater can
# go as by the ring.
DECK = f""".va {VA}
.optical_port src
.optical_port th
.optical_port ad
.optical_port dr
XL src fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xr src th ad dr hp 0 tr ring_nheater
IH 0 hp DC 1e-6
.op
"""
# 0.1 uW stands in for dark: the model has no zero-power path, and it absorbs
# ten thousand times less than the 1 mW trace.
P_DARK_MW = 1e-4
POWERS_MW = (0.25, 0.5, 1.0, 2.0)
I_MA = np.linspace(0.02, 1.25, 110)


def sweep(ckt, p_mw: float):
    """Ramp the heater CURRENT; read the voltage it develops."""
    ckt.set_param("XL", "power_mW", float(p_mw))
    v_h, r_ohm, drop, dT = (np.empty(len(I_MA)) for _ in range(4))
    for i, i_ma in enumerate(I_MA):
        ckt.set_param("IH", "dc", float(i_ma * 1e-3))
        r = ckt.run("op")
        v_h[i] = float(r["V(hp)"][0])
        r_ohm[i] = v_h[i] / (i_ma * 1e-3)
        dT[i] = float(r["V(tr)"][0])
        drop[i] = (float(r["V(dr_re_0)"][0]) ** 2
                   + float(r["V(dr_im_0)"][0]) ** 2)
    return r_ohm, drop, v_h, dT


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    ckt = fc.Circuit()
    ckt.load_str(DECK)

    dark, _, v_dark, _ = sweep(ckt, P_DARK_MW)
    runs = {p: sweep(ckt, p) for p in POWERS_MW}

    print(f"dark sweep: R {dark[0]:.0f} -> {dark[-1]:.0f} ohm, "
          f"V {v_dark[0]:.3f} -> {v_dark[-1]:.3f} V over "
          f"{I_MA[0]:.2f}-{I_MA[-1]:.2f} mA  (paper: 2200 -> ~2800, just over 3 V)")
    print(f"\n{'P_in mW':>8} {'I_res mA':>9} {'dR_res':>8} {'dR_end':>8}"
          f" {'on/off':>7} {'dV_res mV':>10}")
    for p, (r_ohm, drop, v_h, _) in runs.items():
        k = int(np.argmax(drop))
        dr = r_ohm - dark
        print(f"{p:8.2f} {I_MA[k]:9.3f} {dr[k]:8.1f} {dr[-1]:8.1f}"
              f" {dr.min() / dr[-1]:7.1f} {(v_h[k] - v_dark[k]) * 1e3:10.2f}")
    print("paper: light shifts the whole curve DOWN, ~20 ohm off resonance and")
    print("       ~300 ohm on it — a ratio of about 15")
    print("All three land once the coupler is set by the paper's Q of 5900 rather")
    print("than by the FEM gap ratio. It had to be the coupler: this ring's")
    print("linewidth is 31:1 coupling-dominated, so no carrier concentration")
    print("could have got there, and a broader ring both weakens the contrast and")
    print("costs more heat per linewidth of tuning.")

    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6))
    cmap = plt.get_cmap("viridis")
    for n, (p, (r_ohm, drop, _, _)) in enumerate(runs.items()):
        c = cmap(n / max(len(runs) - 1, 1))
        ax[0].plot(I_MA, drop / (p * 1e-3), color=c, label=f"{p:.2f} mW")
        ax[1].plot(I_MA, r_ohm, color=c, label=f"{p:.2f} mW")
        ax[2].plot(I_MA, r_ohm - dark, color=c, label=f"{p:.2f} mW")
    ax[1].plot(I_MA, dark, "k--", lw=1.4, label="dark")
    ax[2].axhline(0, color="k", lw=0.8)

    ax[0].set_ylabel("drop / input")
    ax[0].set_title("the optics: the heater walks the ring\nthrough the laser line")
    ax[1].set_ylabel("heater resistance (ohm)")
    ax[1].set_title("the electrics: self-heating raises R,\nand the light pulls it back down")
    ax[2].set_ylabel("R − R_dark  (ohm)")
    ax[2].set_title("the readout: carriers beat heat,\nso the resonance is a DIP")
    for a in ax:
        a.set_xlabel("heater current (mA)")
        a.legend(fontsize=8)
        a.grid(alpha=0.25)
    fig.suptitle("ring_nheater in current mode — absorbed photons make carriers that "
                 "conduct, so the ring reads its own resonance out as a drop in R",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = RESULTS / "nheater_readout.png"
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
