#!/usr/bin/env python3
"""weight_bank_spacing.py — one 8-ring weight bank, two bus spacings.

Eight `ring_nheater` rings share a thru bus and a drop bus. The thru port of
ring i feeds the in port of ring i+1. The drop port of ring i+1 feeds the add
port of ring i. So the drop bus runs backwards and exits at ring 0. Ring 0 has
a radius of 8 um, and each next ring is 10 nm larger. The heaters carry no
current.

Both buses use the same waveguide length between rings. The script compares
95 um against 100 um. The two buses close a loop through every pair of rings,
so the spacing sets the phase of that loop.

    .venv/bin/python experiments/giona/weight_bank_spacing/weight_bank_spacing.py

Add `--full` to sweep 1450-1600 nm. Writes results/weight_bank_spacing.png
and the raw spectra to results/weight_bank_spacing.npz.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fairchild as fc

HERE = Path(__file__).resolve().parent
VA = HERE.parents[2] / "examples" / "verilog_a" / "models" / "ring_nheater.va"

N_RINGS = 8
R0_M, DR_M = 8.0e-6, 13e-9
SPACINGS_UM = (95.0, 100.0)
# Low power keeps the rings passive: self-heating from absorbed light would
# shift each resonance and make the spectrum depend on the sweep direction.
P_LASER_MW = 0.01
# 10 pm resolves the 257 pm ring linewidth with about 25 points. The default
# 20 nm covers more than one 12 nm FSR in about 100 s. `--full` sweeps
# 1450-1600 nm, which takes about 13 minutes.
WL_LO, WL_HI = (1450.0, 1600.0) if "--full" in sys.argv else (1545.0, 1552.0)
WL_NM = np.arange(WL_LO, WL_HI + 1e-9, 0.02)


def deck(spacing_um: float) -> str:
    ports = ["src"] + [f"{p}{k}" for k in range(N_RINGS) for p in "itad"]
    lines = [
        "* 8-ring N-heater weight bank",
        f".va {VA}",
        *(f".optical_port {p}" for p in ports),
        f"XL src fc_cw_laser power_mW={P_LASER_MW} wavelength_nm=1550",
        "Xin src i0 fc_waveguide L_um=10",
    ]
    for k in range(N_RINGS):
        lines.append(f"Xr{k} i{k} t{k} a{k} d{k} 0 0 th{k} ring_nheater "
                     f"radius={R0_M + k * DR_M:.6e}")
    for k in range(N_RINGS - 1):
        lines.append(f"Xwt{k} t{k} i{k + 1} fc_waveguide L_um={spacing_um}")
        lines.append(f"Xwd{k} d{k + 1} a{k} fc_waveguide L_um={spacing_um}")
    # The bus loop is a cavity, but a DC spectrum needs no photon lifetime.
    lines += [".options optical_delay=0", ".op", ".end"]
    return "\n".join(lines) + "\n"


def power(r, node: str) -> float:
    return float(r[f"V({node}_re_0)"][0]) ** 2 + float(r[f"V({node}_im_0)"][0]) ** 2


def sweep(spacing_um: float) -> tuple[np.ndarray, np.ndarray]:
    ckt = fc.Circuit()
    ckt.load_str(deck(spacing_um))
    p_in = P_LASER_MW * 1e-3
    thru, drop = np.empty(len(WL_NM)), np.empty(len(WL_NM))
    t0 = time.perf_counter()
    for i, wl in enumerate(WL_NM):
        ckt.set_param("XL", "wavelength_nm", float(wl))
        r = ckt.run("op")
        thru[i] = power(r, f"t{N_RINGS - 1}") / p_in
        drop[i] = power(r, "d0") / p_in
    print(f"{spacing_um:g} um: {len(WL_NM)} points in {time.perf_counter() - t0:.1f} s")
    return thru, drop


def db(x: np.ndarray) -> np.ndarray:
    return 10.0 * np.log10(np.maximum(x, 1e-12))


def main() -> None:
    # The sweep takes minutes, so keep it on disk before plotting.
    cache = HERE / "results" / "weight_bank_spacing.npz"
    runs = {s: sweep(s) for s in SPACINGS_UM}
    np.savez(cache, wl_nm=WL_NM, **{f"{k}_{s:g}um": v for s, tv in runs.items()
                                    for k, v in zip(("thru", "drop"), tv)})

    fig, ax = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    for s, (thru, drop) in runs.items():
        ax[0].plot(WL_NM, db(thru), lw=0.7, label=f"{s:g} µm")
        ax[1].plot(WL_NM, db(drop), lw=0.7, label=f"{s:g} µm")
        ax[2].plot(WL_NM, drop - thru, lw=0.7, label=f"{s:g} µm")
    ax[0].set_ylabel(f"ring {N_RINGS - 1} thru (dB)")
    ax[1].set_ylabel("ring 0 drop (dB)")
    ax[2].set_ylabel("drop − thru (linear, of input)")
    ax[2].axhline(0, color="k", lw=0.6)
    ax[2].set_xlabel("wavelength (nm)")
    for a in ax:
        a.legend(fontsize=8, loc="lower right")
        a.grid(alpha=0.25)
    fig.suptitle("8-ring N-heater weight bank: bus spacing 95 µm against 100 µm")
    fig.tight_layout()
    out = HERE / "results" / "weight_bank_spacing.png"
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
