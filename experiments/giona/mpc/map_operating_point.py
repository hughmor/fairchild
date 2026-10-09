#!/usr/bin/env python3
"""map_operating_point.py — where the network can actually be run.

Three things fight each other and only a map shows the winner:

  loop gain      G rises with laser power (kappa = R*a*P0), and the MPC weights
                 need G >= 1.18.
  self-heating   the ring absorbs some of that light and detunes off the
                 resonance the modulation depends on. At 30 mW/channel the ring
                 is already 4.9 K hot with the heater at zero, which is two
                 linewidths.
  estimator      G is only a gain while the probe weight keeps w*G well under 1.
                 Past that the self-coupled neuron is at its own bistability
                 threshold and the number means nothing.

    .venv/bin/python experiments/giona/mpc/map_operating_point.py

Writes results/operating_map.png and results/operating_map.json — the figure is
gitignored, the numbers are not, so the JSON is the durable record.
"""
from __future__ import annotations

import sys
from pathlib import Path

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from characterize import circuit, gain_validated

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
POWERS = (1.0, 3.0, 10.0, 30.0)
IHTS = np.linspace(0.0, 4.0e-3, 6)
PDB = -8.0
NEED = 1.18


def main() -> int:
    RESULTS.mkdir(exist_ok=True)
    G = np.full((len(POWERS), len(IHTS)), np.nan)
    BAD = np.zeros_like(G, dtype=bool)
    DT = np.full_like(G, np.nan)
    for a, p in enumerate(POWERS):
        c = circuit(p)
        for i in range(1, 7):
            c.set_param(f"Vpdb{i}", "dc", PDB)
            c.set_param(f"Iht{i}", "dc", 0.0)
        for b, v in enumerate(IHTS):
            c.set_param("Iht1", "dc", float(v))
            g, dis = gain_validated(c, 1, PDB, float(v))
            G[a, b], BAD[a, b] = g, (not np.isfinite(g)) or dis > 0.15
            DT[a, b] = float(c.run("op")["V(trm1)"][0])
        print(f"{p:5.1f} mW | G " + " ".join(f"{v:6.2f}" for v in G[a])
              + " | dT " + " ".join(f"{v:6.2f}" for v in DT[a]))

    (RESULTS / "operating_map.json").write_text(json.dumps({
        "pdb_V": PDB, "need_G": NEED,
        "powers_mW": list(POWERS), "iht_A": [float(v) for v in IHTS],
        "G": [[None if not np.isfinite(v) else float(v) for v in row] for row in G],
        "probe_out_of_range": BAD.tolist(),
        "dT_K": [[float(v) for v in row] for row in DT],
    }, indent=2) + "\n")
    print(f"wrote {RESULTS / 'operating_map.json'}")

    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.4))
    for a, p in enumerate(POWERS):
        col = plt.get_cmap("viridis")(a / (len(POWERS) - 1))
        ok = ~BAD[a]
        ax[0].plot(IHTS[ok] * 1e3, G[a][ok], "o-", color=col, label=f"{p:.0f} mW")
        ax[0].plot(IHTS[~ok] * 1e3, np.nan_to_num(G[a][~ok]), "x", color=col, ms=9)
        ax[1].plot(IHTS * 1e3, DT[a], "o-", color=col, label=f"{p:.0f} mW")
        ax[2].plot(DT[a], G[a], "o-", color=col, label=f"{p:.0f} mW")
    ax[0].axhline(NEED, color="k", ls="--", lw=1)
    ax[0].annotate("MPC needs 1.18", (0.05, NEED), fontsize=8,
                   textcoords="offset points", xytext=(2, 4))
    ax[0].set_xlabel("ring heater (mA)"); ax[0].set_ylabel("loop gain G")
    ax[0].set_title("gain vs trim\n(x = probe out of range, not a gain)")
    ax[1].set_xlabel("ring heater (mA)"); ax[1].set_ylabel("ring dT (K)")
    ax[1].axhline(0.37, color="k", ls=":", lw=1)
    ax[1].annotate("one linewidth", (0.05, 0.37), fontsize=8,
                   textcoords="offset points", xytext=(2, 4))
    ax[1].set_title("the price of that gain:\nself-heating detunes the ring")
    ax[2].axhline(NEED, color="k", ls="--", lw=1)
    ax[2].set_xlabel("ring dT (K)"); ax[2].set_ylabel("loop gain G")
    ax[2].set_xscale("log")
    ax[2].set_title("what it costs in detuning\nto buy a given gain")
    for a_ in ax:
        a_.legend(fontsize=8); a_.grid(alpha=0.25)
    fig.suptitle("giona MPC operating point — gain against the self-heating it costs",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = RESULTS / "operating_map.png"
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
