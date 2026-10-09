#!/usr/bin/env python3
"""gain_matrix.py — the full neuron-to-neuron gain, measured open loop.

## Why not the closed-loop estimator

`rnn_math.md` (11) gets the loop gain from how a self-weight changes a neuron's
bias sensitivity, `G = (1 - S0/S(w))/w`. It is a ratio of two nearly equal
numbers divided by a small one, and on this deck it does not survive contact:
across a heater sweep it returned 1.64, 1.81, 4.86, 5.12 and once -4255. Two
separate causes, both real and both worth recording:

  * the probe weight has to keep `w*G` well under 1, or the self-coupled neuron
    is at its own bistability threshold and S(w) diverges and changes sign;
  * `reltol` is 1e-3 on a node at 977 mV — a millivolt of slop — and at a 20 mV
    bias step the whole response being differenced is 1.4 mV. Tightening the
    tolerance does not help either: at 1e-10 the solver chases the gmin-floored
    dark optical nets and stops converging.

Even at a 200 mV step and a small probe, four of six neurons still disagreed
between probe weights by more than 15 %.

## What this does instead

`rnn_math.md` §9(b) — the direct route, "no division by a small number". The
coupling is not one scalar anyway: neuron i's response to neuron j is
`G_ij = eta_i * T'_j * R * a`, a sink factor times a source factor, so what the
weights actually need is the whole matrix.

Measure a COLUMN at a time. Put the probe weight on column j only — every row,
so all six banks watch channel j — and perturb neuron j's bias. Nothing feeds
back, because every other column is zero, so each neuron's response is one hop:

    G_ij = (1/w) * d(mc_i) / d(mc_j)

A ratio of two directly-measured responses, no small difference anywhere, and
both neurons sit at the same bias so the mc-to-current conversion cancels. Two
solves per column, twelve for the matrix, against roughly two hundred for the
closed-loop sweep it replaces.

    .venv/bin/python experiments/giona/mpc/gain_matrix.py

Writes results/gain_matrix.json and results/gain_matrix.png.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from characterize import N_NEURON, OPERATING_POINT, circuit

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
W_PROBE = 0.4      # open loop, so this can be large and well resolved
D_PDB = 0.2        # 200 mV: the response has to clear reltol*V ~ 1 mV


def mc(c) -> np.ndarray:
    r = c.run("op")
    return np.array([float(r[f"V(mc{i})"][0]) for i in range(1, N_NEURON + 1)])


def gain_matrix(c, pdb: float, w: float = W_PROBE, d: float = D_PDB) -> np.ndarray:
    G = np.zeros((N_NEURON, N_NEURON))
    for j in range(1, N_NEURON + 1):
        for a in range(1, N_NEURON + 1):
            for b in range(1, N_NEURON + 1):
                c.set_param(f"VW{a}{b}", "dc", w if b == j else 0.0)
        c.set_param(f"Vpdb{j}", "dc", pdb + d)
        hi = mc(c)
        c.set_param(f"Vpdb{j}", "dc", pdb - d)
        lo = mc(c)
        c.set_param(f"Vpdb{j}", "dc", pdb)
        resp = (hi - lo) / (2 * d)
        G[:, j - 1] = resp / resp[j - 1] / w
        # The driven neuron itself has no weight in its own path here — its
        # response is the bias sensitivity, which is the reference, not a gain.
        G[j - 1, j - 1] = np.nan
    return G


def main() -> int:
    RESULTS.mkdir(exist_ok=True)
    op = OPERATING_POINT
    c = circuit(op["p_las_mW"])
    for i in range(1, N_NEURON + 1):
        c.set_param(f"Vpdb{i}", "dc", op["pdb_V"])
        c.set_param(f"Iht{i}", "dc", op["iht_A"])
    G = gain_matrix(c, op["pdb_V"])

    print(f"gain matrix G_ij at {op['p_las_mW']} mW/ch, PDB {op['pdb_V']} V, "
          f"Iht {op['iht_A'] * 1e3:.2f} mA")
    print("(row = responding neuron, column = driven neuron; diagonal not "
          "measurable this way)")
    print(np.array2string(G, precision=3, suppress_small=True))
    off = G[~np.isnan(G)]
    print(f"\noff-diagonal G: {off.min():.3f} to {off.max():.3f}, median {np.median(off):.3f}")
    # Is it rank 1, i.e. eta_i * T'_j? That is what a shared bus should give.
    Gf = np.where(np.isnan(G), np.nanmean(G), G)
    sv = np.linalg.svd(Gf, compute_uv=False)
    print(f"singular values {np.round(sv, 3)} -> rank-1 fraction "
          f"{sv[0] / sv.sum():.1%}")
    print("A rank-1 G means one sink factor per row and one source factor per")
    print("column, which is what a shared bus predicts, and it means the trim")
    print("has two independent knobs per neuron rather than one entangled one.")

    (RESULTS / "gain_matrix.json").write_text(json.dumps({
        **op, "w_probe": W_PROBE, "d_pdb": D_PDB,
        "G": [[None if np.isnan(v) else float(v) for v in row] for row in G],
        "singular_values": [float(v) for v in sv],
    }, indent=2) + "\n")
    print(f"\nwrote {RESULTS / 'gain_matrix.json'}")

    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    im = ax[0].imshow(G, cmap="viridis")
    ax[0].set_xticks(range(N_NEURON), [str(i) for i in range(1, N_NEURON + 1)])
    ax[0].set_yticks(range(N_NEURON), [str(i) for i in range(1, N_NEURON + 1)])
    ax[0].set_xlabel("driven neuron j"); ax[0].set_ylabel("responding neuron i")
    ax[0].set_title("gain matrix $G_{ij}$\n(diagonal not measurable open loop)")
    for i in range(N_NEURON):
        for j in range(N_NEURON):
            if not np.isnan(G[i, j]):
                ax[0].text(j, i, f"{G[i, j]:.2f}", ha="center", va="center",
                           fontsize=7, color="w")
    fig.colorbar(im, ax=ax[0])
    ax[1].bar(range(1, N_NEURON + 1), np.nanmean(G, axis=1), color="#4c72b0",
              label="row mean (sink factor, $\\eta_i$)")
    ax[1].bar(range(1, N_NEURON + 1), np.nanmean(G, axis=0), color="#dd8452",
              width=0.45, label="column mean (source factor, $T'_j$)")
    ax[1].axhline(1.18, color="k", ls="--", lw=1.2)
    ax[1].annotate("MPC needs 1.18", (0.6, 1.18), fontsize=9,
                   textcoords="offset points", xytext=(2, 4))
    ax[1].set_xlabel("neuron"); ax[1].set_ylabel("mean gain")
    ax[1].set_title("where the spread lives:\nresponding side vs driving side")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.25, axis="y")
    fig.suptitle("giona MPC — neuron-to-neuron gain, measured open loop", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(RESULTS / "gain_matrix.png", dpi=130)
    print(f"wrote {RESULTS / 'gain_matrix.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
