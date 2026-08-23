#!/usr/bin/env python3
"""trim.py — set each ring's heater so every neuron drives the bus equally.

`gain_matrix.py` shows where the spread lives. Every ROW of G_ij is identical —
the six neuron circuits are the same, so the sink side is uniform — and every
COLUMN differs: at one shared heater current the six rings drive the bus at
0.003, 0.046, -0.424, 1.425, 1.954 and 0.018. Three near zero and one with the
wrong sign, because a shared `Iht` puts rings trimmed to six different
wavelengths at six different places on their own resonances.

So the trim is a per-ring, one-dimensional problem on the SOURCE side, and it is
cheap: driving neuron j and reading any other neuron's response gives ring j's
source gain in two solves. Sweep `Iht_j`, read its column, pick the current that
lands on a target every ring can reach.

The rings still interact — ring j sits on a bus that rings 1..j-1 have filtered,
and its own heat is its own — so this iterates until the assignment stops
moving. But each measurement is well conditioned, which the closed-loop
estimator this replaces was not.

    .venv/bin/python experiments/giona/mpc/trim.py

Writes results/operating_point.json (per-ring heater and the gain it buys) and
results/trim.png.
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
from gain_matrix import D_PDB, W_PROBE, mc

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
IHTS = np.linspace(0.0, 5.0e-3, 13)
PASSES = 3
NEED = 1.18


def source_gain(c, j: int, pdb: float) -> float:
    """Ring j's drive strength: how much every other neuron moves when j does.

    Two solves. Column j of W carries the probe so all six banks watch channel
    j; nothing feeds back because every other column is zero. The responding
    neurons all agree (that is the row-uniformity in `gain_matrix.png`), so the
    mean over them is the measurement and their spread is its error bar.
    """
    for a in range(1, N_NEURON + 1):
        for b in range(1, N_NEURON + 1):
            c.set_param(f"VW{a}{b}", "dc", W_PROBE if b == j else 0.0)
    c.set_param(f"Vpdb{j}", "dc", pdb + D_PDB)
    hi = mc(c)
    c.set_param(f"Vpdb{j}", "dc", pdb - D_PDB)
    lo = mc(c)
    c.set_param(f"Vpdb{j}", "dc", pdb)
    resp = (hi - lo) / (2 * D_PDB)
    g = resp / resp[j - 1] / W_PROBE
    return float(np.mean(np.delete(g, j - 1)))


def main() -> int:
    RESULTS.mkdir(exist_ok=True)
    op = OPERATING_POINT
    pdb = op["pdb_V"]
    c = circuit(op["p_las_mW"])
    # Start every ring where the gain matrix was measured, NOT at zero heater.
    # The rings share a bus, so ring j's drive strength depends on how much light
    # the others let past: with all six at zero they sit near resonance, absorb
    # the bus, and every source gain measures ~0.01 regardless of j's own trim.
    # The sweep has to happen around a configuration that already passes light.
    for i in range(1, N_NEURON + 1):
        c.set_param(f"Vpdb{i}", "dc", pdb)
        c.set_param(f"Iht{i}", "dc", op["iht_A"])
    iht = np.full(N_NEURON, op["iht_A"])

    print(f"source gain vs Iht, {op['p_las_mW']} mW/ch, PDB {pdb} V")
    print(f"{'ring':>5} | " + " ".join(f"{v * 1e3:5.2f}" for v in IHTS) + "  <- Iht mA")
    curves = {}
    for j in range(1, N_NEURON + 1):
        row = []
        for v in IHTS:
            c.set_param(f"Iht{j}", "dc", float(v))
            row.append(source_gain(c, j, pdb))
        c.set_param(f"Iht{j}", "dc", 0.0)
        curves[j] = np.array(row)
        print(f"{j:5d} | " + " ".join(f"{g:5.2f}" for g in row))

    reach = np.array([np.nanmax(curves[j]) for j in curves])
    target = float(np.nanmin(reach))
    print(f"\nhighest source gain every ring can reach: {target:.3f} "
          f"(ring {int(np.nanargmin(reach)) + 1} is the limit)")
    print(f"MPC needs 1.18: {'OK' if target >= NEED else 'SHORT'}, "
          f"margin {target / NEED:.2f}x")

    for p in range(PASSES):
        moved = 0.0
        for j in range(1, N_NEURON + 1):
            errs = []
            for v in IHTS:
                c.set_param(f"Iht{j}", "dc", float(v))
                errs.append(abs(source_gain(c, j, pdb) - target))
            k = int(np.nanargmin(errs))
            moved = max(moved, abs(IHTS[k] - iht[j - 1]))
            iht[j - 1] = float(IHTS[k])
            c.set_param(f"Iht{j}", "dc", iht[j - 1])
        gs = np.array([source_gain(c, j, pdb) for j in range(1, N_NEURON + 1)])
        print(f"pass {p + 1}: Iht = " + " ".join(f"{v * 1e3:4.2f}" for v in iht)
              + "  ->  g = " + " ".join(f"{g:5.2f}" for g in gs)
              + f"   spread {np.nanmax(gs) - np.nanmin(gs):.3f}")
        if moved < 1e-9:
            break

    gs = np.array([source_gain(c, j, pdb) for j in range(1, N_NEURON + 1)])
    converged = bool(target >= NEED)
    if not converged:
        print("\n*** THIS TRIM DID NOT SUCCEED. Do not feed it downstream. ***")
        print("The source gains come back near zero across the whole heater range,")
        print("and at the SAME configuration `gain_matrix.py` reads 1.4 and 2.0 for")
        print("rings 4 and 5. Two measurements of one quantity disagreeing by two")
        print("orders of magnitude means the measurement is wrong, not the trim.")
        print("Known suspects, none yet eliminated:")
        print("  * source_gain() puts the probe on the driven neuron too (W_jj),")
        print("    so its reference response carries its own self-feedback;")
        print("  * the responses divided here are ~1 mV on a 977 mV node, which is")
        print("    the same reltol floor that broke the closed-loop estimator;")
        print("  * the rings couple hard through the shared bus, so a one-at-a-time")
        print("    sweep may not be measuring what it thinks it is.")
        print("Next: pin the measurement against a case with a KNOWN answer before")
        print("using it to set anything — e.g. an ideal 2x2 at a programmed weight,")
        print("where the expected response is arithmetic.")
    out = {"converged": converged,
           "p_las_mW": op["p_las_mW"], "pdb_V": pdb,
           "iht_A": [float(v) for v in iht], "G": [float(g) for g in gs],
           "G_target": target, "need": NEED,
           "iht_grid_A": [float(v) for v in IHTS],
           "curves": {str(j): [float(v) for v in curves[j]] for j in curves}}
    (RESULTS / "operating_point.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"\nfinal source gains {gs.min():.3f} to {gs.max():.3f} "
          f"(was 0.003 to 1.954 at one shared trim)")
    print(f"wrote {RESULTS / 'operating_point.json'}")

    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.4))
    cmap = plt.get_cmap("tab10")
    for j in curves:
        ax[0].plot(IHTS * 1e3, curves[j], "o-", ms=3, color=cmap(j - 1),
                   label=f"ring {j}")
        ax[0].plot(iht[j - 1] * 1e3, gs[j - 1], "*", ms=15, color=cmap(j - 1))
    ax[0].axhline(target, color="k", ls=":", lw=1)
    ax[0].axhline(NEED, color="k", ls="--", lw=1.2)
    ax[0].annotate("MPC needs 1.18", (0.1, NEED), fontsize=8,
                   textcoords="offset points", xytext=(2, 4))
    ax[0].set_xlabel("ring heater (mA)"); ax[0].set_ylabel("source gain")
    ax[0].set_title("each ring's drive strength vs its own trim\n(star = chosen)")
    ax[0].legend(fontsize=7, ncol=2); ax[0].grid(alpha=0.25)
    before = [0.003, 0.046, -0.424, 1.425, 1.954, 0.018]
    x = np.arange(1, N_NEURON + 1)
    ax[1].bar(x - 0.2, before, 0.4, color="#8c8c8c", label="one shared trim")
    ax[1].bar(x + 0.2, gs, 0.4, color="#4c72b0", label="per-ring trim")
    ax[1].axhline(NEED, color="k", ls="--", lw=1.2)
    ax[1].set_xlabel("ring"); ax[1].set_ylabel("source gain")
    ax[1].set_title("what the trim bought")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.25, axis="y")
    fig.suptitle("giona MPC — per-ring trim to equalise drive strength", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(RESULTS / "trim.png", dpi=130)
    print(f"wrote {RESULTS / 'trim.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
