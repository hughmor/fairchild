#!/usr/bin/env python3
"""trim.py — set each ring's heater so every neuron drives the bus equally.

## Two trims, and the deck only had one

A ring has two knobs and they do different jobs. `n_eff` is the DESIGN trim, set
at layout by the ring's radius; the heater is the OPERATING trim, and it can
only ever add heat. `build_mpc_deck.py` sets `n_eff` to put each ring's COLD
resonance on its channel — which is wrong at the power this network runs at.

At 30 mW/channel the ring self-heats 5.2 K before the heater does anything,
which is 360 pm, and its own notch is 167 pm wide. So every resonance has
already walked two linewidths RED of its channel, and the heater — which only
reddens — cannot bring it back. The coarse sweep duly reports every ring's notch
at the bottom edge of the heater range, because there is no notch inside it.

So this file trims `n_eff` first, to put the resonance on channel AT the
operating power with the heater mid-range, and only then trims the heater. On a
real chip that first step is a layout decision rather than a knob, which is
exactly why it has to be got right in the deck.

## The measurement this uses, and why the earlier ones did not work

`validate_estimator.py` checks the gain chain factor by factor against closed
forms. Three of the four factors are exact and all three reproduce to machine
precision: the weight block and balanced pair give `I = R*w*P_in` to 7e-16, the
tap and tree give exactly 1/16, and the node division gives R_sh = 1394 ohm and
eta = 0.591 at 134 uA, matching `rnn_math.md` (6). Only the ring's slope is
device physics, and it is a ratio of two directly-read numbers.

So the gain is assembled, not fitted:

    G_j = eta * R * dP_win_j / dI_D_j

with `dI_D` from the exact derivative of the model's own diode law. Nothing
divides a small difference of nearly-equal numbers, and **no weight is
involved** — which is what makes this cheap: two solves per point, and the
weight bank can stay at zero throughout.

That measurement says the rings are badly untrimmed at a shared heater current:
G = -0.007, 0.102, -0.002, 8.27, -0.002, -3.88. Four sitting at or near a slope
zero, and the two live ones with opposite signs. A ring at its notch bottom has
no gain and no defined sign, which is precisely `rnn_math.md` §9's warning.

    .venv/bin/python experiments/giona/mpc/trim.py

Writes results/operating_point.json and results/trim.png.
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
from validate_estimator import N_DIODE, R_PD, VT

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
# Coarse-then-fine, because the feature is far narrower than the knob's range.
# At 5 mA the heater tunes this ring about 2 nm while its notch is 167 pm wide,
# so a uniform 17-point sweep samples a sharp resonance every 0.7 linewidths and
# ALIASES: ring 3 came back with eight sign changes, which is not a resonance
# flank, it is undersampling. Find the notch first with a cheap one-solve power
# sweep, then measure the slope only where the slope exists.
IHT_COARSE = np.linspace(0.0, 5.0e-3, 41)   # 1 solve each, find the notch
FINE_SPAN = 0.6e-3                          # +-, around the notch
N_FINE = 15                                 # 2 solves each
PASSES = 3
NEED = 1.18
I_SAT = 5.099e-8
NVT = N_DIODE * VT
R_SH = 1.0 / (1 / 2e3 + 1 / 10e3 + 2 / 17.051e3)
D_PDB = 0.2


def source_gain(c, j: int, pdb: float) -> float:
    """G_j = eta * R * dP_win_j/dI_D_j. Two solves, no weights involved."""
    def read():
        r = c.run("op")
        return (float(r[f"V(mc{j})"][0]),
                float(r[f"V(win1_re_{j - 1})"][0]) ** 2
                + float(r[f"V(win1_im_{j - 1})"][0]) ** 2)

    c.set_param(f"Vpdb{j}", "dc", pdb + D_PDB); mch, pwh = read()
    c.set_param(f"Vpdb{j}", "dc", pdb - D_PDB); mcl, pwl = read()
    c.set_param(f"Vpdb{j}", "dc", pdb)
    mc0 = 0.5 * (mch + mcl)
    i_d = I_SAT * (np.exp(min(-mc0 / NVT, 700.0)) - 1.0)
    eta = R_SH / (R_SH + NVT / i_d)
    di_d = -(i_d + I_SAT) / NVT * (mch - mcl) / (2 * D_PDB)
    dpw = (pwh - pwl) / (2 * D_PDB)
    return float(eta * R_PD * dpw / di_d) if di_d else float("nan")


IHT_PARK = 1.5e-3       # where the heater should sit once n_eff is right


def park_n_eff(c, j: int, n_eff0: float) -> float:
    """Move ring j's DESIGN index until its notch sits at `IHT_PARK`.

    The ring is hot, so its resonance is red of where a cold trim put it; the
    fix is to start it blue by the same amount. Search `n_eff` directly — it is
    an instance parameter — and take the value whose channel is darkest at the
    parked heater current.
    """
    c.set_param(f"Iht{j}", "dc", IHT_PARK)
    # One linewidth is ~167 pm; a nanometre of search is six of them, which is
    # more than the self-heating can have moved it.
    cand = n_eff0 + np.linspace(-1.6e-3, 0.4e-3, 33)
    best, bp = n_eff0, np.inf
    for v in cand:
        c.set_param(f"Xr{j}", "n_eff", float(v))
        r = c.run("op")
        pw = (float(r[f"V(win1_re_{j - 1})"][0]) ** 2
              + float(r[f"V(win1_im_{j - 1})"][0]) ** 2)
        if pw < bp:
            best, bp = float(v), pw
    c.set_param(f"Xr{j}", "n_eff", best)
    return best


def notch_current(c, j: int) -> tuple[float, np.ndarray]:
    """Where ring j's own channel is darkest — one solve per point, no probe."""
    p = []
    for v in IHT_COARSE:
        c.set_param(f"Iht{j}", "dc", float(v))
        r = c.run("op")
        p.append(float(r[f"V(win1_re_{j - 1})"][0]) ** 2
                 + float(r[f"V(win1_im_{j - 1})"][0]) ** 2)
    p = np.array(p)
    return float(IHT_COARSE[int(np.argmin(p))]), p


def main() -> int:
    RESULTS.mkdir(exist_ok=True)
    op = OPERATING_POINT
    pdb = op["pdb_V"]
    c = circuit(op["p_las_mW"])
    for i in range(1, N_NEURON + 1):
        c.set_param(f"Vpdb{i}", "dc", pdb)
        c.set_param(f"Iht{i}", "dc", op["iht_A"])
        for j in range(1, 9):
            c.set_param(f"VW{i}{j}", "dc", 0.0)
    iht = np.full(N_NEURON, IHT_PARK)

    # ── design trim first: put each resonance on channel while HOT ──────────
    from build_mpc_deck import LAMBDAS_NM, N_G_MOD, n_eff_for
    print(f"design trim ({op['p_las_mW']} mW/ch — the rings are hot, so the cold "
          f"trim is wrong by construction)")
    n_eff = {}
    for j in range(1, N_NEURON + 1):
        cold = n_eff_for(LAMBDAS_NM[j - 1], N_G_MOD)
        n_eff[j] = park_n_eff(c, j, cold)
        print(f"  ring {j}: n_eff {cold:.6f} -> {n_eff[j]:.6f} "
              f"({(n_eff[j] - cold) / N_G_MOD * LAMBDAS_NM[j - 1] * 1e3:+.0f} pm of "
              f"pre-compensation)")
    for j in range(1, N_NEURON + 1):
        c.set_param(f"Iht{j}", "dc", IHT_PARK)

    print(f"\nsource gain vs Iht, {op['p_las_mW']} mW/ch, PDB {pdb} V")
    curves, grids, notches, spectra = {}, {}, {}, {}
    for j in range(1, N_NEURON + 1):
        i0, spec = notch_current(c, j)
        notches[j], spectra[j] = i0, spec
        grid = np.clip(np.linspace(i0 - FINE_SPAN, i0 + FINE_SPAN, N_FINE),
                       0.0, 5.0e-3)
        row = []
        for v in grid:
            c.set_param(f"Iht{j}", "dc", float(v))
            row.append(source_gain(c, j, pdb))
        c.set_param(f"Iht{j}", "dc", op["iht_A"])
        curves[j], grids[j] = np.array(row), grid
        print(f"  ring {j}: notch at {i0 * 1e3:.2f} mA, G over +-{FINE_SPAN * 1e3:.1f} mA: "
              + " ".join(f"{g:6.2f}" for g in row))

    # Work in signed gain: a negative column is fine (the weights are bipolar)
    # but the six have to AGREE, or the trim is chasing two different targets.
    reach = np.array([np.nanmax(curves[j]) for j in curves])
    target = float(np.nanmin(reach))
    print(f"\nhighest gain every ring can reach with the same sign: {target:.3f} "
          f"(ring {int(np.nanargmin(reach)) + 1} is the limit)")
    print(f"MPC needs {NEED}: {'OK' if target >= NEED else 'SHORT'}, "
          f"margin {target / NEED:.2f}x")

    for p in range(PASSES):
        moved = 0.0
        for j in range(1, N_NEURON + 1):
            errs = []
            for v in grids[j]:
                c.set_param(f"Iht{j}", "dc", float(v))
                errs.append(abs(source_gain(c, j, pdb) - target))
            k = int(np.nanargmin(errs))
            moved = max(moved, abs(grids[j][k] - iht[j - 1]))
            iht[j - 1] = float(grids[j][k])
            c.set_param(f"Iht{j}", "dc", iht[j - 1])
        gs = np.array([source_gain(c, j, pdb) for j in range(1, N_NEURON + 1)])
        print(f"pass {p + 1}: Iht = " + " ".join(f"{v * 1e3:4.2f}" for v in iht)
              + "  ->  G = " + " ".join(f"{g:6.2f}" for g in gs)
              + f"   spread {np.nanmax(gs) - np.nanmin(gs):.3f}")
        if moved < 1e-9:
            break

    gs = np.array([source_gain(c, j, pdb) for j in range(1, N_NEURON + 1)])
    ok = bool(np.nanmin(gs) >= NEED)
    out = {"converged": ok, "p_las_mW": op["p_las_mW"], "pdb_V": pdb,
           "n_eff": {str(j): float(v) for j, v in n_eff.items()},
           "iht_A": [float(v) for v in iht], "G": [float(g) for g in gs],
           "G_target": target, "need": NEED,
           "notch_A": {str(j): float(v) for j, v in notches.items()},
           "iht_grid_A": {str(j): [float(v) for v in grids[j]] for j in grids},
           "curves": {str(j): [float(v) for v in curves[j]] for j in curves}}
    (RESULTS / "operating_point.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"\nfinal G {np.nanmin(gs):.3f} to {np.nanmax(gs):.3f}  "
          f"{'-- usable' if ok else '-- STILL SHORT, do not use downstream'}")
    print(f"wrote {RESULTS / 'operating_point.json'}")

    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.4))
    cmap = plt.get_cmap("tab10")
    for j in curves:
        ax[0].plot(grids[j] * 1e3, curves[j], "-", color=cmap(j - 1), label=f"ring {j}")
        ax[0].plot(iht[j - 1] * 1e3, gs[j - 1], "*", ms=15, color=cmap(j - 1))
    ax[0].axhline(target, color="k", ls=":", lw=1)
    ax[0].axhline(NEED, color="k", ls="--", lw=1.2)
    ax[0].annotate("MPC needs 1.18", (0.1, NEED), fontsize=8,
                   textcoords="offset points", xytext=(2, 4))
    ax[0].set_xlabel("ring heater (mA)"); ax[0].set_ylabel("source gain G")
    ax[0].set_title("each ring's drive strength vs its own trim\n(star = chosen)")
    ax[0].legend(fontsize=7, ncol=2); ax[0].grid(alpha=0.25)
    before = [-0.007, 0.102, -0.002, 8.268, -0.002, -3.884]
    x = np.arange(1, N_NEURON + 1)
    ax[1].bar(x - 0.2, before, 0.4, color="#8c8c8c", label="one shared trim")
    ax[1].bar(x + 0.2, gs, 0.4, color="#4c72b0", label="per-ring trim")
    ax[1].axhline(NEED, color="k", ls="--", lw=1.2)
    ax[1].axhline(0, color="k", lw=0.8)
    ax[1].set_xlabel("ring"); ax[1].set_ylabel("source gain G")
    ax[1].set_title("what the trim bought")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.25, axis="y")
    fig.suptitle("giona MPC — per-ring trim, on a gain chain validated against "
                 "arithmetic", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(RESULTS / "trim.png", dpi=130)
    print(f"wrote {RESULTS / 'trim.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
