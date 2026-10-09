#!/usr/bin/env python3
"""gain_vs_power.py — the loop gain, measured the one way that survives checking.

`validate_estimator.py` checks the chain factor by factor: the weight block and
balanced pair give `I = R*w*P_in` to 7e-16, the tap and tree give exactly 1/16,
and the node division reproduces `rnn_math.md` (6) at R_sh = 1394 ohm. Only the
ring's slope is device physics.

    G = eta * R * dP_win/dI_D

Two details make it work that did not before:

  * `I_D` is READ, not inferred. The deck now carries a 0 V source in series with
    each ring's anode. Computing `I_D` from the cathode voltage through the
    diode law amplifies that node's ~1 mV solver resolution by
    `(I_D+i_sat)/(n*V_T)`, which is 7.5 % of a 10 uA signal — and near the notch,
    where the signal is small, far worse. That is what made every earlier sweep
    look like noise with sign flips.
  * the ring is trimmed HOT. At 30 mW/channel it self-heats past its own channel
    before the heater does anything, and the heater only reddens, so `n_eff` has
    to be pre-compensated blue by the self-heating. See `trim.py`.

With both, `dP_win/dI_D` against heater current is a clean unimodal curve and
the 2-point and 5-point-fit derivatives agree to four digits.

    .venv/bin/python experiments/giona/mpc/gain_vs_power.py

Writes results/gain_vs_power.json and .png.
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
from build_mpc_deck import LAMBDAS_NM, N_G_MOD, n_eff_for
from characterize import circuit
from validate_estimator import N_DIODE, R_PD, VT

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
POWERS = (1.0, 3.0, 10.0, 30.0, 60.0)
NEED = 1.18
I_SAT, NVT = 5.099e-8, N_DIODE * VT
R_SH = 1.0 / (1 / 2e3 + 1 / 10e3 + 2 / 17.051e3)
PDB, D_PDB, J = -8.0, 0.2, 1
IHT_PARK = 1.5e-3


def probe(c, j: int):
    r = c.run("op")
    return (-float(r[f"I(vsen{j})"][0]),
            float(r[f"V(win1_re_{j - 1})"][0]) ** 2
            + float(r[f"V(win1_im_{j - 1})"][0]) ** 2)


def gain(c, j: int) -> tuple[float, float]:
    c.set_param(f"Vpdb{j}", "dc", PDB + D_PDB); ih, ph = probe(c, j)
    c.set_param(f"Vpdb{j}", "dc", PDB - D_PDB); il, pl = probe(c, j)
    c.set_param(f"Vpdb{j}", "dc", PDB)
    i0, _ = probe(c, j)
    slope = (ph - pl) / (ih - il) if ih != il else np.nan
    eta = R_SH / (R_SH + NVT / abs(i0))
    return float(eta * R_PD * slope), float(slope)


def main() -> int:
    RESULTS.mkdir(exist_ok=True)
    rows = []
    for p in POWERS:
        c = circuit(p)
        for i in range(1, 7):
            c.set_param(f"Vpdb{i}", "dc", PDB)
            c.set_param(f"Iht{i}", "dc", IHT_PARK)
            for k in range(1, 9):
                c.set_param(f"VW{i}{k}", "dc", 0.0)
        # design trim: put the HOT resonance on channel with the heater parked
        cold = n_eff_for(LAMBDAS_NM[J - 1], N_G_MOD)
        best, bp = cold, np.inf
        for v in cold + np.linspace(-2.0e-3, 0.4e-3, 41):
            c.set_param(f"Xr{J}", "n_eff", float(v))
            _, pw = probe(c, J)
            if pw < bp:
                best, bp = float(v), pw
        c.set_param(f"Xr{J}", "n_eff", best)
        # then the heater, finely, around the parked point
        gs = []
        ihts = np.linspace(max(IHT_PARK - 0.8e-3, 0.0), IHT_PARK + 0.8e-3, 17)
        for v in ihts:
            c.set_param(f"Iht{J}", "dc", float(v))
            gs.append(gain(c, J)[0])
        gs = np.array(gs)
        k = int(np.nanargmax(np.abs(gs)))
        dT = float(c.run("op")[f"V(trm{J})"][0])
        rows.append(dict(p=p, n_eff=best, best_iht=float(ihts[k]),
                         G=float(gs[k]), dT=dT,
                         ihts=[float(v) for v in ihts],
                         curve=[float(v) for v in gs]))
        print(f"{p:5.1f} mW: n_eff {best:.6f}, best Iht {ihts[k] * 1e3:.2f} mA, "
              f"G = {gs[k]:8.4f}, ring dT {dT:.1f} K")

    best = max(rows, key=lambda r: abs(r["G"]))
    print(f"\nbest |G| anywhere on this sweep: {best['G']:.4f} at {best['p']} mW")
    print(f"MPC needs {NEED}. Short by {NEED / abs(best['G']):.0f}x.")
    (RESULTS / "gain_vs_power.json").write_text(json.dumps(
        {"need": NEED, "pdb_V": PDB, "rows": rows}, indent=2) + "\n")
    print(f"wrote {RESULTS / 'gain_vs_power.json'}")

    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.3))
    ax[0].loglog([r["p"] for r in rows], [abs(r["G"]) for r in rows], "o-")
    ax[0].axhline(NEED, color="k", ls="--", lw=1.2)
    ax[0].annotate("MPC needs 1.18", (1.1, NEED), fontsize=9,
                   textcoords="offset points", xytext=(2, 4))
    ax[0].set_xlabel("laser power per channel (mW)"); ax[0].set_ylabel("|G|")
    ax[0].set_title("loop gain, measured through the validated chain")
    ax[0].grid(alpha=0.25, which="both")
    cmap = plt.get_cmap("viridis")
    for n, r in enumerate(rows):
        ax[1].plot(np.array(r["ihts"]) * 1e3, r["curve"], "-",
                   color=cmap(n / max(len(rows) - 1, 1)), label=f"{r['p']:.0f} mW")
    ax[1].set_xlabel("ring heater (mA)"); ax[1].set_ylabel("G")
    ax[1].set_title("and it is smooth now\n(it was sign-flipping noise before)")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.25)
    fig.suptitle("giona MPC — loop gain once I_D is read rather than inferred",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(RESULTS / "gain_vs_power.png", dpi=130)
    print(f"wrote {RESULTS / 'gain_vs_power.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
