#!/usr/bin/env python3
"""rnn_plots.py — figures for the giona RNN calibration and programming.

Produces, into `results/`:

  rnn_activation.png     the two candidate activations vs junction current, the
                         diode clamp, the heater knob, the calibrated comb
  rnn_rest_point.png     how to choose the rest current: dT/dI, the AC coupling
                         efficiency, their product (the loop gain), and the
                         measured gain map over (operating point, laser power)
  rnn_dynamics.png       transients across the oscillation threshold

Reads the transient traces that `rnn_hunt` saved as npz; recomputes the
single-ring curves directly (cheap).

Run:  .venv/bin/python experiments/giona/rnn_plots.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rnn_explore as ex  # noqa: E402

RESULTS = HERE / "results"
N_DIODE, V_T = 5.0, 0.025852          # from the fitted card
R_SHUNT_PAR = 1.0 / (1 / 2e3 + 1 / 10e3)   # R1 (on chip) || Rb3 (PCB)

# Loop gain measured on the full deck (rnn_explore --gain), rows PD_B, cols mW/ch.
GAIN_PDB = [-5, -6, -7.2, -8, -9, -10, -12, -15]
GAIN_POW = [12, 30, 60, 120]
GAIN = np.array([
    [-0.0448, -0.1777, -0.4836, -1.0987],
    [-0.0919, -0.4237, -1.2184, -2.8542],
    [+0.0245, -0.0126, -0.4788, -2.3820],
    [+0.1028, +0.4065, +0.7261, -0.2624],
    [+0.1069, +0.5380, +1.5793, +2.4466],
    [+0.0755, +0.4125, +1.5290, +3.8925],
    [+0.0292, +0.1528, +0.6895, +3.2898],
    [-0.0103, -0.0829, -0.2488, -0.0252],
])


def fig_rest_point(out):
    """How to choose the rest current — the two competing effects and the map."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    radii = ex.stage_radii(verbose=False)
    r1, lam = radii[0], ex.LAMBDAS_NM[0]

    I = np.linspace(5e-6, 700e-6, 60)
    T = np.array([ex.probe_i(r1, lam, float(i))[0] for i in I])
    dT_dI = np.gradient(T, I)                      # mW per A

    # AC coupling efficiency: the junction's small-signal conductance against
    # the resistive load it competes with. Only this fraction of the AC
    # photocurrent drives the modulator; the rest leaks off chip.
    g_d = I / (N_DIODE * V_T)
    eta = g_d / (g_d + 1.0 / R_SHUNT_PAR)

    fig, ax = plt.subplots(2, 2, figsize=(12, 8))

    ax[0, 0].plot(I * 1e6, T, lw=1.6)
    ax[0, 0].set_xlabel("rest junction current (µA)")
    ax[0, 0].set_ylabel("transmission (mW)")
    ax[0, 0].set_title("Activation T(I) — notch port, on resonance at rest")
    ax[0, 0].grid(alpha=0.3)

    a = ax[0, 1]
    a.plot(I * 1e6, dT_dI / dT_dI.max(), label="dT/dI  (falls)", lw=1.6)
    a.plot(I * 1e6, eta, label="AC coupling η  (rises)", lw=1.6)
    prod = (dT_dI / dT_dI.max()) * eta
    # Single-ring estimate only: it ignores the 8-ring cascade, where each
    # ring's transmission also changes what every downstream ring sees. The
    # measured map (bottom right) is the authority — and it disagrees, peaking
    # at a much higher rest current.
    a.plot(I * 1e6, prod / prod.max(), "k--", lw=1.8,
           label="product (single-ring estimate)")
    i_best = I[int(np.argmax(prod))]
    a.axvline(i_best * 1e6, color="r", ls=":", lw=1.0,
              label=f"optimum ≈ {i_best * 1e6:.0f} µA")
    a.set_xlabel("rest junction current (µA)")
    a.set_ylabel("normalised")
    a.set_title("Choosing the rest current: two competing effects")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    a = ax[1, 0]
    a.plot(I * 1e6, eta * 100, lw=1.6)
    a.set_xlabel("rest junction current (µA)")
    a.set_ylabel("AC photocurrent reaching the junction (%)")
    a.set_title("Raising the rest point recovers AC gain\n"
                "(g_d rises against the fixed 2k‖10k load)", fontsize=10)
    a.grid(alpha=0.3)

    a = ax[1, 1]
    im = a.imshow(GAIN, aspect="auto", cmap="RdBu_r", origin="upper",
                  vmin=-3, vmax=3,
                  extent=[-0.5, len(GAIN_POW) - 0.5, len(GAIN_PDB) - 0.5, -0.5])
    a.set_xticks(range(len(GAIN_POW)))
    a.set_xticklabels(GAIN_POW)
    a.set_yticks(range(len(GAIN_PDB)))
    a.set_yticklabels(GAIN_PDB)
    a.set_xlabel("laser power (mW/channel)")
    a.set_ylabel("PD_B (V)")
    a.set_title("Measured loop gain G per unit weight\n"
                "|G|·√2 > 1 needed for W=[[1,-1],[1,1]]", fontsize=10)
    for i in range(len(GAIN_PDB)):
        for j in range(len(GAIN_POW)):
            g = GAIN[i, j]
            a.text(j, i, f"{g:+.2f}", ha="center", va="center", fontsize=7,
                   color="white" if abs(g) > 1.6 else "black")
    fig.colorbar(im, ax=a)

    fig.suptitle("giona RNN — picking the rest current and the loop gain it buys")
    fig.tight_layout()
    fig.savefig(out, dpi=115)
    print(f"wrote {out}")


def fig_dynamics(out, traces):
    """Transients across the threshold, from saved npz files."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(traces)
    fig, ax = plt.subplots(2, n, figsize=(4.6 * n, 7), squeeze=False)
    for c, (label, f) in enumerate(traces):
        d = np.load(f)
        t, v, bus = d["t"] * 1e9, d["v"], d["bus"]
        ax[0, c].plot(t, v[0] * 1e3, lw=1.2, label="neuron 1")
        ax[0, c].plot(t, v[1] * 1e3, lw=1.2, label="neuron 2")
        ax[0, c].set_ylabel("V(mod_cathode) (mV)")
        ax[0, c].set_title(label, fontsize=10)
        ax[0, c].legend(fontsize=8)
        ax[0, c].grid(alpha=0.3)
        ax[1, c].plot(t, bus[0], lw=1.2, label="λ1 on bus")
        ax[1, c].plot(t, bus[1], lw=1.2, label="λ2 on bus")
        ax[1, c].set_xlabel("time (ns)")
        ax[1, c].set_ylabel("bus power (mW)")
        ax[1, c].legend(fontsize=8)
        ax[1, c].grid(alpha=0.3)
    fig.suptitle("giona RNN — two-neuron dynamics, W = [[1,−1],[1,1]]")
    fig.tight_layout()
    fig.savefig(out, dpi=115)
    print(f"wrote {out}")


def fig_wta(out, files):
    """Winner-take-all: the latch, the optical outcome, and the decision curve."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    runs = []
    for f in files:
        d = np.load(f)
        runs.append(dict(t=d["t"] * 1e9, v=d["v"], bus=d["bus"],
                         in_w=d["in_w"], p=float(d["p"])))
    runs.sort(key=lambda r: r["in_w"][0] - r["in_w"][1])
    seen, uniq = set(), []
    for r in runs:
        key = round(float(r["in_w"][0] - r["in_w"][1]), 3)
        if key not in seen:
            seen.add(key)
            uniq.append(r)
    runs = uniq

    # The two most decisive runs, for the time-domain panels.
    lo, hi = runs[0], runs[-1]
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))

    for col, r in ((0, lo), (1, hi)):
        w = r["in_w"]
        ax[0, col].plot(r["t"], r["v"][0] * 1e3, lw=1.4, label="neuron 1")
        ax[0, col].plot(r["t"], r["v"][1] * 1e3, lw=1.4, label="neuron 2")
        ax[0, col].axvline(5, color="k", ls=":", lw=0.8)
        ax[0, col].set_title(f"input (w13, w23) = ({w[0]:.2f}, {w[1]:.2f})",
                             fontsize=10)
        ax[0, col].set_ylabel("V(mod_cathode) (mV)")
        ax[0, col].legend(fontsize=8)
        ax[0, col].grid(alpha=0.3)
        ax[1, col].plot(r["t"], r["bus"][0], lw=1.4, label="λ1 on bus")
        ax[1, col].plot(r["t"], r["bus"][1], lw=1.4, label="λ2 on bus")
        ax[1, col].axvline(5, color="k", ls=":", lw=0.8)
        ax[1, col].set_xlabel("time (ns)")
        ax[1, col].set_ylabel("bus power (mW)")
        ax[1, col].legend(fontsize=8)
        ax[1, col].grid(alpha=0.3)

    # Decision curve: final differential vs input asymmetry.
    asym = np.array([r["in_w"][0] - r["in_w"][1] for r in runs])
    split = np.array([(r["v"][0][-1] - r["v"][1][-1]) * 1e3 for r in runs])
    a = ax[0, 2]
    a.plot(asym, split, "o-", lw=1.5)
    a.axhline(0, color="k", lw=0.6)
    a.axvline(0, color="k", lw=0.6)
    a.set_xlabel("input asymmetry  w13 − w23")
    a.set_ylabel("final V1 − V2 (mV)")
    a.set_title("Decision curve: the input picks the basin", fontsize=10)
    a.grid(alpha=0.3)

    a = ax[1, 2]
    contrast = np.array([r["bus"][0][-1] - r["bus"][1][-1] for r in runs])
    a.plot(asym, contrast, "s-", lw=1.5, color="tab:green")
    a.axhline(0, color="k", lw=0.6)
    a.axvline(0, color="k", lw=0.6)
    a.set_xlabel("input asymmetry  w13 − w23")
    a.set_ylabel("final bus λ1 − λ2 (mW)")
    a.set_title("Optical contrast between the two channels", fontsize=10)
    a.grid(alpha=0.3)

    fig.suptitle(f"giona RNN — winner-take-all, W_rec = [[1,−1],[−1,1]], "
                 f"{runs[0]['p']:.0f} mW/channel")
    fig.tight_layout()
    fig.savefig(out, dpi=115)
    print(f"wrote {out}")


def fig_ac_weights(out, files):
    """AC-driven weights: waveforms, the decision variable, the phase portrait.

    Note the geometry. Independent square and triangle weights have a varying
    SUM as well as a varying difference, and the sum is common-mode drive that
    both neurons follow together. So the trajectory travels along the V1 = V2
    diagonal while the decision lives in the (smaller) offset from it — which is
    why the middle row plots the differential explicitly.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    runs = [np.load(f) for f in files]
    n = len(runs)
    fig, ax = plt.subplots(3, n, figsize=(8.4 * n, 12.4), squeeze=False)

    for c, d in enumerate(runs):
        t, v, w, p = d["t"] * 1e9, d["v"], d["w"], float(d["p"])
        asym = w[0] - w[1]
        diff = (v[0] - v[1]) * 1e3
        agree = np.mean(np.sign(diff[asym != 0]) == np.sign(asym[asym != 0])) * 100

        # ── row 0: weight waveforms and both neuron outputs ─────────────────
        a = ax[0, c]
        a.plot(t, w[0], lw=1.0, color="tab:blue", alpha=0.8, label="w13 (square)")
        a.plot(t, w[1], lw=1.0, color="tab:orange", alpha=0.8, label="w23 (triangle, 2×f)")
        a.set_ylabel("weight (V)")
        a2 = a.twinx()
        a2.plot(t, v[0] * 1e3, lw=1.4, color="tab:green", label="neuron 1")
        a2.plot(t, v[1] * 1e3, lw=1.4, color="tab:red", label="neuron 2")
        a2.set_ylabel("V(mod_cathode) (mV)")
        h1, l1 = a.get_legend_handles_labels()
        h2, l2 = a2.get_legend_handles_labels()
        a.legend(h1 + h2, l1 + l2, fontsize=8, loc="lower center", ncol=2)
        a.set_title(f"{p:.0f} mW/channel — weights and neuron outputs", fontsize=11)
        a.set_xlabel("time (ns)")
        a.grid(alpha=0.25)

        # ── row 1: the decision variable against the input asymmetry ────────
        a = ax[1, c]
        a.plot(t, asym, lw=1.1, color="tab:purple", label="w13 − w23 (input)")
        a.axhline(0, color="k", lw=0.6)
        a.set_ylabel("input asymmetry (V)", color="tab:purple")
        a.set_xlabel("time (ns)")
        a2 = a.twinx()
        a2.plot(t, diff, lw=1.4, color="k", label="V1 − V2 (decision)")
        a2.axhline(0, color="grey", lw=0.5, ls=":")
        a2.set_ylabel("V1 − V2 (mV)")
        h1, l1 = a.get_legend_handles_labels()
        h2, l2 = a2.get_legend_handles_labels()
        a.legend(h1 + h2, l1 + l2, fontsize=8, loc="lower center")
        # Where the decision actually errs: a static input-referred offset, not
        # dynamics. Fit diff = k*(asym - offset) and report both.
        k, b = np.polyfit(asym, diff, 1)
        off = -b / k
        a2.axhline(0, color="grey", lw=0.5, ls=":")
        a.axvline(off, color="tab:brown", lw=0.9, ls="--",
                  label=f"offset {off:+.3f}")
        h1, l1 = a.get_legend_handles_labels()
        h2, l2 = a2.get_legend_handles_labels()
        a.legend(h1 + h2, l1 + l2, fontsize=8, loc="lower center")
        a.set_title(f"decision tracks the asymmetry — agreement {agree:.1f}%, "
                    f"slope {k:.0f} mV/unit, offset {off:+.3f}", fontsize=10)
        a.grid(alpha=0.25)

        # ── row 2: phase portrait, coloured by who SHOULD win ───────────────
        # Positive asymmetry suppresses neuron 1 (the neuron path inverts), so
        # w13 > w23 means neuron 2 ought to be the winner.
        a = ax[2, c]
        n2 = asym > 0
        a.scatter(v[0][n2] * 1e3, v[1][n2] * 1e3, s=8, alpha=0.22,
                  color="tab:red", label="w13 > w23  (N2 should win)")
        a.scatter(v[0][~n2] * 1e3, v[1][~n2] * 1e3, s=8, alpha=0.22,
                  color="tab:green", label="w13 < w23  (N1 should win)")
        lo = min(v[0].min(), v[1].min()) * 1e3
        hi = max(v[0].max(), v[1].max()) * 1e3
        pad = 0.04 * (hi - lo)
        a.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k--", lw=0.8,
               alpha=0.7, label="V1 = V2 (no winner)")
        a.set_xlabel("neuron 1  V(mod_cathode) (mV)")
        a.set_ylabel("neuron 2  V(mod_cathode) (mV)")
        a.set_title(f"phase portrait — basins split across the diagonal, "
                    f"max |V1−V2| = {np.abs(diff).max():.1f} mV", fontsize=10)
        a.legend(fontsize=8, loc="upper left")
        a.grid(alpha=0.25)
        a.set_aspect("equal", adjustable="datalim")

    fig.suptitle("giona RNN — winner-take-all under AC weights "
                 "(square w13, triangle w23 at 2× frequency)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, dpi=115)
    print(f"wrote {out}")


def main() -> int:
    RESULTS.mkdir(exist_ok=True)
    fig_rest_point(RESULTS / "rnn_rest_point.png")
    scratch = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    if scratch and scratch.is_dir():
        wta = sorted(scratch.glob("w3_30_*.npz")) + sorted(scratch.glob("c3_*.npz"))
        if len(wta) >= 3:
            fig_wta(RESULTS / "rnn_wta.png", wta)
        ac = sorted(scratch.glob("ac_*.npz"), key=lambda f: float(np.load(f)["p"]))
        if ac:
            fig_ac_weights(RESULTS / "rnn_wta_ac.png", ac)
    return 0


if __name__ == "__main__":
    sys.exit(main())
