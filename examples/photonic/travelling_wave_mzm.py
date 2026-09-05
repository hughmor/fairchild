#!/usr/bin/env python3
"""A travelling-wave MZM, and the velocity mismatch that sets its bandwidth.

`fc_tw_ps` is one arm. An MZM is two of them in an interferometer, driven
push-pull, each electrode fed and terminated by the deck. Nothing in this
example is new machinery — it is the phase-2 element composed the way a real
modulator is built, which is the claim #116 phase 3 makes.

## What is being measured

The electro-optic response, end to end and optically: laser in, detector out,
`.ac` on the RF drive. What limits it here is **walk-off** — the RF and the
light travel at different speeds, so a photon entering at `t` sees a drive that
has moved on by the time it leaves. The closed form is

    H(f) = sinc(π · f · L · (n_m − n_g) / c)

and it is the dashed line on the plot. Nothing in the deck computes it: the
ladder's slices accumulate RF and optical delay at different rates and the sinc
falls out.

## What comes free, and what does not

  * **Segmented electrodes** — cascade several `fc_tw_ps`, each with its own
    drive. Same ladder, more of them.
  * **Distributed drivers** — the electrode is a real transmission line with
    real ends, so a driver sees its impedance and its reflections rather than a
    lumped capacitance.
  * **Travelling-wave detectors do not.** A distributed detector needs the
    photocurrent generated per slice, and `fc_photodetector` is not a segment
    device — it has no length to cut up. That is a new device, not a
    composition.

    MPLBACKEND=Agg python3 examples/photonic/travelling_wave_mzm.py
    python3 examples/photonic/travelling_wave_mzm.py --selftest
"""

import argparse
import math
import os
import sys

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))
import fairchild  # noqa: E402

C0 = 299_792_458.0

L_M = 3e-3          # 3 mm arms
N_G = 4.2           # optical group index
V_PI_L = 0.012      # V·m ⇒ V_pi = 4 V over 3 mm
Z0 = 35.0           # electrode impedance, loaded
F_MAX = 80e9        # sets the slice count
P_LASER_MW = 2.0
RESPONSIVITY = 0.9
R_LOAD = 1.0e3
KL_3DB = math.pi / 4

# Microwave indices to sweep: matched, then progressively faster electrodes.
N_M_SWEEP = [4.2, 3.8, 3.2, 2.5, 1.8]

FREQS = np.concatenate([np.linspace(1e9, 20e9, 20), np.linspace(21e9, 120e9, 34)])

# Sequential ramp: n_m is a magnitude, so one hue light→dark rather than
# categorical hues. Steps taken from the dataviz reference blue ramp.
RAMP = ["#9fc6f0", "#6ea8e6", "#3987e5", "#2a78d6", "#17539c"]
C_FORM = "#52514e"   # the closed form, deliberately recessive ink
INK = "#0b0b0b"
INK_2 = "#52514e"
SURFACE = "#fcfcfb"
GRID = "#e3e2dd"


def deck(n_m: float, v_bias: float = 4.0) -> str:
    """MZI with two travelling-wave arms, driven push-pull.

    The bias reaches the slices through the electrode itself: a DC source
    behind `Z0` into a line terminated in `Z0` is a resistive divider, so each
    slice sits at half the source. That the bias arrives at all is the DC
    through-connection of a lossless line (#111).
    """
    arm = (
        f"fc_tw_ps l_um={L_M * 1e6} v_pi_l={V_PI_L} n_g={N_G} n_m={n_m} "
        f"z0={Z0} f_max={F_MAX:e} alpha_dB_cm=0 pin_at_ref=1"
    )
    return f"""* travelling-wave MZM
.options waveguide_delay=1
.optical_port lin
.optical_port ldark
.optical_port a1
.optical_port a2
.optical_port b1
.optical_port b2
.optical_port obar
.optical_port ocross
Xlas lin fc_cw_laser power_mW={P_LASER_MW} wavelength_nm=1550
Vp rfp 0 DC {v_bias} AC 0.5
Vn rfn 0 DC {-v_bias} AC 0.5 180
Rsp rfp p0 {Z0}
Rsn rfn n0 {Z0}
Rtp pN 0 {Z0}
Rtn nN 0 {Z0}
Xc1 lin ldark a1 a2 fc_dcoupler kappa_L={KL_3DB}
Xarm1 a1 b1 p0 pN {arm}
Xarm2 a2 b2 n0 nN {arm}
Xc2 b1 b2 obar ocross fc_dcoupler kappa_L={KL_3DB}
Xpd obar det 0 fc_photodetector responsivity={RESPONSIVITY} r_shunt=1Meg i_dark_a=0
Rl det 0 {R_LOAD}
"""


def eo_response(n_m: float) -> np.ndarray:
    """|H(f)| at the detector, normalised to its lowest-frequency value."""
    c = fairchild.Circuit()
    c.load_str(deck(n_m))
    r = c.run("ac", fstart=FREQS[0], fstop=FREQS[-1], points=len(FREQS), variation="lin")
    mag = np.abs(np.asarray(r["V(det)"]))
    return np.asarray(r.freq()), mag / mag[0]


def walkoff(f, n_m: float) -> np.ndarray:
    x = np.pi * np.asarray(f) * L_M * (n_m - N_G) / C0
    return np.where(np.abs(x) < 1e-12, 1.0, np.abs(np.sin(x) / np.where(x == 0, 1, x)))


def f_3db(f, h) -> float:
    """First crossing of 1/√2, linearly interpolated."""
    half = 1.0 / math.sqrt(2.0)
    below = np.nonzero(h < half)[0]
    if len(below) == 0:
        return float("nan")
    k = below[0]
    if k == 0:
        return float(f[0])
    return float(np.interp(-half, [-h[k - 1], -h[k]], [f[k - 1], f[k]]))


def figure(results, path):
    import matplotlib.pyplot as plt

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11.2, 4.6))
    fig.patch.set_facecolor(SURFACE)

    for a in (ax, bx):
        a.set_facecolor(SURFACE)
        a.tick_params(colors=INK_2, labelsize=8.5, length=3, width=0.8)
        for side in ("top", "right"):
            a.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            a.spines[side].set_color(GRID)
            a.spines[side].set_linewidth(0.8)
        a.grid(True, color=GRID, linewidth=0.7, alpha=0.9)
        a.set_axisbelow(True)

    ax.set_title("EO response, and the walk-off it comes from", color=INK,
                 fontsize=10.5, loc="left", pad=8, fontweight="semibold")
    ax.set_xlabel("frequency (GHz)", color=INK_2, fontsize=9)
    ax.set_ylabel("|H(f)|, normalised", color=INK_2, fontsize=9)
    for i, ((n_m, f, h), col) in enumerate(zip(results, RAMP)):
        fg = np.asarray(f) / 1e9
        # Closed form underneath and wider, measurement on top and thinner: the
        # point is that one hides the other, and two lines of equal weight would
        # just look like one line.
        ax.plot(fg, walkoff(f, n_m), color=C_FORM, linewidth=3.0,
                linestyle=(0, (3, 2)), alpha=0.5, zorder=2)
        ax.plot(fg, h, color=col, linewidth=1.8, zorder=3)
        # Each curve is labelled where it crosses a *different* level, so the
        # labels walk down the family instead of stacking where the steep ones
        # all pass 0.9 within a few GHz of each other.
        level = 0.95 - 0.13 * i
        below = np.nonzero(h < level)[0]
        k = below[0] if len(below) else len(h) - 1
        ax.annotate(f"n$_m$ = {n_m}", (fg[k], h[k]), textcoords="offset points",
                    xytext=(7, 5), color=col, fontsize=8.5, fontweight="semibold",
                    zorder=6)
    ax.axhline(1 / math.sqrt(2), color=INK_2, linewidth=0.9, linestyle=(0, (2, 3)))
    ax.annotate("−3 dB", (FREQS[-1] / 1e9, 1 / math.sqrt(2)), textcoords="offset points",
                xytext=(-4, 5), ha="right", color=INK_2, fontsize=8.5)
    ax.annotate("dashed: sinc(π f L Δn / c) — the curves sit on it",
                (0.98, 0.90), xycoords="axes fraction", ha="right", color=C_FORM,
                fontsize=8.5)
    ax.set_ylim(0, 1.08)

    bx.set_title("Bandwidth against velocity mismatch", color=INK, fontsize=10.5,
                 loc="left", pad=8, fontweight="semibold")
    bx.set_xlabel("|n$_m$ − n$_g$|", color=INK_2, fontsize=9)
    bx.set_ylabel("3 dB bandwidth (GHz)", color=INK_2, fontsize=9)
    dn = np.array([abs(n_m - N_G) for n_m, _, _ in results])
    meas = np.array([f_3db(f, h) for _, f, h in results]) / 1e9
    # The closed form: sinc drops to 1/√2 at x = 1.3916.
    grid = np.linspace(max(dn.min(), 0.05), dn.max() * 1.05, 100)
    bx.plot(grid, 1.3916 * C0 / (np.pi * L_M * grid) / 1e9, color=C_FORM,
            linewidth=1.4, linestyle=(0, (4, 3)), label="sinc closed form")
    for x, y, col in zip(dn, meas, RAMP):
        bx.plot([x], [y], "o", color=col, markersize=7, markeredgecolor=SURFACE,
                markeredgewidth=1.4, zorder=5)
    bx.set_yscale("log")
    bx.legend(frameon=False, fontsize=8.5, labelcolor=INK_2, loc="upper right")
    # The matched arm has no 3 dB point to plot, which is the whole point of
    # matching, so it is said rather than drawn.
    matched = [n_m for n_m, _, _ in results if abs(n_m - N_G) < 1e-9]
    if matched:
        bx.annotate(
            f"n$_m$ = n$_g$ = {N_G}: flat to {FREQS[-1] / 1e9:.0f} GHz,\nno 3 dB point to plot",
            (0.04, 0.06), xycoords="axes fraction", color=INK, fontsize=8.5,
            fontweight="semibold",
        )
    fig.suptitle("A travelling-wave MZM, built from two fc_tw_ps",
                 color=INK, fontsize=12.5, x=0.008, ha="left",
                 fontweight="semibold", y=0.975)
    fig.text(0.008, 0.915,
             f"{L_M * 1e3:.0f} mm arms, {Z0:.0f} Ω electrodes, push-pull, both ends "
             "terminated. The walk-off is not modelled anywhere — it emerges from "
             "the ladder.",
             color=INK_2, fontsize=9, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(path, dpi=170, facecolor=SURFACE)
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    results = []
    print(f"{'n_m':>6} {'|n_m − n_g|':>12} {'f_3dB (GHz)':>12} {'sinc says':>11}")
    for n_m in N_M_SWEEP:
        f, h = eo_response(n_m)
        results.append((n_m, f, h))
        dn = abs(n_m - N_G)
        got = f_3db(f, h) / 1e9
        want = 1.3916 * C0 / (math.pi * L_M * dn) / 1e9 if dn > 1e-9 else float("inf")
        print(f"{n_m:>6.1f} {dn:>12.2f} {got:>12.1f} {want:>11.1f}")

    if args.selftest:
        for n_m, f, h in results:
            dn = abs(n_m - N_G)
            if dn < 1e-9:
                assert h.min() > 0.99, (
                    f"velocity matched must be flat, dipped to {h.min():.3f}"
                )
                continue
            got = f_3db(f, h)
            want = 1.3916 * C0 / (math.pi * L_M * dn)
            assert abs(got / want - 1.0) < 0.12, (
                f"n_m={n_m}: 3 dB at {got / 1e9:.1f} GHz, sinc says {want / 1e9:.1f} GHz"
            )
        print("\nselftest ok: every arm's bandwidth matches the walk-off sinc")
        return

    out = args.out or os.path.join(os.path.dirname(__file__), "travelling_wave_mzm.png")
    figure(results, out)


if __name__ == "__main__":
    main()
