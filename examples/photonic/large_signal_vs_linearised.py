#!/usr/bin/env python3
"""Where a linearised modulator model stops being trustworthy, and by how much.

A frequency-domain tool answers a modulator by measuring one small-signal
response `H(f)` at the bias point and filtering the drive with it. That is
orders of magnitude cheaper than a transient co-solve and it is the right tool
for a long PRBS. This example is about the case where it is not: a
reverse-biased PN phase shifter whose junction capacitance moves with the very
voltage that drives it.

    C_j(V) = C_j0 / (1 − V/V_bi)^m_j        m_j = 0.5, V_bi = 0.917 V

Over a 4 V swing that capacitance changes by 2.3x, so the RC pole the drive sees
is not a property of the modulator — it is a property of the modulator *and the
bit being sent*. A model holding it fixed does not make a small error. It
answers a different circuit.

## The controlled experiment

Two runs of the same solver, the same optics, the same PRBS, the same driver:

  * **C_j(V)** — `m_j = 0.5`, the real junction.
  * **C_j fixed** — `m_j = 0`, `c_j0` pinned to `C_j(V_bias)`. This is exactly
    the circuit a small-signal extraction sees, so it is exactly what a
    linearised flow is answering.

Only the capacitance law differs. The MZI's optical transfer is `sin²(Δφ/2)` in
both, so the (large) optical nonlinearity is common to the two and cancels out
of the comparison. What is left is the electrical nonlinearity alone, which is
the claim.

    MPLBACKEND=Agg python3 examples/photonic/large_signal_vs_linearised.py
    python3 examples/photonic/large_signal_vs_linearised.py --selftest
"""

import argparse
import math
import os
import sys

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))
import fairchild  # noqa: E402

# ── link design ──────────────────────────────────────────────────────────────
BIT_S = 50e-12                  # 20 Gb/s
STEP_S = 0.5e-12
EDGE_FRAC = 0.2                 # driver rise/fall as a fraction of a bit

P_LASER_MW = 2.0
L_ARM_UM = 3000.0               # 3 mm arms
V_PI_L = 0.012                  # V·m  ⇒ V_pi = 4 V per arm at 3 mm
C_J0 = 750e-15                  # 250 fF/mm
V_BI, M_J = 0.917, 0.5
ALPHA_DB_CM = 2.0
R_DRV = 50.0                    # driver output impedance per arm
V_CENTRE = -3.0                 # both arms reverse-biased

RESPONSIVITY = 0.9
R_LOAD = 1.0e3
C_PD = 5e-15                    # small, so the receiver is not the bottleneck

KL_3DB = math.pi / 4            # kappa·L for a 50/50 coupler
V_PI = V_PI_L / (L_ARM_UM * 1e-6)

PRBS_ORDER = 7
N_SYM = (1 << PRBS_ORDER) - 1
PRBS_TAPS = {7: (7, 6), 9: (9, 5)}

# Swings to sweep, as a fraction of V_pi differential. 1.0 is a full-swing
# driver; the small end is where the two models must agree, and that agreement
# is the control that says the comparison is set up correctly.
SWING_FRACS = [0.1, 0.2, 0.35, 0.5, 0.7, 0.85, 1.0]

# ── palette (validated: see the dataviz reference palette) ───────────────────
C_FULL = "#2a78d6"      # slot 1 — the full co-solve
C_LIN = "#eb6834"       # slot 2 — the linearised capacitance
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_MUTED = "#8a8983"
SURFACE = "#fcfcfb"
GRID = "#e3e2dd"


def c_j(v_pn: float) -> float:
    """The junction capacitance the model actually stamps, at bias `v_pn`."""
    return C_J0 / (1.0 - v_pn / V_BI) ** M_J


def prbs(order: int, n: int) -> list:
    """Maximal-length LFSR sequence, Fibonacci form, all-ones seed."""
    taps = PRBS_TAPS[order]
    reg = (1 << order) - 1
    out = []
    for _ in range(n):
        fb = 0
        for t in taps:
            fb ^= (reg >> (t - 1)) & 1
        out.append(reg & 1)
        reg = ((reg << 1) | fb) & ((1 << order) - 1)
    return out


def arm_dc(sign: int) -> float:
    """This arm's static bias, including the quadrature offset.

    A push-pull MZI biased at `Δφ = 0` sits on a transfer *extremum*, where the
    response is second order in the drive and both rails move the same way — the
    eye never opens, whatever the model. Quadrature is `Δφ = π/2`, which is a
    differential `V_pi/2`, so each arm carries a static `±V_pi/4`.
    """
    return V_CENTRE + sign * 0.25 * V_PI


def pwl(bits: list, sign: int, swing_frac: float) -> str:
    """Push-pull drive for one arm, differential swing = `swing_frac`·V_pi."""
    half = 0.5 * EDGE_FRAC * BIT_S
    quarter = 0.25 * swing_frac * V_PI  # per arm, so the pair swings half of V_pi·frac

    def level(b):
        return arm_dc(sign) + sign * (quarter if b else -quarter)

    pts = [f"0 {level(bits[0]):.6g}"]
    prev = bits[0]
    for k in range(1, len(bits)):
        if bits[k] != prev:
            pts.append(f"{k * BIT_S - half:.6e} {level(prev):.6g}")
            pts.append(f"{k * BIT_S + half:.6e} {level(bits[k]):.6g}")
            prev = bits[k]
    pts.append(f"{len(bits) * BIT_S:.6e} {level(prev):.6g}")
    return f"PWL({' '.join(pts)})"


PORTS = "".join(
    f".optical_port {p}\n" for p in ("lin", "ldark", "a1", "a2", "b1", "b2", "obar", "ocross")
)


def deck(bits: list, swing_frac: float, linearised: bool) -> str:
    """The link, with the junction capacitance either moving or pinned.

    `linearised` swaps `C_j(V)` for the constant a small-signal extraction at
    the bias point would have measured. Nothing else changes.
    """
    def cap(sign: int) -> str:
        # The pinned value is per arm, because the two arms sit at different
        # bias points once the quadrature offset is in. That is what an
        # extraction at the operating point would have measured for each.
        if linearised:
            return f"c_j0={c_j(arm_dc(sign)):.6e} m_j=0"
        return f"c_j0={C_J0:.6e} m_j={M_J}"

    def arm(sign: int) -> str:
        return (
            f"fc_pn_ps_cap l_um={L_ARM_UM} v_pi_l={V_PI_L} {cap(sign)} "
            f"v_bi={V_BI} alpha_dB_cm={ALPHA_DB_CM} pin_at_ref=1"
        )

    return f"""* MZM link, junction capacitance under test
{PORTS}Xlas lin fc_cw_laser power_mW={P_LASER_MW}
Vp pd 0 {pwl(bits, +1, swing_frac)}
Vn nd 0 {pwl(bits, -1, swing_frac)}
Rp pd p {R_DRV}
Rn nd n {R_DRV}
Xc1 lin ldark a1 a2 fc_dcoupler kappa_L={KL_3DB}
Xps1 a1 b1 p 0 {arm(+1)}
Xps2 a2 b2 n 0 {arm(-1)}
Xc2 b1 b2 obar ocross fc_dcoupler kappa_L={KL_3DB}
Xpd obar det 0 fc_photodetector responsivity={RESPONSIVITY} r_shunt=1Meg i_dark_a=0
Rl det 0 {R_LOAD}
Cl det 0 {C_PD}
"""


def run(bits: list, swing_frac: float, linearised: bool) -> dict:
    """One transient. Returns the detector waveform and the junction drive."""
    c = fairchild.Circuit()
    c.load_str(deck(bits, swing_frac, linearised))
    r = c.run("tran", step=STEP_S, stop=len(bits) * BIT_S)
    t = np.asarray(r.time())
    return {
        "t": t,
        "det": np.asarray(r["V(det)"]),
        "vdiff": np.asarray(r["V(p)"]) - np.asarray(r["V(n)"]),
    }


def eye_metrics(t, y, bits) -> dict:
    """Worst-case eye opening and the two rail levels, sampled at bit centres.

    Skips the first eight symbols: the run starts from the operating point and
    the first few bits carry the settling, which is not eye closure.

    Polarity is read from the data rather than assumed. Which MZI output port
    is bright at the bias point depends on the arm phases, and hard-coding
    "ones are the high rail" gives a negative opening on the other port —
    a number that still plots, still trends, and means nothing.
    """
    skip = 8
    centres = [(k + 0.5) * BIT_S for k in range(skip, len(bits))]
    samples = np.interp(centres, t, y)
    ones = np.array([s for s, b in zip(samples, bits[skip:]) if b])
    zeros = np.array([s for s, b in zip(samples, bits[skip:]) if not b])
    hi, lo = (ones, zeros) if ones.mean() > zeros.mean() else (zeros, ones)
    return {
        "opening": float(hi.min() - lo.max()),
        "hi_mean": float(hi.mean()),
        "lo_mean": float(lo.mean()),
    }


def edge_times(t, v, bits) -> tuple:
    """Median 10-90 rise and fall time of the differential junction voltage.

    The qualitative half of the claim. A fixed pole has one time constant, so
    its rise and its fall are the same length by construction — it *cannot*
    produce an asymmetric edge. A junction whose capacitance depends on the
    instantaneous voltage produces one whether you model it or not.
    """
    lo, hi = np.percentile(v, 2), np.percentile(v, 98)
    a, b = lo + 0.1 * (hi - lo), lo + 0.9 * (hi - lo)
    rises, falls = [], []
    for k in range(8, len(bits) - 1):
        if bits[k] == bits[k + 1]:
            continue
        m = (t >= (k + 1) * BIT_S - 0.3 * BIT_S) & (t <= (k + 1) * BIT_S + 0.7 * BIT_S)
        tt, vv = t[m], v[m]
        if len(tt) < 5:
            continue
        if vv[-1] > vv[0]:
            rises.append(np.interp(b, vv, tt) - np.interp(a, vv, tt))
        else:
            rev_v, rev_t = vv[::-1], tt[::-1]
            falls.append(np.interp(a, rev_v, rev_t) - np.interp(b, rev_v, rev_t))
    return (float(np.median(rises)) if rises else float("nan"),
            float(np.median(falls)) if falls else float("nan"))


def sweep(bits: list) -> dict:
    """Both models at every swing, plus the waveforms at full swing."""
    out = {"swing": [], "full": [], "lin": [], "edges_full": [], "edges_lin": []}
    for frac in SWING_FRACS:
        a = run(bits, frac, linearised=False)
        b = run(bits, frac, linearised=True)
        out["swing"].append(frac)
        out["full"].append(eye_metrics(a["t"], a["det"], bits))
        out["lin"].append(eye_metrics(b["t"], b["det"], bits))
        out["edges_full"].append(edge_times(a["t"], a["vdiff"], bits))
        out["edges_lin"].append(edge_times(b["t"], b["vdiff"], bits))
        if frac == SWING_FRACS[-1]:
            out["wave_full"], out["wave_lin"] = a, b
    return out


# ── the figure ───────────────────────────────────────────────────────────────
def fold_eye(t, y, n_ui=2):
    """Fold a waveform into `n_ui` unit intervals, one segment per crossing."""
    span = n_ui * BIT_S
    segs = []
    k = 8
    while (k + n_ui) * BIT_S <= t[-1]:
        t0 = k * BIT_S
        m = (t >= t0) & (t <= t0 + span)
        segs.append(((t[m] - t0) * 1e12, y[m]))
        k += 1
    return segs


def style_axes(ax, title, xlabel, ylabel):
    ax.set_facecolor(SURFACE)
    ax.set_title(title, color=INK, fontsize=10.5, loc="left", pad=8, fontweight="semibold")
    ax.set_xlabel(xlabel, color=INK_2, fontsize=9)
    ax.set_ylabel(ylabel, color=INK_2, fontsize=9)
    ax.tick_params(colors=INK_2, labelsize=8.5, length=3, width=0.8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
        ax.spines[side].set_linewidth(0.8)
    ax.grid(True, color=GRID, linewidth=0.7, alpha=0.9)
    ax.set_axisbelow(True)


def figure(res, bits, path):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.6))
    fig.patch.set_facecolor(SURFACE)
    a, b, c, d = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    # ── A: the cause. C_j over the swing the drive actually covers.
    v = np.linspace(V_CENTRE - 0.5 * V_PI - 0.4, V_CENTRE + 0.5 * V_PI + 0.4, 400)
    style_axes(a, "The cause: the junction capacitance moves with the drive",
               "junction voltage (V)", "C_j (fF)")
    lo, hi = V_CENTRE - 0.5 * V_PI, V_CENTRE + 0.5 * V_PI
    a.axvspan(lo, hi, color=INK_MUTED, alpha=0.10, linewidth=0)
    a.plot(v, [c_j(x) * 1e15 for x in v], color=C_FULL, linewidth=2.0)
    a.axhline(c_j(V_CENTRE) * 1e15, color=C_LIN, linewidth=2.0, linestyle=(0, (5, 3)))
    a.plot([V_CENTRE], [c_j(V_CENTRE) * 1e15], "o", color=INK, markersize=5, zorder=5)
    a.annotate("bias", (V_CENTRE, c_j(V_CENTRE) * 1e15), textcoords="offset points",
               xytext=(6, 8), color=INK, fontsize=8.5)
    a.annotate("C$_j$(V) — what the junction does", (v[10], c_j(v[10]) * 1e15),
               textcoords="offset points", xytext=(8, -14), color=C_FULL, fontsize=8.5)
    a.annotate("C$_j$ fixed — what a small-signal extraction sees",
               (v[-1], c_j(V_CENTRE) * 1e15), textcoords="offset points",
               xytext=(-6, 8), ha="right", color=C_LIN, fontsize=8.5)
    ratio = c_j(hi) / c_j(lo)
    a.annotate(f"{ratio:.1f}× across one full-swing bit", (0.5 * (lo + hi), c_j(lo) * 1e15),
               textcoords="offset points", xytext=(0, -22), ha="center",
               color=INK_2, fontsize=8.5)

    # ── B: the junction voltage, where the asymmetry appears first.
    style_axes(b, "Differential junction voltage, full swing",
               "time (ps)", "V(p) − V(n)  (V)")
    wf, wl = res["wave_full"], res["wave_lin"]
    t0, t1 = 10 * BIT_S, 18 * BIT_S
    for w, col, lab in ((wl, C_LIN, "C$_j$ fixed"), (wf, C_FULL, "C$_j$(V)")):
        m = (w["t"] >= t0) & (w["t"] <= t1)
        b.plot((w["t"][m] - t0) * 1e12, w["vdiff"][m], color=col, linewidth=2.0, label=lab)
    b.legend(frameon=False, fontsize=8.5, labelcolor=INK_2, loc="upper right")

    # ── C: the consequence. Two eyes, overlaid.
    style_axes(c, "Detected eye at full swing", "time (ps)", "V(det)  (V)")
    for w, col, lab in ((wl, C_LIN, "C$_j$ fixed"), (wf, C_FULL, "C$_j$(V)")):
        for i, (tt, yy) in enumerate(fold_eye(w["t"], w["det"])):
            c.plot(tt, yy, color=col, linewidth=0.6, alpha=0.55,
                   label=lab if i == 0 else None)
    c.legend(frameon=False, fontsize=8.5, labelcolor=INK_2, loc="lower right")

    # ── D: the headline. Where the two models part company.
    style_axes(d, "Eye opening vs drive swing", "differential swing (V$_\\pi$)",
               "worst-case eye opening (mV)")
    sw = np.array(res["swing"])
    full = np.array([m["opening"] for m in res["full"]]) * 1e3
    lin = np.array([m["opening"] for m in res["lin"]]) * 1e3
    d.plot(sw, lin, color=C_LIN, linewidth=2.0, marker="o", markersize=5,
           markeredgecolor=SURFACE, markeredgewidth=1.2, label="C$_j$ fixed (linearised)")
    d.plot(sw, full, color=C_FULL, linewidth=2.0, marker="o", markersize=5,
           markeredgecolor=SURFACE, markeredgewidth=1.2, label="C$_j$(V) (co-solved)")
    d.legend(frameon=False, fontsize=8.5, labelcolor=INK_2, loc="upper left")
    err = 100.0 * (lin - full) / np.maximum(full, 1e-12)
    k = int(np.argmax(np.abs(err)))
    d.annotate(f"{err[k]:+.0f}%", (sw[k], 0.5 * (lin[k] + full[k])),
               textcoords="offset points", xytext=(10, 0), color=INK, fontsize=9,
               fontweight="semibold")
    d.annotate("agree at small signal —\nthe control", (sw[0], full[0]),
               textcoords="offset points", xytext=(6, 26), color=INK_2, fontsize=8.5)

    fig.suptitle(
        "A fixed junction capacitance answers a different circuit",
        color=INK, fontsize=13, x=0.011, ha="left", fontweight="semibold", y=0.985,
    )
    fig.text(0.011, 0.945,
             f"50 Gb/s PRBS-{PRBS_ORDER}, {R_DRV:.0f} Ω driver, {L_ARM_UM / 1000:.0f} mm arms. "
             "Both runs share the solver, the optics and the pattern; only the "
             "capacitance law differs.",
             color=INK_2, fontsize=9, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=170, facecolor=SURFACE)
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selftest", action="store_true",
                    help="assert the claims instead of plotting")
    ap.add_argument("--out", default=None, help="figure path")
    args = ap.parse_args()

    bits = prbs(PRBS_ORDER, N_SYM)
    res = sweep(bits)

    print(f"\nC_j at bias {V_CENTRE:+.1f} V: {c_j(V_CENTRE) * 1e15:.1f} fF")
    lo, hi = V_CENTRE - 0.5 * V_PI, V_CENTRE + 0.5 * V_PI
    print(f"C_j over a full-swing bit: {c_j(lo) * 1e15:.1f} fF at {lo:+.1f} V "
          f"→ {c_j(hi) * 1e15:.1f} fF at {hi:+.1f} V  ({c_j(hi) / c_j(lo):.2f}×)\n")
    print(f"{'swing/V_pi':>11}  {'C_j(V) (mV)':>12}  {'C_j fixed (mV)':>15}  {'error':>8}"
          f"  {'rise ps':>8}  {'fall ps':>8}")
    err = []
    for frac, f, l, e_f in zip(res["swing"], res["full"], res["lin"], res["edges_full"]):
        e = 100.0 * (l["opening"] - f["opening"]) / abs(f["opening"])
        err.append(e)
        print(f"{frac:>11.2f}  {f['opening'] * 1e3:>12.2f}  "
              f"{l['opening'] * 1e3:>15.2f}  {e:>7.1f}%"
              f"  {e_f[0] * 1e12:>8.2f}  {e_f[1] * 1e12:>8.2f}")
    r_lin, f_lin = res["edges_lin"][-1]
    r_full, f_full = res["edges_full"][-1]
    print(f"\nfull swing: C_j(V) rise {r_full * 1e12:.2f} ps vs fall {f_full * 1e12:.2f} ps "
          f"({100 * (f_full - r_full) / (0.5 * (f_full + r_full)):+.1f}% apart)")
    print(f"            C_j fixed rise {r_lin * 1e12:.2f} ps vs fall {f_lin * 1e12:.2f} ps "
          f"({100 * (f_lin - r_lin) / (0.5 * (f_lin + r_lin)):+.1f}% apart — a fixed pole "
          "has one time constant)")
    if args.selftest:
        assert abs(err[0]) < 2.0, (
            f"small signal must agree: the two models differ by {err[0]:.1f}% at "
            f"{SWING_FRACS[0]}·V_pi, so the comparison is not controlled"
        )
        assert abs(err[-1]) > 5.0, (
            f"full swing must diverge: only {err[-1]:.1f}% apart, so this demo "
            "does not demonstrate anything"
        )
        print("\nselftest ok: models agree at small signal and part at full swing")
        return

    out = args.out or os.path.join(os.path.dirname(__file__),
                                   "large_signal_vs_linearised.png")
    figure(res, bits, out)


if __name__ == "__main__":
    main()
