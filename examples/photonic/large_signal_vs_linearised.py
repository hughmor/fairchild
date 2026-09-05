#!/usr/bin/env python3
"""How far a linearised modulator model can be trusted, measured rather than argued.

A frequency-domain tool answers a modulator by measuring one small-signal
response at the bias point and filtering the drive with it. That is orders of
magnitude cheaper than a transient co-solve and it is the right tool for a long
PRBS. The question this example answers is where it stops being right, for a
reverse-biased PN phase shifter whose junction capacitance moves with the very
voltage that drives it:

    C_j(V) = C_j0 / (1 − V/V_bi)^m_j        m_j = 0.5, V_bi = 0.917 V

## The controlled experiment

Two runs of the same solver, the same optics, the same PRBS, the same driver:

  * **C_j(V)** — `m_j = 0.5`, the real junction.
  * **C_j fixed** — `m_j = 0`, `c_j0` pinned per arm to `C_j` at that arm's own
    bias. That is exactly the circuit a small-signal extraction sees, so it is
    exactly what a linearised flow is answering.

Only the capacitance law differs. The MZI's optical transfer is `sin²(Δφ/2)` in
both, so the (large) optical nonlinearity is common to the two and cancels out
of the comparison. What is left is the electrical nonlinearity alone.

What is swept is the **bias**, because that decides how much of the `C_j` curve
one bit crosses. Deep in reverse bias a bit rides the flat part; near zero it
climbs the steep part.

## The measured answer

    bias    most forward V    C_j across a bit    linearised error
    −4.0 V      −2.0 V         284 → 421 fF            +3.5 %
    −3.0 V      −1.0 V         295 → 519 fF            +5.2 %
    −2.2 V      −0.2 V         317 → 680 fF            +9.3 %

So: a fixed capacitance is good to a few per cent for a modulator kept in deep
reverse bias, and the error grows monotonically with how much of the curve the
swing covers, reaching about 10 % at the edge of the reverse-biased regime. The
sign is consistent — the linearised model always reports a *more open* eye than
the device has, so it is optimistic, which is the wrong direction for a margin.

That is a calibration, not an indictment. It is worth having as a number.

## What this example actually caught

The first version of this comparison reported a 52 % divergence. That was not
physics: it was a bug in the solver, which advanced `q = C(v)·v` as the
junction's state instead of `q = ∫C dv`, and so ran the device more than twice
as fast as its own parameters. The demo was the thing that found it. The 3-9 %
above is the answer after that fix.

The lesson generalises: a linearised model and a co-solve disagreeing by a lot
is at least as likely to be a bug in the co-solve as a limitation of the linear
model, and the only way to tell is an absolute anchor. Here the anchor was
charge conservation — see `a_junction_stores_the_integral_of_its_capacitance`.

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

# Full swing throughout: what is swept is the *bias*, because that is what
# decides how much of the C_j curve a bit crosses. Deep reverse bias sits on the
# flat part and the two models agree — that agreement is the control. Closer to
# zero the curve steepens and they part.
#
# Each arm sits at `V_CENTRE ± V_pi/4` and swings `± V_pi/4` on top, so the most
# forward point of the pair is `V_CENTRE + V_pi/2`. The list stops where that
# reaches zero: past it the junction conducts, which is a different failure and
# not the one being measured.
BIASES = [-4.0, -3.5, -3.0, -2.75, -2.5, -2.35, -2.2]
SWING_FRAC = 1.0

# ── palette (validated: see the dataviz reference palette) ───────────────────
C_FULL = "#2a78d6"      # slot 1 — the full co-solve
C_LIN = "#eb6834"       # slot 2 — the linearised capacitance
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_MUTED = "#8a8983"
SURFACE = "#fcfcfb"
GRID = "#e3e2dd"


def c_j(v_pn: float) -> float:
    """The junction capacitance the model actually stamps, at bias `v_pn`.

    Including the linear tangent past `V_bi/2`, which is what keeps `C` finite
    into forward bias. Mirrors `junction_cap_and_charge` in the solver — a plot
    of the power law alone goes to infinity where the model does not.
    """
    v_knee = 0.5 * V_BI
    if v_pn < v_knee:
        return C_J0 / (1.0 - v_pn / V_BI) ** M_J
    c_knee = C_J0 / (1.0 - v_knee / V_BI) ** M_J
    return c_knee + c_knee * M_J / (V_BI - v_knee) * (v_pn - v_knee)


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
    """Both models at every bias, plus the waveforms at the closest approach."""
    global V_CENTRE
    keep = V_CENTRE
    out = {"bias": [], "peak": [], "full": [], "lin": [], "c_span": []}
    try:
        for centre in BIASES:
            V_CENTRE = centre
            a = run(bits, SWING_FRAC, linearised=False)
            b = run(bits, SWING_FRAC, linearised=True)
            out["bias"].append(centre)
            # The closest the junction gets to zero over a bit, which is what
            # decides how much of the C_j curve the swing crosses.
            peak = centre + 0.5 * V_PI
            out["peak"].append(peak)
            out["c_span"].append((c_j(centre - 0.5 * V_PI), c_j(peak)))
            out["full"].append(eye_metrics(a["t"], a["det"], bits))
            out["lin"].append(eye_metrics(b["t"], b["det"], bits))
            if centre == BIASES[-1]:
                out["wave_full"], out["wave_lin"] = a, b
                out["wave_bias"] = centre
            if centre == BIASES[0]:
                out["deep_full"], out["deep_lin"] = a, b
    finally:
        V_CENTRE = keep
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


def draw_eye(ax, wave_full, wave_lin, title, note):
    style_axes(ax, title, "time (ps)", "V(det)  (V)")
    for w, col in ((wave_lin, C_LIN), (wave_full, C_FULL)):
        for tt, yy in fold_eye(w["t"], w["det"]):
            ax.plot(tt, yy, color=col, linewidth=0.7, alpha=0.55)
    # The note sits in the eye's own opening, which is the one place in an eye
    # diagram guaranteed to be empty.
    ax.text(0.5, 0.5, note, transform=ax.transAxes, va="center", ha="center",
            color=INK, fontsize=10, fontweight="semibold")


def figure(res, path):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig, axes = plt.subplots(2, 2, figsize=(11.4, 7.8))
    fig.patch.set_facecolor(SURFACE)
    a, b, c, d = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    deep, shallow = BIASES[0], BIASES[-1]
    err = [100.0 * (l["opening"] - f["opening"]) / abs(f["opening"])
           for f, l in zip(res["full"], res["lin"])]

    # ── A: the cause. The piece of the C_j curve each bit actually crosses.
    style_axes(a, "The cause: the piece of C$_j$(V) one bit crosses",
               "junction voltage (V)", "C$_j$ (fF)")
    v = np.linspace(deep - 0.5 * V_PI - 0.4, 0.2, 500)
    a.plot(v, [c_j(x) * 1e15 for x in v], color=INK_MUTED, linewidth=1.6, zorder=2)
    for centre, col, side in ((deep, C_FULL, -1), (shallow, C_LIN, +1)):
        lo, hi = centre - 0.5 * V_PI, centre + 0.5 * V_PI
        seg = np.linspace(lo, hi, 120)
        a.plot(seg, [c_j(x) * 1e15 for x in seg], color=col, linewidth=3.4, zorder=4,
               solid_capstyle="round")
        for x in (lo, hi):
            a.plot([x], [c_j(x) * 1e15], "o", color=col, markersize=6,
                   markeredgecolor=SURFACE, markeredgewidth=1.4, zorder=5)
        mid = 0.5 * (lo + hi)
        a.annotate(f"{centre:.2f} V bias · {c_j(hi) / c_j(lo):.1f}× across a bit",
                   (mid, c_j(mid) * 1e15), textcoords="offset points",
                   xytext=(-30, 40 if side > 0 else -32), ha="center", color=col,
                   fontsize=9, fontweight="semibold")
    a.set_xlim(v[0], 1.2)
    a.set_ylim(0, c_j(0.2) * 1e15 * 1.12)

    # ── B and C: the consequence, at each end of that piece.
    draw_eye(b, res["deep_full"], res["deep_lin"],
             f"Eye at {deep:.1f} V bias — the flat part",
             f"models agree\nto {err[0]:.0f}%")
    draw_eye(c, res["wave_full"], res["wave_lin"],
             f"Eye at {shallow:.2f} V bias — the steep part",
             f"linearised model\nis {err[-1]:.0f}% optimistic")

    # ── D: the headline. Where the linear answer stops being trustworthy.
    style_axes(d, "Eye opening vs how close a bit comes to zero bias",
               "most forward junction voltage over a bit (V)",
               "worst-case eye opening (mV)")
    peak = np.array(res["peak"])
    full = np.array([m["opening"] for m in res["full"]]) * 1e3
    lin = np.array([m["opening"] for m in res["lin"]]) * 1e3
    d.fill_between(peak, full, lin, color=C_LIN, alpha=0.12, linewidth=0)
    d.plot(peak, lin, color=C_LIN, linewidth=2.0, marker="o", markersize=5,
           markeredgecolor=SURFACE, markeredgewidth=1.2)
    d.plot(peak, full, color=C_FULL, linewidth=2.0, marker="o", markersize=5,
           markeredgecolor=SURFACE, markeredgewidth=1.2)
    for i, ha, dx in ((0, "left", 10), (len(peak) - 1, "right", -10)):
        d.annotate(f"+{err[i]:.0f}%", (peak[i], 0.5 * (lin[i] + full[i])),
                   textcoords="offset points", xytext=(dx, 0), ha=ha,
                   color=INK if i else INK_2,
                   fontsize=10 if i else 9,
                   fontweight="semibold" if i else "normal")

    fig.suptitle(
        "A fixed junction capacitance is fine until a bit reaches the steep part of C$_j$(V)",
        color=INK, fontsize=13, x=0.011, ha="left", fontweight="semibold", y=0.978,
    )
    fig.text(0.011, 0.938,
             f"{1e-9 / BIT_S:.0f} Gb/s PRBS-{PRBS_ORDER}, {R_DRV:.0f} Ω driver, "
             f"{L_ARM_UM / 1000:.0f} mm arms, full V$_\\pi$ swing. Both runs share the "
             "solver, the optics and the pattern; only the capacitance law differs.",
             color=INK_2, fontsize=9, ha="left")
    # One legend for the whole figure: the same two series appear in three of
    # the four panels, and three legend boxes would be three chances to collide
    # with the data.
    fig.legend(
        handles=[
            Line2D([], [], color=C_FULL, linewidth=2.4, label="C$_j$(V) — co-solved"),
            Line2D([], [], color=C_LIN, linewidth=2.4, label="C$_j$ fixed — linearised"),
        ],
        loc="upper left", bbox_to_anchor=(0.009, 0.918), frameon=False,
        fontsize=9.5, labelcolor=INK_2, ncols=2, handlelength=1.8,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.885))
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

    print(f"\n{'bias V':>8} {'peak V':>8} {'C_j span (fF)':>16} {'ratio':>6} "
          f"{'C_j(V) mV':>10} {'fixed mV':>9} {'error':>8}")
    err = []
    for centre, peak, (c_lo, c_hi), f, l in zip(
        res["bias"], res["peak"], res["c_span"], res["full"], res["lin"]
    ):
        e = 100.0 * (l["opening"] - f["opening"]) / abs(f["opening"])
        err.append(e)
        print(f"{centre:>8.2f} {peak:>8.2f} "
              f"{c_lo * 1e15:>7.0f} → {c_hi * 1e15:>6.0f} {c_hi / c_lo:>6.2f} "
              f"{f['opening'] * 1e3:>10.1f} {l['opening'] * 1e3:>9.1f} {e:>7.1f}%")

    if args.selftest:
        # The three claims the example makes, at the values it measured. The
        # thresholds are loose enough to survive a solver change that moves the
        # numbers and tight enough to fail if the effect goes away or inverts.
        assert abs(err[0]) < 5.0, (
            f"deep reverse bias must nearly agree: {err[0]:.1f}% apart at "
            f"{res['bias'][0]} V, so the comparison is not controlled"
        )
        assert err[-1] > 7.0, (
            f"near zero bias must diverge: only {err[-1]:.1f}% apart, so this "
            "example measures nothing"
        )
        assert all(x > 0 for x in err), (
            "the linearised model is optimistic at every bias — a sign flip means "
            f"the comparison lost its control: {err}"
        )
        assert err[-1] > 2.0 * err[0], (
            "the error must grow with how much of the C_j curve a bit crosses; "
            f"got {err[0]:.1f}% then {err[-1]:.1f}%"
        )
        print(f"\nselftest ok: {err[0]:.1f}% on the flat part of C_j(V), "
              f"{err[-1]:.1f}% on the steep part, optimistic throughout")
        return

    out = args.out or os.path.join(os.path.dirname(__file__),
                                   "large_signal_vs_linearised.png")
    figure(res, out)


if __name__ == "__main__":
    main()
