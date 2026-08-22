#!/usr/bin/env python3
"""**AC** — design a modulator's electro-optic bandwidth against a target.

A lumped Mach-Zehnder modulator has one central tradeoff and no way around it:
lengthening the arms buys phase efficiency (`V_π ∝ 1/L`) and spends bandwidth
(`f_3dB ∝ 1/L`), because the junction capacitance grows with the same length
that buys the phase. `examples/photonic/noisy_eye_and_ber.py --sweep` shows
that as a table. This solves it instead.

The objective is a least-squares fit of the measured electro-optic response
against a target — "be flat to 12 GHz" — which is what a filter or link
budget actually specifies:

    L = Σ_k [ |H(f_k)|² / |H(f_0)|² − target(f_k) ]²

`Circuit.ac_adjoint` returns `dL/dp` for every parameter from **one** backward
pass over the whole sweep, so adding frequencies costs a forward solve each and
nothing extra in the gradient. That is the property worth having: a filter fit
has hundreds of frequencies and a handful of parameters.

Two parameters, and they pull against each other:

    l_um      arm length      — up: more phase per volt, less bandwidth
    r_drv     driver impedance — down: more bandwidth, more drive current

    python3 examples/optimization/ac_modulator_bandwidth.py [--selftest]
"""

import argparse
import math
import os
import sys
import warnings

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))
import fairchild  # noqa: E402

# `backward` warns when a partial's two finite-difference step sizes disagree,
# and several here legitimately do. Worth understanding rather than hiding: the
# bar is *relative*, and `∂L/∂l_um` and `∂L/∂Rn.r` are individually near zero —
# the length's effect on a DC-normalised response almost cancels, and the two
# drivers are anti-phase. A relative bar on a near-zero number is large however
# small the absolute error is. What this example actually uses is the *sum* of
# the partials, dominated by the well-resolved ones, and that total is checked
# against a full re-solve below and agrees to 1e-8.
#
# So it is reported once, deliberately, rather than left to print sixty times
# with different numbers — `simplefilter("once")` cannot collapse them because
# it keys on the message text and the numbers differ.
SHAKY = []
CALLS = [0]


def quiet_backward(run, cotangent, params):
    """`run.backward`, collecting the resolution warning instead of printing."""
    CALLS[0] += 1
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        g = run.backward(cotangent, params)
    for w in caught:
        if "could not be resolved" in str(w.message):
            SHAKY.append(str(w.message))
    return g

V_PI_L = 0.012          # V·m — 1.2 V·cm
C_PER_UM = 0.25e-15     # 250 fF/mm, so C_j0 tracks the length
V_BI, M_J = 0.917, 0.5
V_CENTRE = -3.0         # both arms reverse-biased
F_TARGET = 12e9         # the bandwidth being designed to

# Log-spaced, because the response is flat over decades and rolls off over one:
# a linear grid spends its samples where nothing happens.
FREQS = np.logspace(8, 11, 40)


def deck(l_um: float, r_drv: float) -> str:
    """MZM from primitives. `c_j0` follows the length — that coupling is the
    whole tradeoff, and hard-coding it would optimise a modulator that cannot
    be built."""
    c_j0 = C_PER_UM * l_um
    arm = (f"fc_pn_ps_cap l_um={l_um} v_pi_l={V_PI_L} c_j0={c_j0} "
           f"v_bi={V_BI} m_j={M_J} alpha_dB_cm=2.0 pin_at_ref=1")
    return f"""* MZM electro-optic response
.optical_port lin
.optical_port dk
.optical_port a1
.optical_port a2
.optical_port b1
.optical_port b2
.optical_port out
.optical_port ou
Xlas lin fc_cw_laser power_mW=1.0 wavelength_nm=1550
Vp pd 0 DC {V_CENTRE + 1.0} AC 0.5
Vn nd 0 DC {V_CENTRE - 1.0} AC 0.5 180
Rp pd p {r_drv}
Rn nd n {r_drv}
Xc1 lin dk a1 a2 fc_dcoupler kappa_L={math.pi / 4}
Xps1 a1 b1 p 0 {arm}
Xps2 a2 b2 n 0 {arm}
Xc2 b1 b2 out ou fc_dcoupler kappa_L={math.pi / 4}
"""


def response(l_um: float, r_drv: float):
    """`(|H|² normalised to DC, the run)` — the electro-optic response."""
    c = fairchild.Circuit()
    c.load_str(deck(l_um, r_drv))
    run = c.ac_adjoint(
        node="out_re_0",
        fstart=float(FREQS[0]),
        fstop=float(FREQS[-1]),
        points=len(FREQS),
        variation="lin",
        src="Vp",
    )
    y = np.asarray(run.response)
    return y, run


def target_curve(freqs: np.ndarray) -> np.ndarray:
    """A single-pole target at `F_TARGET`, normalised to 1 at DC.

    Not a brick wall: an achievable target is the point. Asking a one-pole
    circuit for a response it cannot have makes the least-squares fit converge
    to a compromise that is no longer interpretable.
    """
    return 1.0 / (1.0 + (freqs / F_TARGET) ** 2)


def loss_and_grad(l_um: float, r_drv: float):
    """`L`, `dL/dl_um`, `dL/dr_drv` — one sweep, one backward pass."""
    y, run = response(l_um, r_drv)
    f = np.asarray(run.freqs)
    y0 = y[0]
    tgt = target_curve(f)
    resid = y / y0 - tgt
    loss = float((resid ** 2).sum())

    # dL/dy_k. The normalisation by y[0] couples every frequency to the first
    # sample, so the chain rule picks up a second term there — dropping it is a
    # subtle and very plausible bug, which is why the selftest checks the total
    # against a re-solve rather than each piece.
    dL_dy = 2.0 * resid / y0
    dL_dy[0] -= float((2.0 * resid * (y / y0 ** 2)).sum())

    # **Two things the adjoint cannot know, and both are the designer's.**
    #
    # 1. The arms and the drivers move *together*: a push-pull modulator is
    #    symmetric, and optimising one arm of a pair designs something nobody
    #    builds. So ask for every instance's partial and sum them, rather than
    #    differentiating one and multiplying — the sum is exact whether or not
    #    the two are really identical.
    #
    # 2. `c_j0` is not independent of `l_um` here: `deck()` sets
    #    `c_j0 = C_PER_UM·l_um`, because a longer junction *is* a bigger one and
    #    a design that pretended otherwise would optimise a modulator that
    #    cannot be fabricated. The adjoint honours the netlist, where those are
    #    two separate numbers, so it returns two separate partials. The total
    #    derivative is the chain rule over the design rule:
    #
    #        dL/dL_arm = ∂L/∂l_um + (dc_j0/dl_um)·∂L/∂c_j0
    #
    #    Getting this wrong is quiet: `∂L/∂l_um` alone is ~0 here, because the
    #    optical length scales phase and loss at every frequency alike and the
    #    normalisation divides it straight back out. An optimiser handed that
    #    would conclude length does not matter.
    # `l_um` gets an explicit finite-difference step. The automatic choice is
    # `∛ε·|p|` ≈ 24 nm here, and at ~17 rad/µm of propagation phase that moves
    # the residual far less than the objective's own scale — so the difference
    # is roundoff and `backward` warns it could not resolve it. 0.1 µm is large
    # enough to register and far below any curvature. Without the tuple form
    # there would be no way to act on the warning, which is what makes it worth
    # having rather than suppressing.
    g = quiet_backward(
        run,
        dL_dy,
        [("Xps1.l_um", 0.1), ("Xps2.l_um", 0.1),
         "Xps1.c_j0", "Xps2.c_j0", "Rp.r", "Rn.r"],
    )
    d_len = (g[0] + g[1]) + C_PER_UM * (g[2] + g[3])
    d_r = g[4] + g[5]
    return loss, float(d_len), float(d_r)


def f3db(y: np.ndarray, freqs: np.ndarray) -> float:
    db = 10.0 * np.log10(y / y[0])
    k = int(np.argmax(db < -3.0))
    if k == 0:
        return float("nan")
    return float(np.interp(-3.0, [db[k], db[k - 1]], [freqs[k], freqs[k - 1]]))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--png", default="ac_modulator_bandwidth.png")
    args = ap.parse_args()

    l_um, r_drv = 4000.0, 50.0
    y0, run0 = response(l_um, r_drv)
    f = np.asarray(run0.freqs)
    print(f"start:  L = {l_um:.0f} um, R = {r_drv:.1f} ohm  ->  "
          f"f_3dB = {f3db(y0, f) / 1e9:5.2f} GHz, V_pi = {V_PI_L / (l_um * 1e-6):.2f} V/arm")

    # Gradient check before trusting any of it.
    loss0, gl, gr = loss_and_grad(l_um, r_drv)
    for name, val, g, d in (("l_um", l_um, gl, 1.0), ("r_drv", r_drv, gr, 1e-3)):
        lp = loss_and_grad(l_um + d, r_drv)[0] if name == "l_um" else loss_and_grad(l_um, r_drv + d)[0]
        lm = loss_and_grad(l_um - d, r_drv)[0] if name == "l_um" else loss_and_grad(l_um, r_drv - d)[0]
        fd = (lp - lm) / (2 * d)
        print(f"  dL/d{name:<6} adjoint {g:+.6e}   re-solve {fd:+.6e}   "
              f"rel {abs(g - fd) / max(abs(fd), 1e-30):.2e}")

    # Plain gradient descent with per-parameter scaling: `l_um` and `r_drv`
    # differ by three orders in magnitude, so one step size cannot serve both.
    scale = np.array([l_um, r_drv])
    p = np.array([l_um, r_drv])
    rate, hist = 0.06, []
    for _ in range(60):
        loss, a, b = loss_and_grad(float(p[0]), float(p[1]))
        hist.append((p.copy(), loss))
        step = rate * np.array([a, b]) * scale ** 2
        p = np.clip(p - step, [500.0, 10.0], [8000.0, 200.0])
    l_star, r_star = float(p[0]), float(p[1])
    y_star, run_star = response(l_star, r_star)
    loss_star = loss_and_grad(l_star, r_star)[0]

    print(f"\nfinal:  L = {l_star:.0f} um, R = {r_star:.1f} ohm  ->  "
          f"f_3dB = {f3db(y_star, f) / 1e9:5.2f} GHz, "
          f"V_pi = {V_PI_L / (l_star * 1e-6):.2f} V/arm")
    print(f"  loss {hist[0][1]:.5f} -> {loss_star:.5f} in {len(hist)} sweeps")
    print(f"  target was {F_TARGET / 1e9:.0f} GHz")
    if SHAKY:
        print(f"\n  ({len(SHAKY)} of {CALLS[0]} backward passes reported a partial "
              f"they could not\n   resolve well. Expected here — see the note at the top "
              f"of this file — and the\n   total is checked against a re-solve above.)")

    if args.selftest:
        assert loss_star < 0.35 * hist[0][1], \
            f"loss barely moved: {hist[0][1]:.4f} -> {loss_star:.4f}"
        got = f3db(y_star, f)
        assert abs(got - F_TARGET) / F_TARGET < 0.25, \
            f"f_3dB landed at {got/1e9:.2f} GHz against a {F_TARGET/1e9:.0f} GHz target"
        # The gradient is the whole point, so check it — but **not at the
        # optimum**, where it is zero and agreeing about zero proves nothing.
        # The start point is where it is large and a sign or factor error shows.
        d = 1.0
        fd = (loss_and_grad(4000.0 + d, 50.0)[0] - loss_and_grad(4000.0 - d, 50.0)[0]) / (2 * d)
        ga = loss_and_grad(4000.0, 50.0)[1]
        assert abs(ga - fd) / max(abs(fd), 1e-12) < 1e-3, \
            f"dL/dL_arm disagrees with a re-solve: {ga:.4e} vs {fd:.4e}"
        print("\nselftest OK — target bandwidth reached, gradient matches a re-solve")
        return 0

    import matplotlib.pyplot as plt

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    ax.semilogx(f, 10 * np.log10(y0 / y0[0]), color="0.6", lw=1.6,
                label=f"start — {f3db(y0, f)/1e9:.1f} GHz")
    ax.semilogx(f, 10 * np.log10(y_star / y_star[0]), color="C0", lw=2.0,
                label=f"optimised — {f3db(y_star, f)/1e9:.1f} GHz")
    ax.semilogx(f, 10 * np.log10(target_curve(f)), color="C3", ls="--", lw=1.4,
                label=f"target — {F_TARGET/1e9:.0f} GHz")
    ax.axhline(-3, color="0.7", ls=":", lw=1)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("electro-optic response (dB)")
    ax.set_ylim(-25, 3)
    ax.set_title("AC: fit a modulator's response to a target")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.25)

    # The loss, not the two parameter tracks. Plotting length and resistance on
    # twin axes looked like a single curve here — they happen to fall to almost
    # the same *fraction* of their start (0.69 and 0.70), so the twin scaling
    # superimposed them and the figure read as one parameter moving.
    losses = np.array([h[1] for h in hist])
    ax2.semilogy(np.maximum(losses, 1e-12), color="C1", lw=1.8)
    ax2.set_xlabel("iteration")
    ax2.set_ylabel("loss   Σ (|H|²−target)²")
    ax2.set_title("Both parameters descend together,\n"
                  "one backward pass per iteration")
    ax2.grid(True, which="both", alpha=0.25)
    ax2.annotate(f"{hist[0][0][0]:.0f} µm, {hist[0][0][1]:.0f} Ω", xy=(0, losses[0]),
                 xytext=(4, -14), textcoords="offset points", fontsize=8, color="0.35")
    ax2.annotate(f"{l_star:.0f} µm, {r_star:.0f} Ω",
                 xy=(len(losses) - 1, max(losses[-1], 1e-12)),
                 xytext=(-6, 12), textcoords="offset points", fontsize=8,
                 color="0.35", ha="right")

    out = os.path.join(os.path.dirname(__file__), args.png)
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
