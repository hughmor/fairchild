#!/usr/bin/env python3
"""**DC** — trim a Mach-Zehnder to a target split with a thermal heater.

The job a bias controller does. An MZI comes off the wafer with an arbitrary
phase offset — arm-length mismatch, film thickness, temperature — and a heater
in one arm trims it back. "Hold the through port at quadrature" is the setpoint
a coherent receiver, a switch, or a modulator bias loop actually servos to.

It is a DC problem because a thermo-optic heater is: the thermal time constant
is microseconds against a signal in picoseconds, so the bias loop sees only the
operating point.

    P_bar(V) = P_in·cos²(φ₀/2 + π·V²/(R·P_π)/2)

One parameter, one probe, and the gradient comes from
`Circuit.dc_adjoint` — one solve gives both the value and `dP/dV`, so a
Newton step costs the same as an evaluation.

    python3 examples/optimization/dc_mzi_bias_trim.py [--selftest]
"""

import argparse
import math
import os
import sys

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))
import fairchild  # noqa: E402

P_IN_MW = 1.0
R_HEATER = 500.0        # Ω
P_PI_W = 20e-3          # heater power for π
ARM_UM = 1000.0
# A deliberate arm-length mismatch: this is the offset the trim exists to undo.
# 40 nm of extra path at n_eff ≈ 2.445 is a fraction of a wavelength, which is
# all it takes — and is far tighter than any real fab tolerance.
MISMATCH_UM = 0.04

DECK = f"""* MZI with a thermal trim in one arm
.optical_port lin
.optical_port dk
.optical_port a1
.optical_port a2
.optical_port b1
.optical_port b2
.optical_port bar
.optical_port cross
Xlas lin fc_cw_laser power_mW={P_IN_MW}
Xc1  lin dk a1 a2 fc_dcoupler kappa_L={math.pi / 4}
Xheat a1 b1 htp 0 fc_thermal_ps r_heater={R_HEATER} p_pi={P_PI_W} l_um={ARM_UM}
Xref  a2 b2       fc_waveguide l_um={ARM_UM + MISMATCH_UM} alpha_dB_cm=0
Xc2  b1 b2 bar cross fc_dcoupler kappa_L={math.pi / 4}
Vh   htp 0 DC 0
.op
"""


def probe(v_heater: float):
    """`(P_bar, dP_bar/dV)` at one heater voltage, from a single solve."""
    c = fairchild.Circuit()
    c.load_str(DECK)
    r = c.dc_adjoint(
        probes={"p": ("power", "bar", 0)},
        wrt=["Vh.dc"],
        params={"Vh.dc": v_heater},
    )
    return r.values["p"], float(r.grad["p"][0])


def solve(target_w: float, v0: float = 2.0, iters: int = 40):
    """Newton on `P_bar(V) − target = 0`, which is what the gradient buys.

    A bias servo is a root-find, not a minimisation, so the adjoint's `dP/dV` is
    the Newton denominator directly and each iteration costs one solve.

    **The start point matters, and that is the physics, not a weakness of the
    method.** An MZI transfer is a fringe: `P(V)` is periodic in `V²`, so the
    setpoint has infinitely many solutions and `dP/dV` is zero at every peak
    between them. Start on a peak and Newton has no direction to move; start on
    the wrong flank and it converges to a different fringe than intended. Real
    modulator bias controllers have exactly this problem and solve it the same
    way — bring the servo up on a known flank, then hold. The right-hand panel
    plots `dP/dV` so the dead points are visible rather than implied.
    """
    v, history = v0, []
    for _ in range(iters):
        p, dp = probe(v)
        history.append((v, p))
        err = p - target_w
        if abs(err) < 1e-12:
            break
        if abs(dp) < 1e-6:
            # On a peak: no gradient information at all. Nudge off it rather
            # than dividing by ~0, and say so.
            v += 0.05
            continue
        v -= max(-0.4, min(0.4, err / dp))
        v = max(0.0, v)
    return v, history


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--png", default="dc_mzi_bias_trim.png")
    args = ap.parse_args()

    p_in = P_IN_MW * 1e-3
    quad = 0.5 * p_in

    # The uncorrected part: what the mismatch alone does.
    p0, _ = probe(0.0)
    print(f"untrimmed (V=0):  P_bar = {p0 * 1e3:.4f} mW  "
          f"({100 * p0 / p_in:.1f} % of input)")

    v_star, hist = solve(quad)
    p_star, dp_star = probe(v_star)
    print(f"trimmed to quadrature: V = {v_star:.5f} V "
          f"({v_star ** 2 / R_HEATER * 1e3:.3f} mW of heater power)")
    print(f"  P_bar = {p_star * 1e3:.6f} mW, target {quad * 1e3:.6f} mW, "
          f"error {abs(p_star - quad) / quad * 100:.2e} %")
    print(f"  converged in {len(hist)} solves")
    print(f"  dP/dV at the setpoint = {dp_star * 1e3:.4f} mW/V — the servo's loop gain")

    # Gradient check: the whole example rests on it being right.
    d = 1e-6
    fd = (probe(v_star + d)[0] - probe(v_star - d)[0]) / (2 * d)
    rel = abs(dp_star - fd) / abs(fd)
    print(f"  gradient vs a full re-solve: {rel:.2e} relative")

    if args.selftest:
        assert abs(p_star - quad) / quad < 1e-6, "did not reach quadrature"
        assert rel < 1e-4, f"gradient disagrees with a re-solve by {rel:.2e}"
        assert len(hist) < 25, f"took {len(hist)} solves; Newton should be faster"
        assert abs(dp_star) > 1e-4, "quadrature should be the steepest point, not a null"
        print("\nselftest OK — quadrature reached, gradient matches a re-solve")
        return 0

    import matplotlib.pyplot as plt

    sweep = np.linspace(0.0, 3.6, 121)
    curve = np.array([probe(v)[0] for v in sweep])
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.0), constrained_layout=True)

    ax.plot(sweep, curve * 1e3, color="C0", lw=1.8)
    ax.axhline(quad * 1e3, color="0.5", ls=":", lw=1)
    hv, hp = np.array(hist).T
    ax.plot(hv, hp * 1e3, "o-", color="C3", ms=4, lw=0.9, alpha=0.85,
            label=f"Newton path ({len(hist)} solves)")
    ax.plot([v_star], [p_star * 1e3], "*", color="C3", ms=15, label="setpoint")
    ax.set_xlabel("heater drive (V)")
    ax.set_ylabel("P(bar)  (mW)")
    ax.set_title("DC: trim an MZI to quadrature")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    grads = np.array([probe(v)[1] for v in sweep])
    ax2.plot(sweep, grads * 1e3, color="C2", lw=1.8, label="adjoint dP/dV")
    ax2.axhline(0, color="0.6", lw=0.8)
    ax2.axvline(v_star, color="C3", ls="--", lw=1, label="setpoint")
    ax2.set_xlabel("heater drive (V)")
    ax2.set_ylabel("dP(bar)/dV  (mW/V)")
    ax2.set_title("The gradient, one solve per point\n"
                  "zero at every fringe turn — no servo authority there")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)

    out = os.path.join(os.path.dirname(__file__), args.png)
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
