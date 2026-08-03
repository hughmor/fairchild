#!/usr/bin/env python3
"""Both noise analyses on one receiver, and the trap between them.

`.noise` says how much noise there is per hertz. Transient noise puts it in the
waveform. They agree — but only if you ask them about the same operating point,
and that is the whole point of this example.

**`.noise` linearises about ONE bias.** A modulated link does not have one: a
laser driven above threshold carries RIN and shot noise proportional to its
power, so the `1` rail is far noisier than the `0` rail. Run `.noise` on the
deck as written and you get whatever bias the sources happen to sit at (here,
laser off — thermal noise only), which under-predicts the eye closure by ~17x.
Run it once per rail and it is exact.

That is also why the Q-factor formula has two sigmas in it:

    Q = (mu_1 - mu_0) / (sigma_1 + sigma_0)        BER = erfc(Q/sqrt2)/2

Left panel is the eye — visibly noisier on top than on the bottom, which is
what an RIN-limited link looks like. Right panel is the check: sampled sigma
per rail against the `.noise` integral for that rail, swept over `noisescale`,
which is how you reach a deep-BER number without simulating to 1e-12.

    python3 examples/photonic/noisy_eye_and_ber.py [--selftest]
"""

import argparse
import math
import os
import sys

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))
import fairchild  # noqa: E402

BIT_S = 100e-12                 # 10 Gb/s
STEP_S = 1e-12
N_BITS = 800
V_TH, V_HI = 0.9, 2.4
SLOPE_W_V = 0.6e-3              # W/V — low enough that the eye is noise-limited
RESPONSIVITY = 0.8
R_LOAD = 2.0e3
C_PD = 10e-15                   # tau = 20 ps, so a rail settles well inside a bit
RIN_DB_HZ = -131.0              # a cheap laser: RIN-limited, not thermal-limited

# Run-length-rich, deliberately: only bits whose two predecessors matched are
# usable for a per-rail noise measurement (see `sample_levels`), so a pattern of
# mostly isolated bits would leave ~20 usable samples per rail and a 16 % error
# bar on every sigma. This one leaves ~160 per rail, i.e. ~5 %.
PATTERN = ("0011101111000111001010001111000111110000010011110000111000111110"
           "0111110111001000001111000001110000011100110001110000111110000111"
           "11000001110111110000111000011100" * 6)[:N_BITS]


def netlist(trannoise: bool, seed: int, scale: float, dc_drive=None) -> str:
    """The link. `dc_drive` replaces the pattern with a fixed level, which is
    what `.noise` needs: one rail at a time."""
    noise = (f".options trannoise=1 noiseseed={seed} noisescale={scale}\n"
             if trannoise else "")
    if dc_drive is not None:
        drive = f"Vd  drv 0 DC {dc_drive}"
        return f"""* receiver biased at one rail
{noise}.optical_port opt
{drive}
Xlas opt drv 0 fc_driven_laser slope_w_v={SLOPE_W_V} v_th={V_TH} r_in=1e12 rin_db_hz={RIN_DB_HZ}
Xpd opt det 0 fc_photodetector responsivity={RESPONSIVITY} r_shunt=1Meg i_dark_a=0
Rl  det 0 {R_LOAD}
Cl  det 0 {C_PD}
.end
"""
    pts, t = [], 0.0
    edge = 0.2 * BIT_S
    for bit in PATTERN:
        v = V_HI if bit == "1" else 0.0
        pts += [f"{t + edge:.6e} {v:g}", f"{t + BIT_S:.6e} {v:g}"]
        t += BIT_S
    return f"""* noise-limited 10 Gb/s receiver
{noise}.optical_port opt
Vd  drv 0 PWL(0 0 {' '.join(pts)})
Xlas opt drv 0 fc_driven_laser slope_w_v={SLOPE_W_V} v_th={V_TH} r_in=1e12 rin_db_hz={RIN_DB_HZ}
Xpd opt det 0 fc_photodetector responsivity={RESPONSIVITY} r_shunt=1Meg i_dark_a=0
Rl  det 0 {R_LOAD}
Cl  det 0 {C_PD}
.end
"""


def run_tran(trannoise: bool, seed: int = 1, scale: float = 1.0):
    c = fairchild.Circuit()
    c.load_str(netlist(trannoise, seed, scale))
    r = c.run("tran", step=STEP_S, stop=N_BITS * BIT_S)
    return np.asarray(r.time()), np.asarray(r["V(det)"])


def noise_sigma(v_drive: float) -> float:
    """RMS output noise with the laser held at `v_drive`, from `.noise`.

    The receiver's own pole does the band-limiting, so the integral converges —
    the same reason the time-domain answer does not depend on the timestep.
    """
    c = fairchild.Circuit()
    c.load_str(netlist(False, 1, 1.0, dc_drive=v_drive))
    r = c.run("noise", out="det", src="Vd",
              fstart=1e3, fstop=1e13, points=40, variation="dec")
    return math.sqrt(np.trapezoid(np.asarray(r["onoise"]), np.asarray(r.freq())))


def sample_levels(t, v, settled_only=False):
    """Mid-bit samples, split by the bit that was sent.

    `settled_only` keeps just the bits whose two predecessors matched, which
    matters more than it looks: the rails differ in noise by ~50x here, so a
    zero straight after a one carries the ONE's noise through the receiver's
    RC and swamps the quiet rail's own. An isolated zero is not at the zero
    rail's noise level, and neither simulator nor scope will pretend otherwise.
    """
    skip = 4
    idx = [k for k in range(skip, N_BITS)
           if not settled_only or PATTERN[k - 2:k + 1] in ("000", "111")]
    centres = np.array([(k + 0.6) * BIT_S for k in idx])
    s = np.interp(centres, t, v)
    bits = np.array([PATTERN[k] == "1" for k in idx])
    return s[bits], s[~bits]


SEEDS = (11, 12, 13, 14)


def rail_sigmas(t0, ones0, zeros0, scale: float):
    """Per-rail noise sigma at `scale`, pooled over seeds.

    One transient is ONE realisation of a random process: with ~180 usable
    samples per rail a single run's sigma carries ~5 % of scatter, and the seeds
    above spread over 8 %. Pooling the residuals is both the honest estimator
    and the reason this example can assert a 10 % agreement instead of a 25 %
    one. It is also what you would do on the bench.
    """
    hi, lo = [], []
    for seed in SEEDS:
        _, v = run_tran(True, seed=seed, scale=scale)
        ones, zeros = sample_levels(t0, v, settled_only=True)
        hi.append(ones - ones0)
        lo.append(zeros - zeros0)
    return np.std(np.concatenate(hi)), np.std(np.concatenate(lo))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--png", default="noisy_eye_and_ber.png")
    args = ap.parse_args()

    # `.noise` once per rail — the whole point. The idle-bias number is what you
    # get by running `.noise` on the pattern deck without thinking about it.
    sigma_hi = noise_sigma(V_HI)
    sigma_lo = noise_sigma(0.0)

    # Clean run: the deterministic waveform, for the levels and for subtracting.
    t0, v0 = run_tran(False)
    ones0, zeros0 = sample_levels(t0, v0, settled_only=True)
    separation = ones0.mean() - zeros0.mean()

    scales = np.array([1.0, 2.0, 3.0, 4.0])
    hi, lo = np.array([rail_sigmas(t0, ones0, zeros0, float(s)) for s in scales]).T

    rel_hi = abs(hi[0] - sigma_hi) / sigma_hi
    rel_lo = abs(lo[0] - sigma_lo) / sigma_lo
    # Linearity is measured on the QUIET rail. The loud one is RIN-limited at
    # ~10 % optical-power fluctuation per unit noisescale, so by 4x it swings
    # +/-40 % and the photodiode's square law is no longer small-signal — it
    # comes out ~5 % superlinear, which is physics, not a defect. Worth knowing
    # before extrapolating a BER off a scaled run: the extrapolation is
    # conservative, not optimistic.
    slope = np.polyfit(scales, lo, 1)[0] * scales[0] / lo[0]
    slope_hi = np.polyfit(scales, hi, 1)[0] * scales[0] / hi[0]
    q = separation / (hi + lo)
    ber = 0.5 * np.array([math.erfc(x / math.sqrt(2.0)) for x in q])

    print(f"{'rail':<12}{'.noise sigma':>14}{'.tran sigma':>14}{'error':>9}")
    print(f"{'1 (%.1f V)' % V_HI:<12}{sigma_hi * 1e6:11.1f} uV{hi[0] * 1e6:11.1f} uV"
          f"{100 * rel_hi:+8.1f} %")
    print(f"{'0 (0.0 V)':<12}{sigma_lo * 1e6:11.1f} uV{lo[0] * 1e6:11.1f} uV"
          f"{100 * rel_lo:+8.1f} %")
    print(f"\nlevel dependence  sigma_1 / sigma_0 = {hi[0] / lo[0]:.0f}x — RIN and shot "
          f"follow the optical power,\n  so running .noise at the idle bias predicts "
          f"{sigma_lo * 1e6:.0f} uV for BOTH rails and is "
          f"{hi[0] / sigma_lo:.0f}x low on the 1 rail.\n")
    print(f"eye separation {separation * 1e3:.0f} mV, sampled at 0.6 of the bit\n")
    print(f"{'noisescale':>10}{'sigma_1':>12}{'sigma_0':>10}{'Q':>8}{'BER':>12}")
    for sc, h, l, qq, bb in zip(scales, hi, lo, q, ber):
        print(f"{sc:>10.0f}{h * 1e6:>10.0f} uV{l * 1e6:>8.0f} uV{qq:>8.1f}{bb:>12.1e}")
    print(f"\nsigma is linear in noisescale (quiet rail {slope:.3f}), so Q extrapolates "
          f"as 1/noisescale —\n  which is how you get a 1e-12 BER without simulating "
          f"1e12 bits. The loud rail runs\n  {slope_hi:.3f}: at 4x it swings +/-40 % in "
          f"optical power and the square-law PD stops\n  being small-signal, so the "
          f"extrapolation errs conservative.")

    if args.selftest:
        # ~720 pooled samples per rail: a sigma estimate carries ~2.6 %.
        # 10 % is loose enough for the seed set and far tighter than any
        # factor-of-2 mistake.
        assert rel_hi < 0.10, f"1 rail: .noise vs .tran differ by {100 * rel_hi:.1f} %"
        assert rel_lo < 0.10, f"0 rail: .noise vs .tran differ by {100 * rel_lo:.1f} %"
        assert abs(slope - 1.0) < 0.01, f"quiet rail must be linear in noisescale: {slope:.3f}"
        assert 1.0 <= slope_hi < 1.15, f"loud rail should be mildly superlinear: {slope_hi:.3f}"
        assert hi[0] > 10.0 * lo[0], "the demo should be visibly level-dependent"
        assert 8.0 < q[0] < 40.0, f"Q = {q[0]:.1f}: the demo eye should be open but finite"
        print("\nselftest OK — each rail matches its own .noise integral, "
              "and sigma scales linearly with noisescale")
        return 0

    import matplotlib.pyplot as plt

    # Plot a scaled run: at 1x the fuzz is 3 % of the eye and invisible in
    # print. The right-hand panel is where the 1x number is checked.
    eye_scale = 3.0
    _, v1 = run_tran(True, seed=SEEDS[0], scale=eye_scale)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    for k in range(4, min(N_BITS, 300)):
        m = (t0 >= k * BIT_S) & (t0 < (k + 1) * BIT_S)
        ax.plot((t0[m] - k * BIT_S) / BIT_S, v1[m] * 1e3,
                color="C0", alpha=0.12, lw=0.8)
    ax.axvline(0.6, color="0.4", ls=":", lw=1)
    ax.annotate("sampling instant", (0.59, separation * 1e3 * 0.5),
                fontsize=7, rotation=90, ha="right", va="center", color="0.35")
    ax.set_xlabel("bit phase")
    ax.set_ylabel("V(det)  (mV)")
    ax.set_title(f"10 Gb/s eye at noisescale={eye_scale:.0f} (Q = {q[2]:.0f})\n"
                 "RIN-limited: noise rides the 1 rail, not the 0 rail")
    ax.grid(alpha=0.3)

    ax2.plot(scales, sigma_hi * scales * 1e6, "-", color="C3", lw=5, alpha=0.35,
             label="σ₁ from ∫.noise df at the 1 rail")
    ax2.plot(scales, hi * 1e6, "o--", color="C3", lw=1.3, ms=5,
             label="σ₁ sampled from .tran")
    ax2.plot(scales, sigma_lo * scales * 1e6, "-", color="C0", lw=5, alpha=0.35,
             label="σ₀ from ∫.noise df at the 0 rail")
    ax2.plot(scales, lo * 1e6, "o--", color="C0", lw=1.3, ms=5,
             label="σ₀ sampled from .tran")
    ax2.set_xlabel("noisescale")
    ax2.set_ylabel("output noise σ (µV)")
    ax2.set_yscale("log")
    ax2.set_title("Each rail matches its own operating point")
    ax2.legend(fontsize=7, loc="center right")
    ax2.grid(True, which="both", alpha=0.3)

    out = os.path.join(os.path.dirname(__file__), args.png)
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
