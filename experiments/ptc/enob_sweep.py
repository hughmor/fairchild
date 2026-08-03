"""ENOB vs detected power: does the transient sim land on eq:master_snr?

Sweeps laser power across the floor / shot / ceiling regimes of the universal
curve (manuscript fig:wall_master(b)) and compares the measured MAC residual
against the closed-form ENOB. Noise is isolated by differencing against the
same deck with the noise sources off, so the finite-rise-time ISI term cancels
and what remains is TIA + shot + RIN alone.

    python3 experiments/ptc/enob_sweep.py --png enob_sweep.png
"""
import argparse
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
import ptc  # noqa: E402


def sweep(powers_mW, arch, blocks=10, seed=0):
    rng = np.random.default_rng(seed)
    n = blocks * arch.L
    W = rng.uniform(-1, 1, (arch.N, n))
    X = rng.uniform(-1, 1, (n, arch.S, arch.N))
    rows = []
    for p in powers_mW:
        base = ptc.Hw(p_laser_mW=p, tr_frac=0.02, oversample=400)
        got_q, want = ptc.simulate(arch, base, W, X)
        got_n, _ = ptc.simulate(arch, ptc.Hw(**{**base.__dict__, "noise": True}),
                                W, X)
        i_ph = ptc.calibrate(arch, base, W, X)[0, 0, 0] / arch.L * arch.f_B
        rows.append((i_ph,
                     ptc.enob(got_n - got_q, span=arch.L),
                     ptc.enob_theory(base, arch, i_ph)))
        print(f"  P0 = {p:7.3f} mW  I_ph = {i_ph * 1e6:8.2f} uA"
              f"  ENOB {rows[-1][1]:5.2f} b   theory {rows[-1][2]:5.2f} b",
              flush=True)
    return np.array(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--png", default="enob_sweep.png")
    ap.add_argument("-N", type=int, default=4)
    ap.add_argument("--blocks", type=int, default=10)
    a = ap.parse_args()

    arch = ptc.Arch(N=a.N, S=1, L=a.N, f_B=25e9)
    # Up to the 7 dBm detector cap. Past ~20 mW of source this is not a laser
    # any more, it is the gain a real link would need -- which is exactly the
    # amplifier fairchild does not have yet (see README). Drive it with source
    # power anyway: the receiver model does not care where the photons came
    # from, and only this range exercises the shot knee and the RIN ceiling.
    powers = np.logspace(-2, np.log10(500.0), 12)
    r = sweep(powers, arch, blocks=a.blocks)

    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    fine = np.logspace(np.log10(r[0, 0]), np.log10(r[-1, 0]), 300)
    hw = ptc.Hw()
    ax.plot(fine * 1e6, [ptc.enob_theory(hw, arch, i) for i in fine],
            "-", c="0.35", lw=1.5, label="eq:master_snr (TIA + shot + RIN)")
    ax.plot(r[:, 0] * 1e6, r[:, 1], "o", ms=6, label="fairchild transient")
    ceil = ptc.enob_theory(hw, arch, 1e9)  # RIN-only limit
    ax.axhline(ceil, ls=":", c="C3", lw=1.2)
    ax.text(r[0, 0] * 1e6 * 1.3, ceil + 0.25,
            f"RIN ceiling, {hw.rin_db_hz:.0f} dB/Hz", color="C3", fontsize=8)
    i_star = hw.s_i_tia**2 / (2 * ptc.Q_E)
    ax.axvline(i_star * 1e6, ls="--", c="0.6", lw=1)
    ax.text(i_star * 1e6 * 1.1, -4.4, r"$I^* = A_0/A_1$", fontsize=8, color="0.4")
    # Largest I_ph the passive link can deliver from a 13 dBm source.
    i_cap = r[:, 0].max() * (10 ** (1.3) / 500.0)
    ax.axvline(i_cap * 1e6, ls="-.", c="C2", lw=1)
    ax.text(i_cap * 1e6 * 0.95, -4.4, "13 dBm source,\nno amplifier",
            fontsize=8, color="C2", ha="right")
    ax.set(xscale="log", xlabel=r"$I_{ph}$ at one receiver ($\mu$A)",
           ylabel=f"block ENOB, L = {arch.L} (operand-referenced)",
           title=f"N = {arch.N}, S = 1, $f_B$ = {arch.f_B / 1e9:.0f} GBaud")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(a.png, dpi=140)
    print(f"wrote {a.png}   max |sim - theory| = "
          f"{np.abs(r[:, 1] - r[:, 2]).max():.2f} b")


if __name__ == "__main__":
    main()
