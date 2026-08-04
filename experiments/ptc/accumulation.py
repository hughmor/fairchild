"""Experiment E — how ENOB depends on the accumulation depth L.

The manuscript claims (sec:receiver_model, and the equation flagged for triple
checking at eq:adc_power) that an L-symbol analog accumulation:

  * leaves the *per-MAC* SNR unchanged, because signal and white-noise variance
    both scale by L, and
  * gains the *output-referred* resolution 0.5*log2(L), because the full-scale
    range grows as L while the noise grows as sqrt(L),

so converter power amortizes only as 1/sqrt(L): the readout rate falls by L but
the converter must supply 0.5*log2(L) more bits.

Testing that needs one transient, not several. The integrator's charge is
additive, so a single run blocked at L = 1 gives per-symbol charges that can be
re-summed into any L afterwards -- identical physics, identical noise
realisation, only the bookkeeping changes. Anything that moves is then the
claim, not the seed.

    python3 experiments/ptc/accumulation.py --png results/accumulation.png
"""
import argparse
import pathlib
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import ptc  # noqa: E402


def per_symbol_residual(arch, hw, W, X, seed=0):
    """MAC residual for every single symbol, shape (n_sym, N, S, N)."""
    deck = ptc.build(arch, hw, W, X, seed)
    charge = ptc.readout(deck, ptc.run(deck))
    cal = ptc.calibrate(arch, hw, W, X, seed)
    return ptc.decode(deck, charge, cal, W, X) - ptc.ideal(deck, W, X)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--png", default="results/accumulation.png")
    ap.add_argument("-N", type=int, default=4)
    ap.add_argument("--symbols", type=int, default=512)
    ap.add_argument("--p-mW", type=float, default=10.0)
    a = ap.parse_args()

    # L = 1 so `readout` hands back one charge per symbol; the L sweep is done
    # by re-summing those afterwards.
    arch = ptc.Arch(N=a.N, S=1, L=1, f_B=25e9)
    hw = ptc.Hw(p_laser_mW=a.p_mW, tr_frac=0.1, oversample=40)
    rng = np.random.default_rng(7)
    W = rng.uniform(-1, 1, (arch.N, a.symbols))
    X = rng.uniform(-1, 1, (a.symbols, arch.S, arch.N))

    print(f"N = {arch.N}, {a.symbols} symbols, P0 = {a.p_mW} mW ...", flush=True)
    noisy = per_symbol_residual(arch, ptc.Hw(**{**hw.__dict__, "noise": True}),
                                W, X)
    quiet = per_symbol_residual(arch, hw, W, X)
    resid = noisy - quiet  # ISI and the MZM transfer cancel; noise remains
    i_ph = ptc.calibrate(arch, hw, W, X)[0, 0, 0] * arch.f_B
    print(f"I_ph = {i_ph * 1e6:.1f} uA")

    Ls = [1, 2, 4, 8, 16, 32]
    rows = []
    for L in Ls:
        n = (a.symbols // L) * L
        blocks = resid[:n].reshape(-1, L, *resid.shape[1:]).sum(axis=1)
        # Output-referred: full-scale L against the block's own noise.
        # Per-MAC: the same noise shared out over L MACs, which for a white
        # process means dividing by sqrt(L), not by L.
        rows.append((L,
                     ptc.enob(blocks, span=L),
                     ptc.enob(blocks / np.sqrt(L), span=1),
                     blocks.size))
        print(f"  L = {L:<3} output-referred {rows[-1][1]:5.2f} b   "
              f"per-MAC {rows[-1][2]:5.2f} b   ({rows[-1][3]} blocks x rx)")
    r = np.array(rows, dtype=float)

    per_mac = r[0, 1]
    drift = np.abs(r[:, 2] - per_mac).max()
    slope = np.polyfit(np.log2(r[:, 0]), r[:, 1], 1)[0]
    print(f"per-MAC ENOB drift over L = 1..{Ls[-1]}: {drift:.3f} b")
    print(f"output-referred slope: {slope:.3f} b per octave of L (expect 0.5)")

    fig, ax = plt.subplots(figsize=(6.4, 4.3))
    ax.plot(r[:, 0], r[:, 1], "o-", label="output-referred (block of $L$ MACs)")
    ax.plot(r[:, 0], r[:, 2], "s-", label="per-MAC")
    ax.plot(r[:, 0], per_mac + 0.5 * np.log2(r[:, 0]), ":", c="0.4",
            label=r"$ENOB_1 + \frac{1}{2}\log_2 L$")
    ax.axhline(per_mac, ls="--", c="0.7", lw=1)
    ax.set(xscale="log", xticks=Ls, xlabel="accumulation depth $L$ (symbols)",
           ylabel="ENOB (operand-referenced)",
           title=f"N = {a.N}, S = 1, $f_B$ = 25 GBaud, "
                 rf"$I_{{ph}}$ = {i_ph * 1e6:.0f} $\mu$A")
    ax.set_xticklabels([str(x) for x in Ls])
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    out = pathlib.Path(a.png)
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
