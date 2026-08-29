#!/usr/bin/env python3
"""replicate_deck.py — the two headline measurements, run in the simulator.

    MPLBACKEND=Agg .venv/bin/python experiments/cmim/replicate_deck.py [name ...]

EVERY NUMBER HERE COMES OUT OF A NEWTON SOLVE.  `common/channel.py`'s numpy
model is not used: the channel is `DeckChannel`, which drives
`link/netlists/cmim_link.sp` through `Link.tran` and hands back the receiver's
own waveform.  That is the point of the exercise, and using the fast model for a
single-shot measurement — which is what both of these are — would have been
measuring the wrong thing.

The scale is reduced from the paper's and says so in every result file.  The
paper bifurcates 262,144 spins and multiplies 128 x 128 matrices 500 times; a
transient of that is billions of timesteps.  These use hundreds of spins and
tens of dot products, which is enough to show the behaviour and small enough to
finish.  Nothing else is changed.

    dot_product   accuracy and effective bit precision of the analogue
                  vector-vector product against baud rate, in both the
                  linear-modulator configuration the paper characterises with
                  and the nonlinear one the solver runs in
    bifurcation   spin states against feedback strength, with the machine's own
                  noise switched on
"""
from __future__ import annotations

import json
import pathlib
import sys
import time

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "common"))

import dsp                                   # noqa: E402
from channel import DeckChannel              # noqa: E402
from instruments import AWG_FS, quantise     # noqa: E402

RESULTS = HERE / "results"
RESULTS.mkdir(exist_ok=True)


def save(name, payload):
    (RESULTS / f"{name}.json").write_text(json.dumps(payload, indent=2))
    print(f"  -> results/{name}.json")


def plot(name, draw, size=(10, 3.4)):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=size, dpi=150)
    draw(fig)
    fig.tight_layout()
    fig.savefig(RESULTS / f"{name}.png")
    plt.close(fig)
    print(f"  -> results/{name}.png")


# ── the analogue dot product ────────────────────────────────────────────────
def dot_product(bauds=(4e9, 8e9, 16e9, 32e9, 64e9, 96e9, 106e9, 128e9, 148e9),
                n=32, rows=48, n_pilot=64, oversample=4, seed=0):
    """Accuracy of the accumulated product against baud rate, in the deck.

    Each point is `rows` dot products of length `n`, which is the quantity
    Figure 2(c) plots: its axes are the expected and measured vector-matrix
    product, and its statistics come from 500 sampled 128 x 128 matrices.

    Two configurations, because the paper uses two.  Methods characterises the
    multiplication with "both TFLN MZMs operated in the linear regime", driven
    by "two random vectors, x1 and x2, sampled from a uniform distribution in
    [-1, 1] at 8-bit resolution", against a target of `x1 (*) x2`.  The solver
    instead drives the spin channel across its full half-wave, so its target is
    `w * sin(pi*x/2)`.  The second is the one the machine computes.
    """
    ch = DeckChannel(oversample=oversample)
    rng = np.random.default_rng(seed)
    out = {"config": {"n": n, "rows": rows, "n_pilot": n_pilot,
                      "oversample": oversample,
                      "note": "reduced scale; the paper uses 128x128 and 500 "
                              "matrices"},
           "paper_accuracy_pct": {"4": 98.16, "64": 96.2, "148": 90.7},
           "paper_bits": {"4": 5.03, "32": 4.5, "106": 3.3, "148": 2.79},
           "rows": []}

    for linear in (True, False):
        a_v_x = ch_a = None
        for baud in bauds:
            t0 = time.time()
            x = quantise(rng.uniform(-1, 1, rows * n), 8, 1.0)
            w = quantise(rng.uniform(-1, 1, rows * n), 8, 1.0)
            ix, iw = dsp.interleave(x, w)
            awg_x, awg_w = dsp.tx(ix, iw, baud, 0.2, n_pilot=n_pilot)
            if linear:
                ch.params["a_v_x"] = ch.params.get("a_v_w", 0.4)
            else:
                ch.params.pop("a_v_x", None)
            ch.noise_params = ch.lk.noise_params(**ch.params)
            scope = ch(awg_x, awg_w, baud)
            sym = dsp.rx(scope, baud, n_data=len(ix), n_pilot=n_pilot,
                         clip_sigma=None)
            got = dsp.accumulate_interleaved(sym - np.mean(sym),
                                             np.arange(rows) * n,
                                             dsp.find_parity(sym))
            prod = x * w if linear else w * np.sin(np.pi * x / 2)
            want = prod.reshape(rows, n).sum(axis=1)
            k = min(len(got), len(want))
            a, b = got[:k] - got[:k].mean(), want[:k] - want[:k].mean()
            g = float(a @ b) / float(a @ a) if a @ a else 0.0
            err = np.std(g * a - b) / np.std(b)
            row = {"linear": linear, "baud_GBaud": baud / 1e9,
                   "accuracy_pct": 100 * (1 - err),
                   "bits": float(np.log2(1 / max(err, 1e-9))),
                   "measured": (g * a).tolist(), "expected": b.tolist(),
                   "seconds": time.time() - t0}
            out["rows"].append(row)
            print(f"  {'linear   ' if linear else 'nonlinear'} "
                  f"{baud/1e9:6.0f} GBaud  accuracy {row['accuracy_pct']:6.2f} %"
                  f"  bits {row['bits']:5.2f}   [{row['seconds']:.1f} s]")
    save("deck_dot_product", out)

    def draw(fig):
        ax1 = fig.add_subplot(1, 3, 1)
        ax2 = fig.add_subplot(1, 3, 2)
        ax3 = fig.add_subplot(1, 3, 3)
        for linear, c, lbl in ((True, "C0", "both modulators linear"),
                               (False, "C1", "spin channel nonlinear")):
            r = [x for x in out["rows"] if x["linear"] == linear]
            ax1.plot([x["baud_GBaud"] for x in r],
                     [x["accuracy_pct"] for x in r], "o-", color=c, label=lbl)
            ax2.plot([x["baud_GBaud"] for x in r], [x["bits"] for x in r],
                     "o-", color=c, label=lbl)
        for b, v in out["paper_accuracy_pct"].items():
            ax1.plot(float(b), v, "k*", ms=11)
        ax1.plot([], [], "k*", label="paper")
        for b, v in out["paper_bits"].items():
            ax2.plot(float(b), v, "k*", ms=11)
        ax1.set_xlabel("baud rate (GBaud)")
        ax1.set_ylabel("accuracy (%)")
        ax1.set_title("analogue dot product, in the deck")
        ax1.legend(fontsize=7)
        ax1.grid(alpha=0.3)
        ax2.set_xlabel("baud rate (GBaud)")
        ax2.set_ylabel("effective bit precision")
        ax2.set_title("bit precision")
        ax2.grid(alpha=0.3)
        # Fig 2c style: measured against expected, at the lowest baud rate
        r = [x for x in out["rows"] if not x["linear"]][0]
        ax3.plot(r["expected"], r["measured"], ".", ms=4)
        lim = max(np.abs(r["expected"] + r["measured"]))
        ax3.plot([-lim, lim], [-lim, lim], "k--", lw=1)
        ax3.set_xlabel("expected dot product")
        ax3.set_ylabel("measured")
        ax3.set_title(f"Fig 2c style, {r['baud_GBaud']:.0f} GBaud")
        ax3.grid(alpha=0.3)
    plot("deck_dot_product", draw)


# ── bifurcation, with the machine's own noise ───────────────────────────────
def bifurcation(n_alpha=14, per_alpha=12, iters=20, baud=16e9, n_pilot=64,
                oversample=2, alpha_max=3.0, seed=1, noise=True, block=8):
    """Spin states against feedback strength, iterated through the deck.

    THIS IS THE ONE THE EARLIER RESULT WAS NOT.  `ising.bifurcation` iterates
    `x <- alpha*sin(pi*x/2)` in numpy with no noise of any kind, which is why it
    gave a textbook-clean pitchfork.  Here every iteration is a transient
    through the modulators, the amplifier and the receiver, with
    `.options trannoise=1`, so the spins carry shot noise, amplifier noise and
    signal-spontaneous beat noise, and the states spread the way Figure 2(a)'s
    heatmap does.

    With `J = 0` the weight matrix is `alpha*I`, so the whole alpha sweep fits in
    one transient: spin `i` simply gets its own `alpha`.  That is 30 transients
    for the entire figure rather than 30 per alpha.

    `block` IS LOAD BEARING AND MUST NOT BE 1.  A row of `alpha*I` has one
    non-zero, so pruning the zeros makes every dot product one element long —
    and a block of one is exactly the pairwise cancellation this architecture
    cannot do well, because the two single-channel terms are several times
    larger than the product and nothing else in the block helps cancel them.
    Run it that way and the leaked spin-channel term acts as feedback in its own
    right: the machine bifurcates at every alpha including zero, with no
    threshold at all.  That was the first result out of this function and it
    looked like physics.

    The real machine sends the whole row, zeros included, because the flattened
    weight vector is `N(N+1)` symbols long whether or not the entries are zero.
    So each spin here gets `block` symbols: its own alpha, then `block-1` zero
    weights against other spins.  Those contribute nothing to the product and
    everything to the cancellation.

    The paper uses 262,144 spins at 64 GBaud over 50 iterations.  This uses a
    couple of hundred over 30, which is enough to see both fixed points, the
    threshold, and the spread between them.
    """
    ch = DeckChannel(oversample=oversample, noise=noise)
    rng = np.random.default_rng(seed)
    alphas = np.linspace(0.0, alpha_max, n_alpha)
    a_of_spin = np.repeat(alphas, per_alpha)
    n = len(a_of_spin)
    scale = max(alphas.max(), 1e-9)
    # Block layout: slot 0 carries this spin's own alpha, the rest carry zero
    # weights against other spins.  Fixed for the whole run, so the partners are
    # the same every iteration.
    # Each spin's row is sent as `block` symbols carrying alpha/block against
    # that same spin, so the wanted product accumulates COHERENTLY to alpha*s_i
    # while the receiver's noise, which is independent per symbol, accumulates
    # as its square root.  Signal-to-noise on each spin therefore improves as
    # sqrt(block), and the interleaving still gets `block` pairs to cancel its
    # single-channel terms over.
    #
    # The alternative — one symbol carrying alpha and the rest carrying zero,
    # which is what pruning the zeros out of alpha*I would give — is a
    # single-element dot product.  Its noise is the same and its signal is
    # sqrt(block) smaller, and at this model's per-symbol precision the loop
    # gain then amplifies that noise to saturation: every alpha bifurcates,
    # including zero, and the threshold disappears.  Both encodings compute the
    # same row.  Only one of them can be read back.
    partners = np.repeat(np.arange(n)[:, None], block, axis=1)
    w_block = np.full((n, block), 0.0)
    w_block[:] = (a_of_spin / scale / block)[:, None]
    w_flat = w_block.ravel()

    # ── calibrate the loop gain once, on a DENSE test vector ─────────────
    # The weight vector the sweep itself sends is seven-eighths zeros, so its
    # correlation with the expected value is weak and the SIGN that comes out of
    # it is not reliable — the first attempt calibrated a loop gain of -2.98 and
    # the machine then oscillated instead of bifurcating.  A dense vector gives
    # a strong correlation and an unambiguous sign, and the gain it measures is
    # a property of the hardware, not of the vector.
    def _measure(x_sym, w_sym, starts):
        ix_, iw_ = dsp.interleave(x_sym, w_sym)
        ax_, aw_ = dsp.tx(ix_, iw_, baud, 0.2, n_pilot=n_pilot)
        sc_ = ch(ax_, aw_, baud)
        sy_ = dsp.rx(sc_, baud, n_data=len(ix_), n_pilot=n_pilot,
                     clip_sigma=None)
        return dsp.accumulate_interleaved(sy_ - np.mean(sy_), starts,
                                          dsp.find_parity(sy_))

    cal_rows = max(n // 4, 8)
    cx = rng.uniform(-1, 1, cal_rows * block)
    cw = rng.uniform(-1, 1, cal_rows * block)
    cgot = _measure(np.sin(np.pi * cx / 2), cw, np.arange(cal_rows) * block)
    cwant = (cw * np.sin(np.pi * cx / 2)).reshape(cal_rows, block).sum(axis=1)
    kk = min(len(cgot), len(cwant))
    r_cal = float(np.corrcoef(cgot[:kk], cwant[:kk])[0, 1])
    gain = (float(np.std(cwant[:kk])) / max(float(np.std(cgot[:kk])), 1e-30)
            * np.sign(r_cal) * scale)
    print(f"  loop gain calibrated on a dense vector: {gain:.4g} "
          f"(correlation {r_cal:+.3f})")

    # ── calibrate the spin-channel leak, with every weight at zero ───────
    # A row of alpha*I has ONE non-zero, so the dot product the machine forms
    # for each spin is a single element and the interleaving has almost no block
    # to cancel over.  What survives is a term proportional to the spin's own
    # value, which is feedback: run it uncorrected and the machine bifurcates at
    # every alpha including zero, with no threshold.  That was the first version
    # of this figure and it looked like physics.
    #
    # It is a fixed hardware constant, so it is measured once with the weight
    # channel held at zero and subtracted — which is what a bench does, and what
    # balanced photodetection would do structurally instead.  Methods offers
    # balanced detection as the alternative to "DC filtering and a
    # time-interleaving encoding scheme"; this is the second route with its
    # residual calibrated out.
    s_cal = np.sin(np.pi * np.clip(rng.uniform(-1, 1, n), -1, 1) / 2)
    leak = _measure(s_cal[partners.ravel()], np.zeros(n * block),
                    np.arange(n) * block)
    s_blk = s_cal[partners].sum(axis=1)
    kk = min(len(leak), len(s_blk))
    denom = float(s_blk[:kk] @ s_blk[:kk])
    leak_c = float(leak[:kk] @ s_blk[:kk]) / denom if denom > 0 else 0.0
    print(f"  spin-channel leak coefficient {leak_c:+.4g} "
          f"({100*abs(leak_c)*np.std(s_blk)/max(np.std(leak),1e-30):.0f} % "
          f"of the zero-weight residual)")

    x = rng.uniform(-0.05, 0.05, n)
    hist = np.empty((iters, n))
    t0 = time.time()
    for t in range(iters):
        s = np.sin(np.pi * np.clip(x, -1, 1) / 2)
        x_flat = s[partners.ravel()]
        ix, iw = dsp.interleave(x_flat, w_flat)
        awg_x, awg_w = dsp.tx(ix, iw, baud, 0.2, n_pilot=n_pilot)
        scope = ch(awg_x, awg_w, baud)
        sym = dsp.rx(scope, baud, n_data=len(ix), n_pilot=n_pilot,
                     clip_sigma=None)
        got = dsp.accumulate_interleaved(sym - np.mean(sym),
                                         np.arange(n) * block,
                                         dsp.find_parity(sym))
        got = got - leak_c * s[partners].sum(axis=1)[:len(got)]
        # One loop gain, CALIBRATED ONCE and then held.  It stands for the
        # receiver's transimpedance and the two modulator depths, which are
        # hardware constants and do not belong to the algorithm.
        #
        # Refitting it every iteration by least squares looks equivalent and is
        # not: fitting `want ~ g*got` with a noisy `got` is regression dilution,
        # so `g` comes out systematically small, the loop gain sags below the
        # bifurcation threshold and NOTHING separates at any alpha.  That was
        # the second wrong version of this figure.  A ratio of standard
        # deviations carries the amplitude regardless of how noisy the
        # measurement is.
        want = (w_block * s[partners]).sum(axis=1)
        k = min(len(got), len(want))
        x = np.zeros(n)
        x[:k] = got[:k] * gain
        x = np.clip(x, -1.0, 1.0)
        hist[t] = x
        if t % 10 == 0 or t == iters - 1:
            print(f"  iteration {t:3d}/{iters}  |x| mean {np.mean(np.abs(x)):.4f}"
                  f"   [{time.time()-t0:.0f} s]")

    # "Bifurcated" means pinned at a fixed point, not merely displaced from
    # zero.  With the machine's own noise in the loop a sub-threshold spin
    # wanders, and a loose threshold counts that wandering as a bifurcation:
    # at 0.25 every alpha including zero came out "bifurcated" while the
    # scatter plainly showed a transition around 2/pi.
    frac = np.array([np.mean(np.abs(hist[-1][a_of_spin == a]) > 0.9)
                     for a in alphas])
    mono = np.maximum.accumulate(frac)
    crit = (float(np.interp(0.5, mono, alphas))
            if mono[-1] >= 0.5 > mono[0] else None)
    save("deck_bifurcation", {
        "config": {"n_spins": int(n), "alphas": alphas.tolist(),
                   "per_alpha": per_alpha, "iters": iters,
                   "baud_GBaud": baud / 1e9, "trannoise": noise,
                   "note": "reduced scale; the paper uses 262,144 spins over "
                           "50 iterations at 64 GBaud"},
        "critical_alpha_measured": crit,
        "critical_alpha_closed_form": 2 / np.pi,
        "fraction_bifurcated": frac.tolist(),
        "final_state": hist[-1].tolist(),
        "alpha_of_spin": a_of_spin.tolist(),
        "history": hist.tolist(),
    })
    print(f"  critical alpha {crit}  against 2/pi = {2/np.pi:.4f}")

    def draw(fig):
        a1 = fig.add_subplot(1, 3, 1)
        a1.plot(a_of_spin, hist[-1], ".", ms=3, alpha=0.5)
        a1.axvline(2 / np.pi, color="C3", ls="--", lw=1, label=r"$2/\pi$")
        a1.set_xlabel("feedback strength")
        a1.set_ylabel("spin state")
        a1.set_title(f"bifurcation in the deck, {baud/1e9:.0f} GBaud")
        a1.legend(fontsize=7)
        a2 = fig.add_subplot(1, 3, 2)
        top = a_of_spin >= alphas[-1] - 1e-9
        a2.plot(hist[:, top], lw=0.7, alpha=0.6)
        a2.set_xlabel("iteration")
        a2.set_ylabel("spin state")
        a2.set_title(rf"evolution at $\alpha$ = {alphas[-1]:.2f}")
        a3 = fig.add_subplot(1, 3, 3)
        a3.hist(hist[-1][top], bins=40, color="0.3")
        a3.set_xlabel("final spin state")
        a3.set_title("final histogram")
    plot("deck_bifurcation", draw)


NAMES = {"dot_product": dot_product, "bifurcation": bifurcation}

if __name__ == "__main__":
    for nm in (sys.argv[1:] or list(NAMES)):
        if nm not in NAMES:
            raise SystemExit(f"unknown: {nm}; have {list(NAMES)}")
        print(f"\n=== {nm} (in the simulator) ===")
        t0 = time.time()
        NAMES[nm]()
        print(f"  [{time.time()-t0:.1f} s total]")
