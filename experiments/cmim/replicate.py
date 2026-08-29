#!/usr/bin/env python3
"""replicate.py — work through the paper's results and record what we get.

    MPLBACKEND=Agg .venv/bin/python experiments/cmim/replicate.py [figure ...]

With no argument it runs every figure that currently reproduces.  Results go to
`results/<name>.json` and `results/<name>.png`.  The JSON carries the paper's
own number beside ours for every row, so a regression is visible without
opening a plot.

WHAT RUNS AND WHAT DOES NOT.  See ../README.md for the full status.  In short:
the bifurcation, the lattice scaling, number partitioning and the bit-precision
comparison all run and land near the paper.  The end-to-end analogue
matrix-vector multiplication — Figures 2(c) to 2(e) — does not, and the reason
is a specific under-specified piece of the paper rather than a tuning problem.
`fig2_mvm` runs anyway and prints what it gets, because a number that disagrees
is more useful than a figure that is missing.
"""
from __future__ import annotations

import json
import pathlib
import sys
import time

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "common"))

import ising          # noqa: E402
import problems       # noqa: E402
from instruments import time_to_solution  # noqa: E402

RESULTS = HERE / "results"
RESULTS.mkdir(exist_ok=True)


def save(name: str, payload: dict):
    (RESULTS / f"{name}.json").write_text(json.dumps(payload, indent=2))
    print(f"  -> results/{name}.json")


def plot(name: str, draw):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(9, 3.2), dpi=140)
    draw(fig)
    fig.tight_layout()
    fig.savefig(RESULTS / f"{name}.png")
    plt.close(fig)
    print(f"  -> results/{name}.png")


# ── Figure 2(a), 2(b): bifurcation ──────────────────────────────────────────
def fig2_bifurcation():
    """Feedback strength against final spin state, J = 0.

    Figure 2(a) sweeps alpha with 262,144 uncoupled spins at 64 GBaud and shows
    one fixed point below a critical alpha_0 and two above.  Figure 2(b) shows
    50 iterations at alpha = 3.5.

    The paper gives no closed form for alpha_0.  There is one: with J = 0 the
    update is x <- alpha*sin(pi*x/2), whose fixed point at zero loses stability
    when alpha*pi/2 > 1, so alpha_0 = 2/pi = 0.6366.
    """
    alphas = np.linspace(0.0, 3.5, 71)
    n = 16384
    finals = np.array([ising.bifurcation(a, n, 50, seed=1) for a in alphas])
    frac = (np.abs(finals) > 0.05).mean(axis=1)
    crit = float(np.interp(0.5, frac, alphas))
    # Small initial amplitudes, so the separation is visible over the fifty
    # iterations Figure 2(b) plots rather than complete on the first.
    evo = np.array([ising.bifurcation(3.5, 4096, k, x0_scale=0.02, seed=2)
                    for k in range(1, 51)])

    save("fig2_bifurcation", {
        "critical_alpha_measured": crit,
        "critical_alpha_closed_form": 2 / np.pi,
        "paper": "Fig 2a: one fixed point below alpha_0, two above; "
                 "worked example at alpha = 3.5, 50 iterations",
        "alphas": alphas.tolist(),
        "fraction_bifurcated": frac.tolist(),
    })
    print(f"  critical alpha {crit:.4f} against 2/pi = {2/np.pi:.4f}")

    def draw(fig):
        a1 = fig.add_subplot(1, 3, 1)
        a1.hist2d(np.repeat(alphas, n), finals.ravel(), bins=[71, 80],
                  cmap="magma")
        a1.axvline(2 / np.pi, color="w", ls="--", lw=1)
        a1.set_xlabel("feedback strength")
        a1.set_ylabel("spin state")
        a1.set_title("Fig 2a: bifurcation")
        a2 = fig.add_subplot(1, 3, 2)
        a2.imshow(evo.T, aspect="auto", cmap="magma", origin="lower",
                  extent=[1, 50, 0, 4096])
        a2.set_xlabel("iteration")
        a2.set_ylabel("spin")
        a2.set_title(r"Fig 2b: $\alpha = 3.5$")
        a3 = fig.add_subplot(1, 3, 3)
        a3.hist(evo[-1], bins=60, color="0.3")
        a3.set_xlabel("final spin state")
        a3.set_title("Fig 2b: histogram")
    plot("fig2_bifurcation", draw)


# ── Figure 2(c)-(e): analogue matrix-vector multiplication ──────────────────
def fig2_mvm():
    """Accuracy and effective bit precision of the analogue product.

    THIS ONE DOES NOT REPRODUCE THE PAPER.  It is here so the disagreement is on
    the record with a number attached.  See ../README.md, "What does not work
    yet", for the diagnosis: the time-interleaving scheme is a reconstruction,
    it is algebraically exact but numerically ill conditioned, and it is the
    piece the paper defers to its earlier work for.
    """
    from channel import Channel
    ch = Channel.from_link(seed=0)
    rows = []
    paper = {4e9: 98.16, 64e9: 96.2, 148e9: 90.7}
    paper_bits = {4e9: 5.03, 32e9: 4.5, 106e9: 3.3, 148e9: 2.79}
    for baud in (4e9, 16e9, 32e9, 64e9, 96e9, 106e9, 128e9, 148e9):
        r = ising.calibrate_surrogate(ch, baud=baud, n=128, trials=4)
        rows.append({"baud_GBaud": baud / 1e9,
                     "accuracy_pct": 100 * r["accuracy"],
                     "accuracy_sd_pct": 100 * r["accuracy_sd"],
                     "bits": r["bits"], "bits_sd": r["bits_sd"],
                     "paper_accuracy_pct": paper.get(baud),
                     "paper_bits": paper_bits.get(baud)})
        print(f"  {baud/1e9:6.0f} GBaud  accuracy {100*r['accuracy']:6.2f} %"
              f"  bits {r['bits']:.2f}"
              + (f"   paper {paper[baud]:.2f} %" if baud in paper else ""))
    save("fig2_mvm", {"status": "DOES NOT REPRODUCE — see README",
                      "rows": rows})


# ── Figure 3(a), 3(b) and S7: the square lattice ────────────────────────────
def fig3_lattice(sizes=(10, 14, 20, 30, 40, 50), iters=1000, trials=5):
    """Solution quality against problem size, at the paper's alpha and beta.

    Figure 3(b) runs 100 to 41,209 nodes capped at 1,000 iterations and reports
    97.2 % of the ground state at the largest.  Figure 3(a) reaches the ground
    state of the 101x101 lattice at about iteration 671.

    The sizes here stop well short of 203x203.  That is a runtime choice, not a
    limit: 41,209 spins with 205,233 couplings is a minute per trial in this
    solver, and the shape of the curve is already clear by 2,500 nodes.  Pass
    larger sizes on the command line if you want them.
    """
    out = []
    for L in sizes:
        J, h = problems.square_lattice(L, L)
        gs = problems.lattice_ground_energy(L, L)
        t0 = time.time()
        q, first = [], []
        for s in range(trials):
            r = ising.solve(J, h, alpha=0.86, beta=1.0, iters=iters,
                            mode="surrogate", bits=3.3, seed=s)
            q.append(100 * r.best_energy / gs)
            hit = np.where(r.energy <= gs)[0]
            first.append(int(hit[0]) if len(hit) else None)
        hits = [f for f in first if f is not None]
        out.append({"L": L, "nodes": L * L, "couplings": int(J.nnz),
                    "quality_pct": float(np.mean(q)),
                    "quality_sd_pct": float(np.std(q)),
                    "best_pct": float(np.max(q)),
                    "hit_rate": len(hits) / trials,
                    "iters_to_ground": float(np.mean(hits)) if hits else None,
                    "seconds": time.time() - t0})
        o = out[-1]
        print(f"  {L:3d}x{L:<3d} {L*L:6d} nodes  quality {o['quality_pct']:6.2f} %"
              f"  hit {o['hit_rate']:.0%}"
              + (f"  at iteration {o['iters_to_ground']:.0f}"
                 if o["iters_to_ground"] else "")
              + f"   [{o['seconds']:.1f} s]")
    save("fig3_lattice", {
        "paper": "Fig 3b: 97.2 % at 41,209 nodes, 1,000 iterations; "
                 "Fig 3a: ground state of 101x101 at ~671 iterations",
        "hyperparameters": {"alpha": 0.86, "beta": 1.0, "bits": 3.3},
        "rows": out})

    def draw(fig):
        a1 = fig.add_subplot(1, 2, 1)
        n = [o["nodes"] for o in out]
        a1.errorbar(n, [o["quality_pct"] for o in out],
                    yerr=[o["quality_sd_pct"] for o in out], marker="o")
        a1.axhline(97.2, color="C3", ls="--", lw=1, label="paper, 41,209 nodes")
        a1.set_xscale("log")
        a1.set_xlabel("node count")
        a1.set_ylabel("% of ground state")
        a1.set_title("Fig 3b: solution quality")
        a1.legend(fontsize=7)
        a2 = fig.add_subplot(1, 2, 2)
        L = sizes[-1]
        J, h = problems.square_lattice(L, L)
        r = ising.solve(J, h, alpha=0.86, beta=1.0, iters=iters,
                        mode="surrogate", bits=3.3, seed=0)
        a2.plot(r.energy)
        a2.axhline(problems.lattice_ground_energy(L, L), color="C3", ls="--",
                   lw=1, label="ground state")
        a2.set_xlabel("iteration")
        a2.set_ylabel("Ising energy")
        a2.set_title(f"Fig 3a: {L}x{L} convergence")
        a2.legend(fontsize=7)
    plot("fig3_lattice", draw)


# ── Figure 5: number partitioning ───────────────────────────────────────────
def fig5_partition(sizes=(16, 32, 64, 128, 256), trials=10, iters=1000):
    """Partition N integers drawn from {0..16}, and time the solution.

    The paper reports ground-state solutions for every size from 16 to 256 and
    plots the time to solution against N (Figure 5b).  Number partitioning is
    dense: N spins need N*(N+1) symbols, which is what makes 256 the largest
    the AWG memory holds.
    """
    rng = np.random.default_rng(0)
    out = []
    for N in sizes:
        hits, iters_hit, diffs = 0, [], []
        for s in range(trials):
            v = problems.random_partition_instance(N, 16, rng)
            J, h = problems.number_partition(v)
            r = ising.solve(J, h, alpha=1.0, beta=0.3, iters=iters,
                            mode="surrogate", bits=3.3, seed=s,
                            energy_of=lambda sg, v=v:
                                problems.partition_difference(v, sg))
            diffs.append(r.best_energy)
            if r.best_energy <= (0 if v.sum() % 2 == 0 else 1):
                hits += 1
                iters_hit.append(int(np.argmin(r.energy)))
        n_sym = N * (N + 1)
        out.append({"N": N, "hit_rate": hits / trials,
                    "mean_best_difference": float(np.mean(diffs)),
                    "mean_iters": float(np.mean(iters_hit)) if iters_hit else None,
                    "symbols": n_sym,
                    "tts_s": time_to_solution(int(np.mean(iters_hit)), n_sym,
                                              106e9) if iters_hit else None})
        o = out[-1]
        print(f"  N={N:4d}  ground state {o['hit_rate']:.0%} of runs, "
              f"mean best difference {o['mean_best_difference']:.2f}"
              + (f", TTS {o['tts_s']*1e6:.1f} us" if o["tts_s"] else ""))
    save("fig5_partition", {
        "paper": "ground state for N = 16 to 256; TTS in Fig 5b",
        "hyperparameters": {"alpha": 1.0, "beta": 0.3, "bits": 3.3},
        "rows": out})

    def draw(fig):
        a1 = fig.add_subplot(1, 2, 1)
        a1.bar([str(o["N"]) for o in out], [100 * o["hit_rate"] for o in out],
               color="0.3")
        a1.set_xlabel("set size N")
        a1.set_ylabel("% of runs reaching the ground state")
        a1.set_title("Fig 5: number partitioning")
        a2 = fig.add_subplot(1, 2, 2)
        ok = [o for o in out if o["tts_s"]]
        if ok:
            a2.plot([o["N"] for o in ok], [o["tts_s"] * 1e6 for o in ok],
                    marker="o")
        a2.set_xlabel("set size N")
        a2.set_ylabel("time to solution (us)")
        a2.set_yscale("log")
        a2.set_title("Fig 5b: TTS, pipelined DSP")
    plot("fig5_partition", draw)


# ── Supplementary S3.2: precision against solution quality ──────────────────
def figs14_precision(trials=8, iters=400):
    """How many digital bits it takes to match the analogue machine.

    Supplementary S3.2 and Figure S14: a digital simulation needed 6-bit
    precision to match the hardware, which ran at 3.3 effective bits.  The
    paper attributes the difference to the machine's own noise acting as
    annealing.

    Reproduced here on a 20x20 lattice by sweeping the per-symbol precision
    with and without analogue dither.  The dithered 3.3-bit row against the
    undithered 6-bit row is the comparison.
    """
    J, h = problems.square_lattice(20, 20)
    gs = problems.lattice_ground_energy(20, 20)
    out = []
    for bits in (2.0, 2.5, 3.0, 3.3, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0):
        for dither in (0.0, 0.05):
            q = [100 * ising.solve(J, h, alpha=0.86, beta=1.0, iters=iters,
                                   mode="surrogate", bits=bits,
                                   sigma_rel=dither, seed=s).best_energy / gs
                 for s in range(trials)]
            out.append({"bits": bits, "dither": dither,
                        "quality_pct": float(np.mean(q)),
                        "hit_rate": sum(1 for x in q if x >= 99.99) / trials})
            print(f"  {bits:4.1f} bits, dither {dither:.2f}: "
                  f"{np.mean(q):6.2f} %  hit {out[-1]['hit_rate']:.0%}")
    save("figs14_precision", {
        "paper": "S3.2: digital needs 6 bits to match analogue at 3.3 bits",
        "rows": out})

    def draw(fig):
        ax = fig.add_subplot(1, 1, 1)
        for dither, label in ((0.0, "no dither (a digital solver)"),
                              (0.05, "with analogue dither")):
            r = [o for o in out if o["dither"] == dither]
            ax.plot([o["bits"] for o in r], [o["quality_pct"] for o in r],
                    marker="o", label=label)
        ax.set_xlabel("per-symbol precision (bits)")
        ax.set_ylabel("% of ground state")
        ax.set_title("Fig S14: precision against solution quality, 20x20 lattice")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    plot("figs14_precision", draw)


FIGURES = {
    "fig2_bifurcation": fig2_bifurcation,
    "fig2_mvm": fig2_mvm,
    "fig3_lattice": fig3_lattice,
    "fig5_partition": fig5_partition,
    "figs14_precision": figs14_precision,
}

if __name__ == "__main__":
    names = sys.argv[1:] or list(FIGURES)
    for nm in names:
        if nm not in FIGURES:
            raise SystemExit(f"unknown figure {nm!r}; have {list(FIGURES)}")
        print(f"\n=== {nm} ===")
        t0 = time.time()
        FIGURES[nm]()
        print(f"  [{time.time()-t0:.1f} s]")
