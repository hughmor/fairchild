#!/usr/bin/env python3
"""characterize.py — find the operating point the MPC mapping needs.

The mapping is `Id − G·W = P/s` with the programmed weight |W| <= 1, so the loop
gain `G` is what decides whether the problem fits on the hardware at all. For
the design point in this directory's README that threshold is **G >= 1.18**.

Three knobs set it, and they are not independent:

  laser power   G is proportional to it through kappa = R*a*P0 — but the ring
                absorbs some of that light, warms, and detunes off the very
                resonance the modulation depends on. More power is not simply
                more gain.
  ring trim     `Iht<i>`. The transfer's slope is zero at the notch bottom and
                changes sign across it, so an untrimmed ring has a small G of
                the wrong sign (`rnn_math.md` §9).
  neuron bias   `PDB<i>`. Sets the rest current, which sets both the diode's
                share of the signal current and how much free-carrier
                absorption has already washed out the resonance.

`G` is measured closed-loop, the `rnn_math.md` (11) way: put a self-weight `w`
on one neuron and watch how its bias sensitivity changes,
`G = (1 - S0/S(w))/w`. One neuron at a time — a common-mode bias step moves
every ring on the shared bus and folds crosstalk into the answer.

    .venv/bin/python experiments/giona/mpc/characterize.py            # the operating point
    .venv/bin/python experiments/giona/mpc/characterize.py --scan     # the map it came from

Reports per-neuron `G`, which is what the driver needs: a neuron-to-neuron
spread is a row scaling on W and is calibratable, but only if it is measured.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

import fairchild as fc

HERE = Path(__file__).resolve().parent      # this experiment's own directory
GIONA = HERE.parent                         # experiments/giona
RESULTS = HERE / "results"
DECK = HERE / "netlists" / "giona_mpc_idealW.sp"
N_NEURON = 6

# Where the scan below put it. G = 2.49 against the 1.18 the weights need.
OPERATING_POINT = {"p_las_mW": 30.0, "pdb_V": -8.0, "iht_A": 3.12e-3}


def circuit(p_mw: float) -> fc.Circuit:
    """The deck at a given laser power.

    `p_las` is a `.param`, substituted at parse time, so power is a reload and
    not a `set_param`.
    """
    src = DECK.read_text().replace(".param p_las=30", f".param p_las={p_mw}")
    c = fc.Circuit()
    c.load_str(src)
    for i in range(1, 9):
        for j in range(1, 9):
            c.set_param(f"VW{i}{j}", "dc", 0.0)
    return c


def gain(c, i: int, pdb: float, iht: float, w: float = 0.05, d: float = 0.2) -> float:
    """Closed-loop small-signal gain of neuron `i`.

    `G = (1 - S0/S(w))/w` from `rnn_math.md` (11), where `S(w) = S0/(1 - w*G)`.
    The probe weight has to satisfy two things at once: big enough that S0 and
    S(w) differ resolvably, small enough that `w*G` stays well under 1 — because
    at `w*G -> 1` the self-coupled neuron is at its own bistability threshold and
    S(w) diverges and then changes sign.

    `rnn_math.md` recommends w in [0.1, 0.3] for a deck where G was around 1.
    Here G reaches 5, so w = 0.2 puts w*G at 1 and the estimator returns
    nonsense: measured S(w) ran 0.086, 1.15, then -0.68 across a heater sweep.
    0.05 keeps w*G under 0.3 at the largest G seen and still leaves S(w) a third
    away from S0. `gain_validated` is the check that this is in range.

    `d = 0.2 V`, and not the 0.02 that looks more like a small signal, because
    the SOLVER's resolution sets a floor. `reltol` is 1e-3 on a node sitting at
    977 mV, which is a millivolt of slop — and at d = 0.02 the whole response
    being differenced is 1.4 mV. That is how neuron 5 came back as -4255. At
    d = 0.2 the response is 14 mV and the same measurement gives 5.1 to 6.2.

    Tightening `reltol` instead does NOT work: at 1e-10 the solver starts
    chasing the gmin-floored dark optical nets — `ain1_re_7` asks to move 5e47 V
    — and fails to converge at all. The perturbation has to carry the burden.
    """
    c.set_param(f"Iht{i}", "dc", iht)

    def sens() -> float:
        c.set_param(f"Vpdb{i}", "dc", pdb + d)
        a = float(c.run("op")[f"V(mc{i})"][0])
        c.set_param(f"Vpdb{i}", "dc", pdb - d)
        b = float(c.run("op")[f"V(mc{i})"][0])
        c.set_param(f"Vpdb{i}", "dc", pdb)
        return (a - b) / (2 * d)

    c.set_param(f"VW{i}{i}", "dc", 0.0)
    s0 = sens()
    c.set_param(f"VW{i}{i}", "dc", w)
    sw = sens()
    c.set_param(f"VW{i}{i}", "dc", 0.0)
    if not sw or sw * s0 <= 0:
        # S(w) has crossed zero or flipped sign: w*G >= 1 and the neuron's own
        # self-coupling has gone unstable. The number is not a gain.
        return float("nan")
    return (1.0 - s0 / sw) / w


def gain_validated(c, i: int, pdb: float, iht: float,
                   ws=(0.03, 0.06), d: float = 0.2) -> tuple[float, float]:
    """(G, disagreement) across two probe weights.

    G is a small-signal quantity, so it must not depend on the probe. The spread
    between two `w` is the test that the estimator is in range — this is the
    check whose absence let a heater sweep report 1.64, 1.81, 4.86, 5.12 as if
    those were all gains.
    """
    gs = [gain(c, i, pdb, iht, w=w, d=d) for w in ws]
    if any(not np.isfinite(g) for g in gs):
        return float("nan"), float("inf")
    m = float(np.mean(gs))
    return m, float(abs(gs[1] - gs[0]) / abs(m)) if m else float("inf")


def scan() -> None:
    """G over power, bias and trim — the map the operating point came from."""
    ihts = np.linspace(0, 5e-3, 9)
    print(f"{'P/ch':>6} {'PDB':>6} | " + " ".join(f"{v * 1e3:6.2f}" for v in ihts)
          + "   <- Iht mA")
    for p in (10.0, 30.0):
        c = circuit(p)
        for pdb in (-8.0, -14.0, -20.0):
            row = [gain(c, 1, pdb, float(v)) for v in ihts]
            print(f"{p:6.1f} {pdb:6.1f} | " + " ".join(f"{g:6.2f}" for g in row))
    print("\nMore power is not simply more gain: at 30 mW/channel the rings have")
    print("self-heated 236 pm off resonance, which is 1.4 linewidths, so G grows")
    print("sublinearly — 3x the power from 10 to 30 mW buys 2.8x the gain.")


def per_neuron() -> None:
    op = OPERATING_POINT
    c = circuit(op["p_las_mW"])
    for i in range(1, N_NEURON + 1):
        c.set_param(f"Vpdb{i}", "dc", op["pdb_V"])
        c.set_param(f"Iht{i}", "dc", op["iht_A"])
    print(f"operating point: {op['p_las_mW']} mW/channel, PDB = {op['pdb_V']} V, "
          f"Iht = {op['iht_A'] * 1e3:.2f} mA")
    print(f"\n{'neuron':>7} {'G':>7} {'probe spread':>13} {'mc mV':>10}")
    gs, dis, mcs = [], [], []
    for i in range(1, N_NEURON + 1):
        g, sp = gain_validated(c, i, op["pdb_V"], op["iht_A"])
        mc = float(c.run("op")[f"V(mc{i})"][0]) * 1e3
        gs.append(g); dis.append(sp); mcs.append(mc)
        flag = "" if sp < 0.15 else "  <- probe out of range"
        print(f"{i:7d} {g:7.3f} {sp:12.1%} {mc:10.3f}{flag}")
    gs = np.array(gs)
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "per_neuron_gain.json").write_text(json.dumps({
        **op, "G": [float(v) for v in gs],
        "probe_spread": [float(v) for v in dis],
        "mc_mV": [float(v) for v in mcs],
    }, indent=2) + "\n")
    print(f"\nwrote {RESULTS / 'per_neuron_gain.json'}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    n = np.arange(1, N_NEURON + 1)
    ax.bar(n, gs, color=["#4c72b0" if v < 0.15 else "#c44e52" for v in dis])
    ax.axhline(1.18, color="k", ls="--", lw=1.2)
    ax.annotate("MPC needs 1.18", (0.6, 1.18), fontsize=9,
                textcoords="offset points", xytext=(2, 4))
    for x, (g, sp) in enumerate(zip(gs, dis), start=1):
        ax.errorbar(x, g, yerr=abs(g) * sp, color="k", capsize=4, lw=1)
    ax.set_xlabel("neuron"); ax.set_ylabel("loop gain G")
    ax.set_title(f"per-neuron loop gain at one shared trim\n"
                 f"{op['p_las_mW']:.0f} mW/ch, PDB {op['pdb_V']} V, "
                 f"Iht {op['iht_A'] * 1e3:.2f} mA "
                 f"(bars = probe-weight spread; red = out of range)")
    ax.grid(alpha=0.25, axis="y")
    fig.tight_layout()
    fig.savefig(RESULTS / "per_neuron_gain.png", dpi=130)
    print(f"wrote {RESULTS / 'per_neuron_gain.png'}")
    print(f"G spans {gs.min():.2f}-{gs.max():.2f}, a factor of {gs.max() / gs.min():.1f}")
    print(f"MPC needs G >= 1.18 (see README): "
          f"{'every neuron clears it' if gs.min() >= 1.18 else 'SHORT on some neuron'}, "
          f"margin {gs.min() / 1.18:.2f}x on the worst")
    print("\nThe spread is the thing to fix, not the margin. Row i of W scales by")
    print("1/G_i, so the strong neurons use a small fraction of the weight range and")
    print("their weights land proportionally coarser: at this spread the best-off row")
    print(f"uses {1.18 / gs.max():.0%} of full scale where the worst-off uses "
          f"{1.18 / gs.min():.0%}, so its precision in units of the answer is "
          f"{gs.max() / gs.min():.1f}x worse.")
    print("A single shared trim did this — every ring got ring 1's Iht. Trimming each")
    print("ring to equalise G is the next step, and it is what makes the weight")
    print("precision budget in the README reachable on every row rather than one.")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", action="store_true",
                    help="the power/bias/trim map, not just the chosen point")
    args = ap.parse_args()
    if not DECK.exists():
        print(f"no deck at {DECK} — run build_mpc_deck.py first", file=sys.stderr)
        return 1
    scan() if args.scan else per_neuron()
    return 0


if __name__ == "__main__":
    sys.exit(main())
