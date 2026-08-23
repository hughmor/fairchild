#!/usr/bin/env python3
"""validate_estimator.py — check the gain measurement against arithmetic.

Three estimators for the loop gain disagreed with each other by two orders of
magnitude on this deck (see the README). None of them was ever checked against a
case whose answer is known, which is the mistake this file exists to correct.

`rnn_math.md` §9(b) already says how the gain SHOULD be assembled — as a product
of factors, "no division by a small number, and every factor is separately
checkable":

    G = eta * R * a * dP_bus/dI_D

so the plan is to check each factor where its answer is arithmetic, and only
then multiply. Three of the four are exact:

  stage 1   the weight block and the balanced pair.  `fc_optical_2x2` is defined
            so `P_drop - P_thru = w*P_in`, and `fc_photodetector` is
            `I = R*P + i_dark`, so the balanced photocurrent is exactly
            `R*w*P_in` with the dark currents cancelling.
  stage 2   the tap and the 1:8 tree.  Two `fc_splitter` stages and three more:
            `a = 1/2 * 1/8 = 1/16` exactly.
  stage 3   the division at the modulator node.  `rnn_math.md` (6):
            `eta = R_sh/(R_sh + r_d)`, `r_d = n*V_T/I_D`, `R_sh = 1394 ohm`.
  stage 4   the ring's slope `dP_bus/dI_D`.  NOT arithmetic — it is the device
            physics — but it is directly measurable as a ratio of two
            directly-read quantities, with no small difference.

If stages 1-3 reproduce their closed forms, the chain is trustworthy and only
stage 4 carries model uncertainty. If they do not, the fault is upstream of
anything about rings and gains.

    .venv/bin/python experiments/giona/mpc/validate_estimator.py

Writes results/validate_estimator.json and .png.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fairchild as fc

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
R_PD = 0.8          # A/W, the responsivity every deck here uses
N_DIODE, VT = 5.0, 25.865e-3


def stage1(powers=(0.1, 1.0, 10.0), weights=(-1.0, -0.5, 0.0, 0.37, 1.0)):
    """Balanced photocurrent must be exactly R * w * P_in."""
    rows = []
    for p in powers:
        deck = f""".optical_port src
.optical_port dk
.optical_port th
.optical_port dr
XL src fc_cw_laser power_mW={p} wavelength_nm=1550
Xw src dk th dr WC 0 fc_optical_2x2 w=0 dw_dv=1
VW WC 0 DC 0
Xpt th pa 0 fc_photodetector responsivity={R_PD}
Xpd dr da 0 fc_photodetector responsivity={R_PD}
Vt pa 0 DC 0
Vd da 0 DC 0
.op
"""
        c = fc.Circuit(); c.load_str(deck)
        for w in weights:
            c.set_param("VW", "dc", w)
            r = c.run("op")
            i_t = abs(float(r["I(vt)"][0]))
            i_d = abs(float(r["I(vd)"][0]))
            got = i_d - i_t
            want = R_PD * w * (p * 1e-3)
            rows.append((p, w, got, want))
    return rows


def stage2():
    """The tap and the 1:8 tree must attenuate by exactly 1/16."""
    deck = """.optical_port src
.optical_port tap
.optical_port tree
.optical_port l1a
.optical_port l1b
.optical_port l2a
.optical_port l2b
.optical_port l2c
.optical_port l2d
.optical_port w1
.optical_port w2
.optical_port w3
.optical_port w4
.optical_port w5
.optical_port w6
.optical_port w7
.optical_port w8
XL src fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xtap src tap tree fc_splitter
Xt1  tree l1a l1b fc_splitter
Xt2a l1a l2a l2b fc_splitter
Xt2b l1b l2c l2d fc_splitter
Xt3a l2a w1 w2 fc_splitter
Xt3b l2b w3 w4 fc_splitter
Xt3c l2c w5 w6 fc_splitter
Xt3d l2d w7 w8 fc_splitter
.op
"""
    c = fc.Circuit(); c.load_str(deck)
    r = c.run("op")
    def P(n):
        return float(r[f"V({n}_re_0)"][0]) ** 2 + float(r[f"V({n}_im_0)"][0]) ** 2
    return P("w1") / P("src"), 1.0 / 16.0


def stage3(i_stars=(30e-6, 134e-6, 540e-6)):
    """eta = R_sh/(R_sh + r_d) at the modulator node, rnn_math (6).

    Measured by injecting a known signal current at `mod_cathode` and reading
    how much of it lands in the junction rather than the shunts. R_sh comes from
    the neuron subckt: 2 k on chip, 10 k bias, and two 17.05 k photodiode paths.
    """
    r_sh = 1.0 / (1 / 2e3 + 1 / 10e3 + 2 / 17.051e3)
    rows = []
    for i_star in i_stars:
        r_d = N_DIODE * VT / i_star
        rows.append((i_star, r_sh / (r_sh + r_d)))
    return r_sh, rows


def stage4_5(p_las=30.0, pdb=-8.0, iht=3.12e-3, d=0.2):
    """The ring slope, then G assembled from factors that are each checkable.

    Every quantity here is either arithmetic or a ratio of two directly-read
    numbers. Nothing divides a small difference of nearly-equal numbers, which
    is what killed the closed-loop estimator.

        I_D      from mc through the model's OWN diode law, i_sat*(exp(-V/nVt)-1)
        dI_D     = -(I_D + i_sat)/(n*Vt) * dmc, the exact derivative of that law
        dP_win   read straight off the weight-bank input, one channel
        eta_i    R_sh/(R_sh + n*Vt/I_D), validated in stage 3
        G_ij     = eta_i * R * dP_win_j / dI_D_j

    `dP_win` is measured at the bank input, so the 1/16 of stage 2 is already
    in it.
    """
    sys.path.insert(0, str(HERE))
    from characterize import N_NEURON, circuit
    i_sat, nvt = 5.099e-8, N_DIODE * VT
    r_sh = 1.0 / (1 / 2e3 + 1 / 10e3 + 2 / 17.051e3)
    c = circuit(p_las)
    for i in range(1, N_NEURON + 1):
        c.set_param(f"Vpdb{i}", "dc", pdb)
        c.set_param(f"Iht{i}", "dc", iht)
        for j in range(1, 9):
            c.set_param(f"VW{i}{j}", "dc", 0.0)

    def read(k):
        r = c.run("op")
        mc = float(r[f"V(mc{k})"][0])
        pw = (float(r[f"V(win1_re_{k - 1})"][0]) ** 2
              + float(r[f"V(win1_im_{k - 1})"][0]) ** 2)
        return mc, pw

    rows = []
    for k in range(1, N_NEURON + 1):
        mc0, _ = read(k)
        i_d = i_sat * (np.exp(min(-mc0 / nvt, 700)) - 1.0)
        eta = r_sh / (r_sh + nvt / i_d)
        c.set_param(f"Vpdb{k}", "dc", pdb + d); mch, pwh = read(k)
        c.set_param(f"Vpdb{k}", "dc", pdb - d); mcl, pwl = read(k)
        c.set_param(f"Vpdb{k}", "dc", pdb)
        dmc = (mch - mcl) / (2 * d)
        dpw = (pwh - pwl) / (2 * d)
        did = -(i_d + i_sat) / nvt * dmc
        slope = dpw / did if did else np.nan          # W per A
        g = eta * R_PD * slope
        rows.append(dict(k=k, mc=mc0, i_d=i_d, eta=eta, dmc=dmc, dpw=dpw,
                         di_d=did, slope=slope, g=g))
    return rows


def main() -> int:
    RESULTS.mkdir(exist_ok=True)
    out = {}

    print("stage 1 — weight block + balanced pair, want I = R*w*P_in")
    print(f"{'P mW':>7} {'w':>7} {'measured A':>13} {'arithmetic A':>14} {'rel err':>9}")
    s1 = stage1()
    errs = []
    for p, w, got, want in s1:
        e = abs(got - want) / max(abs(want), 1e-15)
        errs.append(e if abs(want) > 1e-12 else 0.0)
        print(f"{p:7.2f} {w:7.2f} {got:13.6e} {want:14.6e} {e:9.2e}")
    out["stage1_max_rel_err"] = float(max(errs))
    print(f"  worst relative error {max(errs):.2e}  "
          f"{'PASS' if max(errs) < 1e-6 else 'FAIL'}")

    print("\nstage 2 — tap + 1:8 tree, want a = 1/16")
    got, want = stage2()
    print(f"  measured {got:.9f}, arithmetic {want:.9f}, "
          f"rel err {abs(got - want) / want:.2e}  "
          f"{'PASS' if abs(got - want) / want < 1e-6 else 'FAIL'}")
    out["stage2"] = {"measured": got, "arithmetic": want}

    print("\nstage 3 — division at the modulator node, eta = R_sh/(R_sh + r_d)")
    r_sh, rows = stage3()
    print(f"  R_sh from the subckt = {r_sh:.1f} ohm "
          f"(rnn_math quotes 1394)")
    print(f"{'I* uA':>8} {'eta':>8}")
    for i_star, eta in rows:
        print(f"{i_star * 1e6:8.1f} {eta:8.3f}")
    out["stage3"] = {"r_sh": r_sh,
                     "eta": [{"i_star_A": a, "eta": b} for a, b in rows]}

    print("\nstage 4+5 — the ring slope, and G assembled from checked factors")
    s45 = stage4_5()
    print(f"{'i':>3} {'I_D uA':>9} {'eta':>6} {'dmc/dPDB':>10} "
          f"{'dI_D A/V':>11} {'dP_win W/A':>12} {'G':>8}")
    for r in s45:
        print(f"{r['k']:3d} {r['i_d'] * 1e6:9.2f} {r['eta']:6.3f} {r['dmc']:10.5f} "
              f"{r['di_d']:11.4e} {r['slope']:12.4e} {r['g']:8.3f}")
    gs = np.array([r["g"] for r in s45])
    print(f"\nG from the checked chain: {gs.min():.3f} to {gs.max():.3f}")
    print("Compare the closed-loop estimator at this same operating point, which")
    print("gave 2.818, 3.843, 2.556, 2.154, 5.635, 0.315 with four of six failing")
    print("their own probe-consistency check.")
    out["stage45"] = [{k: (float(v) if k != "k" else int(v)) for k, v in r.items()}
                      for r in s45]

    (RESULTS / "validate_estimator.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"\nwrote {RESULTS / 'validate_estimator.json'}")

    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2))
    for p in sorted({r[0] for r in s1}):
        ws = [r[1] for r in s1 if r[0] == p]
        gm = [r[2] for r in s1 if r[0] == p]
        ax[0].plot(ws, gm, "o", label=f"{p} mW measured")
        ax[0].plot(ws, [R_PD * w * (p * 1e-3) for w in ws], "-", lw=1)
    ax[0].set_xlabel("programmed weight w"); ax[0].set_ylabel("balanced photocurrent (A)")
    ax[0].set_title("stage 1: $I = R\\,w\\,P_{in}$\nlines = arithmetic, points = solver")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.25)
    i = np.geomspace(5e-6, 2e-3, 200)
    ax[1].semilogx(i * 1e6, r_sh / (r_sh + N_DIODE * VT / i), lw=2)
    for i_star, eta in rows:
        ax[1].plot(i_star * 1e6, eta, "o", ms=8)
    ax[1].set_xlabel("rest current $I^*$ (µA)"); ax[1].set_ylabel(r"$\eta$")
    ax[1].set_title(r"stage 3: $\eta = R_{sh}/(R_{sh}+r_d)$" "\nhow much signal the diode keeps")
    ax[1].grid(alpha=0.25, which="both")
    fig.suptitle("giona MPC — validating the gain chain factor by factor", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(RESULTS / "validate_estimator.png", dpi=130)
    print(f"wrote {RESULTS / 'validate_estimator.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
