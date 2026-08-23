#!/usr/bin/env python3
"""compare_va_mrm.py — the discrete giona ring vs the Verilog-A add-drop ring.

The giona neurons are built from `examples/photonic/pcells/mrm.sp`: two
`fc_dcoupler`s and two `fc_pn_th_ps LEVEL=4` arcs. The replacement is one
Verilog-A module, `examples/verilog_a/models/mrm_addrop.va`, whose defaults come
from the same fitted card.

Both rings hang off ONE laser in ONE deck and are read on the same solve, so a
difference in the answer is the two models disagreeing and nothing else.

**What is being matched, and what deliberately is not.** The new model is a
superset: it solves for a carrier density and a temperature where the cell
carries fitted coefficients (`dn_di`, `da_di`) and a fitted watts-per-pi
(`p_pi_th`). So the target is agreement *up to the limitations of the cell* —
every static mechanism the cell has, matched at steady state; every mechanism it
does not have, left at a plausible physical value rather than tuned to make the
new model as blind as the old one. Concretely:

  matched (analytic, in the model header)   left alone (new physics)
  ---------------------------------------   ---------------------------------
  coupling, dispersion, background loss     optical self-heating (r_th != 0)
  depletion EO, junction I(V), junction C   junction self-heating
  steady-state dn/dI and dalpha/dI          carrier lifetime as a real pole
  heater watts-per-pi at lambda_ref         thermal time constant c_th

The three sections below follow from that. Section 1 checks the matched half
with the thermal path muted, where the two models are the same equations written
twice and must agree to a picometre. Section 2 runs the shipping defaults, where
they must NOT agree, and the gap should be the size the new physics predicts.
Section 3 is that new physics on its own axis.

    .venv/bin/python experiments/giona/compare_va_mrm.py
    .venv/bin/python experiments/giona/compare_va_mrm.py --refit

`--refit` re-derives `vol_active`/`dalpha_dnc` numerically in the muted-thermal
configuration. It is a check on the algebra, not a fit anyone needs: if the
analytic map is right it comes back at 1.00x.

Writes results/va_mrm_compare.png and results/va_mrm_match.json.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault(
    "FAIRCHILD_OPENVAF",
    "/Users/hugh/Local/src/OpenVAF-Reloaded/target/release/openvaf-r",
)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fairchild as fc

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
RESULTS = HERE / "results"
VA = REPO / "examples" / "verilog_a" / "models" / "mrm_addrop.va"
CELL = REPO / "examples" / "photonic" / "pcells" / "mrm.sp"

# Ring input power. `ringfit.py` puts a 1 mW laser through a 9 dB grating
# coupler, so the captures the card was fitted to saw about this much.
P_MW = 0.126
# `r_th` small enough that no dissipation raises a temperature worth having.
# Only used to isolate the matched half; it is not a value anyone should ship.
R_TH_MUTED = 1e-6

DECK = f""".include {CELL}
.va {VA}
.optical_port src
.optical_port d_th
.optical_port d_ad
.optical_port d_dr
.optical_port v_th
.optical_port v_ad
.optical_port v_dr
XL src fc_cw_laser power_mW={P_MW} wavelength_nm=1550
Xd src d_th d_ad d_dr pn 0 htr 0 mrm
Xv src v_th v_ad v_dr pn 0 htr 0 tv mrm_addrop
VPN pn 0 DC 0
VHT htr 0 DC 0
.op
"""

REVERSE = [(v, 0.0) for v in (0.0, -1.0, -2.0, -3.0, -4.0)]
FORWARD = [(v, 0.0) for v in (0.70, 0.80, 0.85, 0.90)]
HEATER = [(0.0, v) for v in (0.5, 1.0, 1.5, 2.0)]
BIAS = REVERSE + FORWARD + HEATER

R_TH_SHIP = 3139.9          # the analytic map of p_pi_th = 26.4 mW/pi
COLD_NM = 1549.674          # unbiased resonance, the walk starts here
R_HEATER = 368.8            # whole-ring heater resistance (2 x the cell's per-arc)
P_PI_TH = 26.4e-3           # watts per pi, straight off the card
FSR_NM = 11.38              # lambda^2 / (n_g * L_ring)


def power(r, net: str) -> float:
    re = float(r[f"V({net}_re_0)"][0])
    im = float(r[f"V({net}_im_0)"][0])
    return re * re + im * im


def t_db(ckt, wl_nm: float) -> tuple[float, float]:
    """Thru-port transmission of (cell, VA) at one wavelength, in dB."""
    ckt.set_param("XL", "wavelength_nm", float(wl_nm))
    r = ckt.run("op")
    p_in = max(power(r, "src"), 1e-30)
    return (10.0 * np.log10(max(power(r, "d_th"), 1e-30) / p_in),
            10.0 * np.log10(max(power(r, "v_th"), 1e-30) / p_in))


# (half-span nm, samples). Each round's span comfortably exceeds the previous
# round's half-spacing, so the true minimum cannot fall outside the next window;
# round 1's 50 pm spacing cannot step over a 167 pm notch.
ROUNDS = ((0.60, 25), (0.060, 13), (0.012, 13))


class NotchLost(RuntimeError):
    """The bracket did not contain a resonance."""


def locate(ckt, which: int, guess_nm: float) -> tuple[float, float]:
    """(resonance nm, depth dB) by successive bracketing.

    A full spectrum is 500 solves at 33 ms each and all but three get thrown
    away; three narrowing rounds cost 51 and land inside a picometre.

    It raises rather than returning its best guess. An argmin at the edge of the
    window means the window held no notch, and what comes back then is the
    lowest point of a flat baseline — a number that looks like a resonance,
    plots like a resonance, and is off by nanometres. Both models missing the
    same way even agree with each other, which is the worst version of it: this
    script read a 2 V heater point as -0.27 dB and called the two models
    identical there.
    """
    centre = guess_nm
    for span, n in ROUNDS:
        wl = np.linspace(centre - span, centre + span, n)
        y = np.array([t_db(ckt, w)[which] for w in wl])
        i = int(np.argmin(y))
        if i == 0 or i == n - 1:
            raise NotchLost(
                f"no notch within {span:.3f} nm of {centre:.4f} on port {which}: "
                f"the minimum is at the window edge ({y[i]:.2f} dB, baseline "
                f"{np.median(y):.2f} dB)")
        centre = float(wl[i])
    if np.median(y) - y[i] < 0.02:
        raise NotchLost(f"the final bracket at {centre:.4f} is flat "
                        f"({y[i]:.3f} dB vs {np.median(y):.3f} dB baseline)")
    y0, y1, y2 = y[i - 1], y[i], y[i + 1]
    d = 2 * y1 - y0 - y2
    step = wl[1] - wl[0]
    if abs(d) < 1e-12:
        return centre, float(y1)
    frac = float(np.clip((y0 - y2) / (2 * d), -1, 1))
    return centre + frac * step, float(y1 - 0.25 * (y0 - y2) * frac)


def spectrum(ckt, centre_nm: float, half_nm: float = 0.45, n: int = 121):
    wl = np.linspace(centre_nm - half_nm, centre_nm + half_nm, n)
    both = np.array([t_db(ckt, w) for w in wl])
    return wl, both[:, 0], both[:, 1]


def at_bias(ckt, v_pn: float, v_ht: float, r_th: float) -> dict:
    ckt.set_param("VPN", "dc", v_pn)
    ckt.set_param("VHT", "dc", v_ht)
    ckt.set_param("Xv", "r_th", r_th)
    # Where the heater puts the ring, from the card itself: P/p_pi_th pi of
    # phase is that fraction of half an FSR. Guessing this by eye put the 2 V
    # window 300 pm off the notch and the locator returned baseline.
    guess = COLD_NM + (v_ht * v_ht / R_HEATER) / P_PI_TH * FSR_NM / 2.0
    d_res, d_dep = locate(ckt, 0, guess)
    v_res, v_dep = locate(ckt, 1, guess)
    return {"v_pn": v_pn, "v_ht": v_ht, "d_res": d_res, "d_dep": d_dep,
            "v_res": v_res, "v_dep": v_dep}


def grid(ckt, bias, r_th: float) -> list[dict]:
    zero = at_bias(ckt, 0.0, 0.0, r_th)
    rows = [at_bias(ckt, p, h, r_th) for p, h in bias]
    for r in rows:
        r["d_shift"] = (r["d_res"] - zero["d_res"]) * 1e3
        r["v_shift"] = (r["v_res"] - zero["v_res"]) * 1e3
        r["err"] = r["v_shift"] - r["d_shift"]
        r["derr"] = r["v_dep"] - r["d_dep"]
    return rows


def report(title: str, rows) -> None:
    print(f"\n── {title} ──")
    print(f"{'V_pn':>6} {'V_ht':>5} | {'cell shift':>12} {'VA shift':>11} {'err':>9}"
          f" | {'cell depth':>11} {'err':>8}")
    for r in rows:
        print(f"{r['v_pn']:6.2f} {r['v_ht']:5.2f} | {r['d_shift']:9.1f} pm"
              f" {r['v_shift']:8.1f} pm {r['err']:6.1f} pm"
              f" | {r['d_dep']:8.2f} dB {r['derr']:5.2f} dB")
    print(f"  worst shift error {max(abs(r['err']) for r in rows):.1f} pm,"
          f"  worst depth error {max(abs(r['derr']) for r in rows):.2f} dB")


def refit(ckt) -> dict[str, float]:
    """Re-derive vol_active / dalpha_dnc from the cell, thermal path muted.

    Two coefficients against the forward-bias shift and depth, in the one
    configuration where the two models are supposed to be identical. This exists
    to catch an algebra slip in the analytic map, so it starts from that map: a
    result that is not ~1.00x means the header is wrong.
    """
    from scipy.optimize import least_squares
    start = {"vol_active": 2.7531e-17, "dalpha_dnc": 1.0212e-21}
    ref = grid(ckt, FORWARD, R_TH_MUTED)

    def res(x):
        ckt.set_param("Xv", "vol_active", float(x[0] * start["vol_active"]))
        ckt.set_param("Xv", "dalpha_dnc", float(x[1] * start["dalpha_dnc"]))
        rows = grid(ckt, FORWARD, R_TH_MUTED)
        return np.array([v for r, r0 in zip(rows, ref)
                         for v in (r["v_shift"] - r0["d_shift"], r["derr"] * 10.0)])

    sol = least_squares(res, [1.0, 1.0], bounds=([0.1, 0.1], [10.0, 10.0]),
                        diff_step=5e-3, xtol=1e-4, max_nfev=25, verbose=2)
    out = {k: float(v * start[k]) for k, v in zip(start, sol.x)}
    print("\n  analytic map vs numerical re-derivation:")
    for (k, v0), (_, v) in zip(start.items(), out.items()):
        print(f"    {k:12s} {v0:12.5g} -> {v:12.5g}   ({v / v0:.4f}x)")
    for k, v in start.items():          # leave the model on its shipped defaults
        ckt.set_param("Xv", k, v)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refit", action="store_true")
    args = ap.parse_args()
    RESULTS.mkdir(exist_ok=True)

    ckt = fc.Circuit()
    ckt.load_str(DECK)
    out: dict = {"p_mW": P_MW}

    # ── 1. the matched half ─────────────────────────────────────────────────
    print("1. matched half — thermal path muted, so only what the cell also has")
    # Heater rows are absent on purpose: r_th is the heater's ONLY path to the
    # index here, so muting it mutes the heater. They are checked in section 2,
    # which is where r_th carries p_pi_th rather than adding anything.
    muted = grid(ckt, REVERSE + FORWARD, R_TH_MUTED)
    report(f"r_th = {R_TH_MUTED} (muted), {P_MW} mW into the ring — bias only", muted)
    out["muted_worst_pm"] = max(abs(r["err"]) for r in muted)

    if args.refit:
        out["refit"] = refit(ckt)

    # ── 2. shipping defaults ────────────────────────────────────────────────
    print("\n2. shipping defaults — the gap IS the new physics, not an error")
    ship = grid(ckt, BIAS, R_TH_SHIP)
    report(f"r_th = {R_TH_SHIP} K/W (== p_pi_th 26.4 mW/pi), {P_MW} mW", ship)
    out["ship_worst_pm"] = max(abs(r["err"]) for r in ship)

    # The heater rows are the one place r_th is a MATCH, not an addition: they
    # are how p_pi_th was transcribed, so they have to come back right.
    htr = [r for r in ship if r["v_ht"] > 0]
    print(f"\n  heater rows (these are matched, not new): worst"
          f" {max(abs(r['err']) for r in htr):.1f} pm"
          f"  -> p_pi_th reproduced to {max(abs(r['err'] / r['d_shift']) for r in htr):.2%}")
    out["heater_worst_pm"] = max(abs(r["err"]) for r in htr)

    # ── 3. the new physics on its own axis ──────────────────────────────────
    print("\n3. optical self-heating: resonance vs input power")
    powers = np.geomspace(1e-3, 4.0, 11)
    walk_d, walk_v = [], []
    ckt.set_param("VPN", "dc", 0.0)
    ckt.set_param("VHT", "dc", 0.0)
    for p in powers:
        ckt.set_param("XL", "power_mW", float(p))
        ckt.set_param("Xv", "r_th", R_TH_SHIP)
        walk_d.append(locate(ckt, 0, COLD_NM)[0])
        walk_v.append(locate(ckt, 1, COLD_NM + 0.06 * p)[0])
        print(f"  {p:7.3f} mW   cell {walk_d[-1]:10.5f}   VA {walk_v[-1]:10.5f}"
              f"   {(walk_v[-1] - walk_d[-1]) * 1e3:+8.1f} pm")
    ckt.set_param("XL", "power_mW", P_MW)
    out["self_heating_pm_at_1mW"] = float(
        np.interp(1.0, powers, (np.array(walk_v) - np.array(walk_d)) * 1e3))

    (RESULTS / "va_mrm_match.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"\nwrote {RESULTS / 'va_mrm_match.json'}")
    plot(ckt, muted, ship, powers, walk_d, walk_v)


def plot(ckt, muted, ship, powers, walk_d, walk_v) -> None:
    fig = plt.figure(figsize=(15.5, 9.0))
    gs = fig.add_gridspec(2, 3, hspace=0.34, wspace=0.27)

    # (0,0)/(0,1) spectra, reverse and forward, at the shipping defaults.
    for col, (title, bias) in enumerate([("reverse bias", REVERSE),
                                         ("forward bias", [(0.0, 0.0)] + FORWARD)]):
        ax = fig.add_subplot(gs[0, col])
        cmap = plt.get_cmap("viridis")
        for i, (v_pn, v_ht) in enumerate(bias):
            ckt.set_param("VPN", "dc", v_pn)
            ckt.set_param("VHT", "dc", v_ht)
            ckt.set_param("Xv", "r_th", R_TH_SHIP)
            wl, d, v = spectrum(ckt, COLD_NM, 0.5, 81)
            c = cmap(i / max(len(bias) - 1, 1))
            ax.plot(wl, d, "-", color=c, lw=2.4, alpha=0.5)
            ax.plot(wl, v, "--", color=c, lw=1.1)
            ax.plot([], [], "-", color=c, label=f"{v_pn:+.2f} V")
        ax.set_title(f"thru port, {title} @ {P_MW} mW\nsolid = cell, dashed = Verilog-A",
                     fontsize=10)
        ax.set_xlabel("wavelength (nm)")
        ax.set_ylabel("T (dB)")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.25)

    # (0,2) heater tuning — a matched mechanism, so the curves must lie on top.
    ax = fig.add_subplot(gs[0, 2])
    hs = [r for r in ship if r["v_ht"] > 0]
    p_htr = np.array([r["v_ht"] ** 2 / 368.8 * 1e3 for r in hs])
    ax.plot(p_htr, [r["d_shift"] for r in hs], "o-", lw=2.4, alpha=0.6, label="cell")
    ax.plot(p_htr, [r["v_shift"] for r in hs], "s--", label="Verilog-A")
    ax.set_xlabel("heater power (mW)")
    ax.set_ylabel("resonance shift (pm)")
    ax.set_title("heater tuning — matched\n(r_th transcribes p_pi_th = 26.4 mW/pi)",
                 fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)

    # (1,0) matched half.
    ax = fig.add_subplot(gs[1, 0])
    lab_m = [f"{r['v_pn']:+.2f}V" for r in muted]
    ax.bar(range(len(muted)), [r["err"] for r in muted],
           color=["#4c72b0" if r["v_pn"] <= 0 else "#c44e52" for r in muted])
    ax.set_xticks(range(len(muted)))
    ax.set_xticklabels(lab_m, fontsize=7, rotation=60)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_ylim(-1, 1)
    ax.set_ylabel("shift error, VA − cell (pm)")
    ax.set_title("1. matched half, thermal muted\nreverse (blue) | forward (red)",
                 fontsize=10)
    ax.grid(alpha=0.25, axis="y")
    ax.text(0.5, 0.5, f"max |error| = {max(abs(r['err']) for r in muted):.2f} pm\n"
                      f"over {len(muted)} bias points\n"
                      f"depth to {max(abs(r['derr']) for r in muted):.3f} dB",
            transform=ax.transAxes, ha="center", va="center", fontsize=11,
            bbox={"boxstyle": "round", "fc": "#eef4ee", "ec": "#7aa87a"})

    # (1,1) shipping defaults — the same bars, with the new physics switched on.
    ax = fig.add_subplot(gs[1, 1])
    lab_s = [(f"{r['v_pn']:+.2f}V" if r["v_ht"] == 0 else f"htr {r['v_ht']:.1f}V")
             for r in ship]
    ax.bar(range(len(ship)), [r["err"] for r in ship],
           color=["#4c72b0" if r["v_pn"] < 0 else
                  ("#dd8452" if r["v_ht"] > 0 else "#c44e52") for r in ship])
    ax.set_xticks(range(len(ship)))
    ax.set_xticklabels(lab_s, fontsize=7, rotation=60)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_ylabel("shift error, VA − cell (pm)")
    ax.set_title("2. the gap the new physics opens\nblue/red = self-heating, orange = heater (matched)",
                 fontsize=10)
    ax.grid(alpha=0.25, axis="y")

    # (1,2) self-heating walk.
    ax = fig.add_subplot(gs[1, 2])
    ax.semilogx(powers, (np.array(walk_d) - walk_d[0]) * 1e3, "o-", lw=2.4, alpha=0.6,
                label="cell (r_th = 0, blind)")
    ax.semilogx(powers, (np.array(walk_v) - walk_v[0]) * 1e3, "s-",
                label="Verilog-A (solved T)")
    ax.axvline(P_MW, color="k", ls=":", lw=1)
    ax.annotate("giona", (P_MW, 0), textcoords="offset points", xytext=(4, 6), fontsize=8)
    ax.set_xlabel("ring input power (mW)")
    ax.set_ylabel("resonance walk (pm)")
    ax.set_title("3. new physics: absorbed light\nheats the ring", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25, which="both")

    fig.suptitle("mrm_addrop.va vs the discrete mrm.sp cell — matched up to the "
                 "cell's limitations, not beyond them", fontsize=12)
    out = RESULTS / "va_mrm_compare.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
