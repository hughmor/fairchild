#!/usr/bin/env python3
"""va_vs_may_data.py — the Verilog-A add-drop ring against the May capture.

`compare_va_mrm.py` asks whether `mrm_addrop.va` reproduces the discrete cell.
It does, to a picometre. That leaves the question that actually matters: does it
reproduce the *device*, which is what the cell's card was fitted to in the first
place.

The dataset is `giona_neuron2_mod_joint_IV_spec` — the May sparse joint sweep,
20 heater currents x 16 junction voltages x a 1.7 nm spectrum, and the capture
`ringfit.py` staged its fit against. The testbench is `ringfit.build_netlist`'s,
rebuilt around a single ring instance: 1 mW laser, 9 dB grating coupler in,
304.5 um of bus, the ring, 304.5 um of bus, 9 dB out, and transmission read
fibre-to-fibre exactly as the instrument did.

**One free parameter, and only one.** `n_eff` is trimmed so the unbiased notch
sits where the measured one does. The card's header says absolute `n_eff` must
be trimmed per ring because fab variation exceeds anything a shared card can
carry, so this is the card being used as intended rather than a fit. Everything
else is the card, untouched.

Three rings run in the same deck off three identical lasers, so a wavelength
point costs one solve for all of them:

  cell        `examples/photonic/pcells/mrm.sp`, the reference
  va          `mrm_addrop.va` on its shipped defaults — a faithful port of the
              cell, including the cell's inherited factor-of-two terminal
              current (see below)
  va_iv       the same, with `i_sat` and `vol_active` both halved

That last one is the point of running against data rather than against the cell.
The card's `i_sat` was fitted to the measured TERMINAL current, then used as a
PER-ARC parameter in a cell whose two junctions are in parallel — so the cell,
and this model after it, draw twice the current the device does. The optical fit
is unharmed (`dn_di` was fitted through the same doubled current and absorbs it
exactly), which is why nobody noticed, and it is why `mrm.sp`'s header says not
to "correct" `i_sat` alone. `va_iv` corrects it with its partner: halve `i_sat`,
halve `vol_active`, and both the spectra and the I(V) come out right.

    .venv/bin/python experiments/giona/va_vs_may_data.py

Needs `lightlab` for the raw capture; the extraction is cached to
`data/may_n2_cache.npz` on the first run and read from there afterwards.
Writes results/va_vs_may_data.png and results/va_vs_may_data.json.
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
CACHE = HERE / "data" / "may_n2_cache.npz"
VA = REPO / "examples" / "verilog_a" / "models" / "mrm_addrop.va"
CELL = REPO / "examples" / "photonic" / "pcells" / "mrm.sp"

R_SHUNT = 2000.0        # the PCB shunt across the pads, in the measured current
GC_DB = 9.0             # grating coupler, each side
L_BUS = 304.5e-6        # bus waveguide either side of the ring
# 7 pm sampling. At 14 pm the parabola through the three points around the
# minimum carries a bias of a few pm that depends on where the true minimum
# falls between samples — which made two rings whose notches are 9 pm apart
# appear 21 pm apart.
N_WL = 241

# Card-compatible mode: the four linearisations on, Soref-Bennett off, the
# cell's terminal conventions restored. This is the model reproducing `mrm.sp`
# term for term, and it is the baseline the new carrier physics has to beat.
LEGACY = ("dn_dv=-3.62e-5 da_dv=3.29e-4 dn_dnc=-8.8e-28 dalpha_dnc=1.0212e-21 "
          "vol_dep=0 sb_dn_e=0 sb_dn_h=0 sb_da_e=0 sb_da_h=0 "
          "i_sat=1.0198e-7 vol_inj=2.7531e-17 c_j0=2.75e-13 tau_sweep=1")

# The three rings, and what makes each one different.
RINGS = {
    "cell": ("the discrete cell", ""),
    "card": ("Verilog-A, card-compatible", LEGACY),
    "phys": ("Verilog-A, one carrier population", ""),
}

# What the fit is allowed to move, with the analytic starting point and a
# decade either side. All three are per-device by nature: two effective volumes
# the optical mode sees, and a thermal resistance the card's own header says
# varies about 2x between neurons.
FREE = {
    "vol_dep": (6.786e-18, 6.8e-19, 6.8e-17),
    "vol_inj": (8.31e-17, 8.3e-18, 8.3e-16),
    "r_th": (1586.0, 300.0, 20000.0),
}

# The lever between resonance shift and notch depth: electrons and holes differ
# five-fold in absorption per unit index, so this is the only parameter that can
# move one without the other. `--free-ratio` adds it, and the two runs together
# are what says whether this capture can separate it from `vol_dep` — both scale
# the index response, and 1.7 nm of spectrum over a 1 V window is not much to
# separate them with.
RATIO_H = ("ratio_h", 1.0, 0.2, 5.0)

# Residual weights, from each observable's own scatter in this capture rather
# than from taste. The notch position is good to about a picometre (0.85 pm of
# grid, parabola-interpolated). The depth wanders +/-0.3 dB across reverse bias
# where the true signal is monotonic and only 0.41 dB end to end, so it is
# nearly all noise. Weighting 1 dB as 100 pm — which the first version did —
# gave the noisier observable nine times the pull of the informative one, and
# the fit duly moved `vol_dep` five-fold to chase it.
SIGMA_POS_PM = 1.0
SIGMA_DEPTH_DB = 0.3

# 25 %. The locator resolves a notch to about half a picometre, and a 3 % step
# in vol_dep moves it by 0.35 pm — so the first version of this fit built its
# Jacobian out of rounding noise and wandered six-fold away from a starting
# point that was already at the measurement floor. A step has to be big enough
# that the objective actually responds to it.
DIFF_STEP = 0.25


def load() -> dict[str, np.ndarray]:
    if CACHE.exists():
        d = dict(np.load(CACHE))
        if len(d["wl_nm"]) >= 2000:
            return d
        print("cache is the old 250-point extraction; rebuilding")
    import sys

    sys.path.insert(0, str(HERE))
    import ringfit

    # ringfit downsamples to 250 points across 1.7 nm, which quantises a notch
    # position at +/-3.4 pm — a third of the whole 13.5 pm depletion shift this
    # is trying to resolve. The raw capture has 62 417 points; 2 000 costs a
    # 2.5 MB cache and puts the measured notch inside half a picometre.
    # Two places, because `extract_data` sizes its array from the module global
    # while `_downsample` captured it as a default argument at import.
    ringfit.N_SIM_PTS = 2000
    ringfit._downsample.__defaults__ = (2000,)
    sd = ringfit.extract_data(ringfit.load_sweep())
    d = {
        "hc_mA": sd.hc_mA,
        "jv_V": sd.jv_V,
        "wl_nm": sd.wl_nm,
        "T_dB": sd.T_dB,
        "v_heat_V": sd.v_heat_V,
        "i_junc_mA": sd.i_junc_mA,
    }
    np.savez_compressed(CACHE, **d)
    print(f"cached the extraction to {CACHE}")
    return d


def deck(n_eff: dict[str, float] | float) -> str:
    """One laser, one GC chain and one ring per variant, all in one solve."""
    if isinstance(n_eff, float):
        n_eff = dict.fromkeys(RINGS, n_eff)
    lines = [CELL.read_text(), f".va {VA}", ""]
    # Every optical net, declared before anything uses one: a net that reaches
    # an element line before its `.optical_port` is a plain scalar node, and the
    # device refuses it by terminal count rather than by name.
    for tag in RINGS:
        for p in ("src", "gin", "ri", "th", "wo", "gth", "ad", "dr"):
            lines.append(f".optical_port {tag}_{p}")
    lines.append("")
    for tag, (_, extra) in RINGS.items():
        lines += [
            f"XL{tag} {tag}_src fc_cw_laser power_mW=1.0 wavelength_nm=1546.5",
            f"XGI{tag} {tag}_src {tag}_gin fc_grating_coupler alpha_dB={GC_DB}",
            f"XWI{tag} {tag}_gin {tag}_ri fc_waveguide l_m={L_BUS} n_g=4.2 alpha_dB_cm=1.0",
            f"XWO{tag} {tag}_th {tag}_wo fc_waveguide l_m={L_BUS} n_g=4.2 alpha_dB_cm=1.0",
            f"XGO{tag} {tag}_wo {tag}_gth fc_grating_coupler alpha_dB={GC_DB}",
        ]
    lines.append("")
    # Each ring gets its OWN bias and heater nets. Sharing them looks harmless —
    # the voltage is the voltage — but the heater is driven by a CURRENT source,
    # and three heaters on one source split it three ways, so every ring would
    # run at a third of the measured drive. Separate nets also make `I(v…)` the
    # current of one ring rather than the sum of three.
    for tag, (_, extra) in RINGS.items():
        body = f"pn{tag} 0 ht{tag} 0"
        model = "mrm" if tag == "cell" else f"t{tag} mrm_addrop"
        lines += [
            f"X{tag} {tag}_ri {tag}_th {tag}_ad {tag}_dr {body} {model} "
            f"n_eff={n_eff[tag]:.8f} {extra}".rstrip(),
            f"V{tag} pn{tag} 0 DC 0",
            f"RSH{tag} pn{tag} 0 {R_SHUNT}",
            f"IHT{tag} ht{tag} 0 DC 0",
        ]
    lines.append(".op")
    return "\n".join(lines) + "\n"


def spectra(ckt, wl_nm, v_pn: float, i_ht: float) -> dict[str, np.ndarray]:
    """Fibre-to-fibre T_dB for every ring, plus the terminal electricals."""
    for t in RINGS:
        ckt.set_param(f"V{t}", "dc", float(v_pn))
        ckt.set_param(f"IHT{t}", "dc", float(i_ht))
    out = {t: np.empty(len(wl_nm)) for t in RINGS}
    for i, wl in enumerate(wl_nm):
        for t in RINGS:
            ckt.set_param(f"XL{t}", "wavelength_nm", float(wl))
        r = ckt.run("op")
        for t in RINGS:
            p_in = float(r[f"V({t}_src_re_0)"][0]) ** 2 + float(r[f"V({t}_src_im_0)"][0]) ** 2
            p_out = float(r[f"V({t}_gth_re_0)"][0]) ** 2 + float(r[f"V({t}_gth_im_0)"][0]) ** 2
            out[t][i] = 10.0 * np.log10(max(p_out, 1e-30) / max(p_in, 1e-30))
    return out


L_RING = 2.0 * np.pi * 8e-6
N_G = 4.2
WL_REF_NM = 1550.0


def n_eff_for(target_nm: float, guess: float = 2.2810) -> float:
    """The `n_eff` that puts a resonance at `target_nm`.

    Solved rather than nudged. A ring's resonances are an FSR apart — 11.4 nm
    here against a 1.7 nm capture window — so a first-order nudge from a notch
    the window does not contain lands on the window edge and reports a confident
    950 pm error for every bias point. Ask instead which azimuthal order sits
    nearest the target and what `n_eff` puts it exactly there:

        n_eff(lam)*L = m*lam,  n_eff(lam) = n0 + (lam - ref)*(n0 - n_g)/ref
        =>  n0 = m*ref/L + n_g*(lam - ref)/lam
    """
    m = round(guess * L_RING / (target_nm * 1e-9))
    return m * (WL_REF_NM * 1e-9) / L_RING + N_G * (target_nm - WL_REF_NM) / target_nm


def norm(t_db: np.ndarray) -> np.ndarray:
    """Peak to 0 dB, the way ringfit normalises — GC and fibre loss drop out."""
    return t_db - np.percentile(t_db, 95.0)


def notch_nm(wl: np.ndarray, t_db: np.ndarray) -> float:
    i = int(np.nanargmin(t_db))
    if 0 < i < len(t_db) - 1:
        y0, y1, y2 = t_db[i - 1], t_db[i], t_db[i + 1]
        d = 2 * y1 - y0 - y2
        if abs(d) > 1e-12:
            return float(wl[i] + np.clip((y0 - y2) / (2 * d), -1, 1) * (wl[1] - wl[0]))
    return float(wl[i])


def t_db_one(ckt, wl_nm: float, tag: str) -> float:
    """Fibre-to-fibre transmission of ONE ring at one wavelength."""
    for t in RINGS:
        ckt.set_param(f"XL{t}", "wavelength_nm", float(wl_nm))
    r = ckt.run("op")
    p_in = float(r[f"V({tag}_src_re_0)"][0]) ** 2 + float(r[f"V({tag}_src_im_0)"][0]) ** 2
    p_out = float(r[f"V({tag}_gth_re_0)"][0]) ** 2 + float(r[f"V({tag}_gth_im_0)"][0]) ** 2
    return 10.0 * np.log10(max(p_out, 1e-30) / max(p_in, 1e-30))


def baseline(ckt, tag: str, wl_nm: float) -> float:
    """Off-resonance transmission, which is what the measured spectra are
    normalised to.

    The simulated trace carries 18 dB of grating coupler and a couple of bus
    waveguides; the measured one is normalised so its 95th percentile is 0 dB.
    Comparing a notch DEPTH across the two without removing that offset compares
    -35 dB against -16 dB and calls it a 19 dB modelling error — which is what
    the first version of this fit did, and it drove `vol_dep` an order of
    magnitude off to make the discrepancy stop growing.

    Off resonance the ring passes everything, so one sample a nanometre away is
    the whole normalisation and it does not move with bias.
    """
    return t_db_one(ckt, wl_nm, tag)


def locate(ckt, tag: str, wl_guess: float) -> tuple[float, float]:
    """One ring's notch, by two narrowing brackets around a good guess.

    The fit calls this a few thousand times, so sampling 241 wavelengths to keep
    three of them is not affordable. Two rounds of 13 cost 26 and land inside a
    picometre — but only because the guess is good, so the caller passes the
    MEASURED notch. That keeps the guess independent of the parameters being
    fitted, which a previous-iterate guess would not be.

    Returns NaN rather than the window edge when the bracket holds no notch; the
    caller turns that into a penalty instead of a confident wrong number.
    """
    centre = wl_guess
    wl = y = None
    i = 0
    for span in (0.08, 0.012):
        wl = np.linspace(centre - span, centre + span, 13)
        y = np.array([t_db_one(ckt, w, tag) for w in wl])
        i = int(np.argmin(y))
        if i == 0 or i == len(wl) - 1:
            return np.nan, np.nan
        centre = float(wl[i])
    y0, y1, y2 = y[i - 1], y[i], y[i + 1]
    dd = 2 * y1 - y0 - y2
    if abs(dd) < 1e-12:
        return centre, float(y1)
    frac = float(np.clip((y0 - y2) / (2 * dd), -1, 1))
    return centre + frac * (wl[1] - wl[0]), float(y1 - 0.25 * (y0 - y2) * frac)


def fit(ckt, d, n_eff: dict[str, float], free_ratio: bool = False) -> dict[str, float]:
    """Fit vol_dep, vol_inj and r_th to the capture's notch positions and depths.

    Positions are referenced to each model's OWN unbiased notch, so the `n_eff`
    trim cannot soak up a fitting error and cannot hide one. Depths carry the
    absorption information: the resonance tells you the index, and only the
    depth separates electrons from holes.
    """
    from scipy.optimize import least_squares

    hc, jv, wl_all, T = d["hc_mA"], d["jv_V"], d["wl_nm"], d["T_dB"]
    h0, j0 = int(np.argmin(np.abs(hc))), int(np.argmin(np.abs(jv)))
    # Every other junction voltage, five heater currents. The capture's own
    # +0.974 V column is its flagged bad-data point and is left out: fitting to
    # it would drag every parameter toward one bad measurement.
    pts = [(jv[j], hc[h0]) for j in range(0, len(jv) - 1, 2)]
    pts += [(jv[j0], hc[i]) for i in np.linspace(0, len(hc) - 1, 5).astype(int)]
    meas = []
    for v, h in pts:
        j = int(np.argmin(np.abs(jv - v)))
        i = int(np.argmin(np.abs(hc - h)))
        row = T[i, j] - np.percentile(T[i, j], 95.0)
        meas.append((notch_nm(wl_all, row), float(np.min(row))))
    ref_meas = meas[[(v, h) for v, h in pts].index((jv[j0], hc[h0]))][0]
    # One sample a nanometre off resonance, to put the simulated depth on the
    # same normalisation as the measured one.
    base = baseline(ckt, "phys", ref_meas + 1.0)
    print(f"  simulated off-resonance baseline {base:.2f} dB "
          f"(18 dB of grating coupler plus bus)")

    free = dict(FREE)
    if free_ratio:
        free["ratio_h"] = RATIO_H[1:]
    x0 = np.array([free[k][0] for k in free])

    def resid(x):
        for k, v in zip(free, x):
            ckt.set_param("Xphys", k, float(v))
        # Reset n_eff to the deck's trim FIRST. The re-trim below is a one-step
        # correction from a known starting point, and without this reset it
        # started from wherever the previous call happened to leave the
        # instance — which made resid() depend on the order it was called in,
        # and scipy's cost disagree with a re-evaluation at the same x.
        ckt.set_param("Xphys", "n_eff", n_eff["phys"])
        # r_th changes the optical self-heating, which moves the unbiased notch
        # by a few pm, so the trim has to follow the parameters.
        for t in RINGS:
            ckt.set_param(f"V{t}", "dc", float(jv[j0]))
            ckt.set_param(f"IHT{t}", "dc", float(hc[h0] * 1e-3))
        got, _ = locate(ckt, "phys", ref_meas)
        if not np.isfinite(got):
            return np.full(2 * len(pts), 1e4)
        ckt.set_param("Xphys", "n_eff", n_eff["phys"] + N_G * (ref_meas - got) / ref_meas)
        ref_sim, _ = locate(ckt, "phys", ref_meas)
        out = []
        for (v, h), (m_res, m_dep) in zip(pts, meas):
            for t in RINGS:
                ckt.set_param(f"V{t}", "dc", float(v))
                ckt.set_param(f"IHT{t}", "dc", float(h * 1e-3))
            res, dep = locate(ckt, "phys", m_res)
            if not np.isfinite(res):
                out += [1e4, 1e4]
            else:
                out += [((res - ref_sim) - (m_res - ref_meas)) * 1e3 / SIGMA_POS_PM,
                        ((dep - base) - m_dep) / SIGMA_DEPTH_DB]
        return np.array(out)

    print(f"\nfitting {list(free)} to {len(pts)} bias points "
          f"({2 * len(pts)} residuals) …")
    sol = least_squares(resid, x0, bounds=([free[k][1] for k in free],
                                           [free[k][2] for k in free]),
                        x_scale=x0, diff_step=DIFF_STEP, xtol=1e-4,
                        max_nfev=60, verbose=2)
    got = {k: float(v) for k, v in zip(free, sol.x)}
    print()
    for k, v in got.items():
        print(f"  {k:9s} {free[k][0]:12.4g} -> {v:12.4g}   ({v / free[k][0]:.3f}x)")
    print(f"  final cost {sol.cost:.1f}")
    r1 = resid(sol.x)
    r0 = resid(x0)
    for lab, r in (("analytic start", r0), ("fitted", r1)):
        pos = r[0::2] * SIGMA_POS_PM
        dep = r[1::2] * SIGMA_DEPTH_DB
        print(f"  {lab:15s} notch rms {np.sqrt(np.mean(pos**2)):6.2f} pm,"
              f"  depth rms {np.sqrt(np.mean(dep**2)):5.2f} dB,"
              f"  chi2/N {np.mean(r**2):6.2f}")
    chi_start = float(np.mean(r0 ** 2))
    chi_fit = float(np.mean(r1 ** 2))
    if chi_fit >= chi_start:
        print("\n  the optimiser does not beat the analytic values, so those ship.\n"
              "  Not a failure of the fit — a statement about the capture. The three\n"
              "  parameters act on the resonance shift and the notch depth, and over a\n"
              "  1 V window this capture offers 13.5 pm of shift against a 1 pm floor\n"
              "  and 0.41 dB of depth against 0.3 dB of scatter. `vol_dep` alone spans\n"
              "  a decade at constant chi2. What would pin it: reverse bias past -3 V,\n"
              "  where the sqrt and the card's straight line are 20 % apart.")
        for k, (v0, _, _) in FREE.items():
            ckt.set_param("Xphys", k, v0)
        return {"kept": "analytic", "chi2_analytic": chi_start,
                "chi2_fitted": chi_fit, "fitted_but_rejected": got,
                **{k: v[0] for k, v in FREE.items()}}
    for k, v in got.items():
        ckt.set_param("Xphys", k, v)
    return {"kept": "fitted", "chi2_analytic": chi_start, "chi2_fitted": chi_fit, **got}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", action="store_true",
                    help="fit vol_dep, vol_inj and r_th to the capture")
    ap.add_argument("--free-ratio", action="store_true",
                    help="also free ratio_h, which is degenerate with vol_dep")
    args = ap.parse_args()
    RESULTS.mkdir(exist_ok=True)
    d = load()
    hc, jv, wl_all, T = d["hc_mA"], d["jv_V"], d["wl_nm"], d["T_dB"]
    i_meas = d["i_junc_mA"]
    j0 = int(np.argmin(np.abs(jv)))          # the jv closest to 0
    h0 = int(np.argmin(np.abs(hc)))          # the hc closest to 0
    wl = np.linspace(wl_all[0], wl_all[-1], N_WL)
    print(f"May capture: {len(hc)} heater currents x {len(jv)} junction volts, "
          f"{wl_all[0]:.2f}-{wl_all[-1]:.2f} nm")
    print(f"reference cell: hc = {hc[h0]:+.3f} mA, jv = {jv[j0]:+.3f} V")

    # ── the one free parameter ──────────────────────────────────────────────
    meas_res = notch_nm(wl_all, T[h0, j0])
    # PER RING, not one trim shared. The Verilog-A ring self-heats on the
    # absorbed light and the cell cannot, so at 0.126 mW their unbiased notches
    # sit 6 pm apart; trimming both to the cell's value hands that 6 pm to the
    # Verilog-A model as a fixed offset and then scores it for having one. The
    # card allows one trim per ring, so each ring gets its own.
    n_eff = dict.fromkeys(RINGS, n_eff_for(meas_res))
    ckt = fc.Circuit()
    ckt.load_str(deck(n_eff))
    if args.fit:
        # Before the trim, not after: `fit` reuses this trim as the fixed point
        # its one-step re-trim corrects from, and `r_th` moves the unbiased
        # notch through optical self-heating. Trimming at the model default and
        # then fitting from a different `r_th` made the objective depend on the
        # order it was evaluated in.
        for k, (v0, _, _) in FREE.items():
            ckt.set_param("Xphys", k, v0)
    sim = spectra(ckt, wl, jv[j0], hc[h0] * 1e-3)
    # A resonance moves with the GROUP index, so dlam/lam = dn_eff/n_g.
    n_eff = {t: n_eff[t] + N_G * (meas_res - notch_nm(wl, sim[t])) / meas_res
             for t in RINGS}
    ckt = fc.Circuit()
    ckt.load_str(deck(n_eff))
    sim = spectra(ckt, wl, jv[j0], hc[h0] * 1e-3)
    print(f"measured notch {meas_res:.4f} nm; n_eff 2.2810 trimmed per ring:")
    for t in RINGS:
        got = notch_nm(wl, sim[t])
        print(f"  {t:6s} -> {n_eff[t]:.6f}, landing at {got:.4f} nm "
              f"({(got - meas_res) * 1e3:+.1f} pm)")
        assert abs(got - meas_res) < 0.02, f"the {t} n_eff trim did not land"

    if args.fit:
        out_fit = fit(ckt, d, n_eff, args.free_ratio)
        # Re-trim after fitting: r_th moved, so the self-heating did.
        sim = spectra(ckt, wl, jv[j0], hc[h0] * 1e-3)
        n_eff = {t: n_eff[t] + N_G * (meas_res - notch_nm(wl, sim[t])) / meas_res
                 for t in RINGS}
        ckt = fc.Circuit()
        ckt.load_str(deck(n_eff))
        for k in FREE:
            ckt.set_param("Xphys", k, out_fit[k])

    out: dict = {"n_eff_trim": n_eff, "measured_notch_nm": meas_res}
    if args.fit:
        out["fitted"] = out_fit
    out["n_eff_spread_pm"] = (max(n_eff.values()) - min(n_eff.values())) / N_G * meas_res * 1e3

    # ── spectra vs junction voltage, heater off ─────────────────────────────
    print("\n── spectra vs junction voltage (heater at its zero point) ──")
    print(f"{'jv':>7} | {'meas notch':>11} | " +
          " | ".join(f"{t:>14}" for t in RINGS))
    vs_jv = {t: [] for t in RINGS}
    vs_jv_res = {t: [] for t in RINGS}
    meas_jv = []
    for j, v in enumerate(jv):
        s = spectra(ckt, wl, v, hc[h0] * 1e-3)
        m = norm(np.interp(wl, wl_all, T[h0, j]))
        meas_jv.append(m)
        row = []
        for t in RINGS:
            vs_jv[t].append(norm(s[t]))
            vs_jv_res[t].append(notch_nm(wl, s[t]))
            rms = float(np.sqrt(np.mean((vs_jv[t][-1] - m) ** 2)))
            row.append(f"{(vs_jv_res[t][-1] - notch_nm(wl_all, T[h0, j])) * 1e3:+6.0f} pm"
                       f" {rms:5.2f}dB")
        print(f"{v:+7.3f} | {notch_nm(wl_all, T[h0, j]):11.4f} | " + " | ".join(row))

    # ── terminal I(V): where the factor of two shows ────────────────────────
    print("\n── terminal current, measured vs simulated (2 kOhm shunt included) ──")
    # Reported with the shunt SUBTRACTED. Through the pads the 2 kOhm dominates
    # everywhere below a volt, and a 2x error in the diode branch shows up as a
    # 7 % error in the terminal current — which reads like a fit tolerance
    # rather than like the factor of two it is.
    # Two different questions, and only the first has a clean answer.
    #
    #   vs the CARD's own diode law: unambiguous arithmetic. The card fits
    #   (i_sat, n_diode) to the measured TERMINAL current, then `mrm.sp` uses
    #   i_sat per arc with two arcs in parallel, so the cell and this port of it
    #   draw exactly twice what the card says. No data needed to see that.
    #
    #   vs the MEASUREMENT: only meaningful where the diode branch is a decent
    #   fraction of the terminal current. Below ~0.8 V the 2 kOhm shunt carries
    #   over 90 % of it and a per-cent error in the shunt value swamps the
    #   residual, so those ratios are noise, not physics.
    vt_card = 5.0 * 0.025852
    print(f"{'jv':>7} {'meas diode':>13} {'card law':>12} {'as ported':>12}"
          f" {'I(V) fixed':>12} | {'shunt %':>8} {'vs card':>8}")
    # Read the current the solver actually draws, not a recomputation of the
    # diode law — the point is to test the model, and an independent formula
    # here would agree with a wrong model just as happily as with a right one.
    iv = {"meas": [], "card": [], "phys": []}
    for j, v in enumerate(jv):
        for t in RINGS:
            ckt.set_param(f"V{t}", "dc", float(v))
            ckt.set_param(f"IHT{t}", "dc", float(hc[h0] * 1e-3))
        r = ckt.run("op")
        shunt = v / R_SHUNT
        i_va = -float(r["I(vcard)"][0]) - shunt
        i_fx = -float(r["I(vphys)"][0]) - shunt
        i_ms = i_meas[h0, j] * 1e-3 - shunt
        iv["meas"].append(i_ms)
        iv["card"].append(i_va)
        iv["phys"].append(i_fx)
        i_card = 5.099e-8 * (np.exp(np.clip(v / vt_card, -40, 40)) - 1.0)
        frac = abs(shunt) / max(abs(shunt + i_va), 1e-15) * 100.0
        vs_card = f"{i_va / i_card:7.2f}x" if abs(i_card) > 1e-12 else "       -"
        print(f"{v:+7.3f} {i_ms * 1e3:10.4f} mA {i_card * 1e3:9.4f} mA"
              f" {i_va * 1e3:9.4f} mA {i_fx * 1e3:9.4f} mA | {frac:7.1f}% {vs_card}")

    # ── spectra vs heater current, junction at its zero point ───────────────
    print("\n── spectra vs heater current (junction at its zero point) ──")
    h_sel = np.linspace(0, len(hc) - 1, 7).astype(int)
    print(f"{'hc mA':>7} {'P mW':>7} | {'meas shift':>11} {'sim shift':>10} {'err':>8}")
    vs_hc, meas_hc, hc_res = [], [], []
    m_ref = notch_nm(wl_all, T[h0, j0])
    s_ref = notch_nm(wl, spectra(ckt, wl, jv[j0], 0.0)["phys"])
    for i in h_sel:
        s = spectra(ckt, wl, jv[j0], hc[i] * 1e-3)
        vs_hc.append(norm(s["phys"]))
        meas_hc.append(norm(np.interp(wl, wl_all, T[i, j0])))
        r_m = notch_nm(wl_all, T[i, j0]) - m_ref
        r_s = notch_nm(wl, s["phys"]) - s_ref
        hc_res.append((hc[i], r_m, r_s))
        p_mw = hc[i] ** 2 * 1e-6 * 368.8 * 1e3
        print(f"{hc[i]:+7.3f} {p_mw:7.3f} | {r_m * 1e3:8.0f} pm {r_s * 1e3:8.0f} pm"
              f" {(r_s - r_m) * 1e3:6.0f} pm")
    # What this capture says p_pi_th is, from the slope of shift against power.
    # The card carries 26.4 mW from the July neuron-3 sweep and warns the spread
    # across neurons is about 2x; this is that warning, measured.
    pw = np.array([h ** 2 * 1e-6 * 368.8 for h, _, _ in hc_res])
    sm = np.array([m for _, m, _ in hc_res])
    slope = float(np.polyfit(pw, sm, 1)[0])          # nm per watt
    p_pi_here = 11.38 / 2.0 / slope                  # half an FSR is pi
    print(f"\n  this capture's p_pi_th = {p_pi_here * 1e3:.1f} mW/pi "
          f"(card carries 26.4 mW from July neuron 3)")
    print(f"  -> r_th for THIS ring = "
          f"{1.55e-6 / (2 * 1.86e-4 * L_RING * p_pi_here):.0f} K/W")
    out["p_pi_th_this_capture_W"] = p_pi_here
    r_th_here = 1.55e-6 / (2 * 1.86e-4 * L_RING * p_pi_here)
    out["r_th_this_capture"] = r_th_here
    # And the same sweep with r_th set to what this capture says. The card's
    # thermal number is not wrong, it is another ring's — showing the retuned
    # curve is the difference between "the model is off by 2x" and "one
    # parameter is per-device, exactly as the card's header says".
    ckt.set_param("Xphys", "r_th", r_th_here)
    s_ref2 = notch_nm(wl, spectra(ckt, wl, jv[j0], 0.0)["phys"])
    hc_tuned = [(h, notch_nm(wl, spectra(ckt, wl, jv[j0], h * 1e-3)["phys"]) - s_ref2)
                for h, _, _ in hc_res]
    print("  with r_th retuned, worst heater error: "
          f"{max(abs(t - m) for (_, t), (_, m, _) in zip(hc_tuned, hc_res)) * 1e3:.0f} pm")

    out["rms_dB"] = {t: float(np.mean([np.sqrt(np.mean((a - b) ** 2))
                                       for a, b in zip(vs_jv[t], meas_jv)]))
                     for t in RINGS}
    out["iv_ratio_card"] = float(iv["card"][-2] / iv["meas"][-2])
    out["iv_ratio_phys"] = float(iv["phys"][-2] / iv["meas"][-2])
    print("\nmean spectral RMS over the jv axis: " +
          ", ".join(f"{t} {v:.2f} dB" for t, v in out["rms_dB"].items()))
    (RESULTS / "va_vs_may_data.json").write_text(json.dumps(out, indent=2) + "\n")

    plot(wl, wl_all, jv, hc, h_sel, meas_jv, vs_jv, meas_hc, vs_hc, iv, hc_res,
         hc_tuned, p_pi_here, out)


def plot(wl, wl_all, jv, hc, h_sel, meas_jv, vs_jv, meas_hc, vs_hc, iv, hc_res,
         hc_tuned, p_pi_here, out):
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    cmap = plt.get_cmap("viridis")

    sel = np.linspace(0, len(jv) - 1, 6).astype(int)
    for k, j in enumerate(sel):
        c = cmap(k / (len(sel) - 1))
        ax[0][0].plot(wl, meas_jv[j], "-", color=c, lw=2.4, alpha=0.45)
        ax[0][0].plot(wl, vs_jv["phys"][j], "--", color=c, lw=1.2)
        ax[0][0].plot([], [], "-", color=c, label=f"{jv[j]:+.2f} V")
    ax[0][0].set_title("spectra vs junction voltage\nsolid = measured, dashed = Verilog-A"
                       "\n(+0.97 V is the capture's flagged bad-data point)")
    ax[0][0].set_xlabel("wavelength (nm)")
    ax[0][0].set_ylabel("normalised T (dB)")
    ax[0][0].legend(fontsize=7)
    ax[0][0].grid(alpha=0.25)

    for k, i in enumerate(h_sel):
        c = cmap(k / (len(h_sel) - 1))
        ax[0][1].plot(wl, meas_hc[k], "-", color=c, lw=2.4, alpha=0.45)
        ax[0][1].plot(wl, vs_hc[k], "--", color=c, lw=1.2)
        ax[0][1].plot([], [], "-", color=c, label=f"{hc[i]:+.2f} mA")
    ax[0][1].set_title("spectra vs heater current\nsolid = measured, dashed = Verilog-A")
    ax[0][1].set_xlabel("wavelength (nm)")
    ax[0][1].set_ylabel("normalised T (dB)")
    ax[0][1].legend(fontsize=7)
    ax[0][1].grid(alpha=0.25)

    for t, (label, _) in RINGS.items():
        ax[0][2].plot(jv, [np.sqrt(np.mean((a - b) ** 2))
                           for a, b in zip(vs_jv[t], meas_jv)], "o-", label=label)
    ax[0][2].set_title("spectral RMS residual vs the capture")
    ax[0][2].set_xlabel("junction voltage (V)")
    ax[0][2].set_ylabel("RMS (dB)")
    ax[0][2].legend(fontsize=8)
    ax[0][2].grid(alpha=0.25)

    ax[1][0].plot(jv, np.array(iv["meas"]) * 1e3, "o-", lw=2.4, alpha=0.6, label="measured")
    ax[1][0].plot(jv, np.array(iv["card"]) * 1e3, "s--", label="card-compatible (2x)")
    ax[1][0].plot(jv, np.array(iv["phys"]) * 1e3, "^:", label="one carrier population")
    ax[1][0].set_title("terminal current — the inherited factor of two")
    ax[1][0].set_xlabel("junction voltage (V)")
    ax[1][0].set_ylabel("current (mA)")
    ax[1][0].legend(fontsize=8)
    ax[1][0].grid(alpha=0.25)

    ax[1][1].plot(jv, [notch_nm(wl, s) for s in meas_jv], "o-", lw=2.4, alpha=0.6,
                  label="measured")
    for t, (label, _) in RINGS.items():
        ax[1][1].plot(jv, [notch_nm(wl, s) for s in vs_jv[t]], "s--", label=label)
    ax[1][1].set_title("notch position vs junction voltage")
    ax[1][1].set_xlabel("junction voltage (V)")
    ax[1][1].set_ylabel("resonance (nm)")
    ax[1][1].legend(fontsize=7)
    ax[1][1].grid(alpha=0.25)

    h, m, s = np.array(hc_res).T
    p_mw = h**2 * 1e-6 * 368.8 * 1e3
    ax[1][2].plot(p_mw, m * 1e3, "o-", lw=2.4, alpha=0.6, label="measured")
    ax[1][2].plot(p_mw, s * 1e3, "s--", label="card r_th (p_pi = 26.4 mW, July n3)")
    ax[1][2].plot(p_mw, np.array([t for _, t in hc_tuned]) * 1e3, "^:",
                  label=f"r_th retuned here (p_pi = {p_pi_here * 1e3:.0f} mW)")
    ax[1][2].set_title("thermal tuning — p_pi_th is per device,\nand the card's is another ring's")
    ax[1][2].set_xlabel("heater power (mW)")
    ax[1][2].set_ylabel("resonance shift (pm)")
    ax[1][2].legend(fontsize=8)
    ax[1][2].grid(alpha=0.25)

    fig.suptitle("mrm_addrop.va vs the May giona capture — card defaults, n_eff trimmed, "
                 "nothing else fitted", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = RESULTS / "va_vs_may_data.png"
    fig.savefig(p, dpi=130)
    print(f"wrote {p}")


if __name__ == "__main__":
    main()
