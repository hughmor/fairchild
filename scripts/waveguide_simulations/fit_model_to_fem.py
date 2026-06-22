#!/usr/bin/env python3
"""fit_model_to_fem.py — can the circuit PN-shifter FORM reproduce the FEM curve?

SCOPE / WHAT THIS DOES NOT DO: this does **not** diagnose the chip forward-bias
fit failure. The chip-fit question can only be settled by the measured data in
the fitting loop — the FEM curve is itself a model, and agreement (or not)
between two models adds no information about the device. This script asks the
narrower, useful question: how well can each circuit constitutive *form*
reproduce the femwell/Soref-Bennett extraction in `pn_extracted.json`?

Two robust takeaways (which survive the caveats below):
  • The FEM extraction is a coarse STAIRCASE — Δn/Δα are piecewise-constant over
    ~0.1-0.7 V-wide steps (binary depletion mask on a coarse mesh, see
    pn_modulator.py:196) and even non-monotonic. So it is usable for rough
    init/bounds (saturation level, near-0V slopes) but NOT as a clean fit target,
    and a polynomial's apparent "win" here is partly chasing quantization.
  • Deep reverse SATURATES (full depletion of the mode overlap). A linear
    `da_dv_rev·|V|` / `dn_dv_rev·V` term cannot represent saturation; a
    depletion-width / sublinear term (~(1−V/Vbi)^(1−m)) can. That motivates a
    saturating depletion term + the existing exp injection — not a higher poly.

CAVEAT: the V grid is −4..+0.7 V, ~85% reverse, so these fits are
reverse-dominated and say little about forward bias (the fitter drops the
injection term to chase the reverse staircase). Forward bias is the chip data's
job. We score reverse and forward regions separately to make that explicit.

Forms compared (Δn_eff(V) and Δα(V) reproduced; reverse/forward RMS reported):

  F  fc_pn_ps_full's exact form — Δn and Δα share ONE injection factor (e−1)
     with a single ideality n:
        Δn = a_rev·V − a_inj·(e−1)₊
        Δα = b_rev·(−V)₊ + b_inj·(e−1)₊        e = exp(V/(n·V_T))
  P  polynomials Δn=Σc_k V^k (the "just raise the order" idea), fit per-region.
  S  decoupled-injection form — Δn and Δα get INDEPENDENT ideality (n_n, n_a),
     i.e. independent forward curvature (the minimal physics-motivated upgrade
     that mirrors Soref-Bennett's different carrier-density exponents for index
     vs loss).

Run:  .venv/bin/python scripts/waveguide_simulations/fit_model_to_fem.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
JSON = HERE / "pn_modulator" / "pn_extracted.json"

K_OVER_Q = 8.617333262e-5  # V/K


def load_fem():
    d = json.loads(JSON.read_text())
    T = d.get("sim_cfg", {}).get("temperature_K") or d.get("sim_cfg", {}).get("T") or 300.15
    vt = K_OVER_Q * float(T)
    dn = d["delta_neff_vs_v"]
    da = d["delta_alpha_vs_v"]
    V = np.asarray(dn["V"], float)
    dn_eff = np.asarray(dn["dn_eff"], float)
    Va = np.asarray(da["V"], float)
    alpha = np.asarray(da["alpha_per_m"], float)
    assert np.allclose(V, Va), "Δn and Δα are on different V grids"
    return V, dn_eff, alpha, vt, T


def inj_factor(V, n, vt):
    """(e−1)₊ with the same clamp the Rust model uses."""
    arg = np.clip(V / (n * vt), -40.0, 40.0)
    return np.clip(np.exp(arg) - 1.0, 0.0, None)


def region_rms(V, resid):
    """RMS of a residual split into reverse (V<0) and forward (V>0)."""
    rev = V < 0
    fwd = V > 0
    rms = lambda m: float(np.sqrt(np.mean(resid[m] ** 2))) if m.any() else float("nan")
    return rms(rev), rms(fwd), float(np.sqrt(np.mean(resid ** 2)))


def fit_full(V, dn_eff, alpha, vt):
    """Model F: shared injection factor (one ideality n)."""
    dn_s = np.std(dn_eff) or 1.0
    da_s = np.std(alpha) or 1.0

    def resid(p):
        a_rev, a_inj, b_rev, b_inj, n = p
        ij = inj_factor(V, n, vt)
        dn_m = a_rev * V - a_inj * ij
        da_m = b_rev * np.clip(-V, 0, None) + b_inj * ij
        return np.concatenate([(dn_m - dn_eff) / dn_s, (da_m - alpha) / da_s])

    p0 = [5e-5, 1e-4, 30.0, 30.0, 1.0]
    lb = [-1e-2, -1e-2, -1e6, -1e6, 0.5]
    ub = [1e-2, 1e-2, 1e6, 1e6, 4.0]
    r = least_squares(resid, p0, bounds=(lb, ub), max_nfev=20000)
    a_rev, a_inj, b_rev, b_inj, n = r.x
    ij = inj_factor(V, n, vt)
    dn_m = a_rev * V - a_inj * ij
    da_m = b_rev * np.clip(-V, 0, None) + b_inj * ij
    return dict(params=dict(a_rev=a_rev, a_inj=a_inj, b_rev=b_rev, b_inj=b_inj, n=n),
                dn_m=dn_m, da_m=da_m)


def fit_decoupled(V, dn_eff, alpha, vt):
    """Model S: independent ideality for Δn and Δα (decoupled forward curvature)."""
    dn_s = np.std(dn_eff) or 1.0
    da_s = np.std(alpha) or 1.0

    def resid(p):
        a_rev, a_inj, n_n, b_rev, b_inj, n_a = p
        dn_m = a_rev * V - a_inj * inj_factor(V, n_n, vt)
        da_m = b_rev * np.clip(-V, 0, None) + b_inj * inj_factor(V, n_a, vt)
        return np.concatenate([(dn_m - dn_eff) / dn_s, (da_m - alpha) / da_s])

    p0 = [5e-5, 1e-4, 1.0, 30.0, 30.0, 1.0]
    lb = [-1e-2, -1e-2, 0.3, -1e6, -1e6, 0.3]
    ub = [1e-2, 1e-2, 6.0, 1e6, 1e6, 6.0]
    r = least_squares(resid, p0, bounds=(lb, ub), max_nfev=20000)
    a_rev, a_inj, n_n, b_rev, b_inj, n_a = r.x
    dn_m = a_rev * V - a_inj * inj_factor(V, n_n, vt)
    da_m = b_rev * np.clip(-V, 0, None) + b_inj * inj_factor(V, n_a, vt)
    return dict(params=dict(a_rev=a_rev, a_inj=a_inj, n_n=n_n,
                            b_rev=b_rev, b_inj=b_inj, n_a=n_a),
                dn_m=dn_m, da_m=da_m)


def fit_poly(V, y, order):
    c = np.polyfit(V, y, order)
    return np.polyval(c, V)


def main():
    V, dn_eff, alpha, vt, T = load_fem()
    print(f"FEM curve: {len(V)} pts, V∈[{V.min():.2f},{V.max():.2f}] V,  "
          f"T={T} K → V_T={vt*1e3:.3f} mV")
    print(f"  Δn_eff range: [{dn_eff.min():.3e}, {dn_eff.max():.3e}]")
    print(f"  Δα   range: [{alpha.min():.3e}, {alpha.max():.3e}] /m\n")

    def report(name, dn_m, da_m):
        dn_rev, dn_fwd, dn_all = region_rms(V, dn_m - dn_eff)
        da_rev, da_fwd, da_all = region_rms(V, da_m - alpha)
        # forward residual as % of forward data span
        fwd = V > 0
        dn_span = np.ptp(dn_eff[fwd]) or 1.0
        da_span = np.ptp(alpha[fwd]) or 1.0
        print(f"{name:30s}  Δn RMS rev/fwd = {dn_rev:.3e}/{dn_fwd:.3e}"
              f" (fwd {100*dn_fwd/dn_span:5.1f}% span)   "
              f"Δα RMS rev/fwd = {da_rev:.3e}/{da_fwd:.3e}"
              f" (fwd {100*da_fwd/da_span:5.1f}% span)")

    F = fit_full(V, dn_eff, alpha, vt)
    S = fit_decoupled(V, dn_eff, alpha, vt)
    report("F  fc_pn_ps_full (shared n)", F["dn_m"], F["da_m"])
    report("S  decoupled n_n,n_a", S["dn_m"], S["da_m"])
    # Polynomials fit Δn and Δα independently (best case for the poly idea).
    for order in (3, 5, 7):
        dn_p = fit_poly(V, dn_eff, order)
        da_p = fit_poly(V, alpha, order)
        report(f"P  polynomial order {order}", dn_p, da_p)

    print("\nFitted params:")
    print("  F:", {k: f"{v:.4g}" for k, v in F["params"].items()})
    print("  S:", {k: f"{v:.4g}" for k, v in S["params"].items()})

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 2, figsize=(12, 8))
        for j, (y, models, lbl, unit) in enumerate([
            (dn_eff, [("F", F["dn_m"]), ("S", S["dn_m"]),
                      ("poly5", fit_poly(V, dn_eff, 5))], "Δn_eff", ""),
            (alpha, [("F", F["da_m"]), ("S", S["da_m"]),
                     ("poly5", fit_poly(V, alpha, 5))], "Δα", " (/m)")]):
            ax[0, j].plot(V, y, "k.", ms=4, label="FEM (truth)")
            for nm, ym in models:
                ax[0, j].plot(V, ym, label=nm, lw=1.4)
            ax[0, j].axvline(0, color="gray", lw=0.6)
            ax[0, j].set_title(f"{lbl}{unit} vs V"); ax[0, j].legend(fontsize=8)
            ax[0, j].grid(alpha=0.3)
            for nm, ym in models:
                ax[1, j].plot(V, ym - y, label=nm, lw=1.2)
            ax[1, j].axvline(0, color="gray", lw=0.6)
            ax[1, j].set_title(f"{lbl} residual"); ax[1, j].set_xlabel("V (V)")
            ax[1, j].grid(alpha=0.3); ax[1, j].legend(fontsize=8)
        fig.suptitle("PN phase-shifter model vs FEM (Soref-Bennett) ground truth")
        fig.tight_layout()
        out = HERE / "model_expressiveness.png"
        fig.savefig(out, dpi=120)
        print(f"\nwrote {out}")
    except ImportError:
        pass


if __name__ == "__main__":
    main()
