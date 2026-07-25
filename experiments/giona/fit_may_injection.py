#!/usr/bin/env python3
"""fit_may_injection.py — fit the current-parametrized injection (dn_di/da_di)
against the May live-junction dataset (giona_neuron2_mod_joint_IV_spec).

The decoupled workflow the dn_di parametrization enables:
  1. pin (i_sat, n_diode, r_series) from the measured junction IV alone
     (ringfit.prefit_diode_iv, 2 kΩ PCB shunt subtracted),
  2. fit (dn_di, da_di) — Δn = −dn_di·I_fwd, Δα = da_di·I_fwd — as two linear
     coefficients from the forward-bias spectra (0 < JV ≤ 0.9 V, heater off),
     with everything else pinned from the May staged fit (alphas converted to
     post-loss-fix units: fitted-pre-fix / 2).

Run:  .venv/bin/python experiments/giona/fit_may_injection.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ringfit import (  # noqa: E402
    extract_data, load_sweep, prefit_diode_iv, wavelength_sweep, _spectrum_loss,
)

DATA = str(HERE / "data" / "giona_neuron2_mod_joint_IV_spec")
BASE_JSON = HERE / "results" / "giona_neuron7_pn_th_ps_full_fit.json"
OUT_JSON = HERE / "results" / "giona_dn_di_fit.json"
MODEL = "fc_pn_th_ps_full"


def base_params() -> tuple[dict, float]:
    """May staged-fit params, loss-fix-converted; returns (ps_params, kappa_l)."""
    p = json.loads(BASE_JSON.read_text())["params"]
    kappa_l = p.pop("kappa_l")
    p["alpha_db_cm"] = p["alpha_db_cm"] / 2.0   # pre-fix fitted → physical dB/cm
    # Old degenerate injection parametrization off; dn_di path replaces it.
    p["dn_dv_inj"] = 0.0
    p["da_dv_inj"] = 0.0
    return p, kappa_l


def main():
    sweep = load_sweep(DATA)
    sd = extract_data(sweep)
    print(f"sweep: {len(sd.hc_mA)} HC × {len(sd.jv_V)} JV (≤{sd.jv_V.max():.1f} V), "
          f"{sd.T_dB.shape[2]} λ-pts {sd.wl_nm[0]:.2f}–{sd.wl_nm[-1]:.2f} nm")

    # 1. diode electricals from the IV alone (2 kΩ shunt handled inside).
    i_sat, n_diode, r_series = prefit_diode_iv(sd)

    ps, kappa_l = base_params()
    ps.update(i_sat=i_sat, n_diode=n_diode, r_series=r_series)

    # forward-bias slice, heater off
    i0 = int(np.argmin(np.abs(sd.hc_mA)))
    jv_fwd = sd.jv_V[(sd.jv_V > 0) & (sd.jv_V <= 0.9)]
    jv_sel = jv_fwd[np.round(np.linspace(0, len(jv_fwd) - 1,
                                         min(6, len(jv_fwd)))).astype(int)]
    print(f"forward JV points: {np.round(jv_sel, 3)}")

    def spectra_loss(dn_di, da_di):
        p = dict(ps, dn_di=dn_di, da_di=da_di)
        loss = 0.0
        for jv in jv_sel:
            j = int(np.argmin(np.abs(sd.jv_V - jv)))
            sim = wavelength_sweep(MODEL, p, kappa_l, sd.wl_nm, v_pn=float(jv))
            loss += _spectrum_loss(sim, sd.T_dB[i0, j])
        return loss / len(jv_sel)

    # baseline: no injection optics at all
    base_loss = spectra_loss(0.0, 0.0)
    print(f"baseline loss (dn_di=da_di=0): {base_loss:.5f}")

    # 2. DE in log10 space (both coefficients span decades), then polish.
    def cost(x):
        return spectra_loss(10.0 ** x[0], 10.0 ** x[1])

    de = differential_evolution(cost, [(-4.0, 2.0), (-2.0, 6.0)], seed=0,
                                maxiter=30, popsize=8, tol=1e-6, polish=False)
    pol = minimize(cost, de.x, method="Nelder-Mead",
                   options=dict(xatol=1e-3, fatol=1e-6))
    dn_di, da_di = 10.0 ** pol.x[0], 10.0 ** pol.x[1]
    print(f"fitted: dn_di={dn_di:.4g} /A  da_di={da_di:.4g} Np/m/A  "
          f"loss={pol.fun:.5f}  (baseline {base_loss:.5f})")

    out = dict(model=MODEL, source="May neuron2 dataset, post-loss-fix engine",
               i_sat=i_sat, n_diode=n_diode, r_series=r_series,
               dn_di=dn_di, da_di=da_di,
               loss=pol.fun, baseline_loss=base_loss)
    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(f"saved → {OUT_JSON}")

    # validation overlay: forward spectra, measured vs fitted vs no-injection
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, len(jv_sel), figsize=(3.1 * len(jv_sel), 3.6),
                            sharey=True)
    for ax, jv in zip(np.atleast_1d(axs), jv_sel):
        j = int(np.argmin(np.abs(sd.jv_V - jv)))
        ax.plot(sd.wl_nm, sd.T_dB[i0, j], "k.", ms=2, label="meas")
        p = dict(ps, dn_di=dn_di, da_di=da_di)
        ax.plot(sd.wl_nm, wavelength_sweep(MODEL, p, kappa_l, sd.wl_nm,
                                           v_pn=float(jv)) +
                np.median(sd.T_dB[i0, j] - wavelength_sweep(
                    MODEL, p, kappa_l, sd.wl_nm, v_pn=float(jv))),
                "r-", lw=1, label="fit (dn_di)")
        p0 = dict(ps, dn_di=0.0, da_di=0.0)
        s0 = wavelength_sweep(MODEL, p0, kappa_l, sd.wl_nm, v_pn=float(jv))
        ax.plot(sd.wl_nm, s0 + np.median(sd.T_dB[i0, j] - s0), "b--", lw=0.8,
                label="no injection")
        ax.set_title(f"JV={jv:.2f} V", fontsize=9)
        ax.set_xlabel("λ (nm)")
    np.atleast_1d(axs)[0].set_ylabel("T (dB)")
    np.atleast_1d(axs)[0].legend(fontsize=7)
    fig.suptitle("May neuron2 forward-bias spectra — dn_di/da_di fit")
    fig.tight_layout()
    png = HERE / "results" / "may_dn_di_fit.png"
    fig.savefig(png, dpi=120)
    print(f"plot → {png}")


if __name__ == "__main__":
    main()
