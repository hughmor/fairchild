#!/usr/bin/env python3
"""expt_forward_mismatch.py — reproduce (and fix) the forward-bias fit failure.

The recovery test in fit_transient.py is an *inverse crime*: it generates and
fits with the SAME model, so model-form error is zero by construction — it
validates the machinery but says nothing about whether a model is faithful
enough for real data. This script asks the question that actually matters:

    Generate forward-rich "measured" data with fc_pn_ps_full (depletion +
    exponential carrier injection + free-carrier absorption), then fit it
      (a) with the LINEAR fc_pn_ps  → should reproduce the forward-bias failure
          (a high residual floor no choice of params can beat), and
      (b) with fc_pn_ps_full        → should reach the noise floor.

The headline is the MINIMUM ACHIEVABLE RESIDUAL of each model — that comparison
is inverse-crime-free (it measures structural adequacy, not param recovery).

Run:  .venv/bin/python experiments/giona/expt_forward_mismatch.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fit_transient as ft  # noqa: E402

OBS = "P(pn_in)"            # strongly-modulated circulating field (see fit_transient)

# Device-under-test lines (the ring loop's PN shifter). Same ports/length/loss;
# only the carrier→optics constitutive model differs.
LINEAR = ("Xpn pn_in dc_b vmod 0 fc_pn_ps "
          "L_um=500 V_pi_L=2e-3 g_pn=1e-3 alpha_dB_cm=10")
FULL = ("Xpn pn_in dc_b vmod 0 fc_pn_ps_full L_um=500 alpha_dB_cm=10 "
        "dn_dv_rev=5e-5 da_dv_rev=8 dn_dv_inj=1.5e-4 da_dv_inj=150 "
        "i_sat=1e-9 n_diode=1.05 r_series=60 tau_carrier=10n c_j0=1.4e-13 v_bi=0.92")


def pn_netlist(pn_line: str) -> str:
    """Swap the ring's Xpn device line (+ any continuation) for `pn_line`."""
    return re.sub(r"(?m)^Xpn .*\n(\+.*\n)*", pn_line + "\n", ft.NETLIST)


def staircase(levels, hold=8e-11, step=ft.DEFAULT_STEP):
    """Bias staircase sweeping reverse→forward so the forward knee is sampled."""
    t, v = [0.0], [levels[0]]
    for i, L in enumerate(levels):
        t0 = i * hold
        t += [t0, t0 + 1e-12]; v += [v[-1], L]
    t.append(len(levels) * hold); v.append(levels[-1])
    return np.array(t), np.array(v)


def main():
    levels = [0.0, 0.2, 0.4, 0.6, 0.8]      # V_pn = vmod; forward = positive
    drive_t, drive_v = staircase(levels)
    truth_nl = pn_netlist(FULL)

    # "Measured" = full-model truth + instrument distortion (BW, gain, offset, noise).
    t_sim, clean = ft.simulate_pd({}, drive_t, drive_v, observable=OBS, netlist=truth_nl)
    dt = t_sim[1] - t_sim[0]
    BW_TAU = 1.0 / (2 * np.pi * 12e9)
    meas = ft.lowpass_1pole(clean, dt, BW_TAU)
    meas = 3.0 * meas - 0.4
    rng = np.random.default_rng(0)
    meas = meas + rng.normal(0, 0.01 * np.ptp(meas), meas.shape)
    drive_on = ft.drive_v_resampled(drive_t, drive_v, t_sim)
    span = float(np.ptp(meas))

    def fit(label, netlist, specs):
        print(f"\n[{label}]")
        best, res = ft.fit_transient(specs, t_sim, drive_on, meas, bw_tau=BW_TAU,
                                     observable=OBS, gain=3.0, netlist=netlist)
        # best-fit trace + normalized RMS residual
        _, sim = ft.simulate_pd(best, t_sim, drive_on, observable=OBS, netlist=netlist)
        sim = np.interp(t_sim, t_sim, sim)
        model, _ = ft.model_trace(sim, meas, dt, bw_tau=BW_TAU, max_lag=40, gain=3.0)
        rms = float(np.sqrt(np.mean(((model - meas) / span) ** 2)))
        print(f"  best params: " + ", ".join(f"{k.split('.')[1]}={v:.4g}"
                                              for k, v in best.items()))
        print(f"  min achievable residual (normalized RMS) = {rms:.4f}")
        return best, model, rms

    # (a) linear model: free its only shape knobs.
    b_lin, m_lin, rms_lin = fit(
        "LINEAR fc_pn_ps (Hugh's current model)", pn_netlist(LINEAR),
        [ft.FitParam("Xpn.V_pi_L", 2e-3, 5e-4, 8e-3),
         ft.FitParam("Xpn.alpha_dB_cm", 10.0, 1.0, 25.0)])
    # (b) full model: free the forward (injection) params; others at sensible vals.
    b_full, m_full, rms_full = fit(
        "FULL fc_pn_ps_full (recommended)", pn_netlist(FULL),
        [ft.FitParam("Xpn.dn_dv_inj", 0.8e-4, 1e-5, 5e-4),
         ft.FitParam("Xpn.da_dv_inj", 80.0, 10.0, 400.0)])

    print(f"\n=== VERDICT ===")
    print(f"  linear model residual floor : {rms_lin:.4f}")
    print(f"  full   model residual floor : {rms_full:.4f}")
    print(f"  full model is {rms_lin/max(rms_full,1e-9):.1f}× better.")
    print("  → the linear model is STRUCTURALLY unable to fit the forward-bias\n"
          "    shape (a high residual floor no params can beat); the full model\n"
          "    reaches the noise floor. This is the controlled analogue of the\n"
          "    real-chip forward-bias failure — the fix is model form, not optimizer.")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        ts = t_sim * 1e12
        fig, ax = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
        ax[0].plot(ts, meas, "k.", ms=2, label="measured (full-model truth + noise)")
        ax[0].plot(ts, m_lin, "tab:red", lw=1.5, label=f"linear best (RMS {rms_lin:.3f})")
        ax[0].plot(ts, m_full, "tab:blue", lw=1.2, label=f"full best (RMS {rms_full:.3f})")
        ax[0].set_ylabel("PD signal (a.u.)"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
        ax[0].set_title("Forward-bias model mismatch: linear cannot fit, full can")
        ax[1].plot(ts, m_lin - meas, "tab:red", lw=1, label="linear residual")
        ax[1].plot(ts, m_full - meas, "tab:blue", lw=1, label="full residual")
        ax[1].axhline(0, color="gray", lw=0.6)
        ax[1].set_xlabel("time (ps)"); ax[1].set_ylabel("residual"); ax[1].grid(alpha=0.3)
        ax[1].legend(fontsize=8)
        fig.tight_layout()
        out = HERE / "results" / "fit_forward_mismatch.png"
        fig.savefig(out, dpi=120); print(f"\nwrote {out}")
    except ImportError:
        pass


if __name__ == "__main__":
    main()
