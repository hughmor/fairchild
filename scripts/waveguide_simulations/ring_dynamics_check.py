#!/usr/bin/env python3
"""ring_dynamics_check.py — does the MRR carry real cavity dynamics?

The upstream correctness gate for time-domain fitting. A silicon micro-ring has
a photon lifetime τ_ph = Q/ω (tens of ps at Q~1e4) that is a real fraction of a
high-speed bit. If fairchild's ring relaxes *instantly* in transient, an
optimizer fitting dynamics will misattribute the missing cavity lag to carrier
lifetime / RC and return wrong device params.

This drives the all-pass ring in examples/photonic/native_mrr_modulator.sp with
a FAST detuning step at a fine timestep, and checks:
  1. with waveguide_delay ON  → through-port relaxes with τ ≈ τ_ph (analytic),
  2. with waveguide_delay OFF → response is ~instantaneous (no cavity lag).

τ_ph is computed analytically from the netlist parameters (no fit needed to know
the target). Run:  .venv/bin/python scripts/waveguide_simulations/ring_dynamics_check.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import fairchild as fc

C = 2.99792458e8
LAM = 1550e-9
N_G = 4.2
L_RING = 500e-6          # PN-shifter length = ring round trip
ALPHA_DB_CM = 10.0
KAPPA_L = 0.336

def _quiet(nl: str) -> str:
    """Drop the harmless wavelength_nm token from fc_waveguide lines (it's
    ignored by the waveguide model and otherwise spams a warning per run)."""
    out = []
    for ln in nl.splitlines():
        if "fc_waveguide" in ln:
            ln = ln.replace(" wavelength_nm=1550", "")
        out.append(ln)
    return "\n".join(out)


NETLIST = _quiet((Path(__file__).resolve().parents[2]
                  / "examples" / "photonic" / "native_mrr_modulator.sp").read_text())

# The cavity observable is the CIRCULATING field (pn_in), not the (under-coupled,
# nearly flat) through port — the stored ring energy is what relaxes over τ_ph.
CAVITY = "pn_in"


def analytic_tau_ph():
    t_rt = N_G * L_RING / C
    a = 10 ** (-(ALPHA_DB_CM * (L_RING * 100.0)) / 20.0)   # round-trip amplitude (loss)
    t = np.cos(KAPPA_L)                                     # self-coupling amplitude
    r = a * t
    tau = -t_rt / (2.0 * np.log(r))                        # power/energy e-folding
    Q = (2 * np.pi * C / LAM) * tau
    fsr = C / (N_G * L_RING)
    return dict(t_rt=t_rt, a=a, t=t, r=r, tau_ph=tau, Q=Q, fsr=fsr)


# WaveformSource is exported at top level after the __init__ fix; fall back to
# the compiled submodule so this runs against an older install too.
WS = getattr(fc, "WaveformSource", None) or fc.fairchild.WaveformSource


def cavity_power(r):
    re = np.asarray(r[f"V({CAVITY}_re_0)"]); im = np.asarray(r[f"V({CAVITY}_im_0)"])
    return re * re + im * im


def find_resonance():
    """Sweep vmod in op; return (V@max circulating, V@min circulating)."""
    ckt = fc.Circuit(); ckt.load_str(NETLIST)
    Vs = np.linspace(0.0, 4.0, 41)
    P = np.array([cavity_power(
        (ckt.set_param("Vmod", "dc", float(v)), ckt.run("op"))[1])[0] for v in Vs])
    return float(Vs[P.argmax()]), float(Vs[P.argmin()]), Vs, P


def step_response(v_lo, v_hi, delay, step=5e-13, stop=4e-10, t0=4e-11):
    """Fast step v_lo→v_hi at t0; return (t, circulating power).

    NB: waveguide_delay MUST go through run()'s kwarg, not an appended
    `.options` line — the example netlists end with `.end`, after which the
    parser ignores everything (a subtle footgun the fitting harness must avoid).
    """
    ckt = fc.Circuit(); ckt.load_str(NETLIST)
    # Exact step waveform via PWL (1 ps edge), well below τ_ph.
    t = np.array([0.0, t0, t0 + 1e-12, stop])
    v = np.array([v_lo, v_lo, v_hi, v_hi])
    ckt.set_source("Vmod", WS(t, v))
    r = ckt.run("tran", step=step, stop=stop, variable_step=False,
                method="gear", waveguide_delay=bool(delay))
    return np.asarray(r.time()), cavity_power(r)


def fit_tau(t, p, t0):
    """Relaxation time of p(t) after t0, measured as the 1/e-excursion crossing.

    Robust to the ~t_rt round-trip ripple that sits on the decay envelope (a
    log-linear fit would be biased by it) and to the instantaneous edge spike
    (we take the post-edge level a couple ps in, not the spike)."""
    m = t >= t0 + 2e-12
    tt, pp = t[m], p[m]
    if len(pp) < 8:
        return float("nan")
    # Smooth out the ~t_rt round-trip ripple (else a ripple trough trips the
    # crossing early). Window ≈ 2 round trips.
    dt = tt[1] - tt[0]
    win = max(3, int(round(2 * (N_G * L_RING / C) / dt)) | 1)
    if len(pp) > win:
        k = np.ones(win) / win
        pp = np.convolve(pp, k, mode="same")
        edge = win // 2                                  # drop convolution edges
        tt, pp = tt[edge:-edge], pp[edge:-edge]
    p0 = float(np.median(pp[:6]))                       # post-edge level
    p_inf = float(np.median(pp[-max(6, len(pp) // 8):]))  # settled level
    if abs(p0 - p_inf) < 1e-9:
        return float("nan")                              # no relaxation (instant)
    target = p_inf + (p0 - p_inf) / np.e
    s = np.sign(pp - target)
    cross = np.where(np.diff(s) != 0)[0]
    if len(cross) == 0:
        return float("nan")
    i = cross[0]
    frac = (target - pp[i]) / (pp[i + 1] - pp[i]) if pp[i + 1] != pp[i] else 0.0
    return (tt[i] + frac * (tt[i + 1] - tt[i])) - tt[0]


def main():
    A = analytic_tau_ph()
    print(f"Analytic ring: t_rt={A['t_rt']*1e12:.2f} ps  FSR={A['fsr']/1e9:.1f} GHz  "
          f"r={A['r']:.4f}  τ_ph={A['tau_ph']*1e12:.1f} ps  Q={A['Q']:.3g}\n")

    v_on, v_off, Vs, Psw = find_resonance()
    print(f"Resonance sweep ({CAVITY}, circulating): max @ V={v_on:.2f} "
          f"({Psw.max()*1e3:.3f} mW)  min @ V={v_off:.2f} ({Psw.min()*1e3:.3f} mW)  "
          f"build-up {Psw.max()/max(Psw.min(),1e-12):.1f}×\n")

    results = {}
    for delay in (True, False):
        t, p = step_response(v_on, v_off, delay)
        tau = fit_tau(t, p, 4e-11)
        results[delay] = (t, p, tau)
        tag = "ON " if delay else "OFF"
        print(f"waveguide_delay={tag}  measured relaxation τ = "
              f"{tau*1e12:7.2f} ps" if tau == tau else
              f"waveguide_delay={tag}  measured relaxation τ = ~instant (unresolved)")

    tau_on = results[True][2]
    ratio = tau_on / A['tau_ph'] if tau_on == tau_on else float("nan")
    print(f"\nGATE: τ_on/τ_ph = {ratio:.2f}  (≈1 ⇒ cavity dynamics correct)")
    print(f"      delay-OFF should be ≪ τ_ph (near-instant): "
          f"τ_off = {results[False][2]*1e12 if results[False][2]==results[False][2] else 0:.2f} ps")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
        for delay, c in ((True, "tab:blue"), (False, "tab:red")):
            t, p, tau = results[delay]
            ax[0].plot(t * 1e12, p * 1e3, c, lw=1.4,
                       label=f"delay {'ON' if delay else 'OFF'} (τ={tau*1e12:.1f} ps)")
        ax[0].axvline(40, color="gray", lw=0.6, ls="--")
        ax[0].set_xlabel("time (ps)"); ax[0].set_ylabel("circulating power (mW)")
        ax[0].set_title(f"Step detune  (τ_ph={A['tau_ph']*1e12:.0f} ps analytic)")
        ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
        ax[1].plot(Vs, Psw * 1e3, "k.-", lw=1)
        ax[1].set_xlabel("Vmod (V)"); ax[1].set_ylabel("circulating power (mW)")
        ax[1].set_title("Resonance sweep (.op)"); ax[1].grid(alpha=0.3)
        fig.tight_layout()
        out = Path(__file__).resolve().parent / "ring_dynamics_check.png"
        fig.savefig(out, dpi=120); print(f"\nwrote {out}")
    except ImportError:
        pass


if __name__ == "__main__":
    main()
