#!/usr/bin/env python3
"""fit_jul_neuron3.py — staged fit of the neuron3 mod-bank ring to the 2026-07-10 sweep.

Data verdict first (see --report): in this capture the PN junction is
electrically RESISTIVE (IV is a 15.9 kΩ straight line to 0.1 µA RMS over
−6..+4 V — no diode turn-on) and the EO response is tiny and LINEAR
(~+1.8 pm/V, monotonic, not the parabola self-heating would give). There is no
injection regime in this data, so it cannot test the forward-bias model-form
question (that was the May neuron2 dataset). What it does support — cleanly —
is a passive + thermal + linear-EO fit, which the linear fc_pn_th_ps nails.

Stages (each vs the extracted observables from the npz cache):
  1. passive : lineshape T(λ) at hc=0, jv≈0    → n_eff, kappa_L, alpha_dB_cm
  2. thermal : notch λ vs heater current       → p_pi_th   (R_heater fixed; only
               P=I²R/p_pi enters, so R and p_pi are degenerate)
  3. EO      : notch λ vs junction voltage     → dn_dv (linear; the device IS linear)

Cache (extracted once from the 3.6 GB pickle):
  .venv/bin/python experiments/giona/fit_jul_neuron3.py --cache  # rebuild
Run:
  .venv/bin/python experiments/giona/fit_jul_neuron3.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize_scalar

import fairchild as fc

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "giona_mod_neuron3_joint_IV_spec_20260710T151201Z.pkl.gz"
CACHE = HERE / "data" / "neuron3_cache.npz"
OUT_JSON = HERE / "results" / "giona_neuron3_pn_th_ps_fit.json"

# Ring identified in the data: notches at 1541.245 / 1553.438 nm respond to
# this neuron's heater (~400 pm) while the other 10 shift only ~30 pm (thermal
# crosstalk). FSR = 12.193 nm → L = λ²/(FSR·n_g).
LAM0 = 1541.245
FSR_NM = 12.193
N_G = 4.2
L_UM = (LAM0 * 1e-9) ** 2 / (FSR_NM * 1e-9 * N_G) * 1e6   # ≈ 46.4 µm
R_HEATER = 184.37     # fixed (neuron7 electrical value); degenerate with p_pi


# ── cache ────────────────────────────────────────────────────────────────────
def build_cache():
    import gzip, pickle
    with gzip.open(DATA, "rb") as f:
        d = pickle.load(f)
    hc = np.asarray(d["Heater Current (mA)"], float)
    jv = np.asarray(d["Junction Voltage (V)"], float)
    ij = np.asarray(d["Junction Current (mA)"], float)
    sp = d["Spectrum"]
    wl = np.asarray(sp[0, 0].absc, float)
    T = np.empty(sp.shape + (len(wl),), np.float32)
    for i in range(sp.shape[0]):
        for j in range(sp.shape[1]):
            T[i, j] = np.asarray(sp[i, j].ordi, np.float32)
    np.savez_compressed(CACHE, hc_mA=hc[:, 0], jv_V=jv[0, :], i_junc_mA=ij,
                        wl_nm=wl, T_dB=T)
    print(f"cached → {CACHE}")


def load_cache():
    d = np.load(CACHE)
    return d["hc_mA"], d["jv_V"], d["i_junc_mA"], d["wl_nm"], d["T_dB"]


def track_notch(wl, y):
    """Sub-pixel notch minimum via parabolic interpolation."""
    k = int(np.argmin(y))
    if 0 < k < len(y) - 1:
        a, b, c = float(y[k - 1]), float(y[k]), float(y[k + 1])
        denom = a - 2 * b + c
        if denom > 0:
            return wl[k] + 0.5 * (a - c) / denom * (wl[1] - wl[0])
    return wl[k]


def measured_tracks(hc, jv, wl, T):
    """(λ vs HC at jv≈0) and (λ vs JV at hc=0), tracked sequentially so the
    ~1.6 nm thermal excursion never jumps to a neighbouring notch."""
    j0 = int(np.argmin(np.abs(jv)))
    lam_hc, ctr = [], LAM0
    for i in range(len(hc)):
        w = (wl > ctr - 0.25) & (wl < ctr + 0.25)
        ctr = track_notch(wl[w], T[i, j0, w])
        lam_hc.append(ctr)
    lam_jv, ctr = [], LAM0
    for j in range(len(jv)):
        w = (wl > ctr - 0.25) & (wl < ctr + 0.25)
        ctr = track_notch(wl[w], T[0, j, w])
        lam_jv.append(ctr)
    return np.array(lam_hc), np.array(lam_jv), j0


# ── simulator forward model ──────────────────────────────────────────────────
NETLIST = f"""* neuron3 all-pass ring (mod-bank member), linear PN + heater
.optical_port in_opt
.optical_port thru_opt
.optical_port ring_in
.optical_port ring_out
Xlaser in_opt fc_cw_laser power_mW=1.0 wavelength_nm={LAM0}
Xdc in_opt ring_in thru_opt ring_out fc_dcoupler kappa_L=0.3
Xps ring_out ring_in vpn 0 heat 0 fc_pn_th_ps L_um={L_UM:.4f} n_g={N_G}
+ alpha_dB_cm=5 dn_dv=0 g_pn=6.3e-5 R_heater={R_HEATER} p_pi_th=0.05 n_eff=2.42
Vpn vpn 0 DC 0
Iheat 0 heat DC 0
.op
"""


class Ring:
    def __init__(self):
        self.ckt = fc.Circuit()
        self.ckt.load_str(NETLIST)

    def set(self, **kw):
        for k, v in kw.items():
            if "." in k:                       # explicit "Element.param"
                el, p = k.split(".", 1)
            elif k in ("kappa_l", "kappa_L"):  # coupler param
                el, p = "Xdc", k
            else:
                el, p = "Xps", k
            self.ckt.set_param(el, p, float(v))

    def spectrum_dB(self, wl_nm):
        res = self.ckt.sweep("Xlaser.wavelength_nm", list(wl_nm), "op")
        p = np.array([float(r["V(thru_opt_re_0)"][0]) ** 2 +
                      float(r["V(thru_opt_im_0)"][0]) ** 2 for r in res])
        return 10 * np.log10(np.maximum(p, 1e-30))   # dB rel. 1 mW in

    def notch(self, lo, hi, n=60):
        wls = np.linspace(lo, hi, n)
        return track_notch(wls, self.spectrum_dB(wls))


# ── stages ───────────────────────────────────────────────────────────────────
def stage1_passive(wl, T, j0, verbose=True):
    w = (wl > LAM0 - 0.45) & (wl < LAM0 + 0.45)
    wl_m, t_m = wl[w][::52], T[0, j0, w][::52].astype(float)   # ~100 pts
    ring = Ring()
    # n_eff fine-positions the notch: one FSR of n_eff at fixed L is λ/L ≈ 0.033.
    dn_period = LAM0 * 1e-9 / (L_UM * 1e-6)

    def cost(x):
        n_eff, kappa_l, alpha = x
        ring.set(n_eff=n_eff, kappa_l=kappa_l, alpha_dB_cm=alpha)
        sim = ring.spectrum_dB(wl_m)
        off = np.median(t_m - sim)                 # insertion loss (projected)
        return float(np.mean((sim + off - t_m) ** 2))

    res = differential_evolution(
        cost, [(2.42, 2.42 + dn_period), (0.05, 1.2), (0.5, 80.0)],
        seed=0, maxiter=40, popsize=14, tol=1e-7, polish=True)
    n_eff, kappa_l, alpha = res.x
    if verbose:
        print(f"stage1: n_eff={n_eff:.6f} kappa_L={kappa_l:.4f} "
              f"alpha={alpha:.2f} dB/cm  rms={np.sqrt(res.fun):.2f} dB")
    return dict(n_eff=n_eff, kappa_l=kappa_l, alpha_db_cm=alpha, rms_dB=np.sqrt(res.fun))


def stage2_thermal(hc, lam_hc, base, verbose=True):
    ring = Ring()
    ring.set(n_eff=base["n_eff"], kappa_l=base["kappa_l"], alpha_dB_cm=base["alpha_db_cm"])
    sub = np.arange(0, len(hc), 4)                  # 13 bias points is plenty

    def cost(p_pi):
        ring.set(p_pi_th=p_pi)
        err = 0.0
        for i in sub:
            ring.set(**{"Iheat.dc": hc[i] * 1e-3})
            lam = ring.notch(lam_hc[i] - 0.25, lam_hc[i] + 0.25)
            err += (lam - lam_hc[i]) ** 2
        return err / len(sub)

    res = minimize_scalar(cost, bounds=(5e-3, 0.5), method="bounded",
                          options=dict(xatol=1e-4))
    if verbose:
        print(f"stage2: p_pi_th={res.x*1e3:.2f} mW  (R_heater fixed {R_HEATER} Ω) "
              f"rms={np.sqrt(res.fun)*1e3:.1f} pm")
    return dict(p_pi_th=res.x, rms_pm=np.sqrt(res.fun) * 1e3)


def stage3_eo(jv, lam_jv, base, verbose=True):
    ring = Ring()
    ring.set(n_eff=base["n_eff"], kappa_l=base["kappa_l"], alpha_dB_cm=base["alpha_db_cm"],
             **{"Iheat.dc": 0.0})
    sub = np.arange(0, len(jv), 9)                  # 12 bias points

    def cost(dn_dv):
        ring.set(dn_dv=dn_dv)
        err = 0.0
        for j in sub:
            ring.set(**{"Vpn.dc": jv[j]})
            lam = ring.notch(lam_jv[j] - 0.15, lam_jv[j] + 0.15)
            err += (lam - lam_jv[j]) ** 2
        return err / len(sub)

    res = minimize_scalar(cost, bounds=(-5e-5, 5e-5), method="bounded",
                          options=dict(xatol=1e-8))
    if verbose:
        print(f"stage3: dn_dv={res.x:.3e} /V  rms={np.sqrt(res.fun)*1e3:.2f} pm")
    return dict(dn_dv=res.x, rms_pm=np.sqrt(res.fun) * 1e3)


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", action="store_true", help="rebuild npz cache from pickle")
    args = ap.parse_args()
    if args.cache or not CACHE.exists():
        build_cache()

    hc, jv, ij, wl, T = load_cache()
    lam_hc, lam_jv, j0 = measured_tracks(hc, jv, wl, T)
    print(f"ring: L={L_UM:.1f} µm (FSR {FSR_NM} nm)  λ0={LAM0} nm")
    print(f"measured: thermal shift {(lam_hc[-1]-lam_hc[0])*1e3:+.0f} pm over "
          f"{hc[-1]:.2f} mA;  EO {(lam_jv[-1]-lam_jv[0])*1e3:+.1f} pm over "
          f"{jv[0]:.0f}..{jv[-1]:.0f} V (linear, no injection knee)\n")

    p1 = stage1_passive(wl, T, j0)
    p2 = stage2_thermal(hc, lam_hc, p1)
    p3 = stage3_eo(jv, lam_jv, p1)

    best = dict(model="fc_pn_th_ps", l_um=L_UM, n_g=N_G, r_heater=R_HEATER,
                **{k: v for k, v in {**p1, **p2, **p3}.items() if not k.startswith("rms")})
    OUT_JSON.write_text(json.dumps(best, indent=2))
    print(f"\nsaved → {OUT_JSON}")

    # validation overlay
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ring = Ring()
    ring.set(n_eff=p1["n_eff"], kappa_l=p1["kappa_l"], alpha_dB_cm=p1["alpha_db_cm"],
             p_pi_th=p2["p_pi_th"], dn_dv=p3["dn_dv"])
    w = (wl > LAM0 - 0.45) & (wl < LAM0 + 0.45)
    wl_m = wl[w][::52]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for i_hc, c in [(0, "tab:blue"), (10, "tab:orange"), (16, "tab:green")]:
        ring.set(**{"Iheat.dc": hc[i_hc] * 1e-3, "Vpn.dc": jv[j0]})
        win = np.linspace(lam_hc[i_hc] - 0.45, lam_hc[i_hc] + 0.45, 100)
        sim = ring.spectrum_dB(win)
        wm = (wl > win[0]) & (wl < win[-1])
        meas = T[i_hc, j0, wm].astype(float)
        off = np.median(meas[::40] - np.interp(wl[wm][::40], win, sim))
        ax[0].plot(wl[wm], meas, color=c, lw=0.6)
        ax[0].plot(win, sim + off, "--", color=c, lw=1.2,
                   label=f"hc={hc[i_hc]:.2f} mA")
    ax[0].set_title("lineshape: meas (solid) vs fit (dashed)")
    ax[0].set_xlabel("λ (nm)"); ax[0].set_ylabel("T (dB)"); ax[0].legend(fontsize=7)
    sim_hc = []
    for i in range(0, len(hc), 2):
        ring.set(**{"Iheat.dc": hc[i] * 1e-3, "Vpn.dc": jv[j0]})
        sim_hc.append(ring.notch(lam_hc[i] - 0.25, lam_hc[i] + 0.25))
    ax[1].plot(hc, (lam_hc - lam_hc[0]) * 1e3, "k.", ms=4, label="measured")
    ax[1].plot(hc[::2], (np.array(sim_hc) - lam_hc[0]) * 1e3, "r-", lw=1, label="model")
    ax[1].set_title("thermal: Δλ vs heater current"); ax[1].set_xlabel("HC (mA)")
    ax[1].set_ylabel("Δλ (pm)"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    ring.set(**{"Iheat.dc": 0.0})
    sim_jv = []
    for j in range(0, len(jv), 4):
        ring.set(**{"Vpn.dc": jv[j]})
        sim_jv.append(ring.notch(lam_jv[j] - 0.15, lam_jv[j] + 0.15))
    ax[2].plot(jv, (lam_jv - np.mean(lam_jv)) * 1e3, "k.", ms=4, label="measured")
    ax[2].plot(jv[::4], (np.array(sim_jv) - np.mean(lam_jv)) * 1e3, "r-", lw=1, label="model")
    ax[2].set_title("EO: Δλ vs junction V (linear — no injection in device)")
    ax[2].set_xlabel("JV (V)"); ax[2].set_ylabel("Δλ (pm)"); ax[2].legend(fontsize=8)
    ax[2].grid(alpha=0.3)
    fig.tight_layout()
    out = HERE / "results" / "neuron3_fit.png"
    fig.savefig(out, dpi=120)
    print(f"plot → {out}")


if __name__ == "__main__":
    main()
