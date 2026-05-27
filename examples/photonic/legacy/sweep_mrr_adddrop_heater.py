#!/usr/bin/env python3
"""
MRR heater-tuned resonator (add-drop) — wavelength-domain characterization.

Sweeps wavelength at several heater voltages for the L1 and L2 add-drop
heater-integrated ring resonator models.  Each subplot shows through (solid)
and drop (dashed) transmission vs wavelength.  The add port is left dark.

L1: integrated thermal model (R_thermal = 30 kK/W, R_heater = 500 Ω).
L2: external T_node thermal network (Rth connected to ground in netlist).

Requirements: fairchild Python package (maturin develop), numpy, matplotlib.
Compiled OSDI models in legacy/va-models/build/.
"""

import math
import pathlib
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

try:
    import fairchild as fc
except ImportError:
    sys.exit("fairchild not installed — run 'maturin develop' from the repo root first.")

HERE  = pathlib.Path(__file__).resolve().parent
BUILD = HERE.parents[2] / "legacy" / "va-models" / "build"

# ── Parameters ─────────────────────────────────────────────────────────────────

KAPPA_1    = 0.1
KAPPA_2    = 0.1
L_RING_UM  = 100.0
N_G        = 4.2
ALPHA_DCM  = 10.0      # dB/cm — N-doped waveguide
R_HEATER   = 500.0
R_THERMAL  = 30000.0
DN_DT      = 1.86e-4
WL_NM      = 1550.0

_Vpi2 = math.pi / (
    2 * math.pi * DN_DT * (R_THERMAL / R_HEATER) * (L_RING_UM * 1e-6) / (WL_NM * 1e-9)
)
VPI_HEATER = math.sqrt(_Vpi2)   # ≈ 0.83 V

V_HEATERS = [0.0, 0.3, 0.5, 0.7, 0.9]
WL_START  = 1542.0
WL_END    = 1558.0
N_WL      = 501

# ── CMT reference ──────────────────────────────────────────────────────────────

def cmt_adddrop(wl_nm, k1, k2, L_um, n_g, alpha_dB_cm, dphi=0.0):
    """Add-drop ring CMT (mirrors VA model exactly). Returns (T_through, T_drop)."""
    r1 = math.sqrt(1 - k1)
    r2 = math.sqrt(1 - k2)
    t1 = math.sqrt(k1)
    t2 = math.sqrt(k2)
    alpha_Np = alpha_dB_cm * 100 / 8.685889
    a_half   = math.exp(-alpha_Np * L_um * 1e-6 / 2.0)
    phi_rt   = 2 * math.pi * n_g * L_um * 1e-6 / (wl_nm * 1e-9) + dphi
    phi_half = phi_rt / 2.0

    c_half = a_half * math.cos(phi_half)
    s_half = a_half * math.sin(phi_half)
    c_rt   = c_half * c_half - s_half * s_half
    s_rt   = 2.0 * c_half * s_half

    dA  = 1.0 - r1 * r2 * c_rt
    dB  = r1 * r2 * s_rt
    dsq = dA * dA + dB * dB

    N11_re = r1 - r2 * c_rt;    N11_im = -r2 * s_rt
    N12_re = -t1 * t2 * c_half; N12_im = -t1 * t2 * s_half

    T11_re = (N11_re * dA - N11_im * dB) / dsq
    T11_im = (N11_im * dA + N11_re * dB) / dsq
    T12_re = (N12_re * dA - N12_im * dB) / dsq
    T12_im = (N12_im * dA + N12_re * dB) / dsq

    return T11_re**2 + T11_im**2, T12_re**2 + T12_im**2

def heater_dphi(V_heat, r_heater, r_thermal, dn_dt, L_um, wl_nm):
    P  = V_heat * V_heat / r_heater
    dT = P * r_thermal
    return 2 * math.pi * dn_dt * dT * L_um * 1e-6 / (wl_nm * 1e-9)

# ── Netlist builders ───────────────────────────────────────────────────────────

def _osdi(*names):
    return "\n".join(f".osdi {BUILD / n}.osdi" for n in names)

def netlist_l1():
    return f"""\
* MRR heater-tuned add-drop L1 (integrated thermal)
{_osdi("cw_laser", "mrr_heater_l1_adddrop", "photodetector")}
Xlaser  l_re l_im l_wl  cw_laser  power_mW=1.0 wavelength_nm={WL_NM}
Xring   l_re l_im l_wl  th_re th_im l_wl  dp_re dp_im l_wl  ad_re ad_im l_wl \\
        hp hn  mrr_heater_l1_adddrop \\
        kappa_1={KAPPA_1} kappa_2={KAPPA_2} L_ring_um={L_RING_UM} n_g={N_G} \\
        alpha_dB_cm={ALPHA_DCM} R_heater={R_HEATER} R_thermal={int(R_THERMAL)} \\
        dn_dT={DN_DT} wavelength_nm={WL_NM}
Xpd_th  th_re th_im l_wl  ph_th 0  photodetector  responsivity=1.0
Xpd_dp  dp_re dp_im l_wl  ph_dp 0  photodetector  responsivity=1.0
Rad_re  ad_re 0  1G
Rad_im  ad_im 0  1G
Rth     ph_th 0  1k
Rdp     ph_dp 0  1k
Vheat   hp hn  DC 0.0
.optical l_re l_im l_wl th_re th_im dp_re dp_im ad_re ad_im
.op
.end
"""

def netlist_l2():
    return f"""\
* MRR heater-tuned add-drop L2 (external T_node thermal network)
{_osdi("cw_laser", "mrr_heater_l2_adddrop", "photodetector")}
Xlaser  l_re l_im l_wl  cw_laser  power_mW=1.0 wavelength_nm={WL_NM}
Xring   l_re l_im l_wl  th_re th_im l_wl  dp_re dp_im l_wl  ad_re ad_im l_wl \\
        hp hn  tnode  mrr_heater_l2_adddrop \\
        kappa_1={KAPPA_1} kappa_2={KAPPA_2} L_ring_um={L_RING_UM} n_g={N_G} \\
        alpha_dB_cm={ALPHA_DCM} R_heater={R_HEATER} dn_dT={DN_DT} wavelength_nm={WL_NM}
Rth_th  tnode 0  {int(R_THERMAL)}
Xpd_th  th_re th_im l_wl  ph_th 0  photodetector  responsivity=1.0
Xpd_dp  dp_re dp_im l_wl  ph_dp 0  photodetector  responsivity=1.0
Rad_re  ad_re 0  1G
Rad_im  ad_im 0  1G
Rth     ph_th 0  1k
Rdp     ph_dp 0  1k
Vheat   hp hn  DC 0.0
.optical l_re l_im l_wl th_re th_im dp_re dp_im ad_re ad_im
.op
.end
"""

# ── Sweep ──────────────────────────────────────────────────────────────────────

NORM = 1.0e-3 * 1.0 * 1.0e3   # P_in*resp*R_load = 1.0 V at T=1

def run_sweep(netlist_fn, label):
    wl_list = list(np.linspace(WL_START, WL_END, N_WL))
    ckt = fc.Circuit()
    ckt.load_str(netlist_fn())
    data = {}
    for v in V_HEATERS:
        print(f"  {label}  V_heat={v:.2f} V ...", end="", flush=True)
        ckt.set_param("Vheat", "dc", v)
        try:
            results = ckt.sweep("Xlaser.wavelength_nm", wl_list, "op")
            T_th = [r["V(ph_th)"][0] / NORM for r in results]
            T_dp = [r["V(ph_dp)"][0] / NORM for r in results]
        except Exception as exc:
            print(f"  ERROR: {exc}")
            T_th = T_dp = [float("nan")] * len(wl_list)
        data[v] = (T_th, T_dp)
        print(" done")
    return data

# ── Plot ───────────────────────────────────────────────────────────────────────

def plot(wl, data_l1, data_l2):
    colors = cm.Oranges_r(np.linspace(0.15, 0.85, len(V_HEATERS)))
    dphi_list = [heater_dphi(v, R_HEATER, R_THERMAL, DN_DT, L_RING_UM, WL_NM)
                 for v in V_HEATERS]

    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    titles = [
        "L1 — integrated thermal (R_thermal = 30 kK/W, R_heater = 500 Ω)",
        "L2 — external T_node (Rth = 30 kΩ to ground)",
    ]

    for ax, data, title in zip(axes, [data_l1, data_l2], titles):
        for i, (v, color) in enumerate(zip(V_HEATERS, colors)):
            dphi = dphi_list[i]
            T_th, T_dp = data[v]
            lbl = f"{v:.2f} V  (Δφ={dphi/math.pi:.2f}π)"
            ax.plot(wl, T_th, color=color, lw=1.6, label=lbl)
            ax.plot(wl, T_dp, color=color, lw=1.6, ls="--")

            # Analytical CMT reference (thru=solid, drop=dashed, same color, thin)
            T_th_cmt = [cmt_adddrop(w, KAPPA_1, KAPPA_2, L_RING_UM, N_G,
                                    ALPHA_DCM, dphi)[0] for w in wl]
            T_dp_cmt = [cmt_adddrop(w, KAPPA_1, KAPPA_2, L_RING_UM, N_G,
                                    ALPHA_DCM, dphi)[1] for w in wl]
            ax.plot(wl, T_th_cmt, color=color, lw=0.8, ls=":", alpha=0.6)
            ax.plot(wl, T_dp_cmt, color=color, lw=0.8, ls="-.", alpha=0.6)

        ax.set_ylabel("Transmission")
        ax.set_title(title, fontsize=10)
        ax.set_ylim(-0.05, 1.10)
        ax.legend(fontsize=8, loc="center right", ncol=2)
        ax.grid(True, alpha=0.3)
        ax.text(0.02, 0.05, "Thick = simulation   Thin = CMT analytic\n"
                "Solid = through   Dashed = drop",
                transform=ax.transAxes, fontsize=8, color="gray")

    axes[-1].set_xlabel("Wavelength (nm)")
    param_str = (f"κ₁=κ₂={KAPPA_1}  L={L_RING_UM:.0f} µm  n_g={N_G}  "
                 f"α={ALPHA_DCM} dB/cm  Vπ≈{VPI_HEATER:.2f} V")
    fig.suptitle(f"MRR Heater Resonator (add-drop) — Wavelength Sweep\n{param_str}",
                 fontsize=11)
    fig.tight_layout()

    out = HERE / "sweep_mrr_adddrop_heater.png"
    fig.savefig(out, dpi=150)
    print(f"\nFigure saved: {out}")
    plt.close(fig)

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    for name in ("cw_laser", "mrr_heater_l1_adddrop", "mrr_heater_l2_adddrop",
                 "photodetector"):
        if not (BUILD / f"{name}.osdi").exists():
            sys.exit(f"Missing: {BUILD / name}.osdi — compile legacy/va-models first.")

    print(f"Heater Vπ ≈ {VPI_HEATER:.3f} V")
    wl = np.linspace(WL_START, WL_END, N_WL)

    print("Sweeping L1 (integrated thermal, add-drop)...")
    d1 = run_sweep(netlist_l1, "L1")
    print("Sweeping L2 (external T_node, add-drop)...")
    d2 = run_sweep(netlist_l2, "L2")

    plot(wl, d1, d2)


if __name__ == "__main__":
    main()
