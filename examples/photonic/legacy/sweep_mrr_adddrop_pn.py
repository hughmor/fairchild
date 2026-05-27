#!/usr/bin/env python3
"""
MRR modulator (add-drop, PN junction) — wavelength-domain characterization.

Sweeps wavelength at several reverse-bias voltages for the L1, L2, and L3
add-drop ring resonator models.  Each subplot shows the through port (solid)
and drop port (dashed) transmission vs wavelength, with curves for voltages
[0, -3, -6, -9] V.  The add port is left dark (no input signal).

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
ALPHA_DCM  = 2.0
VPI_RT     = 10.0
N_DEP0     = 5e16
VBI        = 1.0
WL_NM      = 1550.0

V_BIASES   = [0.0, -3.0, -6.0, -9.0]
WL_START   = 1542.0
WL_END     = 1558.0
N_WL       = 501

# ── CMT reference ──────────────────────────────────────────────────────────────

def cmt_adddrop(wl_nm, k1, k2, L_um, n_g, alpha_dB_cm, dphi=0.0):
    """
    Add-drop ring CMT — mirrors the VA model formula exactly.
    Returns (T_through, T_drop) with add input = 0.
    """
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

    dA   = 1.0 - r1 * r2 * c_rt
    dB   = r1 * r2 * s_rt
    dsq  = dA * dA + dB * dB

    N11_re = r1 - r2 * c_rt;   N11_im = -r2 * s_rt
    N12_re = -t1 * t2 * c_half; N12_im = -t1 * t2 * s_half

    T11_re = (N11_re * dA - N11_im * dB) / dsq
    T11_im = (N11_im * dA + N11_re * dB) / dsq
    T12_re = (N12_re * dA - N12_im * dB) / dsq
    T12_im = (N12_im * dA + N12_re * dB) / dsq

    T_th = T11_re * T11_re + T11_im * T11_im
    T_dp = T12_re * T12_re + T12_im * T12_im
    return T_th, T_dp

# ── Netlist builders ───────────────────────────────────────────────────────────

def _osdi(*names):
    return "\n".join(f".osdi {BUILD / n}.osdi" for n in names)

def _adddrop_netlist(model_name, extra_params="", extra_elements=""):
    return f"""\
* MRR add-drop ({model_name})
{_osdi("cw_laser", model_name, "photodetector")}
Xlaser  l_re l_im l_wl  cw_laser  power_mW=1.0 wavelength_nm={WL_NM}
Xring   l_re l_im l_wl  th_re th_im l_wl  dp_re dp_im l_wl  ad_re ad_im l_wl \\
        vbias 0  {model_name} \\
        kappa_1={KAPPA_1} kappa_2={KAPPA_2} L_ring_um={L_RING_UM} n_g={N_G} \\
        alpha_dB_cm={ALPHA_DCM} Vpi_rt={VPI_RT} V_ref=0.0 wavelength_nm={WL_NM}{extra_params}
{extra_elements}Xpd_th  th_re th_im l_wl  ph_th 0  photodetector  responsivity=1.0
Xpd_dp  dp_re dp_im l_wl  ph_dp 0  photodetector  responsivity=1.0
Rad_re  ad_re 0  1G
Rad_im  ad_im 0  1G
Rth     ph_th 0  1k
Rdp     ph_dp 0  1k
Vbias   vbias 0  DC 0.0
.optical l_re l_im l_wl th_re th_im dp_re dp_im ad_re ad_im
.op
.end
"""

def netlist_l1():
    return _adddrop_netlist("mrr_modulator_l1_adddrop")

def netlist_l2():
    return _adddrop_netlist(
        "mrr_modulator_l2_adddrop",
        extra_params=f" \\\n        n_dep0={N_DEP0} Vbi={VBI}",
    )

def netlist_l3():
    return f"""\
* MRR add-drop L3 (PN + TPA + T_node)
{_osdi("cw_laser", "mrr_modulator_l3_adddrop", "photodetector")}
Xlaser  l_re l_im l_wl  cw_laser  power_mW=1.0 wavelength_nm={WL_NM}
Xring   l_re l_im l_wl  th_re th_im l_wl  dp_re dp_im l_wl  ad_re ad_im l_wl \\
        vbias 0  tnode  mrr_modulator_l3_adddrop \\
        kappa_1={KAPPA_1} kappa_2={KAPPA_2} L_ring_um={L_RING_UM} n_g={N_G} \\
        alpha_dB_cm={ALPHA_DCM} Vpi_rt={VPI_RT} V_ref=0.0 wavelength_nm={WL_NM} \\
        n_dep0={N_DEP0} Vbi={VBI}
Rtpa    tnode 0  30000
Xpd_th  th_re th_im l_wl  ph_th 0  photodetector  responsivity=1.0
Xpd_dp  dp_re dp_im l_wl  ph_dp 0  photodetector  responsivity=1.0
Rad_re  ad_re 0  1G
Rad_im  ad_im 0  1G
Rth     ph_th 0  1k
Rdp     ph_dp 0  1k
Vbias   vbias 0  DC 0.0
.optical l_re l_im l_wl th_re th_im dp_re dp_im ad_re ad_im
.op
.end
"""

# ── Sweep ──────────────────────────────────────────────────────────────────────

P_IN_MW = 1.0
NORM    = P_IN_MW * 1e-3 * 1.0 * 1e3

def run_sweep(netlist_fn, label):
    wl_list = list(np.linspace(WL_START, WL_END, N_WL))
    ckt = fc.Circuit()
    ckt.load_str(netlist_fn())
    data = {}
    for v in V_BIASES:
        print(f"  {label}  Vbias={v:+.1f} V ...", end="", flush=True)
        ckt.set_param("Vbias", "dc", v)
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

def plot(wl, data_l1, data_l2, data_l3):
    colors = cm.Purples_r(np.linspace(0.2, 0.85, len(V_BIASES)))
    labels = [f"{v:+.0f} V" for v in V_BIASES]

    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
    titles = [
        "L1 — linear PN phase (Vpi_rt = 10 V)",
        "L2 — Soref–Bennett plasma dispersion",
        "L3 — L2 + TPA self-heating",
    ]

    for ax, data, title in zip(axes, [data_l1, data_l2, data_l3], titles):
        for i, (v, color) in enumerate(zip(V_BIASES, colors)):
            T_th, T_dp = data[v]
            ax.plot(wl, T_th, color=color, lw=1.6, label=f"Through {labels[i]}")
            ax.plot(wl, T_dp, color=color, lw=1.6, ls="--")

        # Analytical CMT overlay for L1 (through=solid, drop=dashed, grey)
        if data is data_l1:
            cmt_th = [cmt_adddrop(w, KAPPA_1, KAPPA_2, L_RING_UM, N_G, ALPHA_DCM)[0] for w in wl]
            cmt_dp = [cmt_adddrop(w, KAPPA_1, KAPPA_2, L_RING_UM, N_G, ALPHA_DCM)[1] for w in wl]
            ax.plot(wl, cmt_th, "k-",  lw=0.9, alpha=0.4, label="CMT thru (V=0)")
            ax.plot(wl, cmt_dp, "k--", lw=0.9, alpha=0.4, label="CMT drop (V=0)")

        ax.set_ylabel("Transmission")
        ax.set_title(title, fontsize=10)
        ax.set_ylim(-0.05, 1.10)
        ax.legend(fontsize=7, loc="center right", ncol=2)
        ax.grid(True, alpha=0.3)
        ax.text(0.02, 0.05, "Solid = through port   Dashed = drop port",
                transform=ax.transAxes, fontsize=8, color="gray")

    axes[-1].set_xlabel("Wavelength (nm)")
    param_str = (f"κ₁=κ₂={KAPPA_1}  L={L_RING_UM:.0f} µm  n_g={N_G}  "
                 f"α={ALPHA_DCM} dB/cm  Vpi_rt={VPI_RT} V")
    fig.suptitle(f"MRR Modulator (add-drop, PN) — Wavelength Sweep\n{param_str}", fontsize=11)
    fig.tight_layout()

    out = HERE / "sweep_mrr_adddrop_pn.png"
    fig.savefig(out, dpi=150)
    print(f"\nFigure saved: {out}")
    plt.close(fig)

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    required = ("cw_laser", "mrr_modulator_l1_adddrop", "mrr_modulator_l2_adddrop",
                "mrr_modulator_l3_adddrop", "photodetector")
    for name in required:
        if not (BUILD / f"{name}.osdi").exists():
            sys.exit(f"Missing: {BUILD / name}.osdi — compile legacy/va-models first.")

    wl = np.linspace(WL_START, WL_END, N_WL)

    print("Sweeping L1 (linear PN, add-drop)...")
    d1 = run_sweep(netlist_l1, "L1")
    print("Sweeping L2 (Soref–Bennett, add-drop)...")
    d2 = run_sweep(netlist_l2, "L2")
    print("Sweeping L3 (TPA + T_node, add-drop)...")
    d3 = run_sweep(netlist_l3, "L3")

    plot(wl, d1, d2, d3)


if __name__ == "__main__":
    main()
