#!/usr/bin/env python3
"""
MRR heater-tuned resonator (all-pass) — wavelength-domain characterization.

Sweeps wavelength at several heater voltages for the L1 and L2 all-pass
heater-integrated ring resonator models.

L1: integrated thermal model (R_thermal = 30 kK/W, R_heater = 500 Ω).
    Vπ ≈ 0.83 V (heater voltage for a π round-trip phase shift at 1550 nm).
L2: external T_node thermal network (Rth connected to ground in netlist).
    Same steady-state physics as L1; T_node allows RC thermal dynamics.

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

# ── Ring / heater parameters ───────────────────────────────────────────────────

KAPPA      = 0.1
L_RING_UM  = 100.0
N_G        = 4.2
ALPHA_DCM  = 10.0     # dB/cm — N-doped waveguide (higher than undoped due to FCA)
R_HEATER   = 500.0    # Ω
R_THERMAL  = 30000.0  # K/W
DN_DT      = 1.86e-4  # K⁻¹ (Si thermo-optic coefficient at 1550 nm)
WL_NM      = 1550.0   # design wavelength (nm)

# Vπ = sqrt(π / (2π × dn_dT × (R_thermal/R_heater) × L/λ))
_Vpi2 = (math.pi / (
    2 * math.pi * DN_DT * (R_THERMAL / R_HEATER) * (L_RING_UM * 1e-6) / (WL_NM * 1e-9)
))
VPI_HEATER = math.sqrt(_Vpi2)   # ≈ 0.83 V

# Heater voltages: 0 → ~1.5 × Vπ
V_HEATERS  = [0.0, 0.3, 0.5, 0.7, 0.9]   # V
WL_START   = 1542.0
WL_END     = 1558.0
N_WL       = 501

# ── CMT reference ──────────────────────────────────────────────────────────────

def cmt_allpass(wl_nm, kappa, L_um, n_g, alpha_dB_cm, dphi=0.0):
    """All-pass ring CMT power transmission (matches VA model formula exactly)."""
    r = math.sqrt(1 - kappa)
    alpha_Np = alpha_dB_cm * 100 / 8.685889
    a_rt = math.exp(-alpha_Np * L_um * 1e-6)
    phi = 2 * math.pi * n_g * L_um * 1e-6 / (wl_nm * 1e-9) + dphi
    c = a_rt * math.cos(phi)
    s = a_rt * math.sin(phi)
    kap = 1 - r * r
    dsq = 1 - 2 * r * c + r * r * a_rt * a_rt
    T_re = ((r - c) * (1 - r * c) + r * s * s) / dsq
    T_im = -s * kap / dsq
    return T_re * T_re + T_im * T_im

def heater_dphi(V_heat, r_heater, r_thermal, dn_dt, L_um, wl_nm):
    """Round-trip thermal phase shift for heater voltage V."""
    P = V_heat * V_heat / r_heater
    dT = P * r_thermal
    return 2 * math.pi * dn_dt * dT * L_um * 1e-6 / (wl_nm * 1e-9)

# ── Netlist builders ───────────────────────────────────────────────────────────

def _osdi(*names):
    return "\n".join(f".osdi {BUILD / n}.osdi" for n in names)

def netlist_l1():
    return f"""\
* MRR heater-tuned L1 (all-pass, integrated thermal)
{_osdi("cw_laser", "mrr_heater_l1", "photodetector")}
Xlaser  l_re l_im l_wl  cw_laser  power_mW=1.0 wavelength_nm={WL_NM}
Xring   l_re l_im l_wl  o_re o_im l_wl  hp hn  mrr_heater_l1 \\
        kappa_0={KAPPA} L_ring_um={L_RING_UM} n_g={N_G} alpha_dB_cm={ALPHA_DCM} \\
        R_heater={R_HEATER} R_thermal={R_THERMAL} dn_dT={DN_DT} wavelength_nm={WL_NM}
Xpd     o_re o_im l_wl  ph_out 0  photodetector  responsivity=1.0
Rload   ph_out 0  1k
Vheat   hp hn  DC 0.0
.optical l_re l_im l_wl o_re o_im
.op
.end
"""

def netlist_l2():
    return f"""\
* MRR heater-tuned L2 (all-pass, external T_node thermal network)
{_osdi("cw_laser", "mrr_heater_l2", "photodetector")}
Xlaser  l_re l_im l_wl  cw_laser  power_mW=1.0 wavelength_nm={WL_NM}
Xring   l_re l_im l_wl  o_re o_im l_wl  hp hn  tnode  mrr_heater_l2 \\
        kappa_0={KAPPA} L_ring_um={L_RING_UM} n_g={N_G} alpha_dB_cm={ALPHA_DCM} \\
        R_heater={R_HEATER} dn_dT={DN_DT} wavelength_nm={WL_NM}
Rth     tnode 0  {int(R_THERMAL)}
Xpd     o_re o_im l_wl  ph_out 0  photodetector  responsivity=1.0
Rload   ph_out 0  1k
Vheat   hp hn  DC 0.0
.optical l_re l_im l_wl o_re o_im
.op
.end
"""

# ── Sweep ──────────────────────────────────────────────────────────────────────

P_IN_MW = 1.0
NORM    = P_IN_MW * 1e-3 * 1.0 * 1e3   # 1.0 V at full transmission

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
            T = [r["V(ph_out)"][0] / NORM for r in results]
        except Exception as exc:
            print(f"  ERROR: {exc}")
            T = [float("nan")] * len(wl_list)
        data[v] = T
        print(" done")
    return data

# ── Plot ───────────────────────────────────────────────────────────────────────

def plot(wl, data_l1, data_l2):
    colors = cm.Reds_r(np.linspace(0.15, 0.85, len(V_HEATERS)))
    labels = [f"{v:.2f} V" for v in V_HEATERS]
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
            T_cmt = [cmt_allpass(w, KAPPA, L_RING_UM, N_G, ALPHA_DCM, dphi) for w in wl]
            T_sim = data[v]
            ax.plot(wl, T_sim, color=color, lw=1.8, label=f"{v:.2f} V  (Δφ={dphi/math.pi:.2f}π)")
            ax.plot(wl, T_cmt, color=color, lw=0.9, ls="--", alpha=0.6)

        ax.set_ylabel("Transmission")
        ax.set_title(title, fontsize=10)
        ax.set_ylim(-0.05, 1.10)
        ax.legend(fontsize=8, loc="lower right", ncol=3)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Wavelength (nm)")
    param_str = (f"κ={KAPPA}  L={L_RING_UM:.0f} µm  n_g={N_G}  "
                 f"α={ALPHA_DCM} dB/cm  Vπ≈{VPI_HEATER:.2f} V")
    fig.suptitle(f"MRR Heater Resonator (all-pass) — Wavelength Sweep\n{param_str}", fontsize=11)
    ax.text(0.02, 0.05, "Solid = simulation   Dashed = CMT analytic",
            transform=ax.transAxes, fontsize=8, color="gray")
    fig.tight_layout()

    out = HERE / "sweep_mrr_allpass_heater.png"
    fig.savefig(out, dpi=150)
    print(f"\nFigure saved: {out}")
    plt.close(fig)

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    for name in ("cw_laser", "mrr_heater_l1", "mrr_heater_l2", "photodetector"):
        if not (BUILD / f"{name}.osdi").exists():
            sys.exit(f"Missing OSDI model: {BUILD / name}.osdi\n"
                     "Compile va-models first: cd legacy/va-models && ./build.sh")

    print(f"Heater Vπ ≈ {VPI_HEATER:.3f} V")
    wl = np.linspace(WL_START, WL_END, N_WL)

    print("Sweeping L1 (integrated thermal)...")
    d1 = run_sweep(netlist_l1, "L1")
    print("Sweeping L2 (external T_node)...")
    d2 = run_sweep(netlist_l2, "L2")

    plot(wl, d1, d2)


if __name__ == "__main__":
    main()
