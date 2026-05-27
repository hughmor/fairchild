#!/usr/bin/env python3
"""
MRR modulator (all-pass, PN junction) — wavelength-domain characterization.

Sweeps wavelength at several reverse-bias voltages for the L1, L2, and L3
all-pass ring resonator models.  Each subplot shows transmission vs wavelength
for voltage curves [0, -2, -4, -6, -8] V.  An analytical CMT reference (dashed
grey) is overlaid on the L1 subplot.

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

# ── Ring parameters ────────────────────────────────────────────────────────────

KAPPA      = 0.1
L_RING_UM  = 100.0
N_G        = 4.2
ALPHA_L1   = 2.0      # dB/cm, passive ring
N_DEP0     = 5e16     # cm⁻³, for L2/L3 Soref–Bennett
VBI        = 1.0      # V, built-in potential
VPI_RT     = 10.0     # V, for L1 linear model
WL_NM      = 1550.0   # design wavelength (nm)

V_BIASES   = [0.0, -2.0, -4.0, -6.0, -8.0]   # V (reverse bias)
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

# ── Netlist builders ───────────────────────────────────────────────────────────

def _osdi(*names):
    lines = [f".osdi {BUILD / n}.osdi" for n in names]
    return "\n".join(lines)

def netlist_l1():
    return f"""\
* MRR modulator L1 (all-pass, linear PN)
{_osdi("cw_laser", "mrr_modulator_l1", "photodetector")}
Xlaser  l_re l_im l_wl  cw_laser  power_mW=1.0 wavelength_nm={WL_NM}
Xring   l_re l_im l_wl  o_re o_im l_wl  vbias 0  mrr_modulator_l1 \\
        kappa_0={KAPPA} L_ring_um={L_RING_UM} n_g={N_G} alpha_dB_cm={ALPHA_L1} \\
        Vpi_rt={VPI_RT} V_ref=0.0 wavelength_nm={WL_NM}
Xpd     o_re o_im l_wl  ph_out 0  photodetector  responsivity=1.0
Rload   ph_out 0  1k
Vbias   vbias 0  DC 0.0
.optical l_re l_im l_wl o_re o_im
.op
.end
"""

def netlist_l2():
    return f"""\
* MRR modulator L2 (all-pass, Soref–Bennett)
{_osdi("cw_laser", "mrr_modulator_l2", "photodetector")}
Xlaser  l_re l_im l_wl  cw_laser  power_mW=1.0 wavelength_nm={WL_NM}
Xring   l_re l_im l_wl  o_re o_im l_wl  vbias 0  mrr_modulator_l2 \\
        kappa_0={KAPPA} L_ring_um={L_RING_UM} n_g={N_G} alpha_dB_cm={ALPHA_L1} \\
        n_dep0={N_DEP0} Vbi={VBI} wavelength_nm={WL_NM}
Xpd     o_re o_im l_wl  ph_out 0  photodetector  responsivity=1.0
Rload   ph_out 0  1k
Vbias   vbias 0  DC 0.0
.optical l_re l_im l_wl o_re o_im
.op
.end
"""

def netlist_l3():
    return f"""\
* MRR modulator L3 (all-pass, TPA + T_node)
{_osdi("cw_laser", "mrr_modulator_l3", "photodetector")}
Xlaser  l_re l_im l_wl  cw_laser  power_mW=1.0 wavelength_nm={WL_NM}
Xring   l_re l_im l_wl  o_re o_im l_wl  vbias 0  tnode  mrr_modulator_l3 \\
        kappa_0={KAPPA} L_ring_um={L_RING_UM} n_g={N_G} alpha_dB_cm={ALPHA_L1} \\
        n_dep0={N_DEP0} Vbi={VBI} wavelength_nm={WL_NM}
Rth     tnode 0  30k
Xpd     o_re o_im l_wl  ph_out 0  photodetector  responsivity=1.0
Rload   ph_out 0  1k
Vbias   vbias 0  DC 0.0
.optical l_re l_im l_wl o_re o_im
.op
.end
"""

# ── Sweep ──────────────────────────────────────────────────────────────────────

P_IN_MW  = 1.0
RESP     = 1.0     # A/W
R_LOAD   = 1e3     # Ω
NORM     = P_IN_MW * 1e-3 * RESP * R_LOAD   # = 1.0 V at full transmission

def run_sweep(netlist_fn, label):
    """Sweep wavelength at each bias; return dict {v_bias: [T ...]}."""
    wl_list = list(np.linspace(WL_START, WL_END, N_WL))
    ckt = fc.Circuit()
    ckt.load_str(netlist_fn())

    data = {}
    for v in V_BIASES:
        print(f"  {label}  Vbias={v:+.1f} V ...", end="", flush=True)
        ckt.set_param("Vbias", "dc", v)
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

def plot(wl, data_l1, data_l2, data_l3):
    colors = cm.Blues_r(np.linspace(0.2, 0.85, len(V_BIASES)))
    labels = [f"{v:+.0f} V" for v in V_BIASES]

    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
    titles = [
        "L1 — linear PN phase (Vpi_rt = 10 V)",
        "L2 — Soref–Bennett plasma dispersion",
        "L3 — L2 + TPA self-heating",
    ]

    for ax, data, title in zip(axes, [data_l1, data_l2, data_l3], titles):
        for i, (v, color) in enumerate(zip(V_BIASES, colors)):
            T = data[v]
            ax.plot(wl, T, color=color, lw=1.6, label=labels[i])

        # Analytical CMT reference at V=0 (L1 only)
        if data is data_l1:
            T_cmt = [cmt_allpass(w, KAPPA, L_RING_UM, N_G, ALPHA_L1) for w in wl]
            ax.plot(wl, T_cmt, "k--", lw=1.0, alpha=0.5, label="CMT (V=0)")

        ax.set_ylabel("Transmission")
        ax.set_title(title, fontsize=10)
        ax.set_ylim(-0.05, 1.10)
        ax.legend(fontsize=8, loc="lower right", ncol=3)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Wavelength (nm)")
    param_str = (f"κ={KAPPA}  L={L_RING_UM:.0f} µm  n_g={N_G}  "
                 f"α={ALPHA_L1} dB/cm  P_in={P_IN_MW} mW")
    fig.suptitle(f"MRR Modulator (all-pass, PN) — Wavelength Sweep\n{param_str}", fontsize=11)
    fig.tight_layout()

    out = HERE / "sweep_mrr_allpass_pn.png"
    fig.savefig(out, dpi=150)
    print(f"\nFigure saved: {out}")
    plt.close(fig)

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    for name in ("cw_laser", "mrr_modulator_l1", "mrr_modulator_l2",
                 "mrr_modulator_l3", "photodetector"):
        if not (BUILD / f"{name}.osdi").exists():
            sys.exit(f"Missing OSDI model: {BUILD / name}.osdi\n"
                     "Compile va-models first: cd legacy/va-models && ./build.sh")

    wl = np.linspace(WL_START, WL_END, N_WL)

    print("Sweeping L1 (linear PN)...")
    d1 = run_sweep(netlist_l1, "L1")
    print("Sweeping L2 (Soref–Bennett)...")
    d2 = run_sweep(netlist_l2, "L2")
    print("Sweeping L3 (TPA + T_node)...")
    d3 = run_sweep(netlist_l3, "L3")

    plot(wl, d1, d2, d3)


if __name__ == "__main__":
    main()
