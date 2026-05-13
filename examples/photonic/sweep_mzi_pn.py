#!/usr/bin/env python3
"""
MZI electro-optic modulator (PN junction) — voltage-domain characterization.

Sweeps reverse-bias voltage at λ = 1550 nm for the L1 (single-arm) and L2
(push-pull dual-arm) PN MZI models.  Each subplot shows bar port (solid) and
cross port (dashed) transmission vs bias voltage, with an analytical overlay.

L1: single-arm drive.  Vπ = Vpi_L / L_cm = 2.0 / 0.05 = 40 V.
    At Δφ=0 (V=0): cross full, bar dark.
    At Δφ=π (V=−40 V): bar full, cross dark.

L2: push-pull dual-arm (anode1/cathode1 and anode2/cathode2).
    Topology: V_pn1 = +V, V_pn2 = −V → differential Δφ doubles.
    Effective Vπ ≈ 20 V for Soref–Bennett model at n_dep0 = 5e16 cm⁻³.

Requirements: fairchild Python package (maturin develop), numpy, matplotlib.
Compiled OSDI models in va-models/build/.
"""

import math
import pathlib
import sys

import numpy as np
import matplotlib.pyplot as plt

try:
    import fairchild as fc
except ImportError:
    sys.exit("fairchild not installed — run 'maturin develop' from the repo root first.")

HERE  = pathlib.Path(__file__).resolve().parent
BUILD = HERE.parents[1] / "va-models" / "build"

# ── MZI parameters ─────────────────────────────────────────────────────────────

L_ARM_UM   = 500.0      # µm
N_G        = 4.2
ALPHA_DCM  = 3.0        # dB/cm
VPI_L      = 2.0        # V·cm  (L1 linear model)
V_REF      = 0.0
WL_NM      = 1550.0

L_CM       = L_ARM_UM * 1e-4        # cm
VPI_L1     = VPI_L / L_CM           # ≈ 40 V (single-arm)
VPI_L2_PP  = VPI_L1 / 2             # ≈ 20 V (push-pull ideal)

# Voltage sweep: 0 to −2Vπ (one full sinusoidal period)
V_L1 = np.linspace(0, -2 * VPI_L1, 81)      # 0 to −80 V
V_L2 = np.linspace(0, -2 * VPI_L2_PP, 81)   # 0 to −40 V (approximate)

# ── Analytical MZI transfer function ──────────────────────────────────────────

def mzi_l1_transfer(V_bias, L_um, n_g, alpha_dB_cm, vpi_L, v_ref, wl_nm):
    """
    Analytical MZI L1 bar/cross power transmission (matches VA model).
    Returns (T_bar, T_cross).
    """
    alpha_Np = alpha_dB_cm * 100 / 8.685889
    T_amp    = math.exp(-alpha_Np * L_um * 1e-6 / 2.0)
    L_cm     = L_um * 1e-4
    dphi     = math.pi * (V_bias - v_ref) * L_cm / vpi_L
    T_bar   = T_amp * T_amp * math.sin(dphi / 2) ** 2
    T_cross = T_amp * T_amp * math.cos(dphi / 2) ** 2
    return T_bar, T_cross

# ── Netlist builders ───────────────────────────────────────────────────────────

def _osdi(*names):
    return "\n".join(f".osdi {BUILD / n}.osdi" for n in names)

def netlist_l1():
    return f"""\
* MZI modulator L1 (single-arm PN)
{_osdi("cw_laser", "mzi_modulator_pn_l1", "photodetector")}
Xlaser   l_re l_im l_wl  cw_laser  power_mW=1.0 wavelength_nm={WL_NM}
Xmzi     l_re l_im l_wl  bar_re bar_im l_wl  cross_re cross_im l_wl \\
         vbias 0  mzi_modulator_pn_l1 \\
         L_arm_um={L_ARM_UM} n_g={N_G} alpha_dB_cm={ALPHA_DCM} \\
         Vpi_L={VPI_L} V_ref={V_REF} wavelength_nm={WL_NM}
Xpd_bar   bar_re bar_im l_wl  ph_bar   0  photodetector  responsivity=1.0
Xpd_cross cross_re cross_im l_wl  ph_cross 0  photodetector  responsivity=1.0
Rbar    ph_bar   0  1k
Rcross  ph_cross 0  1k
Vbias   vbias 0  DC 0.0
.optical l_re l_im l_wl bar_re bar_im cross_re cross_im
.op
.end
"""

def netlist_l2():
    # Push-pull: anode1=vbias, cathode1=0, anode2=0, cathode2=vbias
    # → V_pn1 = +V, V_pn2 = −V  (differential drive)
    return f"""\
* MZI modulator L2 (push-pull PN, Soref-Bennett)
{_osdi("cw_laser", "mzi_modulator_pn_l2", "photodetector")}
Xlaser   l_re l_im l_wl  cw_laser  power_mW=1.0 wavelength_nm={WL_NM}
Xmzi     l_re l_im l_wl  bar_re bar_im l_wl  cross_re cross_im l_wl \\
         vbias 0  0 vbias  mzi_modulator_pn_l2 \\
         L_arm_um={L_ARM_UM} n_g={N_G} alpha_dB_cm={ALPHA_DCM} wavelength_nm={WL_NM}
Xpd_bar   bar_re bar_im l_wl  ph_bar   0  photodetector  responsivity=1.0
Xpd_cross cross_re cross_im l_wl  ph_cross 0  photodetector  responsivity=1.0
Rbar    ph_bar   0  1k
Rcross  ph_cross 0  1k
Vbias   vbias 0  DC 0.0
.optical l_re l_im l_wl bar_re bar_im cross_re cross_im
.op
.end
"""

# ── Sweep ──────────────────────────────────────────────────────────────────────

NORM = 1.0e-3 * 1.0 * 1.0e3   # = 1.0 V at full transmission

def run_voltage_sweep(netlist_fn, v_array, label):
    """Sweep Vbias at fixed wavelength; return (T_bar, T_cross) arrays."""
    ckt = fc.Circuit()
    ckt.load_str(netlist_fn())
    ckt.set_param("Xlaser", "wavelength_nm", WL_NM)
    print(f"  {label}  sweeping {len(v_array)} voltage points ...", end="", flush=True)
    try:
        results = ckt.sweep("Vbias.dc", list(v_array), "op")
        T_bar   = np.array([r["V(ph_bar)"][0]   / NORM for r in results])
        T_cross = np.array([r["V(ph_cross)"][0] / NORM for r in results])
    except Exception as exc:
        print(f"  ERROR: {exc}")
        T_bar = T_cross = np.full(len(v_array), float("nan"))
    print(" done")
    return T_bar, T_cross

# ── Plot ───────────────────────────────────────────────────────────────────────

def plot(v_l1, Tb_l1, Tc_l1, v_l2, Tb_l2, Tc_l2):
    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=False)

    # L1 analytical
    T_bar_a  = np.array([mzi_l1_transfer(v, L_ARM_UM, N_G, ALPHA_DCM, VPI_L, V_REF, WL_NM)[0]
                          for v in v_l1])
    T_cross_a = np.array([mzi_l1_transfer(v, L_ARM_UM, N_G, ALPHA_DCM, VPI_L, V_REF, WL_NM)[1]
                          for v in v_l1])

    ax = axes[0]
    ax.plot(v_l1, Tb_l1,   "C0-",  lw=1.8, label="Bar (sim)")
    ax.plot(v_l1, Tc_l1,   "C1-",  lw=1.8, label="Cross (sim)")
    ax.plot(v_l1, T_bar_a,  "C0--", lw=1.0, alpha=0.6, label="Bar (CMT)")
    ax.plot(v_l1, T_cross_a,"C1--", lw=1.0, alpha=0.6, label="Cross (CMT)")
    ax.axvline(-VPI_L1, color="gray", ls=":", lw=1.0, label=f"Vπ = {VPI_L1:.0f} V")
    ax.set_xlabel("Bias voltage (V)")
    ax.set_ylabel("Transmission")
    ax.set_title(f"L1 — single-arm PN  (L={L_ARM_UM:.0f} µm, Vπ·L={VPI_L} V·cm, Vπ≈{VPI_L1:.0f} V)",
                 fontsize=10)
    ax.set_ylim(-0.05, 1.10)
    ax.legend(fontsize=9, ncol=2)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(v_l2, Tb_l2,  "C0-",  lw=1.8, label="Bar (sim)")
    ax.plot(v_l2, Tc_l2,  "C1-",  lw=1.8, label="Cross (sim)")
    ax.axvline(-VPI_L2_PP, color="gray", ls=":", lw=1.0,
               label=f"Vπ_pp ≈ {VPI_L2_PP:.0f} V (ideal)")
    ax.set_xlabel("Bias voltage (V)")
    ax.set_ylabel("Transmission")
    ax.set_title(
        f"L2 — push-pull PN  (Soref–Bennett, n_dep0=5e16 cm⁻³, Vbi=1 V)",
        fontsize=10,
    )
    ax.set_ylim(-0.05, 1.10)
    ax.legend(fontsize=9, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.text(0.98, 0.95,
            "Push-pull topology:\nanode1=V, cathode2=V\n(V_pn1=+V, V_pn2=−V)",
            transform=ax.transAxes, fontsize=8, ha="right", va="top", color="gray")

    fig.suptitle(f"MZI Modulator (PN) — Voltage Sweep @ λ = {WL_NM} nm\n"
                 f"L_arm = {L_ARM_UM:.0f} µm  n_g = {N_G}  α = {ALPHA_DCM} dB/cm",
                 fontsize=11)
    fig.tight_layout()

    out = HERE / "sweep_mzi_pn.png"
    fig.savefig(out, dpi=150)
    print(f"\nFigure saved: {out}")
    plt.close(fig)

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    for name in ("cw_laser", "mzi_modulator_pn_l1", "mzi_modulator_pn_l2",
                 "photodetector"):
        if not (BUILD / f"{name}.osdi").exists():
            sys.exit(f"Missing: {BUILD / name}.osdi — compile va-models first.")

    print(f"MZI L1 Vπ = {VPI_L1:.1f} V  (single-arm)")
    print(f"MZI L2 Vπ ≈ {VPI_L2_PP:.1f} V  (push-pull, ideal linear estimate)")

    print("\nSweeping L1 (single-arm PN)...")
    Tb1, Tc1 = run_voltage_sweep(netlist_l1, V_L1, "L1")

    print("Sweeping L2 (push-pull Soref–Bennett)...")
    Tb2, Tc2 = run_voltage_sweep(netlist_l2, V_L2, "L2")

    plot(V_L1, Tb1, Tc1, V_L2, Tb2, Tc2)


if __name__ == "__main__":
    main()
