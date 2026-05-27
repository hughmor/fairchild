#!/usr/bin/env python3
"""
MZI thermo-optic modulator — heater-voltage-domain characterization.

Sweeps heater voltage at λ = 1550 nm for the L1 (single-arm, integrated
thermal) and L2 (dual-arm, external T_node) thermo-optic MZI models.  Each
subplot shows bar (solid) and cross (dashed) transmission vs heater voltage,
with an analytical CMT overlay.

L1: integrated thermal model.  Vπ ≈ 0.41 V with R_heater=1000 Ω,
    R_thermal=50000 K/W, L_arm=500 µm, dn/dT=1.86×10⁻⁴ K⁻¹.
    At V_heat=0: cross is bright (Δφ=0).
    At V_heat=Vπ: bar is bright (Δφ=π).

L2: external T_node thermal ports (arm 1 driven, arm 2 at 0 V).
    Same steady-state physics as L1 with Rth = 50 kΩ to ground.

Requirements: fairchild Python package (maturin develop), numpy, matplotlib.
Compiled OSDI models in legacy/va-models/build/.
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
BUILD = HERE.parents[2] / "legacy" / "va-models" / "build"

# ── MZI / heater parameters ────────────────────────────────────────────────────

L_ARM_UM   = 500.0       # µm
N_G        = 4.2
ALPHA_DCM  = 3.0         # dB/cm
R_HEATER   = 1000.0      # Ω
R_THERMAL  = 50000.0     # K/W (integrated in L1; external Rth in L2)
DN_DT      = 1.86e-4     # K⁻¹ (Si thermo-optic at 1550 nm)
WL_NM      = 1550.0

# Vπ: solve π = 2π * dn_dT * (Vπ²/R_h) * R_th * L / λ
# Vπ = sqrt(π / (2π × dn_dT × R_th/R_h × L/λ))
_fac   = 2 * math.pi * DN_DT * (R_THERMAL / R_HEATER) * (L_ARM_UM * 1e-6) / (WL_NM * 1e-9)
VPI    = math.sqrt(math.pi / _fac)   # ≈ 0.41 V

# Heater voltage sweep: 0 → 2.5 × Vπ  (shows >1 full period)
V_HEAT = np.linspace(0.0, 2.5 * VPI, 101)

# ── Analytical MZI thermo transfer function ────────────────────────────────────

def mzi_thermo_transfer(V_heat, L_um, n_g, alpha_dB_cm, r_heater, r_thermal,
                        dn_dt, wl_nm):
    """Analytical bar/cross power transmission for thermo MZI (matches VA model)."""
    alpha_Np = alpha_dB_cm * 100 / 8.685889
    T_amp    = math.exp(-alpha_Np * L_um * 1e-6 / 2.0)
    P_heat   = V_heat * V_heat / r_heater
    dT       = P_heat * r_thermal
    dphi     = 2 * math.pi * dn_dt * dT * L_um * 1e-6 / (wl_nm * 1e-9)
    T_bar   = T_amp * T_amp * math.sin(dphi / 2) ** 2
    T_cross = T_amp * T_amp * math.cos(dphi / 2) ** 2
    return T_bar, T_cross

# ── Netlist builders ───────────────────────────────────────────────────────────

def _osdi(*names):
    return "\n".join(f".osdi {BUILD / n}.osdi" for n in names)

def netlist_l1():
    return f"""\
* MZI thermo-optic modulator L1 (integrated thermal, single arm)
{_osdi("cw_laser", "mzi_modulator_thermo_l1", "photodetector")}
Xlaser   l_re l_im l_wl  cw_laser  power_mW=1.0 wavelength_nm={WL_NM}
Xmzi     l_re l_im l_wl  bar_re bar_im l_wl  cross_re cross_im l_wl \\
         hp hn  mzi_modulator_thermo_l1 \\
         L_arm_um={L_ARM_UM} n_g={N_G} alpha_dB_cm={ALPHA_DCM} \\
         R_heater={R_HEATER} R_thermal={int(R_THERMAL)} dn_dT={DN_DT} \\
         wavelength_nm={WL_NM}
Xpd_bar   bar_re bar_im l_wl  ph_bar   0  photodetector  responsivity=1.0
Xpd_cross cross_re cross_im l_wl  ph_cross 0  photodetector  responsivity=1.0
Rbar    ph_bar   0  1k
Rcross  ph_cross 0  1k
Vheat   hp hn  DC 0.0
.optical l_re l_im l_wl bar_re bar_im cross_re cross_im
.op
.end
"""

def netlist_l2():
    # Arm 1 driven by Vheat; arm 2 at 0 V.
    # Both T_nodes connected to external Rth.
    return f"""\
* MZI thermo-optic modulator L2 (external T_nodes; arm 1 driven)
{_osdi("cw_laser", "mzi_modulator_thermo_l2", "photodetector")}
Xlaser   l_re l_im l_wl  cw_laser  power_mW=1.0 wavelength_nm={WL_NM}
Xmzi     l_re l_im l_wl  bar_re bar_im l_wl  cross_re cross_im l_wl \\
         hp1 hn1 T1  hp2 hn2 T2  mzi_modulator_thermo_l2 \\
         L_arm_um={L_ARM_UM} n_g={N_G} alpha_dB_cm={ALPHA_DCM} \\
         R_heater={R_HEATER} dn_dT={DN_DT} wavelength_nm={WL_NM}
Rth1    T1 0  {int(R_THERMAL)}
Rth2    T2 0  {int(R_THERMAL)}
Xpd_bar   bar_re bar_im l_wl  ph_bar   0  photodetector  responsivity=1.0
Xpd_cross cross_re cross_im l_wl  ph_cross 0  photodetector  responsivity=1.0
Rbar    ph_bar   0  1k
Rcross  ph_cross 0  1k
Vheat   hp1 hn1  DC 0.0
Vheat2  hp2 hn2  DC 0.0
.optical l_re l_im l_wl bar_re bar_im cross_re cross_im
.op
.end
"""

# ── Sweep ──────────────────────────────────────────────────────────────────────

NORM = 1.0e-3 * 1.0 * 1.0e3

def run_voltage_sweep(netlist_fn, v_heat_name, label):
    ckt = fc.Circuit()
    ckt.load_str(netlist_fn())
    ckt.set_param("Xlaser", "wavelength_nm", WL_NM)
    print(f"  {label}  sweeping {len(V_HEAT)} heater-voltage points ...", end="", flush=True)
    try:
        results = ckt.sweep(f"{v_heat_name}.dc", list(V_HEAT), "op")
        T_bar   = np.array([r["V(ph_bar)"][0]   / NORM for r in results])
        T_cross = np.array([r["V(ph_cross)"][0] / NORM for r in results])
    except Exception as exc:
        print(f"  ERROR: {exc}")
        T_bar = T_cross = np.full(len(V_HEAT), float("nan"))
    print(" done")
    return T_bar, T_cross

# ── Plot ───────────────────────────────────────────────────────────────────────

def plot(Tb_l1, Tc_l1, Tb_l2, Tc_l2):
    V = V_HEAT

    # Analytical reference
    analytic = [mzi_thermo_transfer(v, L_ARM_UM, N_G, ALPHA_DCM,
                                    R_HEATER, R_THERMAL, DN_DT, WL_NM) for v in V]
    T_bar_a  = np.array([x[0] for x in analytic])
    T_cross_a = np.array([x[1] for x in analytic])

    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    titles = [
        f"L1 — single-arm, integrated thermal  "
        f"(R_thermal={R_THERMAL/1000:.0f} kK/W, R_heater={R_HEATER:.0f} Ω)",
        f"L2 — arm 1 driven, external T_node  "
        f"(Rth={R_THERMAL/1000:.0f} kΩ each arm)",
    ]

    for ax, Tb, Tc, title in zip(axes, [Tb_l1, Tb_l2], [Tc_l1, Tc_l2], titles):
        ax.plot(V, Tb,      "C0-",  lw=1.8, label="Bar (sim)")
        ax.plot(V, Tc,      "C1-",  lw=1.8, label="Cross (sim)")
        ax.plot(V, T_bar_a, "C0--", lw=1.0, alpha=0.6, label="Bar (analytic)")
        ax.plot(V, T_cross_a,"C1--",lw=1.0, alpha=0.6, label="Cross (analytic)")
        ax.axvline(VPI, color="gray", ls=":", lw=1.0, label=f"Vπ ≈ {VPI:.3f} V")
        ax.axvline(2 * VPI, color="gray", ls=":", lw=0.7, alpha=0.5)
        ax.set_ylabel("Transmission")
        ax.set_title(title, fontsize=10)
        ax.set_ylim(-0.05, 1.10)
        ax.legend(fontsize=9, ncol=2)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Heater voltage (V)")
    param_str = (f"L_arm = {L_ARM_UM:.0f} µm  n_g = {N_G}  "
                 f"α = {ALPHA_DCM} dB/cm  dn/dT = {DN_DT:.2e} K⁻¹  Vπ ≈ {VPI:.3f} V")
    fig.suptitle(f"MZI Thermo-Optic Modulator — Heater Voltage Sweep @ λ = {WL_NM} nm\n"
                 f"{param_str}", fontsize=11)
    fig.tight_layout()

    out = HERE / "sweep_mzi_thermo.png"
    fig.savefig(out, dpi=150)
    print(f"\nFigure saved: {out}")
    plt.close(fig)

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    for name in ("cw_laser", "mzi_modulator_thermo_l1", "mzi_modulator_thermo_l2",
                 "photodetector"):
        if not (BUILD / f"{name}.osdi").exists():
            sys.exit(f"Missing: {BUILD / name}.osdi — compile legacy/va-models first.")

    print(f"MZI thermo Vπ ≈ {VPI:.3f} V")

    print("\nSweeping L1 (integrated thermal)...")
    Tb1, Tc1 = run_voltage_sweep(netlist_l1, "Vheat", "L1")

    print("Sweeping L2 (external T_node)...")
    Tb2, Tc2 = run_voltage_sweep(netlist_l2, "Vheat", "L2")

    plot(Tb1, Tc1, Tb2, Tc2)


if __name__ == "__main__":
    main()
