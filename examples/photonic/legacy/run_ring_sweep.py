#!/usr/bin/env python3
"""
Ring resonator wavelength sweep — parametric driver.

Compiles and runs the fairchild ring_sweep example with tunable physical
parameters, then plots V(ph_a) vs wavelength alongside the CMT analytical
transmission curve.

Usage:
    python examples/photonic/run_ring_sweep.py [options]

Examples:
    # Default parameters (100 µm ring, kappa=0.1)
    python examples/photonic/run_ring_sweep.py

    # High-Q ring: larger radius, lower coupling
    python examples/photonic/run_ring_sweep.py --kappa 0.02 --L-um 500 --n-points 201

    # Compare two coupling strengths side by side
    python examples/photonic/run_ring_sweep.py --kappa 0.05 --wl-center 1550
    python examples/photonic/run_ring_sweep.py --kappa 0.20 --wl-center 1550

    # Lossless ring
    python examples/photonic/run_ring_sweep.py --alpha 0.0

Requirements:
    - Compiled va-models (run `cd legacy/va-models && bash build.sh` first)
    - matplotlib  (pip install matplotlib)
"""

import argparse
import math
import csv
import os
import pathlib
import subprocess
import sys
import tempfile

# Repository root (two levels up from scripts/)
REPO_ROOT = pathlib.Path(__file__).parent.parent.parent.parent  # legacy/ -> photonic/ -> examples/ -> repo root

# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="Ring resonator wavelength sweep with tunable parameters.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--kappa",     type=float, default=0.1,   metavar="κ",
                   help="Coupler power cross-coupling fraction (0–1)")
    p.add_argument("--L-um",      type=float, default=100.0, metavar="µm",
                   help="Ring circumference (µm)")
    p.add_argument("--n-g",       type=float, default=4.2,   metavar="n_g",
                   help="Group index of the waveguide")
    p.add_argument("--alpha",     type=float, default=2.0,   metavar="dB/cm",
                   help="Propagation loss (dB/cm)")
    p.add_argument("--power",     type=float, default=1.0,   metavar="mW",
                   help="Laser power (mW)")
    p.add_argument("--r-load",    type=float, default=1e3,   metavar="Ω",
                   help="Photodetector load resistance (Ω)")
    p.add_argument("--wl-center", type=float, default=1551.0, metavar="nm",
                   help="Center wavelength for sweep (nm); sweep spans ±1.5×FSR")
    p.add_argument("--wl-start",  type=float, default=None,  metavar="nm",
                   help="Sweep start wavelength (nm); overrides --wl-center")
    p.add_argument("--wl-end",    type=float, default=None,  metavar="nm",
                   help="Sweep end wavelength (nm); overrides --wl-center")
    p.add_argument("--n-points",  type=int,   default=101,
                   help="Number of wavelength sweep points")
    p.add_argument("--model-dir", type=str,   default=None,
                   help="Directory containing .osdi model files (default: legacy/va-models/build)")
    p.add_argument("--no-sim",    action="store_true",
                   help="Skip simulation; plot CMT only")
    p.add_argument("--no-plot",   action="store_true",
                   help="Skip plot; only run simulation and print summary")
    p.add_argument("--csv",       type=str,   default=None,
                   help="Load existing CSV instead of running simulation")
    return p.parse_args()

# ── CMT helpers ───────────────────────────────────────────────────────────────
def cmt_fsr_nm(n_g, l_ring_um, wl_nm):
    """Free spectral range in nm."""
    return wl_nm**2 / (n_g * l_ring_um * 1e3)  # nm² / (n_g * µm * 1e3 nm/µm)

def cmt_resonance_nearest(wl_center_m, n_g, l_ring_m):
    m = round(n_g * l_ring_m / wl_center_m)
    return n_g * l_ring_m / m

def cmt_fwhm_nm(kappa_0, alpha_db_cm, l_ring_um, n_g, wl_nm):
    """Lorentzian FWHM of the resonance dip in nm (approximate)."""
    l_ring_m = l_ring_um * 1e-6
    alpha_lin = alpha_db_cm * 1e2 / 8.685895
    a = math.exp(-alpha_lin * l_ring_m / 2.0)
    r = math.sqrt(1.0 - kappa_0)
    # FWHM = FSR / (π * F),  F ≈ π*(r*a)^0.5 / (1 - r*a)
    ra = r * a
    finesse = math.pi * math.sqrt(ra) / (1 - ra) if ra < 1 else float("inf")
    fsr = cmt_fsr_nm(n_g, l_ring_um, wl_nm)
    return fsr / finesse if finesse > 0 else float("nan")

def cmt_t_min(kappa_0, alpha_db_cm, l_ring_um):
    l_ring_m = l_ring_um * 1e-6
    alpha_lin = alpha_db_cm * 1e2 / 8.685895
    a = math.exp(-alpha_lin * l_ring_m / 2.0)
    r = math.sqrt(1.0 - kappa_0)
    return (r - a)**2 / (1 - r*a)**2

def cmt_transmission(wl_m, kappa_0, l_ring_m, n_g, alpha_db_cm):
    r = math.sqrt(1.0 - kappa_0)
    alpha_lin = alpha_db_cm * 1e2 / 8.685895
    a = math.exp(-alpha_lin * l_ring_m / 2.0)
    beta = 2.0 * math.pi * n_g / wl_m
    phi = beta * l_ring_m
    return (r*r - 2*r*a*math.cos(phi) + a*a) / (1.0 - 2*r*a*math.cos(phi) + r*r*a*a)

# ── Simulation ────────────────────────────────────────────────────────────────
def run_simulation(args, wl_start, wl_end, csv_path):
    model_dir = args.model_dir or str(REPO_ROOT / "legacy" / "va-models" / "build")
    env = {
        **os.environ,
        "RING_KAPPA_0":     str(args.kappa),
        "RING_L_UM":        str(args.L_um),
        "RING_N_G":         str(args.n_g),
        "RING_ALPHA_DB_CM": str(args.alpha),
        "RING_POWER_MW":    str(args.power),
        "RING_R_LOAD":      str(args.r_load),
        "RING_WL_START_NM": str(wl_start),
        "RING_WL_END_NM":   str(wl_end),
        "RING_N_POINTS":    str(args.n_points),
        "RING_CSV_OUT":     str(csv_path),
        "RING_MODEL_DIR":   model_dir,
    }
    cmd = [
        "cargo", "run",
        "--example", "ring_sweep",
        "--manifest-path", str(REPO_ROOT / "crates" / "fairchild-osdi" / "Cargo.toml"),
        "--quiet",
    ]
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, env=env, cwd=str(REPO_ROOT))
    if result.returncode != 0:
        print("Simulation failed.", file=sys.stderr)
        sys.exit(1)

# ── Plot ──────────────────────────────────────────────────────────────────────
def plot(args, wl_nm_list, v_sim_list, t_cmt_list, wl_start, wl_end):
    try:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    except ImportError:
        print("matplotlib not installed — skipping plot (pip install matplotlib)")
        return

    l_ring_m  = args.L_um * 1e-6
    fsr       = cmt_fsr_nm(args.n_g, args.L_um, (wl_start + wl_end) / 2)
    fwhm      = cmt_fwhm_nm(args.kappa, args.alpha, args.L_um, args.n_g, (wl_start + wl_end) / 2)
    t_min_cmt = cmt_t_min(args.kappa, args.alpha, args.L_um)
    finesse   = fsr / fwhm if fwhm > 0 else float("nan")

    # Find simulated resonance(s) — annotate the deepest one
    v_max = max(v_sim_list)
    v_min = min(v_sim_list)
    sim_res_nm = wl_nm_list[v_sim_list.index(v_min)]
    cmt_res_nm = cmt_resonance_nearest(sim_res_nm * 1e-9, args.n_g, l_ring_m) * 1e9
    dip_pct = (1 - v_min / v_max) * 100

    # Scale CMT to match V axis
    t_max  = max(t_cmt_list)
    v_scale = v_max / t_max if t_max > 0 else 1.0
    t_scaled = [t * v_scale for t in t_cmt_list]

    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.plot(wl_nm_list, v_sim_list, color="C0", lw=1.8, label="Simulated V(ph_a)")
    ax.plot(wl_nm_list, t_scaled,   color="C1", lw=1.5, ls="--", alpha=0.85, label="CMT (scaled)")

    # Annotate simulated resonance
    ax.axvline(sim_res_nm, color="C0", lw=0.9, ls=":")
    ax.annotate(
        f"{sim_res_nm:.3f} nm\n(sim, Δ={abs(sim_res_nm-cmt_res_nm):.3f} nm from CMT)",
        xy=(sim_res_nm, v_min), xytext=(sim_res_nm + fsr * 0.15, v_min + (v_max - v_min) * 0.25),
        arrowprops=dict(arrowstyle="->", color="C0"), color="C0", fontsize=8,
    )

    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("V(ph_a)  [V]")
    ax.set_xlim(wl_start, wl_end)
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.grid(True, which="major", alpha=0.3)
    ax.grid(True, which="minor", alpha=0.1)
    ax.legend(loc="upper right")

    param_str = (
        f"L={args.L_um:.0f} µm  n_g={args.n_g}  κ={args.kappa}  "
        f"α={args.alpha} dB/cm  P={args.power} mW\n"
        f"FSR={fsr:.2f} nm  FWHM={fwhm:.3f} nm  F={finesse:.0f}  "
        f"T_min(CMT)={t_min_cmt:.3f}  dip(sim)={dip_pct:.1f}%"
    )
    ax.set_title(f"Ring resonator — {param_str}", fontsize=9)
    fig.tight_layout()

    out_dir = REPO_ROOT / "docs" / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    slug = f"ring_kappa{args.kappa}_L{args.L_um:.0f}um_ng{args.n_g}_alpha{args.alpha}"
    out_path = out_dir / f"{slug}.png"
    fig.savefig(out_path, dpi=150)
    print(f"Plot saved to {out_path}")
    plt.show()

# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()

    # Determine sweep range
    l_ring_m = args.L_um * 1e-6
    if args.wl_start is not None and args.wl_end is not None:
        wl_start, wl_end = args.wl_start, args.wl_end
    else:
        # Auto-range: ±1.5 FSR around the nearest resonance to wl_center
        wl_center = args.wl_center
        fsr = cmt_fsr_nm(args.n_g, args.L_um, wl_center)
        res = cmt_resonance_nearest(wl_center * 1e-9, args.n_g, l_ring_m) * 1e9
        wl_start = res - 1.5 * fsr
        wl_end   = res + 1.5 * fsr

    # Print analytical predictions
    fsr  = cmt_fsr_nm(args.n_g, args.L_um, (wl_start + wl_end) / 2)
    fwhm = cmt_fwhm_nm(args.kappa, args.alpha, args.L_um, args.n_g, (wl_start + wl_end) / 2)
    print(f"\nAnalytical (CMT):")
    print(f"  FSR   = {fsr:.3f} nm")
    print(f"  FWHM  = {fwhm:.4f} nm")
    print(f"  F     = {fsr/fwhm:.1f}" if fwhm > 0 else "  F     = ∞")
    print(f"  T_min = {cmt_t_min(args.kappa, args.alpha, args.L_um):.4f}")
    print(f"  Sweep: {wl_start:.3f}–{wl_end:.3f} nm  ({args.n_points} points)\n")

    # Load or run simulation
    wl_nm_list, v_sim_list, t_cmt_list = [], [], []

    if args.csv:
        csv_path = pathlib.Path(args.csv)
        print(f"Loading existing CSV: {csv_path}")
    elif args.no_sim:
        csv_path = None
    else:
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
            csv_path = pathlib.Path(tmp.name)
        run_simulation(args, wl_start, wl_end, csv_path)

    if csv_path and csv_path.exists():
        with open(csv_path) as f:
            for row in csv.DictReader(f):
                wl_nm_list.append(float(row["wavelength_nm"]))
                v_sim_list.append(float(row["V_ph_a_V"]))
                t_cmt_list.append(float(row["T_cmt"]))

    if not wl_nm_list:
        # CMT-only mode: generate analytical curve
        step = (wl_end - wl_start) / (args.n_points - 1)
        wl_nm_list  = [wl_start + i * step for i in range(args.n_points)]
        t_cmt_list  = [cmt_transmission(wl * 1e-9, args.kappa, l_ring_m, args.n_g, args.alpha)
                       for wl in wl_nm_list]
        # Scale CMT to V axis for consistent y-axis
        p_in_w = args.power * 1e-3
        v_sim_list = [t * p_in_w * args.r_load for t in t_cmt_list]

    if not args.no_plot:
        plot(args, wl_nm_list, v_sim_list, t_cmt_list, wl_start, wl_end)

if __name__ == "__main__":
    main()
