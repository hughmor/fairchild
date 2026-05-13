#!/usr/bin/env python3
"""
Plot ring resonator wavelength sweep: simulated V(ph_a) vs CMT transmission.

Usage:
    python scripts/plot_ring_sweep.py [sweep.csv]

The CSV is produced by the ring_resonator_wavelength_sweep test:
    cargo test -p fairchild-osdi --test ring_resonator ring_resonator_wavelength_sweep

Default input path: /tmp/ring_resonator_sweep.csv
"""

import sys
import math
import csv
import pathlib

# ── physical parameters (must match ring_resonator.rs constants) ────────────
L_RING_UM   = 100.0
N_G         = 4.2
KAPPA_0     = 0.1
ALPHA_DB_CM = 2.0
POWER_MW    = 1.0
R_LOAD      = 1e3   # Ω

def cmt_transmission(wavelength_m: float) -> float:
    r = math.sqrt(1.0 - KAPPA_0)
    alpha_lin = ALPHA_DB_CM * 1e2 / 8.685895   # dB/cm → Np/m
    l_ring_m  = L_RING_UM * 1e-6
    a   = math.exp(-alpha_lin * l_ring_m / 2.0)
    beta = 2.0 * math.pi * N_G / wavelength_m
    phi  = beta * l_ring_m
    return (r*r - 2*r*a*math.cos(phi) + a*a) / (1.0 - 2*r*a*math.cos(phi) + r*r*a*a)

def cmt_resonance_nearest(lambda_center_m: float) -> float:
    l_ring_m = L_RING_UM * 1e-6
    m = round(N_G * l_ring_m / lambda_center_m)
    return N_G * l_ring_m / m

# ── load CSV ─────────────────────────────────────────────────────────────────
csv_path = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path("/tmp/ring_resonator_sweep.csv")

if not csv_path.exists():
    print(f"CSV not found: {csv_path}")
    print("Run: cargo test -p fairchild-osdi --test ring_resonator ring_resonator_wavelength_sweep")
    sys.exit(1)

wl_nm, v_sim, t_cmt = [], [], []
with open(csv_path) as f:
    for row in csv.DictReader(f):
        wl = float(row["wavelength_nm"])
        wl_nm.append(wl)
        v_sim.append(float(row["V_ph_a_V"]))
        t_cmt.append(cmt_transmission(wl * 1e-9))

# ── derived quantities ────────────────────────────────────────────────────────
v_max      = max(v_sim)
v_min      = min(v_sim)
sim_res_nm = wl_nm[v_sim.index(v_min)]
cmt_res_nm = cmt_resonance_nearest(sim_res_nm * 1e-9) * 1e9

# Normalise CMT to match simulated V scale (V_max ≈ T_max * P_in * R_load)
t_max  = max(t_cmt)
v_scale = v_max / t_max if t_max > 0 else 1.0
t_scaled = [t * v_scale for t in t_cmt]

print(f"Sweep: {wl_nm[0]:.1f}–{wl_nm[-1]:.1f} nm  ({len(wl_nm)} points)")
print(f"V_max = {v_max:.4f} V  V_min = {v_min:.4f} V  dip = {(1-v_min/v_max)*100:.1f}%")
print(f"Simulated resonance: {sim_res_nm:.3f} nm")
print(f"CMT resonance:       {cmt_res_nm:.3f} nm")
print(f"Δλ = {abs(sim_res_nm - cmt_res_nm):.4f} nm  (tolerance: 0.1 nm)")

# ── plot ──────────────────────────────────────────────────────────────────────
try:
    import matplotlib.pyplot as plt
except ImportError:
    print("\nmatplotlib not installed — skipping plot. Install with: pip install matplotlib")
    sys.exit(0)

fig, ax = plt.subplots(figsize=(9, 4))

ax.plot(wl_nm, v_sim,    color="C0", lw=1.5, label="Simulated V(ph_a)")
ax.plot(wl_nm, t_scaled, color="C1", lw=1.5, linestyle="--", label="CMT (scaled)")

ax.axvline(sim_res_nm, color="C0", lw=0.8, linestyle=":")
ax.axvline(cmt_res_nm, color="C1", lw=0.8, linestyle=":")

ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel("V(ph_a)  [V]")
ax.set_title(
    f"Ring resonator (L={L_RING_UM:.0f} µm, n_g={N_G}, κ={KAPPA_0}, α={ALPHA_DB_CM} dB/cm)\n"
    f"Simulated resonance {sim_res_nm:.3f} nm  |  CMT {cmt_res_nm:.3f} nm  |  Δλ = {abs(sim_res_nm-cmt_res_nm):.4f} nm"
)
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()

_repo_root = pathlib.Path(__file__).parent.parent.parent
out = _repo_root / "docs" / "plots" / "ring_resonator_sweep.png"
fig.savefig(out, dpi=150)
print(f"\nPlot saved to {out}")
plt.show()
