"""Two-waveguide supermode coupling for a circular MRR.

Computes κ(g) for two parallel rib waveguides via the supermode method:
    κ(g) = π · (n_sym − n_asym) / λ        (per unit length)
Then integrates κ(g(z)) along the circular ring's approach/recession to
get the total coupling angle κ·L_eff for a "lumped" directional-coupler
model.

For a ring of radius R coupled to a straight bus with minimum gap g_0,
the centre-to-centre gap along the bus is
    g(z) = g_0 + R · (1 − cos(z/R))   ≈   g_0 + z²/(2R)   for |z| ≪ R.

Because κ falls off roughly exponentially with gap (κ ≈ κ_0·exp(−g/L_d)),
the integral has an approximately Gaussian profile in z; we just do it
numerically from the fitted κ(g).

Output: prints κ·L_eff (the value fc_dcoupler needs as `kappa_L`) and
an effective coupler length suggestion.
"""
import numpy as np
from collections import OrderedDict
from shapely.geometry import box
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio
from femwell.mesh import mesh_from_OrderedDict
from femwell.maxwell.waveguide import compute_modes
from wg_sim import n_silicon, n_sio2


# ------------------------------------------------------------------ geometry
def build_two_ribs(gap_um, core_w=0.5, core_h=0.22, slab_w=8.0, slab_h=0.09,
                   sim_w=14.0, box_h=1.0, clad_h=1.0):
    """Two identical 500x220 rib waveguides separated edge-to-edge by `gap_um`.

    Slab is continuous under both cores (standard rib).  Mesh is fine in
    the gap region (resolution 10 nm) where the mode tails matter most.
    """
    # Centre-to-centre of the two cores:  d = gap + core_w
    d = gap_um + core_w
    x_left  = -d / 2.0
    x_right = +d / 2.0
    slab = box(-slab_w/2, 0, slab_w/2, slab_h)
    core_left  = box(x_left - core_w/2,  slab_h, x_left + core_w/2,  slab_h + core_h)
    core_right = box(x_right - core_w/2, slab_h, x_right + core_w/2, slab_h + core_h)
    si = slab.union(core_left).union(core_right)
    box_region = box(-sim_w/2, -box_h, sim_w/2, 0)
    clad = box(-sim_w/2, -box_h, sim_w/2, slab_h + core_h + clad_h)
    clad = clad.difference(si).difference(box_region)
    polys = OrderedDict(core=si, box=box_region, clad=clad)
    # Finer resolution in the gap so the evanescent overlap is captured.
    resolutions = dict(
        core={"resolution": 0.02, "distance": 0.4},
        box ={"resolution": 0.15, "distance": 0.4},
        clad={"resolution": 0.15, "distance": 0.4},
    )
    return polys, resolutions


def supermodes(gap_um, wavelength_um=1.55, n_guess=2.76):
    """Return (n_sym, n_asym) for the lowest two TE-like supermodes."""
    polys, resolutions = build_two_ribs(gap_um)
    mesh = from_meshio(mesh_from_OrderedDict(polys, resolutions,
                                              default_resolution_max=0.25))
    basis0 = Basis(mesh, ElementTriP0())
    epsilon = basis0.zeros(dtype=complex)
    n_si = n_silicon(wavelength_um)
    n_ox = n_sio2(wavelength_um)
    for subdomain, n_val in [("core", n_si), ("box", n_ox), ("clad", n_ox)]:
        epsilon[basis0.get_dofs(elements=subdomain)] = n_val**2
    modes = compute_modes(basis0, epsilon, wavelength=wavelength_um,
                          num_modes=2, order=1, n_guess=n_guess)
    # Sort descending in real(n_eff): symmetric mode has higher n_eff.
    neffs = sorted([complex(m.n_eff) for m in modes], key=lambda c: -c.real)
    return neffs[0].real, neffs[1].real


# ------------------------------------------------------------------ kappa(g)
def kappa_from_supermodes(n_sym, n_asym, wavelength_um=1.55):
    """κ = π · Δn / λ  (in 1/µm if λ is in µm)."""
    return np.pi * (n_sym - n_asym) / wavelength_um


# ------------------------------------------------------------------ ring integration
def kappa_L_circular(kappa_g_fn, R_um, g0_um, span_um=None, n_pts=801):
    """Integrate κ(g(z)) along z for a ring of radius R_um, min gap g0_um.

    g(z) = g0 + R · (1 − cos(z/R))
    The integration window auto-sizes to where κ has dropped to 1 % of peak.
    """
    if span_um is None:
        # Find z where g(z) - g0 has caused κ to drop ~100×; an order-of-mag
        # bigger gap typically reduces κ by ~100×.  Take z where g = g0 + 1 µm.
        # g0 + R(1-cos(z/R)) = g0 + 1  →  z = R·arccos(1 − 1/R).
        # For R = 8 µm that's z ≈ 4.0 µm; double it as a safety margin.
        span_um = 2.0 * R_um * np.arccos(max(1.0 - 1.0/R_um, -1.0))
    z = np.linspace(-span_um, span_um, n_pts)
    g = g0_um + R_um * (1.0 - np.cos(z / R_um))
    k = kappa_g_fn(g)
    kL = np.trapz(k, z)
    return kL, z, g, k


# ------------------------------------------------------------------ main
def main():
    WL = 1.55
    print(f"Two-rib supermode sweep at λ = {WL} µm")
    print("="*70)

    # Gap sweep: tight cluster near 300 nm (the design point), plus
    # wider gaps so we can fit κ(g) = κ0·exp(−g/L_d) cleanly.
    gaps_um = np.array([0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.80, 1.00])
    kappas  = np.zeros_like(gaps_um)
    for i, g in enumerate(gaps_um):
        n_sym, n_asym = supermodes(g)
        kappas[i] = kappa_from_supermodes(n_sym, n_asym, WL)
        print(f"  gap = {g*1000:6.1f} nm:  n_s = {n_sym:.5f}  n_a = {n_asym:.5f}  "
              f"κ = {kappas[i]:.4f} /µm = {kappas[i]*1e6:.3e} /m")

    # Fit κ(g) = κ0 · exp(−g/L_d) on the log domain.
    log_k = np.log(kappas)
    slope, intercept = np.polyfit(gaps_um, log_k, 1)
    L_d = -1.0 / slope
    kappa_0 = np.exp(intercept)
    print()
    print(f"Exponential fit:  κ(g) = {kappa_0:.4f} · exp(−g / {L_d:.4f} µm)  (g in µm)")
    print(f"  κ(300 nm) measured = {kappas[gaps_um==0.30][0]:.4f}  "
          f"vs fit = {kappa_0*np.exp(-0.30/L_d):.4f}")

    def kappa_fn(g_um):
        return kappa_0 * np.exp(-g_um / L_d)

    # ----- Integrate over the circular ring approach -----
    # From the giona schematic: PS_PN_TH L = 25.13 µm, two arms in the ring
    # → ring circumference C = 2 · 25.13 = 50.26 µm → R = C / (2π) ≈ 8 µm.
    R_um = 8.0
    g0_um = 0.30
    kL, z, g, k_of_z = kappa_L_circular(kappa_fn, R_um, g0_um)

    # Recommended effective coupler length: take where κ has dropped to 1 % of
    # the peak.  This is purely a documentation choice (kappa·L is what matters
    # in the cos/sin transfer).
    peak = kappa_fn(g0_um)
    mask = k_of_z >= 0.01 * peak
    if mask.any():
        L_eff_um = float(z[mask][-1] - z[mask][0])
    else:
        L_eff_um = 2.0 * np.sqrt(2.0 * R_um * L_d)  # Gaussian σ rule

    print()
    print("="*70)
    print(f"Ring R = {R_um} µm, minimum gap g₀ = {g0_um*1000:.0f} nm")
    print("="*70)
    print(f"  κ·L_total (point-coupler equivalent) = {kL:.5f} rad")
    print(f"    → power cross-coupling at resonance: |sin(κL)|² = "
          f"{np.sin(kL)**2:.5f}")
    print(f"    → through-coupling:                |cos(κL)|² = "
          f"{np.cos(kL)**2:.5f}")
    print(f"  L_eff (1 % κ-floor)                 = {L_eff_um:.3f} µm")
    print(f"  kappa_per_m (= κL/L_eff)            = "
          f"{kL/(L_eff_um*1e-6):.3e} /m")
    print()
    print("→ fc_dcoupler defaults to set:")
    print(f"     kappa_L = {kL:.5f}")
    print(f"     L_um    = {L_eff_um:.3f}     (so kappa = {kL/(L_eff_um*1e-6):.4g} /m)")

    # Save CSV of κ(g) table for reproducibility / regression.
    import csv, os
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "coupler_kappa_vs_gap.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gap_um", "kappa_per_um", "kappa_per_m"])
        for gi, ki in zip(gaps_um, kappas):
            w.writerow([f"{gi:.4f}", f"{ki:.6e}", f"{ki*1e6:.6e}"])
        w.writerow([])
        w.writerow(["# ring sim:"])
        w.writerow(["R_um", f"{R_um}"])
        w.writerow(["g0_um", f"{g0_um}"])
        w.writerow(["kappa_L_rad", f"{kL:.6f}"])
        w.writerow(["L_eff_um", f"{L_eff_um:.3f}"])
        w.writerow(["kappa_0_per_um", f"{kappa_0:.6f}"])
        w.writerow(["L_d_um", f"{L_d:.6f}"])
    print(f"\nSaved κ(g) table + ring integration to coupler_kappa_vs_gap.csv")

    # Plot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    g_fine = np.linspace(gaps_um.min(), gaps_um.max(), 200)
    axes[0].semilogy(gaps_um*1000, kappas, "o", label="femwell")
    axes[0].semilogy(g_fine*1000, kappa_fn(g_fine), "-",
                     label=f"fit: {kappa_0:.3f}·exp(−g/{L_d*1000:.0f}nm)")
    axes[0].set_xlabel("Gap [nm]")
    axes[0].set_ylabel("κ [rad/µm]")
    axes[0].set_title("Coupling κ(g) — two parallel rib WG")
    axes[0].legend(); axes[0].grid(alpha=0.3, which="both")

    axes[1].plot(z, k_of_z*1000, label=f"κ(z), R={R_um} µm, g₀={g0_um*1000:.0f} nm")
    axes[1].axhline(0.01*peak*1000, color="grey", ls="--", lw=0.5,
                    label="1 % of peak")
    axes[1].set_xlabel("z along bus [µm]")
    axes[1].set_ylabel("κ(z) [mrad/µm]")
    axes[1].set_title(f"κ·L_total = {kL:.4f} rad   "
                      f"→ |sin(κL)|² = {np.sin(kL)**2:.4f}")
    axes[1].legend(); axes[1].grid(alpha=0.3)

    plt.tight_layout()
    out = os.path.join(here, "coupler_sim.png")
    plt.savefig(out, dpi=130)
    print(f"Plot saved to {out}")


if __name__ == "__main__":
    main()
