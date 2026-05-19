"""
Femwell simulations for SOI waveguides:
1. Strip waveguide 500x220 nm at 1550 nm
2. Rib waveguide (500x220 nm core + 4500x90 nm slab) at 1550 nm
3. C-band sweep for both -> n_eff(lambda), n_g
4. Rib waveguide bent at 8 um radius -> n_eff(lambda)
"""
import numpy as np
from collections import OrderedDict
from shapely.geometry import box
from shapely.affinity import translate
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio
from femwell.mesh import mesh_from_OrderedDict
from femwell.maxwell.waveguide import compute_modes


# -------- Material dispersion ---------------------------------------------
def n_silicon(wl_um):
    """Salzberg & Villa Sellmeier for crystalline Si (valid 1.36-11 um)."""
    l = np.asarray(wl_um, dtype=float)
    n2 = 1 + (10.6684293 * l**2 / (l**2 - 0.301516485**2)
              + 0.003043475 * l**2 / (l**2 - 1.13475115**2)
              + 1.54133408 * l**2 / (l**2 - 1104.0**2))
    return np.sqrt(n2)

def n_sio2(wl_um):
    """Malitson Sellmeier for fused silica."""
    l = np.asarray(wl_um, dtype=float)
    n2 = 1 + (0.6961663 * l**2 / (l**2 - 0.0684043**2)
              + 0.4079426 * l**2 / (l**2 - 0.1162414**2)
              + 0.8974794 * l**2 / (l**2 - 9.896161**2))
    return np.sqrt(n2)


# -------- Geometry builders -----------------------------------------------
def build_strip(core_w=0.5, core_h=0.22, sim_w=4.0, sim_h=3.0,
                box_h=1.0, clad_h=1.0):
    """500x220 nm Si strip, fully clad by SiO2."""
    # center the core at x=0, y=0 (core bottom on the BOX interface at y=0)
    core = box(-core_w/2, 0, core_w/2, core_h)
    # BOX (buried oxide) below
    box_region = box(-sim_w/2, -box_h, sim_w/2, 0)
    # cladding everywhere else above BOX
    clad = box(-sim_w/2, -box_h, sim_w/2, core_h + clad_h)
    # subtract core from cladding
    clad = clad.difference(core).difference(box_region)
    polys = OrderedDict(core=core, box=box_region, clad=clad)
    resolutions = dict(core={"resolution": 0.02, "distance": 0.5},
                       box={"resolution": 0.1, "distance": 0.5},
                       clad={"resolution": 0.1, "distance": 0.5})
    return polys, resolutions

def build_rib(core_w=0.5, core_h=0.22, slab_w=4.5, slab_h=0.09,
              sim_w=6.0, box_h=1.0, clad_h=1.0):
    """Rib: 500x220 nm Si on top of 4500x90 nm Si slab."""
    # core sits on top of slab. slab base at y=0, slab top at y=slab_h.
    # core base at y=slab_h, top at y=slab_h+core_h
    slab = box(-slab_w/2, 0, slab_w/2, slab_h)
    core_only = box(-core_w/2, slab_h, core_w/2, slab_h + core_h)
    # the slab "wings" (slab minus the footprint under the core) - but we
    # actually want core+slab as separate regions, with core being Si and
    # slab being Si too (same material). It's simplest to define one Si region.
    si_region = slab.union(core_only)
    box_region = box(-sim_w/2, -box_h, sim_w/2, 0)
    clad = box(-sim_w/2, -box_h, sim_w/2, slab_h + core_h + clad_h)
    clad = clad.difference(si_region).difference(box_region)
    polys = OrderedDict(core=si_region, box=box_region, clad=clad)
    resolutions = dict(core={"resolution": 0.02, "distance": 0.5},
                       box={"resolution": 0.1, "distance": 0.5},
                       clad={"resolution": 0.1, "distance": 0.5})
    return polys, resolutions


def solve_neff(polys, resolutions, wavelength_um, num_modes=2, radius=np.inf,
               n_guess=None):
    """Mesh + compute modes. Returns list of n_eff (complex)."""
    mesh = from_meshio(mesh_from_OrderedDict(polys, resolutions,
                                              default_resolution_max=0.2))
    basis0 = Basis(mesh, ElementTriP0())
    epsilon = basis0.zeros(dtype=complex)
    n_si = n_silicon(wavelength_um)
    n_ox = n_sio2(wavelength_um)
    for subdomain, n_val in [("core", n_si), ("box", n_ox), ("clad", n_ox)]:
        epsilon[basis0.get_dofs(elements=subdomain)] = n_val**2
    modes = compute_modes(basis0, epsilon, wavelength=wavelength_um,
                          num_modes=num_modes, order=1, radius=radius,
                          n_guess=n_guess)
    return modes


if __name__ == "__main__":
    WL = 1.55
    print("="*70)
    print(f"Material indices at {WL} um:")
    print(f"  n_Si  = {n_silicon(WL):.4f}")
    print(f"  n_SiO2 = {n_sio2(WL):.4f}")
    print("="*70)

    print("\n-- Strip 500x220 nm --")
    polys, res = build_strip()
    modes = solve_neff(polys, res, WL, num_modes=2)
    for i, m in enumerate(modes):
        print(f"  mode {i}: n_eff = {np.real(m.n_eff):.4f} + "
              f"{np.imag(m.n_eff):.2e}j   "
              f"(TE frac = {m.te_fraction:.2f})")

    print("\n-- Rib (500x220 core, 4500x90 slab) --")
    polys, res = build_rib()
    modes = solve_neff(polys, res, WL, num_modes=2)
    for i, m in enumerate(modes):
        print(f"  mode {i}: n_eff = {np.real(m.n_eff):.4f} + "
              f"{np.imag(m.n_eff):.2e}j   "
              f"(TE frac = {m.te_fraction:.2f})")
