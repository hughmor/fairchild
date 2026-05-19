"""
Template for simulating doped (N or P) and PN-junction SOI waveguides
in femwell, with index perturbation + loss via the Soref-Bennett model.

Fill in:
  - N_e_cm3, N_h_cm3: free electron / hole concentrations in cm^-3
  - For PN: define geometry of N-region, P-region, depletion region

Soref & Bennett (1987), updated by Nedeljkovic et al. (2011) for 1550 nm:
    Delta_n_e(N_e) = -8.8e-22 * N_e
    Delta_n_h(N_h) = -8.5e-18 * N_h^0.8        (N in cm^-3)
    Delta_alpha_e(N_e) = 8.5e-18 * N_e         [1/cm]
    Delta_alpha_h(N_h) = 6.0e-18 * N_h         [1/cm]

Total complex perm. of doped Si:
    n_doped = n_Si + Delta_n_e + Delta_n_h
    alpha   = Delta_alpha_e + Delta_alpha_h    [1/cm]   (intensity attenuation)
    k       = alpha * lambda / (4*pi)   with lambda in cm  -> imag part of n
    eps_r   = (n_doped + 1j*k)**2
"""
import numpy as np
from collections import OrderedDict
from shapely.geometry import box
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio
from femwell.mesh import mesh_from_OrderedDict
from femwell.maxwell.waveguide import compute_modes
from wg_sim import n_silicon, n_sio2


# ----- Soref-Bennett carrier-induced index/loss at 1550 nm ----------------
def soref_bennett(N_e_cm3, N_h_cm3, wavelength_um=1.55):
    """Returns (dn_e, dn_h, alpha_e, alpha_h).
    alpha_* in 1/cm (intensity).  Values for 1550 nm; use Nedeljkovic 2011
    for other wavelengths."""
    dn_e = -8.8e-22 * N_e_cm3
    dn_h = -8.5e-18 * (N_h_cm3 ** 0.8) if N_h_cm3 > 0 else 0.0
    a_e  = 8.5e-18 * N_e_cm3
    a_h  = 6.0e-18 * N_h_cm3
    return dn_e, dn_h, a_e, a_h


def doped_silicon_eps(N_e, N_h, wavelength_um=1.55):
    """Complex permittivity of doped silicon."""
    n_si = n_silicon(wavelength_um)
    dn_e, dn_h, a_e, a_h = soref_bennett(N_e, N_h, wavelength_um)
    n_real = n_si + dn_e + dn_h
    alpha_total = a_e + a_h                          # 1/cm
    # k from intensity attenuation: alpha[1/cm] = 4*pi*k/lambda[cm]
    lam_cm = wavelength_um * 1e-4
    k = alpha_total * lam_cm / (4 * np.pi)
    return (n_real + 1j*k) ** 2


# ----- Geometry for a PN-junction rib waveguide ---------------------------
def build_pn_rib(N_donor=1e18, N_acceptor=1e18,
                 core_w=0.5, core_h=0.22,
                 slab_w=4.5, slab_h=0.09,
                 junction_x=0.0,        # x-position of metallurgical junction
                 W_depl_n=0.05, W_depl_p=0.05,  # depletion widths (um)
                 sim_w=6.0, box_h=1.0, clad_h=1.0):
    """Rib waveguide split into N, P, and intrinsic depletion sub-regions.
    Each subdomain has its own (epsilon) set later via doped_silicon_eps.
    """
    # geometry: slab on bottom (y=[0, slab_h]), core on top (y=[slab_h, slab_h+core_h])
    y_top = slab_h + core_h
    # Vertical full-height silicon column split horizontally:
    # P region: x < junction_x - W_depl_p
    # depletion: junction_x - W_depl_p < x < junction_x + W_depl_n
    # N region: x > junction_x + W_depl_n
    p_x_max = junction_x - W_depl_p
    n_x_min = junction_x + W_depl_n
    # Union of slab + core for silicon outline
    si_outline = box(-slab_w/2, 0, slab_w/2, slab_h).union(
                    box(-core_w/2, slab_h, core_w/2, y_top))
    # Split into 3 regions via box intersection
    region_p   = si_outline.intersection(box(-slab_w/2, 0, p_x_max, y_top))
    region_dep = si_outline.intersection(box(p_x_max,    0, n_x_min, y_top))
    region_n   = si_outline.intersection(box(n_x_min,    0, slab_w/2, y_top))

    box_region = box(-sim_w/2, -box_h, sim_w/2, 0)
    clad = box(-sim_w/2, -box_h, sim_w/2, y_top + clad_h)
    clad = clad.difference(si_outline).difference(box_region)

    polys = OrderedDict(p_region=region_p,
                        depletion=region_dep,
                        n_region=region_n,
                        box=box_region,
                        clad=clad)
    resolutions = dict(p_region={"resolution": 0.02, "distance": 0.5},
                       depletion={"resolution": 0.005, "distance": 0.3},
                       n_region={"resolution": 0.02, "distance": 0.5},
                       box={"resolution": 0.1, "distance": 0.5},
                       clad={"resolution": 0.1, "distance": 0.5})
    return polys, resolutions


def solve_pn(N_donor, N_acceptor, wavelength_um=1.55, **kw):
    polys, res = build_pn_rib(N_donor=N_donor, N_acceptor=N_acceptor, **kw)
    mesh = from_meshio(mesh_from_OrderedDict(polys, res,
                                              default_resolution_max=0.15))
    basis0 = Basis(mesh, ElementTriP0())
    epsilon = basis0.zeros(dtype=complex)

    eps_p   = doped_silicon_eps(N_e=0,        N_h=N_acceptor, wavelength_um=wavelength_um)
    eps_n   = doped_silicon_eps(N_e=N_donor,  N_h=0,          wavelength_um=wavelength_um)
    eps_dep = n_silicon(wavelength_um)**2     # intrinsic in depletion region
    eps_ox  = n_sio2(wavelength_um)**2

    for name, e in [("p_region", eps_p),  ("n_region", eps_n),
                    ("depletion", eps_dep), ("box", eps_ox), ("clad", eps_ox)]:
        epsilon[basis0.get_dofs(elements=name)] = e

    modes = compute_modes(basis0, epsilon, wavelength=wavelength_um,
                          num_modes=1, order=1, n_guess=2.76)
    return modes


# ----- Loss extraction helper ---------------------------------------------
def loss_dB_per_cm(n_eff_complex, wavelength_um):
    """alpha [dB/cm] from Im(n_eff)."""
    alpha_per_um = 2*np.pi * abs(n_eff_complex.imag) / wavelength_um
    return alpha_per_um * 1e4 * 20 / np.log(10)


if __name__ == "__main__":
    # ---- Example 1: uniformly N-doped rib (no junction) -----------------
    # Easiest version: rebuild the rib geometry but assign doped_silicon_eps
    # to the "core" subdomain instead of pure silicon.
    print("="*70)
    print("Example: Uniformly N-doped rib  (N_e = 1e18 cm^-3)")
    print("="*70)
    from wg_sim import build_rib
    polys, res = build_rib()
    mesh = from_meshio(mesh_from_OrderedDict(polys, res,
                                              default_resolution_max=0.15))
    basis0 = Basis(mesh, ElementTriP0())
    eps = basis0.zeros(dtype=complex)
    WL = 1.55
    eps_si_doped = doped_silicon_eps(N_e=1e18, N_h=0, wavelength_um=WL)
    eps_ox = n_sio2(WL)**2
    for name, e in [("core", eps_si_doped), ("box", eps_ox), ("clad", eps_ox)]:
        eps[basis0.get_dofs(elements=name)] = e
    modes = compute_modes(basis0, eps, wavelength=WL, num_modes=1, n_guess=2.76)
    m = modes[0]
    print(f"  n_eff = {m.n_eff.real:.5f} + {m.n_eff.imag:.3e}j")
    print(f"  loss  = {loss_dB_per_cm(m.n_eff, WL):.3f} dB/cm")
    # compare to undoped:
    dn_e, dn_h, a_e, a_h = soref_bennett(1e18, 0, WL)
    print(f"  (Soref-Bennett: dn = {dn_e:.2e}, alpha_bulk = {a_e:.3f} /cm "
          f"= {a_e * 10 / np.log(10):.3f} dB/cm bulk material)")

    # ---- Example 2: PN-junction rib ------------------------------------
    print("\n" + "="*70)
    print("Example: PN junction rib  (N_donor=N_acceptor=1e18 cm^-3, "
          "depletion ~100 nm)")
    print("="*70)
    modes = solve_pn(N_donor=1e18, N_acceptor=1e18, wavelength_um=WL)
    m = modes[0]
    print(f"  n_eff = {m.n_eff.real:.5f} + {m.n_eff.imag:.3e}j")
    print(f"  loss  = {loss_dB_per_cm(m.n_eff, WL):.3f} dB/cm")

    # To get d(neff)/dV (modulator efficiency V*pi*L):  sweep the depletion
    # widths W_depl_n, W_depl_p (functions of applied reverse bias from a
    # device simulator like DEVSIM/Sentaurus) and re-solve.
