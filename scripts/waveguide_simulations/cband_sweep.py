"""C-band sweeps: n_eff(lambda) for strip, rib (straight), and rib (R=8 um).
Computes n_g from finite differences at the center wavelength."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from wg_sim import build_strip, build_rib, solve_neff, n_silicon, n_sio2


def sweep(builder, wls, label, radius=np.inf, n_guess=2.4):
    neffs = []
    for wl in wls:
        polys, res = builder()
        modes = solve_neff(polys, res, wl, num_modes=1, radius=radius,
                           n_guess=n_guess)
        n_eff = complex(modes[0].n_eff)
        neffs.append(n_eff)
        print(f"  {label} @ {wl*1000:.1f} nm:  n_eff = "
              f"{n_eff.real:.5f}   loss(imag) = {n_eff.imag:.2e}")
    return np.array(neffs)


def compute_ng(wls, neffs_real):
    """n_g = n_eff - lambda * d(n_eff)/d(lambda) via central differences."""
    n_eff_c = neffs_real[len(wls)//2]
    wl_c = wls[len(wls)//2]
    # use 1st-order polynomial fit for robust dn/dlambda at the center
    p = np.polyfit(wls, neffs_real, 1)   # p[0] = dn/dl
    dn_dl = p[0]
    n_g = n_eff_c - wl_c * dn_dl
    return n_g, dn_dl


# C-band: 1530-1565 nm, 8 points
wls = np.linspace(1.530, 1.565, 8)

print("="*70)
print("Strip 500x220 nm  (TE, straight)")
print("="*70)
neff_strip = sweep(build_strip, wls, "strip", n_guess=2.44)

print("\n" + "="*70)
print("Rib (500x220 core + 4500x90 slab) (TE, straight)")
print("="*70)
neff_rib = sweep(build_rib, wls, "rib", n_guess=2.76)

print("\n" + "="*70)
print("Rib (500x220 core + 4500x90 slab) (TE, R = 8 um)")
print("="*70)
neff_rib_bent = sweep(build_rib, wls, "rib_R8", radius=8.0, n_guess=2.76)

# n_g
ng_strip, dn_dl_strip = compute_ng(wls, neff_strip.real)
ng_rib, dn_dl_rib = compute_ng(wls, neff_rib.real)
ng_rib_bent, dn_dl_rib_bent = compute_ng(wls, neff_rib_bent.real)

# bulk Si group index for reference
ng_si_bulk = n_silicon(1.55) - 1.55 * np.polyfit(wls, n_silicon(wls), 1)[0]

print("\n" + "="*70)
print("Summary at 1550 nm (interpolated to center of sweep)")
print("="*70)
print(f"  Strip:        n_eff = {np.interp(1.5475, wls, neff_strip.real):.4f}  "
      f"n_g = {ng_strip:.4f}  dn/dlambda = {dn_dl_strip:.4f}/um")
print(f"  Rib straight: n_eff = {np.interp(1.5475, wls, neff_rib.real):.4f}  "
      f"n_g = {ng_rib:.4f}  dn/dlambda = {dn_dl_rib:.4f}/um")
print(f"  Rib R=8um:    n_eff = {np.interp(1.5475, wls, neff_rib_bent.real):.4f}  "
      f"n_g = {ng_rib_bent:.4f}  dn/dlambda = {dn_dl_rib_bent:.4f}/um")
print(f"  (bulk Si n_g for reference: {ng_si_bulk:.4f})")

# also report imaginary part = radiation loss for the bend
print("\nBend losses (rib R=8um):")
for wl, n in zip(wls, neff_rib_bent):
    # alpha [dB/cm] = 20 / ln(10) * (2*pi/lambda) * Im(n_eff)  with lambda, length in same units
    # Convert: 2*pi*Im(neff)/wavelength_um -> [1/um], then *1e4 -> [1/cm], then *(20/ln10) -> dB/cm
    alpha_per_um = 2*np.pi*abs(n.imag)/wl   # nepers per um
    alpha_db_cm = alpha_per_um * 1e4 * 20/np.log(10)
    print(f"   {wl*1000:.1f} nm: Im(n_eff) = {n.imag:.2e}  "
          f"-> loss = {alpha_db_cm:.3f} dB/cm  (radiation only)")

# save data to CSV
import csv
with open("/home/claude/cband_sweep.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["wavelength_um", "neff_strip", "neff_rib_straight",
                "neff_rib_R8um", "loss_imag_rib_R8um"])
    for i, wl in enumerate(wls):
        w.writerow([wl, neff_strip[i].real, neff_rib[i].real,
                    neff_rib_bent[i].real, neff_rib_bent[i].imag])
print("\nData saved to cband_sweep.csv")

# Plot
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(wls*1000, neff_strip.real, "o-", label=f"Strip (n_g={ng_strip:.3f})")
ax[0].plot(wls*1000, neff_rib.real, "s-", label=f"Rib straight (n_g={ng_rib:.3f})")
ax[0].plot(wls*1000, neff_rib_bent.real, "^--",
           label=f"Rib R=8um (n_g={ng_rib_bent:.3f})")
ax[0].set_xlabel("Wavelength (nm)")
ax[0].set_ylabel("n_eff (real)")
ax[0].set_title("n_eff across the C band")
ax[0].legend()
ax[0].grid(alpha=0.3)

ax[1].semilogy(wls*1000, np.abs(neff_rib_bent.imag), "^-", color="C2",
               label="Rib R=8um")
ax[1].set_xlabel("Wavelength (nm)")
ax[1].set_ylabel("|Im(n_eff)|  (bend radiation)")
ax[1].set_title("Radiation loss vs wavelength")
ax[1].legend()
ax[1].grid(alpha=0.3, which="both")

plt.tight_layout()
plt.savefig("/home/claude/cband_sweep.png", dpi=130)
print("Plot saved to cband_sweep.png")
