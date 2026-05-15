#!/usr/bin/env python3
"""
Native MRR wavelength sweep — characterise the ring's optical transfer
function by sweeping the laser wavelength and recording the through-port
power at each step.  This is the photonic equivalent of an .ac sweep:
we run a DC OP at every wavelength point with V_pn = 0 (on-resonance
condition) and at V_pn = Vπ (off-resonance) and plot both.

Reveals the ring's free spectral range (FSR), extinction ratio, and how
the resonance migrates with applied PN bias — a one-shot view of an
electro-optic modulator's static response.

Uses only B-phase native devices; companion to native_mrr_modulator.{sp,py}.
"""
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np

try:
    import fairchild
except ImportError as e:
    sys.exit(f"fairchild Python package not installed: {e}\n"
             "Build with: cd crates/fairchild-py && maturin develop --release")

HERE = pathlib.Path(__file__).resolve().parent


def mrr_netlist(wavelength_nm: float, v_pn: float) -> str:
    """Build the same circuit as native_mrr_modulator.sp, parameterised.

    SPICE convention: the first logical line is the title (consumed
    verbatim, even if it looks like a directive).  The `* …` comment
    below is mandatory — without it, `.optical_port laser_out` would
    be eaten as the title and the bundle wouldn't expand.
    """
    return f"""\
* Native MRR wavelength sweep (parameterised), λ={wavelength_nm} nm, V_pn={v_pn} V
.optical_port laser_out
.optical_port wg1_out
.optical_port dc_b
.optical_port dc_c
.optical_port pn_in
.optical_port pd_in

Xlaser laser_out fc_cw_laser power_mW=1.0 wavelength_nm={wavelength_nm}

Xwg1 laser_out wg1_out fc_waveguide L_um=50 n_g=4.2 alpha_dB_cm=2.0 wavelength_nm=1550

Xdc wg1_out dc_b dc_c pn_in fc_dcoupler kappa_L=0.336

Xpn pn_in dc_b vmod 0 fc_pn_ps L_um=500 V_pi_L=2e-3 g_pn=1e-3 alpha_dB_cm=10

Xwg2 dc_c pd_in fc_waveguide L_um=50 n_g=4.2 alpha_dB_cm=2.0 wavelength_nm=1550

Xpd pd_in pd_anode 0 fc_photodetector responsivity=0.8 i_dark_a=1e-9 r_shunt=1Meg

Vbias bias 0 DC 1.0
Rload pd_anode bias 1k
Vmod vmod 0 DC {v_pn}
.op
.end
"""


def sweep(wavelengths_nm, v_pn):
    """Return (Re, Im, V_pd) arrays — one element per swept wavelength."""
    re  = np.zeros_like(wavelengths_nm)
    im  = np.zeros_like(wavelengths_nm)
    v_pd = np.zeros_like(wavelengths_nm)
    c = fairchild.Circuit()
    for i, wl in enumerate(wavelengths_nm):
        c.load_str(mrr_netlist(float(wl), v_pn))
        r = c.run("op")
        re[i]  = r["V(pd_in_re_0)"][0]
        im[i]  = r["V(pd_in_im_0)"][0]
        v_pd[i] = r["V(pd_anode)"][0]
    return re, im, v_pd


def main():
    # Sweep wavelength across a small range that covers one full FSR
    # of a ~500 µm silicon ring (FSR ≈ λ²/(n_g·L_ring) ≈ 1.14 nm for
    # n_g=4.2, L=500 µm at λ=1.55 µm).
    wavelengths = np.linspace(1549.0, 1551.0, 200)

    print("Sweeping V_pn = 0 ...")
    re0, im0, pd0 = sweep(wavelengths, v_pn=0.0)
    print("Sweeping V_pn = 4 V (≈ Vπ) ...")
    re4, im4, pd4 = sweep(wavelengths, v_pn=4.0)

    p0 = re0 ** 2 + im0 ** 2
    p4 = re4 ** 2 + im4 ** 2

    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 6))

    axes[0].plot(wavelengths, p0 * 1e3, label="V(pn) = 0",  color="tab:blue")
    axes[0].plot(wavelengths, p4 * 1e3, label="V(pn) = 4 V", color="tab:orange")
    axes[0].set_ylabel("Optical power at PD (mW)")
    axes[0].set_title("Native MRR wavelength sweep — through-port transmission")
    axes[0].grid(alpha=0.3)
    axes[0].legend(loc="lower left")

    axes[1].plot(wavelengths, pd0, label="V(pn) = 0",  color="tab:blue")
    axes[1].plot(wavelengths, pd4, label="V(pn) = 4 V", color="tab:orange")
    axes[1].set_xlabel("wavelength (nm)")
    axes[1].set_ylabel("V(pd_anode) (V)")
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    out = HERE / "native_mrr_wavelength_sweep.png"
    fig.savefig(out, dpi=120)
    print(f"wrote {out}")
    plt.show()


if __name__ == "__main__":
    main()
