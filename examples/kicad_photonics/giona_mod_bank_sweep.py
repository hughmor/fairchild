#%%
#!/usr/bin/env python3
"""Wavelength sweep on giona_mod_bank.sp — 8-MRM cascade spectrum.

Sweeps the CW laser wavelength across a small band, solves DC operating
point at each wavelength, and plots V(V_THRU) and V(V_DROP) — proportional
to thru-port and drop-port optical power respectively.

Usage:  python giona_mod_bank_sweep.py
"""
import os
import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
NETLIST = os.path.join(HERE, "giona_mod_bank.sp")
FAIRCHILD = os.path.abspath(os.path.join(HERE, "..", "..", "target", "release", "fairchild"))


def run_one(wavelength_nm: float, v_mzm: float = 0.0) -> dict:
    """Run a single DC OP at the given laser wavelength.

    `v_mzm` overrides the MZM modulation voltage so the bank can see real
    optical power (V_mzm=0 → MZM bright; V_mzm=Vπ=1 → MZM dark, default in
    the netlist gives near-zero light into the bank).
    """
    cmd = [
        FAIRCHILD,
        "--file", NETLIST,
        "--param", f"XCWL1.wavelength_nm={wavelength_nm}",
        "--param", f"V1.dc={v_mzm}",
        "--probe", "V(v_thru),V(v_drop)",
    ]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    # Last non-empty line is the data row.
    lines = [ln for ln in out.strip().split("\n") if ln.strip()]
    header = lines[-2].split(",")
    # First column is the analysis label ("dc_op"); skip it.
    values = lines[-1].split(",")[1:]
    cols = header[1:]
    return dict(zip(cols, [float(v) for v in values]))


def main():
    # Sweep over ~3 nm around 1550 nm.  With L_ring ≈ 50 µm and n_g = 4.2,
    # FSR(λ) ≈ λ²/(n_g·L) ≈ 11.4 nm, so we sample within one FSR.
    wavelengths = 1550.0 + np.linspace(-6, 6, 1000)  # 2.5 pm steps
    v_thru = np.zeros_like(wavelengths)
    v_drop = np.zeros_like(wavelengths)

    for i, wl in enumerate(wavelengths):
        try:
            r = run_one(wl, v_mzm=0.0)
            v_thru[i] = r["V(v_thru)"]
            v_drop[i] = r["V(v_drop)"]
            if i % 20 == 0:
                print(f"  λ={wl:.3f} nm  V_thru={v_thru[i]*1e3:.3f} mV  "
                      f"V_drop={v_drop[i]*1e3:.3f} mV")
        except subprocess.CalledProcessError as e:
            print(f"  λ={wl:.3f} nm  FAILED: {e.stderr.splitlines()[-1]}")
            v_thru[i] = np.nan
            v_drop[i] = np.nan

    # Convert V_pd → optical power (V_pd = I_ph·R_load = R_pd·P_opt·R_load).
    # With R_load = 1 kΩ and responsivity = 0.7 A/W, 1 V = 1.43 mW.
    p_thru_mw = -v_thru / (0.7 * 1.0)     # mW (R_load already in kΩ → V/I → mA → mW @ 1V cathode)
    p_drop_mw = -v_drop / (0.7 * 1.0)
    p_thru_dbm = 10.0 * np.log10(np.clip(p_thru_mw, 1e-9, None))
    p_drop_dbm = 10.0 * np.log10(np.clip(p_drop_mw, 1e-9, None))

    fig, (ax_lin, ax_dB) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    ax_lin.plot(wavelengths, p_thru_mw, label="Thru-port",  color="tab:blue")
    ax_lin.plot(wavelengths, p_drop_mw, label="Drop-port",  color="tab:orange")
    ax_lin.set_ylabel("PD output [mW equiv.]")
    ax_lin.legend(loc="best")
    ax_lin.grid(alpha=0.3)
    ax_lin.set_title("8-MRM cascade — DC wavelength sweep")

    ax_dB.plot(wavelengths, p_thru_dbm, color="tab:blue")
    ax_dB.plot(wavelengths, p_drop_dbm, color="tab:orange")
    ax_dB.set_ylabel("PD output [dB]")
    ax_dB.set_xlabel("Wavelength [nm]")
    ax_dB.grid(alpha=0.3)
    ax_dB.set_ylim(-40, 5)

    out_png = os.path.join(HERE, "giona_mod_bank_spectrum.png")
    plt.tight_layout()
    plt.savefig(out_png, dpi=120)
    print(f"Saved {out_png}")


if __name__ == "__main__":
    main()

# %%
