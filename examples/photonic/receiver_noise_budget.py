#!/usr/bin/env python3
"""Direct-detection receiver noise budget: thermal vs shot vs RIN.

Sweeps the laser power of a laser → PIN → load-resistor link and reads the
output noise PSD out of `.noise`, then compares each term against the closed
form:

    S_V = ( 4kT/R_L  +  2q·I  +  RIN·I² ) · Z²        I = responsivity · P

The three terms scale as P⁰, P¹ and P². That ordering is the whole point of
the plot: turning the laser up buys SNR until RIN takes over, and after that
it buys nothing at all — the SNR ceiling is 1/(RIN·B), independent of power,
responsivity and load.

    python3 examples/photonic/receiver_noise_budget.py [--selftest]
"""

import argparse
import os
import sys

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))
import fairchild  # noqa: E402

Q = 1.602176634e-19
KB = 1.380649e-23
TEMP_K = 300.15

R_LOAD = 1.0e3
R_SHUNT = 1.0e9
RESPONSIVITY = 0.8
RIN_DB_HZ = -155.0
FREQ_HZ = 1.0e9
BANDWIDTH_HZ = 10.0e9  # for the SNR ceiling

Z = 1.0 / (1.0 / R_LOAD + 1.0 / R_SHUNT)
RIN = 10.0 ** (RIN_DB_HZ / 10.0)


def netlist(power_mw: float, rin: bool) -> str:
    rin_param = f" rin_db_hz={RIN_DB_HZ}" if rin else ""
    return f"""* direct-detection receiver
.optical_port opt
Xlas opt fc_cw_laser power_mw={power_mw}{rin_param}
Xpd  opt det bias fc_photodetector responsivity={RESPONSIVITY} r_shunt={R_SHUNT} i_dark_a=0
Rl   det bias {R_LOAD}
Vb   bias 0 DC 0
"""


def onoise(power_mw: float, rin: bool) -> float:
    """Output-referred noise PSD at FREQ_HZ, V²/Hz."""
    c = fairchild.Circuit()
    c.load_str(netlist(power_mw, rin))
    r = c.run("noise", out="det", out_neg="bias", src="Vb",
              fstart=FREQ_HZ, fstop=FREQ_HZ, points=1, variation="lin")
    return float(np.asarray(r["onoise"])[0])


def sweep(powers_mw):
    total = np.array([onoise(p, rin=True) for p in powers_mw])
    no_rin = np.array([onoise(p, rin=False) for p in powers_mw])
    i_ph = RESPONSIVITY * powers_mw * 1e-3
    thermal = np.full_like(powers_mw, 4.0 * KB * TEMP_K / R_LOAD * Z * Z)
    shot = 2.0 * Q * i_ph * Z * Z
    rin = RIN * i_ph**2 * Z * Z
    return total, no_rin, thermal, shot, rin, i_ph


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true",
                    help="assert the simulated budget matches the closed form")
    ap.add_argument("--png", default="receiver_noise_budget.png")
    args = ap.parse_args()

    dbm = np.linspace(-30.0, 13.0, 40)
    powers_mw = 10.0 ** (dbm / 10.0)
    total, no_rin, thermal, shot, rin, i_ph = sweep(powers_mw)

    # Every term isolated: the sim's RIN contribution is (total − no_rin), its
    # shot contribution is (no_rin − thermal-only), and thermal is the P→0
    # limit. Comparing sums would let two errors cancel.
    err_total = np.max(np.abs(total - (thermal + shot + rin)) / (thermal + shot + rin))
    err_rin = np.max(np.abs((total - no_rin) - rin) / rin)
    err_shot = np.max(np.abs((no_rin - thermal) - shot) / shot)

    snr = (i_ph * Z) ** 2 / (total * BANDWIDTH_HZ)
    ceiling = 1.0 / (RIN * BANDWIDTH_HZ)
    # Crossover powers: where each pair of terms is equal.
    p_shot_thermal = 4.0 * KB * TEMP_K / R_LOAD / (2.0 * Q) / RESPONSIVITY * 1e3
    p_rin_shot = 2.0 * Q / RIN / RESPONSIVITY * 1e3

    print(f"max rel error, total budget : {err_total:.2e}")
    print(f"max rel error, RIN term     : {err_rin:.2e}")
    print(f"max rel error, shot term    : {err_shot:.2e}")
    print(f"shot = thermal at {10 * np.log10(p_shot_thermal):+.1f} dBm")
    print(f"RIN  = shot    at {10 * np.log10(p_rin_shot):+.1f} dBm")
    print(f"SNR ceiling in {BANDWIDTH_HZ / 1e9:.0f} GHz: "
          f"{10 * np.log10(ceiling):.1f} dB (reached {10 * np.log10(snr[-1]):.1f} dB)")

    if args.selftest:
        assert err_total < 1e-6, err_total
        assert err_rin < 1e-6, err_rin
        assert err_shot < 1e-6, err_shot
        assert snr[-1] < ceiling, "SNR cannot exceed the RIN ceiling"
        assert snr[-1] > 0.7 * ceiling, "sweep should reach the RIN-limited regime"
        print("selftest OK — simulated budget matches the closed form term by term")
        return 0

    import matplotlib.pyplot as plt

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    ax.loglog(powers_mw, thermal, "--", label=f"thermal  4kT/R  ({R_LOAD / 1e3:g} kΩ)")
    ax.loglog(powers_mw, shot, "--", label="shot  2qI")
    ax.loglog(powers_mw, rin, "--", label=f"RIN  ({RIN_DB_HZ:g} dB/Hz)")
    ax.loglog(powers_mw, total, "k", lw=2, label="fairchild .noise")
    for p, tag in ((p_shot_thermal, "shot=thermal"), (p_rin_shot, "RIN=shot")):
        ax.axvline(p, color="0.7", lw=0.8)
        ax.annotate(tag, (p, total[0]), rotation=90, fontsize=7,
                    ha="right", va="bottom", color="0.4")
    ax.set_xlabel("optical power (mW)")
    ax.set_ylabel(f"output noise PSD at {FREQ_HZ / 1e9:g} GHz (V²/Hz)")
    ax.set_title("Receiver noise budget")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    ax2.semilogx(powers_mw, 10 * np.log10(snr), "k", lw=2, label="SNR")
    ax2.axhline(10 * np.log10(ceiling), color="C2", ls="--",
                label=f"RIN ceiling  1/(RIN·B) = {10 * np.log10(ceiling):.0f} dB")
    ax2.set_xlabel("optical power (mW)")
    ax2.set_ylabel(f"SNR in {BANDWIDTH_HZ / 1e9:g} GHz (dB)")
    ax2.set_title("More power stops helping")
    ax2.legend(fontsize=8)
    ax2.grid(True, which="both", alpha=0.3)

    out = os.path.join(os.path.dirname(__file__), args.png)
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
