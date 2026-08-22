#!/usr/bin/env python3
"""The two new optical models, side by side.

**Left** — `fc_driven_laser`: a 5 Gb/s directly-modulated transmitter. One
`PULSE` source is the whole modulator; the optical power follows the drive
above threshold and sits on the floor below it, and a PIN + load resistor
turns that back into volts. The eye opens without an `fc_mzm` anywhere.

**Right** — `fc_facet`: return loss through a length of waveguide. The light
pays the propagation loss twice and the reflectance once, so a facet with a
known `reflectance` is a calibrated back-reflection to test a link against.

    python3 examples/photonic/dml_link_and_facet.py [--selftest]
"""

import argparse
import os
import sys

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))
import fairchild  # noqa: E402

# ── Transmitter ────────────────────────────────────────────────────────────
BIT_S = 200e-12           # 5 Gb/s
V_TH = 0.9                # laser threshold
V_HI = 2.4
SLOPE_W_V = 4e-3          # W/V above threshold
RESPONSIVITY = 0.8
R_LOAD = 200.0
C_PD = 150e-15            # receiver pole: tau = 30 ps, ~5.3 GHz with R_LOAD
STEP_S = 2e-12
N_BITS = 24
PATTERN = "010011000111010110100101"

# ── Facet sweep ────────────────────────────────────────────────────────────
WG_UM = 2000.0
WG_ALPHA_DB_CM = 2.0
A_IN = 0.031_622_776_601_683_79   # √(1 mW)


def dml_netlist() -> str:
    """PRBS-ish pattern as a PWL drive on the laser's own electrical port."""
    pts, t = [], 0.0
    edge = 0.15 * BIT_S
    for bit in PATTERN:
        v = V_HI if bit == "1" else 0.0
        pts += [f"{t + edge:.6e} {v:g}", f"{t + BIT_S:.6e} {v:g}"]
        t += BIT_S
    return f"""* directly-modulated 5 Gb/s link
.optical_port opt
.optical_port det_in
Vd  drv 0 PWL(0 0 {' '.join(pts)})
Xlas opt drv 0 fc_driven_laser slope_w_v={SLOPE_W_V} v_th={V_TH} r_in=1e12
Xwg opt det_in fc_waveguide L_um=500 n_g=4.2 alpha_dB_cm=2.0
Xpd det_in pa 0 fc_photodetector responsivity={RESPONSIVITY} r_shunt=1Meg i_dark_a=0
Rl  pa 0 {R_LOAD}
Cl  pa 0 {C_PD}
.tran {STEP_S} {N_BITS * BIT_S}
"""


def facet_netlist(reflectance: float) -> str:
    """Launch on the forward wires with plain sources — a laser here would
    also be driving the backward wires, to zero, and fight the reflection."""
    return f"""* facet return loss
.options enable_bidirectional=1
Xwg s_re_fw s_im_fw s_re_bw s_im_bw s_wl f_re_fw f_im_fw f_re_bw f_im_bw f_wl
+   fc_waveguide L_um={WG_UM} n_eff=2.445 n_g=4.2 alpha_dB_cm={WG_ALPHA_DB_CM}
Xf  f_re_fw f_im_fw f_re_bw f_im_bw f_wl fc_facet reflectance={reflectance}
V_re s_re_fw 0 DC {A_IN}
V_im s_im_fw 0 DC 0
V_wl s_wl 0 DC 1.55e-6
.op
"""


def run_dml():
    c = fairchild.Circuit()
    c.load_str(dml_netlist())
    r = c.run("tran", step=STEP_S, stop=N_BITS * BIT_S)
    return np.asarray(r.time()), np.asarray(r["V(pa)"])


def return_loss(reflectance: float) -> float:
    """Returned power at the launch port, W."""
    c = fairchild.Circuit()
    c.load_str(facet_netlist(reflectance))
    r = c.run("op")
    re, im = float(r["V(s_re_bw)"][0]), float(r["V(s_im_bw)"][0])
    return re * re + im * im


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--png", default="dml_link_and_facet.png")
    args = ap.parse_args()

    t, v_pa = run_dml()
    p_hi = SLOPE_W_V * (V_HI - V_TH)
    dc_hi = RESPONSIVITY * p_hi * R_LOAD * 10 ** (-2.0 * 0.05 / 10.0)  # 500 µm at 2 dB/cm

    # Sample mid-bit, after the first couple of bits have settled.
    centres = np.array([(k + 0.5) * BIT_S for k in range(2, N_BITS)])
    sampled = np.interp(centres, t, v_pa)
    ones = sampled[[b == "1" for b in PATTERN[2:]]]
    zeros = sampled[[b == "0" for b in PATTERN[2:]]]
    eye = ones.min() - zeros.max()
    er_db = 10 * np.log10(ones.mean() / max(zeros.mean(), 1e-18))

    rs = np.logspace(-4, 0, 25)
    p_ret = np.array([return_loss(r) for r in rs])
    one_way = 10 ** (-WG_ALPHA_DB_CM * (WG_UM * 1e-4) / 10.0)
    p_ret_analytic = A_IN**2 * one_way**2 * rs

    print(f"DML: high level {ones.mean() * 1e3:.1f} mV (DC limit {dc_hi * 1e3:.1f} mV), "
          f"eye {eye * 1e3:.1f} mV, ER {er_db:.1f} dB")
    print(f"facet: round-trip loss {-10 * np.log10(one_way**2):.2f} dB of waveguide "
          f"+ reflectance; max rel error "
          f"{np.max(np.abs(p_ret - p_ret_analytic) / p_ret_analytic):.2e}")

    if args.selftest:
        assert eye > 0, "the eye must be open"
        assert er_db > 10, f"received extinction {er_db:.1f} dB is too small"
        # Bandwidth-limited, so the peak must approach the DC level but not pass it.
        assert 0.5 * dc_hi < ones.mean() < 1.001 * dc_hi, ones.mean()
        assert np.max(np.abs(p_ret - p_ret_analytic) / p_ret_analytic) < 1e-9
        print("selftest OK — eye open, received extinction > 10 dB, "
              "return loss matches the closed form")
        return 0

    import matplotlib.pyplot as plt

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    for k in range(2, N_BITS):
        m = (t >= k * BIT_S) & (t < (k + 1) * BIT_S)
        ax.plot((t[m] - k * BIT_S) / BIT_S, v_pa[m] * 1e3, color="C0", alpha=0.5, lw=1)
    ax.axhline(dc_hi * 1e3, color="0.6", ls="--", lw=0.8, label="DC high level")
    ax.set_xlabel("bit phase")
    ax.set_ylabel("V(pa)  (mV)")
    ax.set_title(f"fc_driven_laser — {1e-9 / BIT_S:.0f} Gb/s eye, "
                 f"{eye * 1e3:.0f} mV open")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax2.loglog(rs, p_ret_analytic * 1e3, color="C2", lw=6, alpha=0.45,
               label="closed form  P·T²·R")
    ax2.loglog(rs, p_ret * 1e3, "k--", lw=1.4, label="fairchild .op")
    ax2.set_xlabel("facet reflectance (power)")
    ax2.set_ylabel("returned power at the launch port (mW)")
    ax2.set_title(f"fc_facet — {WG_UM:.0f} µm round trip")
    ax2.legend(fontsize=8)
    ax2.grid(True, which="both", alpha=0.3)

    out = os.path.join(os.path.dirname(__file__), args.png)
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
