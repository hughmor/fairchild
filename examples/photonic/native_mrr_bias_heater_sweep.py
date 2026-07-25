#!/usr/bin/env python3
"""
Add-drop micro-ring modulator — static bias and heater characterisation.

The two sweeps you actually run on a fabricated MRM, on a model whose defaults
came from fitting a real silicon device (see `experiments/giona/` for the fits
and `experiments/giona/giona_pn_th_ps.inc` for the card's provenance):

  1. **PN bias, −1 → +1 V** at heater off. Reverse bias depletes the junction
     (Δn ∝ V, a small *red* shift); forward bias injects carriers, and there the
     index follows the *current*, not the voltage (Δn = −dn_di·I) — so the shift
     is exponential in V and *blue*, an order of magnitude larger by +1 V.
  2. **Heater current, 0 → 1 mA** at zero bias. Thermo-optic, so the shift is
     linear in dissipated power P = I²R and therefore *quadratic* in current.

Topology (add-drop: two couplers, the ring split into two heated PN arms whose
heaters are wired in series, so one current drives both):

      IN ──► CPL2 ─────────────────────────► THRU
              │  ▲
            PS2  PS1        (each = ½ ring: PN junction + heater)
              ▼  │
      DROP ◄── CPL1 ◄─── ADD

Both ports are plotted: the neuron this device came from reads through and drop
into a balanced photodetector pair, and the drop port carries the deeper,
higher-contrast resonance.

Run:      .venv/bin/python examples/photonic/native_mrr_bias_heater_sweep.py
Selftest: same, with --selftest (asserts the physics, no plotting)
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np

try:
    import fairchild as fc
except ImportError as e:
    sys.exit(f"fairchild Python package not installed: {e}\n"
             "Build with: maturin develop --release -m crates/fairchild-py/Cargo.toml")

HERE = pathlib.Path(__file__).resolve().parent

# ── device card ──────────────────────────────────────────────────────────────
# LEVEL=4 dispatches the fc_pn_th_ps family to its full model: depletion
# (dn_dv, C_j(V)) + carrier injection driven by the diode current (dn_di, da_di)
# + TPA/self-heating. Values are the fitted defaults for a real 25 µm half-ring
# arm; l_m/r_heater/p_pi_th are therefore PER ARM (two arms = one ring).
CARD = """\
* Add-drop MRM bias + heater sweep
.model mrm_ps fc_pn_th_ps LEVEL=4
+ l_m=2.5378e-5 n_g=4.2 n_eff=2.2810 pin_at_ref=0 alpha_db_cm=10.7
+ r_heater=184.4 p_pi_th=26.4e-3
+ dn_dv=-3.62e-5 da_dv=3.29e-4 c_j0=1.375e-13 v_bi=0.917 m_j=0.5
+ i_sat=5.099e-8 n_diode=5.0 r_series=0 tau_carrier=10e-9
+ dn_dv_inj=0 da_dv_inj=0 dn_di=3.99 da_di=4.63e6
+ beta_tpa=7.9e-12 a_eff_m2=1.257e-13 dn_dt=1.86e-4 r_th=0
"""

# fc_pn_th_ps LEVEL=4 terminals: [in, out, anode, cathode, heat_p, heat_n].
# anode=vpn, cathode=GND ⇒ vpn > 0 forward-biases (injection), vpn < 0 depletes.
NETLIST = """\
.optical_port in_opt
.optical_port thru_opt
.optical_port drop_opt
.optical_port add_opt
.optical_port ring_a
.optical_port ring_b
.optical_port ring_c
.optical_port ring_d
Xlaser in_opt fc_cw_laser power_mW=1.0 wavelength_nm=1546.5
XCPL2  in_opt ring_a thru_opt ring_b fc_dcoupler kappa_L=0.183
XCPL1  ring_c add_opt ring_d drop_opt fc_dcoupler kappa_L=0.183
XPS1   ring_d ring_a vpn 0 hbias hmid mrm_ps
XPS2   ring_b ring_c vpn 0 hmid 0     mrm_ps
Vpn  vpn   0 DC 0
Iheat hbias 0 DC 0
.op
"""

R_HEATER_ARM = 184.4          # Ω per arm; series pair ⇒ 2× this total
P_PI_ARM = 26.4e-3            # W per arm for π of phase


class Ring:
    """One loaded circuit; every sweep point is a set_param + .op re-solve."""

    def __init__(self):
        self.ckt = fc.Circuit()
        self.ckt.load_str(CARD + NETLIST)

    def op(self, wl_nm: float, v_pn: float, i_heat: float):
        self.ckt.set_param("Xlaser", "wavelength_nm", float(wl_nm))
        self.ckt.set_param("Vpn", "dc", float(v_pn))
        self.ckt.set_param("Iheat", "dc", float(i_heat))
        return self.ckt.run("op")

    def spectra(self, wl_nm, v_pn=0.0, i_heat=0.0):
        """Through- and drop-port transmission (dB, relative to input power)."""
        thru = np.empty(len(wl_nm))
        drop = np.empty(len(wl_nm))
        for i, wl in enumerate(wl_nm):
            r = self.op(wl, v_pn, i_heat)
            p_in = _power(r, "in_opt")
            thru[i] = 10 * np.log10(max(_power(r, "thru_opt"), 1e-30) / p_in)
            drop[i] = 10 * np.log10(max(_power(r, "drop_opt"), 1e-30) / p_in)
        return thru, drop

    def junction_current(self, v_pn: float) -> float:
        """|I| through the PN pair at this bias (both arms, off resonance)."""
        return abs(float(self.op(1540.0, v_pn, 0.0)["I(vpn)"][0]))


def _power(r, port: str) -> float:
    """Optical power on a bundle port: |E|² from its expanded re/im nodes."""
    re = float(r[f"V({port}_re_0)"][0])
    im = float(r[f"V({port}_im_0)"][0])
    return max(re * re + im * im, 1e-30)


def track_notch(wl, t_dB) -> float:
    """Resonance wavelength to sub-sample precision: parabola through the
    minimum and its two neighbours. Needed because the depletion shift over
    the whole ±1 V range is only ~10 pm — coarser than the notch itself."""
    i = int(np.argmin(t_dB))
    if i == 0 or i == len(t_dB) - 1:
        return float(wl[i])
    y0, y1, y2 = t_dB[i - 1], t_dB[i], t_dB[i + 1]
    denom = y0 - 2 * y1 + y2
    shift = 0.0 if denom == 0 else 0.5 * (y0 - y2) / denom
    return float(wl[i] + shift * (wl[1] - wl[0]))


def locate_resonance(ring: Ring, lo=1540.0, hi=1552.0, n=241) -> float:
    """Find the unbiased resonance by coarse scan, so the fine sweep windows
    self-centre instead of hard-coding a wavelength that fab variation moves."""
    wl = np.linspace(lo, hi, n)
    thru, _ = ring.spectra(wl)
    return track_notch(wl, thru)


# ── sweeps ───────────────────────────────────────────────────────────────────
V_BIAS = np.array([-1.0, -0.75, -0.5, -0.25, 0.0, 0.4, 0.6, 0.8, 1.0])
I_HEAT = np.linspace(0.0, 1e-3, 6)
WINDOW_NM, WINDOW_PTS = 0.5, 301


def run_sweeps(verbose=True):
    ring = Ring()
    # Forward bias blue-shifts hard, so hang the window asymmetrically below λ0.
    lam0 = locate_resonance(ring)
    wl = np.linspace(lam0 - WINDOW_NM, lam0 + 0.15, WINDOW_PTS)
    # Re-read λ0 on the fine grid: the coarse scan's 50 pm sampling leaves a
    # ~0.1 pm bias, which is not negligible against a 13 pm depletion shift.
    lam0 = track_notch(wl, ring.spectra(wl)[0])
    if verbose:
        print(f"unbiased resonance: {lam0:.4f} nm")

    bias = {"v": V_BIAS, "wl": wl, "thru": [], "drop": [], "dlam": [], "i_pn": []}
    for v in V_BIAS:
        thru, drop = ring.spectra(wl, v_pn=float(v))
        bias["thru"].append(thru)
        bias["drop"].append(drop)
        bias["dlam"].append((track_notch(wl, thru) - lam0) * 1e3)   # pm
        bias["i_pn"].append(ring.junction_current(float(v)))
        if verbose:
            print(f"  V={v:+.2f} V  Δλ={bias['dlam'][-1]:+8.2f} pm  "
                  f"I={bias['i_pn'][-1] * 1e6:9.3f} µA")

    heat = {"i": I_HEAT, "wl": wl, "thru": [], "drop": [], "dlam": []}
    for i_h in I_HEAT:
        thru, drop = ring.spectra(wl, i_heat=float(i_h))
        heat["thru"].append(thru)
        heat["drop"].append(drop)
        heat["dlam"].append((track_notch(wl, thru) - lam0) * 1e3)
        if verbose:
            print(f"  I={i_h * 1e3:.2f} mA  P={i_h ** 2 * 2 * R_HEATER_ARM * 1e3:6.3f} mW"
                  f"  Δλ={heat['dlam'][-1]:+8.2f} pm")

    for d in (bias, heat):
        d["dlam"] = np.array(d["dlam"])
    bias["i_pn"] = np.array(bias["i_pn"])
    return lam0, bias, heat


# ── plotting ─────────────────────────────────────────────────────────────────
def plot(lam0, bias, heat, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(3, 2, figsize=(11, 10.5))

    # Rows 0/1: the spectra families, through port then drop port. Separate axes
    # per port — on a shared one the drop peak lands on the through baseline.
    sweeps = [
        (0, bias, "v", "coolwarm", lambda x: f"{x:+.2f} V", "PN bias (heater off)"),
        (1, heat, "i", "inferno", lambda x: f"{x * 1e3:.1f} mA",
         "heater current (0 V bias)"),
    ]
    for col, d, key, cname, fmt, what in sweeps:
        cmap = plt.get_cmap(cname)
        n = len(d[key])
        for row, port in ((0, "thru"), (1, "drop")):
            a = ax[row, col]
            for k, val in enumerate(d[key]):
                frac = k / (n - 1)
                a.plot(d["wl"], d[port][k], lw=1.1, label=fmt(val),
                       color=cmap(frac if cname == "coolwarm"
                                  else 0.15 + 0.7 * frac))
            a.set_title(f"{port.capitalize()}-port transmission vs {what}")
            a.set_xlabel("wavelength (nm)")
            a.set_ylabel("transmission (dB)")
            a.grid(alpha=0.3)
            a.legend(fontsize=7, ncol=2 if n > 6 else 1,
                     loc="lower left" if port == "thru" else "upper left")

    ax[2, 0].plot(bias["v"], bias["dlam"], "o-", color="tab:red")
    ax[2, 0].axhline(0, color="k", lw=0.6)
    ax[2, 0].axvline(0, color="k", lw=0.6)
    ax[2, 0].set_xlabel("PN bias (V)")
    ax[2, 0].set_ylabel("Δλ_res (pm)")
    ax[2, 0].set_title("Resonance shift vs bias — depletion (V) vs injection (I)")
    ax[2, 0].grid(alpha=0.3)
    axi = ax[2, 0].twinx()
    axi.semilogy(bias["v"], np.maximum(bias["i_pn"] * 1e6, 1e-6), "s:",
                 color="tab:gray", ms=3.5)
    axi.set_ylabel("|I_pn| (µA)", color="tab:gray")
    axi.tick_params(axis="y", colors="tab:gray")

    p_mW = heat["i"] ** 2 * 2 * R_HEATER_ARM * 1e3
    ax[2, 1].plot(heat["i"] * 1e3, heat["dlam"], "o-", color="tab:orange",
                  label="simulated")
    if heat["dlam"][-1] != 0:
        ax[2, 1].plot(heat["i"] * 1e3, heat["dlam"][-1] * p_mW / p_mW[-1],
                      "k--", lw=0.8, label="∝ P = I²R")
    ax[2, 1].set_xlabel("heater current (mA)")
    ax[2, 1].set_ylabel("Δλ_res (pm)")
    ax[2, 1].set_title("Resonance shift vs heater current (quadratic in I)")
    ax[2, 1].legend(fontsize=8)
    ax[2, 1].grid(alpha=0.3)

    fig.suptitle(f"Add-drop MRM static characterisation — λ_res(0 V, 0 mA) "
                 f"= {lam0:.3f} nm", fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    print(f"wrote {out}")


# ── selftest ─────────────────────────────────────────────────────────────────
def selftest(lam0, bias, heat) -> int:
    v, dl, ip = bias["v"], bias["dlam"], bias["i_pn"]
    rev, fwd = v < 0, v > 0

    # Depletion: reverse bias red-shifts (dn_dv < 0 with v_junc < 0 ⇒ Δn > 0),
    # monotonically in |V|, and stays small.
    assert np.all(dl[rev] > 0), f"reverse bias should red-shift, got {dl[rev]}"
    assert np.all(np.diff(dl[rev][::-1]) > 0), f"non-monotonic depletion: {dl[rev]}"
    assert dl[rev].max() < 50, f"depletion shift implausibly large: {dl[rev].max()} pm"

    # Injection: forward bias blue-shifts, and much harder — carrier density
    # follows the diode current, which is exponential in V.
    assert np.all(dl[fwd] < 0), f"forward bias should blue-shift, got {dl[fwd]}"
    assert abs(dl[fwd].min()) > 5 * dl[rev].max(), (
        f"injection ({dl[fwd].min():.1f} pm) should dominate depletion "
        f"({dl[rev].max():.1f} pm)")
    # The dn_dv term applies at any bias, so forward Δλ = (linear in V) +
    # (injection). Subtract the reverse-bias slope and what's left must track
    # the *current* — that is the dn_di claim, and it is what makes the forward
    # branch exponential in V. Only well-resolved points (>3 pm) qualify.
    slope = np.polyfit(v[rev], dl[rev], 1)[0]          # pm/V, from depletion only
    excess = dl[fwd] - slope * v[fwd]
    big = np.abs(excess) > 3.0
    ratio = excess[big] / ip[fwd][big]
    assert np.all(excess < 0), f"injection excess should blue-shift: {excess}"
    assert ratio.max() / ratio.min() < 1.5, (
        f"injection Δλ should track I; Δλ_exc/I spread "
        f"{ratio.max() / ratio.min():.2f} over {ip[fwd][big] * 1e6} µA")

    # Thermal: shift is linear in P = I²R ⇒ quadratic in I, and same sign as
    # dn_dt > 0 (heating raises n ⇒ red shift).
    p = heat["i"] ** 2
    dlh = heat["dlam"]
    assert np.all(dlh >= 0), f"heating should red-shift, got {dlh}"
    assert np.all(np.diff(dlh) > -1e-9), f"non-monotonic thermal tuning: {dlh}"
    nz = p > 0
    k = dlh[nz] / p[nz]
    assert k.max() / k.min() < 1.05, f"Δλ not ∝ I²: k spread {k.max()/k.min():.3f}"

    print(f"selftest OK — depletion +{dl[rev].max():.1f} pm @ -1 V, "
          f"injection {dl[fwd].min():.1f} pm @ +1 V, "
          f"thermal +{dlh[-1]:.1f} pm @ 1 mA "
          f"(π needs {np.sqrt(P_PI_ARM / 2 / R_HEATER_ARM) * 1e3:.1f} mA)")
    return 0


def main() -> int:
    st = "--selftest" in sys.argv
    lam0, bias, heat = run_sweeps(verbose=not st)
    if st:
        return selftest(lam0, bias, heat)
    selftest(lam0, bias, heat)
    plot(lam0, bias, heat, HERE / "native_mrr_bias_heater_sweep.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
