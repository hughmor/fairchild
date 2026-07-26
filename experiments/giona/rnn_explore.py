#!/usr/bin/env python3
"""rnn_explore.py — calibrate and program the giona 8-neuron WDM recurrent net.

Works on `netlists/giona_rnn_perfectW.sp` (hand-written; this script never
rewrites it — it drives it through `set_param`, and prints the `.param` edits it
recommends).

Stages, each runnable on its own:

  --radii       Solve each ring's radius so ring i is resonant with channel i.
  --activation  Characterise the ring activation vs junction bias, on both the
                peak port (ReLU-like) and the notch port (sigmoid-like), and
                show what the heater does to it.
  --bias        Derive the DC bias each neuron needs to sit at a chosen
                operating point, in the recurrent case.
  --osc         Two-neuron oscillator: W = [[1,-1],[1,1]], transient.

Run:  .venv/bin/python experiments/giona/rnn_explore.py --radii --activation
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

try:
    import fairchild as fc
except ImportError as e:
    sys.exit(f"fairchild not installed: {e}")

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
NETLIST = HERE / "netlists" / "giona_rnn_perfectW.sp"
MRM_CELL = REPO / "examples" / "photonic" / "pcells" / "mrm.sp"
RESULTS = HERE / "results"

N = 8
LAMBDAS_NM = [1546.12, 1546.92, 1547.72, 1548.51,
              1549.32, 1550.12, 1550.92, 1551.72]
N_EFF = 2.302216932          # shared by every ring in the hand-written deck
R0 = 8e-6

# The chip's optical path from a ring's bus output to a weight block input:
# 2 mm waveguide, a 1:1 tap, then a 1:8 log tree.
TAP_TREE_ATT = 0.5 / 8.0
RESPONSIVITY = 0.8


# ── single-ring probe ────────────────────────────────────────────────────────
def probe_deck(radius: float, wl_nm: float, v_c: float, i_ht: float,
               feed: str) -> str:
    """One ring, one wavelength. `feed` is 'in' (bus) or 'add' (drop bus).

    Hugh's deck wires pn_a=0 and pn_c=mod_cathode, so v_junction = -V(pn_c):
    a NEGATIVE cathode voltage forward-biases the junction into injection.
    """
    src = "pin" if feed == "in" else "pad"
    # The unused optical input is left undriven: its wires settle to exactly 0,
    # which is the physical "no light" we want. (Verified on the full deck.)
    return (
        f"* ring probe\n.include {MRM_CELL}\n"
        ".optical_port pin\n.optical_port pth\n"
        ".optical_port pad\n.optical_port pdr\n"
        f"Xl {src} fc_cw_laser power_mW=1.0 wavelength_nm={wl_nm:.6f}\n"
        f"Xr pin pth pad pdr 0 vc ht 0 mrm radius={radius:.9e} n_eff={N_EFF}\n"
        f"Vc vc 0 DC {v_c:.9g}\nIht 0 ht DC {i_ht:.9g}\n.op\n"
    )


def probe(radius: float, wl_nm: float, v_c=0.0, i_ht=0.0, feed="in"):
    """Return (P_thru mW, P_drop mW, I_junction A)."""
    c = fc.Circuit()
    c.load_str(probe_deck(radius, wl_nm, v_c, i_ht, feed))
    r = c.run("op")

    def p(port):
        return (float(r[f"V({port}_re_0)"][0]) ** 2
                + float(r[f"V({port}_im_0)"][0]) ** 2) / 1e-3

    return p("pth"), p("pdr"), float(r["I(vc)"][0])


def resonance_nm(radius: float, near_nm: float, span=0.9, n=181) -> float:
    """Drop-port resonance of `radius`, the interior peak nearest `near_nm`.

    Scans a window around a *predicted* wavelength rather than around the
    target, and refuses an edge maximum — an edge means the peak is outside the
    window, and taking it silently is how a solve ends up chasing its own tail.
    """
    for _ in range(4):
        wl = np.linspace(near_nm - span, near_nm + span, n)
        d = np.array([probe(radius, float(w))[1] for w in wl])
        peaks = [i for i in range(1, n - 1) if d[i] > d[i - 1] and d[i] > d[i + 1]]
        if peaks:
            i = min(peaks, key=lambda j: abs(wl[j] - near_nm))
            y0, y1, y2 = d[i - 1], d[i], d[i + 1]
            den = y0 - 2 * y1 + y2
            sh = 0.0 if den == 0 else 0.5 * (y0 - y2) / den
            return float(wl[i] + sh * (wl[1] - wl[0]))
        span *= 2.0          # nothing interior — widen and retry
    raise RuntimeError(f"no resonance found near {near_nm} nm for r={radius:e}")


def dlambda_dr(radius: float, near_nm: float) -> float:
    """Local sensitivity dλ_res/dr, measured.

    Not λ/r: the segment's first-order dispersion (n_eff walking toward n_g)
    cancels about half of the geometric shift, so the naive λ/r estimate
    overshoots by ~2x. Measure it instead of trusting the algebra.
    """
    dr = 20e-9
    a = resonance_nm(radius - dr, near_nm)
    b = resonance_nm(radius + dr, near_nm)
    return (b - a) / (2 * dr)


# ── stage 1: radii ───────────────────────────────────────────────────────────
def stage_radii(verbose=True) -> list[float]:
    """Solve radius_i so ring i is resonant with channel i."""
    if verbose:
        print("── stage 1: radius calibration ──")
    lam0 = resonance_nm(R0, LAMBDAS_NM[0])
    slope = dlambda_dr(R0, lam0)
    if verbose:
        print(f"  r = {R0 * 1e6:.4f} µm → {lam0:.4f} nm")
        naive = lam0 / R0 * 1e-9          # nm of λ per nm of radius
        print(f"  measured dλ_res/dr = {slope * 1e-9 * 1e3:.1f} pm per nm of "
              f"radius (naive λ/r says {naive * 1e3:.1f} pm/nm — dispersion "
              f"cancels {100 * (1 - slope * 1e-9 / naive):.0f}% of it)")
        print(f"  → {0.8 / slope * 1e9:.3f} nm of radius per 0.80 nm channel step")

    radii = []
    for k, target in enumerate(LAMBDAS_NM):
        r = R0 + (target - lam0) / slope       # first-order, measured slope
        got = resonance_nm(r, target)
        for _ in range(4):                     # secant on the real measurement
            if abs(got - target) < 5e-4:
                break
            r += (target - got) / slope
            got = resonance_nm(r, target)
        radii.append(r)
        if verbose:
            print(f"  ring {k + 1}: r = {r * 1e6:.6f} µm → {got:.4f} nm "
                  f"(err {(got - target) * 1e3:+.2f} pm)")
    if verbose:
        print("\n  paste into giona_rnn_perfectW.sp:")
        for k, r in enumerate(radii, start=1):
            print(f"  .param radius{k}={r:.7e}")
        step = np.diff(radii).mean()
        print(f"\n  step {step * 1e9:.3f} nm/channel — the deck's 13 nm step is "
              f"{13 / (step * 1e9):.2f}x too coarse (it spread the rings "
              f"{13 / (step * 1e9) * 0.8:.2f} nm apart, not 0.80 nm)")
    return radii


# ── stage 2: activation ──────────────────────────────────────────────────────
# Port map of the add-drop ring (verified against mrm.sp's coupler wiring):
#   in  --(straight)--> thru        NOTCH: dips at resonance
#   in  --(ring)------> drop        PEAK
#   add --(straight)--> drop        NOTCH
#   add --(ring)------> thru        PEAK   <- the one we want for a ReLU
# So feeding the source into the add/drop bus and reading the BUS output gives a
# peak response: no light off resonance, light when the ring is pulled onto the
# channel. Carrier injection blue-shifts (dn = -dn_di*I), so at rest the ring
# must sit RED of its channel and injection pulls it in — the heater sets that
# rest offset, since it red-shifts.
PORTS = {
    "peak": ("add", 0, "add→thru (peak) — ReLU-like"),
    "notch": ("in", 0, "in→thru (notch) — sigmoid-like"),
}


def activation(radius: float, wl_nm: float, v_list, shape="peak", i_ht=0.0):
    """Transmission (mW) and junction current (A) versus cathode voltage."""
    feed, port, _ = PORTS[shape]
    t, i = [], []
    for v in v_list:
        pt, pd, ij = probe(radius, wl_nm, v_c=float(v), i_ht=i_ht, feed=feed)
        t.append((pt, pd)[port])
        i.append(ij)
    return np.array(t), np.array(i)


def heater_shift_pm(radius: float, lam_ref: float, i_ht: float) -> float:
    """Resonance shift from heater current, in pm (positive = red)."""
    c = fc.Circuit()
    # Measure by finding the resonance with the heater on.
    lo = resonance_nm_heated(radius, lam_ref, i_ht)
    return (lo - lam_ref) * 1e3


def resonance_nm_heated(radius: float, near_nm: float, i_ht: float,
                        span=1.2, n=241) -> float:
    wl = np.linspace(near_nm - span, near_nm + span, n)
    d = np.array([probe(radius, float(w), i_ht=i_ht)[1] for w in wl])
    i = int(np.argmax(d))
    if 0 < i < n - 1:
        y0, y1, y2 = d[i - 1], d[i], d[i + 1]
        den = y0 - 2 * y1 + y2
        sh = 0.0 if den == 0 else 0.5 * (y0 - y2) / den
        return float(wl[i] + sh * (wl[1] - wl[0]))
    return float(wl[i])


def stage_activation(radii: list[float]):
    print("\n── stage 2: activation function ──")
    r1, lam = radii[0], LAMBDAS_NM[0]

    for v in (-0.9, +0.9):
        _, _, i = probe(r1, lam, v_c=v)
        print(f"  V(pn_c) = {v:+.1f} V → I_junction = {i * 1e6:+8.2f} µA "
              f"({'INJECTION (forward)' if abs(i) > 1e-6 else 'depletion (reverse)'})")
    print("  ⇒ the neuron node must go NEGATIVE to activate a ring.")

    print("\n  heater tuning (this is what sets the rest offset / threshold):")
    for i_ht in (0.0, 2e-3, 4e-3, 6e-3, 8e-3):
        sh = resonance_nm_heated(r1, lam, i_ht) - lam
        print(f"    {i_ht * 1e3:.1f} mA → {sh * 1e3:+7.1f} pm  "
              f"(P = {i_ht ** 2 * 2 * 184.4 * 1e3:.2f} mW)")

    v = np.linspace(0.0, -1.1, 111)
    for shape in ("peak", "notch"):
        print(f"\n  {PORTS[shape][2]}")
        print("    rest offset      T(0 V)    T(-0.7)   T(-0.9)   T(-1.1)   swing")
        for det_pm in (0, 100, 200, 300, 400):
            # Channel sits `det_pm` BLUE of the ring's resonance at rest, which
            # is what a heater does; probing below resonance is equivalent.
            t, _ = activation(r1, lam - det_pm * 1e-3, v, shape=shape)
            g = lambda vv: t[int(np.argmin(abs(v - vv)))]
            print(f"    {det_pm:4d} pm blue   {t[0]:7.4f}   {g(-0.7):7.4f}   "
                  f"{g(-0.9):7.4f}   {t[-1]:7.4f}   {t.max() - t.min():6.4f}")
    return v


# ── driving the hand-written deck ────────────────────────────────────────────
# `.param` values are substituted at parse time, so they cannot be reached with
# set_param. Overriding the text is the honest way to sweep them — and it keeps
# this script read-only with respect to Hugh's netlist.
def load_deck(**overrides):
    src = NETLIST.read_text()
    for k, v in overrides.items():
        pat = re.compile(rf"^(\.param\s+{re.escape(k)}\s*=\s*)\S+", re.M)
        src, n = pat.subn(lambda m: f"{m.group(1)}{v:.9g}", src)
        if n == 0:
            raise KeyError(f".param {k} not found in {NETLIST.name}")
    c = fc.Circuit()
    c.load_str(src)
    return c


CAL_RADII = {f"radius{k + 1}": r for k, r in enumerate(
    [7.9998737e-06, 8.0074122e-06, 8.0149580e-06, 8.0224153e-06,
     8.0300681e-06, 8.0376334e-06, 8.0452040e-06, 8.0527814e-06])}


def read_state(r):
    """Per-neuron junction voltage/current and per-channel powers of interest."""
    P = lambda p, k: (float(r[f"V({p}_re_{k})"][0]) ** 2
                      + float(r[f"V({p}_im_{k})"][0]) ** 2) / 1e-3
    v = np.array([float(r[f"V(mod_cathode{i})"][0]) for i in range(1, N + 1)])
    return {
        "v": v,
        "bus": np.array([P("bout8", k) for k in range(N)]),
        "win": np.array([[P(f"win{i}", k) for k in range(N)]
                         for i in range(1, N + 1)]),
    }


# ── stage 3: DC bias rule ────────────────────────────────────────────────────
def stage_bias(i_star_uA=90.0, powers_mW=None, weights=None):
    """Set every neuron's rest operating point, in the recurrent case.

    The rule. Neuron i's junction current is
        I_i = I_bias(PD_B_i) + g * sum_j w_ij * T_j(I_j)
    with g collecting responsivity and the tap/tree attenuation. In the
    feedforward case you know the incoming powers because they are upstream. In
    the recurrent case you do not — but you do not need to: CHOOSE the rest
    state I* for every neuron at once, and the required bias follows in closed
    form,
        I_bias(PD_B_i) = I* - g * sum_j w_ij * T_j(I*),
    because T_j(I*) is then a known constant. That makes I* a fixed point by
    construction; no iteration on the recurrence is involved. Stability of that
    fixed point is a separate question (for an oscillator you want it unstable).

    In practice I_bias(PD_B) is only piecewise-linear (the 2k/10k network plus a
    clamping diode), so this does the closed-form estimate and then two secant
    corrections against the measured node voltage.
    """
    print("\n── stage 3: DC operating point ──")
    powers = powers_mW if powers_mW is not None else [1.0] * N
    W = np.zeros((N, N)) if weights is None else np.asarray(weights, float)

    ov = dict(CAL_RADII)
    for k in range(N):
        ov[f"p_{k + 1}"] = powers[k]
    for i in range(N):
        for j in range(N):
            ov[f"w_{i + 1}{j + 1}"] = W[i, j]

    # Target node voltage for the chosen current, read off the diode clamp.
    r1 = CAL_RADII["radius1"]
    v_target = None
    for v in np.linspace(-0.5, -1.3, 81):
        _, _, ij = probe(r1, LAMBDAS_NM[0], v_c=float(v))
        if ij * 1e6 >= i_star_uA:
            v_target = float(v)
            break
    print(f"  target: I_junction = {i_star_uA:.0f} µA ⇒ V(mod_cathode) ≈ "
          f"{v_target:+.4f} V")

    bias = np.zeros(N)
    for it in range(8):
        c = load_deck(**ov, **{f"PD_B{i + 1}": bias[i] for i in range(N)})
        st = read_state(c.run("op"))
        err = st["v"] - v_target
        print(f"  iter {it}: PD_B = {bias[0]:+7.3f} V, V = {st['v'][0]:+.4f} V, "
              f"max |V-V*| = {np.abs(err).max() * 1e3:6.2f} mV")
        if np.abs(err).max() < 2e-3:
            break
        # Re-measure dV/dPD_B every step: the diode's incremental resistance
        # collapses as it turns on, so a slope taken at V=0 badly overshoots.
        c2 = load_deck(**ov, **{f"PD_B{i + 1}": bias[i] - 0.25 for i in range(N)})
        st2 = read_state(c2.run("op"))
        slope = (st["v"] - st2["v"]) / 0.25
        slope = np.where(np.abs(slope) > 1e-4, slope, 1 / 6)
        bias = bias - err / slope
    print("\n  paste into giona_rnn_perfectW.sp:")
    for i in range(N):
        print(f"  .param PD_B{i + 1}={bias[i]:.6f}")
    if np.abs(bias).max() > 3.0:
        print(f"\n  !! {np.abs(bias).max():.2f} V exceeds the +/-3 V rails in the deck.")
        print("     The 2k shunt to ground (R1) sinks ~440 uA at this operating")
        print("     point, so 10k (Rb3) needs several volts to push 90 uA of")
        print("     injection through the junction. Other reachable knobs:")
        for name, lo, hi in (("PD_P1", 3.0, 1.0), ("PD_N1", -3.0, -1.0)):
            a = load_deck(**ov, **{name: lo}); va = read_state(a.run("op"))["v"][0]
            b = load_deck(**ov, **{name: hi}); vb = read_state(b.run("op"))["v"][0]
            print(f"       d V(mod_cathode) / d {name} = "
                  f"{(vb - va) / (hi - lo) * 1e3:+7.2f} mV/V   "
                  f"(via the 17k PD shunt)")
    return bias, v_target, ov


# ── plots ────────────────────────────────────────────────────────────────────
def stage_plots(radii):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    RESULTS.mkdir(exist_ok=True)

    fig, ax = plt.subplots(2, 3, figsize=(16, 9))

    # (0,0) resonance comb after radius calibration
    wl = np.linspace(1545.6, 1552.3, 400)
    for k, (r, lam) in enumerate(zip(radii, LAMBDAS_NM)):
        d = np.array([probe(r, float(w))[1] for w in wl])
        ax[0, 0].plot(wl, d, lw=1.0, label=f"ring {k + 1}")
        ax[0, 0].axvline(lam, color="k", lw=0.4, ls=":")
    ax[0, 0].set_xlabel("wavelength (nm)"); ax[0, 0].set_ylabel("drop (mW)")
    ax[0, 0].set_title("Calibrated comb: ring i on channel i")
    ax[0, 0].legend(fontsize=6, ncol=2); ax[0, 0].grid(alpha=0.3)

    # (0,1) heater tuning
    ih = np.linspace(0, 8e-3, 17)
    sh = [ (resonance_nm(radii[0], LAMBDAS_NM[0] + 1.6 * (i / 8e-3) ** 2 * 1.0
                         if i > 0 else LAMBDAS_NM[0], span=1.0)
            if False else 0) for i in ih ]
    sh = []
    for i in ih:
        pred = LAMBDAS_NM[0] + 11.4 * (np.pi * i ** 2 * 2 * 184.4 / 26.4e-3) / (2 * np.pi)
        sh.append(resonance_nm_heated_at(radii[0], pred, float(i)) - LAMBDAS_NM[0])
    sh = np.array(sh)
    ax[0, 1].plot(ih * 1e3, sh * 1e3, "o-")
    ax[0, 1].axhline(800, color="r", ls="--", lw=0.8, label="one channel (800 pm)")
    ax[0, 1].set_xlabel("heater current (mA)"); ax[0, 1].set_ylabel("Δλ_res (pm)")
    ax[0, 1].set_title("Heater: coarse lock / channel placement")
    ax[0, 1].legend(fontsize=8); ax[0, 1].grid(alpha=0.3)

    # (0,2) node voltage vs injected current — the diode clamp
    I = np.linspace(0, 800e-6, 41)
    vn = [probe_i(radii[0], LAMBDAS_NM[0], float(i))[2] for i in I]
    ax[0, 2].plot(I * 1e6, vn)
    ax[0, 2].set_xlabel("junction current (µA)")
    ax[0, 2].set_ylabel("V(mod_cathode) (V)")
    ax[0, 2].set_title("The diode clamps the node: current is the loop variable")
    ax[0, 2].grid(alpha=0.3)

    # (1,0)/(1,1) activation vs current, both ports
    I = np.concatenate([np.linspace(0, 200e-6, 21), np.linspace(220e-6, 900e-6, 18)])
    for col, (shape, title) in enumerate(
            (("notch", "NOTCH (in→thru): sigmoid / saturating ReLU"),
             ("peak", "PEAK (add→thru): bump, turns over"))):
        a = ax[1, col]
        for det in (0, 100, 200, 300):
            t = [probe_i(radii[0], LAMBDAS_NM[0] - det * 1e-3, float(i),
                         feed=PORTS[shape][0])[PORTS[shape][1]] for i in I]
            a.plot(I * 1e6, t, "-", lw=1.3, label=f"{det} pm blue")
        a.set_xlabel("junction current (µA)"); a.set_ylabel("transmission (mW)")
        a.set_title(title, fontsize=10); a.legend(fontsize=8); a.grid(alpha=0.3)

    # (1,2) available photocurrent vs laser power — the loop-gain budget
    pw = np.array([1, 2, 5, 10, 20])
    avail = []
    for P0 in pw:
        c = load_deck(**CAL_RADII, **{f"p_{k + 1}": float(P0) for k in range(N)})
        st = read_state(c.run("op"))
        avail.append(st["win"][0].sum() * 1e-3 * RESPONSIVITY * 1e6)
    ax[1, 2].plot(pw, avail, "o-", label="Σ_j |photocurrent| at a block")
    ax[1, 2].axhline(90, color="r", ls="--", lw=0.8, label="sigmoid half-max (90 µA)")
    ax[1, 2].axhline(400, color="darkred", ls=":", lw=0.8, label="saturation (400 µA)")
    ax[1, 2].set_xlabel("laser power per channel (mW)")
    ax[1, 2].set_ylabel("available photocurrent (µA)")
    ax[1, 2].set_title("Loop-gain budget: tap+tree costs 16×")
    ax[1, 2].legend(fontsize=8); ax[1, 2].grid(alpha=0.3); ax[1, 2].set_xscale("log")

    fig.suptitle("giona RNN calibration — activation, tuning knobs, link budget")
    fig.tight_layout()
    out = RESULTS / "rnn_activation.png"
    fig.savefig(out, dpi=110)
    print(f"wrote {out}")
    return out


def resonance_nm_heated_at(radius, near_nm, i_ht, span=1.0, n=161):
    wl = np.linspace(near_nm - span, near_nm + span, n)
    d = np.array([probe(radius, float(w), i_ht=i_ht)[1] for w in wl])
    i = int(np.argmax(d))
    if 0 < i < n - 1:
        y0, y1, y2 = d[i - 1], d[i], d[i + 1]
        den = y0 - 2 * y1 + y2
        sh = 0.0 if den == 0 else 0.5 * (y0 - y2) / den
        return float(wl[i] + sh * (wl[1] - wl[0]))
    return float(wl[i])


def probe_i(radius, wl_nm, i_inj, i_ht=0.0, feed="in"):
    """Ring driven by a junction CURRENT source. Returns (thru, drop, V_node)."""
    src = "pin" if feed == "in" else "pad"
    deck = (f"* ring, current-driven\n.include {MRM_CELL}\n"
            ".optical_port pin\n.optical_port pth\n"
            ".optical_port pad\n.optical_port pdr\n"
            f"Xl {src} fc_cw_laser power_mW=1.0 wavelength_nm={wl_nm:.6f}\n")
    if feed == "add":
        # The coupler routes BOTH outputs' lambda from port a1, so an add-fed
        # ring has no wavelength label unless the bus port carries one. Drive
        # just that wire: it is a label, not a signal.
        deck += f"Vwl pin_wl_0 0 DC {wl_nm * 1e-9:.9e}\n"
    deck += (f"Xr pin pth pad pdr 0 vc ht 0 mrm radius={radius:.9e}"
             f" n_eff={N_EFF}\n"
             f"Ij vc 0 DC {i_inj:.9e}\nRl vc 0 1e9\nIht 0 ht DC {i_ht:.9e}\n.op\n")
    c = fc.Circuit(); c.load_str(deck); r = c.run("op")
    P = lambda p: (float(r[f"V({p}_re_0)"][0]) ** 2
                   + float(r[f"V({p}_im_0)"][0]) ** 2) / 1e-3
    return P("pth"), P("pdr"), float(r["V(vc)"][0])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--radii", action="store_true")
    ap.add_argument("--activation", action="store_true")
    ap.add_argument("--bias", action="store_true")
    ap.add_argument("--plots", action="store_true")
    args = ap.parse_args()

    if args.radii or args.activation or args.plots:
        radii = stage_radii(verbose=args.radii or args.activation)
        if args.activation:
            stage_activation(radii)
        if args.plots:
            stage_plots(radii)
    if args.bias:
        stage_bias()
    return 0


if __name__ == "__main__":
    sys.exit(main())
