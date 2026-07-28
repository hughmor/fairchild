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
  --space       Explore the activation space the way the chip is wired: the
                real 2k/10k bias network, heater current for detuning, and a
                signed input photocurrent.
  --insitu      Neuron transconductance vs bias, measured in the full 8-neuron
                deck rather than on a single-ring bench.
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
    "drop": ("add", 0, "add→thru (drop) — ReLU-like"),
    "thru": ("in", 0, "in→thru (thru) — sigmoid-like"),
}


def activation(radius: float, wl_nm: float, v_list, shape="drop", i_ht=0.0):
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
    for shape in ("drop", "thru"):
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


# ── stage 2b: activation space, with the real bias network ───────────────────
# `probe_i` above drives the junction from a bare current source and detunes by
# moving the laser. Neither is what the chip does. This probe reproduces the DC
# environment of one neuron exactly as `neuron_junction_wdm8.sp` builds it:
#
#      mod_cathode ──┬── 2k ──────────────── gnd          (on-chip shunt)
#                    ├── 10k ─────────────── PD_B         (off-chip bias R)
#                    ├── 17k+1+50 ────────── PD_N = -3 V  (thru-PD shunt path)
#                    ├── 17k+1+50 ────────── PD_P = +3 V  (drop-PD shunt path)
#                    ├── ring PN junction ── gnd          (anode grounded)
#                    └── I_sig                            (stands in for the
#                                                          balanced photocurrent)
#
# The two ±3 V paths matter: they are a balanced pair that cancels only at
# V = 0, so at a forward-biased node they act as a further ~8.5k to ground.
# Ignoring them (as the old probe did) understates the shunt by ~20 %.
#
# Validated against the real deck with the lasers off: this probe reproduces
# V(mod_cathode1) to 5 decimals at PD_B = -13.4 V and -17.8 V.
#
# Sign convention for I_sig follows the hardware. The thru PD sits with its
# CATHODE on mod_cathode, so its photocurrent pulls the node down (more forward
# bias, more injection); the drop PD sits anode-on-node and pushes it up. So
# I_sig > 0 means "net thru-side light" = activating, and the sweep runs through
# zero to negative to cover a negative weighted sum.
R_ON_CHIP = 2e3
R_BIAS = 10e3
R_PD_SHUNT = 17e3 + 1.0 + 50.0
PD_P, PD_N = 3.0, -3.0
# Small-signal shunt seen by I_sig: everything to a stiff supply, in parallel.
R_SHUNT_AC = 1.0 / (1.0 / R_ON_CHIP + 1.0 / R_BIAS + 2.0 / R_PD_SHUNT)
N_DIODE, V_T = 5.0, 1.380649e-23 * 300.15 / 1.602176634e-19


def bias_deck(i_sig: float, i_ht: float, pd_b: float, radius: float,
              wl_nm: float, feed: str, power_mW: float) -> str:
    src = "pin" if feed == "in" else "pad"
    return (
        f"* ring activation with the real DC bias network\n.include {MRM_CELL}\n"
        ".optical_port pin\n.optical_port pth\n"
        ".optical_port pad\n.optical_port pdr\n"
        f"Xl {src} fc_cw_laser power_mW={power_mW:.6g} wavelength_nm={wl_nm:.6f}\n"
        f"Xr pin pth pad pdr 0 pn_c ht 0 mrm radius={radius:.9e} n_eff={N_EFF}\n"
        # 0 V source purely to read the junction current.
        "Vsense pn_c mc DC 0\n"
        f"R1 mc 0 {R_ON_CHIP:g}\n"
        f"Rb3 mc bias_com {R_BIAS:g}\n"
        f"Vb bias_com 0 DC {pd_b:.9g}\n"
        f"Rsh1 mc th_a 17e3\nRse1 th_a th_a1 1\nRb1 th_a1 bias_neg 50\n"
        f"Vpdn bias_neg 0 DC {PD_N:g}\n"
        f"Rsh2 mc dr_c 17e3\nRse2 dr_c dr_c1 1\nRb2 dr_c1 bias_pos 50\n"
        f"Vpdp bias_pos 0 DC {PD_P:g}\n"
        f"Isig mc 0 DC {i_sig:.9e}\n"
        f"Iht 0 ht DC {i_ht:.9e}\n.op\n"
    )


def bias_probe(i_sig=0.0, i_ht=0.0, pd_b=-13.4, radius=None, wl_nm=None,
               feed="add", power_mW=1.0) -> dict:
    """One ring in its real DC bias network. Returns transmissions, node V, I_j."""
    radius = CAL_R[0] if radius is None else radius
    wl_nm = LAMBDAS_NM[0] if wl_nm is None else wl_nm
    c = fc.Circuit()
    c.load_str(bias_deck(i_sig, i_ht, pd_b, radius, wl_nm, feed, power_mW))
    r = c.run("op")
    p = lambda q: (float(r[f"V({q}_re_0)"][0]) ** 2
                   + float(r[f"V({q}_im_0)"][0]) ** 2) / 1e-3
    return {"thru": p("pth"), "drop": p("pdr"),
            "v": float(r["V(mc)"][0]), "i_j": float(r["I(vsense)"][0])}


CAL_R = [7.9998737e-06, 8.0074122e-06, 8.0149580e-06, 8.0224153e-06,
         8.0300681e-06, 8.0376334e-06, 8.0452040e-06, 8.0527814e-06]

# Peak port (add→thru): dark off resonance, bright when the ring is pulled onto
# the channel. Injection blue-shifts, the heater red-shifts, so the heater sets
# WHERE ALONG THE CURRENT AXIS the ring crosses the channel — it is the
# threshold knob, the hardware's version of a bias term.
I_SIG = np.linspace(-800e-6, 900e-6, 70)


# Which (feed, read) pair gives which lineshape. Both useful shapes read the
# thru port; what distinguishes them is which bus the light entered on, so the
# feed cannot be inferred from the port name.
SHAPES = {
    "peak": ("add", "thru", "peak port (add→thru) — ReLU-like"),
    "notch": ("in", "thru", "notch port (in→thru) — sigmoid-like"),
}


def activation_curve(i_ht: float, pd_b: float, shape="peak", i_sig=None):
    """Transmission, node voltage and junction current over the I_sig sweep."""
    i_sig = I_SIG if i_sig is None else i_sig
    feed, port, _ = SHAPES[shape]
    out = [bias_probe(float(i), i_ht=i_ht, pd_b=pd_b, feed=feed) for i in i_sig]
    return (np.array([o[port] for o in out]),
            np.array([o["v"] for o in out]),
            np.array([o["i_j"] for o in out]))


def eta_measured(pd_b: float, i_ht=0.0, h=20e-6) -> tuple[float, float]:
    """(rest junction current, dI_j/dI_sig) at the operating point.

    The fraction of signal photocurrent that reaches the junction instead of
    leaking into the shunts — the AC coupling efficiency, measured rather than
    inferred from a small-signal model.
    """
    a = bias_probe(-h, i_ht=i_ht, pd_b=pd_b)
    b = bias_probe(+h, i_ht=i_ht, pd_b=pd_b)
    rest = bias_probe(0.0, i_ht=i_ht, pd_b=pd_b)["i_j"]
    return rest, (b["i_j"] - a["i_j"]) / (2 * h)


def centre_heater(pd_b: float, lo=0.0, hi=6e-3, coarse=25, fine=13) -> float:
    """Heater current that puts the ring's resonance on its channel at I_sig = 0.

    Coarse scan then a local refine. Not a root-find on a formula: the injection
    that comes with the bias already blue-shifts the ring, so how much heat is
    needed depends on the bias — the two knobs are not independent.
    """
    feed, port, _ = SHAPES["peak"]
    grid = np.linspace(lo, hi, coarse)
    t = [bias_probe(0.0, i_ht=float(i), pd_b=pd_b, feed=feed)[port] for i in grid]
    k = int(np.argmax(t))
    step = grid[1] - grid[0]
    g2 = np.linspace(max(lo, grid[k] - step), min(hi, grid[k] + step), fine)
    t2 = [bias_probe(0.0, i_ht=float(i), pd_b=pd_b, feed=feed)[port] for i in g2]
    return float(g2[int(np.argmax(t2))])


def slope_at(pd_b: float, i_ht: float, h=15e-6) -> float:
    """dT/dI_sig (mW per µA) at I_sig = 0, on the peak port.

    Per unit *input photocurrent*, not per unit junction current, so it already
    folds in the shunt loss — this is the number that closes the optical loop.
    """
    feed, port, _ = SHAPES["peak"]
    a = bias_probe(-h, i_ht=i_ht, pd_b=pd_b, feed=feed)[port]
    b = bias_probe(+h, i_ht=i_ht, pd_b=pd_b, feed=feed)[port]
    return (b - a) / (2 * h * 1e6)


def flank_heater(pd_b: float, lo=0.0, hi=6e-3, coarse=31, fine=15) -> float:
    """Heater current that puts the STEEPEST FLANK at I_sig = 0.

    Not `centre_heater`: centring the resonance on the channel parks the
    operating point at the peak, where dT/dI is zero by construction. A neuron
    wants maximum sensitivity, which is halfway down the flank.
    """
    grid = np.linspace(lo, hi, coarse)
    sl = [abs(slope_at(pd_b, float(i))) for i in grid]
    k = int(np.argmax(sl))
    step = grid[1] - grid[0]
    g2 = np.linspace(max(lo, grid[k] - step), min(hi, grid[k] + step), fine)
    s2 = [abs(slope_at(pd_b, float(i))) for i in g2]
    best = float(g2[int(np.argmax(s2))])
    # An optimum at the range floor is not an optimum: the heater can only add
    # red shift, so it means the ring is already past its steepest flank at this
    # bias and no heater setting recovers it. Flag it rather than reporting the
    # boundary value as if it were a choice. (Same trap `resonance_nm` refuses.)
    return best, best <= lo + 1e-12 or best >= hi - 1e-12


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
    # ax[0, 0].set_title("Calibrated comb: ring i on channel i")
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
    # ax[0, 1].set_title("Heater: coarse lock / channel placement")
    ax[0, 1].legend(fontsize=8); ax[0, 1].grid(alpha=0.3)

    # (0,2) node voltage vs injected current — the diode clamp
    I = np.linspace(0, 800e-6, 41)
    vn = [probe_i(radii[0], LAMBDAS_NM[0], float(i))[2] for i in I]
    ax[0, 2].plot(I * 1e6, vn)
    ax[0, 2].set_xlabel("mod junction current (µA)")
    ax[0, 2].set_ylabel("mof cathode voltage (V)")
    # ax[0, 2].set_title("The diode clamps the node: current is the loop variable")
    ax[0, 2].grid(alpha=0.3)

    # (1,0)/(1,1) activation vs current, both ports
    I = np.concatenate([np.linspace(0, 200e-6, 21), np.linspace(220e-6, 900e-6, 18)])
    for col, (shape, title) in enumerate(
            (("thru", "NOTCH (in→thru): sigmoid / saturating ReLU"),
             ("drop", "PEAK (add→thru): bump, turns over"))):
        a = ax[1, col]
        for det in (0, 100, 200, 300):
            t = [probe_i(radii[0], LAMBDAS_NM[0] - det * 1e-3, float(i),
                         feed=PORTS[shape][0])[PORTS[shape][1]] for i in I]
            a.plot(I * 1e6, t, "-", lw=1.3, label=f"{det} pm blue")
        a.set_xlabel("mod junction current (µA)"); a.set_ylabel(f"{shape} transmission (mW)")
        # a.set_title(title, fontsize=10)
        a.legend(fontsize=8)
        a.grid(alpha=0.3)

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
    ax[1, 2].set_ylabel("max photocurrent (µA)")
    # ax[1, 2].set_title("Loop-gain budget: tap+tree costs 16×")
    ax[1, 2].legend(fontsize=8); ax[1, 2].grid(alpha=0.3); ax[1, 2].set_xscale("log")

    fig.suptitle("giona simulated calibration (`fairchild`)")
    fig.tight_layout()
    out = RESULTS / "rnn_activation.png"
    fig.savefig(out, dpi=110)
    print(f"wrote {out}")
    return out


def insitu_transfer(pd_b: float, power_mW=30.0, ws=(0.0, 0.10, 0.20, 0.30)):
    """Neuron 1's transfer measured in the full 8-neuron deck, not on a bench.

    Sweeps neuron 1's own weight and watches its channel on the bus output. The
    input photocurrent is w · P_at_block · R, both of which the deck reports, so
    the slope is directly comparable to the single-ring |dT/dI_sig|.

    Deliberately a real weight range rather than a tiny probe: the closed-loop
    sensitivity difference at w = 0.05 lives in the 6th decimal, and dividing
    that by w amplifies solver tolerance into a meaningless "loop gain".
    """
    i_in, p_bus, v0 = [], [], None
    for w in ws:
        ov = dict(CAL_RADII)
        ov |= {f"p_{k + 1}": power_mW for k in range(N)}
        ov |= {f"PD_B{k + 1}": pd_b for k in range(N)}
        ov["w_11"] = w
        st = read_state(load_deck(**ov).run("op"))
        if v0 is None:
            v0 = st["v"][0]
        i_in.append(w * st["win"][0][0] * 1e-3 * RESPONSIVITY)
        p_bus.append(st["bus"][0])
    i_in, p_bus = np.array(i_in) * 1e6, np.array(p_bus)
    slope = float(np.polyfit(i_in, p_bus, 1)[0])
    return {"v_rest": v0, "i_in": i_in, "p_bus": p_bus, "slope": abs(slope)}


def stage_insitu(power_mW=30.0,
                 biases=(-5.0, -8.0, -10.0, -13.4, -17.8, -22.0)):
    """Does the bench measurement hold in the real circuit?

    Same transconductance question asked of the full 8-neuron deck at its real
    laser power. This is what confirmed the rest-point choice.
    """
    print(f"\n── in situ, full deck at {power_mW:.0f} mW/channel "
          "(neuron 1's own weight swept) ──")
    print("    PD_B    V_rest    |dP_bus/dI_in|   I_in at w=0.3")
    ins = {}
    for pb in biases:
        d = ins[pb] = insitu_transfer(pb, power_mW=power_mW)
        print(f"   {pb:6.1f}  {d['v_rest']:+8.5f}   {d['slope']:11.5f}   "
              f"{d['i_in'][-1]:7.1f} µA")
    ref = ins.get(-13.4, {"slope": float("nan")})["slope"]
    best = max(biases, key=lambda pb: ins[pb]["slope"])
    print(f"\n  → best at PD_B = {best:.1f} V: "
          f"{ins[best]['slope'] / ref:.1f}× the transconductance of the "
          f"-13.4 V the WTA runs used.")
    return ins


def stage_space():
    """Explore the activation space the way the chip is actually wired.

    Three knobs, and the point of the figure is that they are not orthogonal:
      I_sig  — the balanced photocurrent, positive OR negative (the input)
      I_ht   — heater current: sets WHERE on the current axis the ring crosses
               its channel, i.e. the threshold
      PD_B   — bias supply through the real 10k: sets the rest injection
               current, which fixes both the electrical coupling efficiency and
               (because injection blue-shifts) how much heat the threshold costs
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    RESULTS.mkdir(exist_ok=True)

    pd_ref = -13.4                      # the bias the WTA runs at
    # Four biases across the top row. Six on the bottom row, where the panels
    # are one point per bias and can afford the extra resolution.
    top_biases = [-5.0, -8.0, -13.4, -22.0]
    biases = [-5.0, -8.0, -10.0, -13.4, -17.8, -22.0]
    # Heater window relative to each bias's own centring current. Absolute
    # values would not be comparable across panels: the injection that comes
    # with the bias already blue-shifts the ring, so -5 V needs ~0.5 mA to hold
    # its channel and -22 V needs ~3.5 mA (that is finding 3 below).
    # Spanned rather than offset, so clamping at 0 mA still leaves four distinct
    # curves instead of two coincident ones.
    HT_BELOW, HT_ABOVE, HT_N = 1.0e-3, 0.5e-3, 4

    print("\n── activation space, real bias network ──")
    rest = bias_probe(0.0, pd_b=pd_ref)
    print(f"  PD_B = {pd_ref} V, I_sig = 0  →  V(node) = {rest['v']:+.5f} V, "
          f"I_junction = {rest['i_j'] * 1e6:.1f} µA")
    print(f"  small-signal shunt (2k ‖ 10k ‖ 2×17.05k) = {R_SHUNT_AC:.0f} Ω")

    fig, ax = plt.subplots(2, 4, figsize=(21, 9))
    # Top row shares y: the peak height genuinely falls as the bias rises, and
    # per-panel autoscale would hide exactly that.
    for c in range(1, 4):
        ax[0, c].sharey(ax[0, 0])

    # Top row: one panel per bias. Both lineshapes together — solid for the peak
    # port (add→thru, ReLU-like), dotted for the notch port (in→thru,
    # sigmoid-like) — since they are the same ring read two ways and the useful
    # comparison is between them at a given operating point.
    print("\n  top row: heater currents used, relative to each bias's centring "
          "current")
    for col, pb in enumerate(top_biases):
        ih_c = centre_heater(pb)
        hts = list(np.linspace(max(0.0, ih_c - HT_BELOW), ih_c + HT_ABOVE, HT_N))
        print(f"    PD_B = {pb:6.1f} V: centring {ih_c * 1e3:.2f} mA, using "
              + ", ".join(f"{h * 1e3:.2f}" for h in hts) + " mA")
        a = ax[0, col]
        for i, ih in enumerate(hts):
            c = f"C{i}"
            for sh, ls in (("peak", "-"), ("notch", ":")):
                t = activation_curve(ih, pb, shape=sh)[0]
                a.plot(I_SIG * 1e6, t, ls, color=c, lw=1.5,
                       label=f"{ih * 1e3:.2f} mA" if sh == "peak" else None)
        a.axvline(0, color="k", lw=0.5, ls=":")
        a.set_xlabel("signal photocurrent I_sig (µA)   → activating")
        if col == 0:
            a.set_ylabel("transmission (mW)")
        a.set_title(f"PD_B = {pb:.1f} V   (rest "
                    f"{bias_probe(0.0, pd_b=pb)['i_j'] * 1e6:.0f} µA)\n"
                    "solid: add→thru (peak/ReLU)   dotted: in→thru "
                    "(notch/sigmoid)", fontsize=9)
        a.legend(fontsize=7, title="heater", title_fontsize=7)
        a.grid(alpha=0.3)

    # Removed from the figure but kept for reference: what the heater knob does
    # on its own — threshold position and peak height versus heater current, at
    # one bias. Superseded by the top row, which shows the same translation
    # happening inside each panel.
    #
    # To bring it back, widen the grid — ax[0, 2] now holds a bias panel.
    #
    # curves = {ih: {sh: activation_curve(ih, pd_ref, shape=sh)
    #                for sh in SHAPES} for ih in [0.0, 1.5e-3, 2.0e-3, 2.5e-3,
    #                                             3.0e-3, 3.5e-3]}
    # thr = [I_SIG[int(np.argmax(curves[ih]["peak"][0]))] * 1e6 for ih in curves]
    # pk = [curves[ih]["peak"][0].max() for ih in curves]
    # a = ax[0, 2]
    # a.plot(np.array(list(curves)) * 1e3, thr, "o-", color="C0")
    # a.set_xlabel("heater current (mA)")
    # a.set_ylabel("threshold: I_sig at peak (µA)", color="C0")
    # a.tick_params(axis="y", labelcolor="C0")
    # a.axhline(0, color="k", lw=0.5, ls=":")
    # a2 = a.twinx()
    # a2.plot(np.array(list(curves)) * 1e3, pk, "s--", color="C3")
    # a2.set_ylabel("peak transmission (mW)", color="C3")
    # a2.tick_params(axis="y", labelcolor="C3")
    # a.set_title("The heater is the threshold knob\n"
    #             "…but a threshold at higher injection costs contrast",
    #             fontsize=10)
    # a.grid(alpha=0.3)

    # (0,3)/(1,0) the node itself: the diode clamps the voltage, so the loop
    # variable is current — and only part of I_sig gets past the shunts.
    node = {pb: activation_curve(0.0, pb, shape="peak") for pb in biases}
    a = ax[1, 0]
    for pb in biases:
        a.plot(I_SIG * 1e6, node[pb][1] * 1e3, lw=1.3, label=f"{pb:.1f} V")
    a.set_xlabel("signal photocurrent I_sig (µA)")
    a.set_ylabel("node voltage V(mod_cathode) (mV)")
    a.set_title("The junction clamps the node\n"
                "so current, not voltage, is the loop variable", fontsize=10)
    a.legend(fontsize=7, title="PD_B", title_fontsize=7); a.grid(alpha=0.3)

    a = ax[1, 1]
    for pb in biases:
        a.plot(I_SIG * 1e6, node[pb][2] * 1e6, lw=1.3, label=f"{pb:.1f} V")
    a.plot(I_SIG * 1e6, I_SIG * 1e6, "k:", lw=0.8, label="slope 1 (no loss)")
    a.set_xlabel("signal photocurrent I_sig (µA)")
    a.set_ylabel("junction current I_j (µA)")
    a.set_title("Only part of I_sig reaches the junction\n"
                "the rest leaks into 2 kΩ ‖ 10 kΩ ‖ 2×17 kΩ", fontsize=10)
    a.legend(fontsize=7, title="PD_B", title_fontsize=7); a.grid(alpha=0.3)

    # coupling efficiency vs rest current: measured against the divider.
    a = ax[1, 2]
    meas = [eta_measured(pb) for pb in biases]
    ir = np.array([m[0] for m in meas]) * 1e6
    et = np.array([m[1] for m in meas])
    a.plot(ir, et * 100, "o", color="C0", label="measured  dI_j/dI_sig")
    grid = np.logspace(np.log10(max(ir.min(), 1.0)), np.log10(ir.max()), 100)
    r_d = N_DIODE * V_T / (grid * 1e-6)
    a.plot(grid, 100 * R_SHUNT_AC / (R_SHUNT_AC + r_d), "-", color="0.5",
           label=f"R_sh/(R_sh+r_d), R_sh={R_SHUNT_AC:.0f} Ω")
    for x, y, pb in zip(ir, et, biases):
        a.annotate(f"{pb:.0f} V", (x, y * 100), fontsize=6,
                   textcoords="offset points", xytext=(3, -8))
    a.set_xscale("log")
    a.set_xlabel("rest junction current (µA)")
    a.set_ylabel("fraction of I_sig reaching the junction (%)")
    a.set_title("Bias sets the electrical coupling\n"
                "r_d = n·V_T/I falls, so more current couples better", fontsize=10)
    a.legend(fontsize=7); a.grid(alpha=0.3, which="both")

    # (1,2) the actual rest-point decision. At each bias, re-trim the heater so
    # the steepest flank — not the peak — sits at I_sig = 0, then ask what that
    # operating point is worth. dT/dI_sig is the number that closes the loop: it
    # already contains the shunt loss, since it is per unit input photocurrent.
    a = ax[1, 3]
    rows = []
    for pb in biases:
        ih, at_edge = flank_heater(pb)
        t, _, _ = activation_curve(ih, pb, shape="peak")
        rest_i, eta = eta_measured(pb, i_ht=ih)
        rows.append((pb, ih, rest_i, eta, t.max(), abs(slope_at(pb, ih)), at_edge))
    print("\n  rest-point choice — heater re-trimmed per bias so the STEEPEST "
          "FLANK sits at I_sig = 0:")
    print("    PD_B    heater   I_rest    eta    peak T   |dT/dI_sig|")
    print("     (V)     (mA)     (µA)     (%)     (mW)      (mW/µA)")
    for pb, ih, ri, e, pt, sl, edge in rows:
        print(f"   {pb:6.1f}  {ih * 1e3:6.2f}  {ri * 1e6:7.1f}  {e * 100:5.1f}  "
              f"{pt:7.4f}  {sl:10.5f}"
              + ("   <- heater at range floor: already past the flank" if edge else ""))
    ri = np.array([r[2] for r in rows]) * 1e6
    sl = np.array([r[5] for r in rows])
    a.plot(ri, [r[3] * 100 for r in rows], "o-", color="C0",
           label="coupling η (%)")
    a.plot(ri, [r[4] / max(r[4] for r in rows) * 100 for r in rows], "s-",
           color="C3", label="resonance contrast (norm.)")
    a.plot(ri, sl / sl.max() * 100, "^-", color="C2", lw=2.2,
           label="|dT/dI_sig| — the loop gain (norm.)")
    best = rows[int(np.argmax(sl))]
    a.axvline(best[2] * 1e6, color="C2", ls="--", lw=0.8)
    a.annotate(f"best {best[5]:.4f} mW/µA\nat {best[2] * 1e6:.0f} µA\n"
               f"(PD_B {best[0]:.1f} V, heater {best[1] * 1e3:.2f} mA)",
               (0.52, 0.86), xycoords="axes fraction", fontsize=7, color="C2")
    a.set_xscale("log")
    a.set_xlabel("rest junction current (µA)")
    a.set_ylabel("normalised (%)")
    a.set_title("Choosing the rest current\n"
                "coupling wants it high, free-carrier loss wants it low",
                fontsize=10)
    a.legend(fontsize=7, loc="lower left"); a.grid(alpha=0.3, which="both")

    # Removed from the figure but kept for reference: the same question asked of
    # the full 8-neuron deck at its real laser power, which is what confirmed
    # the rest-point choice (-8 V gave 9.3x the transconductance of -13.4 V).
    # Still runnable on its own — `rnn_explore.py --insitu` prints the table.
    #
    # To bring it back, widen the grid — ax[1, 3] now holds the rest-point panel.
    #
    # a = ax[1, 3]
    # ins = {pb: insitu_transfer(pb) for pb in biases}
    # for pb in biases:
    #     d = ins[pb]
    #     a.plot(d["i_in"], d["p_bus"], "o-", lw=1.3,
    #            label=f"{pb:.1f} V  ({d['slope']:.4f})")
    # a.set_xlabel("input photocurrent to neuron 1 (µA)")
    # a.set_ylabel("ring 1's channel on the bus output (mW)")
    # a.set_title("In situ: full 8-neuron deck at 30 mW/channel\n"
    #             "legend shows |dP/dI_in| in mW/µA", fontsize=10)
    # a.legend(fontsize=6.5, title="PD_B", title_fontsize=6.5); a.grid(alpha=0.3)

    fig.suptitle("giona ring activation — swept through the real bias network "
                 "(2 kΩ on-chip ‖ 10 kΩ off-chip, heater detuning, signed input)")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = RESULTS / "rnn_activation_space.png"
    fig.savefig(out, dpi=110)
    print(f"\nwrote {out}")
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
    ap.add_argument("--space", action="store_true",
                    help="activation space through the real bias network")
    ap.add_argument("--insitu", action="store_true",
                    help="transconductance vs bias measured on the full deck")
    args = ap.parse_args()

    if args.space:
        stage_space()
    if args.insitu:
        stage_insitu()

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
