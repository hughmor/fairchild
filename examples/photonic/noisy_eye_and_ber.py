#!/usr/bin/env python3
"""A 10 Gb/s electro-optic link, built from primitives, with the noise in it.

CW laser → Mach-Zehnder modulator → photodiode → load. The modulator is not a
behavioural block: it is two directional couplers and two reverse-biased PN
phase shifters, the same devices you would place in a layout, driven push-pull
by a PRBS. Everything below is measured from the circuit rather than asserted
about it.

**Can primitives reach 10 GHz with realistic numbers?** Yes, and the tradeoff is
the classic lumped-modulator one — `.ac` measures it:

    arms    C_j0     driver    f_3dB (modulator)   V_pi/arm
    1 mm    250 fF    25 ohm         52.0 GHz         12.0 V
    2 mm    500 fF    25 ohm         26.0 GHz          6.0 V
    2 mm    500 fF    50 ohm         13.2 GHz          6.0 V
    3 mm    750 fF    25 ohm         17.3 GHz          4.0 V
    5 mm   1250 fF    25 ohm         10.4 GHz          2.4 V

Bandwidth and drive voltage trade against each other through the arm length —
both `f_3dB` and `V_pi` go as `1/L`, because the junction capacitance grows with
length while the phase efficiency does too. Reproduce the table with `--sweep`.
This example takes 3 mm: 17.3 GHz for the modulator, 8.3 GHz for the link once
the receiver's own pole is in series, against a 5 GHz Nyquist.

**Is the noise real?** It is the point of the example, so it is checked three
ways, not asserted: the sampled sigma from `.tran` against the `.noise` integral
at that rail, both against the closed form `sigma/mu = sqrt(RIN·B)` in the
RIN-limited regime, and linearity in `noisescale`.

**Why is this link RIN-limited, when published eyes usually are not?** Because
of how much light lands on the detector. The budget printed below comes out

    RIN·I²      93.7 %       I_ph = 1.57 mA
    shot 2qI     6.1 %
    thermal      0.2 %

which is the correct answer for these numbers and *not* how a real receiver is
built. RIN scales as `I²` while shot goes as `I` and thermal not at all, so RIN
wins whenever the photocurrent is large — and a milliamp into a kilohm is a
1.5 V swing, which no real front end wants. A load resistor big enough to give a
readable voltage is what puts a link in the RIN-limited corner, and it is slow
for the same reason: `1/(2πRC)` and `4kT/R` pull opposite ways.

**`--tia` is the answer to that**, and it is worth running:

    receiver              sigma_1    sigma_0    dominant term
    load resistor        11839 uV     547 uV    RIN, 94 %
    Verilog-A TIA         4114 uV    4001 uV    amplifier input noise, 93 %

A transimpedance amplifier holds the summing node at a virtual ground, so the
diode's capacitance no longer sets the bandwidth and a real receiver can detect
50 µW instead of 2 mW. Its own input-referred current noise — 15 pA/√Hz here —
then sets the floor. Note the second column: with the TIA the two rails have
*the same* noise, because the dominant generator is in the amplifier and does
not care how much light arrived. That is what leaving the RIN-limited regime
looks like, and it is why published eye diagrams do not show noise riding the
one rail the way the resistor version does.

`--tia` needs `examples/verilog_a/build/va_tia.osdi`, which means OpenVAF; the
default front end is a resistor so the example runs with no toolchain at all.

Note also that the modulator is perfectly balanced, so the extinction ratio is
numerically infinite (the 70+ dB below is the solver's floor, not a spec). A real
MZM reaches 20-30 dB, and the 0 rail's noise and the eye closure both depend on
that.

    python3 examples/photonic/noisy_eye_and_ber.py [--selftest] [--pam4]
"""

import argparse
import math
import os
import sys

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))
import fairchild  # noqa: E402

# ── link design ──────────────────────────────────────────────────────────────
BIT_S = 100e-12                 # 10 Gb/s (and 10 GBd for PAM-4)
STEP_S = 1e-12
EDGE_FRAC = 0.25                # driver rise/fall as a fraction of a bit

P_LASER_MW = 2.0             # resistor front end
P_LASER_MW_TIA = 0.05        # a real receiver detects far less
RIN_DB_HZ = -145.0              # a decent DFB; -155 is a good one, -130 a cheap one

L_ARM_UM = 3000.0               # 3 mm arms
V_PI_L = 0.012                  # V*m (1.2 V*cm) -> V_pi = 4 V per arm at 3 mm
C_J0 = 750e-15                  # 250 fF/mm
V_BI, M_J = 0.917, 0.5
ALPHA_DB_CM = 2.0
R_DRV = 25.0                    # driver output impedance per arm
V_CENTRE = -3.0                 # both arms reverse-biased

RESPONSIVITY = 0.9
R_LOAD = 1.0e3
C_PD = 15e-15                   # tau = 15 ps

# ── receiver front end ───────────────────────────────────────────────────────
# Two of them, and the difference is the point. The resistor is what you can
# run with no toolchain; the TIA is what a real link uses. See `--tia`.
TIA_OSDI = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "verilog_a", "build", "va_tia.osdi")
)
TIA_ZT = 2.0e3                  # V/A
TIA_RIN = 50.0                  # ohm — the virtual ground the diode sees
TIA_F3DB = 12e9                 # Hz
TIA_IN = 15e-12                 # A/sqrt(Hz), input-referred
TIA_ROUT = 50.0
TIA_RLOAD = 1.0e6               # high-Z probe, so z_t is not divided

V_PI = V_PI_L / (L_ARM_UM * 1e-6)      # per arm
SWING = V_PI / 4.0                     # per arm; differential swings V_pi

KL_3DB = math.pi / 4            # kappa*L for a 50/50 coupler

PRBS_ORDER = 9                  # 511 symbols, the short-pattern standard
N_SYM = (1 << PRBS_ORDER) - 1

# ── pattern generation ───────────────────────────────────────────────────────
# A real maximal-length PRBS, not a hand-typed pattern: the run-length
# distribution (up to `order` consecutive ones) is what loads a receiver's
# baseline and closes an eye, and it is the reason a short repeating pattern
# flatters a link.
PRBS_TAPS = {7: (7, 6), 9: (9, 5), 11: (11, 9), 15: (15, 14)}


def prbs(order: int, n: int) -> list:
    """Maximal-length LFSR sequence, Fibonacci form, all-ones seed."""
    taps = PRBS_TAPS[order]
    reg = (1 << order) - 1
    out = []
    for _ in range(n):
        fb = 0
        for t in taps:
            fb ^= (reg >> (t - 1)) & 1
        out.append(reg & 1)
        reg = ((reg << 1) | fb) & ((1 << order) - 1)
    return out


def pam4_levels(bits: list) -> list:
    """Gray-coded PAM-4 symbols in 0..3 from a bit stream (MSB first)."""
    gray = {(0, 0): 0, (0, 1): 1, (1, 1): 2, (1, 0): 3}
    return [gray[(bits[2 * i], bits[2 * i + 1])] for i in range(len(bits) // 2)]


def drive_fractions(levels: list, n_levels: int) -> list:
    """Map symbol -> differential drive as a fraction of V_pi.

    The MZM's power transfer is `sin²(Δφ/2)` with `Δφ = π·d`, so equally spaced
    *optical* levels need unequally spaced *drive* levels. Pre-distorting here
    is what a real PAM-4 transmitter's DAC does, and skipping it is why an
    un-corrected MZM PAM-4 eye has three unequal openings.
    """
    return [(2.0 / math.pi) * math.asin(math.sqrt(s / (n_levels - 1.0))) for s in levels]


def pwl(fracs: list, sign: int) -> str:
    """Push-pull drive for one arm.

    arm(d) = V_CENTRE + sign·2·SWING·d, so the *differential* arm voltage is
    4·SWING·d = V_pi·d and the phase difference is π·d.
    """
    half = 0.5 * EDGE_FRAC * BIT_S

    def level(d):
        return V_CENTRE + sign * 2.0 * SWING * d

    # Points only at transitions, and each transition straddles the symbol
    # boundary: hold the old level to `k·T − half`, reach the new one by
    # `k·T + half`. Ramping *after* the boundary instead puts every crossing a
    # quarter-bit late, which shows up as an eye whose opening is not where the
    # sampler is.
    pts = [f"0 {level(fracs[0]):.6g}"]
    prev = fracs[0]
    for k in range(1, len(fracs)):
        if fracs[k] != prev:
            pts.append(f"{k * BIT_S - half:.6e} {level(prev):.6g}")
            pts.append(f"{k * BIT_S + half:.6e} {level(fracs[k]):.6g}")
            prev = fracs[k]
    pts.append(f"{len(fracs) * BIT_S:.6e} {level(prev):.6g}")
    return f"PWL({' '.join(pts)})"


# ── the circuit ──────────────────────────────────────────────────────────────
PORTS = "".join(
    f".optical_port {p}\n" for p in ("lin", "ldark", "a1", "a2", "b1", "b2", "obar", "ocross")
)


def receiver(use_tia: bool) -> tuple:
    """`(preamble, netlist_fragment, probe_node)` for the chosen front end."""
    if not use_tia:
        return "", (f"Rl det 0 {R_LOAD}\nCl det 0 {C_PD}\n"), "det"
    return (
        f".osdi {TIA_OSDI}\n",
        f"Cpd det 0 {C_PD}\n"
        f"Xtia det tout 0 va_tia z_t={TIA_ZT} r_in={TIA_RIN} f_3db={TIA_F3DB} "
        f"i_n_in={TIA_IN} v_out_dc=0 v_swing=1.5 r_out={TIA_ROUT}\n"
        f"Rl tout 0 {TIA_RLOAD}\n",
        "tout",
    )


def link(drive_p: str, drive_n: str, trannoise=False, seed=1, scale=1.0,
         use_tia=False) -> str:
    """CW laser → 50/50 coupler → two PN arms → 50/50 coupler → photodiode."""
    noise = (
        f".options trannoise=1 noiseseed={seed} noisescale={scale} variable_step=0\n"
        if trannoise
        else ""
    )
    arm = (
        f"fc_pn_ps_cap l_um={L_ARM_UM} v_pi_l={V_PI_L} c_j0={C_J0} "
        f"v_bi={V_BI} m_j={M_J} alpha_dB_cm={ALPHA_DB_CM} pin_at_ref=1"
    )
    pre, rx, _ = receiver(use_tia)
    p_mw = P_LASER_MW_TIA if use_tia else P_LASER_MW
    return f"""* 10 Gb/s MZM link
{noise}{pre}{PORTS}Xlas lin fc_cw_laser power_mW={p_mw} rin_db_hz={RIN_DB_HZ}
Vp pd 0 {drive_p}
Vn nd 0 {drive_n}
Rp pd p {R_DRV}
Rn nd n {R_DRV}
Xc1 lin ldark a1 a2 fc_dcoupler kappa_L={KL_3DB}
Xps1 a1 b1 p 0 {arm}
Xps2 a2 b2 n 0 {arm}
Xc2 b1 b2 obar ocross fc_dcoupler kappa_L={KL_3DB}
Xpd obar det 0 fc_photodetector responsivity={RESPONSIVITY} r_shunt=1Meg i_dark_a=0
{rx}.end
"""


def run_tran(fracs, trannoise=False, seed=1, scale=1.0, use_tia=False):
    stop = len(fracs) * BIT_S
    c = fairchild.Circuit()
    c.load_str(link(pwl(fracs, +1), pwl(fracs, -1), trannoise, seed, scale, use_tia))
    r = c.run("tran", step=STEP_S, stop=stop)
    return np.asarray(r.time()), np.asarray(r[f"V({receiver(use_tia)[2]})"])


def bandwidth(l_um=L_ARM_UM, c_j0=C_J0, r_drv=R_DRV, c_pd=C_PD):
    """Small-signal electro-optic response at quadrature, and its 3 dB point.

    Measured optically and end to end — the photodiode's voltage is what a
    receiver actually sees. Pass a tiny `c_pd` to take the receiver's own pole
    out of the way and isolate the modulator.
    """
    q = 0.5  # quadrature: half of V_pi differential
    dp = f"DC {V_CENTRE + 2 * SWING * q} AC 0.5"
    dn = f"DC {V_CENTRE - 2 * SWING * q} AC 0.5 180"
    arm = (
        f"fc_pn_ps_cap l_um={l_um} v_pi_l={V_PI_L} c_j0={c_j0} "
        f"v_bi={V_BI} m_j={M_J} alpha_dB_cm={ALPHA_DB_CM} pin_at_ref=1"
    )
    deck = f"""* MZM small-signal response
{PORTS}Xlas lin fc_cw_laser power_mW={P_LASER_MW}
Vp pd 0 {dp}
Vn nd 0 {dn}
Rp pd p {r_drv}
Rn nd n {r_drv}
Xc1 lin ldark a1 a2 fc_dcoupler kappa_L={KL_3DB}
Xps1 a1 b1 p 0 {arm}
Xps2 a2 b2 n 0 {arm}
Xc2 b1 b2 obar ocross fc_dcoupler kappa_L={KL_3DB}
Xpd obar det 0 fc_photodetector responsivity={RESPONSIVITY} r_shunt=1Meg i_dark_a=0
Rl det 0 {R_LOAD}
Cl det 0 {c_pd}
"""
    c = fairchild.Circuit()
    c.load_str(deck)
    r = c.run("ac", fstart=1e7, fstop=1e12, points=40, variation="dec")
    f = np.asarray(r.freq())
    db = 20.0 * np.log10(np.abs(np.asarray(r["V(det)"])) / abs(r["V(det)"][0]))
    k = int(np.argmax(db < -3.0))
    f3 = np.interp(-3.0, [db[k], db[k - 1]], [f[k], f[k - 1]]) if k else float("nan")
    return f, db, f3


def noise_sigma(frac: float, use_tia: bool = False) -> float:
    """RMS output noise with the modulator parked at drive fraction `frac`.

    `.noise` linearises about ONE bias. A modulated link does not have one — RIN
    and shot noise both follow the optical power — so this is called once per
    rail rather than once for the deck.
    """
    dp = f"DC {V_CENTRE + 2 * SWING * frac}"
    dn = f"DC {V_CENTRE - 2 * SWING * frac}"
    c = fairchild.Circuit()
    c.load_str(link(dp, dn, use_tia=use_tia))
    r = c.run("noise", out=receiver(use_tia)[2], src="Vp",
              fstart=1e3, fstop=1e13, points=40, variation="dec")
    return math.sqrt(np.trapezoid(np.asarray(r["onoise"]), np.asarray(r.freq())))


# ── sampling ─────────────────────────────────────────────────────────────────
SKIP = 8  # let the first few symbols settle before measuring anything


def sample_at(t, v, n_sym, phase=0.5, settled=None):
    """One sample per symbol at `phase` through the symbol.

    `settled` keeps only symbols whose two predecessors carried the same value,
    and it matters more than it looks. The rails differ in noise by ~100x here,
    so a zero straight after a one is still discharging the one's noise through
    the receiver's RC: measured on the 0 rail, sampling every zero reports 3.5 mV
    where the rail's own noise is 0.5 mV — a 570 % error that looks like a broken
    noise model and is really intersymbol interference. An isolated zero is not
    at the zero rail's noise level, and no scope would pretend otherwise.
    """
    idx = np.arange(SKIP, n_sym)
    if settled is not None:
        idx = np.array([k for k in idx if settled[k - 2] == settled[k - 1] == settled[k]])
    return idx, np.interp((idx + phase) * BIT_S, t, v)


def eye_traces(t, v, n_sym, n_show=220):
    """2-UI windows centred on a symbol boundary, which is what puts a crossing
    at each edge of the plot and the eye in the middle."""
    out = []
    for k in range(SKIP, min(n_sym - 1, SKIP + n_show)):
        t0 = (k - 0.5) * BIT_S
        m = (t >= t0) & (t < t0 + 2 * BIT_S)
        if m.sum() > 4:
            out.append(((t[m] - t0) / BIT_S - 0.5, v[m] * 1e3))
    return out


# ── main ─────────────────────────────────────────────────────────────────────
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--tia", action="store_true",
                    help="read the diode with the Verilog-A TIA instead of a load "
                         "resistor (needs examples/verilog_a/build/va_tia.osdi)")
    ap.add_argument("--sweep", action="store_true",
                    help="measure the length/bandwidth/V_pi tradeoff and exit")
    ap.add_argument("--png", default="noisy_eye_and_ber.png")
    args = ap.parse_args()

    if args.sweep:
        print(f"{'arms':>8}{'C_j0':>9}{'driver':>9}{'f_3dB (mod)':>14}{'V_pi/arm':>11}")
        for l_um, c_j0, r in ((1000.0, 250e-15, 25.0), (2000.0, 500e-15, 25.0),
                              (2000.0, 500e-15, 50.0), (3000.0, 750e-15, 25.0),
                              (5000.0, 1250e-15, 25.0)):
            _, _, f3m = bandwidth(l_um=l_um, c_j0=c_j0, r_drv=r, c_pd=1e-18)
            print(f"{l_um / 1000:>6.0f} mm{c_j0 * 1e15:>7.0f} fF{r:>7.0f} Ω"
                  f"{f3m / 1e9:>11.1f} GHz{V_PI_L / (l_um * 1e-6):>9.1f} V")
        return 0

    use_tia = args.tia
    if use_tia and not os.path.exists(TIA_OSDI):
        sys.exit(f"--tia needs {TIA_OSDI}\n"
                 "  build it:  cd examples/verilog_a && OPENVAF=... ./build.sh")
    rx_name = "Verilog-A TIA" if use_tia else "load resistor"
    probe = receiver(use_tia)[2]

    bits = prbs(PRBS_ORDER, N_SYM)
    nrz = [float(b) for b in bits]
    sym4 = pam4_levels(prbs(PRBS_ORDER, 2 * (N_SYM // 2)))
    pam = drive_fractions(sym4, 4)

    # The modulator on its own (receiver pole pushed far out), then the link.
    _, _, f3_mod = bandwidth(c_pd=1e-18)
    f_ac, db_ac, f3 = bandwidth()
    print(f"receiver:  {rx_name}, probing V({probe})")
    print(f"modulator: V_pi = {V_PI:.1f} V/arm, f_3dB = {f3_mod / 1e9:.1f} GHz "
          f"({L_ARM_UM / 1000:.0f} mm arms, {C_J0 * 1e15:.0f} fF, {R_DRV:.0f} ohm driver)")
    print(f"link (modulator + receiver): f_3dB = {f3 / 1e9:.1f} GHz, against a "
          f"{1 / (2 * BIT_S) / 1e9:.0f} GHz Nyquist")

    # Clean run: the deterministic waveform, for levels and for subtracting.
    t0, v0 = run_tran(nrz, use_tia=use_tia)
    idx, s0 = sample_at(t0, v0, N_SYM, settled=bits)
    ones0 = s0[[bits[k] == 1 for k in idx]]
    zeros0 = s0[[bits[k] == 0 for k in idx]]
    mu1, mu0 = ones0.mean(), zeros0.mean()
    print(f"NRZ levels: mu1 = {mu1 * 1e3:.0f} mV, mu0 = {mu0 * 1e3:.1f} mV, "
          f"ER = {10 * math.log10(mu1 / max(mu0, 1e-12)):.1f} dB")

    # `.noise` per rail.
    sig1_ac, sig0_ac = noise_sigma(1.0, use_tia), noise_sigma(0.0, use_tia)

    # Transient noise, pooled over seeds, at several amplitudes.
    scales = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
    seeds = (11, 12)
    sig1, sig0 = [], []
    for sc in scales:
        hi, lo = [], []
        for sd in seeds:
            _, v = run_tran(nrz, trannoise=True, seed=sd, scale=float(sc), use_tia=use_tia)
            _, s = sample_at(t0, v, N_SYM, settled=bits)
            hi.append(s[[bits[k] == 1 for k in idx]] - ones0)
            lo.append(s[[bits[k] == 0 for k in idx]] - zeros0)
        sig1.append(np.std(np.concatenate(hi)))
        sig0.append(np.std(np.concatenate(lo)))
    sig1, sig0 = np.array(sig1), np.array(sig0)

    # The closed form, so the number defends itself rather than matching a
    # second measurement of the same thing. In the RIN-limited regime
    # sigma/mu = sqrt(RIN·B_n) with B_n the receiver's noise bandwidth.
    rin = 10 ** (RIN_DB_HZ / 10.0)
    b_n = ((math.pi / 2) * TIA_F3DB if use_tia
           else (math.pi / 2) / (2 * math.pi * R_LOAD * C_PD))
    sig1_cf = 0.0   # filled once the budget below is assembled

    # Which term dominates, and by how much. This is the part that decides
    # whether an eye looks RIN-limited, and it is set by the photocurrent.
    q_e, k_b, t_k = 1.602176634e-19, 1.380649e-23, 300.0
    if use_tia:
        z_eff = TIA_ZT * TIA_RLOAD / (TIA_ROUT + TIA_RLOAD)
        i_ph, b_n = mu1 / z_eff, (math.pi / 2) * TIA_F3DB
        terms = {"TIA input-referred": TIA_IN ** 2}
    else:
        z_eff = R_LOAD
        i_ph = mu1 / R_LOAD
        terms = {"thermal 4kT/R_L": 4 * k_b * t_k / R_LOAD}
    terms["shot 2qI"] = 2 * q_e * i_ph
    terms["RIN·I²"] = rin * i_ph * i_ph
    tot = sum(terms.values())
    sig1_cf = math.sqrt(tot * b_n) * z_eff
    print(f"\nnoise budget at the 1 rail — {rx_name} "
          f"(I_ph = {i_ph * 1e6:.1f} uA, B_n = {b_n / 1e9:.1f} GHz):")
    for name, psd in sorted(terms.items(), key=lambda kv: -kv[1]):
        print(f"  {name:<20}{psd:9.2e} A²/Hz  {100 * psd / tot:5.1f} %   "
              f"σ = {math.sqrt(psd * b_n) * z_eff * 1e6:7.0f} uV")

    rel1 = abs(sig1[0] - sig1_ac) / sig1_ac
    rel0 = abs(sig0[0] - sig0_ac) / sig0_ac
    rel_cf = abs(sig1[0] - sig1_cf) / sig1_cf
    slope0 = np.polyfit(scales, sig0, 1)[0] * scales[0] / sig0[0]

    print(f"\n{'':12}{'.noise':>12}{'.tran':>12}{'closed form':>14}")
    print(f"{'sigma_1':12}{sig1_ac * 1e6:9.0f} uV{sig1[0] * 1e6:9.0f} uV{sig1_cf * 1e6:11.0f} uV")
    print(f"{'sigma_0':12}{sig0_ac * 1e6:9.0f} uV{sig0[0] * 1e6:9.0f} uV{'--':>14}")
    print(f"\n  .tran vs .noise: 1 rail {100 * rel1:+.1f} %, 0 rail {100 * rel0:+.1f} %")
    print(f"  .tran vs closed form on the 1 rail: {100 * rel_cf:+.1f} %")
    print(f"  RIN-limited SNR ceiling 1/(RIN·B) = {10 * math.log10(1 / (rin * b_n)):.1f} dB "
          f"in B_n = {b_n / 1e9:.1f} GHz")
    print(f"  sigma_1/mu_1 = {100 * sig1[0] / mu1:.2f} %,  sigma_1/sigma_0 = "
          f"{sig1[0] / sig0[0]:.1f}x")
    if use_tia:
        print("  sigma_1 ~ sigma_0 because the dominant generator is the amplifier's own\n"
              "  input noise, which does not care how much light arrived. That is what\n"
              "  leaving the RIN-limited regime looks like.")

    q = (mu1 - mu0) / (sig1 + sig0)
    ber = 0.5 * np.array([math.erfc(x / math.sqrt(2.0)) for x in q])
    print(f"\n{'noisescale':>11}{'sigma_1':>11}{'sigma_0':>10}{'Q':>8}{'BER':>11}")
    for sc, h, l, qq, bb in zip(scales, sig1, sig0, q, ber):
        print(f"{sc:>11.0f}{h * 1e6:>9.0f} uV{l * 1e6:>8.0f} uV{qq:>8.1f}{bb:>11.1e}")

    # PAM-4.
    tp, vp = run_tran(pam, use_tia=use_tia)
    _, sp = sample_at(tp, vp, len(pam))
    lv = [sp[[sym4[k] == L for k in range(SKIP, len(pam))]].mean() for L in range(4)]
    print(f"\nPAM-4 levels (mV): " + ", ".join(f"{x * 1e3:.0f}" for x in lv))
    gaps = np.diff(lv)
    print(f"  eye openings: " + ", ".join(f"{g * 1e3:.0f}" for g in gaps) +
          f"  (spread {100 * (gaps.max() - gaps.min()) / gaps.mean():.0f} % — "
          "pre-distortion is what keeps them equal)")

    if args.selftest:
        assert f3_mod > 10e9, f"the modulator itself must clear 10 GHz: {f3_mod / 1e9:.1f} GHz"
        # A link carries NRZ at ~0.7x the bit rate; 7 GHz for 10 Gb/s.
        assert f3 > 0.7 / BIT_S, f"the link must carry 10 Gb/s: f_3dB = {f3 / 1e9:.1f} GHz"
        assert rel1 < 0.12, f"1 rail: .noise vs .tran differ by {100 * rel1:.1f} %"
        assert rel0 < 0.12, f"0 rail: .noise vs .tran differ by {100 * rel0:.1f} %"
        assert rel_cf < 0.10, f"1 rail vs closed form: {100 * rel_cf:.1f} %"
        # The quiet rail is the linearity check: it is thermal-dominated, so it
        # scales exactly. The loud one is RIN-dominated and bends slightly by
        # 16x, which is the square-law detector and not a defect.
        assert abs(slope0 - 1.0) < 0.02, f"quiet rail must be linear in noisescale: {slope0:.3f}"
        slope1 = np.polyfit(scales, sig1, 1)[0] * scales[0] / sig1[0]
        assert 0.98 < slope1 < 1.10, f"loud rail should be near-linear: {slope1:.3f}"
        assert 10 * math.log10(mu1 / max(mu0, 1e-12)) > 15.0, "extinction ratio too low"
        assert q[0] > 6.0, f"the 1x eye should be open: Q = {q[0]:.1f}"
        if use_tia:
            # The whole point of the TIA: the noise stops tracking the signal.
            assert sig1[0] / sig0[0] < 1.5, \
                f"a TIA-limited receiver's rails should have similar noise: " \
                f"{sig1[0] / sig0[0]:.1f}x"
            assert terms["TIA input-referred"] / tot > 0.8, \
                "the TIA should dominate its own receiver"
        else:
            assert sig1[0] > 10.0 * sig0[0], \
                "a RIN-limited receiver should be visibly level-dependent"
        assert (gaps.max() - gaps.min()) / gaps.mean() < 0.20, \
            "pre-distorted PAM-4 openings should be within 20 % of each other"
        print("\nselftest OK — modulator carries 10 Gb/s, both noise analyses agree "
              "with each other and with the closed form, PAM-4 levels are equalised")
        return 0

    import matplotlib.pyplot as plt

    # Pick the eye's noisescale from the measured noise rather than hard-coding
    # it: enough fuzz to see in print, not enough to close the eye. The two
    # front ends are an order of magnitude apart in sigma/separation — 0.76 %
    # for the resistor, 5.3 % for the TIA — so one constant cannot serve both,
    # and the one tuned for the resistor buries the TIA's eye completely.
    target = 0.06
    eye_scale = float(np.clip(round(target / (sig1[0] / (mu1 - mu0))), 1.0, 16.0))
    q_at_eye = float(np.interp(eye_scale, scales, q))
    _, v_nrz = run_tran(nrz, trannoise=True, seed=seeds[0], scale=eye_scale, use_tia=use_tia)
    _, v_pam = run_tran(pam, trannoise=True, seed=seeds[0], scale=eye_scale, use_tia=use_tia)

    fig, ax = plt.subplots(2, 2, figsize=(11.5, 7.4), constrained_layout=True)

    for a, v, n, title in (
        (ax[0][0], v_nrz, N_SYM, f"NRZ, 10 Gb/s — Q = {q_at_eye:.0f}"),
        (ax[0][1], v_pam, len(pam), "PAM-4, 10 GBd (20 Gb/s)"),
    ):
        for x, y in eye_traces(t0[: len(v)], v, n):
            a.plot(x, y, color="C0", alpha=0.10, lw=0.8)
        a.axvline(0.0, color="0.55", ls=":", lw=0.9)
        a.set_xlabel("symbol phase (UI)")
        a.set_ylabel("V(det)  (mV)")
        a.set_xlim(-0.5, 1.5)
        a.set_title(f"{title}   (noisescale={eye_scale:.0f})")
        a.grid(alpha=0.25)

    a = ax[1][0]
    a.semilogx(f_ac, db_ac, color="C2", lw=1.8)
    a.axhline(-3, color="0.6", ls=":", lw=1)
    a.axvline(f3, color="C3", ls="--", lw=1.2)
    a.annotate(f"{f3 / 1e9:.1f} GHz", (f3, -14), fontsize=9, color="C3",
               ha="right", rotation=90, va="bottom")
    a.axvline(1 / (2 * BIT_S), color="0.4", ls="-.", lw=1)
    a.annotate("Nyquist, 10 Gb/s", (1 / (2 * BIT_S), -28), fontsize=8, color="0.35",
               ha="right", rotation=90, va="bottom")
    a.set_xlim(2e8, 1e11)
    a.set_ylim(-30, 3)
    a.set_xlabel("frequency (Hz)")
    a.set_ylabel("electro-optic response (dB)")
    a.set_title(f"Link response — modulator alone: {f3_mod / 1e9:.0f} GHz\n"
                f"{L_ARM_UM / 1000:.0f} mm arms, {R_DRV:.0f} Ω driver, V$_\\pi$ = {V_PI:.0f} V/arm")
    a.grid(True, which="both", alpha=0.25)

    a = ax[1][1]
    a.plot(scales, sig1_ac * scales * 1e6, "-", color="C3", lw=5, alpha=0.30,
           label="σ₁  ∫.noise df at the 1 rail")
    a.plot(scales, sig1 * 1e6, "o--", color="C3", lw=1.3, ms=5, label="σ₁  sampled from .tran")
    a.plot(scales, sig1_cf * scales * 1e6, ":", color="0.25", lw=1.6,
           label="σ₁  μ·√(RIN·B), closed form")
    a.plot(scales, sig0_ac * scales * 1e6, "-", color="C0", lw=5, alpha=0.30,
           label="σ₀  ∫.noise df at the 0 rail")
    a.plot(scales, sig0 * 1e6, "o--", color="C0", lw=1.3, ms=5, label="σ₀  sampled from .tran")
    a.set_xlabel("noisescale")
    a.set_ylabel("output noise σ (µV)")
    a.set_yscale("log")
    a.set_title("Each rail against its own operating point")
    a.legend(fontsize=7.5, loc="center right")
    a.grid(True, which="both", alpha=0.25)

    out = os.path.join(os.path.dirname(__file__), args.png)
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
