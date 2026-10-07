#!/usr/bin/env python3
"""A ring modulator driven faster than its cavity can respond.

One deck. One option. Two answers that disagree by more than the signal.

A micro-ring stores light. Detune it and the stored field does not vanish — it
beats against the new drive and decays over the cavity's photon lifetime. Every
frequency-domain tool, and most behavioural compact models, carry a *static*
resonance: they evaluate the ring's steady-state transfer at the instantaneous
voltage and move on. That is right while the bit period is long compared with
the photon lifetime, and it is wrong the moment it is not.

Here both answers come out of the same netlist and the same device models. The
only difference is one line:

    .options optical_delay=0      # the round trip is instantaneous
    .options waveguide_delay=1    # the round trip takes n_g·L/c

fairchild says so itself before either run finishes::

    warning: this deck has a closed optical path (through rd), and optical
    group delays are off — so its round trip is instantaneous and the cavity
    has no photon lifetime. Wavelength sweeps are unaffected; a transient
    response is not.

**Nothing below is fitted.** The ring's ringing period and its decay are each
predicted by a closed form that appears nowhere in any device:

    beat period   1 / (Δφ/2π · FSR),   FSR = c/(n_g·L)
    amplitude τ   −T_rt / ln(a·t),     T_rt = n_g·L/c

where ``a`` is the round-trip amplitude transmission from ``alpha_dB_cm`` and
``t = cos(kappa_L)`` is the coupler. The measurements sit on both.

**The trade the static model cannot see.** Weak coupling buys extinction and
costs photon lifetime, and the two are the same parameter. A designer choosing
``kappa_L`` from a static sweep is reading half the answer. Panel D is the
other half: where in (coupling × bit rate) the static model is still usable,
and where it stops being.

    python3 examples/photonic/ring_photon_lifetime.py [--selftest] [--png FILE]
"""

import argparse
import math
import os
import sys

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))
import fairchild  # noqa: E402

C = 2.99792458e8

# ── the ring ─────────────────────────────────────────────────────────────────
# A 10 µm-radius ring: 62.8 µm round trip, which is one FSR of 1.14 THz.
L_RING_UM = 62.8
N_G = 4.2
ALPHA_DB_CM = 2.0               # a good silicon waveguide
V_PI_L = 2e-3                   # V·m -> V_pi = 31.8 V over this ring
LAMBDA_NM = 1550.0
P_LASER_MW = 1.0
RESPONSIVITY = 0.9
R_LOAD = 1.0e3

L_RING_M = L_RING_UM * 1e-6
T_RT = N_G * L_RING_M / C                       # round-trip time, 0.88 ps
FSR_HZ = 1.0 / T_RT
A_RT = 10.0 ** (-ALPHA_DB_CM * (L_RING_M * 1e2) / 20.0)   # amplitude per trip
V_PI = V_PI_L / L_RING_M

# The headline case. kappa_L = 0.2 gives a 41 ps photon lifetime, so a 100 ps
# bit is marginal and a 25 ps bit is hopeless — which is the point.
KAPPA_L = 0.20
V_STEP = 0.5                    # a step of ~2 linewidths, not a swing to V_pi


def tau_amplitude(kappa_l: float) -> float:
    """Field amplitude decay time of the loaded cavity, from the closed form.

    `a·t` is the amplitude surviving one round trip, so the field falls by that
    factor every `T_rt` — an exponential whose time constant is this. Energy
    decays twice as fast; the beat against a steady drive is a field cross-term,
    so it is this one that shows in the waveform.
    """
    return -T_RT / math.log(A_RT * math.cos(kappa_l))


def loaded_q(kappa_l: float) -> float:
    """`Q = ω₀·τ_energy`, and the energy decay is half the amplitude decay."""
    omega_0 = 2.0 * math.pi * C / (LAMBDA_NM * 1e-9)
    return omega_0 * tau_amplitude(kappa_l) / 2.0


def detuning_hz(dv: float) -> float:
    """Resonance shift for an applied `dv`, as a frequency.

    A round-trip phase of `Δφ` moves the resonance by `Δφ/2π` of a free spectral
    range. Nothing in the deck computes this; it is what the ringing must beat at.
    """
    return (math.pi * dv / V_PI) / (2.0 * math.pi) * FSR_HZ


def deck(kappa_l: float, drive: str, delays: bool) -> str:
    """The ring, as a netlist. `delays` is the entire difference between models."""
    option = "waveguide_delay=1" if delays else "optical_delay=0"
    return f"""* micro-ring modulator
.options {option}
.optical_port lin
.optical_port rb
.optical_port thru
.optical_port rd
Xlas lin fc_cw_laser power_mW={P_LASER_MW} wavelength_nm={LAMBDA_NM}
Xdc  lin rb thru rd fc_dcoupler kappa_L={kappa_l}
Xring rd rb vmod 0 fc_pn_ps
+    L_um={L_RING_UM} V_pi_L={V_PI_L} g_pn=1e-3
+    alpha_dB_cm={ALPHA_DB_CM} n_g={N_G} wl_ref_nm={LAMBDA_NM} pin_at_ref=1
Xpd  thru det 0 fc_photodetector responsivity={RESPONSIVITY}
Rl   det 0 {R_LOAD}
{drive}
"""


def run(kappa_l, drive, delays, step, stop):
    c = fairchild.Circuit()
    c.load_str(deck(kappa_l, drive, delays))
    r = c.run("tran", step=step, stop=stop)
    return np.asarray(r.time()), np.asarray(r["V(det)"])


# Full transmission puts this much on the load, so every trace divides by it and
# becomes a transmission rather than a volt.
V_FULL = RESPONSIVITY * P_LASER_MW * 1e-3 * R_LOAD


def static_transfer(kappa_l, n=401, v_lo=-2.0, v_hi=2.0):
    """Transmission against bias, with the cavity instantaneous.

    A slow ramp with delays off *is* the quasi-static model: the optical loop
    closes algebraically inside one Newton solve, so every point is the ring's
    steady state at that instant's voltage. This is the curve an S-matrix tool
    gives, and the lookup the behavioural models apply.
    """
    span = 20e-9
    t, v = run(kappa_l, f"Vmod vmod 0 PWL(0 {v_lo} {span:.3e} {v_hi})",
               False, span / n, span)
    return v_lo + (v_hi - v_lo) * (t / span), v / V_FULL


def quasi_static(kappa_l, t_drive, v_drive):
    """Map a drive through the static transfer, which is what the fast tools do."""
    vs, ts = static_transfer(kappa_l)
    order = np.argsort(vs)
    return np.interp(v_drive, vs[order], ts[order])


def ringdown(t, dev):
    """Beat period and amplitude decay of a cavity's transient, from the trace.

    `dev` is what the photon lifetime adds: the co-solved response minus the
    quasi-static one. It carries two periodicities. The one under test is the
    cavity beating against the detuned drive. The other is the round trip
    itself, because a ring is a sampled system and emits a copy every `T_rt`.

    A boxcar of exactly one round trip removes the second — the comb's lines sit
    at multiples of `1/T_rt`, which is where that filter's nulls are — and leaves
    the envelope. The extrema of what remains give the period from their spacing
    and the decay from their heights.

    Returned as (period, tau_amplitude). Both are predicted by closed forms that
    appear nowhere in any device, which is why they are worth measuring.
    """
    h = t[1] - t[0]
    w = max(int(round(T_RT / h)), 1)
    sm = np.convolve(dev, np.ones(w) / w, mode="same")
    ext = [i for i in range(1, len(sm) - 1)
           if (sm[i] - sm[i - 1]) * (sm[i + 1] - sm[i]) < 0]
    # The ringdown starts once the cavity has dumped, so the fit starts after
    # the largest excursion. What comes before it is the cavity *filling* — a
    # different process, and including it biases the decay long by a quarter
    # (50.4 ps against a true 40.8 on the headline case).
    top = max(ext, key=lambda i: abs(sm[i]))
    peaks = [(t[i], abs(sm[i])) for i in ext if t[i] > t[top]]
    # Drop the tail, where the deviation is at the solver's floor rather than
    # the cavity's: fitting noise would bias it the other way.
    peaks = [p for p in peaks if p[1] > 1e-3 * peaks[0][1]]
    if len(peaks) < 3:
        return float("nan"), float("nan")
    period = 2.0 * float(np.mean(np.diff([p[0] for p in peaks])))
    slope = np.polyfit([p[0] for p in peaks], np.log([p[1] for p in peaks]), 1)[0]
    return period, -1.0 / slope


def prbs(order: int, n: int) -> list:
    """A maximal-length shift register, so the pattern is not a square wave."""
    taps = {7: (7, 6), 9: (9, 5), 11: (11, 9)}[order]
    reg = [1] * order
    out = []
    for _ in range(n):
        bit = reg[taps[0] - 1] ^ reg[taps[1] - 1]
        reg = [bit] + reg[:-1]
        out.append(bit)
    return out


def pair_prbs(kappa_l, rate, n_bits=48):
    """The same PRBS through both models, on one time grid.

    Edges take a tenth of a bit, so the drive itself is not what limits either
    answer — the difference between the traces is the cavity and nothing else.
    """
    bits = prbs(7, n_bits)
    bit_s = 1.0 / rate
    edge = bit_s / 10.0
    pts = ["0 0"]
    for k, b in enumerate(bits):
        t0 = k * bit_s
        pts.append(f"{t0:.6e} {V_STEP * (bits[k - 1] if k else 0)}")
        pts.append(f"{t0 + edge:.6e} {V_STEP * b}")
    stop = n_bits * bit_s
    drive = "Vmod vmod 0 PWL(" + " ".join(pts) + f" {stop:.6e} {V_STEP * bits[-1]})"
    h = min(T_RT / 3.0, bit_s / 60.0)
    tc, vc = run(kappa_l, drive, True, h, stop)
    ts, vs = run(kappa_l, drive, False, h, stop)
    n = min(len(tc), len(ts))
    return tc[:n], vc[:n] / V_FULL, ts[:n], vs[:n] / V_FULL, bits


def eye_opening(t, v, bit_s, bits, skip=4):
    """Vertical eye opening at the best sampling phase, as a transmission.

    The margin a decision circuit actually has: the worst `1` minus the best
    `0`, maximised over where in the bit you sample. Negative means the eye is
    closed and no threshold separates the rails.

    This is bounded and it is what a designer reads, which an RMS error is not —
    a cavity dumping overshoots to several times the swing, so an RMS metric
    saturates long before the eye does and stops ranking anything.
    """
    best = -1e9
    for frac in np.linspace(0.05, 0.95, 19):
        ones, zeros = [], []
        for k in range(skip, len(bits)):
            s = np.interp((k + frac) * bit_s, t, v)
            (ones if bits[k] else zeros).append(s)
        if ones and zeros:
            best = max(best, min(ones) - max(zeros))
    return best


def eye_segments(t, v, bit_s, skip=4):
    """Fold the trace onto one bit period, dropping the pattern's first bits."""
    out = []
    k = skip
    while (k + 2) * bit_s < t[-1]:
        m = (t >= k * bit_s) & (t < (k + 2) * bit_s)
        if m.sum() > 2:
            out.append((t[m] - k * bit_s, v[m]))
        k += 1
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true",
                    help="assert the physics instead of plotting it")
    ap.add_argument("--png", default="ring_photon_lifetime.png")
    args = ap.parse_args()

    tau = tau_amplitude(KAPPA_L)
    print(f"ring: {L_RING_UM:.1f} um round trip, FSR {FSR_HZ / 1e12:.3f} THz, "
          f"V_pi {V_PI:.1f} V")
    print(f"kappa_L {KAPPA_L}: loaded Q {loaded_q(KAPPA_L):.3g}, "
          f"amplitude tau {tau * 1e12:.2f} ps")

    # ── B: the step response, where the two models part company ──────────────
    t_step, stop = 50e-12, 400e-12
    drive = f"Vmod vmod 0 PWL(0 0 {t_step:.3e} 0 {t_step + 1e-15:.3e} {V_STEP})"
    t_co, v_co = run(KAPPA_L, drive, True, 0.1e-12, stop)
    t_qs, v_qs = run(KAPPA_L, drive, False, 0.1e-12, stop)
    v_co, v_qs = v_co / V_FULL, v_qs / V_FULL

    f_beat = detuning_hz(V_STEP)
    print(f"step {V_STEP} V -> detuning {f_beat / 1e9:.2f} GHz, "
          f"predicted beat period {1e12 / f_beat:.1f} ps")

    after = t_co > t_step
    ta, va = t_co[after], v_co[after]
    period_meas, tau_meas = ringdown(ta, va - v_qs[after])
    overshoot = va.max()
    print(f"measured: beat {period_meas * 1e12:.1f} ps (predicted "
          f"{1e12 / f_beat:.1f}), amplitude tau {tau_meas * 1e12:.1f} ps "
          f"(predicted {tau * 1e12:.1f})")
    print(f"          overshoot {overshoot:.2f} x the laser's own power, "
          f"quasi-static settles instantly at {v_qs[-1]:.3f}")

    if args.selftest:
        # The photon lifetime has two terms, `a` and `t`, and one operating
        # point cannot test both. At kappa_L = 0.2 the coupler carries 93 % of
        # the round-trip decay, so doubling `alpha_dB_cm` moves tau by 6 % and
        # a loss bug hides inside the tolerance. At kappa_L = 0.05 the split is
        # roughly even and the same bug moves tau by 35 %. Both are checked.
        for kl in (0.05, 0.10, KAPPA_L):
            closed = tau_amplitude(kl)
            span = max(12 * closed, 400e-12)
            d = (f"Vmod vmod 0 PWL(0 0 {span * 0.1:.3e} 0 "
                 f"{span * 0.1 + 1e-15:.3e} {V_STEP})")
            h = closed / 400.0
            tc, vc = run(kl, d, True, h, span)
            tq, vq = run(kl, d, False, h, span)
            # The two runs can land one point apart: a delay shortens the
            # internal step, and the last interval is clamped to `stop`.
            n = min(len(tc), len(tq))
            tc, vc, vq = tc[:n], vc[:n], vq[:n]
            m = tc > span * 0.1
            _, tm = ringdown(tc[m], ((vc - vq) / V_FULL)[m])
            print(f"  kappa_L {kl:.2f}: tau {tm * 1e12:7.1f} ps vs closed form "
                  f"{closed * 1e12:7.1f} ps  ({100 * (tm - closed) / closed:+.1f} %)")
            assert abs(tm - closed) / closed < 0.10, (
                f"kappa_L {kl}: tau {tm * 1e12:.1f} ps vs {closed * 1e12:.1f} ps")

        # The cavity empties over its photon lifetime. This is the claim.
        assert abs(tau_meas - tau) / tau < 0.10, (
            f"tau {tau_meas * 1e12:.1f} ps vs closed form {tau * 1e12:.1f} ps")
        # And it rings at the detuning it was given.
        assert abs(period_meas - 1.0 / f_beat) / (1.0 / f_beat) < 0.10, (
            f"beat {period_meas * 1e12:.1f} ps vs predicted {1e12 / f_beat:.1f} ps")
        # The stored field leaves through the through port, so it transiently
        # carries more than the laser put in. A static model cannot exceed 1.
        assert overshoot > 1.5, f"overshoot {overshoot:.2f} is not a cavity dumping"
        assert v_qs.max() <= 1.001, f"quasi-static reached {v_qs.max():.3f} > 1"
        # Both models agree on the endpoints, which is why the disagreement in
        # between is dynamics rather than a different circuit.
        assert abs(v_co[0] - v_qs[0]) < 1e-3, "models disagree before the step"
        assert abs(v_co[-1] - v_qs[-1]) < 5e-3, "models disagree once settled"
        print("selftest: ok")
        return 0

    return plot(args.png, t_co, v_co, t_qs, v_qs, t_step, f_beat, tau)


def plot(png, t_co, v_co, t_qs, v_qs, t_step, f_beat, tau):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("A ring modulator driven faster than its cavity can respond",
                 fontsize=13)

    # A — the static transfer both models share
    vs, ts = static_transfer(KAPPA_L)
    ax[0, 0].plot(vs, ts, color="#1f6feb")
    ax[0, 0].axvline(0.0, color="#999", lw=0.8, ls=":")
    ax[0, 0].axvline(V_STEP, color="#999", lw=0.8, ls=":")
    ax[0, 0].set(xlabel="bias (V)", ylabel="transmission",
                 title="A — the steady-state transfer, identical in both models")
    ax[0, 0].grid(alpha=0.3)

    # B — the step response
    ax[0, 1].plot(t_qs * 1e12, v_qs, color="#999", lw=2,
                  label="static resonance (what S-matrix tools give)")
    ax[0, 1].plot(t_co * 1e12, v_co, color="#d1242f", lw=1.4,
                  label="co-solved, with the photon lifetime")
    after_step = t_co > t_step
    t_first = t_co[after_step][int(np.argmax(v_co[after_step]))]
    for k in range(3):
        ax[0, 1].axvline((t_first + k / f_beat) * 1e12, color="#1f6feb",
                         lw=0.8, ls="--", alpha=0.7)
    ax[0, 1].axhline(1.0, color="#333", lw=0.8, ls=":")
    ax[0, 1].text(0.98, 0.95, f"dashed: 1/detuning = {1e12 / f_beat:.0f} ps\n"
                              "dotted: the laser's own power",
                  transform=ax[0, 1].transAxes, ha="right", va="top", fontsize=8)
    ax[0, 1].set(xlabel="time (ps)", ylabel="transmission",
                 title="B — one 0.5 V step, two answers")
    ax[0, 1].legend(fontsize=8, loc="center right")
    ax[0, 1].grid(alpha=0.3)

    # C — the eye, at a bit rate the static model says is fine
    rate = 25e9
    t_c, v_c, t_s, v_s, bits_c = pair_prbs(KAPPA_L, rate, n_bits=64)
    bit_s = 1.0 / rate
    for tt, vv, colour, name, z in ((t_c, v_c, "#d1242f", "co-solved", 1),
                                    (t_s, v_s, "#444", "static resonance", 2)):
        for seg_t, seg_v in eye_segments(tt, vv, bit_s):
            ax[1, 0].plot(seg_t * 1e12, seg_v, color=colour, lw=0.5,
                          alpha=0.35, zorder=z)
        ax[1, 0].plot([], [], color=colour, lw=1.5, label=name)
    o_c = eye_opening(t_c, v_c, bit_s, bits_c)
    o_s = eye_opening(t_s, v_s, bit_s, bits_c)
    ax[1, 0].set(xlabel="time within a bit (ps)", ylabel="transmission",
                 title=f"C — the eye at {rate / 1e9:.0f} Gb/s: "
                       f"static says {o_s:.2f}, the cavity gives {o_c:.2f}")
    ax[1, 0].legend(fontsize=8, loc="upper right")
    ax[1, 0].grid(alpha=0.3)

    # D — where the static model is still usable
    # A ring coupled harder than this has almost no static eye to begin with,
    # so a ratio against it divides by nearly nothing and says more about the
    # denominator than about the cavity.
    kls = [0.05, 0.08, 0.12, 0.20, 0.32]
    rates = [10e9, 25e9, 50e9, 100e9]
    grid = np.zeros((len(kls), len(rates)))
    for i, kl in enumerate(kls):
        for j, r in enumerate(rates):
            tc, vc, ts, vs, bb = pair_prbs(kl, r, n_bits=40)
            o_s = eye_opening(ts, vs, 1.0 / r, bb)
            o_c = eye_opening(tc, vc, 1.0 / r, bb)
            # 100 = the fast model is right. Below = it promised an eye that is
            # not there. Above = the cavity's overshoot opens the rails wider
            # than it predicted, which is also a wrong answer, just a flattering
            # one. Clipped so one extreme cell does not set the whole scale.
            grid[i, j] = min(200.0, 100.0 * max(o_c, 0.0) / max(o_s, 1e-9))
    im = ax[1, 1].imshow(grid, origin="lower", aspect="auto", cmap="RdBu_r",
                         vmin=0, vmax=200,
                         extent=(-0.5, len(rates) - 0.5, -0.5, len(kls) - 0.5))
    ax[1, 1].set_xticks(range(len(rates)),
                        [f"{r / 1e9:.0f}" for r in rates])
    ax[1, 1].set_yticks(range(len(kls)),
                        [f"{tau_amplitude(k) * 1e12:.0f}" for k in kls])
    for i in range(len(kls)):
        for j in range(len(rates)):
            ax[1, 1].text(j, i, f"{grid[i, j]:.0f}", ha="center", va="center",
                          fontsize=8,
                          color="w" if grid[i, j] < 35 or grid[i, j] > 165 else "k")
    fig.colorbar(im, ax=ax[1, 1],
                 label="real eye, as % of the static model's prediction")
    ax[1, 1].set(xlabel="bit rate (Gb/s)", ylabel="photon lifetime (ps)",
                 title="D — 100 means the fast model is right")

    fig.tight_layout()
    out = png if os.path.isabs(png) else os.path.join(os.path.dirname(__file__), png)
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
