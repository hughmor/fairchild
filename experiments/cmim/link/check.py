#!/usr/bin/env python3
"""check.py — assert that the CMIM analogue deck does what the paper says.

Every check compares the circuit against an ABSOLUTE anchor: a closed-form
power budget, the sine transfer function of an interferometer, the Gaussian
gain fit, the discrete-time variance of a filtered white sequence.  None of
them compares one part of the deck against another, because a fault shared by
both sides of such a comparison is invisible.

Run it after any change to a model or to the deck:

    .venv/bin/python experiments/cmim/link/check.py

Each check names, in its comment, the sabotage that makes it fail.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "common"))
from link import Link  # noqa: E402

Q = 1.602176634e-19
H_PLANCK = 6.62607015e-34
C_LIGHT = 299792458.0
BW_OPT = 4.77e12          # the receiver's optical noise bandwidth, Hz
# Output noise, in volts rms, per amp per root hertz injected at the detector.
# Measured against the analytic discrete-time variance of the two backward-Euler
# pole sections at a 2 ps step, and agreeing with it to 0.6 %.
UNIT_2PS = 6.5486e-05 / 1e-12

FAILED: list[str] = []


def check(name, got, want, tol, unit=""):
    ok = abs(got - want) <= tol
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:44s} "
          f"{got:12.6g} vs {want:12.6g} {unit}")
    if not ok:
        FAILED.append(name)


def check_true(name, cond, note=""):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name:44s} {note}")
    if not cond:
        FAILED.append(name)


def power(op: dict, port: str) -> float:
    lo = {k.lower(): v for k, v in op.items()}
    return lo[f"v({port}_re_0)"] ** 2 + lo[f"v({port}_im_0)"] ** 2


def soft_clip(v, swing, n=4.0):
    u = v / swing
    return swing * u / (1.0 + abs(u) ** n) ** (1.0 / n)


def main() -> int:
    lk = Link()
    print("CMIM analogue link — checks against closed form\n")

    ports = ["las", "m1o", "m2o", "comb", "soao"]
    probes = ([f"v({p}_{q}_0)" for p in ports for q in ("re", "im")]
              + ["v(rxout)", "v(asesoa)"])
    op = lk.op(probes)
    lo = {k.lower(): v for k, v in op.items()}
    p = {q: power(op, q) for q in ports}

    # ── 1. power budget, stage by stage ──────────────────────────────────
    # Anchor: arithmetic on the paper's loss numbers.
    # Sabotage: change il_dB in va_tfln_mzm; both modulator rows fail.
    print("power budget:")
    check("laser output (mW)", p["las"] * 1e3, 4.47, 1e-3, "mW")
    er, il = 10 ** (42.0 / 10.0), 10 ** (-6.0 / 10.0)
    t_quad = il * ((1 - 1 / er) * 0.5 + 1 / er)
    check("MZM1 transmission at quadrature", p["m1o"] / p["las"], t_quad, 1e-6)
    check("MZM2 transmission at quadrature", p["m2o"] / p["m1o"], t_quad, 1e-6)
    # 1e-6, not 1e-12: `kappa_L` is stored by scaling `kappa_per_m` against a
    # fixed length, so the round trip costs about 4e-8.  That is the device, not
    # the deck's decimal expansion of pi/4.
    check("50:50 combiner split", p["comb"] / p["m2o"], 0.5, 1e-6)

    # QD SOA gain: the Gaussian fit evaluated by hand at 20 C and 1310 nm,
    # divided by the input-saturation term.  Sabotage: change bw_3dB_nm.
    g_pk, lam_pk, bw = 18.2, 1300.0, 27.3
    two_sig2 = (bw / 2) ** 2 / np.log(g_pk / (g_pk - 3.0))
    g0_dB = g_pk * np.exp(-((1310.0 - lam_pk) ** 2) / two_sig2)
    sat = 1.0 + p["comb"] / (1e-3 * 10 ** (4.1 / 10))
    g_dB = g0_dB - 10 * np.log10(sat)
    check("QD SOA gain (dB)", 10 * np.log10(p["soao"] / p["comb"]), g_dB, 2e-3, "dB")
    check("QD SOA input power (dBm)",
          10 * np.log10(p["comb"] / 1e-3), -14.5, 0.3, "dBm")

    # The exported ASE density, against Supplementary equation 20 by hand.
    # Sabotage: change nf_dB or f_pol; this row moves and nothing else does.
    hnu = H_PLANCK * C_LIGHT / 1310e-9
    check("exported ASE density (W/Hz)", lo["v(asesoa)"],
          10 ** (5.6 / 10) * 10 ** (g_dB / 10) * hnu / 2.0, 1e-21, "W/Hz")

    # Receiver: responsivity times power times transimpedance, through the soft
    # clip and the output divider.  Sabotage: change clip_n and this row alone
    # moves, by 0.03 %, which is why the tolerance is tight.
    i_ph = 0.8 * (p["soao"] + lo["v(asesoa)"] * BW_OPT) + 10e-9
    check("receiver output (V)", lo["v(rxout)"],
          soft_clip(300.0 * i_ph, 2.0) * 1e6 / (1e6 + 50.0), 1e-6, "V")

    # ── 2. the modulator transfer really is a bounded sine ───────────────
    # This is the nonlinearity the algorithm depends on.  Anchor: the analytic
    # transfer function.  Sabotage: move phi_bias_deg off -90 and the residual
    # explodes.  a_v = 1 with a wide swing keeps the driver out of its own clip,
    # so what is measured is the modulator.  The driver's 50 ohm source into the
    # modulator's 50 ohm termination halves the drive, so the half-wave voltage
    # referred to the AWG is 2*V_pi/a_v = 3.0 V.
    print("\nspin-channel transfer function:")
    vin = np.linspace(-3.0, 3.0, 41)
    t_meas = np.array([power(lk.op(probes, v_x=v, a_v_x=1.0, sw_x=200.0), "m1o")
                       for v in vin]) / p["las"]
    t_model = il * ((1 - 1 / er)
                    * (1 + np.cos(-np.pi / 2 + np.pi * vin / 3.0)) / 2 + 1 / er)
    check("sine transfer, worst residual",
          float(np.max(np.abs(t_meas - t_model))), 0.0, 1e-6)
    check("half-wave voltage at the AWG (V)",
          float(vin[int(np.argmax(t_meas))] - vin[int(np.argmin(t_meas))]),
          3.0, 0.16, "V")
    check("extinction ratio (dB)",
          float(10 * np.log10(t_meas.max() / t_meas.min())), 42.0, 0.6, "dB")

    # ── 3. bandwidth comes out of the poles, not out of a parameter ──────
    # Anchor: the analytic cascade of every pole ON THE SPIN PATH.  There are
    # five, and which five is the point of this check:
    #
    #   65.0 GHz   the driver, which stands in for the AWG's analogue bandwidth
    #  127.3 GHz   MZM1's electrode capacitance against the driver's 50 ohm
    #              source in parallel with its own 50 ohm termination.  Nobody
    #              writes this one down and it costs 0.6 dB at 111 GHz.
    #  110.0 GHz   MZM1's electro-optic pole
    #  100.0 GHz   the photodetector
    #   70.0 GHz   the transimpedance amplifier
    #
    # MZM2's electro-optic pole is NOT here, and that is physics rather than an
    # omission: with its drive held at DC the second modulator is a constant
    # attenuator, and its bandwidth applies to its own electrode, not to light
    # passing through it.  Adding it moves the fit from 0.03 dB to 7.6 dB.
    # Sabotage: change f_eo, f_tia or c_elec and the measured curve follows.
    print("\nsmall-signal response:")
    f, h = lk.ac(1e8, 3e11, 30)
    mag = 20 * np.log10(np.abs(h) / np.abs(h[0]))
    f_elec = 1.0 / (2 * np.pi * 25.0 * 50e-15)
    poles = [65e9, f_elec, 110e9, 100e9, 70e9]
    ideal = -10 * np.sum([np.log10(1 + (f / fp) ** 2) for fp in poles], axis=0)
    check("cascade matches the declared poles (dB)",
          float(np.max(np.abs(mag - ideal)[f < 1.5e11])), 0.0, 0.05, "dB")
    f3 = float(np.interp(-3.0, mag[::-1], f[::-1]))
    f6 = float(np.interp(-6.0, mag[::-1], f[::-1]))
    print(f"         raw link 3 dB {f3/1e9:.1f} GHz, 6 dB {f6/1e9:.1f} GHz."
          f"  Paper, WITH pre-emphasis: 55 and 75 GHz (S1.5.1).")

    # ── 4. noise has the right size and the right dependence ─────────────
    # Anchor: the analytic discrete-time variance of the injected sequence
    # through the two backward-Euler pole sections.
    # Sabotage: drop the sig_spont term and the measured rows fall by 8x.
    print("\nnoise:")
    n, dt = 8192, 2e-12
    t = np.arange(n) * dt
    z = np.zeros(n)

    _, v_off = lk.tran(t, z, z, tstep=dt, trannoise=False)
    check("no noise when trannoise is off (V)",
          float(np.std(v_off[n // 4:])), 0.0, 1e-9, "V")

    stds = []
    for i_mA in (0, 60, 100):
        npar = lk.noise_params(i_noise_mA=i_mA, p_ase_ref=-30.0)
        _, v = lk.tran(t, z, z, tstep=dt, trannoise=True,
                       i_noise_mA=i_mA, p_ase_ref=-30.0, **npar)
        got = float(np.std(v[n // 4:]))
        s = (2 * Q * npar["i_ph_n"]
             + 4 * 0.64 * npar["p_sig_n"] * npar["rho_ase_n"]
             + 2 * 0.64 * npar["rho_ase_n"] ** 2 * BW_OPT
             + 3e-12 ** 2)
        want = UNIT_2PS * np.sqrt(s)
        stds.append(got)
        # One transient is one realisation.  The docs put several per cent of
        # scatter between seeds, so 8 % is the band, not the expectation.
        check(f"noise at {i_mA:3d} mA bulk SOA (mV)", got * 1e3, want * 1e3,
              0.08 * want * 1e3, "mV")
    check_true("bulk SOA bias raises the noise monotonically",
               stds[0] < stds[1] < stds[2],
               f"{stds[0]*1e3:.2f} < {stds[1]*1e3:.2f} < {stds[2]*1e3:.2f} mV")

    print()
    if FAILED:
        print(f"{len(FAILED)} FAILED: " + ", ".join(FAILED))
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
