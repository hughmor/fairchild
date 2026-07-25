#!/usr/bin/env python3
"""fit_transient.py — time-domain parameter fitting for a ring modulator.

The transient counterpart to ringfit.py (which fits CW transmission
spectra). Here we drive the modulator with a measured numpy waveform and fit
device params to a measured time-domain photodetector trace — the workflow for
calibrating against a lightlab AWG-drive + off-chip-PD-on-scope capture (both
sides arrive as numpy arrays).

Pipeline per objective evaluation:
    set_param(fit params) → set_source(measured drive into Vmod)
      → run("tran", waveguide_delay=ON, sub-round-trip step)
      → simulated PD trace
      → measurement model: 1-pole PD/TIA bandwidth, time-lag alignment,
        and gain/offset projected out in closed form (nuisance params)
      → residual vs measured trace

Two requirements proven upstream (ring_dynamics_check.py) and enforced here:
  • waveguide_delay MUST be ON (else the ring's photon-lifetime lag is invisible
    and the optimizer misattributes it to carrier τ / RC),
  • the timestep MUST resolve the round trip (τ_g = n_g·L/c), else same problem.

Validation without chip data: `--selftest` runs a synthetic parameter-recovery
— simulate with known params, distort like a real instrument (bandwidth, gain,
offset, lag, noise), then fit and check the params come back. This also probes
identifiability (which params the time-domain data can actually constrain).

Run:  .venv/bin/python examples/kicad_photonics/fit_transient.py --selftest
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, least_squares

import fairchild as fc

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
RING = REPO / "examples" / "photonic" / "native_mrr_modulator.sp"

# Ring round trip → the timestep ceiling and the delay requirement.
N_G, L_RING, C = 4.2, 500e-6, 2.99792458e8
T_RT = N_G * L_RING / C            # ~7 ps
DEFAULT_STEP = T_RT / 8.0          # resolve the round trip


def _quiet(nl: str) -> str:
    return "\n".join(l.replace(" wavelength_nm=1550", "") if "fc_waveguide" in l else l
                     for l in nl.splitlines())


NETLIST = _quiet(RING.read_text())
WS = getattr(fc, "WaveformSource", None) or fc.fairchild.WaveformSource


# ── parameter spec (mirrors ringfit.ParamSpec, standalone) ──────────────
@dataclass
class FitParam:
    name: str          # "<element>.<param>", e.g. "Xpn.V_pi_L"
    init: float
    lo: float
    hi: float
    free: bool = True


# ── the simulator forward model ──────────────────────────────────────────────
def simulate_pd(overrides: dict[str, float], drive_t, drive_v, *,
                step=DEFAULT_STEP, waveguide_delay=True, observable="V(pd_anode)",
                netlist=NETLIST):
    """Drive Vmod with (drive_t, drive_v); return (t, observable).

    observable is either a raw signal name ("V(pd_anode)", "I(...)") or an
    optical-power readout "P(<port>)" = |A|² = re²+im² for an optical bundle.
    `netlist` lets a caller swap the device under test (e.g. a different PN model).
    """
    ckt = fc.Circuit()
    ckt.load_str(netlist)
    for key, val in overrides.items():
        el, _, p = key.partition(".")
        ckt.set_param(el, p, float(val))
    ckt.set_source("Vmod", WS(np.asarray(drive_t, float), np.asarray(drive_v, float)))
    stop = float(drive_t[-1])
    r = ckt.run("tran", step=step, stop=stop, variable_step=False,
                method="gear", waveguide_delay=waveguide_delay)
    t = np.asarray(r.time())
    if observable.startswith("P(") and observable.endswith(")"):
        base = observable[2:-1]
        re = np.asarray(r[f"V({base}_re_0)"]); im = np.asarray(r[f"V({base}_im_0)"])
        return t, re * re + im * im
    return t, np.asarray(r[observable])


# ── measurement model ────────────────────────────────────────────────────────
def lowpass_1pole(y, dt, tau):
    """First-order PD/TIA bandwidth: dy/dt = (x − y)/τ, forward-Euler IIR."""
    if tau <= 0:
        return y
    a = dt / (tau + dt)
    out = np.empty_like(y)
    acc = y[0]
    for i, x in enumerate(y):
        acc += a * (x - acc)
        out[i] = acc
    return out


def best_lag(sim, meas, max_lag):
    """Integer-sample lag that best aligns sim to meas (cross-correlation)."""
    s = sim - sim.mean()
    m = meas - meas.mean()
    lags = np.arange(-max_lag, max_lag + 1)
    cc = [np.dot(np.roll(s, k), m) for k in lags]
    return int(lags[int(np.argmax(cc))])


def project_gain_offset(sim, meas, gain=None):
    """meas ≈ a·sim + b. If `gain` is given (calibrated TIA×responsivity), fit
    only the DC offset b; else fit both in closed form.

    Projecting the gain out is robust but makes amplitude-like params (e.g.
    V_pi_L) degenerate with it — identifiable only through the response
    NONLINEARITY. With a calibrated gain those params become identifiable."""
    if gain is None:
        A = np.vstack([sim, np.ones_like(sim)]).T
        (a, b), *_ = np.linalg.lstsq(A, meas, rcond=None)
    else:
        a = gain
        b = float(np.mean(meas - a * sim))
    return a * sim + b, a, b


def model_trace(sim, meas, dt, *, bw_tau, max_lag, gain=None):
    """Apply instrument model to a raw sim trace so it is comparable to meas:
    bandwidth low-pass → lag-align → gain/offset (gain known or projected)."""
    y = lowpass_1pole(sim, dt, bw_tau)
    lag = best_lag(y, meas, max_lag)
    y = np.roll(y, lag)
    fit, a, b = project_gain_offset(y, meas, gain=gain)
    return fit, dict(lag=lag, gain=a, offset=b)


# ── fitting ───────────────────────────────────────────────────────────────────
def fit_transient(specs: list[FitParam], drive_t, drive_v, meas, *,
                  bw_tau, step=DEFAULT_STEP, max_lag=40, observable="V(pd_anode)",
                  gain=None, seed=0, netlist=NETLIST, verbose=True):
    """Fit free params to a measured trace. offset/lag (and optionally gain) are
    handled per-eval. Uses differential_evolution + an LM polish: resonant
    time-domain cost surfaces are spiky/multi-modal (the resonance position is
    very sensitive to phase params), so a local optimizer alone gets trapped on
    the plateau around the narrow global minimum. Pass gain=<calibrated> to make
    amplitude-like params (V_pi_L) identifiable; gain=None projects it out."""
    free = [s for s in specs if s.free]
    fixed = {s.name: s.init for s in specs if not s.free}
    lo = np.array([s.lo for s in free]); hi = np.array([s.hi for s in free])
    span = float(np.ptp(meas)) or 1.0
    ncall = [0]

    def residuals(x):
        ncall[0] += 1
        ov = dict(fixed)
        ov.update({s.name: xi for s, xi in zip(free, x)})
        t_s, sim = simulate_pd(ov, drive_t, drive_v, step=step,
                               observable=observable, netlist=netlist)
        sim = np.interp(drive_t, t_s, sim)   # solver grid → measurement grid
        dt = drive_t[1] - drive_t[0]
        fitm, _ = model_trace(sim, meas, dt, bw_tau=bw_tau, max_lag=max_lag, gain=gain)
        return (fitm - meas) / span

    def cost(x):
        r = residuals(x)
        return float(0.5 * np.dot(r, r))

    de = differential_evolution(cost, bounds=list(zip(lo, hi)), seed=seed,
                                maxiter=40, popsize=12, tol=1e-4, mutation=(0.4, 1.2),
                                polish=False)
    pol = least_squares(residuals, de.x, bounds=(lo, hi),
                        diff_step=0.03, xtol=1e-8, ftol=1e-8, verbose=0)
    best = dict(fixed); best.update({s.name: v for s, v in zip(free, pol.x)})
    if verbose:
        print(f"  DE+polish: {ncall[0]} sims, cost={pol.cost:.3e}")
    return best, pol


# ── synthetic parameter-recovery self-test ────────────────────────────────────
def _make_drive(stop=6e-10, step=DEFAULT_STEP, bit=6e-11, v_lo=0.0, v_hi=1.5, seed=1):
    """A short pseudo-random NRZ drive on Vmod (bias+signal)."""
    t = np.arange(0.0, stop + step, step)
    rng = np.random.default_rng(seed)
    nbits = int(stop / bit) + 1
    bits = rng.integers(0, 2, nbits)
    v = np.where(bits[np.clip((t / bit).astype(int), 0, nbits - 1)] > 0, v_hi, v_lo)
    return t, v


def selftest():
    print(f"Ring: t_rt={T_RT*1e12:.1f} ps, step={DEFAULT_STEP*1e12:.2f} ps "
          f"(τ_g/8), waveguide_delay=ON\n")
    drive_t, drive_v = _make_drive()

    # Observable: the circulating field (strongly modulated). The realistic
    # through-port (V(pd_anode)) is weakly modulated on this under-coupled ring,
    # so its trace barely constrains the params — an identifiability limit, not
    # a harness bug (a real modulator needs adequate through/drop extinction).
    OBS = "P(pn_in)"

    # Ground-truth device params.
    true = {"Xpn.V_pi_L": 2.0e-3, "Xpn.alpha_dB_cm": 10.0}
    t_sim, clean = simulate_pd(true, drive_t, drive_v, observable=OBS)
    # The solver may use its own grid; resample both drive and clean onto t_sim.
    meas = clean.copy()

    # Distort like a real instrument: 1-pole PD/TIA BW, gain, offset, lag, noise.
    dt = t_sim[1] - t_sim[0]
    BW_TAU = 1.0 / (2 * np.pi * 12e9)          # ~12 GHz PD+TIA
    GAIN, OFFSET, LAG = 3.0, -0.4, 7
    meas = lowpass_1pole(meas, dt, BW_TAU)
    meas = GAIN * meas + OFFSET
    meas = np.roll(meas, LAG)
    rng = np.random.default_rng(0)
    meas = meas + rng.normal(0, 0.01 * np.ptp(meas), meas.shape)
    drive_on_sim = drive_v_resampled(drive_t, drive_v, t_sim)

    def run_case(label, gain):
        specs = [
            FitParam("Xpn.V_pi_L", init=3.2e-3, lo=5e-4, hi=6e-3),
            FitParam("Xpn.alpha_dB_cm", init=5.0, lo=1.0, hi=25.0),
        ]
        print(f"\n[{label}]")
        best, _ = fit_transient(specs, t_sim, drive_on_sim, meas,
                                bw_tau=BW_TAU, step=DEFAULT_STEP, observable=OBS,
                                gain=gain)
        print("  param                 true        init        recovered     err%")
        errs = {}
        for s in specs:
            tv, rv = true[s.name], best[s.name]
            errs[s.name] = 100 * (rv - tv) / tv
            print(f"    {s.name:18s} {tv:10.4g} {s.init:10.4g} {rv:12.5g} {errs[s.name]:8.1f}")
        return errs

    # Positive control: calibrated gain (Hugh will know TIA×responsivity).
    e_cal = run_case("calibrated gain", gain=GAIN)
    # Robustness: gain unknown/projected. With a GLOBAL optimizer both params
    # still recover — the resonance-shape traversal encodes V_pi_L nonlinearly,
    # not as a pure scale, so it is not actually degenerate with gain.
    e_proj = run_case("gain projected (unknown)", gain=None)

    ok = (all(abs(e) < 10.0 for e in e_cal.values())
          and all(abs(e) < 10.0 for e in e_proj.values()))
    print(f"\nRECOVERY {'PASS' if ok else 'FAIL'} — V_pi_L and alpha within 10% "
          f"(calibrated AND projected gain)")
    print("Lessons baked in: (1) global optimization (differential_evolution) is\n"
          "  REQUIRED — the resonant cost surface is spiky/multi-modal and local\n"
          "  least-squares gets trapped on the plateau around the narrow minimum;\n"
          "  (2) the binding identifiability constraint is the OBSERVABLE — we fit\n"
          "  the strongly-modulated circulating field; a weakly-coupled through\n"
          "  port (V(pd_anode) here varies ~4%) barely constrains the params.")
    return 0 if ok else 1


def drive_v_resampled(drive_t, drive_v, t_sim):
    """Resample the NRZ drive onto the solver's output grid (zero-order hold)."""
    idx = np.clip(np.searchsorted(drive_t, t_sim, side="right") - 1, 0, len(drive_v) - 1)
    return drive_v[idx]


def main():
    ap = argparse.ArgumentParser(description="Time-domain ring-modulator fitting")
    ap.add_argument("--selftest", action="store_true",
                    help="synthetic parameter-recovery (no chip data needed)")
    args = ap.parse_args()
    if args.selftest:
        raise SystemExit(selftest())
    ap.error("real-data fitting entry point TBD; use --selftest for now")


if __name__ == "__main__":
    main()
