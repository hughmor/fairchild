#!/usr/bin/env python3
"""
WDM weight bank with `fc_optical_2x2` — a behavioural stand-in for a cascade of
ring modulators sharing a through bus and a drop bus.

This is the block you reach for when you care about what comes *after* the
weights: the balanced photodetector pair, the O/E/O nonlinearity, the receiver
electronics. One instance replaces N rings, so the weights become N numbers you
set directly instead of N rings' worth of resonance, coupling, and bias
parameters — and because there is no resonance left, the transient timestep is
set by your electronics rather than by an optical round trip.

```
   λ0…λ3 ──►┌─────────────────┐──► thru ──►┐
            │  fc_optical_2x2 │            ├─► balanced PD pair ──► V_out
   (dark) ──►│  w_k per λ      │──► drop ──►┘   (drop − thru)
            └────────┬────────┘
                     │  wctl_0 … wctl_3   (one control wire per channel)
```

The weight is defined so a balanced pair reads it directly:
`P_drop − P_thru = w · P_in`, with `w = −1` all-through, `+1` all-drop.

Two things this demonstrates that a static netlist can't:
  1. **Per-channel control.** `.electrical_port wctl N` makes one control wire
     per wavelength, so the netlist scales with N by changing one number.
  2. **Weights that move during a transient** — `dw_dv` ties each weight to its
     own control voltage, so a `.tran` can sweep a weight while the optical
     signal is flowing. `set_param` only works between runs.

Run:      .venv/bin/python examples/photonic/native_weight_bank.py
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

N_CH = 4
LAMBDAS_NM = [1549.0, 1550.0, 1551.0, 1552.0]
P_IN_MW = 1.0
RESPONSIVITY = 0.8      # A/W
R_TIA = 1e3             # Ω — transimpedance across the balanced pair


def netlist(n=N_CH, *, dw_dv=0.0, w0=0.0, tau_s=0.0, tran=None) -> str:
    """Weight bank + balanced PD pair. `tran` is a PWL/PULSE spec for wctl."""
    lines = [f"* WDM weight bank, {n} channels"]
    for port in ("bus", "dark", "thru", "drop"):
        lines.append(f".optical_port {port} {n}")
    # One control wire per wavelength channel. The parser rejects a width here
    # that disagrees with the optical ports, since a bundle-aware device can't
    # mean anything sensible when they differ.
    lines.append(f".electrical_port wctl {n}")

    for k, wl in enumerate(LAMBDAS_NM[:n]):
        lines.append(f"Xl{k} bus_re_{k} bus_im_{k} bus_wl_{k} fc_cw_laser "
                     f"power_mW={P_IN_MW} wavelength_nm={wl}")
        # Second input port unused — it is an input, so its wires need drivers.
        lines += [f"Vdr{k} dark_re_{k} 0 DC 0",
                  f"Vdi{k} dark_im_{k} 0 DC 0",
                  f"Vdw{k} dark_wl_{k} 0 DC {wl * 1e-9:.6e}"]
        drive = tran if tran else "DC 0"
        lines.append(f"Vc{k} wctl_{k} 0 {drive}")

    lines.append(f"Xwb bus dark thru drop wctl 0 fc_optical_2x2 "
                 f"w={w0} dw_dv={dw_dv} tau_s={tau_s}")

    # Balanced pair: both PDs push photocurrent into the same node with opposite
    # sign, so V(bal) ∝ (P_drop − P_thru) = Σ w_k·P_k — the weighted sum.
    for k in range(n):
        lines.append(f"Xpd_d{k} drop_re_{k} drop_im_{k} drop_wl_{k} bal 0 "
                     f"fc_photodetector responsivity={RESPONSIVITY} "
                     f"i_dark_a=0 r_shunt=1e12")
        lines.append(f"Xpd_t{k} thru_re_{k} thru_im_{k} thru_wl_{k} 0 bal "
                     f"fc_photodetector responsivity={RESPONSIVITY} "
                     f"i_dark_a=0 r_shunt=1e12")
    lines.append(f"Rtia bal 0 {R_TIA}")
    return "\n".join(lines) + "\n"


def optical_power(res, port: str, k: int) -> float:
    re = float(res[f"V({port}_re_{k})"][0])
    im = float(res[f"V({port}_im_{k})"][0])
    return re * re + im * im


# ── 1. static weights, set per channel ───────────────────────────────────────
def static_weights(weights) -> tuple[np.ndarray, float]:
    """Return (measured per-channel weights, balanced-PD voltage)."""
    c = fc.Circuit()
    c.load_str(netlist() + ".op\n")
    for k, w in enumerate(weights):
        c.set_param("Xwb", f"w_{k}", float(w))
    r = c.run("op")
    meas = np.array([
        (optical_power(r, "drop", k) - optical_power(r, "thru", k)) / (P_IN_MW * 1e-3)
        for k in range(N_CH)
    ])
    return meas, float(r["V(bal)"][0])


# ── 2. weights swept by their control voltages, in one transient ─────────────
def dynamic_weights():
    """Ramp every channel's control voltage on a different schedule and watch
    the balanced output follow the weighted sum in real time."""
    # w_k = 0 + 1.0·V(wctl_k); each channel gets its own PWL ramp so the
    # weighted sum is a moving target rather than a scaled copy of one signal.
    c = fc.Circuit()
    c.load_str(netlist(dw_dv=1.0, w0=0.0,
                       tran="PWL(0 -1 20n 1 40n -1)") + ".tran 200p 40n\n")
    # Stagger the channels: channel k's ramp is delayed by k·4 ns, so the
    # weighted sum is a moving target rather than a scaled copy of one signal.
    WS = getattr(fc, "WaveformSource", None) or fc.fairchild.WaveformSource
    for k in range(N_CH):
        t0 = k * 4e-9
        t = np.array([0.0, t0, t0 + 12e-9, t0 + 24e-9, 40e-9])
        v = np.array([-1.0, -1.0, 1.0, -1.0, -1.0])
        c.set_source(f"Vc{k}", WS(t, v))
    r = c.run("tran", step=200e-12, stop=40e-9)
    t = np.asarray(r.time())
    v_bal = np.asarray(r["V(bal)"])
    ctl = np.array([np.asarray(r[f"V(wctl_{k})"]) for k in range(N_CH)])
    return t, v_bal, ctl


def expected_bal(weights) -> float:
    """V(bal) = R·R_λ·Σ w_k·P_k for a balanced pair into a transimpedance."""
    return R_TIA * RESPONSIVITY * float(np.sum(weights)) * P_IN_MW * 1e-3


def selftest() -> int:
    # Bipolar weights land exactly, and the balanced pair reads their sum.
    want = np.array([1.0, -1.0, 0.0, 0.5])
    meas, v_bal = static_weights(want)
    assert np.allclose(meas, want, atol=1e-9), f"weights {meas} != {want}"
    want_v = expected_bal(want)
    # Relative: Newton's convergence tolerance scales with the node voltage,
    # and V(bal) reaches volts when every channel is fully on one port.
    assert np.isclose(v_bal, want_v, rtol=1e-8, atol=1e-12), \
        f"V(bal)={v_bal} expected {want_v}"

    # All-through and all-drop are the rails, symmetric about zero.
    lo, _ = static_weights([-1.0] * N_CH)
    hi, _ = static_weights([1.0] * N_CH)
    assert np.allclose(lo, -1.0) and np.allclose(hi, 1.0), f"rails {lo} {hi}"

    # A weight beyond the rails clamps rather than producing gain.
    over, v_over = static_weights([4.0] * N_CH)
    assert np.allclose(over, 1.0, atol=1e-9), f"clamp failed: {over}"
    assert np.isclose(v_over, expected_bal([1.0] * N_CH), rtol=1e-8, atol=1e-12)

    # Dynamic: the balanced output must track R·R_λ·P·Σ w_k(t) with w_k = V_ctl,
    # clamped — that is the whole claim of the voltage-controlled path.
    t, v_bal_t, ctl = dynamic_weights()
    pred = R_TIA * RESPONSIVITY * P_IN_MW * 1e-3 * np.clip(ctl, -1, 1).sum(axis=0)
    err = np.max(np.abs(v_bal_t - pred))
    assert err < 1e-7, f"transient mismatch {err:.3e} V"
    assert np.ptp(v_bal_t) > 0.5 * abs(expected_bal([1.0] * N_CH)), \
        "control ramps produced no meaningful swing"

    print(f"selftest OK — {N_CH} channels, weights exact to 1e-9, "
          f"V(bal) rails ±{abs(expected_bal([1.0] * N_CH)) * 1e3:.1f} mV, "
          f"transient tracks Σw_k(t) to {err:.1e} V")
    return 0


def plot(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(11, 7))

    # (0,0) weight transfer curve, per channel, swept one channel at a time.
    sweep = np.linspace(-1.4, 1.4, 29)
    for k in range(N_CH):
        got = []
        for w in sweep:
            weights = np.zeros(N_CH)
            weights[k] = w
            got.append(static_weights(weights)[0][k])
        ax[0, 0].plot(sweep, got, lw=1.2, label=f"λ{k} = {LAMBDAS_NM[k]:.0f} nm")
    ax[0, 0].plot(sweep, np.clip(sweep, -1, 1), "k--", lw=0.8, label="clamp")
    ax[0, 0].set_xlabel("requested weight")
    ax[0, 0].set_ylabel("measured (P_drop − P_thru) / P_in")
    ax[0, 0].set_title("Per-channel weight is exact, and clamps at ±1")
    ax[0, 0].legend(fontsize=7)
    ax[0, 0].grid(alpha=0.3)

    # (0,1) power split vs weight — passivity: the two ports always sum to P_in.
    c = fc.Circuit()
    c.load_str(netlist() + ".op\n")
    thru, drop = [], []
    for w in sweep:
        c.set_param("Xwb", "w_0", float(w))
        r = c.run("op")
        thru.append(optical_power(r, "thru", 0) / (P_IN_MW * 1e-3))
        drop.append(optical_power(r, "drop", 0) / (P_IN_MW * 1e-3))
    thru, drop = np.array(thru), np.array(drop)
    ax[0, 1].plot(sweep, thru, label="thru")
    ax[0, 1].plot(sweep, drop, label="drop")
    ax[0, 1].plot(sweep, thru + drop, "k--", lw=0.8, label="sum (passivity)")
    ax[0, 1].set_xlabel("weight")
    ax[0, 1].set_ylabel("fraction of P_in")
    ax[0, 1].set_title("Lossless split — no weight can make gain")
    ax[0, 1].legend(fontsize=8)
    ax[0, 1].grid(alpha=0.3)

    # (1,0)+(1,1) dynamic: control ramps and the balanced sum they produce.
    t, v_bal, ctl = dynamic_weights()
    for k in range(N_CH):
        ax[1, 0].plot(t * 1e9, np.clip(ctl[k], -1, 1), lw=1.1, label=f"w_{k}(t)")
    ax[1, 0].set_xlabel("time (ns)")
    ax[1, 0].set_ylabel("weight")
    ax[1, 0].set_title("Weights swept by their own control wires, in one .tran")
    ax[1, 0].legend(fontsize=7, ncol=2)
    ax[1, 0].grid(alpha=0.3)

    pred = R_TIA * RESPONSIVITY * P_IN_MW * 1e-3 * np.clip(ctl, -1, 1).sum(axis=0)
    ax[1, 1].plot(t * 1e9, v_bal * 1e3, lw=1.4, label="V(bal) simulated")
    ax[1, 1].plot(t * 1e9, pred * 1e3, "k--", lw=0.8, label="R·ℜ·P·Σ w_k(t)")
    ax[1, 1].set_xlabel("time (ns)")
    ax[1, 1].set_ylabel("balanced output (mV)")
    ax[1, 1].set_title("Balanced PD pair reads the weighted sum")
    ax[1, 1].legend(fontsize=8)
    ax[1, 1].grid(alpha=0.3)

    fig.suptitle("WDM weight bank — fc_optical_2x2, one instance, "
                 f"{N_CH} wavelength channels", fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    print(f"wrote {out}")


def main() -> int:
    if "--selftest" in sys.argv:
        return selftest()
    selftest()
    plot(HERE / "native_weight_bank.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
