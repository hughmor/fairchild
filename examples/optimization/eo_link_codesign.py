#!/usr/bin/env python3
"""
Co-design a 10 Gb/s electro-optic link: modulator length and receiver load, together.

The electrical and the optical halves of a link are usually optimised by two
people in two tools, trading spreadsheets across the boundary.  They should not
be, because the trade-off *is* the boundary:

      laser ──► MZM ──► photodiode ──┬── v(pout)
                 ▲                   │
              [ Cmod ]              Rl ∥ Cl
                 ▲
      Vdrv ──[ Rd ]

  * A **longer modulator** is more efficient — V_pi falls as 1/L — so the same
    driver swing buys more optical extinction.
  * A longer modulator is also **slower**: its capacitance rises as L, and
    `Rd·Cmod` starts eating the bit period, so less of that swing arrives.
  * A **larger receiver load** turns more photocurrent into volts, and also
    slows the receiver through `Rl·Cl`.

Each of those is a real optimum in the interior, and neither can be found from
one side of the boundary alone.  Here they are found together, from a single
gradient that crosses the domain twice — driver → optics → detector → load.

**The design variables are not netlist parameters, and that is the point.**  The
simulator differentiates `V_pi`, `Cmod` and `Rl`; the *design* has one knob, `L`,
which moves `V_pi` and `Cmod` together because they are two consequences of one
length.  Writing that relation in JAX and calling `jax.grad` composes it with the
adjoint automatically — the chain rule from `dL/d(netlist)` to `dL/d(design)` is
never written down.  That is what the `jax.custom_vjp` adapter is for; do it by
hand with `run.backward()` and numpy if you prefer (see `step_response_match.py`).

Bounds matter here, and physically so.  A lumped MZM whose V_pi falls below the
drive swing wraps onto the next fringe of its transfer function, where the
objective has spurious optima that are real features of the model and not of any
useful design.  The length is therefore bounded, smoothly, by optimising through
a sigmoid — so the gradient stays defined everywhere and the optimiser cannot
walk out of the physics.

Requires JAX:  uv run --with jax python examples/optimization/eo_link_codesign.py

Run:      .venv/bin/python examples/optimization/eo_link_codesign.py
Selftest: same, with --selftest (asserts the optimum + gradient, no plotting)
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np

try:
    import fairchild as fc
except ImportError as e:  # pragma: no cover - environment guard
    sys.exit(f"fairchild Python package not installed: {e}\n"
             "Build with: maturin develop --release -m crates/fairchild-py/Cargo.toml")

HERE = pathlib.Path(__file__).resolve().parent

TB = 100e-12          # 10 Gb/s
STEP, STOP = 2e-12, 3.5 * TB

NETLIST = """* 10 Gb/s electro-optic link
.optical_port in0
.optical_port out0
Xl0 in0 fc_cw_laser power_mW=1.0 wavelength_nm=1550
Vdrv drv 0 PULSE(0.5 2.5 0 5p 5p 45p 100p)
Rd drv vsig 50
Cmod vsig 0 300f
Xmzm in0 out0 vsig 0 fc_mzm V_pi=3.0 alpha=1.0 e_r=1000
Xpd out0 pout 0 fc_photodetector responsivity=0.8
Rl pout 0 900
Cl pout 0 30f
.tran 2p 400p
.end
"""

# Modulator scaling: both of these follow from one length, which is exactly why
# the two cannot be optimised separately.
L_REF, VPI_REF = 1000.0, 3.0        # µm, V   — V_pi = VPI_REF·L_REF/L
C_PER_UM = 0.3e-15                  # F/µm    — Cmod = C_PER_UM·L

L_MIN, L_MAX = 300.0, 1600.0        # µm  — upper bound keeps V_pi above the drive swing
R_MIN, R_MAX = 150.0, 3000.0        # Ω

PARAMS = ["Xmzm.V_pi", "Cmod.c", "Rl.r"]
# A deliberately poor starting design: a short, inefficient modulator into a
# low-impedance load — fast, but with almost no signal to be fast about.
START = np.array([-1.5, -2.0])      # sigmoid → L ≈ 540 um, Rl ≈ 490 Ohm


def netlist_values(design):
    """The three netlist parameters implied by a design `(L_um, Rl)`."""
    length, load = design
    return {"Xmzm.V_pi": VPI_REF * L_REF / length,
            "Cmod.c": C_PER_UM * length,
            "Rl.r": load}


def simulate(ckt, design):
    return ckt.tran_adjoint({"v": "pout"}, step=STEP, stop=STOP,
                            method="tr", reltol=1e-11, params=netlist_values(design))


def sampling_instants(ckt):
    """Recover a clock phase once, then hold it — as a real receiver does.

    The two samples are half a bit apart, so one lands on each drive level, and
    the phase is the one that opens the nominal eye widest.  That is clock
    recovery in one line, and it matters for more than realism:

      * **Hold it fixed** so the objective stays smooth.  Re-optimising the phase
        every iteration would let the sampler chase the waveform, and a sampling
        index that jumps makes the gradient jump with it.
      * **Sample mid-bit**, not at the peaks.  On the peaks a slow receiver has
        already settled, so more gain always wins and `Rl` runs away to its
        bound — the bandwidth half of the trade-off simply disappears from the
        objective.  Where you sample decides which physics you are optimising.
    """
    run = simulate(ckt, (L_REF, 900.0))
    t, v = run.time, run.probes["v"]
    half = int(round(0.5 * TB / STEP))
    candidates = np.flatnonzero((t >= STOP - TB) & (t + 0.5 * TB <= STOP))
    phase = max(candidates, key=lambda i: abs(v[i] - v[i + half]))
    a, b = phase, phase + half
    return (a, b) if v[a] > v[b] else (b, a)


# ── the design problem, in JAX ────────────────────────────────────────────────
def build_objective(ckt, i_one, i_zero):
    import jax
    import jax.numpy as jnp
    from fairchild.jax import differentiable

    sim = differentiable(ckt, probes={"v": "pout"}, params=PARAMS,
                         step=STEP, stop=STOP, method="tr", reltol=1e-11)

    def design(u):
        """Unconstrained `u` → bounded (L, Rl).  Smooth, so grad flows through."""
        s = jax.nn.sigmoid(u)
        return (L_MIN + (L_MAX - L_MIN) * s[0], R_MIN + (R_MAX - R_MIN) * s[1])

    def loss(u):
        length, load = design(u)
        # The one place the physics of the design is written down.  jax.grad
        # composes it with the adjoint; nothing here differentiates the circuit.
        values = jnp.array([VPI_REF * L_REF / length, C_PER_UM * length, load])
        v = sim(values)["v"]
        return -(v[i_one] - v[i_zero])       # maximise eye height

    return loss, design, jax.jit(jax.value_and_grad(loss))


def optimise(value_and_grad, design, iters=90, rate=0.25, verbose=True):
    """Adam on the unconstrained variables."""
    u = START.copy()
    m = np.zeros_like(u)
    s = np.zeros_like(u)
    path, history = [], []

    for i in range(1, iters + 1):
        loss, g = value_and_grad(u)
        loss, g = float(loss), np.asarray(g, dtype=float)
        length, load = (float(x) for x in design(u))
        path.append((length, load))
        history.append(-loss)
        m = 0.9 * m + 0.1 * g
        s = 0.999 * s + 0.001 * g * g
        u -= rate * (m / (1 - 0.9 ** i)) / (np.sqrt(s / (1 - 0.999 ** i)) + 1e-12)
        if verbose and (i % 15 == 0 or i == 1):
            print(f"  iter {i:3d}   eye {-loss * 1e3:6.1f} mV   "
                  f"L {length:6.0f} um (V_pi {VPI_REF * L_REF / length:4.2f} V, "
                  f"Cmod {C_PER_UM * length * 1e15:4.0f} fF)   Rl {load:6.0f} Ohm")

    length, load = (float(x) for x in design(u))
    return (length, load), np.array(path), np.array(history)


def landscape(ckt, i_one, i_zero, n=13):
    lengths = np.linspace(L_MIN + 60, L_MAX - 60, n)
    loads = np.linspace(R_MIN + 60, R_MAX - 400, n)
    eye = np.empty((n, n))
    for a, length in enumerate(lengths):
        for b, load in enumerate(loads):
            v = simulate(ckt, (length, load)).probes["v"]
            eye[a, b] = v[i_one] - v[i_zero]
    return lengths, loads, eye


def plot(ckt, best, path, history, i_one, i_zero, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sig = 1.0 / (1.0 + np.exp(-START))
    start = (L_MIN + (L_MAX - L_MIN) * sig[0], R_MIN + (R_MAX - R_MIN) * sig[1])
    r0, r1 = simulate(ckt, start), simulate(ckt, best)
    t = r0.time * 1e12

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.3))

    ax = axes[0]
    ax.plot(t, r0.probes["v"] * 1e3, "--", color="tab:red", label="start")
    ax.plot(t, r1.probes["v"] * 1e3, color="tab:green", label="optimised")
    for i, style in ((i_one, "^"), (i_zero, "v")):
        ax.plot(t[i], r1.probes["v"][i] * 1e3, style, color="k", ms=9, zorder=5)
    ax.set(xlabel="time (ps)", ylabel="v(pout)  (mV)",
           title="Received waveform (markers = sampling instants)")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)

    bx = axes[1]
    lengths, loads, eye = landscape(ckt, i_one, i_zero)
    cs = bx.contourf(loads, lengths, eye * 1e3, levels=18, cmap="viridis")
    fig.colorbar(cs, ax=bx, label="eye height (mV)")
    bx.plot(path[:, 1], path[:, 0], "-o", color="w", ms=3, lw=1.4)
    bx.plot(best[1], best[0], "*", color="tab:red", ms=16)
    bx.set(xlabel="receiver load Rl (Ω)", ylabel="modulator length L (µm)",
           title="Co-design landscape and the optimiser's path")

    cx = axes[2]
    cx.plot(history * 1e3, color="tab:blue")
    cx.set(xlabel="iteration", ylabel="eye height (mV)",
           title="Adam on the adjoint gradient")
    cx.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out, dpi=120)
    print(f"wrote {out}")


# ── selftest ─────────────────────────────────────────────────────────────────
def selftest(ckt, best, history, i_one, i_zero) -> int:
    length, load = best

    # 1. The optimiser improved the link, and landed inside the bounds rather
    #    than against them — an optimum on a bound would mean the trade-off this
    #    example claims to show is not actually in the model.
    assert history[-1] > history[0], f"eye did not improve: {history[0]} → {history[-1]}"
    assert L_MIN + 20 < length < L_MAX - 20, f"L hit a bound at {length:.0f} um"
    assert R_MIN + 20 < load < R_MAX - 20, f"Rl hit a bound at {load:.0f} Ohm"

    # 2. It is the actual optimum: no neighbour on a coarse grid does better.
    best_eye = history[-1]
    for dl in (-0.25, 0.25):
        for dr in (-0.25, 0.25):
            v = simulate(ckt, (length * (1 + dl), load * (1 + dr))).probes["v"]
            neighbour = v[i_one] - v[i_zero]
            assert neighbour <= best_eye + 1e-6, (
                f"L={length * (1 + dl):.0f}, Rl={load * (1 + dr):.0f} beats the "
                f"optimum: {neighbour * 1e3:.2f} mV > {best_eye * 1e3:.2f} mV")

    # 3. The gradient crossing both domains is the real one — checked against a
    #    full re-solve of the transient, which is the only reference that proves
    #    the adjoint differentiates the system the integrator actually solved.
    #
    #    Deliberately *not* at the optimum: that is where the gradient is zero,
    #    and agreeing about zero to a few percent proves nothing.  Check it where
    #    the slope is steep enough to be worth getting right.
    at = (L_MIN + (L_MAX - L_MIN) * 0.5, R_MIN + (R_MAX - R_MIN) * 0.5)  # mid-range
    run = simulate(ckt, at)
    cot = np.zeros(len(run.time))
    cot[i_one], cot[i_zero] = 1.0, -1.0
    g = run.backward({"v": cot}, PARAMS)

    base = netlist_values(at)

    def eye_at(values):
        v = ckt.tran_adjoint({"v": "pout"}, step=STEP, stop=STOP, method="tr",
                             reltol=1e-11, params=values).probes["v"]
        return v[i_one] - v[i_zero]

    for k, (name, delta) in enumerate(zip(PARAMS, (1e-6, 1e-19, 1e-4))):
        hi, lo = dict(base), dict(base)
        hi[name] += delta
        lo[name] -= delta
        fd = (eye_at(hi) - eye_at(lo)) / (2 * delta)
        err = abs(g[k] - fd) / abs(fd)
        assert err < 1e-5, f"d(eye)/d{name}: adjoint {g[k]:.8e} vs re-solve {fd:.8e} ({err:.2e})"

    print(f"selftest OK — eye {history[0] * 1e3:.1f} → {history[-1] * 1e3:.1f} mV at "
          f"L={length:.0f} um (V_pi {VPI_REF * L_REF / length:.2f} V, "
          f"Cmod {C_PER_UM * length * 1e15:.0f} fF), Rl={load:.0f} Ohm; "
          f"gradient matches a full re-solve to <1e-5")
    return 0


def main() -> int:
    st = "--selftest" in sys.argv
    try:
        import jax
    except ImportError:
        sys.exit("this example needs JAX:\n"
                 "  uv run --with jax python examples/optimization/eo_link_codesign.py\n"
                 "For the same optimisation without a framework, see "
                 "examples/optimization/step_response_match.py")

    # Float64, before anything else touches JAX.  At float32 there are ~7 digits
    # to spend, which is less than the spread between a well-chosen
    # finite-difference step and a bad one — the gradient would be reporting the
    # dtype rather than the circuit.  `fairchild.jax` refuses to build without it.
    jax.config.update("jax_enable_x64", True)

    ckt = fc.Circuit()
    ckt.load_str(NETLIST)
    i_one, i_zero = sampling_instants(ckt)

    if not st:
        print("co-designing modulator length and receiver load:")
    _, design, value_and_grad = build_objective(ckt, i_one, i_zero)
    best, path, history = optimise(value_and_grad, design, verbose=not st)
    rc = selftest(ckt, best, history, i_one, i_zero)
    if not st:
        plot(ckt, best, path, history, i_one, i_zero, HERE / "eo_link_codesign.png")
    return rc


if __name__ == "__main__":
    sys.exit(main())
