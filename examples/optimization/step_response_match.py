#!/usr/bin/env python3
"""
Recover a filter's component values from its step response, by gradient descent.

The plainest possible use of the transient adjoint, on a circuit where you can
check every number by hand.  A two-pole RC ladder is simulated at known
component values to make a *target* waveform; the values are then thrown away,
the optimiser is started somewhere wrong, and it has to find its way back using
nothing but `dL/dp` from the adjoint.

      in ──[ R1 ]──┬──[ R2 ]──┬── out
                   │          │
                  C1         C2
                   │          │
                  gnd        gnd

Recovering known values is the point: any optimiser can make a loss go down, but
only a *correct* gradient walks to the answer that generated the data.  A sign
error, or a dropped history term in the co-state recursion, still descends —
just to somewhere else.

Why the adjoint rather than differencing the simulation:

  * **Cost.** One co-state solve covers every parameter, where forward
    differences pay a full transient re-solve each.  Two parameters here; a real
    filter has thirty.
  * **Accuracy, which matters more.** Differencing a converged transient
    differences a quantity known only to `reltol`, and that error compounds step
    over step.  The adjoint differences the *residual*, which is an explicit
    function evaluated to machine precision.  Below, the gradient is checked
    against a re-solve and agrees to better than 1e-5 — and that reference is
    itself the noisy side.

This example uses only numpy — no autodiff framework.  `run.backward()` returns
`dL/dp` for the netlist parameters, and the chain rule onto whatever design
variables you actually optimise is yours to apply; here it is the one line that
turns a log-space step into a component value, which keeps R and C positive
without any bounds machinery.  See `eo_link_codesign.py` for the same thing
done automatically by `jax.grad`.

Run:      .venv/bin/python examples/optimization/step_response_match.py
Selftest: same, with --selftest (asserts recovery + gradient, no plotting)
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

NETLIST = """* two-pole RC ladder, step excited
V1 in 0 PULSE(0 1 0 10p 10p 1 2)
R1 in mid 1k
C1 mid 0 100p
R2 mid out 1k
C2 out 0 220p
.tran 1n 1u
"""

# The values the target waveform is made from, and which the optimiser must find.
#
# Two, not three.  With the topology fixed the ladder's response has exactly two
# independent coefficients — `R1C1 + R1C2 + R2C2` and `R1R2C1C2` — so fitting
# three components is degenerate: a whole one-parameter family produces the
# identical waveform, and the optimiser drifts along it forever with the loss
# still falling.  Identifiability is a property of the problem, not of the
# gradient, and no amount of correct `dL/dp` fixes an under-determined fit.
TRUTH = {"R1.r": 2.2e3, "C1.c": 47e-12}
START = {"R1.r": 1.0e3, "C1.c": 150e-12}

STEP, STOP = 5e-9, 1.5e-6


def simulate(ckt: "fc.Circuit", values: dict):
    """One differentiable transient at `values`."""
    return ckt.tran_adjoint({"v": "out"}, step=STEP, stop=STOP,
                            method="tr", reltol=1e-11, params=values)


def loss_and_grad(ckt, names, u, target):
    """Sum-of-squares against `target`, and its gradient w.r.t. log-parameters.

    The optimiser works in `u = log(p)` so components stay positive and a step
    means the same thing to a kilohm and to a picofarad.  That reparameterisation
    is a chain rule the caller owns: `dL/du = dL/dp · p`.
    """
    values = {n: float(np.exp(x)) for n, x in zip(names, u)}
    run = simulate(ckt, values)
    y = run.probes["v"]
    residual = y - target
    # dL/dy at every timepoint, which is all the adjoint needs to know about
    # the loss — it never sees the loss itself.
    grad_p = run.backward({"v": 2.0 * residual}, names)
    return float(np.sum(residual ** 2)), grad_p * np.exp(u), y


def optimise(ckt, target, iters=500, rate=0.15, verbose=True):
    """Adam on the log-parameters.  Ten lines, so the gradient stays visible."""
    names = list(TRUTH)
    u = np.log([START[n] for n in names])
    m = np.zeros_like(u)
    v = np.zeros_like(u)
    history = []

    for i in range(1, iters + 1):
        loss, g, _ = loss_and_grad(ckt, names, u, target)
        history.append(loss)
        m = 0.9 * m + 0.1 * g
        v = 0.999 * v + 0.001 * g * g
        u -= rate * (m / (1 - 0.9 ** i)) / (np.sqrt(v / (1 - 0.999 ** i)) + 1e-12)
        if verbose and (i % 100 == 0 or i == 1):
            got = ", ".join(f"{n}={np.exp(x):.4g}" for n, x in zip(names, u))
            print(f"  iter {i:4d}   loss {loss:.3e}   {got}")

    return names, np.exp(u), np.array(history)


def build():
    ckt = fc.Circuit()
    ckt.load_str(NETLIST)
    target = simulate(ckt, TRUTH).probes["v"]
    return ckt, target


def plot(ckt, target, names, found, history, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = simulate(ckt, TRUTH).time * 1e9
    before = simulate(ckt, START).probes["v"]
    after = simulate(ckt, dict(zip(names, found))).probes["v"]

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax.plot(t, target, "k", lw=2.4, label="target")
    ax.plot(t, before, "--", color="tab:red", label="start")
    ax.plot(t, after, color="tab:green", label="recovered")
    ax.set(xlabel="time (ns)", ylabel="v(out)  (V)",
           title="Step response: start → recovered")
    ax.legend(loc="lower right")
    ax.grid(alpha=0.3)

    bx.semilogy(history, color="tab:blue")
    bx.set(xlabel="iteration", ylabel="Σ (v − target)²",
           title="Adam on the adjoint gradient")
    bx.grid(alpha=0.3, which="both")

    fig.tight_layout()
    fig.savefig(out, dpi=120)
    print(f"wrote {out}")


# ── selftest ─────────────────────────────────────────────────────────────────
def selftest(ckt, target, names, found) -> int:
    # 1. The optimiser found the values that made the data.
    for name, got in zip(names, found):
        want = TRUTH[name]
        err = abs(got - want) / want
        assert err < 2e-2, f"{name}: recovered {got:.5g}, truth {want:.5g} ({err:.1%} off)"

    # 2. The gradient is the real one.  Check it against a full re-solve of the
    #    transient either side of nominal — the only reference that proves the
    #    adjoint differentiates the system the integrator actually solved.
    at = {"R1.r": 1.5e3, "C1.c": 80e-12}
    run = simulate(ckt, at)
    y = run.probes["v"]
    keys = list(at)
    g = run.backward({"v": 2.0 * (y - target)}, keys)

    def loss_at(values):
        v = simulate(ckt, values).probes["v"]
        return float(np.sum((v - target) ** 2))

    for k, (key, delta) in enumerate([("R1.r", 1e-3), ("C1.c", 1e-17)]):
        hi, lo = dict(at), dict(at)
        hi[key] += delta
        lo[key] -= delta
        fd = (loss_at(hi) - loss_at(lo)) / (2 * delta)
        err = abs(g[keys.index(key)] - fd) / abs(fd)
        assert err < 1e-5, f"d(loss)/d{key}: adjoint {g[k]:.8e} vs re-solve {fd:.8e} ({err:.2e})"

    # 3. A parameter that reaches nothing is an error, not a zero.
    try:
        run.backward({"v": np.zeros_like(y)}, ["Rnope.r"])
    except RuntimeError as e:
        assert "reach nothing" in str(e), f"wrong error: {e}"
    else:  # pragma: no cover - guard
        raise AssertionError("a missing element should have been reported")

    print("selftest OK — recovered " +
          ", ".join(f"{n}={v:.4g} (truth {TRUTH[n]:.4g})" for n, v in zip(names, found)) +
          "; gradient matches a full re-solve to <1e-5")
    return 0


def main() -> int:
    st = "--selftest" in sys.argv
    ckt, target = build()
    if not st:
        print("recovering R1 and C1 from a step response:")
    names, found, history = optimise(ckt, target, verbose=not st)
    rc = selftest(ckt, target, names, found)
    if not st:
        plot(ckt, target, names, found, history, HERE / "step_response_match.png")
    return rc


if __name__ == "__main__":
    sys.exit(main())
