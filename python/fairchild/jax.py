"""JAX adapter — make a fairchild transient run something `jax.grad` can see through.

The simulator is not traceable: it is a Newton solve in Rust, and JAX cannot
differentiate it by tracing.  It does not need to.  `jax.custom_vjp` lets you
hand JAX a forward rule and a vector-Jacobian product and be believed, and the
adjoint *is* a vector-Jacobian product — `backward()` takes `dL/d(probe)` at
every timepoint and returns `dL/dp`, which is exactly the signature a VJP rule
has to satisfy.

So the whole adapter is a shape declaration and two callbacks::

    import jax, jax.numpy as jnp
    from fairchild.jax import differentiable

    sim = differentiable(ckt, probes={"v": "pout"},
                         params=["Xmzm.V_pi", "Rd.r"],
                         step=10e-12, stop=2e-9)

    def loss(p):
        return jnp.sum((sim(p)["v"] - target) ** 2)

    g = jax.grad(loss)(p0)          # one adjoint solve, not one per parameter

Everything JAX gives you on top — `jit`, `vmap` over a batch of designs,
composing the circuit loss with a neural network's — then works, because from
JAX's side this is an ordinary differentiable primitive.

**Float64 is not optional here.**  JAX defaults to float32, which leaves ~7
digits — less than the spread between a good finite-difference step and a bad
one, and the gradient would be noise dressed as a number.  This module refuses
to run without ``jax.config.update("jax_enable_x64", True)``.
"""

from __future__ import annotations

from typing import Callable, Dict, Mapping, Sequence

import numpy as np


def differentiable(
    circuit,
    probes: Mapping[str, object],
    params: Sequence[str],
    **run_kwargs,
) -> Callable:
    """Wrap `circuit` as a differentiable `f(params) -> {probe: waveform}`.

    `probes` and `run_kwargs` are passed straight to
    :meth:`fairchild.Circuit.tran_adjoint`; `params` are `"element.param"`
    strings, and the argument to the returned function is an array of their
    values in that order.
    """
    import jax

    if not jax.config.read("jax_enable_x64"):
        raise RuntimeError(
            "fairchild.jax needs 64-bit JAX: call "
            'jax.config.update("jax_enable_x64", True) before building the '
            "simulator.  At float32 a circuit gradient is dominated by the "
            "dtype, not by the circuit."
        )

    names = list(params)
    probes = dict(probes)

    def _run(values: np.ndarray):
        overrides = {n: float(v) for n, v in zip(names, np.asarray(values, dtype=float))}
        return circuit.tran_adjoint(probes, params=overrides, **run_kwargs)

    # One run at the netlist's own values, to settle the output shapes JAX needs
    # declared up front — and to fail loudly here, rather than inside a traced
    # callback, if a probe name or a solver option is wrong.  Not at the caller's
    # first parameter vector, because there is no such thing yet, and not at
    # zeros, which is a short circuit and a zero-farad capacitor.
    probe_len = len(circuit.tran_adjoint(probes, **run_kwargs).time)
    shapes = {
        name: jax.ShapeDtypeStruct((probe_len,), np.float64) for name in probes
    }

    # ponytail: a one-entry cache, so the backward pass reuses the forward run
    # instead of solving the circuit twice.  Keyed on the parameter bytes, so a
    # miss is a re-run and never a wrong answer.  Grow it only if someone
    # genuinely interleaves grads at different points.
    cache: Dict[bytes, object] = {}

    def _cached(values: np.ndarray):
        key = np.asarray(values, dtype=float).tobytes()
        run = cache.get(key)
        if run is None:
            run = _run(values)
            cache.clear()
            cache[key] = run
        return run

    def _forward_np(values):
        run = _cached(values)
        out = run.probes
        return {name: np.asarray(out[name]) for name in probes}

    def _backward_np(values, cotangents):
        run = _cached(values)
        given = {
            name: np.ascontiguousarray(np.asarray(c), dtype=float)
            for name, c in cotangents.items()
        }
        return np.asarray(run.backward(given, names))

    @jax.custom_vjp
    def simulate(values):
        return jax.pure_callback(_forward_np, shapes, values)

    def _fwd(values):
        return simulate(values), values

    def _bwd(values, cotangents):
        grads = jax.pure_callback(
            _backward_np,
            jax.ShapeDtypeStruct((len(names),), np.float64),
            values,
            cotangents,
        )
        return (grads,)

    simulate.defvjp(_fwd, _bwd)
    return simulate
