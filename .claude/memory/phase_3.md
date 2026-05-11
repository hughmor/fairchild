# Phase 3 — Python Bindings

**Goal**: Jupyter-native workflow — import fairchild as a Python package, build/run circuits programmatically, get numpy arrays back.

**Milestone**: Jupyter notebook demonstrating SiPh transceiver simulation with eye diagram.

**Status**: 📋 Not started (after Phase 2)

---

## Core API

```python
import fairchild as fc

ckt = fc.Circuit()
ckt.load("transceiver.sp")
result = ckt.run(analysis="tran", stop=10e-9, step=0.01e-9)

t = result.time()           # numpy array
v_out = result["V(out)"]    # numpy array
```

## Numpy waveform input

```python
t = np.linspace(0, 10e-9, 1000)
v_in = np.where((t % 1e-9) < 0.5e-9, 1.8, 0.0)

ckt.set_source("VIN", fc.WaveformSource(t, v_in))
result = ckt.run(analysis="tran", stop=10e-9)
```

`WaveformSource(t, v)` converts to internal PWL at the Rust boundary — the Phase 1.5 PWL infrastructure handles the rest. No Python-only solver code path.

## Implementation

Crate: `fairchild-py` (PyO3 + pyo3-numpy)
- `Circuit`, `Netlist`, `SimResults` Python classes
- Numpy array bridge for all output data
- Async: `.run_async()` returning a Future
- JAX compatibility: expose as `jax.pure_callback` (enables Optax, scipy.optimize)
- Packaging: `maturin` + `pyproject.toml` at repo root
