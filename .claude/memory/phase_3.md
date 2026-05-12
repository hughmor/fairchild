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

## Parametric sweep API

Run a parameter sweep with parallel independent simulations:

```python
results = ckt.sweep("R1.resistance", [1e3, 2e3, 5e3, 10e3], analysis="tran", stop=10e-9)
# returns list of SimResults, one per parameter value
```

Rust-side design:
- `run_sweep(netlist, param_name, values: Vec<f64>) -> Vec<TranResult>` in fairchild-core
- Parallelised with `rayon` (independent per-point; `Device: Send + Sync` already guaranteed)
- Sweepable: element values (R, L, C), `.model` params (IS, N, VTO, KP), source amplitudes
- CLI: `--sweep "R1.resistance=1k,2k,5k,10k"` → one CSV column-group per sweep point

This is also a prerequisite for Phase 4 inverse design: the sweep loop feeds the adjoint.

## Implementation

Crate: `fairchild-py` (PyO3 + pyo3-numpy)
- `Circuit`, `Netlist`, `SimResults` Python classes
- Numpy array bridge for all output data
- Async: `.run_async()` returning a Future
- JAX compatibility: expose as `jax.pure_callback` (enables Optax, scipy.optimize)
- Packaging: `maturin` + `pyproject.toml` at repo root
