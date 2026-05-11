# Phase 4 — Differentiable Simulation

**Goal**: `∂L/∂p` flows through the photonic circuit. Gradient-based inverse design via adjoint method — cost proportional to one forward simulation, regardless of number of parameters.

**Milestone**: Optimize ring resonator coupling gap to hit target center wavelength; convergence in <100 gradient steps from random init.

**Status**: 📋 Not started (after Phase 3)

---

## Discrete adjoint for GEAR integrator

Given forward solution `x(t₀...tₙ)`, the adjoint `λ(t)` satisfies:
```
-dλ/dt = (∂f/∂x)ᵀ λ + (∂L/∂x)ᵀ
```
Integrated backward in time using the same GEAR/BE/TR coefficients as the forward pass.

## OSDI adjoint extension

OSDI already provides `∂f/∂x` for Newton-Raphson. Also need `∂f/∂p` (parameter Jacobian):
- Option A: add a call to the OSDI interface
- Option B: finite differences (simpler, slower — acceptable for initial implementation)

## Python API

```python
with fc.GradientContext(params=["ring.kappa", "ring.L"]) as ctx:
    result = ckt.run(analysis="tran", stop=10e-9)
    loss = result.target_spectrum(measured_s21)
    loss.backward()

grads = ctx.grads  # {"ring.kappa": float, "ring.L": float}
```

## Validation

Test: perturb each parameter by ε, run two forward simulations, compute FD gradient. Compare against adjoint gradient. Must agree to within numerical precision.

## Device trait extension

Add `load_jacobian_param(&self, mat: &mut ParamJacMatrix, p_idx: usize)` to `Device` trait so every stamp records symbolic dependency on parameters for the adjoint pass.
