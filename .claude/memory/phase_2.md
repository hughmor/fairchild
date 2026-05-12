# Phase 2 — Photonic Discipline

**Goal**: First working optical circuit: CW laser → waveguide → directional coupler → photodetector, co-simulated with electronics in the same Newton-Raphson loop.

**Milestone**: Transient simulation of ring resonator modulation; resonance wavelength within 0.1 nm of coupled-mode theory analytical solution.

**Status**: 🔄 In progress — 3 of 4 ring resonator tests pass; sweep test assertions need one more fix

---

## Steps

### Step 1: OpenVAF extension validation (critical path) ✅

Verify OpenVAF-Reloaded handles a custom `optical` discipline without modification.

**Hypothesis**: OpenVAF treats all ports as generic analog nodes; discipline info is used only for connection checking (our elaborator), not the compiled OSDI code. If true, no compiler modification needed.

**Fallback**: fork OpenVAF-Reloaded and add optical nature/discipline to the front-end.

Action: compile a minimal `.va` with `optical` discipline and verify the OSDI output loads.

### Step 2: Optical discipline definition ✅

File: `va-models/disciplines/optical.vams`

```verilog-ams
nature Optical_Amplitude
  units = "sqrt(W)";  access = Oamp;  abstol = 1e-12;
endnature
nature Optical_Phase
  units = "rad";  access = Ophase;  abstol = 1e-9;
endnature
discipline optical
  potential Optical_Phase;
  flow Optical_Amplitude;
enddiscipline
```

### Step 3: Port keikawa model library ✅

Source: `keikawa/Verilog-A-photonic-model-library`

Validate/adapt these for our optical discipline:
- `waveguide.va` — propagation loss, group index
- `directional_coupler.va` — wavelength-dependent coupling ratio
- `cw_laser.va` — CW output with power, RIN noise, linewidth
- `photodetector.va` — responsivity, bandwidth, shot noise

### Step 4: Discipline checking in elaborator ✅

Add to `fairchild-parser`: discipline tracking per net; emit error on optical↔electrical mismatch. Mixed-domain components (modulator, photodetector) declare both.

### Step 5: Verilog-A model golden tests ✅

ngspice cannot simulate `.va` natively (ADMS is fragile and not our target).
**Reference simulator for `.va` models: VACASK** (or Xyce with ADMS-compiled models).

Test strategy for OSDI-compiled VA models:
- Compile `.va` with OpenVAF-Reloaded → load via OSDI → run DC sweep / transient
- Compare voltage/current output against VACASK running the same `.va` netlist
- For photonic models where no VACASK reference exists: validate against coupled-mode theory
  or waveguide analytical expressions (document this explicitly — these are not "golden" tests
  in the regression sense, but physics checks)

Each VA model in `va-models/` should have a paired test in `crates/fairchild-osdi/tests/`
following the pattern in `osdi_dc_op.rs`.

### Step 6: `@cross` event detection (deferred from Phase 1.5)

Zero-crossing step refinement. Implement after first optical circuit is working.

---

## Implementation notes (learned during Phase 2)

**Norton-equivalent only**: VA models MUST use `Oamp(p) <+` flow contributions, never `Ophase(p) <+` potential contributions. Potential contributions add internal branch-variable nodes in OpenVAF's OSDI output (`num_nodes > num_terminals`), which the runtime doesn't allocate. Use Norton equivalent: `Oamp(out) <+ G*(target - Ophase(out))` with `G = 1e6`. This also means input ports that only read their voltage (high-impedance) should have a tiny shunt: `Oamp(in) <+ 1e-12 * (-Ophase(in))` to ensure the node appears in the Jacobian.

**XOsdi in SPICE**: `X<name> net0 ... netN model_name [key=val ...]`. Nets map to OSDI terminal indices in declaration order. Discipline checking via `.optical net1 net2 ...` + `check_disciplines()`. Build the XOsdi model library path with `.osdi <path>` before the X instance.

---

## OSDI runtime bugs fixed (session 2026-05-12)

### Bug 1: `set_real_param` used wrong access() id
OSDI `access()` uses the **absolute** `param_opvar` index, not the relative
model-param index. The `param_opvar` array is laid out as:
`[inst_params (0..n_inst) | model_params (n_inst..n_total) | opvars]`
OpenVAF always inserts at least `$mfactor` as inst_params[0], so for cw_laser
`n_inst=1`: power_mW is at absolute index 1, not 0. Using the wrong index
mapped the write to a different offset in the model buffer (no visible effect on sim).

Fix: iterate `params[0..n_inst]` with `id = PARA_KIND_INST | i`, and
`params[n_inst..n_total]` with `id = PARA_KIND_MODEL | i` (absolute `i`).

### Bug 2: model writes not propagated to instance
OSDI `eval()` reads model-param-derived values from the **instance struct**
(pre-computed by `setup_instance`), not the model struct directly. After
`access(SET)` updates the model buffer, a `refresh_instance()` call (re-runs
`setup_instance`) is required to push the new value into the instance cache.

Fix: `refresh_instance()` added to `OsdiDevice`; called from `set_real_param`
after every model-param write.

### Bug 3: `build_devices` silently dropped XOsdi element params
`Element::XOsdi { nets, model_name, .. }` destructured `params` as `..`,
so netlist params like `kappa_0=0.1` were never applied (kappa_0 stayed at
default 0.5, giving only 1.3% resonance dip instead of 8.4%).

Fix: destructure as `Element::XOsdi { nets, model_name, params, .. }` and
call `dev.set_real_param(name, *value)` for each `(name, value)` in `params`.

---

## Sweep test status and remaining fix

Three tests pass: `access_ptr_diagnostic`, `set_real_param_verification`,
`ring_single_point_diagnostic`.

The sweep test `ring_resonator_wavelength_sweep` fails due to two wrong assertions:

**Issue 1 — wrong CMT reference resonance**:
`cmt_resonance_nearest(1551e-9)` always returns the mode-271 resonance at
1549.82 nm. But the simulation's global minimum falls at 1544.14 nm, which is
the mode-272 resonance (CMT: 1544.12 nm, Δ = 0.02 nm — within tolerance).
The fix is to use `cmt_resonance_nearest(sim_res_nm * 1e-9)` so the CMT lookup
is relative to wherever the simulation found the dip.

**Issue 2 — wrong off-resonance voltage**:
`v_off_resonance = p_in * (1 - kappa_0) * R_load = 0.9 V` is wrong. At ring
anti-resonance the through-port transmission approaches 1.0 (not 1-kappa).
Actual V_max from sweep ≈ 1.0 V.  The `v_min < v_off_resonance * 0.8`
check (threshold = 0.72 V) can never pass when V_min ≈ 0.93 V.
Fix: compute V_off as `sweep_results.iter().map(|(_,v)| *v).fold(NEG_INFINITY, f64::max)`,
and tighten threshold to `v_min < v_off_resonance * 0.975` (≥2.5% dip).
With nearest-sample offset of 0.02 nm and FWHM ≈ 0.098 nm, expected
apparent T_min ≈ 0.930 (7% dip), so 2.5% is a safe conservative threshold.

---

## Next steps for next session

1. Fix `ring_resonator_wavelength_sweep` assertions as described above.
2. Verify sweep test passes end-to-end.
3. Add `examples/ring_resonator_sweep.sp` SPICE netlist for the sweep.
4. Add `scripts/plot_ring_sweep.py` to plot V(ph_a) vs wavelength alongside
   CMT transmission curve. The sweep test already writes CSV to a temp path —
   the example should write it to `examples/ring_resonator_sweep.csv`.
5. Commit and update phase status to ✅.

---

## Key design: carrier frequency (ω₀)

`ω₀` is a **connection attribute**, not a MNA node:
- Set by laser model parameter `wavelength_nm`
- Propagated at elaboration time through passive connections
- Passed to device `eval()` via `SimContext.omega_0` (constant across NR iterations)
- For OSDI models: inserted into `OsdiSimParas` before each `eval` call

WDM: each channel = separate `(o_re, o_im)` port pair with its own `ω₀_k`.

## Signal representation

Per optical port pair: two MNA nodes `(o_re, o_im)` = complex SVEA envelope `A(t)` in √W.
Photonic models stamp amplitude/phase contributions into those rows.
Carrier `ω₀` is passed via `SimContext`; models compute `β(ω₀)`, coupling ratios, etc.
