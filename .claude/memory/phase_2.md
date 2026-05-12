# Phase 2 — Photonic Discipline

**Goal**: First working optical circuit: CW laser → waveguide → directional coupler → photodetector, co-simulated with electronics in the same Newton-Raphson loop.

**Milestone**: Transient simulation of ring resonator modulation; resonance wavelength within 0.1 nm of coupled-mode theory analytical solution.

**Status**: 🔜 Not started

---

## Steps

### Step 1: OpenVAF extension validation (critical path)

Verify OpenVAF-Reloaded handles a custom `optical` discipline without modification.

**Hypothesis**: OpenVAF treats all ports as generic analog nodes; discipline info is used only for connection checking (our elaborator), not the compiled OSDI code. If true, no compiler modification needed.

**Fallback**: fork OpenVAF-Reloaded and add optical nature/discipline to the front-end.

Action: compile a minimal `.va` with `optical` discipline and verify the OSDI output loads.

### Step 2: Optical discipline definition

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

### Step 3: Port keikawa model library

Source: `keikawa/Verilog-A-photonic-model-library`

Validate/adapt these for our optical discipline:
- `waveguide.va` — propagation loss, group index
- `directional_coupler.va` — wavelength-dependent coupling ratio
- `cw_laser.va` — CW output with power, RIN noise, linewidth
- `photodetector.va` — responsivity, bandwidth, shot noise

### Step 4: Discipline checking in elaborator

Add to `fairchild-parser`: discipline tracking per net; emit error on optical↔electrical mismatch. Mixed-domain components (modulator, photodetector) declare both.

### Step 5: Verilog-A model golden tests

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
