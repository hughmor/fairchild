# Fairchild Photonic Model Library

## Overview

Fairchild ships a library of Verilog-A photonic device models that compile to
OSDI shared libraries via OpenVAF.  Each family has three abstraction levels:

| Level | Name | Physics | Use case |
|-------|------|---------|----------|
| L1 | Quasi-static heuristic | Linear Δφ(V) or ΔT, no carrier transport | Fast sweeps, system-level design |
| L2 | Physics-based | Plasma dispersion, Shockley I–V, depletion Cj | Accurate modulator eye diagrams |
| L3 | Nonlinear effects | TPA self-heating, photocarrier generation | High-power, pulsed operation |

All models use the **3-wire optical discipline** (`optical` net type with
`Ophase`/`Oamp` potentials) and Norton-form current contributions for NR
stability.

---

## Building the Models

```bash
cd va-models
./build.sh               # compile all models
MODEL_FILTER="mrr pn" ./build.sh  # compile only matching models
```

Output: `va-models/build/*.osdi`

---

## Port Convention

Every optical port uses three wires: `re`, `im`, `lambda`.

```
Xlaser  lre lim wl  cw_laser  power_mW=1.0 wavelength_nm=1550.0
         ^^  ^^  ^^
         │   │   └── lambda (wavelength in µm, shared across circuit)
         │   └────── imaginary amplitude (√W SVEA envelope)
         └────────── real amplitude
```

The `wl` (lambda) node must be connected to **every** optical element's lambda
port.  A single `wl` node can serve the whole circuit.  The CW laser drives
this node to the design wavelength.

```spice
.optical  lre lim wl ore oim   ; declare all optical nodes
```

---

## CW Laser (`cw_laser`)

```spice
.osdi ../../va-models/build/cw_laser.osdi
Xlaser  out_re out_im out_lambda  cw_laser  power_mW=1.0 wavelength_nm=1550.0
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `power_mW` | 1.0 | Output power (mW) |
| `wavelength_nm` | 1550.0 | Carrier wavelength (nm) |

---

## PN Junction Phase Shifter

### L1 — Linear heuristic

```spice
.osdi ../../va-models/build/pn_phase_shifter_l1.osdi
Xpnps  in_re in_im wl  out_re out_im wl  anode cathode  pn_phase_shifter_l1 \
       L_um=500.0 n_g=4.2 alpha_dB_cm=3.0 Vpi_L=2.0 V_ref=0.0 wavelength_nm=1550.0
```

Port order: 3 optical in + 3 optical out + 2 electrical (anode, cathode).

| Parameter | Default | Description |
|-----------|---------|-------------|
| `L_um` | 500.0 | Waveguide length (µm) |
| `n_g` | 4.2 | Group index |
| `alpha_dB_cm` | 3.0 | Propagation loss (dB/cm) |
| `Vpi_L` | 2.0 | Voltage–length product for π phase shift (V·cm) |
| `V_ref` | 0.0 | Reference voltage for Δφ=0 (V) |
| `wavelength_nm` | 1550.0 | Design wavelength (nm) |

**Physics:** Δφ = π × (V − V_ref) × L_cm / Vpi_L

### L2 — Soref-Bennett + Shockley + depletion Cj

```spice
.osdi ../../va-models/build/pn_phase_shifter_l2.osdi
Xpnps  in_re in_im wl  out_re out_im wl  anode cathode  pn_phase_shifter_l2 \
       L_um=500.0 n_g=4.2 alpha_dB_cm=3.0 Vbi=0.9 n_dep0=1e17 \
       IS=1e-14 n_id=1.5 Cj0=50e-15 mj=0.5 wavelength_nm=1550.0
```

**Additional parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Vbi` | 0.9 | Built-in PN junction voltage (V) |
| `n_dep0` | 1e17 | Peak depletion carrier density at V=0 (cm⁻³) |
| `IS` | 1e-14 | Saturation current for Shockley I–V (A) |
| `n_id` | 1.5 | Ideality factor |
| `Cj0` | 50e-15 | Zero-bias junction capacitance (F) |
| `mj` | 0.5 | Grading coefficient (0.5=abrupt) |

---

## Thermo-optic Phase Shifter

### L1 — Integrated thermal resistance

```spice
.osdi ../../va-models/build/thermo_phase_shifter_l1.osdi
Xheat  in_re in_im wl  out_re out_im wl  heat_p heat_n  thermo_phase_shifter_l1 \
       L_um=200.0 n_g=4.2 alpha_dB_cm=2.5 \
       R_heater=1000.0 R_thermal=50000.0 dn_dT=1.86e-4 wavelength_nm=1550.0
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `R_heater` | 1000.0 | Heater electrical resistance (Ω) |
| `R_thermal` | 50000.0 | Thermal resistance (K/W) |
| `dn_dT` | 1.86e-4 | Thermo-optic coefficient (K⁻¹); Si at 1550 nm |
| `L_um` | 200.0 | Waveguide length (µm) |

**Physics:** P = V²/R\_heater; ΔT = P × R\_thermal; Δφ = 2π·n\_g·dn\_dT·ΔT·L/λ

### L2 — External T_node electro-thermal port

Allows user to connect any thermal RC network to `T_node` (V = ΔT in K).

```spice
.osdi ../../va-models/build/thermo_phase_shifter_l2.osdi
Xheat  in_re in_im wl  out_re out_im wl  heat_p heat_n  T_node  thermo_phase_shifter_l2 \
       L_um=200.0 n_g=4.2 alpha_dB_cm=2.5 \
       R_heater=1000.0 dn_dT=1.86e-4 wavelength_nm=1550.0

* External thermal circuit: τ = Rth × Cth = 50 ms
Rth  T_node  0  50000   ; 50 kΩ → 50 kK/W
Cth  T_node  0  1e-6    ; 1 µF  → 1 µJ/K
```

Port order: 3 optical in + 3 optical out + 2 electrical (heat_p, heat_n) + 1 electrical (T_node).

---

## Photodetector

### L1 — Basic PIN (`photodetector`)

```spice
.osdi ../../va-models/build/photodetector.osdi
Xpd  in_re in_im wl  anode cathode  photodetector  responsivity=1.0
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `responsivity` | 1.0 | A/W at operating wavelength |
| `I_dark_A` | 1e-9 | Dark current (A) |
| `R_shunt` | 1e6 | Junction shunt resistance (Ω) |

### L2 — Bias-dependent junction capacitance (`photodetector_l2`)

```spice
.osdi ../../va-models/build/photodetector_l2.osdi
Xpd  in_re in_im wl  anode cathode  photodetector_l2 \
     responsivity=1.0 Cj0=50e-15 Vbi=0.9 mj=0.5
```

Adds `Cj(V) = Cj0 / (1 - V/Vbi)^mj` with transient charge contribution
`I_cap = d/dt(Cj·V)`. Critical for accurate bandwidth simulations.

---

## MRR Modulator

All-pass ring resonator with embedded PN junction.  Transfer function:
```
T = (r − a·e^{jφ}) / (1 − r·a·e^{jφ})
where r = √(1−κ), a = exp(−α·L/2), φ = 2π·n_g·L/λ + Δφ_PN
```

### L1 — Quasi-static CMT + linear PN

```spice
.osdi ../../va-models/build/mrr_modulator_l1.osdi
Xmod  in_re in_im wl  out_re out_im wl  anode cathode  mrr_modulator_l1 \
      kappa_0=0.1 L_ring_um=100.0 n_g=4.2 alpha_dB_cm=2.0 \
      Vpi_rt=10.0 V_ref=0.0 wavelength_nm=1550.0
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `kappa_0` | 0.1 | Power coupling coefficient (0–1) |
| `L_ring_um` | 100.0 | Ring circumference (µm) |
| `n_g` | 4.2 | Group index |
| `alpha_dB_cm` | 2.0 | Round-trip propagation loss (dB/cm) |
| `Vpi_rt` | 10.0 | Voltage for π round-trip phase shift (V) |
| `V_ref` | 0.0 | PN reference voltage for Δφ=0 (V) |
| `wavelength_nm` | 1550.0 | Design wavelength / bootstrap guard (nm) |

**FSR:** λ²/(n\_g·L\_ring) ≈ 5.72 nm (for L=100 µm at 1550 nm, n\_g=4.2)

### L2 — Soref-Bennett + Shockley + FCA

Adds Soref-Bennett plasma dispersion (Δn, Δα from ΔN), Shockley diode I–V,
and free-carrier absorption in the round-trip loss.

```spice
.osdi ../../va-models/build/mrr_modulator_l2.osdi
Xmod  in_re in_im wl  out_re out_im wl  anode cathode  mrr_modulator_l2 \
      kappa_0=0.1 L_ring_um=100.0 n_g=4.2 alpha_dB_cm=2.0 \
      Vbi=0.9 n_dep0=1e17 IS=1e-14 n_id=1.5 Cj0=50e-15 wavelength_nm=1550.0
```

### L3 — TPA self-heating + photocarrier generation

Adds two-photon absorption (TPA: β_TPA ≈ 0.7 cm/GW for Si), intra-ring
power buildup, free-carrier absorption from TPA-generated carriers, and
a thermal T_node port for electro-thermal co-simulation.

```spice
.osdi ../../va-models/build/mrr_modulator_l3.osdi
Xmod  in_re in_im wl  out_re out_im wl  anode cathode  T_node  mrr_modulator_l3 \
      kappa_0=0.1 L_ring_um=100.0 n_g=4.2 alpha_dB_cm=2.0 \
      Vpi_rt=10.0 beta_tpa=0.7 A_eff_um2=0.1 tau_fc_ns=10.0 wavelength_nm=1550.0
Rth  T_node  0  5000   ; external thermal impedance
Cth  T_node  0  1e-9   ; external thermal capacitance
```

---

## N-doped Heater MRR

Ring resonator with N-doped waveguide heater (higher propagation loss from FCA
but no separate optical phase shifter required).

### L1 — Integrated thermal model

```spice
.osdi ../../va-models/build/mrr_heater_l1.osdi
Xring  in_re in_im wl  out_re out_im wl  heat_p heat_n  mrr_heater_l1 \
       kappa_0=0.1 L_ring_um=100.0 n_g=4.2 alpha_dB_cm=10.0 \
       R_heater=500.0 R_thermal=30000.0 dn_dT=1.86e-4 wavelength_nm=1550.0
```

Typical: R\_heater=500 Ω, alpha=10 dB/cm (N-doped FCA), dn/dT=1.86e-4 K⁻¹.

### L2 — External T_node

Same as L1 but with external thermal port for full transient thermal simulation.

```spice
.osdi ../../va-models/build/mrr_heater_l2.osdi
Xring  in_re in_im wl  out_re out_im wl  heat_p heat_n  T_node  mrr_heater_l2 \
       kappa_0=0.1 L_ring_um=100.0 n_g=4.2 alpha_dB_cm=10.0 \
       R_heater=500.0 dn_dT=1.86e-4 wavelength_nm=1550.0
Rth  T_node  0  5000
Cth  T_node  0  1e-9
```

---

## CLI Usage

### Basic simulation

```bash
fairchild -f netlist.sp                          # DC op, print all signals
fairchild -f netlist.sp --probe "V(ph_a)"        # filter to one signal
fairchild -f netlist.sp --verbose                # show iteration counts
fairchild -f netlist.sp --format nutmeg -o out.raw
```

### Parameter overrides

```bash
# Override element parameters without editing the netlist
fairchild -f mrr_modulator_dc.sp --param "Vbias.dc=-2.0" --probe "V(ph_a)"
fairchild -f ring_sweep.sp --param "Xlaser.wavelength_nm=1549.0"
fairchild -f circuit.sp --param "Rload.resistance=2e3" --param "Xcoupler.kappa_0=0.05"
```

Supported element types and param names:
- `VoltageSource`: `dc`, `value`, `v`
- `CurrentSource`: `dc`, `value`, `i`
- `Resistor`: `resistance`, `value`, `r`
- `Capacitor`: `capacitance`, `value`, `c`
- `Inductor`: `inductance`, `value`, `l`
- `XOsdi`: any model parameter name

### Bias voltage sweep (shell loop)

```bash
for V in 0 -1 -2 -3 -4 -5; do
  echo -n "V=$V: "
  fairchild -f mrr_modulator_dc.sp --param "Vbias.dc=$V" --probe "V(ph_a)" -q \
    | tail -1
done
```

### Validation

```bash
fairchild -f netlist.sp --check         # discipline-check only, exit 0/1
fairchild -f netlist.sp --list-nodes    # enumerate circuit nodes
fairchild -f netlist.sp --list-models   # show .model cards + .osdi paths
```

---

## Python API

After `maturin develop` (or `pip install`):

```python
import fairchild as fc
import numpy as np

ckt = fc.Circuit()
ckt.load("examples/photonic/mrr_modulator_dc.sp")

# Single DC operating point
result = ckt.run("op")
print(f"V(ph_a) = {result['V(ph_a)'][0]:.4f} V")

# Parametric voltage sweep
v_values = np.linspace(0, -5, 21).tolist()
results = ckt.sweep("vbias.dc", v_values, "op")
p_through = [r["V(ph_a)"][0] / 1e3 for r in results]   # mW

# Transient simulation
ckt.load("examples/photonic/thermo_heater_tran.sp")
result = ckt.run("tran", stop=0.2, step=1e-3)
time   = result.time()           # numpy array, seconds
v_pha  = result["V(ph_a)"]       # numpy array, volts
t_dev  = result["V(t_dev)"]      # thermal node temperature (K above ambient)
```

### Parameter overrides in Python

```python
ckt = fc.Circuit()
ckt.load("mrr_modulator_dc.sp")
ckt.set_param("Xlaser", "power_mW", 2.0)     # 2 mW laser
ckt.set_param("Xmod", "kappa_0", 0.05)       # tighter coupling
result = ckt.run("op")
```

---

## Common Pitfalls

**Lambda wire missing from laser:**  The CW laser has 3 ports — `out_re`,
`out_im`, `out_lambda`.  Always connect all three.  Omitting the lambda node
causes the lambda Norton source to drive the wrong node, corrupting the
amplitude signal:

```spice
* WRONG — missing lambda port, lambda drives out_re accidentally:
Xlaser  lre lim  cw_laser  power_mW=1.0 wavelength_nm=1550.0

* CORRECT:
Xlaser  lre lim wl  cw_laser  power_mW=1.0 wavelength_nm=1550.0
```

**Electrical ports as key=val:**  XOsdi instances specify nets positionally.
Electrical ports (anode, cathode, heat_p, heat_n, T_node) must appear in the
positional net list, not as `anode=vbias` key=val parameters:

```spice
* WRONG — electrical ports listed as params, they get dropped:
Xmod  lre lim wl  ore oim wl  mrr_modulator_l1  anode=vbias cathode=0

* CORRECT — all nets positional, params are key=val only:
Xmod  lre lim wl  ore oim wl  vbias 0  mrr_modulator_l1 \
      kappa_0=0.1 ...
```

**Missing .optical declaration:**  All nodes using the `optical` discipline
must be listed in a `.optical` statement for discipline checking to pass.

**T_node port needs Rth:**  L2 thermo and L3 MRR models have an external
`T_node` port.  Connect at least a resistor to ground (`Rth T_node 0 R`) or
the node will be floating.
