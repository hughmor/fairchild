# fairchild User Guide

---

## Table of Contents

1. [SPICE Netlist Syntax](#1-spice-netlist-syntax)
2. [Elements Reference](#2-elements-reference)
3. [Waveform Sources](#3-waveform-sources)
4. [Model Cards](#4-model-cards)
5. [Analysis Directives](#5-analysis-directives)
6. [CLI Reference](#6-cli-reference)
7. [Output Formats](#7-output-formats)
8. [Solver Theory](#8-solver-theory)
9. [Convergence Knobs](#9-convergence-knobs)
10. [OSDI Model Loading](#10-osdi-model-loading)
11. [Photonic Model Library](#11-photonic-model-library)

---

## 1. SPICE Netlist Syntax

A fairchild netlist is a text file following standard SPICE conventions:

- **First line** is always the title (ignored by the parser, used in output headers).
- Lines beginning with `*` are comments.
- Keywords are case-insensitive (`NMOS` = `nmos`).
- Numbers accept SI suffixes: `k` = 1e3, `m` = 1e-3, `u` = 1e-6, `n` = 1e-9, `p` = 1e-12, `f` = 1e-15, `g` = 1e9, `t` = 1e12.
- Node `0` (also `gnd`, `GND`) is ground — always zero volts.

```spice
* My circuit title
V1  in 0  DC 5
R1  in out 1k
C1  out 0  1u
.tran 10u 5m
.end
```

---

## 2. Elements Reference

### Resistor
```
R<name>  <pos>  <neg>  <resistance>
R1  in out  1k
```

### Capacitor
```
C<name>  <pos>  <neg>  <capacitance>
C1  out 0  1u
```

### Inductor
```
L<name>  <pos>  <neg>  <inductance>
L1  a b  1m
```

### Voltage Source
```
V<name>  <pos>  <neg>  DC <value>
V<name>  <pos>  <neg>  PULSE(<v0> <v1> <td> <tr> <tf> <pw> <per>)
V1  vdd 0  DC 3.3
V2  in  0  PULSE(0 1 10n 1n 1n 40n 100n)
```

### Current Source
```
I<name>  <pos>  <neg>  DC <value>
I<name>  <pos>  <neg>  PULSE(<v0> <v1> <td> <tr> <tf> <pw> <per>)
I1  a 0  DC 1m
```

### Diode
```
D<name>  <anode>  <cathode>  <model_name>
D1  in  out  D1N4148
```

Requires a `.model` card with `D` type.

### MOSFET (Level 1)
```
M<name>  <drain>  <gate>  <source>  <bulk>  <model_name>  [W=<w>]  [L=<l>]
MN  out in 0 0  nm  W=10u L=1u
MP  out in vdd vdd  pm  W=20u L=1u
```

Requires a `.model` card with `NMOS` or `PMOS` type.

---

## 3. Waveform Sources

### DC
```
DC <value>
```
Constant value for all time. Used for bias supplies.

### PULSE
```
PULSE(<v0> <v1> <td> <tr> <tf> <pw> <per>)
```

| Parameter | Description |
|-----------|-------------|
| `v0` | Initial value (V or A) |
| `v1` | Pulsed value |
| `td` | Delay before first edge (s) |
| `tr` | Rise time v0 → v1 (s) |
| `tf` | Fall time v1 → v0 (s) |
| `pw` | Pulse width at v1 (s) |
| `per` | Period (s); 0 = single pulse |

The waveform is piecewise-linear. The DC operating point uses `v0`.

### PWL
```
PWL(<t0> <v0> <t1> <v1> ...)
```

Piecewise-linear waveform defined by (time, value) pairs. Time values must be in ascending
order. Before `t0` the waveform holds `v0`; after the last point it holds the last value.

```spice
V1  in 0  PWL(0 0  1u 5  2u 5  3u 0)
```

The DC operating point uses the value at `t = t0`.

---

## 4. Model Cards

```
.model  <name>  <type>  (<param>=<value>  ...)
```

### Diode (Shockley model)
```
.model D1N4148 D (IS=2.52n N=1.752)
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `IS` | Saturation current (A) | 1e-14 |
| `N` | Ideality factor | 1.0 |

### MOSFET Level 1 (Shichman-Hodges)
```
.model nm NMOS (VTO=0.7 KP=100u LAMBDA=0.01 GAMMA=0.5 PHI=0.6)
.model pm PMOS (VTO=-0.7 KP=50u)
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `VTO` | Threshold voltage (V) | 0.5 (NMOS), −0.5 (PMOS) |
| `KP` | Transconductance parameter (A/V²) | 20µ |
| `LAMBDA` | Channel-length modulation (1/V) | 0.0 |
| `GAMMA` | Body-effect coefficient (V^0.5) | 0.0 |
| `PHI` | Surface potential (V) | 0.6 |

Instance parameters `W` (width) and `L` (length) override model defaults. The effective
transconductance is `β = KP × W/L`.

---

## 5. Analysis Directives

### DC Operating Point
```
.op
```
Solves the nonlinear DC equations via Newton-Raphson. Output: one row of node voltages
and branch currents.

### Transient
```
.tran  <step>  <stop>
```
Integrate from t=0 to t=`stop` with timestep `step`. Uses Backward Euler for the first
step then Trapezoidal Rule for all subsequent steps (BE+TR, order 2).

Example:
```
.tran 10n 1u     * 10 ns step, 1 µs total
```

### AC Sweep
```
.ac  DEC|OCT|LIN  <points>  <fstart>  <fstop>
```

Small-signal AC sweep around the DC operating point.

| Variation | Meaning |
|-----------|---------|
| `DEC` | Logarithmic; `points` per decade |
| `OCT` | Logarithmic; `points` per octave |
| `LIN` | Linear; `points` total across `[fstart, fstop]` |

Example:
```
.ac dec 20 1 100k    * 20 points/decade from 1 Hz to 100 kHz
```

AC sweep can also be specified entirely via CLI flags (see Section 6), which override any
`.ac` directive in the netlist.

---

## 6. CLI Reference

```
fairchild [OPTIONS] --file <FILE>

Options:
  -f, --file <FILE>                   Input SPICE netlist
      --format <FORMAT>               csv (default) or nutmeg
  -o, --output <FILE>                 Output file (default: stdout)
      --probe <SIGNAL,...>            Comma-separated signals to include in CSV output
                                      e.g. --probe "V(out),I(V1)"
      --param <ELEMENT.PARAM=VALUE>   Override a circuit parameter (repeatable)
                                      e.g. --param "Vbias.dc=-2.0" --param "Rload.value=2e3"
      --check                         Parse + discipline-check only; exit 0 if valid
      --list-nodes                    Print all node names then exit
      --list-models                   Print .model cards + .osdi paths then exit
  -v, --verbose                       Show iteration counts and timing
  -q, --quiet                         Suppress all warnings
  -h, --help                          Print help
```

**--param** supported element types:

| Element type | Accepted param names |
|---|---|
| `VoltageSource` (Vxxx) | `dc`, `value`, `v` |
| `CurrentSource` (Ixxx) | `dc`, `value`, `i` |
| `Resistor` (Rxxx) | `resistance`, `value`, `r` |
| `Capacitor` (Cxxx) | `capacitance`, `value`, `c` |
| `Inductor` (Lxxx) | `inductance`, `value`, `l` |
| `XOsdi` (Xxxx) | any model parameter name |

### Examples

```bash
# DC only
fairchild -f netlist.sp

# Transient to file
fairchild -f netlist.sp -o waveforms.csv

# Nutmeg format (read back with ngspice or Python)
fairchild -f netlist.sp --format nutmeg -o waveforms.raw

# Filter output to selected signals
fairchild -f netlist.sp --probe "V(ph_a),V(T_dev)"

# Override bias voltage and coupling, verbose output
fairchild -f mrr_modulator_dc.sp --param "Vbias.dc=-2.0" --param "Xmod.kappa_0=0.05" -v

# Validate netlist (no simulation)
fairchild -f netlist.sp --check

# Voltage sweep in a shell loop
for V in 0 -1 -2 -3 -4 -5; do
  fairchild -f mrr_modulator_dc.sp --param "Vbias.dc=$V" --probe "V(ph_a)" -q | tail -1
done
```

---

## 7. Output Formats

### CSV

Default. One row per timepoint (transient), one row per frequency (AC), or one data row (DC).

```
time,V(out),V(in),I(v1)
0.000000e0,0.000000e0,0.000000e0,0.000000e0
1.000000e-8,3.297000e0,...
```

### Nutmeg rawfile

ASCII Nutmeg format compatible with ngspice's `rawread` command and the `spyci` Python
library. Contains a header block (title, variable names, types) followed by tabular data.

```bash
# Read back in ngspice:
ngspice
> rawread waveforms.raw
> plot V(out)
```

---

## 8. Solver Theory

### Modified Nodal Analysis (MNA)

fairchild assembles the circuit equations in MNA form: `A·x = b`, where:
- `x` is the unknown vector: `[V(n₁), V(n₂), ..., I(V₁), I(V₂), ...]`
- Node voltages occupy the first `n_nodes` entries.
- Voltage-source branch currents occupy the remaining entries (KVL constraints).
- Ground node is eliminated (always 0 V).

Each element stamps its conductance/current contributions into `A` and `b` following
standard MNA stamp rules.

### DC Newton-Raphson

For nonlinear circuits (diodes, MOSFETs), the DC solve uses Newton-Raphson iteration:

```
x_{k+1} = x_k − J(x_k)⁻¹ · f(x_k)
```

where `f(x) = A·x − b(x)` is the residual and `J = ∂f/∂x` is the Jacobian.

Each nonlinear device provides analytic `g_m, g_ds, g_mbs` (conductances) and a Norton
equivalent current (`J_eq`) that linearize around the current operating point.

**Convergence criteria** (SPICE standard):
- `|ΔV| < VNTOL + RELTOL · |V|` for all node voltages
- `|ΔI| < ABSTOL + RELTOL · |I|` for all branch currents

**Damping**: Per-iteration voltage steps are clamped to ±`VMAX = 0.5 V` to prevent
divergence in strongly nonlinear regions.

**Homotopy** (convergence aids):
- *GMIN stepping*: add a small conductance `GMIN = 1e-12 S` across every nonlinear element,
  then ramp it to zero. Makes the initial Jacobian well-conditioned.
- *Source stepping*: ramp all voltage/current sources from 0 to their final values in steps.
  Provides a continuous path from the linear (zero-source) solution.

### Transient Integration

#### Backward Euler (BE, order 1)

For a capacitor with voltage V and capacitance C, the companion model at timestep h is:

```
G_eq = C/h
I_hist = −(C/h) · V(t − h)
```

This stamps as an equivalent conductance + current source, integrating the MNA system at each
time step. BE is unconditionally stable and first-order accurate: LTE ∝ h².

#### Trapezoidal Rule (TR, order 2)

The TR companion model uses the average of the current and previous derivatives:

```
G_eq = 2C/h
I_hist = (2C/h) · V(t − h) + I_L(t − h)
```

TR is second-order accurate (LTE ∝ h³) but can exhibit numerical ringing on stiff circuits
near sharp transitions. fairchild uses BE for the first step then switches to TR.

#### Variable-step BE + LTE control

The variable-step solver (`tran_nr_var`) estimates the local truncation error at each step
using a predictor-corrector approach:

1. **Predict**: `x_pred = x + (h/h_prev) · (x − x_prev)` (linear extrapolation)
2. **Correct**: run Newton-Raphson from `x_pred` to get `x_corr`
3. **LTE estimate**: `lte = max_i |x_corr[i] − x_pred[i]| · 0.5 / (VNTOL + RELTOL · |x_corr[i]|)`
4. **Accept/reject**:
   - Accept if `lte ≤ 1` or `h ≤ h_min`
   - Reject and retry with `h_new = h · (0.9 / lte)^0.5`
5. **Step control**: `h_new = h · (0.9 / lte)^0.5`, clamped to `[0.1h, 4h]` and `h ≤ step`

### AC Small-Signal Analysis

After finding the DC operating point x₀, the solver:

1. Evaluates each nonlinear device's linearized conductances at x₀.
2. Builds the frequency-dependent admittance matrix `Y(jω) = G + jωC − jL/ω`.
3. For each frequency ω, solves `Y·V = I_ac` for complex node voltages.

The 2N×2N real block system avoids complex arithmetic:

```
[ G  −ωC ] [V_re]   [I_re]
[ ωC  G  ] [V_im] = [I_im]
```

---

## 9. Convergence Knobs

| Constant | Value | Description |
|----------|-------|-------------|
| `VNTOL` | 1 µV | Absolute voltage tolerance |
| `RELTOL` | 0.1% | Relative tolerance |
| `GMIN` | 1 pS | Minimum conductance (diagonal regularization) |
| `VMAX` | 0.5 V | Max ΔV per Newton iteration |
| `MAX_ITER` | 150 | Maximum Newton iterations per timepoint |

If NR fails to converge within `MAX_ITER` iterations, the solver tries source stepping and
GMIN stepping homotopy before returning an error.

---

## 10. OSDI Model Loading

fairchild can load compact models compiled by
[OpenVAF-Reloaded](https://codeberg.org/arpadbuermen/OpenVAF-Reloaded) as OSDI v0.4
shared libraries (`.osdi` files).

### Using OSDI models

```spice
.osdi /path/to/bsim4.osdi
.model nmos_bsim4 nmos4 (tox=3n vth0=0.4 ...)
```

The `.osdi` directive loads the shared library; the `.model` card name must match the
module name exported by the compiled model.

### Compiling Verilog-A models

```bash
# Install OpenVAF-Reloaded (see their README)
openvaf bsim4.va -o bsim4.osdi
```

### OSDI transient simulations

OSDI device models (e.g., BSIM4, PSP compiled via OpenVAF-Reloaded) are fully supported
in both DC and transient analyses. Reactive (capacitive) Jacobian entries are stamped
correctly via the `write_jacobian_array_react` copy path with `alpha = 1/h` scaling.
The variable-step solver respects OSDI reactive contributions at realistic step sizes.

---

## 11. Photonic Model Library

Fairchild includes a full library of Verilog-A photonic device models spanning
five device families at L1/L2/L3 abstraction levels:

- **CW laser** (`cw_laser`) — drives 3-wire optical ports (re, im, lambda)
- **PN junction phase shifters** (`pn_phase_shifter_l1/l2`) — carrier depletion modulation
- **Thermo-optic phase shifters** (`thermo_phase_shifter_l1/l2`) — resistive heating with external thermal RC
- **Photodetectors** (`photodetector`, `photodetector_l2`) — responsivity + bias-dependent Cj
- **MRR modulators** (`mrr_modulator_l1/l2/l3`) — all-pass CMT ring + PN, Soref-Bennett, TPA
- **N-doped heater MRR** (`mrr_heater_l1/l2`) — thermo-optic ring tuning

For detailed usage, port conventions, parameter tables, and common pitfalls, see:

**[docs/photonic_models.md](photonic_models.md)**

### Quick start

```bash
cd va-models && ./build.sh          # compile all models to build/*.osdi
cd examples/photonic
fairchild -f mrr_modulator_dc.sp --probe "V(ph_a)" --verbose
```

See `examples/photonic/mrr_voltage_sweep.py` for a parametric voltage sweep.
