# fairchild User Guide

---

## Contents

1. [Netlist syntax](#1-netlist-syntax)
2. [Elements reference](#2-elements-reference)
3. [Waveform sources](#3-waveform-sources)
4. [Model cards](#4-model-cards)
5. [Analyses](#5-analyses)
6. [Directives](#6-directives)
7. [SimOptions and convergence knobs](#7-simoptions-and-convergence-knobs)
8. [CLI reference](#8-cli-reference)
9. [Python bindings](#9-python-bindings)
10. [Output formats](#10-output-formats)
11. [Solver theory](#11-solver-theory)
12. [Photonic devices](#12-photonic-devices)
13. [OSDI model loading](#13-osdi-model-loading)

---

## 1. Netlist syntax

A fairchild netlist follows standard SPICE conventions:

- First line is the title (used in output headers).
- `*` and `;` start comments. `;` may follow content on a line.
- Continuation lines start with `+`.
- Keywords are case-insensitive (`NMOS` = `nmos`).
- SI suffixes: `k`=1e3, `meg`=1e6, `g`=1e9, `t`=1e12, `m`=1e-3, `u`=1e-6,
  `n`=1e-9, `p`=1e-12, `f`=1e-15.
- Node `0` (also `gnd`, `GND`) is ground.

```spice
* My circuit title
V1  in 0  DC 5
R1  in out 1k
C1  out 0  1u
.tran 10u 5m
.end
```

---

## 2. Elements reference

### Passive elements

```
R<name>  <pos> <neg>  <resistance>
C<name>  <pos> <neg>  <capacitance>     [IC=<v0>]
L<name>  <pos> <neg>  <inductance>      [IC=<i0>]
```

### Independent sources

```
V<name>  <pos> <neg>  <DC|PULSE|PWL|SIN|EXP|SFFM|AM>(...)
I<name>  <pos> <neg>  <DC|PULSE|PWL|SIN|EXP|SFFM|AM>(...)
```

See §3 for waveform shapes.

### Diode

```
D<name>  <anode> <cathode>  <model_name>
```

Requires a `.model … D (…)` card.

### MOSFET (Level 1, Shichman-Hodges)

```
M<name>  <drain> <gate> <source> <bulk>  <model_name>  [W=<w>]  [L=<l>]
```

Requires a `.model … NMOS|PMOS (…)` card.

### Behavioral source (B-element)

```
B<name>  <pos> <neg>  V=<expression>
B<name>  <pos> <neg>  I=<expression>
```

Expressions reference node voltages `V(n)`, branch currents `I(V…)`, time
`time`, and parameters `{name}`. Operators: `+ - * / ^`, unary minus,
parentheses. Functions: `abs sgn sqrt exp ln log10 sin cos tan asin acos atan
sinh cosh tanh min max if`. The expression is differentiated symbolically for
the Jacobian; supports static and dynamic cases.

```spice
B1  out 0  V=V(in)*V(in) + 0.5*sin(2*3.14159*1k*time)
```

### Subcircuit instantiation

```
X<name>  <node1> <node2> ...  <subckt_or_model>  [key=val ...]
```

Used for `.subckt` hierarchical blocks, OSDI device instances, and native
photonic devices.

---

## 3. Waveform sources

| Shape | Form |
|---|---|
| Constant | `DC <value>` |
| Pulse | `PULSE(<v0> <v1> <td> <tr> <tf> <pw> <per>)` |
| Piecewise linear | `PWL(<t0> <v0> <t1> <v1> ...)` |
| Sinusoid | `SIN(<v0> <va> <freq> [<td> <theta> [<phi>]])` |
| Exponential | `EXP(<v0> <v1> <td1> <tau1> <td2> <tau2>)` |
| Single-frequency FM | `SFFM(<v0> <va> <fc> <mod_index> <fs>)` |
| AM | `AM(<sa> <fc> <fm> <off> [<td>])` |

The DC operating point uses `v0` (or the value at `t = t0` for PWL).

---

## 4. Model cards

```
.model  <name>  <type>  (<param>=<value>  ...)
```

### Diode (`D`)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `IS` | Saturation current (A) | 1e-14 |
| `N`  | Ideality factor | 1.0 |

### MOSFET Level 1 (`NMOS` / `PMOS`)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `VTO` | Threshold voltage (V) | 0.5 / −0.5 |
| `KP` | Transconductance (A/V²) | 20µ |
| `LAMBDA` | Channel-length modulation (1/V) | 0 |
| `GAMMA` | Body-effect coefficient (V^0.5) | 0 |
| `PHI` | Surface potential (V) | 0.6 |

Instance `W` / `L` override model defaults.

---

## 5. Analyses

### DC operating point

```
.op
```

Solves the nonlinear DC equations via Newton-Raphson. Returns one row of node
voltages + branch currents.

### DC sweep

```
.dc  <src>  <start>  <stop>  <step>  [<src2> <start2> <stop2> <step2>]
```

Sweeps `src` (a voltage or current source) over the range, optionally nested
with a second source.

### Transient

```
.tran  <step>  <stop>  [UIC]
```

Integrates from 0 to `stop`. `UIC` (or `.options uic`) skips DC and uses
`.ic` / element `IC=` values as the initial state.

Integration method is selected by `.options method=be|tr|gear` (default TR
with BE first-step bootstrap; GEAR/BDF-2 recommended for stiff circuits or
where ringing on TR is unacceptable).

### AC sweep

```
.ac  DEC|OCT|LIN  <points>  <fstart>  <fstop>
```

Small-signal sweep around the DC operating point.

### Noise

```
.noise  V(<out_pos>[,<out_neg>])  <input_src>  DEC|OCT|LIN  <points>  <fstart>  <fstop>
```

Adjoint-method noise analysis. Device noise sources:

| Device | Source |
|---|---|
| Resistor | 4kT/R thermal |
| Diode | 2qI_d shot |
| MOSFET | 8kTg_m/3 channel |

Output is `onoise` (V²/Hz at the output port) and `inoise` (equivalent input
PSD referred to `input_src`). OSDI devices can plug in via the `Device` trait's
`noise_sources()` hook.

---

## 6. Directives

### Initial conditions

```
.ic V(<node>)=<value> ...
.nodeset V(<node>)=<value> ...
```

`.ic` seeds transient when `UIC` is set; `.nodeset` is a soft hint for the DC
solver (cleared after the first iteration converges).

### Measurements

```
.measure tran  <name>  FIND  V(<n>)  AT=<t>
.measure tran  <name>  MAX|MIN|AVG|RMS|PP   V(<n>)  [FROM=<t>] [TO=<t>]
.measure tran  <name>  INTEG|DERIV  V(<n>)  [FROM=<t>] [TO=<t>]
.measure tran  <name>  TRIG  V(<n>)=<v> [CROSS=<k>]  TARG  V(<n>)=<v> [CROSS=<k>]
```

Post-processed after a `.tran` run; emitted in the output and exposed in
Python as `result.measurements`.

### Libraries and includes

```
.lib  "<file>"  <section>
... (section body) ...
.endl  <section>

.include "<file>"
```

`.lib` reads the named section out of a `.lib`/`.endl`-delimited file.
`.include` substitutes a file inline (depth ≤ 16; relative to the referencing
file).

### Parametrisation

```
.param  <name>=<value>  [<name2>=<value2> ...]
```

Parameters substitute into element values, model-card values, and B-element
expressions as `{name}`.

### Corner and temperature sweeps

```
.temp <T1> [<T2> ...]
.alter <label>
   ... overrides ...
.endalter
```

`.temp` re-runs every analysis once per listed temperature (°C). `.alter`
blocks describe deltas from the base netlist; each block produces a full
re-run with overrides applied.

### Solver options

```
.options  <key>=<value>  ...
```

Any field of `SimOptions` (§7) can be set this way.

### Subcircuits

```
.subckt <name>  <p1> <p2> ...
...
.ends [name]
```

Two-pass parse-time flattening; supports nested `.subckt` and `.param`.

### OSDI

```
.osdi  <path/to/model.osdi>
.model <card_name> <module_name> (<params>)
```

Loads a compiled OpenVAF-Reloaded shared library; `.model` names the SPICE
card and matches the module exported in the library.

### Bus vectors and optical ports

```
.optical_port  <name>  [<N>]
```

Declares an N-channel optical bundle. Each `name` becomes three (or `3·N`)
underlying wires `name_re[_k]`, `name_im[_k]`, `name_wl[_k]` for SVEA (re, im)
amplitude and per-channel wavelength. Photonic device instances on the bundle
are auto-replicated per channel. See §12.

```
net[M..N]    expands to net_M, net_{M+1}, ..., net_N
```

Bus vector expansion can appear anywhere an X-element net or `.optical_port`
name does.

---

## 7. SimOptions and convergence knobs

Every solver entry point takes a `SimOptions` struct. The same fields are
addressable from the netlist (`.options key=val`), CLI (`--opt key=val` or
convenience flags), and Python (`Circuit.run("…", key=val)`).

| Field | Default | Description |
|---|---|---|
| `reltol` | 1e-3 | NR relative tolerance |
| `abstol` | 1e-12 | NR current absolute tolerance (A) |
| `vntol` | 1e-6 | NR voltage absolute tolerance (V) |
| `vmax` | 0.5 | Per-iteration |ΔV| clamp |
| `gmin` | 1e-12 | Diagonal regularising conductance (S) |
| `gminmax` | 1.0 | GMIN-stepping starting value |
| `itl1` | 150 | DC max NR iterations |
| `itl4` | 50 | Transient per-step max NR iterations |
| `max_rejections` | 30 | Var-step max step rejections |
| `method` | `tr` | `be` / `tr` / `gear` |
| `max_step` | ∞ | Transient max step (s) |
| `srcsteps` | 11 | Source-stepping homotopy resolution |
| `temp_k` | 300.15 | Operating temperature (K) |
| `uic` | false | Use `.ic` / element `IC=` instead of DC |
| `pnjlim` | true | Diode / MOSFET junction limiting in NR |
| `solver` | `auto` | `auto` / `dense` / `sparse` linear backend |

Setting any of these from the netlist:

```spice
.options reltol=1e-5  method=gear  maxstep=1n  solver=sparse
```

From the CLI:

```bash
fairchild -f circuit.sp --opt reltol=1e-5 --opt method=gear \
                       --maxstep 1n --solver sparse
```

From Python:

```python
circuit.run("tran", step=1e-9, stop=100e-9,
            reltol=1e-5, method="gear", solver="sparse")
```

Convergence aids that don't appear as fields but always run:

- **GMIN stepping**: starts at `gminmax`, ramps to `gmin`.
- **Source stepping**: ramps sources 0 → final in `srcsteps` steps.
- **Pseudo-transient continuation**: not yet implemented.

---

## 8. CLI reference

```
fairchild [OPTIONS] --file <FILE>
```

| Flag | Description |
|---|---|
| `-f, --file <FILE>` | Input SPICE netlist |
| `--format csv\|nutmeg` | Output format (default csv) |
| `-o, --output <FILE>` | Output destination (default stdout) |
| `--probe <SIG,…>` | Comma-separated CSV signal filter |
| `--param ELEM.PARAM=VAL` | Override a circuit parameter (repeatable) |
| `--opt KEY=VAL` | Override a SimOptions field (repeatable) |
| `--reltol`, `--gmin`, `--method`, `--maxstep`, `--solver` | Convenience flags |
| `--no-pnjlim` | Disable junction limiting |
| `--check` | Parse + discipline-check only |
| `--list-nodes` / `--list-models` | Inspect parsed netlist, then exit |
| `-v` / `-q` | Verbose / quiet |

Examples:

```bash
# DC operating point
fairchild -f netlist.sp

# Transient with GEAR + tightened tolerance + sparse LU, into Nutmeg
fairchild -f netlist.sp \
   --opt method=gear --opt reltol=1e-5 --solver sparse \
   --format nutmeg -o out.raw

# AC sweep specified entirely via netlist; CSV-filter on two probes
fairchild -f rlc.sp --probe "V(out),I(V1)"

# Photonic transient with parameter override
fairchild -f examples/photonic/native_mrr_modulator.sp \
   --param "Xpn.V_pi_L=1.5e-3" --probe "V(pd_anode)"

# Parse-only check
fairchild -f circuit.sp --check
```

---

## 9. Python bindings

```bash
pip install fairchild     # once published; until then, see crates/fairchild-py/README
```

```python
import fairchild

c = fairchild.Circuit()
c.load("rc_step.sp")                     # load from file
# or: c.load_str(netlist_string)

# Override scalar element parameters before running:
c.set_param("Rload", "resistance", 2e3)

# Run any analysis with SimOptions kwargs:
op    = c.run("op")
tran  = c.run("tran", step=1e-9, stop=1e-6,
              method="gear", reltol=1e-5, variable_step=True)
ac    = c.run("ac", points="dec", n=20, fstart=1, fstop=1e6)
noise = c.run("noise", points="dec", n=20, fstart=1, fstop=1e6,
              v_out_pos="out", input_src="V1")

# Parametric sweep — equivalent of Monte Carlo / corner runs.
results = c.sweep("Rload.resistance", [1e3, 2e3, 5e3], "tran",
                  step=1e-9, stop=1e-6)

# Access:
t  = tran.time()
vo = tran["V(out)"]
meas = tran.measurements        # dict of name → value from .measure
```

`Circuit.set_source(name, WaveformSource.pulse(0, 1, 0, 1e-9, 1e-9, 10e-9, 20e-9))`
replaces a source's waveform without re-parsing.

---

## 10. Output formats

### CSV

Default. One row per timepoint (transient), frequency (AC/noise), or sweep
point (DC sweep). Header row is comma-separated signal names. `--probe`
filters the column set.

### Nutmeg rawfile

ASCII-only format compatible with ngspice's `rawread` and the `spyci` Python
library. Always emits the full signal set.

```bash
fairchild -f circuit.sp --format nutmeg -o waveforms.raw
ngspice -c "rawread waveforms.raw; plot V(out)"
```

---

## 11. Solver theory

### Modified Nodal Analysis (MNA)

fairchild assembles `A·x = b`, where `x` is `[V(n₁), …, V(nₖ), I_branch(V₁), …]`
and ground (node 0) is eliminated. Each device stamps its conductance and
current contributions following standard MNA stamp rules.

### DC Newton-Raphson

Linearise around the current iterate: `J(xₖ) · Δx = −f(xₖ)`, where
`f(x) = A(x)·x − b(x)`. Diodes and MOSFETs provide analytic `g_m, g_ds, g_mb`
and Norton equivalents (`J_eq`).

Convergence criteria:

- `|ΔV| < vntol + reltol · |V|` for all node voltages.
- `|ΔI| < abstol + reltol · |I|` for all branch currents.

Convergence aids: per-iteration `|ΔV|` clamp (`vmax`), junction limiting
(`pnjlim` for diodes, `fetlim` for MOSFETs), GMIN stepping (ramp diagonal
conductance from `gminmax` to `gmin`), source stepping (ramp sources 0 →
final in `srcsteps` steps).

### Transient integration

Three integrators, selected via `method`:

- **BE** (BDF-1): unconditionally stable, first-order. LTE ∝ h².
- **TR** (trapezoidal): second-order, A-stable. LTE ∝ h³. Can ring on near-
  step inputs.
- **GEAR** (BDF-2): second-order, L-stable. No ringing. Variable-coefficient
  stamp for non-uniform steps; demotes to BE on the first step, after a
  rejection, or on extreme step-ratio transitions.

The variable-step controller uses divided-difference LTE estimation,
accept/reject with `lte ≤ 1`, and `h_new = h · (0.9/lte)^0.5`, clamped to
`[0.1h, 4h] ∩ [h_min, max_step]`.

### AC small-signal

After DC, the solver linearises each device at the operating point and forms
the complex admittance matrix `Y(jω) = G + jωC − j/(ωL)`. Solved per
frequency via the 2N×2N real block representation.

### Noise

Adjoint-method: solve the adjoint system at each frequency, propagate device
noise sources (resistor 4kT/R, diode 2qI_d, MOSFET 8kTg_m/3, OSDI hook) to
the output. Equivalent input noise referred to the named input source.

---

## 12. Photonic devices

### Conventions

Every "optical port" is a **3-wire bundle** `(re, im, λ)`:

- `re`, `im` — slowly-varying-envelope (SVEA) complex amplitude in √W. The
  optical power at the port is `|A|² = V(re)² + V(im)²`. The carrier
  frequency is implicit; only the envelope is solved.
- `λ` — propagation wavelength in metres. A device-local wire that allows
  wavelength-dependent physics (e.g. waveguide propagation phase) without
  forcing a global parameter.

`.optical_port NAME [N]` declares a bundle. For `N > 1`, each bundle-using
device instance is auto-replicated into `N` parallel single-channel
instances by the parser. Electrical nets and scalar nets shared across the
replicas (e.g. one `V_pn` driving every replica of a PN-loaded ring) make
single-modulator-multi-wavelength topologies trivial — see
`examples/photonic/native_wdm_mrr_modulator.sp`.

All values internally are SI. Convenience aliases like `L_um` and
`wavelength_nm` accept the named unit and convert. **Beware SPICE SI
prefixes**: writing `L_um=8u` parses as `L_um = 8e-6` (because `u` =
1e-6 in SPICE), then the device interprets that as `8e-6 µm = 8 pm` —
almost certainly not what you meant. Either drop the suffix (`L_um=8`,
read as 8 µm) or use the SI alias (`L_m=8e-6`, read as 8 µm).

Direction of energy flow is fixed by the device — there are no
back-propagating modes in this version, and reflections / standing waves
are out of scope.

### `fc_cw_laser` — constant-wave laser

```
X<name>  out  fc_cw_laser  [param=val …]
```

| Port | Role |
|---|---|
| `out` | bundle, optical output |

| Parameter | Default | Description |
|---|---|---|
| `power_mW` | 1.0 | Output optical power. Sets `re² + im² = power_mW × 1e−3`. |
| `power_W` | — | Alternative spelling; overrides `power_mW`. |
| `phase_deg` | 0 | Initial phase of the SVEA carrier. |
| `wavelength_nm` | 1550 | Output wavelength (drives the `λ` wire). |
| `re_amp` / `im_amp` | derived | Direct override of the SVEA components. |

`power_mW`, `power_W`, and `phase_deg` are not orthogonal: `power_*` set
the magnitude of `(re, im)` while preserving phase; `phase_deg` rotates the
current magnitude. Setting `re_amp` / `im_amp` directly bypasses both.

**Physics.** Three direct-potential equations fix `V(out_re) = √P · cos(φ₀)`,
`V(out_im) = √P · sin(φ₀)`, `V(out_λ) = λ`. No electrical input; no noise;
no spectral linewidth.

### `fc_waveguide` — lossy waveguide

```
X<name>  in  out  fc_waveguide  [param=val …]
```

| Port | Role |
|---|---|
| `in`  | bundle, optical input |
| `out` | bundle, optical output |

| Parameter | Default | Description |
|---|---|---|
| `L_um` | 100 | Length, µm. |
| `L_m` / `length` | — | Length, m (overrides `L_um`). |
| `n_g` | 4.2 | Group index. |
| `alpha_dB_cm` | 2.0 | Power loss (dB/cm). |
| `wavelength_nm` | 1550 | **Bootstrap** wavelength used when the input `λ` wire is unbound (initial NR iterate). Otherwise the wire value wins. |
| `wavelength_m` | — | Same, in metres. |

**Physics.** `A_out = A_in · exp(−α L / 2) · exp(−j β L)` with `β = 2π n_g / λ`
and `α` in nepers/m (the `alpha_dB_cm` value is converted internally). The
`λ` wire is read at evaluation time, so wavelength-dependent propagation
phase is captured automatically — this is what makes the ring resonator
example see a true resonance dip when you sweep the laser wavelength.

### `fc_dcoupler` — 2×2 directional coupler

```
X<name>  a  b  c  d  fc_dcoupler  [param=val …]
```

| Port | Role |
|---|---|
| `a` | bundle, input arm 1 |
| `b` | bundle, input arm 2 |
| `c` | bundle, through output (paired with `a`) |
| `d` | bundle, cross output (paired with `a`) |

| Parameter | Default | Description |
|---|---|---|
| `kappa_per_m` / `kappa` | 100 | Coupling rate (rad/m). |
| `L_um` / `L_m` / `length` | 5e-3 m | Interaction length. |
| `kappa_L` / `kappaL` | — | Override for `κ·L` directly (preserves length, scales `kappa_per_m`). |

`kappa_L=0` gives a perfect through, `kappa_L=π/2` gives full cross. The
mass-action numbers in `examples/photonic/native_mrr_modulator.sp` use
`kappa_L=0.336` (≈ 11% cross-coupled power) for a critically-coupled ring.

**Physics.** Lossless coupling matrix `[c; d] = [cos(κL), j sin(κL); j sin(κL), cos(κL)] · [a; b]`.
In SVEA `(re, im)` form this is six direct-potential equations. `c_λ = a_λ`,
and `d_λ = a_λ` too — the second was originally `d_λ = b_λ` but that creates
a feedback loop with no driving source in closed-loop topologies (e.g. ring
resonators) and causes the PN PS to read `λ ≈ 0`. Both outputs now route
the input-arm-`a` wavelength. For asymmetric WDM topologies where you
genuinely want different wavelengths on the two arms, that's a future
device — file an issue.

### `fc_splitter` — 1×2 Y-junction (3 dB)

```
X<name>  in  out_a  out_b  fc_splitter
```

| Port | Role |
|---|---|
| `in` | bundle, optical input |
| `out_a`, `out_b` | bundles, optical outputs |

No parameters. Lossless equal-power split: `out_a = out_b = in / √2`.
Wavelength duplicated to both outputs. Use for combining MZI arms in
reverse (with a second `fc_splitter` instance) or for any branching
topology where you want half-power outputs.

### `fc_photodetector` — PIN photodetector

```
X<name>  in  anode  cathode  fc_photodetector  [param=val …]
```

| Port | Role |
|---|---|
| `in` | bundle, optical input |
| `anode` / `cathode` | scalar electrical nodes |

| Parameter | Default | Description |
|---|---|---|
| `responsivity` | 1.0 | A/W. |
| `i_dark` / `i_dark_a` | 1e-9 | Dark current (A). |
| `r_shunt` | 1e6 | Shunt resistance (Ω) — junction non-idealness. |

**Physics.** Photocurrent `I_ph = responsivity × (V(in_re)² + V(in_im)²) + i_dark`
flows from cathode to anode internally (reverse-biased convention).
Externally, the anode sources current. A linear shunt `1/r_shunt` is
stamped between anode and cathode. The photocurrent is nonlinear in the
optical amplitudes, so a Norton equivalent linearised at the current
operating point is contributed to the NR loop; `∂I_ph/∂V_re = 2R·V_re`,
`∂I_ph/∂V_im = 2R·V_im`.

No bandwidth limit, no transit-time delay, no avalanche gain. A linear
capacitance and finite responsivity bandwidth are future parameters.

### `fc_pn_ps` — PN-junction phase shifter

```
X<name>  in  out  anode  cathode  fc_pn_ps  [param=val …]
```

| Port | Role |
|---|---|
| `in` / `out` | bundles, optical pass-through |
| `anode` / `cathode` | scalar electrical nodes (PN-junction terminals) |

| Parameter | Default | Description |
|---|---|---|
| `L_um` | 1000 (1 mm) | PN-section length, µm. |
| `L_m` / `length` | — | Same, in metres. |
| `n_g` | 4.2 | Group index (sets the wavelength-dependent propagation phase). |
| `wavelength_nm` | 1550 | Reference wavelength: propagation phase is zero at `λ = wavelength_nm`. Pin this to your laser's λ so the device is "on resonance" by default. |
| `dn_dv` | 1e-4 | Effective-index change per applied volt (small-signal). |
| `V_pi_L` | — | Convenience override: `Vπ·L` in V·m. Setting this overrides `dn_dv` so that `φ = π` when `V = Vπ`. |
| `g_pn` | 1e-3 | Linearised PN-junction conductance (S). Connects anode and cathode through `1/g_pn`. |
| `alpha_dB_cm` | 0 | Propagation loss along the PN section. For a closed-loop ring this loss sets the extinction ratio of the resonance dip — without it the ring is all-pass. |

**Physics.** Optical: `φ = φ_prop + φ_eo` where

- `φ_prop = 2π n_g L (1/λ − 1/λ_ref)` — wavelength-dependent propagation
  phase, zeroed at the reference wavelength.
- `φ_eo = 2π L (dn/dV) V_pn / λ` — the electro-optic shift.

The transfer is `A_out = A_in · exp(−α L / 2) · exp(−j φ)`. Wavelength
passes through unchanged.

Electrical: a single linear conductance `g_pn` between anode and cathode.
The model does not yet capture the diode I-V or the bias-dependent
depletion-region capacitance — for a forward-biased carrier-injection
device, replace `g_pn` with a more complete junction model in a future
revision.

`V_pi_L` and `dn_dv` are not orthogonal: setting `V_pi_L` recomputes
`dn_dv` from the current `wavelength_nm` and `L`. If you set `V_pi_L`
and later change `wavelength_nm`, **re-set `V_pi_L`** so the EO calibration
tracks the new λ_ref.

### `fc_mux` — N → 1 WDM multiplexer

```
X<name>  bus  ch_0  ch_1  ...  ch_{N-1}  fc_mux
```

| Port | Role |
|---|---|
| `bus` | bundle, N-channel optical output |
| `ch_k` (k = 0..N-1) | bundle, single-channel optical input |

No parameters; channel count `N` is inferred from instance arity (number of
positional nets minus 1).

**Physics.** Identity routing per channel: `V(bus_k.re) = V(ch_k.re)`,
`V(bus_k.im) = V(ch_k.im)`, `V(bus_k.λ) = V(ch_k.λ)` for k = 0..N-1.
There's no wavelength selectivity — this is a topology marker, not an
AWG. For wavelength-selective combining you'd write a different device
that conditionally couples based on input wavelength; not in scope yet.

The parser special-cases `fc_mux` so that (a) the bus and channel
bundles can have different channel widths without erroring, and (b)
the device is NOT replicated per channel — one instance handles all N
channels at once.

### `fc_demux` — 1 → N WDM demultiplexer

```
X<name>  bus  ch_0  ch_1  ...  ch_{N-1}  fc_demux
```

| Port | Role |
|---|---|
| `bus` | bundle, N-channel optical input |
| `ch_k` (k = 0..N-1) | bundle, single-channel optical output |

Symmetric counterpart to `fc_mux`: `V(ch_k.*) = V(bus_k.*)`. Same
no-wavelength-selectivity caveat.

Typical WDM topology:

```spice
.optical_port ch0
.optical_port ch1
.optical_port bus 2
.optical_port out_bus 2
.optical_port d0
.optical_port d1

Xl0 ch0 fc_cw_laser power_mW=1.0 wavelength_nm=1549.9
Xl1 ch1 fc_cw_laser power_mW=1.0 wavelength_nm=1550.1
Xmux bus ch0 ch1 fc_mux
Xwg  bus out_bus fc_waveguide L_um=500 alpha_dB_cm=2 wavelength_nm=1550
Xdemux out_bus d0 d1 fc_demux
Xpd0 d0 v_pd0 0 fc_photodetector responsivity=0.8
Xpd1 d1 v_pd1 0 fc_photodetector responsivity=0.8
```

The `bus` and `out_bus` bundles each carry 2 channels; `Xwg` replicates
into 2 parallel single-channel waveguides automatically because both
its input and output bundles are 2-channel.

### `fc_thermal_ps` — thermo-optic phase shifter

```
X<name>  in  out  heat_p  heat_n  fc_thermal_ps  [param=val …]
```

| Port | Role |
|---|---|
| `in` / `out` | bundles, optical pass-through |
| `heat_p` / `heat_n` | scalar electrical nodes (heater drive) |

| Parameter | Default | Description |
|---|---|---|
| `r_heater` / `r` | 1000 | Heater resistance (Ω). |
| `p_pi` / `p_pi_w` | 10e-3 | Heater power for π phase shift (W). |

**Physics.** Electrical: linear resistor `R_heater` between `heat_p` and
`heat_n`. Joule power `P = V² / R_heater` is converted instantaneously into
an optical phase shift `φ = π · P / P_pi`. No thermal RC — the conversion
is algebraic, no transient lag. For a circuit where thermal time constants
matter (most real heaters), build an external RC between the drive node
and an intermediate `T_dev` net and feed `T_dev` to the heater pin pair
with a B-element converting temperature to equivalent voltage.

The transfer is `A_out = A_in · exp(−j φ)` — lossless. Wavelength passes
through unchanged.

### Registering PDK-specific aliases

Native devices can be aliased into a PDK-specific taxonomy without
duplicating the device. `DeviceRegistry::register_alias("pdk_foo_waveguide",
"fc_waveguide", remap)` registers a new card name that delegates to
`fc_waveguide` with a parameter-remapping closure. This is the PDK
private-fork extension point — it lets foundry naming and parameter
conventions live downstream without leaking into the upstream master
branch. See `crates/fairchild-core/src/device_registry.rs:300` for an
example.

### Legacy OSDI Verilog-A models

Pre-Phase-B Verilog-A models (MRR, MZI, PN-PS, thermo PS, photodetector) are
still loadable via `.osdi`. They use a different discipline (Norton-equivalent
flow contributions) and the 12-pin underlying-wire syntax. The CLI prints a
one-shot hint pointing at the native devices when a photonic `.osdi` library
is loaded. New work should use native devices; legacy models survive for
back-compat and third-party clear-text Verilog-A. See `docs/photonic_models.md`.

---

## 13. OSDI model loading

fairchild loads compact models compiled by [OpenVAF-Reloaded][openvaf] as
OSDI v0.4 shared libraries (`.osdi`).

```spice
.osdi /path/to/bsim4.osdi
.model nmos_bsim4 nmos4 (tox=3n vth0=0.4 ...)
```

The `.osdi` directive loads the shared library; the `.model` second token
must match the module name exported by the compiled model. Parameter names
match case-insensitively (VA preserves case; fairchild lowercases).

Reactive Jacobian contributions are stamped via the `write_jacobian_array_react`
copy path with `α = 1/h` scaling. The aliasing path (`load_jacobian_resist`)
is broken in OpenVAF-Reloaded and not used.

Encrypted PDK Verilog-A (IEEE-1735 / Cadence NCPROTECT) is fundamentally
unsupported by OpenVAF — that's an upstream Cadence-key problem, not a
fairchild limitation.

```bash
# Build
openvaf bsim4.va -o bsim4.osdi

# Use
fairchild -f my_circuit.sp
```

[openvaf]: https://codeberg.org/arpadbuermen/OpenVAF-Reloaded
