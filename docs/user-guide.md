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
13. [Writing custom devices](#13-writing-custom-devices)
14. [OSDI model loading](#14-osdi-model-loading)

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
| `lambda_center_m` | 1.55e-6 | Photonic band-centre default (laser λ, PN-PS reference, waveguide bootstrap). Set via `lambda_center_nm` for nm units. |

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

`.optical_port NAME [N]` declares a bundle. How the parser handles a
multi-channel bundle depends on the device on it:

- **Pure-optical devices** (`fc_cw_laser`, `fc_waveguide`, `fc_dcoupler`,
  `fc_splitter`): the parser replicates the X-element into `N` parallel
  single-channel instances. Each instance handles one wavelength; they're
  independent.
- **Bundle-aware devices with shared electrical state**
  (`fc_pn_ps`, `fc_thermal_ps`, `fc_photodetector`): the parser does NOT
  replicate. One device instance handles all `N` optical channels
  internally while keeping one shared electrical interface
  (anode/cathode or heat_p/heat_n). This is what makes "single physical
  modulator, multiple wavelengths" work correctly — the V_pn supply sees
  one PN junction, not N parallel copies; the photodetector sums
  photocurrents across channels into one anode current and presents one
  dark current and one shunt.
- **Bundle-bridges** (`fc_mux`, `fc_demux`): the parser does NOT replicate
  and intentionally allows mismatched channel counts on different pins
  (N-channel bus side + N single-channel pins on the other side).

Electrical nets and scalar nets (`vmod`, `0`, etc.) wired into a bundle-
aware device are shared across all channels; the device sees one shared
voltage / current per electrical pin. See
`examples/photonic/native_wdm_mrr_modulator.sp` for a full topology.

### Band centre

Photonic devices need a default wavelength for things like the initial
NR iterate (when no laser has yet driven the λ wire) and for any
device parameter that defaults to "the design wavelength". One global
option sets this band-wide default:

```spice
.options lambda_center_nm=1310    * O-band
.options lambda_center_nm=1550    * C-band (default)
.options lambda_center_m=1.31e-6  * same, in metres
```

Or via the CLI: `--opt lambda_center_nm=1310`. Or in Python:
`Circuit.run("op", lambda_center_nm=1310)`.

Devices with their own `wavelength_nm` parameter (the laser's output
wavelength, the PN-PS's reference wavelength) override the band centre
when set explicitly. The waveguide doesn't have a `wavelength_nm`
parameter at all — its λ comes entirely from the input wire and the
band-centre is only a bootstrap fallback.

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

The `wavelength_nm` parameter is accepted for backward compatibility but
no longer does anything — the waveguide reads λ directly from the input
bundle's λ wire, and the laser drives that wire to whatever wavelength it
was configured with. A hard-coded 1.55 µm is used only to seed the very
first NR iterate (where the wire is still at 0 V); after iteration 1 the
laser's value wins.

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

**WDM behaviour.** Bundle-aware. On an N-channel optical input, ONE
`fc_photodetector` instance handles all N channels: it sums the photocurrents
(`I_ph = responsivity · Σ_k (V(re_k)² + V(im_k)²) + i_dark`) into one anode
current, presents one shared dark current, and stamps one shunt resistance
— not N copies. Responsivity is the same on every channel in this first-
pass model; per-channel responsivity (for wavelength-selective detection)
is a future parameter.

**Physics.** Photocurrent flows from cathode to anode internally
(reverse-biased convention). Externally, the anode sources current. A
linear shunt `1/r_shunt` is stamped between anode and cathode. The
photocurrent is nonlinear in the optical amplitudes, so a Norton
equivalent linearised at the current operating point is contributed to
the NR loop; `∂I_ph/∂V(re_k) = 2R·V(re_k)`, `∂I_ph/∂V(im_k) = 2R·V(im_k)`
for every channel `k`.

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
| `g_pn` | 1e-3 | Linearised PN-junction conductance (S). Connects anode and cathode through `1/g_pn`. ONE conductance shared across all N optical channels — see WDM note below. |
| `alpha_dB_cm` | 0 | Propagation loss along the PN section. For a closed-loop ring this loss sets the extinction ratio of the resonance dip — without it the ring is all-pass. |

**WDM behaviour.** `fc_pn_ps` is bundle-aware. On an N-channel optical
bus, ONE `fc_pn_ps` instance handles all N optical paths and presents ONE
shared electrical interface — the V_pn supply sees exactly `g_pn`, not
`N · g_pn`. Per-channel wavelength is read independently from each
channel's λ wire, so a single Vπ-driven modulator naturally produces
wavelength-diverse phase shifts.

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

**WDM behaviour.** Like `fc_pn_ps`, this is bundle-aware. One device, one
shared heater, all N optical channels see the same phase shift (no
wavelength dependence in this first-pass model).

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

## 13. Writing custom devices

This is the path for adding a new SPICE primitive to fairchild — say, a
better PN-junction model with proper diode I-V, or a new photonic device
class. You need a working Rust toolchain (`rustup install stable`) but
you don't need deep Rust experience: every existing device is a copy-
paste template, and the trait you implement has small methods that
follow the same MNA stamp patterns SPICE has used for 50 years.

The work has three parts:

1. Write a `struct` for your device + an `impl Device` block.
2. Register it in `device_registry.rs` so SPICE element names dispatch
   to your factory.
3. Add a regression test under `crates/fairchild-core/tests/`.

### 13.1 Where files live

| File | Purpose |
|---|---|
| `crates/fairchild-core/src/models/<your_family>.rs` | The device struct + `impl Device`. Put diodes in `diode.rs`, MOSFETs in `mosfet1.rs`, photonics in `photonic.rs`, or create a new file for a new family. |
| `crates/fairchild-core/src/models/mod.rs` | One-line `pub use` re-export so the type is visible from `crate::models`. |
| `crates/fairchild-core/src/device_registry.rs` | The factory call that maps `"my_device"` (SPICE card name) to a constructor. |
| `crates/fairchild-core/tests/<your_test>.rs` | Regression test — exercise DC OP / transient / AC against a closed-form expectation. |

### 13.2 Minimum Rust syntax you'll touch

```rust
use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

pub struct MyDevice {
    // Parameters as plain f64 fields. Defaults set in new().
    my_param: f64,
    // Terminal node indices. Use [NodeId; N] for fixed arity, Vec<NodeId> for variable.
    // NodeId is Option<usize>: None means "tied to ground", Some(i) is row i in the MNA matrix.
    nodes: [NodeId; 2],
    // Cached values from eval() that load_jacobian / load_residual will read.
    // Anything that depends on the current iterate goes here.
    g_cached: f64,
}

impl MyDevice {
    pub fn new() -> Self {
        Self { my_param: 1.0, nodes: [None; 2], g_cached: 0.0 }
    }
}

impl Device for MyDevice {
    fn num_terminals(&self) -> usize { 2 }
    fn setup_model(&mut self, _ctx: &SimContext) {}
    fn setup_instance(&mut self, terminals: &[NodeId], _ctx: &SimContext) {
        for i in 0..2 { self.nodes[i] = terminals[i]; }
    }
    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "my_param" => { self.my_param = value; true }
            _ => false,
        }
    }
    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        // Compute anything that depends on the current solution.  x[i] is
        // the voltage / current at MNA row i.  Store results in self.*_cached.
    }
    fn load_residual(&self, b: &mut [f64]) {
        // Add current contributions to the RHS vector b (one entry per MNA row).
    }
    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        // Add conductance contributions to mat.a[row][col].
    }
    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) { self.load_residual(b); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) { self.load_jacobian(mat); }
}
```

That's the minimum. The hard part is what to put inside `eval`,
`load_residual`, `load_jacobian` — that's the device physics. The
sections below give you the standard patterns.

### 13.3 Stamp patterns by physics type

#### Linear resistance (Ohm's law)

A resistor `R` between nodes `a` and `b` stamps a conductance `g = 1/R`:

```rust
fn load_jacobian(&self, mat: &mut MnaMatrix) {
    let g = 1.0 / self.r;
    if let Some(a) = self.nodes[0] {
        mat.a[a][a] += g;
        if let Some(b) = self.nodes[1] { mat.a[a][b] -= g; }
    }
    if let Some(b) = self.nodes[1] {
        mat.a[b][b] += g;
        if let Some(a) = self.nodes[0] { mat.a[b][a] -= g; }
    }
}
```

`load_residual` does nothing — a linear resistor has no current source
to contribute to the RHS, only diagonal/off-diagonal entries in the
Jacobian. The `if let Some(a) = self.nodes[i]` guards skip stamps to
ground (NodeId is None = ground).

#### Reactive element (capacitor / inductor) between nodes

A capacitor `C dv/dt = i` is integrated by the time-stepping rule. The
companion model at step h is `i = (C/h) · (v − v_prev)`. In MNA terms
that's a conductance `g_eq = C·α` (α = 1/h for Backward Euler, 2/h for
Trapezoidal — passed by the solver into `load_*_tran`) plus a history
current `−g_eq · v_prev` on the RHS.

```rust
pub struct MyCap {
    c: f64,
    nodes: [NodeId; 2],
    v_prev: f64,    // committed voltage from the previous accepted step
}

impl Device for MyCap {
    fn num_terminals(&self) -> usize { 2 }
    fn setup_instance(&mut self, t: &[NodeId], _: &SimContext) {
        self.nodes = [t[0], t[1]];
    }
    fn eval(&mut self, _x: &[f64], _: EvalFlags, _: &SimContext) {}

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, alpha: f64) {
        let g_eq = self.c * alpha;
        if let Some(a) = self.nodes[0] {
            mat.a[a][a] += g_eq;
            if let Some(b) = self.nodes[1] { mat.a[a][b] -= g_eq; }
        }
        if let Some(b) = self.nodes[1] {
            mat.a[b][b] += g_eq;
            if let Some(a) = self.nodes[0] { mat.a[b][a] -= g_eq; }
        }
    }
    fn load_residual_tran(&self, b: &mut [f64], alpha: f64) {
        let g_eq = self.c * alpha;
        let i_hist = g_eq * self.v_prev;
        if let Some(a) = self.nodes[0] { b[a] += i_hist; }
        if let Some(c) = self.nodes[1] { b[c] -= i_hist; }
    }
    fn commit_timestep(&mut self, x: &[f64]) {
        let va = self.nodes[0].map_or(0.0, |i| x[i]);
        let vb = self.nodes[1].map_or(0.0, |i| x[i]);
        self.v_prev = va - vb;
    }

    fn load_residual(&self, _b: &mut [f64]) {}              // DC: open circuit
    fn load_jacobian(&self, _mat: &mut MnaMatrix) {}        // DC: no current
}
```

An inductor is the dual: `L di/dt = v`, integrate to get `v = L·α·(i−i_prev)`.
You'd introduce an auxiliary branch row for the current (via
`num_extra_nodes` returning 1 and `bind_extra_nodes` recording the row),
then stamp `L·α` on the diagonal of that row and `-L·α·i_prev` to the RHS.
See built-in `Inductor` in fairchild-core for the worked example.

The takeaway: **time derivatives are α-scaled stamps + a `commit_timestep`
that snapshots the new state.** That's the same pattern Verilog-A's
`ddt(x)` compiles to.

#### Nonlinear current (diode-like)

For a current `I = f(V)` that's nonlinear in V (diode I-V curve, MOSFET
drain current, photodetector responsivity), Newton-Raphson linearises:

```
I_real(V) ≈ I(V_op) + (dI/dV)·(V − V_op)
```

You stamp `dI/dV` like a conductance in the Jacobian, and an equivalent
current source `I_eq = I(V_op) − (dI/dV)·V_op` on the RHS. This is the
standard "Norton equivalent" companion model. See `ShockleyDiode` in
`crates/fairchild-core/src/models/diode.rs` for a textbook example, or
`NativePhotodetector` in `photonic.rs` for a 2-input version where the
current depends on `V(in_re)² + V(in_im)²`.

The key insight: `eval()` evaluates the I-V curve at the current
operating point and caches `I_op`, `dI/dV`, and `V_op`. `load_residual`
contributes `±I_eq` to the two terminal rows. `load_jacobian`
contributes `dI/dV` like a conductance.

#### Integrals (Verilog-A `idt(x)`)

An integral `y = idt(x)` is the dual of `ddt`: you need an auxiliary state
variable that accumulates `h·x` per step. Introduce one auxiliary MNA row
via `num_extra_nodes() = 1` and `bind_extra_nodes`, then for that row stamp
`+1` against `y` and `-h·α^{-1}` against `x`, with the history term on the
RHS. In practice this is rare in SPICE-style devices — circuits usually
prefer differential equations that the solver integrates implicitly.

#### Voltage-source contributions (Verilog-A `V(p,n) <+ expr`)

For `V(out) = expr`, you need an auxiliary branch row to enforce the
voltage relation as an extra equation, plus an entry coupling the
branch's "current" back to the KCL at `out`. The pattern is:

```
[ 1   ... |  k_i ... ]    branch row: V(out) − Σ k_i · V(in_i) = expr_rhs
[ ...     |  ...     ]
[ ...   1 |  ...     ]    KCL at out: + branch current
```

`fc_waveguide.stamp_potential_eq` is the canonical worked example. The
inputs are `(out_node, [(in_node, coefficient)...])` and the function
takes care of both stamps. You request the branch rows via
`num_extra_nodes()` (returning N branches) and bind them via
`bind_extra_nodes`.

#### Bundle ports (the `(re, im, λ)` photonic convention)

A bundle is just 3 ordinary scalar nets named with `_re`, `_im`, `_λ`
suffixes. The parser's `.optical_port NAME` directive lets the user
write a single token and the parser expands it to 3 wires. Your device
just sees those 3 wires as positions in its terminal list. There's
nothing special at the device layer.

For variable-arity (e.g. an N-channel waveguide), use `Vec<NodeId>` for
nodes, derive N from `terminals.len()` in `setup_instance`, and loop
over channels in `load_jacobian`. `NativeMux` / `NativeDemux` are the
worked examples. The parser side requires adding the device's model
name to the "don't replicate" list in `expand_optical_ports`.

### 13.4 Verilog-A ↔ Rust cheat sheet

| Verilog-A | Rust equivalent in fairchild |
|---|---|
| `parameter real x = 1.0;` | `x: f64` field on the struct, `"x" => { self.x = value; true }` in `set_real_param`. |
| `module foo(p, n);` ports | `nodes: [NodeId; N]` array in the struct, indexed positionally. |
| `analog begin … end` body | `eval()` (read x → cache), `load_residual()` (RHS stamps), `load_jacobian()` (Jacobian stamps). |
| `I(p, n) <+ expr` | Compute I and dI/dV in `eval`, cache; stamp dI/dV like a conductance in `load_jacobian` and `I_eq` in `load_residual`. |
| `V(p, n) <+ expr` | Use a direct-potential branch row (`fc_waveguide.stamp_potential_eq`). |
| `ddt(x)` | `α · (x − x_prev)`. Stamp `α·coef` in `load_jacobian_tran`, history on RHS in `load_residual_tran`, snapshot in `commit_timestep`. |
| `idt(x)` | Auxiliary state variable accumulating `h·x` per step. Rare; use a differential reformulation if you can. |
| `cross(expr, dir)` event | Add a breakpoint via `SimContext` (not yet supported — file an issue). |
| `@(initial_step)` | `setup_instance` runs once; the first `eval` is the equivalent of initial_step. |
| `node_temperature` / `$temperature` | `ctx.temp_k` in `SimContext`. |

### 13.5 Registering the device

Once your `impl Device` compiles, register the SPICE card name in
`crates/fairchild-core/src/device_registry.rs`. Add to `mod.rs`
re-exports first:

```rust
// crates/fairchild-core/src/models/mod.rs
pub use my_family::MyDevice;
```

Then in `register_native_photonics` (or wherever your family belongs):

```rust
self.register("my_device", |terminals, ctx| {
    let mut d = MyDevice::new();
    d.setup_model(ctx);
    d.setup_instance(terminals, ctx);
    Box::new(d) as Box<dyn Device>
});
```

The string `"my_device"` is what the user writes in the netlist:
`X1 a b my_device my_param=2.0`.

### 13.6 Testing

Write a regression test at `crates/fairchild-core/tests/my_device.rs`:

```rust
use fairchild_core::{DeviceRegistry, dc_op_nr_with_registry};
use fairchild_parser::parse_spice;

#[test]
fn my_device_dc_op_matches_closed_form() {
    let netlist = parse_spice(
        "* test\n\
         V1 in 0 DC 1.0\n\
         X1 in out my_device my_param=2.0\n\
         R1 out 0 1k\n\
         .op\n.end\n"
    ).unwrap();
    let r = dc_op_nr_with_registry(&netlist, &DeviceRegistry::new()).unwrap();
    let v = r.node_voltage("out").unwrap();
    let expected = /* your closed-form expression */;
    assert!((v - expected).abs() < 1e-6,
        "V(out) = {v}, expected {expected}");
}
```

Run with `cargo test --release --test my_device`. Once it passes, the
device is part of the standard build and your custom netlists can use it.

### 13.7 Where the change propagates

After adding a device, the only files you touch are:

1. `crates/fairchild-core/src/models/<your_file>.rs` (the device itself)
2. `crates/fairchild-core/src/models/mod.rs` (re-export)
3. `crates/fairchild-core/src/device_registry.rs` (factory)
4. `crates/fairchild-core/tests/<your_test>.rs` (regression)

You do NOT need to touch the parser, the CLI, the Python bindings, or
the netlist analyses — those all dispatch through `DeviceRegistry` and
will pick up your new device automatically once it's registered.

The parser-side `.optical_port` bundle handling DOES need a touch if
your device is variable-arity (treats the bundle as multi-channel
internally) — add its name to the "don't replicate" list in
`expand_optical_ports`. Pure scalar devices don't need this.

If you find yourself wanting a more general extension point (e.g. a
plugin registry that loads compiled `.so` files from outside the
fairchild source tree), that's the OpenVAF / OSDI path described in the
next section. For most new devices, writing a small Rust struct is
simpler and produces faster code.

---

## 14. OSDI model loading

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
