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
14. [Verilog-A models (OSDI)](#14-verilog-a-models-osdi)

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

### BJT (Gummel-Poon Level 1, NPN / PNP)

```
Q<name>  <collector> <base> <emitter> [<substrate>]  <model_name>
```

Substrate node is optional; if absent it is tied to ground internally.

Requires a `.model … NPN (…)` or `.model … PNP (…)` card.

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

### Transmission line (lossless)

```
T<name>  <A+> <A-> <B+> <B->  Z0=<ohms>  TD=<seconds>
T<name>  <A+> <A-> <B+> <B->  Z0=<ohms>  F=<hz> [NL=<wavelengths>]
```

Ideal lossless two-port delay line (Branin's method), characteristic impedance
`Z0` and one-way delay `TD`. Instead of `TD` you may give a frequency `F` with
optional normalised length `NL` (default `0.25`); then `TD = NL / F`. The delay
is intrinsic and always modelled in transient analysis; at DC the line is an
ideal through-connection. No `.model` card — parameters are on the element line.

```spice
T1  in 0 out 0  Z0=50 TD=1n        ; 50 Ω, 1 ns one-way delay
```

(Lossy lines with LTRA-style loss/dispersion are not yet supported.)

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

### BJT Gummel-Poon Level 1 (`NPN` / `PNP`)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `IS` | Transport saturation current (A) | 1e-16 |
| `BF` (`HFE`) | Forward current gain β | 100 |
| `BR` (`HRC`) | Reverse current gain | 1 |
| `NF` | Forward emission coefficient | 1.0 |
| `NR` | Reverse emission coefficient | 1.0 |
| `VAF` (`VA`) | Forward Early voltage (V); `∞` = no Early effect | ∞ |
| `VAR` (`VB`) | Reverse Early voltage (V) | ∞ |
| `TF` | Forward transit time (s) — B-E diffusion charge | 0 |
| `TR` | Reverse transit time (s) — B-C diffusion charge | 0 |
| `CJE` | Zero-bias B-E depletion capacitance (F) | 0 |
| `VJE` | B-E built-in junction potential (V) | 0.75 |
| `MJE` | B-E grading coefficient | 0.33 |
| `CJC` | Zero-bias B-C depletion capacitance (F) | 0 |
| `VJC` | B-C built-in junction potential (V) | 0.75 |
| `MJC` | B-C grading coefficient | 0.33 |
| `FC` | Forward-bias depletion cap linearisation coefficient | 0.5 |
| `RB` | Base ohmic series resistance (Ω) | 0 |
| `RC` | Collector ohmic series resistance (Ω) | 0 |
| `RE` | Emitter ohmic series resistance (Ω) | 0 |

Both transit-time diffusion charges (TF·IF, TR·IR) and depletion junction
capacitances (CJE, CJC) are stamped as companion models in transient analysis.
Non-zero `RB`/`RC`/`RE` add internal collector/base/emitter nodes with the
series resistance to the external terminal (the intrinsic transistor operates on
the internal nodes); the operating point matches ngspice to 6 significant
figures. Only constant `RB` is modelled (no current-dependent `RBM`/`IRB`).
Series resistances RB/RC/RE are accepted and silently ignored (Tier-2 gap).

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

Small-signal sweep around the DC operating point. The reactance matrix includes
both netlist capacitors/inductors **and** device-internal small-signal
reactances at the operating point — diode junction capacitance `Cj(V)`, MOSFET
Meyer gate caps and depletion junction caps, and photonic parasitics — so a
reverse-biased varactor or a transistor's high-frequency rolloff is modelled
correctly. (Through 2026-05-30 these device caps were applied in transient but
omitted from `.ac`/`.noise`; they are now consistent across all analyses.)

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
`noise_sources()` hook. Like `.ac`, the noise small-signal network now includes
device-internal capacitances, so high-frequency noise shaping from device caps
is captured.

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

### Subcircuits (and PCells)

```
.subckt <name>  <port1> <port2> …  [param=default …]
…
.ends [name]

X<inst>  <net1> <net2> …  <name>  [param=value …]
```

Two-pass parse-time flattening. Instances may nest (with cycle detection);
nested *definitions* — a `.subckt` inside a `.subckt` — are rejected. Internal
nets are namespaced `<inst>.<net>`; `0`/`gnd` stays global.

Parameters resolve global `.param` < subckt defaults < call-site overrides, and
are referenced as `{…}`, which may hold **arithmetic** over parameters, numeric
literals, `pi`, and the usual functions:

```spice
.subckt ring in out radius=8e-6 n_g=4.2
Xwg in out fc_waveguide l_m={2*pi*radius} n_g={n_g}
.ends
```

An undefined name, or an expression that evaluates non-finite, is an error —
never a silent zero.

**`.model` cards inside a `.subckt` are per-instance.** The card is name-mangled
to `<inst>.<card>` and references from inside that instance are retargeted, so
each instance carries its own model built from its own parameters. That is what
makes a `.subckt` a real parameterized cell for photonics, where `LEVEL` is only
ever read from a card:

```spice
.subckt mrm_arc a b vpn gnd radius=8e-6 alpha_db_cm=10.7 dn_di=3.99
.model arc_ps fc_pn_th_ps LEVEL=4 l_m={pi*radius} alpha_db_cm={alpha_db_cm}
+ dn_di={dn_di}
Xps a b vpn gnd arc_ps
.ends
```

**Separate files per cell** work through `.include`, resolved relative to the
including deck (`Circuit.load`) or the working directory (`Circuit.load_str`).
See `examples/photonic/pcells/` for `mrm.sp` and `source_bank.sp`, and
`examples/photonic/native_pcell_link.py` for a deck that includes both.

**Optical bundles across the boundary.** A subcircuit's port list is a flat list
of nets, so an optical port costs 3 (or 5) tokens — `in_re in_im in_wl`. When an
instance references a declared `.optical_port`/`.electrical_port`, the port count
of the subckt picks the semantics:

| Declared ports match | Meaning |
|---|---|
| the **flattened** width | one instance carries the whole N-channel bus (a WDM block: N lasers + a mux) |
| the **per-channel** width | replicate, one instance per channel (a single-channel cell along a bus) |

Anything else is a port-count error naming both candidate widths. The two
coincide at N=1, so there is nothing to disambiguate there.

Note the consequence for replication: every element inside is duplicated,
*including electrical ones*. A single-channel ring cell handed an 8-channel bus
becomes 8 rings — right for a ring bank, but if they share an electrical net
they also share 8 parallel junctions. Give each replica its own drive with an
`.electrical_port`, or use a bundle-aware device (`fc_optical_2x2`) when one
electrical interface must serve every channel.

### OSDI (Verilog-A)

```
.osdi  <path/to/model.osdi>
.model <card_name> <module_name> (<params>)
```

`.osdi` loads a compiled OpenVAF-Reloaded shared library and registers every
module in it under its own name. The optional `.model` card binds a second
name to one of those modules with default parameters, which is how foundry
PDKs are written. Instantiate with `X`, or with `D`/`M`/`Q` for a model of
matching arity. Full treatment in §14.

### Bus vectors and optical ports

```
.optical_port  <name>  [<N>]
```

Declares an N-channel optical bundle. Each `name` becomes three (or `3·N`)
underlying wires `name_re[_k]`, `name_im[_k]`, `name_wl[_k]` for SVEA (re, im)
amplitude and per-channel wavelength. Photonic device instances on the bundle
are auto-replicated per channel. See §12.

```
.electrical_port  <name>  [<N>]
```

The electrical sibling: one plain wire per channel, `name_0 … name_{N-1}`.
Like `.optical_port` it makes `name` usable as a *single net token* on an
X-element line, which the parser expands in place. This is how a bundle-aware
device takes one control signal per WDM channel without needing a per-N variant
of the device — see `fc_optical_2x2`. These wires stay in the **electrical**
discipline, so ordinary `V`/`I`/`R`/`C` may drive them.

Both directives create a *port*. Contrast `.optical <net> …` and
`.optical_bus <N> <re_base> <im_base> <wl_base>`, which are **discipline
annotations only** — they tag already-named nets as optical and create nothing
you can reference by a single name. `.optical_port NAME N` does strictly more
than `.optical_bus`; the latter is kept for older netlists.

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
| `solver` | `auto` | `auto` / `dense` / `sparse` / `klu` linear backend (`klu` needs the `klu` build feature) |
| `equilibrate` | false | Two-sided (Ruiz) matrix scaling before LU; improves conditioning of badly-scaled systems, transparent to the solution |
| `cond_estimate` | false | Print a 2-norm condition-number estimate κ(A) of the MNA matrix at the DC operating point |
| `lambda_center_m` | 1.55e-6 | Photonic band-centre default (laser λ, PN-PS reference, waveguide bootstrap). Set via `lambda_center_nm` for nm units. |
| `waveguide_delay` | false | Model optical group delay τ_g = L·n_g/c as a true delay line on every segment-based device — the waveguide **and** the active phase shifters/modulators (default: instantaneous transmission). Aliases: `optical_delay`, `wg_delay`. See §12 `fc_waveguide`. |

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

`.optical_port NAME [N]` declares a bundle. **Every photonic device is
bundle-aware**: a single device instance handles all N optical channels.
This is the rule, not an exception — WDM operation comes from connecting
an N-channel `.optical_port` to a device, not from any per-device opt-in.
The parser dispatches per the centralised `BundleArity` table in
`fairchild-parser`:

- **`BundleArity::Aware`** (default for photonics): `fc_waveguide`,
  `fc_splitter`, `fc_dcoupler`, `fc_grating_coupler`, `fc_pn_ps`,
  `fc_thermal_ps`, `fc_photodetector`, `fc_optical_2x2`, `fc_awgr`. The parser
  flattens every bundle
  into its underlying wires and emits ONE X-element with the combined
  terminal vector; the device's `setup_instance` derives the channel
  count from `terminals.len()`. Pure-optical devices run independent
  per-channel propagation. Devices with electrical state (`fc_pn_ps`,
  `fc_thermal_ps`, `fc_photodetector`) keep ONE shared electrical
  interface — anode/cathode, heat_p/heat_n — so the V_pn supply sees
  one PN junction (not N parallel ones), the photodetector sums
  photocurrents into one anode current with one shared dark current
  and shunt, etc.
- **`BundleArity::Bridge`** (`fc_mux`, `fc_demux`): also flattens, but
  also bypasses the channel-count matching check (N-channel bus side
  ↔ N single-channel pins on the other side).
- **`BundleArity::Scalar`**: the laser (`fc_cw_laser`) and every non-
  photonic device. The parser replicates the X-element into N parallel
  instances when bundles are connected. A laser is fundamentally a
  single-wavelength source — to drive a WDM bus, instantiate one laser
  per channel and combine them through `fc_mux`.

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
| `n_eff` | 2.445 | Effective index at `wl_ref_nm` (sets the accumulated phase per unit length). |
| `n_g` | 4.2 | Group index at `wl_ref_nm` (sets the dispersion slope dn_eff/dλ). |
| `wl_ref_nm` | 1550 | Reference wavelength at which `n_eff` and `n_g` are quoted. Alias: `wl_ref_m` in metres. Defaults from `.options lambda_center_nm`. |
| `alpha_dB_cm` | 2.0 | Power loss (dB/cm). |

The `wavelength_nm` parameter is accepted for backward compatibility but
no longer does anything — the waveguide reads λ directly from the input
bundle's λ wire, and the laser drives that wire to whatever wavelength it
was configured with. A hard-coded 1.55 µm is used only to seed the very
first NR iterate (where the wire is still at 0 V); after iteration 1 the
laser's value wins.

**Physics.** `A_out = A_in · exp(−α L / 2) · exp(−j β L)` with
`β = 2π · n_eff(λ) / λ` and `α` in nepers/m (the `alpha_dB_cm` value is
converted internally). The effective index is first-order dispersion-
corrected from the (`n_eff`, `n_g`) pair at `wl_ref_nm`:

```
n_eff(λ) = n_eff(λ_0) + (λ − λ_0) · (n_eff(λ_0) − n_g(λ_0)) / λ_0
```

so that `n_g(λ_0) = n_eff − λ · dn_eff/dλ` reproduces by construction.
This is the correct physics for ring resonator FSR / Q calculations and
for any wavelength sweep where the propagation phase is the observable.
The `λ` wire is read at evaluation time, so wavelength-dependent
propagation phase is captured automatically — this is what makes the
ring resonator example see a true resonance dip when you sweep the
laser wavelength.

**Group delay (optional).** By default the transmission is instantaneous —
correct for DC and steady-state spectra. With `.options waveguide_delay=1` the
waveguide instead delays its output envelope by the group delay
`τ_g = L · n_g / c`, reconstructed from a per-channel history buffer. Enable it
when the optical modulation bandwidth approaches `1/τ_g` (high-speed links, long
delay lines); leave it off otherwise (cheaper, and the delay is negligible at
low modulation rates). See `examples/photonic/waveguide_delay_demo.sp`.

The corresponding group delay `τ_g = L · n_g / c` is computed at setup
time and stored on the device. It is informational at this tier — the
waveguide currently stamps an instantaneous envelope transfer (no time-
domain delay line); τ_g matters only when modulation bandwidth is
comparable to 1/τ_g (typically tens to hundreds of GHz on chip), which
this device's first-pass model does not yet reproduce. A future
transmission-line device will use τ_g directly.

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

### `fc_optical_2x2` — behavioural per-channel 2×2 transfer block

```
X<name>  in1 in2 thru drop  wctl  ctl_ret  fc_optical_2x2  [param=val …]
```

A 2-in/2-out block whose response you *specify* rather than derive, with an
independent matrix per wavelength channel:

```
[ thru ]   [ s11  s12 ] [ in1 ]
[ drop ] = [ s21  s22 ] [ in2 ]
```

The motivating case is a weight bank — a cascade of ring modulators sharing a
through bus and a drop bus — collapsed to one instance with N weights. That
drops the rings' free parameters and, more importantly, their resonance: with no
resonance the transient timestep is set by your electronics rather than by a
sub-round-trip optical constraint.

Terminals are `4·wpc·N + N + 1`: the four optical bundle ports (all N channels
of each, in port order), then N control wires, then one shared control return.
Declare it with vector ports and the netlist scales by changing one number:

```
.optical_port     in1  4
.optical_port     in2  4
.optical_port     thru 4
.optical_port     drop 4
.electrical_port  wctl 4
Xwb in1 in2 thru drop wctl 0 fc_optical_2x2 w=0 dw_dv=0.4
```

A control bus whose width disagrees with the optical ports is a parse error —
that check lives in the parser, the only layer that still knows each port's
declared width.

| Parameter | Default | Description |
|---|---|---|
| `w`, `w_<k>` | 0 | Bipolar weight, clamped to [−1, 1]. Defined so `P_drop − P_thru = w·P_in`: `−1` all-through, `+1` all-drop, `0` a 50/50 split. Passivity is automatic. |
| `dw_dv`, `dw_dv_<k>` | 0 | Weight per volt on that channel's control wire. This is what moves a weight *during* a transient — `set_param` only works between runs. |
| `s11_mag_<k>` … `s22_deg_<k>` | identity-ish | Explicit matrix entries (magnitude, phase in degrees). Setting any switches that channel out of weight mode. |
| `il_db` | 0 | Extra insertion loss, power dB, applied to every entry. |
| `tau_s` | 0 | Latency (s). Engages a delay line in transient. |
| `allow_gain` | 0 | Permit a matrix with largest singular value > 1. |

Unindexed parameter names broadcast to every channel; the `_<k>` form overrides
one. So `w=0 w_2=0.8` sets channel 2 only.

**Weight mode** builds the lossless coupler-form matrix `s11 = s22 = cos θ`,
`s12 = s21 = −j sin θ` with `θ = ½·acos(−w)` — exactly the `fc_dcoupler` matrix
at `κL = θ`, but with `w` (the number a balanced photodetector pair reads)
as the knob instead of a coupling length.

**Latency caveat.** What `tau_s` delays is the *field*: the output is
`S(t)·in(t − τ)`, matching `OpticalSegment`. A step on a control voltage
therefore reaches the output immediately; only a step on the input field is
delayed. Also note resolving a latency needs a timestep of order `tau_s`, so
leave it at 0 when you want the speed.

**Passivity guard.** An explicit matrix with gain is rejected unless
`allow_gain=1`, because a gain block inside a feedback path diverges silently.
Weight mode is unitary by construction and its clamp keeps it that way under
any control voltage.

**Not yet supported.** Bidirectional mode (`.options enable_bidirectional=1`):
the backward-travelling fields would need their own matrix, and reflection
entries only become meaningful there. The device rejects it outright rather
than leaving the backward wires undriven.

Worked example: `examples/photonic/native_weight_bank.py` (4 channels, balanced
PD pair, weights swept by their control wires inside one `.tran`).

### `fc_splitter` — 1×2 Y-junction (configurable loss + asymmetry)

```
X<name>  in  out_a  out_b  fc_splitter  [param=val …]
```

| Port | Role |
|---|---|
| `in` | bundle, optical input |
| `out_a`, `out_b` | bundles, optical outputs |

| Parameter | Default | Description |
|---|---|---|
| `alpha` | 1.0 | Total intensity transmission (lossless when 1.0). |
| `alpha_dB` (or `il_dB`) | — | Insertion loss in dB; sets `alpha` = 10^(−`alpha_dB`/10). |
| `r` (or `split_ratio`) | 0.5 | Fraction of intensity routed to `out_a`. `out_b` receives `alpha − r`. |

Wavelength duplicated to both outputs. Amplitude transmission to each
arm is √r (for `out_a`) and √(α − r) (for `out_b`). Defaults reproduce
the original 3 dB lossless equal-power split.

### `fc_grating_coupler` — fibre ↔ chip grating coupler

```
X<name>  in  out  fc_grating_coupler  [param=val …]
```

| Port | Role |
|---|---|
| `in` | bundle, optical input |
| `out` | bundle, optical output |

| Parameter | Default | Description |
|---|---|---|
| `alpha_dB` (or `il_dB`) | 3.0 | Insertion loss in dB (amplitude transmission = 10^(−`alpha_dB`/20)). |
| `alpha` | — | Linear amplitude transmission (sets `alpha_dB` = −20·log₁₀(`alpha`)). |

Models a zero-length waveguide with flat amplitude attenuation. No phase
accumulation and no wavelength dependence at this tier; suitable for
testbench inputs/outputs where coupling efficiency is the only physics
that matters. Bundle-aware via parser per-channel replication (pure
optical, no shared electrical).

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
| `r_series` (or `r_s`) | 0.0 | Series resistance (Ω) in line with the anode. Allocates an internal junction node and routes the photocurrent through it. |
| `c_par` (or `c_j0`) | — | Accepted for forward-compatibility; deferred to the L2 PD model with bias-dependent junction capacitance. No effect at this tier. |

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
| `n_eff` | 2.445 | Effective index at `wl_ref_nm` (sets the propagation phase). |
| `n_g` | 4.2 | Group index at `wl_ref_nm` (sets the dispersion slope so `n_eff(λ) = n_eff + (λ−λ_0)·(n_eff−n_g)/λ_0`). |
| `wavelength_nm` | 1550 | Reference wavelength: propagation phase is zero at `λ = wavelength_nm`. Pin this to your laser's λ so the device is "on resonance" by default. (Alias: `wl_ref_nm`.) |
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

Channel count `N` is inferred from instance arity (number of positional nets
minus 1). All parameters are optional.

**Physics.** By default, identity routing per channel: `V(bus_k.re) =
V(ch_k.re)`, `V(bus_k.im) = V(ch_k.im)`, `V(bus_k.λ) = V(ch_k.λ)` for
k = 0..N-1 — a topology marker, not a filter.

Setting any parameter below switches on a **diagonal** spectral response:
channel `k` is scaled by its own passband, evaluated at the wavelength that
channel actually carries. λ labels are never attenuated.

| Param | Default | Meaning |
|---|---|---|
| `il_db` | 0 | insertion loss (flat, if no `fwhm_ghz`) |
| `lambda0_nm` | `.options lambda_center_nm` | grid anchor (channel 0's centre) |
| `df_ghz` | 100 | channel spacing (a **frequency** grid) |
| `fwhm_ghz` | — | passband FWHM; setting it gives each channel a passband |
| `shape_p` | 1 | 1 = Gaussian, 2–4 = flat-top |
| `dlambda_dt_pm_per_k` | 0 | thermal grid drift (silica ≈ 11, SOI ≈ 80) |
| `t_nom_k` | 300.15 | reference temperature for the drift |

```spice
Xmux bus ch0 ch1 fc_mux il_db=3 fwhm_ghz=40 df_ghz=100
```

That gets you insertion loss and the penalty a detuned laser pays on the
passband skirt. What it deliberately does **not** get you is cross-channel
crosstalk — and for a mux that is not a limitation: the N inputs land in N
distinct channel slots, so leakage has nowhere to go. See `fc_demux` for the
case where it *is* a limitation, and `fc_awgr` for the fix.

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

Symmetric counterpart to `fc_mux`: `V(ch_k.*) = V(bus_k.*)`, with the same
optional diagonal response and the same parameter list.

**Why a demux has no crosstalk parameter.** A real demux does leak channel `k`
into output port `j ≠ k` — but an `fc_demux` output port carries a *single*
channel, so representing that leak would mean adding two different carriers
into one complex envelope. That is not allowed: `|a_j + a_k|²` would contribute
a static beat term where the physical ~100 GHz beat is filtered out by any real
photodiode, injecting a spurious DC offset into every downstream detector.
Fields may only be summed within one channel slot.

For a demux **with** crosstalk, use `fc_awgr` with `N−1` input ports left dark.
Its output ports are N-channel buses, which is exactly the somewhere the
leakage needs to live.

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

### `fc_awgr` — N×N arrayed-waveguide grating router

```
X<name>  in_0 … in_{N-1}  out_0 … out_{N-1}  fc_awgr  [params]
```

| Port | Role |
|---|---|
| `in_i` (i = 0..N-1) | bundle, **N-channel** optical input |
| `out_j` (j = 0..N-1) | bundle, **N-channel** optical output |

Every one of the `2N` ports must be declared with `N` channels, giving
`2·wpc·N²` terminals (`wpc` = 3, or 5 bidirectional — which is rejected, see
below). `N` is recovered as `√(len / (2·wpc))` and must come out exact.

**Routing.** Input port `i` channel `k` leaves on output port `(i + k) mod N`,
still in channel slot `k` — the cyclic wavelength shift that makes an AWGR an
all-to-all interconnect: every output receives exactly one wavelength from
every input, with no switching. Taking channel slot index ≡ wavelength index,
the whole device is one complex matrix per slot:

```
out_j[k] = Σ_i  t_ij(λ_{i,k}) · in_i[k]
```

Both crosstalk mechanisms live inside that single form, which is why this
device can model crosstalk and `fc_demux` cannot: wrong-*port* crosstalk
arrives at the same wavelength, so it lands in the same slot and coherently
sums with the wanted signal; wrong-*wavelength* crosstalk lands in its own slot
and stays separate. Neither path ever mixes carriers.

**Modes**, chosen by which parameters are present rather than a mode string:

- **ideal** — nothing set. The exact cyclic permutation, lossless. Routes by
  slot index and deliberately does *not* consult the grid, so an off-grid comb
  still routes rather than going silently dark.
- **gauss** — `fwhm_ghz > 0`. Super-Gaussian passbands on a periodic frequency
  grid, floored by the crosstalk spec.
- **table** — measured `N×N` spectra from CSV, via a `.model` card.

| Param | Default | Meaning |
|---|---|---|
| `lambda0_nm` | `.options lambda_center_nm` | grid anchor: centre of the `(j−i) mod N == 0` pairs |
| `df_ghz` | 100 | channel spacing (a **frequency** grid — an AWG is periodic in f) |
| `fsr_ghz` | `N·df_ghz` | free spectral range; the default *is* the cyclic condition |
| `fwhm_ghz` | — | passband FWHM; positive selects gauss mode, 0 stays ideal |
| `shape_p` | 1 | 1 = Gaussian, 2–4 = flat-top (MMI-input AWGs) |
| `il_db` | 3 (gauss) / 0 (ideal) | peak insertion loss |
| `il_tilt_db` | 0 | extra loss at the outermost channel vs the centre one |
| `xt_adj_db` | −30 | adjacent-channel crosstalk floor |
| `xt_bg_db` | −40 | non-adjacent crosstalk floor |
| `dlambda_dt_pm_per_k` | 0 | thermal grid drift (silica ≈ 11, SOI ≈ 80) |
| `t_nom_k` | 300.15 | reference temperature for the drift |
| `lambda_src` | 0 | which input port's λ tags the outputs mirror |

The crosstalk floors are not decoration: a pure Gaussian tail three channels
out is below −1000 dB, whereas a fabricated AWG floors at −25…−35 dB from phase
errors in the array. `max(gaussian, floor)` reproduces the datasheet shape from
the two numbers datasheets actually quote.

```spice
.optical_port in0 8
* … in1 … in7, out0 … out7 …
Xr in0 in1 in2 in3 in4 in5 in6 in7  out0 out1 out2 out3 out4 out5 out6 out7
+   fc_awgr df_ghz=100 fwhm_ghz=40 il_db=3 xt_adj_db=-30
.options vntol=1e-14 reltol=1e-12
```

**Set a tight `vntol`.** The λ wires carry ~1.55e-6, so the default
`vntol = 1e-6` is the same order as the entire quantity: Newton's step test can
be satisfied while λ is still ~10 pm out, and 10 pm is a real detuning for a
40 GHz passband. The router then reports the transmission for the wrong
wavelength with no error. At N ≤ 5 the first step lands accurately enough that
it never shows, which is what makes it a trap — measured at N = 8, the routed
output reads 0 instead of 1.109 under defaults. Put `.options vntol=1e-14
reltol=1e-12` in any deck containing an `fc_awgr`.

**Measured spectra.** The file path is a string and X-line instance params are
numeric, so a measured router comes in through a `.model` card:

```spice
.model awg8 fc_awgr sfile="awgr8.csv"
Xr in0 … in7 out0 … out7 awg8
```

CSV layout is `wavelength_nm` then `t_<in>_<out>_db` per pair, with optional
`t_<in>_<out>_deg`; rows may be unordered and missing pairs read as dark.

**Also a mux and a demux.** A demux *is* this device with `N−1` input ports
left dark; a mux is it with `N−1` output ports left dark. Not a workaround —
the same star-coupler-plus-array silicon used three ways, which is why an AWG
is cyclic in the first place. Dark ports contribute nothing regardless of their
λ tags.

**Modelling scope.** Transmission is a static coefficient evaluated at each
channel's carrier — the exact narrowband limit, with relative field error
`≈ 2·ln2·(B/FWHM)²` for modulation bandwidth `B`. Insertion loss, detuning
penalty, crosstalk (including its coherent accumulation) and thermal drift are
exact; sideband shaping/ISI, PM→AM off a detuned slope, and channel skew are
not modelled. Analytic modes are purely real, so crosstalk terms all add in
phase — the pessimistic bound. Bidirectional propagation and `tau_s` latency
are rejected/absent. Full discussion, plus the measured cost table, in
`docs/photonic_awgr.md`; worked example in
`examples/photonic/native_awgr_router.py`.

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

### `fc_pn_th_ps` — combined PN + thermal phase shifter

```
X<name>  in  out  anode  cathode  heat_p  heat_n  fc_pn_th_ps  [param=val …]
```

| Port | Role |
|---|---|
| `in` / `out` | bundles, optical pass-through |
| `anode` / `cathode` | scalar nodes — PN-junction terminals |
| `heat_p` / `heat_n` | scalar nodes — heater terminals |

| Parameter | Default | Description |
|---|---|---|
| `L_um` / `L_m` / `length` | 1 mm | Section length. |
| `n_g` | 4.2 | Group index. |
| `dn_dv` | 1e-4 | Effective-index change per applied PN volt. |
| `V_pi_L` | — | Convenience override (V·m) — sets `dn_dv` to give `Vπ·L`. |
| `g_pn` | 1e-3 | PN-junction conductance (S) between anode and cathode. |
| `r_heater` / `r` | 1000 | Heater resistance (Ω) between heat_p and heat_n. |
| `p_pi` / `p_pi_th` | 10e-3 | Heater power for π phase shift (W). |
| `alpha_dB_cm` | 0 | Propagation loss along the section. |

**Physics.** φ_total = φ_prop + φ_eo_pn + φ_th_heater. The two electrical
interfaces are independent — driving only the PN gives `fc_pn_ps`
behaviour, driving only the heater gives `fc_thermal_ps`, driving both
sums the phase shifts. Bundle-aware: ONE physical device handles all N
optical channels with one shared PN junction AND one shared heater
resistor.

### `fc_mzm` — idealised testbench Mach-Zehnder modulator

```
X<name>  in  out  sig  gnd  fc_mzm  [param=val …]
```

| Port | Role |
|---|---|
| `in` / `out` | bundles, optical pass-through |
| `sig` / `gnd` | scalar electrical nodes (modulation drive: `V_mod = V(sig) − V(gnd)`) |

| Parameter | Default | Description |
|---|---|---|
| `V_pi` | 3.0 | Half-wave voltage (DC + AC). |
| `alpha` | 1.0 | Intensity transmission at the bright point (V_mod=0). |
| `alpha_dB` (or `il_dB`) | — | Insertion loss in dB; sets `alpha` = 10^(−`alpha_dB`/10). |
| `e_r` (or `er`) | 1000 | Extinction ratio (linear). |
| `e_r_dB` | — | Extinction ratio in dB; sets `e_r` = 10^(`e_r_dB`/10). |
| `f_c` | 1e10 | First-order EO cutoff frequency (Hz). **Accepted but not yet active** — the V_sig path is instantaneous; f_c lands when device-internal reactive state is wired up. |

**Physics.** Intensity transmission:
```
T(V_mod) = α · [(1 − 1/E_r) · (1 + cos(π V_mod / V_π)) / 2  +  1/E_r]
```
ranges from `α` (bright, V_mod=0) to `α/E_r` (dark, V_mod=V_π). Amplitude
transmission `t_amp = √T`. Use this as the source-side MZM in a
testbench schematic; for a chip-level MZI you'd combine
`fc_splitter` + two `fc_pn_th_ps` arms + a second `fc_splitter`.

### `fc_pn_ps_cap` — depletion-mode PN phase shifter with C_j(V)

```
X<name>  in  out  anode  cathode  fc_pn_ps_cap  [param=val …]
```

Same pin layout as `fc_pn_ps` (bundle-aware optical in/out plus PN
anode/cathode).  Adds bias-dependent junction capacitance and an
optional linear loss-vs-bias coefficient on top of the small-signal
`dn_dv` from `fc_pn_ps`.

| Parameter | Default | Description |
|---|---|---|
| (all `fc_pn_ps` params) | — | `L_um`, `n_g`, `dn_dv`, `g_pn`, `V_pi_L`, `alpha_dB_cm` carry over. |
| `c_j0` | 20 fF | Zero-bias junction capacitance (F). |
| `v_bi` | 0.7 | Built-in voltage (V) — knee at V_pn = V_bi/2. |
| `m_j` | 0.5 | Grading coefficient (0.5 = abrupt junction). |
| `da_dv` | 0 | Linear loss-vs-bias coefficient (Np/m per V); adds extra propagation absorption in reverse bias. |

**Physics.** `C_j(V_pn) = C_j0 / (1 − V_pn/V_bi)^m_j` for V_pn ≤ V_bi/2,
linearly continued above the knee for NR stability when the user drives
the junction into forward bias.  The integrator owns the companion-model
state for this junction capacitance via the new `Device::reactive_branches`
hook (a single Capacitor branch between anode and cathode, value
re-queried per NR iteration).  Forward-injection physics (high dn/dV,
carrier recombination time constants, large da/dV) is reserved for a
future `fc_pn_ps_inj` class.

### `fc_thermal_ps_rc` — thermal phase shifter with τ_th

```
X<name>  in  out  heat_p  heat_n  fc_thermal_ps_rc  [param=val …]
```

Same pin layout as `fc_thermal_ps`.  Adds a first-order thermal RC: the
optical phase shift tracks the *filtered* heater dissipation rather than
the instantaneous Joule power, so transient warm-up / cool-down on the
thermo-optic time scale shows up in the simulation.

| Parameter | Default | Description |
|---|---|---|
| (all `fc_thermal_ps` params) | — | `r_heater` / `r`, `p_pi` carry over. |
| `tau_th` (or `tau`) | 10 µs | Thermal time constant (s). |

**Physics.** `dT/dt = (P − T) / tau_th`, with T in normalised "power-
equivalent" units so steady-state φ = π · T / P_pi matches the L1
model.  In transient, an abrupt change in V_h propagates to T (and
hence to φ) through the LPF.  At DC the state equation reduces to
T = P and the optical output is identical to `fc_thermal_ps`.

This is the canonical "path B" device — T(t) lives as an MNA state row
the device allocates via `num_extra_nodes` and stamps a discretised
state equation in `load_jacobian_tran`.  The previous-timestep T_old
is captured via `commit_timestep`.

### Authoring a custom active photonic device

All the active phase-shifter / modulator classes above share one
architecture (see `crates/fairchild-core/src/models/photonic/`): an
**`OpticalSegment`** (the optical core — per-channel re/im/λ propagation,
bundle stamping, group delay, λ bootstrap) driven by a swappable
**`PhotonicActiveModel`** (the electrical/thermal physics). A device is
`ActiveOpticalDevice { segment, model }`; new physics is a new
`PhotonicActiveModel`, not a new `Device` impl. Its `eval` receives the
node voltages **and the per-channel optical intensity** (so optical→
electrical back-action — photoconductive guides, detectors, TPA self-
heating — is expressible), and returns an `OpticalPerturbation`
(`dn_eff`, `dphi`, `dalpha`) the segment applies.

- **Bias-dependent reactance** (junction `C_j(V)`, electrode cap): return
  it from `PhotonicActiveModel::reactive_branches`; the integrator does
  the BE / TR / BDF-2 companion handling and state advance. `PnDrive`
  (`fc_pn_ps_cap`) is the canonical example, and the same caps now also
  appear in `.ac` / `.noise`.
- **Device-owned discretised state** (a thermal-RC node, carrier density):
  allocate it via `num_internal_nodes`, bind it in `bind_internal`, and
  override `stamp_tran` / `stamp_residual_tran` for the BE state equation.
  `HeaterRc` (`fc_thermal_ps_rc`) is the canonical example.
- **Bolt a heater onto any drive**: wrap it in `WithHeater` (the
  `fc_pn_th_ps*` family) — no new struct needed.

The L3 `FullPnDrive` (`fc_pn_ps_full`) exercises the full surface at once:
implicit junction-voltage solve (series resistance), TPA + self-heating
read from the optical intensity, and depletion + injection regimes.

### Declarative models, no recompile — `fc_phase_shifter_expr`

For the common case — a phase shifter whose constitutive map is a closed-form
function of bias — you can define the physics **on the `.model` card**, no Rust
and no rebuild. The base type `fc_phase_shifter_expr` reads expression-valued
params over the variables `V` (anode−cathode bias), `T` (temperature, K) and
`lambda` (centre wavelength, m), using the same grammar as B-sources (with SPICE
suffixes):

```spice
.model myps fc_phase_shifter_expr
+   dneff  = "-3.1e-5*V - 1.2e-5*V*V"    ; Δn_eff(V)  (free-carrier dispersion)
+   dalpha = "8.0*exp(V/0.7)"            ; excess loss Δα(V)  (Np/m)
+   g_pn   = 1m                          ; junction conductance (numeric)
+   L_um   = 480  n_eff = 2.76           ; optics (numeric, as usual)
Xps in0 out0 a 0 myps pin_at_ref=1
```

`dneff`/`dalpha` are parsed once and evaluated per Newton iterate. Quote the
expression (`"…"` or `{…}`) so spaces and parentheses survive parsing. This is
"Tier 1" of the runtime-loadable-models plan; stateful physics (carrier ODEs,
lookup tables) is future work (scripting / a plugin ABI). Internally it is just
another `PhotonicActiveModel` (`ExprDrive`) on the shared `OpticalSegment`, so it
composes with the rest (e.g. WDM bundles, the `pin_at_ref` convention).

### `fc_circulator` — 3-port bidirectional circulator

```
X<name>  p1  p2  p3  fc_circulator
```

Routes light cyclically: incoming at port 1 exits port 2, incoming at
port 2 exits port 3, incoming at port 3 exits port 1. **Requires
`enable_bidirectional=1`** — without it the device errors out on
instantiation (a circulator is meaningless in unidirectional mode).

Wire convention at every port: `re_fw` / `im_fw` flow INWARD (incoming
to the device); `re_bw` / `im_bw` flow OUTWARD (leaving toward whatever
is connected). λ is tied across all three ports.

Typical use — round-trip monitoring of a device-under-test (DUT):

```
[laser] → port 1 (drive p1_re_fw_*)
port 2 ↔ [DUT]   (forward stimulus + reflected return)
port 3 → [PD]    (PD reads p3_re_bw_*: only the reflection back from p2)
```

No insertion loss or isolation parameters at this tier — the model is
an ideal 3-port circulator.

### Tiered photonic models (separate device classes)

`fc_pn_ps`, `fc_thermal_ps`, and `fc_pn_th_ps` ship today as the basic
linearised first-pass models. Richer physics is exposed by *separate
model names* — not by a `level=` switch. Reasons:

- Tiers aren't strictly linear: a future model may add carrier
  dynamics without adding C_j(V), or vice versa.
- Side variants (forward-only vs. reverse-only PN, photoconductive vs.
  ohmic heater) deserve their own names, not a flag on a base type.
- Selection happens in one place — the KiCad symbol's `model=` field —
  so users still edit one property to choose physics.

Planned upcoming classes (subject to refinement as L2/L3 land):

- `fc_pn_ps_cap` — small-signal `dn/dV` + bias-dependent depletion C_j(V).
- `fc_pn_ps_carrier` — full carrier dynamics, TPA, self-heating from
  linear + nonlinear loss, photocurrent at the electrodes.
- `fc_thermal_ps_rc` — heater with thermal time constant tau_th.
- `fc_thermal_ps_photo` — N-doped photoconductive heater (a different
  physics regime, not just a tier increment of `fc_thermal_ps`).

**SPICE `.model` cards** are supported as a sharing mechanism — drop
a text element on the schematic like

```
.model fast_pn fc_pn_ps_cap c_j0=20f v_bi=0.7
```

and reference it on instances via a `params=fast_pn` symbol field
(separate from `model=`, which still names the device class). This is
optional — instance-side params alone are the common case.

### Bidirectional propagation (`enable_bidirectional`)

By default fairchild's photonic discipline is unidirectional: each
`.optical_port` channel carries 3 wires (`re`, `im`, `λ`) and light
only flows in the direction the device's port topology implies.

Enable bidirectional propagation with

```spice
.options enable_bidirectional=1
```

(or `bidirectional=1`, or via `--opt enable_bidirectional=1` / Python
kwarg). Each channel then carries 5 wires (`re_fw`, `im_fw`, `re_bw`,
`im_bw`, `λ`) and every bundle-aware photonic device stamps an
independent forward path and backward path. The wavelength wire is
shared between the two directions.

Wire-name convention under bidir:

| Direction | Wires |
|---|---|
| Unidirectional (default) | `<port>_re_<k>`, `<port>_im_<k>`, `<port>_wl_<k>` |
| Bidirectional            | `<port>_re_fw_<k>`, `<port>_im_fw_<k>`, `<port>_re_bw_<k>`, `<port>_im_bw_<k>`, `<port>_wl_<k>` |

**Backscattering is not modelled at this tier.** Even with bidir on,
no device couples the forward and backward paths through a scattering
matrix — they're independent. That means you can drive a ring in
either direction and read it out from the appropriate side without
re-wiring the schematic, but reflections off endpoints or coupling
from forward into backward inside a device aren't represented.
`fc_circulator` lets you route light through a directional path on a
testbench; an unconnected port acts as a perfect terminator (no
device drives the floating wires).

### Combined-physics rule for shared-state devices

A bundle-aware device like `fc_pn_ps` represents ONE physical PN
junction interacting with multiple optical modes (channels +, once
implemented, directions). The shared physics MUST integrate across all
modes:

- The single V_pn drives one current — N parallel channels do NOT
  create N parallel conductances. (Already enforced.)
- When future tiers add thermal state from absorbed power, the heat
  generated is the SUM of absorbed power across every channel and
  direction, and the resulting Δn applies to EVERY channel/direction.
- When future tiers add carrier dynamics, the free-carrier density is
  driven by absorption from ALL modes; the resulting Δn and Δα affect
  ALL modes.

This is the natural consequence of "one physical device handles all
the modes" — it shouldn't require special user awareness, just correct
implementation. Bundle-aware stamps already sum currents across
channels into the shared electrical nodes; the same pattern extends to
shared thermal and carrier states.

### Registering PDK-specific aliases

Native devices can be aliased into a PDK-specific taxonomy without
duplicating the device. `DeviceRegistry::register_alias("pdk_foo_waveguide",
"fc_waveguide", remap)` registers a new card name that delegates to
`fc_waveguide` with a parameter-remapping closure. This is the PDK
private-fork extension point — it lets foundry naming and parameter
conventions live downstream without leaking into the upstream master
branch. See `crates/fairchild-core/src/device_registry.rs:300` for an
example.

### OSDI Verilog-A models

OSDI (`.osdi` shared objects compiled by OpenVAF) is the **supported path for
electrical device models distributed as Verilog-A** — foundry transistor models
(BSIM, PSP, HiCUM, …). fairchild does not hand-write BSIM in Rust; load it via
`.osdi <path>` and instantiate with an `X` element. The loader is verified in CI
by the `osdi-mock` fixture.

**Photonics can be written in Verilog-A too** — the complex-envelope
representation is three ordinary real unknowns per channel, so a custom
`optical_field` / `optical_lambda` discipline needs no compiler support and
interoperates exactly with the native devices. See §14.3 and
`examples/verilog_a/`.

What a Verilog-A optical model *cannot* reach is the rest of the abstraction
layer in §12: WDM bundle awareness, bidirectional propagation, `DelayLine`
group delay, and `PhotonicActiveModel` composition. It is single-channel and
forward-only. So prefer the native devices for anything needing those — and in
particular do not start from the pre-Phase-B models under `legacy/`, which are
on the superseded Norton-drive discipline and carry a factor-of-two loss bug
(see `legacy/README.md`). The CLI prints a one-shot hint when one is loaded.

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
name to the `BundleArity` table read by `expand_bundle_ports`.

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
`bundle_arity_for` (read by `expand_bundle_ports`). Pure scalar
devices don't need this.

If you find yourself wanting a more general extension point (e.g. a
plugin registry that loads compiled `.so` files from outside the
fairchild source tree), that's the OpenVAF / OSDI path described in the
next section. For most new devices, writing a small Rust struct is
simpler and produces faster code.

---

## 14. Verilog-A models (OSDI)

fairchild does not parse Verilog-A. You compile it with
[OpenVAF-Reloaded][openvaf] into an OSDI v0.4 shared library (`.osdi`), and
`crates/fairchild-osdi` `dlopen`s it and drives it through the OSDI ABI. This
is the supported route to foundry compact models — BSIM, PSP, HiCUM — which
fairchild will never hand-write in Rust, and it is equally the route to your
own models, electrical or optical.

Worked examples and eight ready models live in `examples/verilog_a/`.

### 14.1 Writing a model

An ordinary Verilog-A module. What the fairchild runtime supports:

| Verilog-A | supported | notes |
|---|---|---|
| `I(a,b) <+ expr` | yes | the resistive branch |
| `V(a,b) <+ expr` | yes | adds an internal branch unknown |
| `ddt(q)` | yes | transient **and** `.ac`/`.noise` |
| internal nodes | yes | `num_nodes > num_terminals`; fairchild allocates the MNA rows |
| `analog function` | yes | |
| `$abstime` | yes | reads the transient clock; 0 in DC and AC |
| `$temperature` | yes | from `.options temp` |
| custom disciplines | yes | OSDI treats them as metadata — see §14.3 |
| `$limit(v, "pnjlim", …)` | yes | installed into the library's `OSDI_LIM_TABLE` at load |
| other `$limit` functions | degrades | identity limiter + a warning; runs unlimited, never crashes |
| `$limexp` | **no** | OpenVAF 23.5 rejects it — a compile error, not a silent no-op |
| `$strobe`, `$finish` | no | |

`ddt` is integrated with whatever `.options method` selects, through the same
`crate::reactive` engine as a discrete `C`. A Verilog-A `ddt(C*V)` and a
netlist `C` of the same value are the same circuit element to machine
precision, under `be`, `tr` and `gear`, fixed step and variable step alike.

That requires the charge itself, not just its Jacobian, so a model must expose
`load_residual_react` — OpenVAF always emits it. A hand-written library that
declares reactive Jacobian entries without it falls back to Backward Euler
rather than stamping a companion with no history behind it.

Junction limiting works the way a compact model expects:

```verilog
Vcrit = Vt * ln(Vt / (sqrt(2.0) * Is));
Vd    = $limit(V(internal, cathode), "pnjlim", Vt, Vcrit);
```

OpenVAF compiles that into a call through the library's exported
`OSDI_LIM_TABLE`, whose entries ship **null** for the simulator to fill in.
fairchild fills every one of them at load time.

The limiter name is not a fixed vocabulary — OpenVAF does not validate it,
it forwards whatever string literal the model wrote — so no simulator can
implement them all. fairchild implements `pnjlim`; anything else (or `pnjlim`
called with an unexpected number of arguments, which would be an ABI mismatch)
gets an **identity limiter** and a warning naming it. That model then runs
*without* limiting for that call: convergence may be slower or need a smaller
step, but the answer is unaffected and the process does not die. Adding a
limiter is one row in `LIMITERS` (`fairchild-osdi/src/loader.rs`).

Limiting changes the Newton path, not the solution: the `$limit` diode in
`examples/verilog_a/` converges to 0.6333213 V against 0.6333214 V for the same
model without it.

A bare `exp()` is fine too — the Newton loop's Armijo line search carries it,
checked to a 500 V drive — but limiting is what scales to stiffer circuits.

### 14.2 Compiling

```bash
openvaf-r -I models models/va_diode.va -o build/va_diode.osdi
```

`-I` sets the include path for `\`include "optical.vams"`. On macOS
`openvaf-r` needs LLVM 18 on the loader path:

```bash
export DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib
```

`examples/verilog_a/build.sh` wraps both. Compiled `.osdi` files are platform
binaries and are gitignored.

Encrypted PDK Verilog-A (IEEE-1735 / Cadence NCPROTECT) is unsupported by
OpenVAF — an upstream Cadence-key problem, not a fairchild limitation.

### 14.3 Optical models

fairchild carries an optical signal on ordinary real-valued MNA unknowns —
three per channel — so a custom optical discipline needs no compiler support
at all; OSDI passes it through as metadata.

| wire | carries | units |
|---|---|---|
| `<port>_re` | real part of the field envelope | sqrt(W) |
| `<port>_im` | imaginary part | sqrt(W) |
| `<port>_wl` | carrier wavelength | m |

Optical power is `re² + im²`. These are exactly the wires a native
`.optical_port p` expands to (`p_re_0`, `p_im_0`, `p_wl_0`), so a Verilog-A
model drops straight into a link built from native `fc_*` devices — address the
underlying wire names on the `X` line. `examples/verilog_a/models/optical.vams`
is the discipline; `wg_compare.sp` pins a Verilog-A waveguide against the
native one in a single deck.

Two rules, both of which cost real debugging time to learn:

**Take the wavelength from a parameter, never off the λ wire.** The wire
exists so a chain stays self-consistent and native devices can read it. Do not
put `OWL(...)` inside an expression you contribute from. Propagation phase is
thousands of radians — 400 µm at n_g = 4.2 is about 6800 rad — so OpenVAF
differentiating it against the λ unknown puts ∂φ/∂λ = φ/λ ≈ 1e9 per metre into
the Jacobian, and Newton fails to converge at some wavelengths and not others.
Native devices dodge this by freezing λ at the previous iterate, which
Verilog-A has no way to express. Sweep wavelength with `--param` or
`set_param` instead.

**`alpha_dB_cm` is power dB, so the amplitude factor divides by 20** —
`10^(−dB/20)`. Everything under `legacy/` predates commit `0f689cb` and is a
factor of two out.

Verilog-A optical models are single-channel and forward-only. WDM bundle
awareness, bidirectional propagation, `DelayLine` group delay and
`PhotonicActiveModel` composition are native-Rust features a Verilog-A model
cannot reach; see §12.

### 14.4 Instantiating

```spice
.osdi build/va_diode.osdi                  ; relative to the netlist file
Xd1  a out  va_diode  Is=1e-14 Rs=0.5      ; model name == module name
```

`.osdi` registers every descriptor in the library under its module name. From
there it resolves like any other model: `X` takes an arbitrary terminal list,
and `D`, `M`, `Q` work for two-, four- and three-terminal models respectively.
All four carry instance parameters.

The foundry idiom puts process parameters on a card and geometry per instance:

```spice
.osdi  bsim4.osdi
.model nch  bsim4 (tox=3n vth0=0.4 …)
M1  d g s b  nch  W=1u L=100n
```

`.model` binds a card name to a loaded descriptor, with the card's parameters
as defaults; an instance parameter overrides the card. Parameter names match
case-insensitively (Verilog-A preserves case, fairchild lowercases). A card
whose second token names no loaded descriptor and no built-in type is left
alone, so a typo surfaces as `unknown model` at the element that uses it.

`--param ELEMENT.PARAM=VALUE` reaches `X`, `R`, `C` and `L` elements only; to
sweep a Verilog-A transistor use the Python bindings' `set_param`, or a
`.param` in the netlist referenced as `{name}`.

### 14.5 Implementation notes

Reactive Jacobian contributions are stamped through the
`write_jacobian_array_react` copy path — `α·∂q/∂x` in transient, `ω·∂q/∂x` in
`.ac` and `.noise`, the same entries in the same positions. The aliasing path
(`load_jacobian_resist`) is broken in OpenVAF-Reloaded and not used.

`OsdiDevice` deliberately does not report `small_signal_reactances`: a
Verilog-A charge is a general matrix, and ∂q_i/∂v_j ≠ ∂q_j/∂v_i —
transcapacitance is exactly what a BSIM-class model is made of, so reciprocal
two-terminal branches would be silently wrong. It overrides
`Device::load_reactive_jacobian` instead.

[openvaf]: https://codeberg.org/arpadbuermen/OpenVAF-Reloaded
