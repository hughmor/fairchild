<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="logos/logo_dark.svg">
    <source media="(prefers-color-scheme: light)" srcset="logos/logo.svg">
    <img alt="fairchild" src="logos/logo.svg" width="360">
  </picture>
</p>

# User guide

Everything needed to write a netlist, run it, and read the result. If you know
SPICE, most of this will be familiar and you can skim to
[§5 Analyses](#5-analyses) and [§7 SimOptions](#7-simoptions-and-convergence-knobs);
the photonics have [their own guide](photonic-models.md).

## Contents

**Writing a netlist**

1. [Netlist syntax](#1-netlist-syntax)
2. [Elements reference](#2-elements-reference)
3. [Waveform sources](#3-waveform-sources)
4. [Model cards](#4-model-cards)

**Running it**

5. [Analyses](#5-analyses) — `.op`, `.dc`, `.tran`, `.ac`, `.noise`, `.tf`, `.sens`, `.pz`
6. [Directives](#6-directives)
7. [SimOptions and convergence knobs](#7-simoptions-and-convergence-knobs)
8. [CLI reference](#8-cli-reference)
9. [Python bindings](#9-python-bindings)
10. [Output formats](#10-output-formats)

**Going deeper**

11. [Solver theory](#11-solver-theory)
12. [Photonic devices](#12-photonic-devices) → [full reference](photonic-models.md)
13. [Writing custom devices](#13-writing-custom-devices)
14. [Verilog-A models (OSDI)](#14-verilog-a-models-osdi)

### Related documents

| | |
|---|---|
| [Photonic models](photonic-models.md) | The optical discipline and every `fc_*` device |
| [SPICE support](spice_support.md) | Every ngspice construct: supported, or how it fails |
| [Model status](model_status.md) | Per parameter: parsed, stamped, validated |
| [Benchmarks](benchmarks.md) | Accuracy and speed against ngspice |
| [KiCad integration](../kicad_integration.md) | Schematic capture |

---

## 1. Netlist syntax

A fairchild netlist follows standard SPICE conventions:

- First line is the title (used in output headers).
- `*` and `;` start comments. `;` may follow content on a line.
- Continuation lines start with `+`.
- Keywords are case-insensitive (`NMOS` = `nmos`).
- SI suffixes: `k`=1e3, `meg`=1e6, `g`=1e9, `t`=1e12, `m`=1e-3, `u`=1e-6,
  `n`=1e-9, `p`=1e-12, `f`=1e-15.
- **RKM is accepted**: the suffix may stand in for the decimal point, so `4k7`
  is 4700 and `2n2` is 2.2 nF — the notation on the part, without transcribing
  it. `m` is the one exception and is a hard error: SPICE reads `m` as milli
  and RKM reads `M` as mega, so `4M7` is 4.7 mΩ or 4.7 MΩ depending on who is
  reading. Write `4.7m` or `4.7meg`. Note that **ngspice does not support RKM**
  — it reads `4k7` as 4000 and drops the `7` without a word — so an RKM deck is
  a fairchild deck, and a deck you intend to run under both should not use it.
- Node `0` (also `gnd`, `GND`) is ground.
- `.end` is optional — end-of-file ends a deck. It is accepted for
  compatibility, but nothing may follow it: a trailing line is an error rather
  than input read off the end of the file and dropped. That is what makes
  "append a line to a working deck and re-run" a workflow you can trust. An
  `.include`d file's own `.end` ends *that file*, not the deck including it.

```spice
* My circuit title
V1  in 0  DC 5
R1  in out 1k
C1  out 0  1u
.tran 10u 5m
```

---

## 2. Elements reference

> For a per-parameter breakdown of what is actually *stamped* versus merely
> accepted, and what each is validated against, see
> [`model_status.md`](model_status.md). Several parameters here are parsed
> for compatibility and do nothing.

### Passive elements

```
R<name>  <pos> <neg>  <resistance>      [m=<n>] [cpar=<F>]
C<name>  <pos> <neg>  <capacitance>     [IC=<v0>] [m=<n>] [esr=<Ω>] [esl=<H>] [rpar=<Ω>]
L<name>  <pos> <neg>  <inductance>      [IC=<i0>] [m=<n>] [rser=<Ω>] [cpar=<F>]
```

`m=<n>` is the instance multiplier — *n of this element in parallel* — applied
exactly: a resistance or inductance divides, a capacitance multiplies, and the
parasitics scale with the copies they belong to. Any other `key=value` on these
lines is an **error** naming the key: ignoring one would leave the element at its
bare value and give a clean answer for a different component.

### Independent sources

```
V<name>  <pos> <neg>  [DC <value>]  [AC <mag> [phase]]  [PULSE|PWL|SIN|EXP|SFFM|AM(...)]
I<name>  <pos> <neg>  [DC <value>]  [AC <mag> [phase]]  [PULSE|PWL|SIN|EXP|SFFM|AM(...)]
```

See [§3](#3-waveform-sources) for waveform shapes.

A source can carry a DC value, an AC small-signal spec, and a transient function
at once — each analysis reads the part that concerns it. The `AC` spec may
appear anywhere after the nodes, and **only** sources carrying one are excited
in `.ac`; see [AC sweep](#ac-sweep).

```spice
V1  in  0  DC 0.7  AC 1        ; biased at 0.7 V, swept at unit amplitude
V2  clk 0  PULSE(0 1.8 0 1n 1n 10n 20n)
```

### Diode

```
D<name>  <anode> <cathode>  <model_name>  [AREA=<n>]
```

Requires a `.model … D (…)` card. `AREA` scales the junction — `IS` and `CJO`
with it, `RS` against it — so `AREA=2` is exactly two of the same diode in
parallel. Any other instance parameter is named on stderr rather than dropped.

### BJT (Gummel-Poon Level 1, NPN / PNP)

```
Q<name>  <collector> <base> <emitter> [<substrate>]  <model_name>  [AREA=<n>]
```

Substrate node is optional; if absent it is tied to ground internally.

Requires a `.model … NPN (…)` or `.model … PNP (…)` card. `AREA` scales the
device: saturation and knee currents and junction capacitances with it, ohmic
series resistances against it.

### MOSFET (Level 1, Shichman-Hodges)

```
M<name>  <drain> <gate> <source> <bulk>  <model_name>  [W=<w>]  [L=<l>]
+        [AS=<a>] [AD=<a>] [PS=<p>] [PD=<p>]
```

Requires a `.model … NMOS|PMOS (…)` card. `AS`/`AD`/`PS`/`PD` are the source and
drain junction areas and perimeters, which is what `CJ` and `CJSW` multiply.

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

### Controlled sources (`E`, `F`, `G`, `H`)

The four linear dependent sources. Two are controlled by a voltage across a node
pair, two by the current through a named voltage source; two produce a voltage,
two produce a current.

```
E<name>  <pos> <neg>  <nc+> <nc->  <gain>     V = gain · (V(nc+) − V(nc−))    VCVS
G<name>  <pos> <neg>  <nc+> <nc->  <gain>     I = gain · (V(nc+) − V(nc−))    VCCS
H<name>  <pos> <neg>  <Vctrl>      <gain>     V = gain · I(Vctrl)             CCVS
F<name>  <pos> <neg>  <Vctrl>      <gain>     I = gain · I(Vctrl)             CCCS
```

```spice
E1  out 0  in 0   3.0        ; ×3 voltage amplifier
G1  out 0  in 0   1m         ; 1 mS transconductance
Vsense  a b  DC 0            ; a current probe …
F1  out 0  Vsense  10        ; … whose current is mirrored ×10
```

`gain` is dimensionless for `E` and `F`, siemens for `G`, ohms for `H`. Each is
stamped every Newton iteration rather than folded into the operating point, so
they work in `.ac` as well as `.op`/`.tran`.

The current sources follow SPICE's convention: positive output current flows
*into* `pos`, through the source, and *out of* `neg`. So `G1 out 0 in 0 1m` with
`V(in) = 1 V` into a 1 kΩ load at `out` gives `V(out) = −1 V`. Swap the output
nodes for the other sign. A quick way to confirm you have it right: a `G` wired
across its own control nodes with `gain = 1/R` must behave exactly as that
resistor.

These are desugared onto the `B` element, which already owns everything they
need. For anything nonlinear, write the `B` element directly.

> `POLY(n)`, `VALUE={…}` and `TABLE` are different elements wearing the same
> letters and mean something quite different. They are refused by name, with a
> message pointing at `B` — rather than being read as a node called `POLY(1)`.

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

### Switches (`S` voltage-controlled, `W` current-controlled)

```
S<name>  <N+> <N-> <NC+> <NC->  <model> [ON|OFF]
W<name>  <N+> <N-> <vsource>    <model> [ON|OFF]
```

A resistor whose value is `RON` or `ROFF` depending on a control quantity —
`V(NC+,NC-)` for `S`, the branch current of a named voltage source for `W`.
The trailing `ON`/`OFF` keyword sets the initial state (default `OFF`).

```spice
.model swmod SW  (VT=2.5 VH=0 RON=10  ROFF=1e9)
.model cs    CSW (IT=1m   IH=0 RON=0.1 ROFF=1e9)
S1  in out  clk 0  swmod OFF     ; sample-and-hold gate
W1  a  b    Vsense cs            ; trips on current through Vsense
```

Switching is a **hard step**, matching ngspice:

```
ctrl > threshold + hysteresis  → ON
ctrl < threshold − hysteresis  → OFF
otherwise                      → hold the previous state
```

Two practical notes, both consequences of that discontinuity:

- **Resolve the timestep.** Pick `h` small enough that a hold capacitor's
  companion conductance (`2C/h` under TR) dominates `1/RON`. Otherwise the
  switched node can move far enough in one step to re-cross the threshold, and
  Newton oscillates between the two states.
- **Use `VH` on feedback paths.** Any switch whose own output can reach its
  control input will chatter with `VH=0`; the hysteresis band is the fix, and
  it is worth reaching for before raising `itl1`.

Inside the hysteresis band the state is held from the **last accepted
timestep**, not the last Newton iterate, so a switch cannot flip-flop within
one NR loop. One consequence: `.dc` sweep points run in parallel and therefore
do not inherit each other's state, so a sweep through the band reports the
instance's `ON`/`OFF` keyword rather than the path-dependent value ngspice
would give. With the default `VH=0` there is no band and no difference.

See `examples/electronic/switched_capacitor.sp`.

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
| `RS` | Series resistance (Ω) | 0 |
| `CJO` (`CJ0`) | Zero-bias junction capacitance (F) | 0 |
| `VJ` | Built-in junction potential (V) | 1.0 |
| `M` (`MJ`) | Grading coefficient | 0.5 |
| `FC` | Forward-bias depletion cap linearisation coefficient | 0.5 |
| `TT` | Transit time (s) — diffusion charge | 0 |

Reverse breakdown (`BV`, `IBV`) is accepted so a foundry card loads, and warned
about: a Zener or an ESD clamp simulates as an ordinary diode with no knee.

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
| `IKF` | Forward high-injection knee current (A); 0 = no roll-off | 0 |
| `IKR` | Reverse high-injection knee current (A) | 0 |
| `ISE` (`C2`) | B-E leakage saturation current (A) | 0 |
| `NE` | B-E leakage emission coefficient | 1.5 |
| `ISC` (`C4`) | B-C leakage saturation current (A) | 0 |
| `NC` | B-C leakage emission coefficient | 2.0 |
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

The base charge carries the Early effect (`VAF`/`VAR`) and the high-injection
knee (`IKF`/`IKR`) together, as in SPICE, so beta rolls off above the knee;
`ISE`/`NE` and `ISC`/`NC` add the non-ideal leakage that pulls beta down at low
current. All of it is checked against ngspice over VBE = 0.4 → 0.9 V.

A parameter this simulator accepts and does not model — `CJS`, `XCJC`, `RBM`,
the transit-time bias modulation, the temperature coefficients, `KF`/`AF` —
produces one warning per card saying what the deck loses.
[Model status](model_status.md) §4 lists every BJT parameter with its true
status; check it before trusting a vendor card.

### MOSFET Level 1 (`NMOS` / `PMOS`)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `VTO` | Threshold voltage (V) | 0.5 / −0.5 |
| `KP` | Transconductance (A/V²) | 20µ |
| `LAMBDA` | Channel-length modulation (1/V) | 0 |
| `GAMMA` | Body-effect coefficient (V^0.5) | 0 |
| `PHI` | Surface potential (V) | 0.6 |
| `CGSO`, `CGDO`, `CGBO` | Gate overlap caps per unit W (or L, for CGBO) | 0 |
| `COX` / `TOX` | Oxide cap density (F/m²) / oxide thickness (m) | 0 |
| `CJ`, `CJSW` | Bulk junction cap per unit area / perimeter | 0 |
| `PB` | Junction potential (V) | 0.8 |
| `MJ`, `MJSW` | Grading coefficient, bottom / sidewall | 0.5 / 0.33 |
| `FC` | Forward-bias depletion cap linearisation coefficient | 0.5 |

Instance `W` / `L` override model defaults.

Only Level 1 exists; a card asking for another `LEVEL` says so and is simulated
as Level 1. The drain and source ohmic resistances (`RD`, `RS`, `RSH`), the
bulk diodes (`IS`, `JS`), velocity saturation (`VMAX`) and the short-channel
corrections are accepted and **not** stamped — each warns, naming what is
missing. For a foundry PDK the answer is the OSDI/Verilog-A path (§14).

### Switch (`SW` / `CSW`)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `RON` | On-resistance (Ω); must be > 0 | 1 |
| `ROFF` | Off-resistance (Ω); must be > 0 | 1e12 |
| `VT` (`SW`) / `IT` (`CSW`) | Threshold on the control quantity | 0 |
| `VH` (`SW`) / `IH` (`CSW`) | Hysteresis half-width; magnitude is used | 0 |

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
with a second source. Names are case-insensitive.

A name that matches no source in the netlist is an error naming the sources you
could have meant, on either axis. It is worth saying explicitly because the
alternative is not obviously wrong: sweeping nothing still produces a table of
the right shape, one operating point repeated down every column.

### Transient

```
.tran  <step>  <stop>  [tstart]  [tmax]  [UIC]
```

Integrates from 0 to `stop`. `UIC` (or `.options uic`) skips DC and uses
`.ic` / element `IC=` values as the initial state; it may occupy any trailing
slot.

| Argument | Meaning |
|---|---|
| `step` | printing/initial step |
| `stop` | end of integration |
| `tstart` | trim the *saved* result to `t ≥ tstart`. Integration still begins at 0 — SPICE's `tstart` selects what is kept, not where the maths starts |
| `tmax` | ceiling on the internal step. Clamps `max_step` with `min`, so it can never loosen a limit set by `.options maxstep` |

The CLI runs every `.tran` card in the deck, each with its own `tstart`, `tmax`
and `UIC`. From Python the card applies only to a `run("tran")` that asked for
no timing of its own — see [Who owns the
run](#who-owns-the-run--the-deck-or-the-caller).

Integration method is selected by `.options method=be|tr|gear` (default TR
with BE first-step bootstrap; GEAR/BDF-2 recommended for stiff circuits or
where ringing on TR is unacceptable). Device-internal capacitances honour this
setting — a diode's `Cj`, a MOSFET's Meyer caps and a BJT's junction caps
integrate the same way a netlist `C` beside them does.

Set `.options variable_step=1` for LTE-controlled stepping. Note that transient
noise requires a fixed step and refuses the variable-step path — see
[Noise](#noise).

### AC sweep

```
.ac  DEC|OCT|LIN  <points>  <fstart>  <fstop>
```

Small-signal sweep around the DC operating point.

**Which sources are driven.** Put an `AC` spec on the source you mean to sweep:

```spice
V1 in 0 DC 0 AC 1        ; magnitude 1, phase 0
V2 in 0 AC 2 90          ; magnitude 2, phase +90°
Vbias b 0 DC 3.3         ; no AC spec — not an AC source, contributes nothing
```

The spec may sit anywhere after the nodes, before or after a DC value or a
transient function. The rule is ngspice's and it is strict: **a source without an
`AC` spec is not an AC source.** A deck with no `AC` spec anywhere is a hard
`SimError::NoAcSource` rather than a quiet zero, because the alternative —
driving everything at unit amplitude — excites DC bias rails as though they were
signal generators, and in a multi-rail circuit that is wrong in a way no single
number reveals.

> Decks written before 0.3.0 relied on the old permissive default. Add `AC 1` to
> the source you intended to sweep.

**What is in the reactance matrix.** Netlist capacitors and inductors, **and**
device-internal small-signal reactances at the operating point — diode junction
capacitance `Cj(V)`, MOSFET Meyer gate caps and depletion junction caps, and
photonic parasitics — so a reverse-biased varactor or a transistor's
high-frequency rolloff is modelled correctly.

An inductor is a short at DC but reactive in `.ac` and `.noise`; both behaviours
come from one assembler told which one it is, so an LC resonance is not flattened
by the DC rule.

### Noise

fairchild has **two** noise analyses over one set of physical sources:

| | What you get | When to reach for it |
|---|---|---|
| **`.noise`** — small-signal | PSDs vs frequency, and a per-source budget | Sensitivity, noise figure, "which device dominates" |
| **`.options trannoise=1`** — time domain | Random currents injected into `.tran` | Eye closure, jitter, comparator/PLL behaviour, anything nonlinear |

They are not alternatives with different answers. `crate::noise::NoiseSources`
enumerates the generators once and both analyses consume that list, so the
time-domain variance of an unfiltered node is exactly the `.noise` PSD
integrated over the resolved band. `transient_noise_agrees_with_the_noise_
analysis` asserts it on three circuits, each biased so a different source
dominates.

Use `.noise` unless the thing you care about is nonlinear or a threshold
crossing. It costs two solves per frequency; a transient-noise run costs a
transient, and you need many of them to estimate a tail.

#### `.noise` — small-signal

```
.noise  V(<out_pos>[,<out_neg>])  <input_src>  DEC|OCT|LIN  <points>  <fstart>  <fstop>
```

```bash
fairchild -f receiver.sp -o out.csv             # onoise / inoise columns
```

```python
r = c.run("noise", out="det", out_neg="bias", src="Vb",
          fstart=1e6, fstop=40e9, points=20, variation="dec")
psd   = r["onoise"]           # V²/Hz at the output port
vrthz = r["onoise_vrthz"]     # V/√Hz, the same thing datasheets quote
inoise = r["inoise"]          # referred back to `src`
```

Adjoint-method: one transposed solve per frequency gives the transfer impedance
from every internal injection to the output, so the cost does not grow with the
number of noise sources. Device noise sources:

| Device | Source |
|---|---|
| Resistor | 4kT/R thermal |
| Diode | 2qI_d shot |
| MOSFET | 8kTg_m/3 channel |
| `fc_photodetector` | 2q(I_ph + I_dark) shot |
| `fc_cw_laser`, `fc_driven_laser` | RIN, `S_P = 10^(rin_dB_Hz/10) · P²` — off unless `rin_db_hz` is set |

Output is `onoise` (V²/Hz at the output port) and `inoise` (equivalent input
PSD referred to `input_src`). OSDI devices can plug in via the `Device` trait's
`noise_sources()` hook. Like `.ac`, the noise small-signal network now includes
device-internal capacitances, so high-frequency noise shaping from device caps
is captured.

> **`.noise` linearises about ONE operating point, and a modulated link does
> not have one.** Shot noise follows the current and RIN follows the optical
> power, so a link's `1` rail can be tens of times noisier than its `0` rail.
> Run `.noise` on a deck whose sources idle at zero and you get the idle
> answer — in `examples/photonic/noisy_eye_and_ber.py` the two rails are 22×
> apart. Bias the deck at the level you care about, one run per rail. This is
> also why the Q-factor formula has two sigmas in it:
> `Q = (µ₁ − µ₀)/(σ₁ + σ₀)`.

**Optical noise.** A laser + PIN + load resistor gives the textbook receiver
budget, and fairchild reproduces it term by term:

```
S_V(f) = ( 4kT/R_L  +  2q·I  +  RIN·I² ) · |Z(f)|²        I = responsivity · P
```

Thermal is flat in power, shot is linear, RIN is quadratic — so RIN is the
floor a link cannot buy its way out of by turning the laser up. Both optical
sources are flat with frequency; the receiver's own poles shape them through
`Z(f)`, which is why device capacitances belong in the small-signal network.

RIN reaches the circuit through *both* field wires at once (one intensity
fluctuation, split by the emission phase), so it is stamped as a single
correlated generator via `Device::correlated_noise_sources`. Treating the two
wires as independent sources would under-report by up to 2× depending on
`phase_deg` — a bug that hides completely at the default 0°.

#### Transient noise — time domain

```
.options trannoise=1 [noiseseed=<n>] [noisescale=<x>]
```

```bash
fairchild -f receiver.sp -o out.csv              # `.options` line lives in the deck
```

```python
r = c.run("tran", step=2e-12, stop=20e-9, trannoise=True, noiseseed=7)
v = r["V(det)"]                                    # now a noisy waveform
```

| Option | Default | Meaning |
|---|---|---|
| `trannoise` | `0` | Inject the noise sources into `.tran`. |
| `noiseseed` | `1` | RNG seed. Same seed ⇒ same waveform; sweep it for independent trials. |
| `noisescale` | `1.0` | Multiplier on the noise **amplitude**, so power goes as its square. `noisescale=3` gives 9× the noise power — the usual way to pull a deep-BER eye closure into a simulation short enough to run, then extrapolate. |

Every generator in the `.noise` table above is drawn once per timestep, held
across the step, and injected as a current at the same terminals the PSD is
defined between. A generator with one-sided PSD `S` is realised as

```
i_n = √(S / 2h) · N(0, 1)          h = timestep
```

because a zero-order-held sequence of variance `σ²` at interval `h` has PSD
`2σ²h` below its Nyquist frequency. Integrating `S` over the resolved band
`[0, 1/2h]` gives back `σ²`, which is the consistency the two analyses share.

**Off by default, and deterministic when off** — a `.tran` is expected to be
reproducible, and every golden in the tree depends on it. A noisy run is still
reproducible: fix `noiseseed` and you get the same waveform every time.

**Fixed step only.** `.options variable_step=1` with `trannoise=1` is an error,
not an approximation: the LTE controller reads a fresh random sample as a fast
signal and shrinks the step to chase it, and the step size then becomes
correlated with the noise, biasing the very spectrum it was meant to reproduce.
Fixed steps are what SDE solvers use, for this reason.

**Bandwidth is set by the timestep, and that is usually fine.** The injected
noise is white up to `1/2h` and absent above it. That truncation does not bias
anything you measure through a circuit, because a transient that resolves your
circuit has its Nyquist frequency well above the circuit's bandwidth and the
circuit filters the difference away — an RC low-pass settles at `kT/C` for any
`h ≪ RC`, independent of `h`, which is the test that pins this. What *does*
depend on the step is a node with no bandwidth limit of its own: a bare resistor
divider has variance `S_V/2h` and always will, in any simulator, because
white noise of unbounded bandwidth has unbounded power. If you are reading a
number that has no filter in front of it, put one there.

**How long to run.** The variance of a variance estimate over `N` independent
samples is `2/N`, and samples are independent only about once per circuit time
constant — so a 1 % variance estimate needs ~20 000 time constants, not
20 000 timesteps. One transient is one realisation of a random process: expect
several per cent of scatter between seeds and pool a few before believing a
number.

**`noisescale` is only exactly linear while the circuit is.** A doubled
amplitude doubles σ as long as the circuit still responds linearly to it. Push
it far enough and it will not — once the injected amplitude is a large fraction
of the signal, a square-law detector stops being small-signal about it. In
`noisy_eye_and_ber.py` the quiet rail tracks `noisescale` to 0.03 % out to 16×
and the loud one to 2 %, because even at 16× the noise is only a few per cent of
that rail. A BER
extrapolated from a scaled run errs conservative when the scaling does bend,
which is the direction you want, but check it rather than assuming it: the
example asserts the slope.

#### Worked examples

| Script | Shows |
|---|---|
| `examples/photonic/receiver_noise_budget.py` | `.noise` only — thermal / shot / RIN vs optical power, and the SNR ceiling `1/(RIN·B)` |
| `examples/photonic/noisy_eye_and_ber.py` | Both — NRZ and PAM-4 eyes through an MZM built from primitives, Q and BER, and each rail checked against its own `.noise` integral *and* the closed form |

#### What neither analysis models

Flicker (1/f) and RTS noise — no device implements them, in either domain.
`.noise` reports PSDs and never injects; transient noise injects and never
reports a PSD. Nothing here is correlated between separate generators, so a
device with physically correlated drain and gate noise would need a single
multi-tap `CorrelatedNoise` (the mechanism exists — laser RIN uses it).

---

### Small-signal reports: `.tf`, `.sens`, `.pz`

These three answer a question about the circuit rather than producing a
waveform, so they return a table instead of a sweep. All three linearise about
the DC operating point.

```
.tf   <v(node[,ref])|i(vsrc)>  <input_source>
.sens <v(node[,ref])|i(vsrc)>  [<element>[.<param>] …]
.pz   <n1> <n2> <n3> <n4>  cur|vol  pol|zer|pz
```

```bash
fairchild -f amp.sp                    # cards run in deck order, like any other
```

```python
c.tf()                                  # the deck's .tf card, whole
c.tf(out="v(out)", src="Vin")           # explicit; needs no card
# {'gain': 0.75, 'r_in': 4000.0, 'r_out': 750.0, 'out_value': 0.75}

c.sens(out="v(out)")                    # every R/C/L/V/I value in the deck
c.sens(out="v(out)", wrt=["r1", "m1.w"])
# [{'param': 'r1.value', 'nominal': 1000.0, 'sensitivity': -1.875e-4,
#   'normalised': -0.1875, 'reached': True, 'fd_error': 3.6e-11}, …]

c.pz(in_pos="in", out_pos="out")        # in_neg/out_neg default to ground
# {'poles': [(-50000+998749.2j), (-50000-998749.2j)], 'zeros': [], …}
```

Pass no analysis arguments and the deck's card is adopted whole; pass any and
the card is not used at all — the same rule `run("tran")` follows, so the
numbers in one result always come from one place.

`run("tf")`, `run("sens")` and `run("pz")` work too and are the *same call* —
`run` delegates to the method, so the two spellings cannot drift apart. What
they do not do is return a `SimResult`: that class is arrays indexed by signal
name, and a report has neither an axis nor signals, so every accessor would
have to answer with an empty array — which is also what a `SimResult` returns
when something went wrong. These return a dict (or a list of dicts) instead.

**`.tf`** gives the gain, the resistance the input source sees, and the
resistance the output port presents. Signs follow ngspice: a branch current
counts positive into a source's `+` terminal, so a driving source reads negative
and a sense source in the return path reads positive.

**`.sens`** is the adjoint, not ngspice's per-parameter re-solve — every
parameter costs one transposed solve between them, and the result is good to
~1e-10 relative rather than to `reltol`. **Read the `reached` flag.** A
parameter the adjoint could not perturb reports `0.0` with `reached = False`,
and a genuine insensitivity reports `0.0` with `reached = True`; a gradient
descent that cannot tell them apart stalls at what looks like a stationary
point. Which model parameters are reachable is `docs/model_status.md`.

**`.pz`** reports roots in rad/s. `vol` drives the input port from a voltage
source (using the deck's own, if one is already there) and `cur` injects current
into an open port — the two give different pole sets, and that is the physics,
not an inconsistency. The eigensolve is dense: past 400 unknowns it is a hard
error naming the limit rather than an unbounded wait.

> **A pole-zero listing is only as linear as the operating point it was taken
> at.** Like `.noise` above, `.pz` and `.tf` describe the circuit *at one bias*.
> A deck idling at zero reports the poles of the idle circuit, which for
> anything with a nonlinear device is not the circuit you care about. Bias the
> deck where you want the answer.

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

An included file may be written in the **Spectre** dialect (`.scs`) rather than
SPICE, and so may the top-level deck — the dialect is detected from the content of
each file, so you never choose a mode and a SPICE deck may pull in a Spectre model
library. Subcircuits (`subckt`/`inline subckt`, with their `parameters` hoisted to
the header), `model` cards, `if`/`else` blocks and `real f(){return …;}` functions
all read across. Spectre statements are transliterated to their SPICE equivalents, which is
why everything else in this guide still applies unchanged; the surface currently
read, and what it refuses, is tabulated in
[`spice_support.md` §5](spice_support.md).

### Parametrisation

```
.param  <name>=<value>  [<name2>=<value2> ...]
.func   <name>(<arg> [, <arg> ...]) = <expression>
```

Parameters substitute into element values, model-card values and B-element
expressions as `{name}`, and a value may be any expression rather than only a
number:

```spice
.param w=2u l=0.18u
.param area={w*l}          ← braces
.param perim='2*(w+l)'     ← single quotes: the HSPICE spelling, same meaning
.func  ratio(a,b) = a/b
.param aspect={ratio(w,l)}
M1 d g s b nm W={w} L={l}
.model nm NMOS (VTO={0.7*1.05} KP=100u)
```

Values resolve in file order over the parameters already defined, including
earlier on the same line, so define before you use. A value may contain spaces
only inside braces or quotes — `.param a = 1 + 2 b = 3` has no unambiguous
reading. A SPICE suffix stays a number: `1k` is 1000, not an expression over an
undefined `k`.

`.func` is expanded where it is called, at parse time, so it works anywhere an
expression does — a `{…}` value, another `.param`, a `.model` value, a B-source,
a `.measure`. A definition may come after its first use. An argument shadows a
`.param` of the same name inside the body. Recursion, a repeated argument name, a
name that shadows a built-in function such as `sin`, and a call with the wrong
number of arguments are all errors.

Double quotes mean something different and are not substituted: on a `.model`
line `"…"` is a device constitutive map over the device's own bias
(`dneff="5.0e-5*V"`), which is not a parse-time value. See
[Photonic devices](#12-photonic-devices).

### Conditional netlists

```
.if (<condition>)  …  .elseif (<condition>)  …  .else  …  .endif
```

The condition is a parse-time expression over `.param` values and `.func` calls,
with the comparison and logical operators of the expression grammar. Only the
taken branch is collected: a `.model`, `.subckt`, `.param` or element in a branch
that was not taken does not exist afterwards. Blocks nest, and a condition in a
branch that cannot run is not evaluated, so a dead branch may reference names that
were never defined.

```spice
.param corner=2
.if (corner==1)
.include models_tt.lib
.elseif (corner==2)
.include models_ff.lib
.else
.include models_ss.lib
.endif
```

Inside a `.subckt` the condition is evaluated **per instance**, against that
instance's parameters — so a switch on a subcircuit parameter selects for each
instance, and a dead branch is dropped whole for that instance:

```spice
.subckt rsel a b mode=0
.if (mode == 1)
.param r=2k
.else
.param r=1k
.endif
R1 a b {r}
.ends
Xa in 0 rsel mode=1     ← 2 kΩ
Xb in 0 rsel            ← 1 kΩ
```

One thing is refused rather than guessed:

- **A condition over an undefined name.** `nope==1` would compare NaN against 1
  and yield a perfectly ordinary `false`, so a misspelled corner variable would
  quietly select the other branch. Names are checked before the condition is
  evaluated.

### Corner and temperature sweeps

```
.temp <T1> [<T2> ...]
.alter <label>
   ... overrides ...
.endalter
```

`.temp` re-runs every analysis once per listed temperature (°C). `.alter`
blocks describe deltas from the base netlist; each block produces a full
re-run with overrides applied. The two cross: a deck with two `.alter` blocks
and three temperatures has **nine** corners (base plus two blocks, times three),
and every analysis in the deck runs at each of them.

A single `.temp 75` is not a sweep — it just sets the temperature, so the deck
still has one corner.

The CLI runs the whole grid, one output file per corner (see
[§8](#8-cli-reference)). From Python it is the same grid, expanded by the same
code, reachable two ways:

```python
c.corners              # [{'alter': 'base', 'temp_c': -40.0}, …] — what the deck declares
c.run_all()            # every analysis at every corner; one row per (corner, analysis)

for corner in c.corners:                      # or drive them yourself
    r = c.run("tran", alter=corner["alter"], temp=corner["temp_c"])
```

`run_all()` returns rows of `{'alter', 'temp_c', 'kind', 'result'}`. Each row
names its corner because the results are otherwise indistinguishable — two
`.tran` rows from a two-corner deck are the same analysis at different
temperatures, and nothing inside a `SimResult` says which. Corners are
independent, so they run in parallel.

`run_all(alter=…)` and `run_all(temp=…)` are **errors**, not filters: asking the
run-everything verb for one corner is a contradiction. Use `run()` for that.
Solver options (`reltol`, `method`, …) do apply to every corner.

### Solver options

```
.options  <key>=<value>  ...
```

Any field of `SimOptions` ([§7](#7-simoptions-and-convergence-knobs)) can be set
this way, under the key the §7 table lists — which is not always the field's
own name: temperature is `temp`, in °C, not the kelvin field `temp_k`. An
unrecognised key warns and has no effect.

### Subcircuits (and PCells)

```
.subckt <name>  <port1> <port2> …  [param=default …]
…
.ends [name]

X<inst>  <net1> <net2> …  <name>  [param=value …]
```

Two-pass parse-time flattening. Instances may nest (with cycle detection);
nested *definitions* — a `.subckt` inside a `.subckt` — are rejected. Internal
nets are namespaced `<inst>.<net>`; `0`/`gnd` stays global. Node and branch
references inside an expression are in the same scope as the element that holds
them, so a `B`/`E`/`F`/`G`/`H` inside a subcircuit reads its own instance's nets
and sources.

```
.global  <net1> <net2> …
```

A `.global` net is the same node in every scope, whether or not it appears in a
port list — the supplies of a CDL or foundry deck, which no layout tool threads
through every port list. Nesting is unlimited, and the declaration may come after
the instance that needs it. A net that is **both a port and global is refused**:
the port would take the caller's net while every reference inside took the global
one, and picking either silently is wrong for the deck that meant the other.

```spice
.global vdd vss
.subckt inv a y          ← no supply ports
Rpull vdd y 1k
Rdown y vss 1k
.ends
Vsup vdd 0 1.8
Xi in out inv            ← Rpull connects to the top-level vdd
```

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

**A subcircuit's parameters resolve per instance**, in this order: the enclosing
scope, then the header defaults with the call's overrides already in place, then
the body's `.param` lines. So a default may be an expression over an earlier
parameter, a body `.param` follows an override rather than the default, and two
instances of one definition resolve independently:

```spice
.subckt rdiv a b n=1 rsh='1000/n'   ← a default over another parameter
.param rtot={rsh*n}                 ← computed per instance, not per definition
R1 a b {rtot}
.ends
X1 in 0 rdiv n=2                    ← rsh=500, rtot=1000
```

**`m=` on the instance** means m of the whole subcircuit in parallel, and scales
everything the body flattened to (see [Passive elements](#passive-elements)).
`m` is the simulator's parameter, so it needs no declaration — but a definition
that *declares* `m` owns it, and nothing is scaled here, because a wrapper that
forwards `m` to the device inside is already doing the scaling. An element with no
exact scaling — a diode, a MOSFET, a switch — is an error naming it, since a factor
quietly ignored is a wrong answer exactly the size of the factor.

Two more things are refused rather than guessed:

- **A parameter the definition does not declare.** `X1 in 0 rdiv nn=2` is an error
  naming `nn` and listing what `rdiv` declares — a typo that left the default in
  place would report a clean answer for a different circuit.
- **Overriding a body `.param`.** It is computed, not an interface: overriding it
  and recomputing it are different circuits. Move it to the header to make it
  overridable.

An instance-parameter value must also be readable — a number, or `{…}` that
resolves to one. `n=2*3` unbraced is an error, not a dropped assignment.

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
matching arity. Full treatment in [§14](#14-verilog-a-models-osdi).

### Bus vectors and optical ports

```
.optical_port  <name>  [<N>]
```

Declares an N-channel optical bundle. Each `name` becomes three (or `3·N`)
underlying wires `name_re[_k]`, `name_im[_k]`, `name_wl[_k]` for SVEA (re, im)
amplitude and per-channel wavelength. Photonic device instances on the bundle
are auto-replicated per channel. See [Photonic models](photonic-models.md).

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
| `temptol` | 1e-3 | NR temperature absolute tolerance (K), on thermal rows |
| `vmax` | 0.5 | Per-iteration |ΔV| clamp |
| `gmin` | 1e-12 | Diagonal regularising conductance (S) |
| `gminmax` | 1.0 | GMIN-stepping starting value |
| `itl1` | 150 | DC max NR iterations |
| `itl4` | 150 | Transient per-step max NR iterations |
| `max_rejections` | 30 | Var-step max step rejections |
| `method` | `tr` | `be` / `tr` / `gear` |
| `max_step` | ∞ | Transient max step (s) |
| `tstart` | 0 | Discard transient output before this time (s); the run still integrates from 0. Also `.tran`'s third positional argument |
| `variable_step` | false | LTE-controlled variable-step transient; `step` becomes the initial/maximum stride. Alias: `variablestep`. CLI `--variable-step` |
| `srcsteps` | 10 | Source-stepping homotopy resolution. Alias: `srcmax` |
| `temp` | 27 | Operating temperature in **°C** (stored internally as kelvin; the struct field is `temp_k`, but `temp_k` is not an accepted key). `tnom` takes the same units |
| `uic` | false | Use `.ic` / element `IC=` instead of DC |
| `pnjlim` | true | Diode / MOSFET junction limiting in NR |
| `solver` | `auto` | `auto` / `dense` / `sparse` / `klu` linear backend (`klu` needs the `klu` build feature; `auto` picks dense below ~20 nodes, then `klu` when built in, `sparse` otherwise) |
| `equilibrate` | false | Two-sided (Ruiz) matrix scaling before LU; improves conditioning of badly-scaled systems, transparent to the solution |
| `cond_estimate` | false | Print a 2-norm condition-number estimate κ(A) of the MNA matrix at the DC operating point |
| `lambda_center_m` | 1.55e-6 | Photonic band-centre default (laser λ, PN-PS reference, waveguide bootstrap). Set via `lambda_center_nm` for nm units. |
| `trannoise` | false | Inject the `.noise` sources into `.tran` as random currents. Fixed step only. Aliases: `tran_noise`, `transient_noise`. See [Noise](#noise). |
| `noiseseed` | 1 | Transient-noise RNG seed. Same seed ⇒ same waveform. |
| `noisescale` | 1.0 | Multiplier on transient-noise **amplitude** (power goes as its square). |
| `enable_bidirectional` | false | Bundles carry forward **and** backward fields (5 wires/channel instead of 3); reflective devices become meaningful. Aliases: `bidirectional`, `bidirectional_propagation`. See [§14](#14-verilog-a-models-osdi) |
| `sanity_check` | true | Netlist preflight pass (R=0, duplicate refdes, zero-parameter `fc_*`, …) warning before analysis. Disable with `nosanitycheck` / `sanity_check=0` |
| `verbose` | false | Solver progress notes to stderr: matrix size/NNZ, which convergence phase ran, top residual rows on NR failure. CLI `-v` |
| `waveguide_delay` | false | Model optical group delay τ_g = L·n_g/c as a true delay line on every segment-based device — the waveguide **and** the active phase shifters/modulators (default: instantaneous transmission). Aliases: `optical_delay`, `wg_delay`. See [`fc_waveguide`](photonic-models.md#fc_waveguide--lossy-waveguide). |

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
| `-o, --output <FILE>` | Output destination (default stdout). With more than one `.alter` × `.temp` corner, the path becomes a *stem* and one file per corner is written (e.g. `out.alter_pvtfast.temp_-40c.csv`), the corners running in parallel — unless `--single-output` or `--verbose` |
| `--single-output` | Bundle all corner outputs into the single `--output` file (with `# alter=…` / `# temp_c=…` header lines), run serially in deterministic order |
| `--probe <SIG,…>` | Comma-separated CSV signal filter. An unmatched name is an error; `V(node)` selects an AC sweep's `mag_`/`phase_deg_` pair |
| `--param ELEM.PARAM=VAL` | Override a circuit parameter (repeatable) |
| `--opt KEY=VAL` | Override a SimOptions field (repeatable) |
| `--reltol`, `--gmin`, `--method`, `--maxstep`, `--solver` | Convenience flags |
| `--variable-step` | LTE-controlled variable-step transient (same as `.options variable_step=1`) |
| `--no-pnjlim` | Disable junction limiting |
| `--check` | Parse + discipline-check only |
| `--list-nodes` / `--list-models` | Inspect parsed netlist, then exit |
| `-v` / `-q` | Verbose / quiet. `-q` silences every warning, including the ones raised inside the parser and the solver — a skipped `.control` block, an ignored `.print`, an unrecognised `.options` key, a MOSFET card asking for an unimplemented level. Errors are unaffected |

Examples:

```bash
# DC operating point
fairchild -f netlist.sp

# Transient with GEAR + tightened tolerance + sparse LU, into Nutmeg
fairchild -f netlist.sp \
   --opt method=gear --opt reltol=1e-5 --solver sparse \
   --format nutmeg -o out.raw

# AC sweep specified entirely via netlist; CSV-filter on one probe
# (selects mag_V(out) and phase_deg_V(out); an AC CSV carries no currents,
# so probing I(V1) here would be refused by name)
fairchild -f rlc.sp --probe "V(out)"

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
ac    = c.run("ac", variation="dec", points=20, fstart=1, fstop=1e6)
noise = c.run("noise", variation="dec", points=20, fstart=1, fstop=1e6,
              out_pos="out", src="V1")

# The small-signal reports return tables, not waveforms.  Either spelling works
# — c.run("tf") delegates to c.tf() — but both return a dict, not a SimResult.
# Each takes the deck's card when given no arguments of its own.
gain  = c.tf(out="v(out)", src="Vin")["gain"]
grads = c.sens(out="v(out)")                      # check each row's ['reached']
poles = c.pz(in_pos="in", out_pos="out")["poles"] # rad/s, complex

# What the deck declares, without running any of it:
c.analyses          # [{'kind': 'tran', …}, {'kind': 'pz', …}]
c.corners           # [{'alter': 'base', 'temp_c': 27.0}, …]

# Run the whole deck — every analysis at every corner, like the CLI does.
for row in c.run_all():
    print(row["alter"], row["temp_c"], row["kind"], row["result"])

# Or drive the corners yourself, one at a time:
for corner in c.corners:
    r = c.run("tran", alter=corner["alter"], temp=corner["temp_c"])

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

#### Who owns the run — the deck, or the caller

A deck *declares* what could be run; from Python it never runs anything by
itself. `Circuit.analyses` is that declaration, and `run(kind)` chooses from it:

```python
c.analyses
# [{'kind': 'op'},
#  {'kind': 'tran', 'step': 1e-12, 'stop': 5e-09, 'tstart': 0.0,
#   'tmax': 1e-13, 'uic': False}]

c.run("tran")                          # the deck's .tran card, whole
c.run("tran", step=1e-12, stop=5e-9)   # your numbers; the card is not read

for a in c.analyses:                   # what the CLI does, in deck order
    c.run("dc_sweep" if a["kind"] == "dc" else a["kind"])
```

A card is taken **as a unit or not at all**: `run("tran")` with no timing takes
`step`, `stop`, `tstart`, `tmax` and `UIC` off the card; passing either `step` or
`stop` means you own the timing and none of the card applies. Half a card is
never applied, so every number in a run comes from one place. `.ac`, `.noise` and
`.dc` (through `run("dc_sweep")`) work the same way.

Two errors rather than a guess: no card *and* no kwargs, and a deck declaring two
cards of the same kind. `"dc"` remains an alias for `"op"` — a `.dc` card is
reached with `run("dc_sweep")`.

Directives outside that list divide the same way. `.model`, `.param`, `.subckt`,
`.include`, `.ic`, `.osdi` and friends define *what the circuit is* and are
always honoured identically in both frontends. `.print`, `.plot`, `.probe`, `.save` and
`.width` select *what to report*, which the frontend owns — use `--probe` from
the CLI, or index the returned result. They load and warn rather than narrowing
anything: every signal is available either way. `.control` is imperative script —
`run`, `let`, `write`, loops — and control flow belongs in this API, so the block
is skipped with a warning naming the commands in it rather than interpreted. That
is permanent, not pending: if the deck's only analysis was a `tran` command inside
the block, declare a `.tran` card or pass the timing here instead.

### Gradients — `dc_adjoint`, `ac_adjoint`, `tran_adjoint`

Each of the three analyses can report **`dL/dp`**: how a scalar you care about
moves when a design parameter moves. The cost is one extra solve regardless of
how many parameters you ask about, which is what makes it usable inside an
optimiser.

```python
# DC — value and gradient from a single solve.
r = c.dc_adjoint(probes={"p": ("power", "bar", 0)},
                 wrt=["Vh.dc"],
                 params={"Vh.dc": 2.5})       # per-run overrides
r.values["p"], r.grad["p"][0]

# AC — one backward pass covers the whole sweep.
run = c.ac_adjoint("out", quantity="mag2", variation="dec", points=8,
                   fstart=1e8, fstop=1e11, src="V1")
y = np.asarray(run.response)                   # per frequency; a property
g = run.backward(2 * (y - target),             # dL/dy, one weight per point
                 ["R1.r", ("Xps.l_um", 0.1)])  # (name, step) pins the FD step

# Transient — a functional of the waveform.
run = c.tran_adjoint({"v": "out",                  # node voltage
                      "p": ("power", "out0", 0)},  # optical power, W
                     step=1e-11, stop=2e-9)
y = run.probes["v"]
g = run.backward({"v": 2 * (y - target)},          # dL/dy per timepoint
                 ["R1.r", "C1.c"])
```

`wrt` and `params` are separate on purpose: `params` overrides values for this
run, `wrt` names what to differentiate with respect to.

Three things to know before you rely on the numbers:

- **The result is a partial with respect to a netlist parameter.** If your
  design variable moves several netlist parameters at once — an optical length
  that sets both phase and junction capacitance — the chain rule onto the design
  variable is yours to apply, or JAX's if you wrap it in `custom_vjp`.
- **`reached` is reported, not assumed.** A misspelled parameter, or a device
  that ignores the write, comes back flagged rather than as a confident zero.
- **Check the gradient, and not at the optimum.** The gradient is zero there, so
  agreeing about zero proves nothing.
- **`tran_adjoint` needs a fixed timestep.** `variable_step=1` is refused, not
  downgraded: the co-state recursion replays the forward pass's step sequence,
  and an adaptive controller would re-decide that sequence under a perturbed
  parameter — differentiating a schedule nobody solved. Compare against a finite
  difference from a fixed-step `run("tran")`, not an adaptive one.

Worked examples for all three, each checking itself against a full re-solve:
[`examples/optimization/`](../examples/optimization).

---

## 10. Output formats

### CSV

Default. One row per timepoint (transient), frequency (AC/noise), or sweep
point (DC sweep). Header row is comma-separated signal names. `--probe`
filters the column set; a probe that matches no column of an analysis's
output is an error naming it, never a silently thinner CSV. For an AC sweep,
`V(node)` selects that node's `mag_`/`phase_deg_` column pair.

### Nutmeg rawfile

ASCII-only format compatible with ngspice's `rawread` and the `spyci` Python
library. Always emits the full signal set.

```bash
fairchild -f circuit.sp --format nutmeg -o waveforms.raw
ngspice -c "rawread waveforms.raw; plot V(out)"
```

Every analysis in the deck writes one plot, appended to the same file in deck
order, using ngspice's plot names so a reader can classify them:

| Analysis | `Plotname:` | `Flags:` | Sweep variable |
|---|---|---|---|
| `.op` | `Operating Point` | `real` | — (one point) |
| `.dc` | `DC transfer characteristic` | `real` | sweep source name |
| `.ac` | `AC Analysis` | `complex` | `frequency` |
| `.tran` | `Transient Analysis` | `real` | `time` |
| `.noise` | `Noise Spectral Density Curves` | `real` | `frequency` |

**`.noise` units differ between the two formats, deliberately.** The rawfile
emits `onoise_spectrum` / `inoise_spectrum` as amplitude density in **V/√Hz** —
that is what ngspice puts under those names, and what a reader assumes. The CSV
gives both, with the units in the column names (`onoise_v2hz`, `onoise_vrthz`).
Where the input-referred PSD is not computable (the transfer function is too
small to invert) both formats write `NaN` rather than `0`.

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
- `|ΔT| < temptol + reltol · |T|` for thermal rows (see below).
  The tolerance is per *quantity*, not one number for the whole solution
  vector — see `crates/fairchild-core/src/tolerance.rs`. There used to be a
  third class, `lambdatol`, for optical wavelength wires; λ is now resolved
  before the solve rather than solved for, so there is no λ row to converge and
  `.options lambdatol=…` warns that it no longer applies.

### Thermal nodes

A Verilog-A model that declares a node `thermal` gets an unknown whose potential
is a **temperature rise above ambient in kelvin** and whose flow is a **power in
watts**:

```verilog
`include "disciplines.vams"

module self_heated_r(p, n, h);
    inout p, n;  electrical p, n;
    inout h;     thermal h;
    parameter real r0 = 1000.0, alpha_r = 4.0e-3;
    real r_now;
    analog begin
        r_now = r0 * (1.0 + alpha_r * Temp(h));
        I(p, n) <+ V(p, n) / r_now;
        Pwr(h)  <+ -V(p, n) * V(p, n) / r_now;   // watts out of the device
    end
endmodule
```

Nothing is declared in the deck. OSDI carries the discipline's units through to
the descriptor, so fairchild reads which rows are thermal off the model itself
and bounds them with `temptol` (1 mK) instead of `vntol` (1 µV) — a microvolt is
a nonsense convergence bound on a temperature. This works for a thermal port and
for a purely internal self-heating node alike.

**The thermal network is written in ordinary SPICE primitives**, because on a
thermal node they already mean the right thing:

```
X1  p 0 h  self_heated_r
Rth h 0 500      ; 500 K/W to ambient — R stamps ΔT/R watts
Cth h 0 1u       ; 1 J/K of heat capacity
```

`R` is a thermal resistance in K/W, `C` a heat capacity in J/K, `I` a heat
source in watts, and `V` a clamped temperature. Two devices sharing one thermal
net are thermally coupled; an `R` between two thermal nets is crosstalk. There is
deliberately **no** discipline check refusing electrical elements here, unlike on
optical wires — the electrical primitives *are* the thermal primitives.

`Temp` is a rise, not an absolute temperature. Ambient comes from `$temperature`
(set by `.options temp=`, in °C), so absolute temperature is
`$temperature + Temp(h)`. Keeping them apart is what lets a deck sweep ambient
while self-heating stays solved on top.

`examples/verilog_a/models/mrm_wdm.va` is a worked example: a microring whose
resonance moves with a solved temperature, with `th` exposed so a deck can wire
ring-to-ring thermal crosstalk.

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

The photonic discipline has its own reference: **[Photonic
models](photonic-models.md)**. It covers how an optical field is represented as
MNA unknowns, bundle ports and WDM, bidirectional propagation, every `fc_*`
device with its parameters, the phase-shifter tiers, optical noise, and what is
deliberately not modelled.

The rest of this guide applies unchanged — a photonic device is an `X` element,
it is solved by the same Newton loop, and every analysis, option and output
format below works on a circuit containing one.

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
         .op\n"
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

fairchild does not parse Verilog-A. It is compiled by
[OpenVAF-Reloaded][openvaf] into an OSDI v0.4 shared library (`.osdi`), which
`crates/fairchild-osdi` `dlopen`s and drives through the OSDI ABI. This is the
supported route to foundry compact models — BSIM, PSP, HiCUM — which fairchild
will never hand-write in Rust, and it is equally the route to your own models,
electrical or optical.

**You do not have to run the compiler yourself.** Name the source and fairchild
compiles it, caching the result:

```spice
.va models/va_diode.va          ; Verilog-A source — compiled on demand
.model nch va_diode (Is=1e-14)
M1 d g s b nch
```

Spectre's `ahdl_include "va_diode.va"` does the same thing, so a foundry PDK
loads as written. The explicit two-step route also still works, and is what
belongs in CI:

```spice
.osdi build/va_diode.osdi       ; a compiled artefact — no toolchain needed
```

Which keyword you use does not matter: a path ending `.va` or `.vams` is
treated as source and a path ending `.osdi` as an artefact, whichever of `.va`
/ `.osdi` named it.

Worked examples and nine ready models live in `examples/verilog_a/`.

### 14.1 Writing a model

An ordinary Verilog-A module. What the fairchild runtime supports:

| Verilog-A | supported | notes |
|---|---|---|
| `I(a,b) <+ expr` | yes | the resistive branch |
| `V(a,b) <+ expr` | yes | adds an internal branch unknown |
| `ddt(q)` | yes | transient **and** `.ac`/`.noise` |
| `white_noise(pwr, "name")` | yes | reaches `.noise` **and** `.options trannoise=1` |
| `flicker_noise(k, af, …)` | ⚠️ | exact in `.noise`; transient injects white samples and can only realise one density, so it probes mid-band and warns |
| internal nodes | yes | `num_nodes > num_terminals`; fairchild allocates the MNA rows |
| `analog function` | yes | |
| `$abstime` | yes | reads the transient clock; 0 in DC and AC |
| `$temperature` | yes | from `.options temp` |
| custom disciplines | yes | OSDI treats them as metadata — see [§14.3](#143-optical-models) |
| `$limit(v, "pnjlim", …)` | yes | installed into the library's `OSDI_LIM_TABLE` at load |
| other `$limit` functions | degrades | identity limiter + a warning; runs unlimited, never crashes |
| `$limexp` | **no** | OpenVAF 23.5 rejects it — a compile error, not a silent no-op |
| `@(timer(…))` | **no** | OpenVAF 23.5 parses only `initial_step`/`final_step`, so a model cannot clock itself |
| `$strobe`, `$display`, `$warning`, `$error`, … | routed, except on macOS/aarch64 | a model's own messages come out through the same switch as fairchild's warnings, so `--quiet` governs them. On macOS/aarch64 they are **removed from the source before compiling** — see below |
| `$finish` | no | |

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

### 14.1.1 A model's diagnostics, and one platform where they are removed

A model talking to the simulator — `$strobe`, `$display`, `$warning`, `$error`,
`$fatal` — arrives through the library's `osdi_log` hook and is printed like any
other fairchild warning, so `--quiet` silences it. Its parameter-range
complaints arrive separately, through the OSDI init-error interface, and are
reported with the parameter they name.

**On macOS/aarch64 those calls are removed from the source before it is
compiled**, and fairchild says so:

```
warning: removed 429 call(s) to $strobe from 'bsim4.va' before compiling: the
Verilog-A compiler miscompiles it on macos-aarch64 and the model would crash the
run instead of printing. The circuit is unchanged, but any condition the model
would have reported — including a parameter it considers out of range — is now
invisible. Set FAIRCHILD_VA_KEEP_UNSUPPORTED=1 to compile the source as written.
```

This is an upstream code-generation bug, not a fairchild limitation: OpenVAF
lowers `snprintf`'s three *fixed* arguments onto the stack, which AArch64 Apple's
variadic convention passes in registers, so the model writes its message to
whatever is in `x0` and the process dies inside `vsnprintf`. On x86-64 the two
conventions coincide, which is why it does not show up on Linux.

It matters because *every* real compact model talks: `bsim4.va` carries 429
`$strobe` calls about its own parameters, and before this transform none of them
could be simulated on an Apple-silicon machine at all — `SIGSEGV`, no output.
The constitutive equations are untouched, so the circuit is the circuit; what is
lost is the model's own commentary, and that is what the warning is for.

`FAIRCHILD_VA_KEEP_UNSUPPORTED=1` compiles the source exactly as written, which
is how to check whether an upstream fix has landed. The transformed source is
kept next to the artefact cache (`--emit-generated <dir>` puts it somewhere you
choose) so you can read what was actually compiled.

### 14.2 Compiling

Either fairchild drives the compile or you do. Both need `openvaf-r` installed;
neither links it in, which is why fairchild stays Apache-2.0 while
OpenVAF-Reloaded is GPL-3.0.

**Fairchild drives it.** Put `.va` (or Spectre `ahdl_include`) in the deck and
run normally. Relevant flags:

| flag | effect |
|---|---|
| `--openvaf <path>` | the compiler to use. Default: `openvaf-r`, then `openvaf`, from `PATH` |
| `--va-include <dir>` | an OpenVAF `-I` search directory. Repeatable, order preserved. Each source's own directory is always searched last |
| `--no-va-compile` | never compile: any `.va` in the deck becomes an error |

`FAIRCHILD_OPENVAF` and `FAIRCHILD_VA_CACHE` set the compiler and the cache
directory for a caller with no command line — the Python binding, mainly.
The cache defaults to `$XDG_CACHE_HOME/fairchild/va`, else `~/.cache/fairchild/va`.

A cached model is reused only when the compiler would produce the same thing:
the key is a hash of OpenVAF's own fully preprocessed source (its
`--print-expansion`) plus the compiler version, so editing any file in the
`` `include `` closure — not just the top one — recompiles. A missing compiler
is an error naming what was looked for, never a skipped device.

**You drive it.** The explicit route, for CI and for anyone without a
toolchain on the machine that runs the simulation:

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

fairchild carries an optical signal on ordinary real-valued MNA unknowns — two
per channel, plus a wavelength label that is resolved before the solve rather
than solved for — so a custom optical discipline needs no compiler support at
all; OSDI passes it through as metadata.

| wire | carries | units |
|---|---|---|
| `<port>_re` | real part of the field envelope | sqrt(W) |
| `<port>_im` | imaginary part | sqrt(W) |
| `<port>_wl` | carrier wavelength — a label, not an unknown | m |

Optical power is `re² + im²`. These are exactly the wires a native
`.optical_port p` expands to (`p_re_0`, `p_im_0`, `p_wl_0`), so a Verilog-A
model drops straight into a link built from native `fc_*` devices — address the
underlying wire names on the `X` line. `examples/verilog_a/models/optical.vams`
is the discipline; `wg_compare.sp` pins a Verilog-A waveguide against the
native one in a single deck.

Two rules, both of which cost real debugging time to learn:

**λ is not a solver unknown, so never contribute to a λ wire.** The wire is
declared — the port list is positional and a bundle is `wpc·N` wires wide either
way — but it carries no equation: a wavelength is routed from its source, never
computed, so fairchild resolves every λ net before the solve. Take the
wavelength from a parameter and sweep it with `--param` or `set_param`. Writing
`OWL(...) <+ …` in a fixed-port model contributes into a node the matrix does
not have, which is a model that compiles, runs, and propagates nothing. (Before
λ was resolved, the objection was the derivative: propagation phase is
thousands of radians — 400 µm at n_g = 4.2 is about 6800 rad — so OpenVAF
differentiating it against a λ unknown put ∂φ/∂λ = φ/λ ≈ 1e9 per metre into the
Jacobian and Newton failed to converge at some wavelengths and not others.)

**`alpha_dB_cm` is power dB, so the amplitude factor divides by 20** —
`10^(−dB/20)`. Everything under `legacy/` predates commit `0f689cb` and is a
factor of two out.

**A Verilog-A model can carry a WDM bundle**, as long as it declares the ports
for it. Dispatch is decided by terminal count: write out the `3·N` wires a
channel bundle expands to (`5·N` under `enable_bidirectional`), and an
`.optical_port bus N` connects to it as one instance serving the whole bus. The
order is positional and per channel — `[re, im, λ]` for channel 0, then channel
1, and so on, input bundle before output bundle — so a wrong guess is a wrong
answer rather than an error. `crates/fairchild-osdi/tests/models/wg_wdm2.va` is
a worked two-channel example.

### One source, any channel count

Writing out `3·N` ports by hand serves exactly one N. To write a model once and
run it at any width, declare the port as a **bundle** and let fairchild generate
the ports for whatever the deck asks for:

| construct | meaning |
|---|---|
| `optical_bundle p, q;` | ports whose width comes from the deck |
| `N(p)` | that width — usable as a loop bound |
| `E_RE(p,k)` / `E_IM(p,k)` | channel `k`'s forward field wires |
| `E_RE_BW(p,k)` / `E_IM_BW(p,k)` | its backward pair — only under `enable_bidirectional` |
| `LAMBDA(p,k)` | channel `k`'s wavelength, filled from the deck's sources |

```verilog
module wg_bundle(a, b);
    optical_bundle a, b;
    parameter real l_um = 1000.0;
    integer k;
    analog begin
        for (k = 0; k < N(a); k = k + 1) begin
            ph = 2.0 * `M_PI * n_g * l_um * 1e-6 / LAMBDA(a, k);
            OF(E_RE(b, k)) <+ amp * (OF(E_RE(a, k)) * cos(ph) + OF(E_IM(a, k)) * sin(ph));
        end
    end
endmodule
```

The channel count appears nowhere. `.optical_port bus 8` in the deck is the only
place a width is written, and the model is compiled for that width and cached
like any other source. Change the deck to 100 channels and the model does not
change.

**There is nothing to write about λ.** `LAMBDA(p, k)` becomes a generated
`parameter real wl_k`, so the physics reads a wavelength with no derivative
attached — and fairchild fills it from whatever source the deck routes to that
port, which it also knows because a bundle model declares its λ terminals and
its slot-for-slot routing. Nothing has to be passed along by hand.

There used to be a second accessor, `WL(p, k)`, for the λ *wire*, and the author
wrote `OWL(WL(out,k)) <+ OWL(WL(in,k));` to carry the tag through the matrix.
λ is no longer in the matrix; a source still using `WL` is refused by name, with
the fix in the message.

`wl_0=1546.12n wl_1=1546.92n` on the instance still works, and applies to a
channel no source reaches. Where a source does reach it and disagrees, the
resolved wavelength wins and says so — two answers for one wavelength is a deck
bug, not a preference.

The loop form is fixed at `for (k = 0; k < N(p); k = k + 1) begin … end`.
Anything else is refused by name rather than expanded into something plausible,
because a mis-generated model is a silently wrong device. `--emit-generated DIR`
writes the expanded source so you can read exactly what was compiled — the first
thing worth looking at when a model behaves at one channel and not at eight.

Still native-only: `DelayLine` group delay and `PhotonicActiveModel`
composition. See [Photonic models](photonic-models.md).

Bidirectional is a middle case, and worth stating precisely rather than leaving
to be discovered. Under `enable_bidirectional` a bundle port expands to `5·N`
wires and `E_RE_BW` / `E_IM_BW` resolve to the backward pair, so a generated
model *can* address light travelling the other way. What does not exist is any
end-to-end deck exercising one: the expansion is unit-tested, the physics is
not. Treat it as unproven rather than unsupported, and anchor a new
bidirectional model on a hand-computed budget rather than on a native device.

### 14.4 Instantiating

```spice
.va   models/va_diode.va                   ; relative to the file that names it
Xd1  a out  va_diode  Is=1e-14 Rs=0.5      ; model name == module name
```

Either directive registers every descriptor in the library under its module
name. From there it resolves like any other model: `X` takes an arbitrary
terminal list, and `D`, `M`, `Q` work for two-, four- and three-terminal models
respectively. All four carry instance parameters.

A relative path resolves against the file the directive is *written in*, not
the top-level deck — so a PDK library in a subdirectory can name its own
sources as siblings, and `.include`-ing it works from anywhere. `.va` sources
load before `.osdi` artefacts, both in deck order.

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

`--param ELEMENT.PARAM=VALUE` overrides an instance parameter on the element
line, including `M`, `Q` and `D`, and takes the same engineering suffixes as a
netlist (`--param "M1.W=1u"`). It does **not** reach `.model` card parameters —
for those, edit the card, or use a `.param` referenced as `{name}`.

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
