# SPICE surface — what a netlist may contain

*Audited against the source and against ngspice 46 on 2026-08-04, by running
every case rather than reading the parser. If this table and the simulator
disagree, the simulator is the bug.*

`docs/model_status.md` covers the **device** dimension: which model-card and
instance parameters are parsed, stamped and validated. This document covers the
**syntax** dimension: which element letters, dot-commands and source functions a
netlist may use at all — in either the SPICE spelling or the Spectre one (§5) —
and, the part that matters, **how fairchild fails when it cannot honour
something.**

Three failure modes, and only one of them is dangerous:

| | meaning |
|---|---|
| **error** | Parse or setup fails, non-zero exit, named directive/letter and line number. Safe: you cannot get a wrong number out of it. |
| **warn** | Runs, prints to **stderr**. Safe only if you are reading stderr — note the CSV goes to stdout, and `--quiet` silences every one of them. |
| **SILENT** | Runs, says nothing, and the answer may be wrong. **These are bugs.** |

The good news up front: every unimplemented element letter and every
unimplemented dot-command is a hard **error**. The silent set is recorded in §4 —
**every item in it is now fixed** — and the section is kept as the record of what
they were, because each one is a shape of bug worth recognising again.

---

## 1. Element letters

The full ngspice set, from `devhelp` on the installed build. "fail" is what
fairchild does with a syntactically plausible line using that letter.

| | ngspice device | fairchild | fail mode |
|---|---|---|---|
| `A` | XSPICE code model | ❌ | error |
| `B` | behavioural (arbitrary) source | ✅ `V=`/`I=` expression (unbraced only) | — |
| `C` | capacitor | ✅ + ESR/ESL parasitics | — |
| `D` | diode | ✅ | — |
| `E` | VCVS — linear voltage-controlled voltage source | ✅ | — |
| `F` | CCCS — current-controlled current source | ✅ | — |
| `G` | VCCS — voltage-controlled current source | ✅ | — |
| `H` | CCVS — current-controlled voltage source | ✅ | — |
| `I` | independent current source | ✅ | — |
| `J` | JFET | ❌ | error |
| `K` | mutual inductance | ✅ | — |
| `L` | inductor | ✅ + ESR parasitic | — |
| `M` | MOSFET | ✅ Level 1 only — `LEVEL≠1` warns loudly | — |
| `N` | numerical device (NUMD, NBJT) | ❌ | error |
| `O` | lossy transmission line (LTRA) | ❌ | error |
| `P` | coupled multiconductor line | ❌ | error |
| `Q` | BJT | ✅ Gummel-Poon L1 | — |
| `R` | resistor | ✅ + ESR/ESL parasitics | — |
| `S` | voltage-controlled switch | ✅ | — |
| `T` | lossless transmission line | ✅ | — |
| `U` | uniform RC line (URC) | ❌ | error |
| `V` | independent voltage source | ✅ incl. `AC <mag> [phase]` | — |
| `W` | current-controlled switch | ✅ | — |
| `X` | subcircuit instance | ✅ + OSDI/Verilog-A instances | — |
| `Y` | simple lossy line (TransLine/txl) | ❌ | error |
| `Z` | MESFET / HFET | ❌ | error |

**`E`/`F`/`G`/`H` are supported in their linear form only:**

```spice
E<n> p n nc+ nc- <gain>     V = gain·(V(nc+) − V(nc-))
G<n> p n nc+ nc- <gain>     I = gain·(V(nc+) − V(nc-))
H<n> p n <Vctrl>  <gain>    V = gain·I(Vctrl)
F<n> p n <Vctrl>  <gain>    I = gain·I(Vctrl)
```

All four are desugared onto the B-element rather than given their own stamps — a
VCVS *is* `B… V=gain*(V(cp)-V(cn))` — so they inherit its auxiliary branch row,
its per-reference Jacobian columns, and its tested sign conventions. Values are
pinned against ngspice in `crates/fairchild-core/tests/controlled_sources.rs`,
including the one that is easy to get backwards: a `G` pushes current *out of*
`n+`, so a VCCS wired across its own output nodes with `gm = 1/R` is exactly that
resistor.

`POLY(n)`, `VALUE={…}` and `TABLE` — the other SPICE spellings of these four —
are refused by name rather than mis-read as a node. Write a polynomial or
expression source as a `B` element instead.

`R`/`C` referencing a `.model` card — ngspice's semiconductor resistor
(`R1 n1 n2 rmod L=… W=…`) — is an error, not a silent mis-parse.

## 2. Dot-commands

| | ngspice meaning | fairchild | fail mode |
|---|---|---|---|
| `.title` | title line | ✅ (first line, implicit) | — |
| `.end` | end of deck | ✅ | — |
| `.subckt` / `.ends` | subcircuit definition | ✅ nested instantiation, `{}` arithmetic | — |
| `.include` | include a file | ✅ 16-deep | — |
| `.lib` / `.endl` | library section | ✅ | — |
| `.model` | model card | ✅ | — |
| `.param` | parameter definition | ✅ values may be expressions: `{…}`, `'…'`, or bare | §4.8 |
| `.func` | user function | ✅ expanded at parse time | §4.8 |
| `.global` | global nets | ✅ (a port of the same name is refused) | §4.10 |
| `.csparam` | constant → control variable | ❌ | error |
| `.table` | lookup table | ❌ | error |
| `.if` / `.elseif` / `.else` / `.endif` | netlist conditionals | ✅ resolved at parse time (not inside `.subckt`) | §4.8 |
| `.control` / `.endc` | interactive control block | ⚠️ skipped, **warns once** | §4.7 |
| `.op` | operating point | ✅ | §4.4 |
| `.dc` | DC sweep | ✅ nested, parallel | §4.4 |
| `.ac` | AC small-signal | ✅ magnitude and phase honoured | §4.4 |
| `.tran` | transient | ✅ incl. `tstart`, `tmax`, `UIC` | §4.4 |
| `.noise` | noise analysis | ✅ | §4.4 |
| `.disto` | small-signal distortion | ❌ | error |
| `.pz` | pole-zero | ❌ | error |
| `.sens` | sensitivity | ❌ | error |
| `.tf` | transfer function | ❌ | error |
| `.four` / `.fourier` | Fourier analysis | ❌ | error |
| `.sp` | S-parameter | ❌ | error (verified) |
| `.ic` | initial conditions | ✅ (honoured with `UIC`) | — |
| `.nodeset` | DC solution hint | ✅ | — |
| `.save` | select saved vectors | ⚠️ ignored, **warns once** | §4.6 |
| `.print` | tabular output | ⚠️ ignored, **warns once** | §4.6 |
| `.plot` | line-printer plot | ⚠️ ignored, **warns once** | §4.6 |
| `.probe` | select probes | ⚠️ ignored, **warns once** | §4.6 |
| `.width` | output width | ⚠️ ignored, **warns once** | §4.6 |
| `.measure` / `.meas` | measurements | ✅ (tran; see model_status §10) | — |
| `.options` / `.option` | simulator options | ✅ known keys; unknown keys **warn** | — |
| `.temp` | temperature (incl. sweep) | ✅ | — |
| `.backanno` | (LTspice) back-annotation | ⚠️ silently ignored | §4.6 |

`.control` blocks deserve a note: a great many ngspice decks in the wild wrap the
whole run in `.control … .endc`. The block is now skipped with a warning rather
than refusing the file — see §4.7 for what that does and does not buy you. It is
**not** interpreted, and will not be.

**fairchild extensions** (not ngspice): `.osdi`, `.optical`, `.optical_port`,
`.optical_bus`, `.electrical_port`, and `.alter` (which is HSPICE's).

## 3. Independent-source functions

| | fairchild | fail mode |
|---|---|---|
| `DC <v>` | ✅ | — |
| `PULSE(…)` | ✅ | — |
| `SIN(…)` | ✅ | — |
| `EXP(…)` | ✅ | — |
| `PWL(…)` | ✅ | — |
| `SFFM(…)` | ✅ | — |
| `AM(…)` | ✅ | — |
| `TRNOISE(…)` | ❌ | error |
| `TRRANDOM(…)` | ❌ | error |
| `AC <mag> [phase]` | ✅ | — |

Time-domain noise exists, but as `.options trannoise=1` over the `.noise` source
list rather than as a `TRNOISE` source function.

---

## 4. The silent set

Kept as the record. Each was found by running a deck, not by reading code.
§4.1–§4.5 were the original five; §4.9 and the first half of §4.8 are two more
found later, in the same way. §4.6 and §4.7 were never silent — `.save` and
`.control` were loud errors — but they are recorded here because the same audit
turned them up and one rule decided them: who owns the decision, the deck or the
caller.

### 4.1 `AC <mag> [phase]` on a source line was discarded — FIXED

```spice
V1 in 0 DC 0 AC 2      ← the "AC 2" is parsed away and never reaches the solver
```

`.ac` always drives with unit amplitude and zero phase. Measured on an RC
divider at 1 kHz:

| | `AC 1` | `AC 2` | `AC 5` | no AC spec |
|---|---|---|---|---|
| ngspice | 0.99996 | 1.99992 | 4.99980 | — |
| fairchild | 0.99998 | 0.99998 | 0.99998 | 0.99998 |

So an ngspice deck using `AC 2` gives results **2× too small, with no
diagnostic**. This is the worst item in this document: supported analysis,
supported syntax, wrong number, silence. Phase is likewise unavailable, so
multi-source AC (anything needing a 90° drive) cannot be expressed at all.

**Fixed.** `AcSpec { mag, phase_deg }` on the source element, split off the line
by `split_ac_spec` (the spec may sit anywhere after the nodes), and used by
`build_ac_rhs`, which now returns both a real and an imaginary RHS. Verified
against ngspice: magnitude scales exactly, and `AC 1 90` on an RC at its corner
gives +45° where `AC 1 0` gives −45°.

**SPICE semantics, and only those.** A source without an `AC` spec is not an AC
source and contributes nothing; a deck with no AC source at all is a hard
`no AC source` error rather than a quiet zero.

There was briefly a compatibility fallback here — a deck declaring no spec
anywhere kept the old unit drive — and it was removed on purpose. It made the
rule depend on the deck's contents, and it preserved the worst part of the
original bug: with no spec anywhere, *every* source in the circuit is still
driven at unit amplitude, so DC bias rails get excited as though they were
signal generators. In a multi-rail circuit that is wrong in a way no single
number reveals. Requiring the spec costs one token per deck and removes the
whole class.

### 4.2 `.model x D(IS=…)` with no space before `(` lost the first parameter — FIXED

```spice
.model xx D (IS=1e-16)   → I(V1) = 5.670e-5 A     correct
.model xx D(IS=1e-16)    → I(V1) = 5.670e-3 A     IS silently defaulted to 1e-14
```

Two orders of magnitude apart. The kind token becomes `d(is=1e-16`, which still
`starts_with('d')`, so it dispatches as a diode and the parameter vanishes.
MOSFET and BJT cards escape it — their dispatch is exact, so they raise
`unknown model` — which makes the diode the one device that fails quietly.
ngspice accepts the no-space form.

**Fixed.** `parse_model` splits the kind token on `(` before lowercasing it, and
folds the glued remainder back into the parameter list. Both spellings now give
−5.670295e-5 A, against ngspice's −5.67035e-05.

### 4.3 Unknown `.options` keys were accepted in silence — FIXED

```spice
.options reltol=1e-4 banana=7 trtol=7 chgtol=1e-14   ← all four accepted, no output
```

`trtol` and `chgtol` are real ngspice options someone would reasonably set and
expect to matter. `banana` shows there is no validation at all.

**Fixed.** `SimOptions::set` already returned `false` for an unrecognised key —
`from_netlist` was discarding it. It now warns. Known keys stay quiet.

### 4.4 `.tran`'s third and fourth arguments were discarded, then half-applied — FIXED

```spice
.tran 1n 20n 10n        tstart — ngspice suppresses output before 10 ns
.tran 1n 20n 0 0.1n     tmax   — ngspice caps the step at 0.1 ns
```

Measured originally: with `tstart=10n`, ngspice's first output row is at 1e-8;
fairchild's was at 0 and the row count unchanged. With `tmax=0.1n`, fairchild
returned the same 21 rows as without it — the step cap ignored, so a deck that
asked for finer resolution silently did not get it.

(`UIC` *is* honoured — verified, `.ic v(out)=0.9` survives to t=0 with `UIC` and
is discarded without it.)

**Fixed, twice.** Both are parsed (skipping `UIC`, which may occupy any trailing
slot). The first fix picked them up in `SimOptions::from_netlist`, which every
frontend calls whether or not it is running the deck's analyses — so
`Circuit.run("tran", step=…, stop=…)`, which takes its timing entirely from the
caller, still inherited the card's `tmax`. Half the card applied and nobody chose
that. A deck with two `.tran` lines was worse: both runs got the tightest `tmax`
of the two.

An analysis card is now honoured **as a unit**, by whoever is about to run it
(`SimOptions::apply_tran_card`):

- The CLI runs every card in deck order, each with its own `tstart`, `tmax` and
  `UIC` — unchanged behaviour for a single-card deck, correct now for several.
- `Circuit.run("tran")` with no timing kwargs adopts the deck's card whole:
  `step`, `stop`, `tstart`, `tmax`, `UIC`. Pass any of `step`/`stop` and the card
  is not consulted at all. The same rule covers `.ac`, `.dc` (through
  `run("dc_sweep")`) and `.noise`; `Circuit.analyses` lists what a deck declares.
- With neither a card nor kwargs, `run("tran")` is an error naming both fixes,
  and a deck declaring two cards of one kind is an error rather than a silent
  pick of the first.

Verified: `tstart=10n` starts at 1e-8 like ngspice; `tmax=0.1n` turns 21 rows
into 202 when the card is adopted, and leaves a caller-timed run at 21.

**Behaviour change worth noting:** a deck's `tstart`/`tmax` no longer reach a
Python or C-API run that supplied its own `step` and `stop`. A deck that was
relying on the leak gets a coarser `max_step` than before — pass `maxstep=` (or
no timing at all, and take the card) to get it back.

### 4.5 `LEVEL=` on a MOSFET card was only warned about generically — FIXED

```spice
.model nm NMOS (LEVEL=3 VTO=0.7 KP=100u KAPPA=0.2 THETA=0.1)
warning: MOSFET model 'nm' params not yet implemented (using defaults): level, kappa, theta
```

Not silent, so not a bug — but the message treats `LEVEL` as one defaulted
parameter among several, when what actually happened is that a Level 3 model was
simulated as Level 1. That is a different device, not a defaulted coefficient,
and it deserves to say so.

**Fixed.** `LEVEL != 1` gets its own warning saying the card is being simulated
as Level 1 and that currents and capacitances will differ from the intended
model — not merely in the unset parameters. `LEVEL=1` stays quiet.

### 4.6 Output-selection directives were ignored in silence, or not at all — FIXED

```spice
.print tran V(out)      ← ignored, silently; you get every node
.probe V(out)           ← ignored, silently
.save V(out)            ← hard error
.width out=80           ← hard error
```

Output selection is `--probe` (CLI) or indexing the returned result (Python), so
a deck's version has nowhere to go — every signal is available either way. That
part is by design. Two things were not: it happened in silence, so a deck whose
`.print tran V(out)` you expect to narrow the output gets every node instead;
and `.save`/`.width`, which are the *same class of directive*, failed a
different way — refusing to load a deck the others accept. Nobody chose that
split either.

**Fixed.** All five load and warn, once per directive however many lines of it a
deck carries:

```
warning: .print is ignored — output selection belongs to the frontend, not the
deck: use --probe (CLI) or index the returned result (Python). Every signal is
available either way (see docs/spice_support.md §4.6)
```

`.backanno` (LTspice) stays silent: it selects nothing and there is no fairchild
mechanism to point at.

This is the Select class of the directive rule — see
[who owns the run](user-guide.md#who-owns-the-run--the-deck-or-the-caller).

### 4.7 `.control` refused to load the deck — FIXED, and it will never be interpreted

```spice
.control
run
let vpk = maximum(v(out))
write rc.raw v(out)
plot v(out)
.endc
```

A great many ngspice decks in the wild wrap the whole run in a block like this.
Refusing the file over it blocked loading a large share of third-party decks for
a construct that, in the large majority of real blocks, contains only
`run`/`write`/`plot` — all three of which the frontend already does.

**Fixed by skipping, not by interpreting.** The block is consumed in pass 1,
never reaches the parser proper, and the commands found in it are named once on
stderr:

```
warning: .control block skipped — its commands are not interpreted (run, let,
write, plot). fairchild is not an ngspice shell: control flow belongs in Python
(fairchild.Circuit) or in CLI flags, and output selection is --probe. An analysis
that existed only inside the block will not run: give the deck a
.tran/.ac/.dc/.op card, or drive the run from Python
(see docs/spice_support.md §4.7)
```

A `.control` with no `.endc`, and an `.endc` with no `.control`, are both hard
errors. Guessing where an unterminated block ended would discard the rest of the
deck silently, which is how a circuit loses half its elements and still runs.

**What this does not buy you.** If the deck's only analysis was a `tran`/`ac`
command *inside* the block, nothing runs — you get this warning and then "no
analyses found in netlist". Add a real `.tran` card, or drive the analysis from
Python. Nothing is silently substituted.

**Why it will not be interpreted.** `.control` is imperative script: `run`,
`let`, `write`, `alter`, loops, conditionals. The Python bindings exist so that
control flow lives in Python. A second scripting language inside the simulator
would be the largest single item on this list *and* would compete with the
feature it duplicates. This is the Script class of the directive rule — the one
class fairchild declines rather than honours; see
[who owns the run](user-guide.md#who-owns-the-run--the-deck-or-the-caller).

### 4.8 A `.model` value written as a parameter was silently dropped — FIXED

```spice
.param vt=0.7
.model nm NMOS (VTO={vt} KP=100u)     ← VTO defaulted; nothing said so
```

`.model` lines were the one place parameter substitution did not run. `{vt}` was
not evaluated, so `vto` landed in the card's *expression* params, which the
MOSFET path never reads — the card was simulated with the default threshold. A
different transistor than the deck asked for, with no warning. Subcircuit-local
cards were already substituted at instantiation, which is what made the top-level
omission look deliberate.

**Fixed**, together with the two gaps that made it hard to hit and easy to
mis-diagnose:

- **`.model` values are substituted** like any element value. `{…}` and `'…'` are
  parse-time expressions everywhere in a deck now, including on a card.
  `"…"` on a card still means a device constitutive map over the device's own bias
  (`dneff="5.0e-5*V"`) and is left untouched — that is now the only spelling for
  one.
- **`.param` values may be expressions.** `.param b={a*3}` was `invalid number
  '{a*3}'`. Values may be braced, single-quoted (the HSPICE spelling), or bare,
  and resolve in file order over the parameters already defined — including
  earlier on the same line. A value may contain spaces only inside braces or
  quotes, which is HSPICE's rule for HSPICE's reason: `.param a = 1 + 2 b = 3` has
  no unambiguous reading. `1k` stays a number.
- **`.func name(args) = body`** is expanded at parse time, over the syntax tree
  rather than the source text — `f(x)=x+1` called as `2*f(3)` is 8, and textual
  expansion gets that wrong without parenthesising every substitution. It works
  anywhere an expression does: `{…}`, `.param`, a `.model` value, a B-source,
  `.measure`. Definitions may follow their first use. Refused: a name that shadows
  a built-in, a repeated formal, recursion (there is no finite expansion),
  and the wrong argument count.

Also fixed in passing: an **undefined parameter is now named**. `{2*nope}` said
only that the result was not finite.

**`.if` / `.elseif` / `.else` / `.endif`** land here too, as netlist
preprocessing: the condition is an ordinary parse-time expression over `.param`
values and `.func` calls, and only the taken branch is collected — a `.model`,
`.subckt` or `.param` in a dead branch does not exist afterwards. A condition in a
branch that cannot run is not evaluated, so a dead branch may reference names that
were never defined.

Two refusals worth stating:

- **`.if` inside a `.subckt` is an error.** Its condition would be evaluated once,
  against the subcircuit's *default* parameters, and then apply to every instance
  — a wrong answer per instance rather than a missing feature.
- **A condition over an undefined name is an error**, not `false`. `nope == 1`
  compares NaN against 1 and yields a perfectly finite `false`, so a misspelled
  corner variable would have silently selected the other branch. Names are checked
  before the expression is evaluated.

### 4.9 An unknown function read as zero — FIXED

```spice
.param a=2
R1 in 0 {frobnicate(a)}     ← a 0 Ω resistor, silently
```

The expression evaluator resolved an unknown function name — or a known one with
the wrong argument count — to `0.0`. Every other undefined thing in that grammar
resolves to NaN precisely so the caller's finiteness check can refuse it; calls
were the hole.

**Fixed.** Unresolvable calls evaluate to NaN, and the parser refuses the
expression by name (`unknown function 'frobnicate'`) rather than reporting that
something in it was not finite.

### 4.10 A controlled source inside a subcircuit read zero — FIXED

```spice
.subckt amp inp outp
R1 inp mid 1k
R2 mid 0 1k
B1 outp 0 V=v(mid)*2
.ends
X1 a y amp        ← V(y) was 0.0 V; it is 2 × 0.5 = 1.0 V
```

Subcircuit flattening renames nodes and elements. It renamed the element's own
terminals and left the references *inside* its expression alone, so `v(mid)` still
named a top-level `mid` that does not exist — and an unknown node reads as zero.
`E`, `F`, `G` and `H` desugar onto the B-element, so all four were silently dead
inside any subcircuit, in DC, AC and transient alike. A current-controlled pair
missed twice over: `F1 … Vsense …` references a source name that flattening had
prefixed.

**Fixed.** `Expr::rename_refs` applies the same two maps flattening already uses —
the node map (port → call-site net, ground stays ground, `.global` passes through,
everything else namespaced) and the element prefix for a branch reference — at the
point that was already remapping the terminals.

**`.global` lands here**, being the other half of the same resolver:

```spice
.global vdd vss        ← the same node in every scope, port list or not
```

Declarations are collected before any instance is expanded, so a `.global` may
follow the instance that needs it. Nesting is unlimited: a supply reaches a
subcircuit two levels down without appearing in either port list, which is what
CDL and foundry decks written by a layout tool expect.

A net that is **both a port and global is refused**. The port would take the
caller's net while every reference inside took the global one, and picking either
silently is wrong for the deck that meant the other. `.global 0` warns instead —
ground is already global in every scope, so the declaration is redundant rather
than wrong.

---

### 4.11 A subcircuit parameter was resolved once, for the default — FIXED

```spice
.subckt rdiv a b n=1
.param rtot={1000*n}      ← evaluated when the DEFINITION was read, with n=1
R1 a b {rtot}
.ends
X1 in 0 rdiv n=2          ← 1 kΩ, not 2 kΩ. No warning, no error, wrong current.
```

A `.subckt` is the only construct whose parameters have more than one answer: one
per instance. Both were collected as **numbers** when the definition was read —
header defaults through `parse_value`, body `.param` lines evaluated against the
defaults then and there — so every expression froze at the default. An instance
that overrode `n` got its own `n` and everybody else's `rtot`.

The same freeze made the honest half of the seam impossible: a header default
*written* as an expression (`.subckt r a b w=1u rsh='100/w'`) could not be a number
at collection time, so it was a hard error — which is what a foundry wrapper looks
like from the first line.

**Fixed.** A definition now keeps its parameters as **source text**, and each
instance resolves them in order: enclosing scope, then header defaults with the
call's overrides in place *before* anything reads them, then the body's `.param`
assignments. Two instances of one definition resolve independently.

Two refusals came with it, both at the same seam:

- **An undeclared instance parameter is an error.** `X1 p q rdiv nn=2` used to add
  `nn` to the scope, leave `n` at its default, and report a clean answer for a
  circuit nobody described. The error names the parameter and lists what the
  definition declares.
- **Overriding a body `.param` is an error**, not a silent choice between
  overriding and recomputing. The message says which line computes it and that
  moving it to the header makes it overridable.

```spice
X1 p q rdiv n=2*3         ← also now an error: an unreadable instance-parameter
                            value used to be dropped, leaving the default in place
```

---

## 5. The Spectre dialect

Foundry model libraries are written in Spectre (`.scs`), so fairchild reads that
spelling too. It is a **front end, not a second parser**: every statement is
transliterated into the equivalent SPICE statement and the existing passes do the
rest, so exactly one place still decides what a resistor is and how a subcircuit
flattens. Transliteration is line-aligned, so an error still names the line you
wrote.

The dialect is detected **per file, from the content** — `simulator lang=` or a
leading `//` comment — not from the extension. A SPICE deck may therefore
`.include` a Spectre library and vice versa, and neither caller chooses a mode.

| Spectre | becomes |
|---|---|
| `parameters a=1 b=2*a` | `.param a=1 b={2*a}` |
| `R1 (in out) resistor r=1k` | `R1 in out 1k` |
| `C1 (out 0) capacitor c=1n`, `L… inductor l=` | `C1 out 0 1n`, `L… …` |
| `V1 (in 0) vsource dc=1.8`, `isource` | `V1 in 0 DC 1.8`, `I…` |
| `X1 (a b) mycell w=1u` | `X1 a b mycell w=1u` |
| `X1 a b mycell` (bare nodes) | same — both spellings are read |
| `dc1 dc`, `tr1 tran stop=1n`, `ac1 ac`, `n1 noise` | `.op`, `.tran … 1n`, `.ac …`, `.noise …` |
| `include "f"` / `include "f" section=s` | `.include "f"` / `.lib "f" s` |
| `vdd!`, `global vdd!` | a net, plus `.global vdd!` |
| `opt1 options temp=85` | `.temp 85` |
| a trailing `\\` or leading `+` continuation, `//` and column-1 `*` comments | joined / stripped |

Failure modes, same three as above:

| construct | mode |
|---|---|
| `subckt` / `inline subckt` / `ends`, `model`, `if`/`else`, `function` | **error** naming the construct — not read yet, so a flat deck is the current surface |
| Any other unreadable statement | **error** with the line and the two forms it accepts |
| `ahdl_include "m.va"` | **warn**, skipped — Verilog-A needs an out-of-band OpenVAF compile and a `.osdi` line, so the devices it defines would otherwise be unknown models |
| `save`, `assert`, `statistics`, `montecarlo`, `sweep`, `alter`, `altergroup`, `check`, `info`, `shell` | **warn**, skipped — they do not change a single solve |
| `options` keys other than `temp` | **warn** (as `.options` above) |

Nothing here is SILENT: a statement is either transliterated, reported as
skipped, or refused.

## 6. What is left

§4 is done, and `E`/`F`/`G`/`H` are in. Remaining, by (value × cheapness):

1. **Spectre `subckt` / `model` / `if`** — a foundry library is nothing but these,
   so the flat surface in §5 does not yet read one. Each has an exact SPICE
   equivalent (`.subckt`, `.model`, `.if`), which is what the front end will emit.
2. **Accept `{…}` on a `B` line** — let the B-element claim its own braces before
   `.param` substitution runs. Small, and it is the difference between an
   ngspice behavioural deck loading and not.
3. **Model binning and `level=` routing** — what actually stands between this and
   a foundry deck; neither is implemented, and picking a wrong W/L bin would be a
   silent wrong answer, so geometry outside every bin must be a hard error.
4. Everything else in §1–§2 is a clean error and can wait for a use case.

## How to update this document

Re-run the probes rather than reading the parser — every "SILENT" above was
found by running a deck and comparing against ngspice, and two of them
contradicted what the code appeared to do.
