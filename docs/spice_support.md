# SPICE surface — what a netlist may contain

*Audited against the source and against ngspice 46 on 2026-08-04, by running
every case rather than reading the parser. If this table and the simulator
disagree, the simulator is the bug.*

`docs/model_status.md` covers the **device** dimension: which model-card and
instance parameters are parsed, stamped and validated. This document covers the
**syntax** dimension: which element letters, dot-commands and source functions a
netlist may use at all, and — the part that matters — **how fairchild fails when
it cannot honour something.**

Three failure modes, and only one of them is dangerous:

| | meaning |
|---|---|
| **error** | Parse or setup fails, non-zero exit, named directive/letter and line number. Safe: you cannot get a wrong number out of it. |
| **warn** | Runs, prints to **stderr**. Safe only if you are reading stderr — note the CSV goes to stdout. |
| **SILENT** | Runs, says nothing, and the answer may be wrong. **These are bugs.** |

The good news up front: every unimplemented element letter and every
unimplemented dot-command is a hard **error**. The silent set was small, and §4
records it — **all five silent items are now fixed**; the section is kept as the
record of what they were, because each one is a shape of bug worth recognising
again.

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
| `E` | **VCVS** — linear voltage-controlled voltage source | ❌ | error |
| `F` | **CCCS** — current-controlled current source | ❌ | error |
| `G` | **VCCS** — voltage-controlled current source | ❌ | error |
| `H` | **CCVS** — current-controlled voltage source | ❌ | error |
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

**`E`/`F`/`G`/`H` are the notable gap.** Linear controlled sources are core SPICE
and appear in most real op-amp, feedback and macromodel decks. `B` reaches the
same physics — `B1 out 0 V=2*V(in)` gives 3.0 V from a 1.5 V input, verified —
but a deck written for any other simulator will not load.

Note the brace form does **not** work, though ngspice accepts it:

```spice
B1 out 0 V=2*V(in)      ✅ 3.0 V
B1 out 0 V={2*V(in)}    ❌ error — ngspice gives 3.0 V
```

`{…}` is claimed by `.param` substitution before the B-element sees it, so the
expression is evaluated as a parameter and fails on the node reference. The error
message says so clearly, so this is a compatibility gap rather than a trap —
but it is a one-character difference away from a deck that runs.

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
| `.param` | parameter definition | ✅ (`{expr}` references) | — |
| `.func` | user function | ❌ | error |
| `.global` | global nets | ❌ | error |
| `.csparam` | constant → control variable | ❌ | error |
| `.table` | lookup table | ❌ | error |
| `.if` / `.elseif` / `.else` / `.endif` | netlist conditionals | ❌ | error |
| `.control` / `.endc` | interactive control block | ❌ | error |
| `.op` | operating point | ✅ | — |
| `.dc` | DC sweep | ✅ nested, parallel | — |
| `.ac` | AC small-signal | ✅ magnitude and phase honoured | — |
| `.tran` | transient | ✅ incl. `tstart`, `tmax`, `UIC` | — |
| `.noise` | noise analysis | ✅ | — |
| `.disto` | small-signal distortion | ❌ | error |
| `.pz` | pole-zero | ❌ | error |
| `.sens` | sensitivity | ❌ | error |
| `.tf` | transfer function | ❌ | error |
| `.four` / `.fourier` | Fourier analysis | ❌ | error |
| `.sp` | S-parameter | ❌ | error (verified) |
| `.ic` | initial conditions | ✅ (honoured with `UIC`) | — |
| `.nodeset` | DC solution hint | ✅ | — |
| `.save` | select saved vectors | ❌ | error |
| `.print` | tabular output | ⚠️ **silently ignored** | §4.6 |
| `.plot` | line-printer plot | ⚠️ **silently ignored** | §4.6 |
| `.probe` | select probes | ⚠️ **silently ignored** | §4.6 |
| `.width` | output width | ❌ | error |
| `.measure` / `.meas` | measurements | ✅ (tran; see model_status §10) | — |
| `.options` / `.option` | simulator options | ✅ known keys; unknown keys **warn** | — |
| `.temp` | temperature (incl. sweep) | ✅ | — |
| `.backanno` | (LTspice) back-annotation | ⚠️ silently ignored | §4.6 |

`.control` blocks deserve a note: a great many ngspice decks in the wild wrap the
whole run in `.control … .endc`, so this single gap blocks loading those decks
outright. It is an error rather than a silent skip, which is correct — but if
loading third-party decks matters, skipping the block with a warning would get
further than refusing the file.

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

## 4. The silent set — all five now fixed

Kept as the record. Each was found by running a deck, not by reading code.

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

One deliberate wrinkle: a deck that declares **no** `AC` spec anywhere keeps the
old unit drive, because there is no stated intent to contradict and the one
in-tree `.ac` example relies on it. As soon as any source declares a spec the
deck is read strictly, ngspice-style — which also fixes a second silent bug
nobody had noticed: pairing an `AC 1` drive with a plain DC bias source used to
excite the bias source too.

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

### 4.4 `.tran`'s third and fourth arguments were discarded — FIXED

```spice
.tran 1n 20n 10n        tstart — ngspice suppresses output before 10 ns
.tran 1n 20n 0 0.1n     tmax   — ngspice caps the step at 0.1 ns
```

Measured: with `tstart=10n`, ngspice's first output row is at 1e-8; fairchild's
is at 0 and the row count is unchanged. With `tmax=0.1n`, fairchild returns the
same 21 rows as without it — the step cap is ignored, so a deck that asked for
finer resolution silently does not get it.

(`UIC` *is* honoured — verified, `.ic v(out)=0.9` survives to t=0 with `UIC` and
is discarded without it.)

**Fixed.** Both are parsed (skipping `UIC`, which may occupy any trailing slot)
and picked up in `SimOptions::from_netlist`, so the CLI and the Python bindings
get them from one place. `tstart` trims the finished result — SPICE's tstart
selects what is *saved*, not where integration begins. `tmax` clamps
`max_step` with `min`, so it can never loosen a ceiling set via
`.options maxstep`. Verified: `tstart=10n` now starts at 1e-8 like ngspice, and
`tmax=0.1n` turns 21 rows into 202.

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

### 4.6 `.print` / `.plot` / `.probe` / `.backanno` are ignored by design — no change

Output selection is `--probe` instead, so these have nowhere to go. Benign, but
undocumented until now: a deck whose `.print tran V(out)` you expect to narrow
the output gets every node instead.

---

## 5. What is left

§4 is done. Remaining, by (value × cheapness):

1. **`E`/`F`/`G`/`H` controlled sources** — the largest *capability* gap, and the
   most likely reason a third-party deck will not load. All four are linear
   stamps; `B` already proves the machinery.
2. **`.control` block skip-with-warning** — unblocks a large fraction of
   real-world ngspice decks without implementing the shell.
3. **Accept `{…}` on a `B` line** — let the B-element claim its own braces before
   `.param` substitution runs. Small, and it is the difference between an
   ngspice behavioural deck loading and not.
4. Everything else in §1–§2 is a clean error and can wait for a use case.

## How to update this document

Re-run the probes rather than reading the parser — every "SILENT" above was
found by running a deck and comparing against ngspice, and two of them
contradicted what the code appeared to do.
