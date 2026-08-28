# Model status — what is parsed, what is stamped, what is validated

*Audited against the source on 2026-08-22. If you find a disagreement between
this table and the simulator, the simulator is the bug — and for the
"accepted, not modelled" rows, a test fails.*

This document exists because "supported" is not a binary. A parameter can be
accepted by the parser, silently discarded before it reaches the matrix, and
never covered by a test — and from the outside all three look identical: the
netlist runs and you get a number.

Three columns, and they mean exactly this:

| Column | Meaning |
|---|---|
| **Parsed** | The parameter is read off the `.model` card or element line without an error. |
| **Stamped** | It actually changes the matrix or the residual. A parameter that is parsed but not stamped **changes nothing**. |
| **Validated** | A test pins it — `ngspice` where the column says so, otherwise an analytic or equivalence test in-tree. |

Legend: ✅ yes · ⚠️ partial (see the note) · ❌ no.

**How a card value may be written.** Any value on a `.model` card may be a
parse-time expression over `.param` values and `.func` calls, braced (`VTO={vt}`)
or single-quoted (`VTO='vt*1.05'`); both are evaluated before the card is built,
so the tables below apply unchanged. Until 2026-08-18 they were *not* evaluated on
a top-level card: the value landed in the card's expression params, which only the
`fc_phase_shifter_expr` and `fc_awgr` kinds read, and every other card defaulted
the parameter in silence. Double quotes still mean a device constitutive map over
the device's own bias (`dneff="5.0e-5*V"`) and are left symbolic — a card whose
kind cannot read one now warns instead of dropping it.

**The short version of what to watch out for:**

- Every parameter is now either stamped or named. A model-card parameter this
  simulator accepts and does not model produces one warning per card saying what
  the deck loses (`IKF ignored: high-injection roll-off is not modelled, so the
  forward current keeps its exponential slope past the knee`), and an instance
  parameter a device cannot honour is named per instance. The lists live in
  `crate::unmodelled`, and `the_unmodelled_tables_match_the_audit_document`
  fails if they and the tables below disagree.
- Rows below marked **⚠️ accepted, not modelled** are exactly those warnings.
  They are the ones to read before trusting a foundry card.
- `.noise` has no ngspice comparison. Optical noise (PD shot, laser RIN) is
  checked against the analytic receiver budget, and `.options trannoise=1`
  injects the same generators into `.tran`. Neither domain models flicker (1/f)
  or RTS noise.
- The photonic models are validated against analytic forms and against
  themselves, never against an external simulator. That is the biggest gap in
  this document — see §9.

---

## 1. Passive elements

| Element | Parsed | Stamped | Validated |
|---|---|---|---|
| `R<name> n+ n- <value>` | ✅ | ✅ | ✅ ngspice (`rc_step`, dividers) |
| `C<name> n+ n- <value>` | ✅ | ✅ | ✅ ngspice transient + AC |
| `L<name> n+ n- <value>` | ✅ | ✅ | ✅ ngspice (`rl_step`) |
| `K<name> L1 L2 <k>` | ✅ | ✅ | ✅ ngspice 1:1 transformer step, 2 % |

### Passive parasitics (instance parameters)

Each expands into an equivalent sub-network of ordinary elements at parse time,
so they are exactly as validated as R/L/C are. **Which parameter is accepted on
which element is not uniform** — an `esl=` on a resistor is silently no-ops,
because the R arm never looks at it:

| Parameter | Accepted on | Expands to | Validated |
|---|---|---|---|
| `CPAR` | **R, L** | a parallel `C` across the element | ⚠️ in-tree only (`passive_parasitics.rs`) |
| `ESR` | **L** (as series R), **C** (as series R) | a series `R` through an internal node | ⚠️ in-tree only |
| `ESL` | **C** | a series `L` through an internal node | ⚠️ in-tree only |
| `RPAR` | **C** | a parallel `R` across the element | ⚠️ in-tree only |

No ngspice comparison exists — ngspice has no equivalent syntax, so the
expansion is checked in-tree against the hand-written equivalent network.

---

## 2. Independent sources

`V` and `I` accept `DC`, `PULSE`, `PWL`, `SIN`, `EXP`, `SFFM`, `AM`. All are
parsed, evaluated, and exercised by the transient goldens.

✅ **`AC <mag> [phase]` is honoured**, magnitude and phase both, and the spec may
sit anywhere after the nodes — before or after a DC value or a transient
function. Verified against ngspice: on an RC at its corner, `AC 1 90` gives
+45.0009° where `AC 1 0` gives −44.9991°.

The rule is ngspice's, and it is strict on purpose: **a source with no `AC` spec
is not an AC source and contributes nothing** to `.ac`. A deck with no `AC` spec
anywhere is `SimError::NoAcSource` rather than a quiet zero. Before 0.3.0 an
unspecified deck drove *every* source at unit amplitude, which excited DC bias
rails as though they were signal generators — wrong in a way no single number
reveals.

---

## 3. Diode (`D` / `.model … D`)

| Parameter | Parsed | Stamped | Validated |
|---|---|---|---|
| `IS` | ✅ | ✅ | ✅ ngspice DC |
| `N` | ✅ | ✅ | ✅ ngspice DC |
| `RS` | ✅ | ✅ | ✅ ngspice DC (`diode_series_rd`), and equal to an external resistor of the same value at 0.4/0.7/1.0 V (`rs_in_the_model_equals_an_external_resistor`) |
| `CJO` / `CJ0` | ✅ | ✅ | ✅ ngspice transient; equivalence test vs a discrete `C` |
| `VJ` | ✅ | ✅ | ⚠️ transitively via `CJO` |
| `M` / `MJ` | ✅ | ✅ | ⚠️ transitively; `M=0` is exercised directly by the integrator equivalence test |
| `FC` | ✅ | ✅ | ⚠️ transitively |
| `TT` | ✅ | ✅ | ⚠️ transitively (transit-time charge) |
| `BV` | ✅ | ✅ reverse breakdown, knee adjusted so `I(-BV) = IBV` | ✅ ngspice across the knee, and a Zener shunt regulator |
| `IBV` | ✅ | ✅ default 1 mA | ✅ ngspice at three values, plus the clamp below the leakage floor |
| `TNOM`, `EG`, `XTI` | ✅ | ✅ `IS(T)`, `VJ(T)`, `CJO(T)` — see *Temperature* below | ✅ ngspice at −40/27/75/125 °C, at three (EG, XTI) pairs, and `Cj` at two `M` |
| `KF`, `AF` | ✅ | ✅ flicker noise, `KF·|Id|^AF / f` across the junction, in `.noise` and transient noise | ✅ ngspice at two `AF` values, slope and magnitude |
| `ISR`, `NR` | ✅ | ✅ recombination current, with SPICE's generation factor | ✅ ngspice at four `VJ`/`M` pairs and three biases; `NR` is a **deliberate divergence**, see below |
| `IKF` | ✅ | ✅ `Id/(1 + sqrt(Id/IKF))`, on the **total** forward current | ✅ ngspice at twelve points, residual 1e−6 |
| `TRS1`, `TRS2` | ✅ | ✅ `RS(T) = RS·(1 + TRS1·dT + TRS2·dT²)` | ✅ ngspice at four temperatures and three cards, residual ~1e−6 |
| `CTA`, `VPT` | ⚠️ accepted, not modelled — **and ngspice ignores both**, see below | ❌ | ✅ measured inert in the reference |

Anything not on either list is an unknown parameter and is warned about as one.

### Reverse breakdown

`BV` and `IBV` are modelled. Past `-BV` the current is the Shockley exponential
mirrored about an adjusted knee voltage, so a Zener or an ESD clamp clamps instead
of blocking. Every constant was back-solved from ngspice rather than read from the
SPICE source, and the tests compare at runtime
(`ngspice_diode_breakdown_golden.rs`).

The knee is **adjusted**: the card gives a voltage *and* a current at that
voltage, and both hold at once only if the exponential's offset solves
`IS·(exp((BV − bv_adj)/vte) − 1 + bv_adj/vte) = IBV`. With `BV=5, IBV=1 mA,
IS=10 fA` that is 4.3449 V, which predicts ngspice's current at 4.5 V
(−4.02313e−12) exactly. The exponential's slope is `1/(N·vt)`, fitted from
ngspice rather than assumed. When `IBV` falls below the leakage the card already
has at `-BV` there is no offset to solve for, and `bv_adj` is `BV` unshifted —
ngspice does the same.

**One divergence, on purpose.** `AREA` scales the breakdown current here, so
`area=N` equals N diodes in parallel. ngspice's breakdown branch is *exactly
independent* of `area` — measured at 4.8, 5.0, 5.1 and 5.3 V, ratio 1.0000, while
its forward current doubles correctly — because deriving the knee offset from
`IS·AREA` doubles the prefactor and lifts the offset by `vte·ln 2`, and the two
cancel. ngspice then disagrees with itself: two diodes in parallel give exactly
twice the breakdown current of one `area=2` diode. This tree's rule is already
written down in `area_scales_the_diode_exactly` — "AREA=2 *is* two devices" — and
an `area=10` Zener silently carrying a tenth of its knee current is the failure
this codebase refuses.

**Mild reverse also diverges, and here fairchild is the exact one.** ngspice
smooths its reverse saturation with a cubic fit, `-IS·(1 + (3·vte/(vd·e))³)`,
which reads 1.86e-4 low against the Shockley law it is fitting at −0.5 V.
fairchild evaluates the law. The test asserts both halves — fairchild against the
closed form, and the gap to ngspice against ngspice's own fit — so the difference
stays the one term we know about.

### Instance parameters

| Parameter | Parsed | Stamped | Validated |
|---|---|---|---|
| `AREA` | ✅ | ✅ scales IS and CJO, divides RS | ✅ exact agreement with N parallel diodes, incl. with RS |

Any other instance parameter is named on stderr per instance rather than
dropped. A value that cannot be read (`area=2x`) is a parse error — it used to
be discarded, which read as "this simulator ignores AREA".

Shot noise (`2q·|Id|`) is stamped for `.noise`.

### Temperature

`crate::temperature` owns every temperature law, because the same bandgap
expression appears in the diode's saturation current, the BJT's, and the MOSFET's
threshold, and three copies are three chances to disagree. Devices hold a
*factor* derived once per solve, derived idempotently from the nominal value so a
second `setup_model` cannot compound it and so `AREA` still multiplies on top.

Every law was **back-solved from ngspice**, not read from the SPICE source, which
carries several variants of each:

| law | form | how it was measured |
|---|---|---|
| diode `IS(T)` | `IS·exp((T/TNOM−1)·EG/(N·vt))·(T/TNOM)^(XTI/N)` | reverse saturation at −1 V *is* `−IS(T)`; fitted `XTI` 2.9999, `EG` 1.1100 against defaults of 3 and 1.11 |
| BJT `IS(T)` | `IS·(T/TNOM)^XTI·exp(EG·(T/TNOM−1)/vt)` | `IC` at fixed `V_BE`, dividing out `exp(V_BE/vt(T))`; residual 1.6e−6 |
| BJT `BF(T)` | `BF·(T/TNOM)^XTB` | `IC/IB`, which cancels `IS(T)`; 1.5278 measured against 1.5278 |
| MOSFET `KP(T)` | `KP·(T/TNOM)^−1.5` | slope of `sqrt(Id)` against two `Vgs`, which separates it from the threshold |
| MOSFET `PHI(T)`, `VTO(T)` | SPICE3's `pbfact`/bandgap form | the same fit's intercept; residual 0.04 mV |
| junction potential `VJ(T)`, `PB(T)`, `VJE`/`VJC` | the same law as `PHI(T)` | reused, not rewritten |
| junction capacitance `CJO(T)`, `CJ`/`CJSW`, `CJE`/`CJC` | `cjfact·cjfact1` about `pbo` | 1.2e−6…8.2e−5 across −40/27/75/125 °C at two `M` |

Note the **`/N`**: the emission coefficient divides both temperature terms in the
diode law and neither in the BJT's. The two coincide at `N = 1`, which is why
`the_diode_law_divides_by_n` exists — it is the only test that separates them, and
sabotage confirms it is the only one that catches the difference.

`k` and `q` in `temperature.rs` are SPICE3's values rather than the more accurate
ones used elsewhere in this tree. The bandgap expressions are curve fits whose
coefficients were published against those constants, and a better `k` moves
`PHI(T)` in the fourth decimal and stops matching ngspice.

### The junction potential and capacitance

`VJ(T)` is **the same function** as `PHI(T)` — a MOSFET's surface potential and a
diode's junction potential are the same quantity in SPICE — so
`scaled_junction_potential` delegates to `scaled_phi` rather than restating it.

The capacitance law has two halves and they are not redundant:

```text
pbo     = (VJ − pbfact(TNOM)) / (TNOM/300.15)
cjfact  = 1 / (1 + M·(4e-4·(TNOM−300.15) − (VJ    − pbo)/pbo))
cjfact1 =      1 + M·(4e-4·(T   −300.15) − (VJ(T) − pbo)/pbo)
CJO(T)  = CJO · cjfact · cjfact1
```

`cjfact` un-references the card from its own `TNOM` and `cjfact1` re-references it
to `T`, so a card extracted at `TNOM` and run at `TNOM` comes back unshifted. That
identity is the property no cross-simulator comparison can check — both would
share an offset — and it is what would move every existing transient golden at
once if the two halves failed to cancel, so it is asserted directly.

`M` appears three times (both corrections and the depletion exponent), which is
why the tests sweep it. A MOSFET carries **two** factors rather than one, because
`CJ` grades by `MJ` and `CJSW` by `MJSW` and SPICE derives each correction from
its own coefficient.

Measured by reading a reverse-biased diode's own `Cj` as the C of a 1 MΩ
single-pole RC. The probe frequency sits near the pole deliberately: a first
attempt at 1 kHz against a ~450 kHz pole attenuated by 6 ppm, where six printed
digits carry no information about `C`.

**`TNOM` is now honoured on all three families**, so it is off the
accepted-but-not-modelled lists. What a card still cannot get is a *temperature
coefficient of its own* — `TRS1`/`TRS2` for a diode's `RS`, and `CTA` — which are
separate parameters and remain listed.

### `gmin`, `RS`, and step limiting

`.options gmin` is a real conductance across the junction — it is in `Id` as well
as in `∂Id/∂V`, so a reverse-biased diode carries `IS + gmin·V` and the leakage
follows the option. It used to be in the slope only, which the Norton form
cancelled out exactly; see `docs/spice_support.md`.

`RS` is solved, not lagged. The junction voltage satisfies
`Vd_j + RS·Id(Vd_j) = Vd_terminal` at every eval, by a scalar Newton inside
`ShockleyDiode::eval`. It used to be one fixed-point step per outer iteration
using `Id` from the previous eval — and `Vd_j` is internal state the outer Newton
cannot see, so its convergence test was satisfied with the lag still open. That
happens immediately whenever a voltage source pins the diode's terminals, since
the visible unknowns then stop moving on the first iteration. Measured against
ngspice, which gives the junction a real internal node: 2.7% low at 0.7 V, and
100% wrong at 1.0 V once `gmin` gave the reverse branch a conductance to converge
onto. `RS = 0` skips the solve entirely, which is every diode that does not set
the parameter.

**There is no mirrored `pnjlim` for the breakdown exponential**, and that is a
decision rather than an omission. Mirroring the forward limiter about `-bv_adj`
looks obviously right — the knee is as steep as forward conduction — and it was
written and then removed, because it produced a silent wrong answer. `vd_prev` is
state the outer Newton cannot see: the mirror compressed the walk into the knee
while a free node jumped straight to the supply, so the stamp kept reading "barely
conducting", and in that region the terminal current is under `abstol`, so the
visible unknowns stopped moving and Newton reported success. Measured at
`.options vmax=1e6`: a 12 V / 1 kΩ Zener regulator read `out = 12 V` with the
mirror and the correct 5.0501 V without it. At the default `vmax` it changed no
answer on any deck, including 1 kV into 10 Ω, which is why it took a non-default
setting to find.

The forward `pnjlim` is safe for the reason the mirror was not: while it is active
the current changes by orders of magnitude per iteration, so a stalled walk cannot
pass the convergence test. Reverse breakdown has a flat plateau under `abstol`
instead. What bounds a step into breakdown now is the trust region
(`vmax + reltol·|v|`, #90), which covers both exponentials and keeps no state.

---

## 4. BJT (`Q` / `.model … NPN|PNP`) — Gummel-Poon Level 1

| Parameter | Parsed | Stamped | Validated |
|---|---|---|---|
| `IS` | ✅ | ✅ | ✅ ngspice DC (3 op-point tests) |
| `BF`, `BR` | ✅ | ✅ | ✅ ngspice DC (forward-active + saturation) |
| `NF`, `NR` | ✅ | ✅ | ⚠️ transitively |
| `VAF` / `VA`, `VAR` / `VB` | ✅ | ✅ | ✅ ngspice IC(VCE) at four bias points — **the sign was inverted until #63, and no golden set VAF** |
| `IKF`, `IKR` | ✅ | ✅ high-injection knee, via the base charge | ✅ ngspice, VBE 0.4→0.9 V |
| `ISE`, `NE`, `ISC`, `NC` | ✅ | ✅ non-ideal junction leakage | ✅ ngspice, same sweep |
| `TF`, `TR` | ✅ | ✅ | ⚠️ transitively (transit-time charge) |
| `RB`, `RC`, `RE` | ✅ | ✅ | ✅ ngspice op-point to 6 significant figures |
| `CJE`, `VJE`, `MJE` | ✅ | ✅ | ✅ ngspice transient (CE stage, 5 %) |
| `CJC`, `VJC`, `MJC` | ✅ | ✅ | ✅ ngspice transient |
| `FC` | ✅ | ✅ | ⚠️ transitively |
| `TNOM`, `EG`, `XTI` | ✅ | ✅ `IS(T)`, and `VJE`/`VJC`/`CJE`/`CJC` at temperature | ✅ ngspice at −40/27/75/125 °C |
| `XTB` | ✅ | ✅ `BF(T)`/`BR(T)` | ✅ ngspice, and `XTB=0` pinned as the control |
| `KF`, `AF` | ✅ | ✅ flicker noise, `KF·|Ib|^AF / f` across base-emitter — the **base** current, as SPICE does | ✅ ngspice at two `AF` values |
| `ISS` | ✅ | ✅ substrate junction DC branch, Shockley | ✅ ngspice DC, six biases over five decades |
| `CJS`, `VJS`, `MJS` | ✅ | ✅ collector-substrate depletion capacitance | ✅ ngspice AC at eight biases, including past `VJS` |
| `FCS` | ⚠️ accepted, not modelled — **and ngspice ignores it too**, see below | ❌ | ✅ measured inert in the reference |
| `XCJC` | ✅ | ✅ splits `CJC` across `RB`, internal share to the internal base | ✅ ngspice AC at six values, monotonic and matching |
| `XTF`, `VTF`, `ITF` | ✅ | ✅ `TF_eff = TF·(1 + XTF·(IF/(IF+ITF))²·exp(VBC/(1.44·VTF)))`, with the capacitance as the charge's analytic derivative | ✅ ngspice AC at nine cards, plus a `VCE` sweep for the transcapacitance and a transient for its linearisation |
| `RBM`, `IRB` | ✅ | ✅ the base resistance falls with base current, two laws selected by `IRB` | ✅ ngspice DC at 16 cards and biases, plus both closed-form limits |
| `PTF` | ⚠️ accepted, not modelled — a frequency-dependent transconductance, see below | ❌ | ✅ measured real in the reference, phase only |

### Instance parameters

| Parameter | Parsed | Stamped | Validated |
|---|---|---|---|
| `AREA` | ✅ | ✅ scales IS, IKF/IKR, ISE/ISC, CJE/CJC; divides RB/RC/RE | ✅ ngspice, and IC exactly doubled |

Any other instance parameter is named on stderr per instance. `build_devices`'
BJT arm did not destructure `params` at all before #26, so nothing on a `Q` line
could reach the device.

`RB` is a constant resistance. The current-dependent base resistance (`RBM`,
`IRB`) is not modelled — and is warned about.

Shot noise on both junctions is stamped for `.noise`, from the terminal currents
of the last `eval`. It used to read the *Norton offsets* instead, which equal
those currents only when every terminal sits at 0 V.

### `gmin`

`.options gmin` is a real conductance across each *modelled* junction, added at
the terminal pair (`gpi`, `gmu`) and carried in `ic`/`ib`. It was previously
folded into `gbe`/`gbc`, which are transport quantities divided by `BF`/`BR` on
the way out — so it arrived as `gmin/100` and was not a conductance across
anything — and the Norton form then cancelled it out of the terminal currents
entirely. A reverse-biased BJT carried `2·IS = 2e-16` where ngspice carries
`gmin·V ≈ 1e-12`.

### `XCJC`, and why the split matters

`XCJC·CJC` connects to the **internal** base node, so `RB` sits in series with it.
The rest hangs off the external base pin, where it does not. All of `CJC` used to
sit outside `RB`, which is the `XCJC = 0` case and 0.575 of the right answer on a
default card. Measured against ngspice, AC `|V(b)|` at 100 MHz behind 1 kΩ with
`RB = 1k` and `CJC = 10p`:

| card | ngspice | relative |
|---|---|---|
| absent | 5.174797e−01 | 1.0000 |
| `XCJC=1.0` | 5.174797e−01 | 1.0000 |
| `XCJC=0.5` | 3.733585e−01 | 0.7215 |
| `XCJC=0.0` | 2.975813e−01 | 0.5751 |

Absent and `XCJC=1.0` being identical is how the default was read off rather than
assumed. With `RB = 0` the two base nodes alias and the split is invisible, which
is every card that gives no base resistance.

### The transit time, and the only transcapacitance in the tree

`TF` used to be constant, so `fT` did not fall at high current and did not rise
with `VCE` — the two things `XTF`/`VTF`/`ITF` exist to describe.

```text
TF_eff = TF·(1 + xf),   xf = XTF·(IF/(IF+ITF))²·exp(VBC/(1.44·VTF))
```

`ITF = 0` and `VTF = 0` each disable their own factor, which is SPICE's convention
and is why both default to zero rather than to infinity. Measured: `TF_eff/TF` is
2.0000, 6.0000 and 11.0000 for `XTF` of 1, 5 and 10 with the other two absent.

**The capacitance is not `TF_eff·gm`.** It is the derivative of the charge, and the
two factors differ:

```text
q  = TF·(1 + xf)·IF/qb
∂q/∂vbe = TF/qb·[gbe·(1 + xf·(3 − 2·tmp)) − (1 + xf)·i_diff·∂qb/∂vbe]
∂q/∂vbc = TF·i_diff·[xf/(1.44·VTF) − (1 + xf)·∂qb/∂vbc/qb]
```

The `(3 − 2·tmp)` comes out of `IF·∂xf/∂vbe = 2·xf·(1−tmp)·gbe`. Extracting a
capacitance from an ngspice divider gave 9.70 where the derivative form predicts
9.93 and the charge form predicts 7.35 — enough to say which form is right, not
enough to pin it, so the tests compare the AC response of a whole deck instead.

`∂q/∂vbc` is a **transcapacitance**: the base-emitter charge varies with the
base-collector voltage. It is the only asymmetric reactance in this tree, so it
cannot be a `ReactiveBranchSpec` and goes through `load_reactive_jacobian`. The
`qb` half of it was missing before, on any card with a finite `VAF` or `IKF`.

It has to appear in three places, and leaving it out of any one of them is a
different fault. Out of the transient Jacobian: a wrong Newton step. Out of the AC
matrix: a wrong bandwidth. Out of the companion's `cv`: a spurious
`scale·cbe_x·vbc` in the converged transient answer, which no AC test can see. The
last one is worth 83% on a switching deck and has its own transient golden.

### The diode's forward characteristic below and above the ideal region

Three laws, all measured from ngspice rather than read.

**Recombination.** Below about 0.4 V the depletion region's recombination current
dominates, which is why a real diode's low-current ideality is nearer two than one.
It was absent, so the low-bias current was the ideal exponential and nothing else.

```text
Irec = ISR·(exp(V/(NR·vt)) − 1)·((1 − V/VJ)² + 0.005)^(M/2)
```

Matched to 1.1e−4…2.7e−4 across `(VJ, M)` of (1.0, 0.5), (0.75, 0.33), (0.6, 0.5)
and (1.0, 0.0), at 0.2, 0.35 and 0.5 V. The `0.005` is what keeps the generation
factor finite at `V = VJ`, where the bare square is zero and its `M/2` power has an
infinite slope.

**High injection.** Above `IKF` the current bends from exponential towards `sqrt`,
because the injected carrier density reaches the doping.

```text
Id = Id_total/(1 + sqrt(Id_total/IKF))
```

Matched to 1e−6 at twelve points. The knee applies to the **total** forward current,
ideal plus recombination, not to the ideal alone: at 0.5 V with `ISR = 1e-8` and
`IKF = 1e-3` ngspice reads 8.555990e−05, against 8.555977e−05 for knee-on-total and
1.143950e−04 for knee-on-ideal. Forward only, because in reverse the square root is
not real and there is no high injection to describe.

**`RS` with temperature.** `RS(T) = RS·(1 + TRS1·dT + TRS2·dT²)`, matched to ~1e−6
at −40, 0, 75 and 125 °C across three cards, extracted by binary search on a fixed
`RS` at the same temperature.

`ISR` scales with `AREA` like `IS` and takes `IS`'s temperature factor. Its own law
would need its own `EG`/`XTI` pair, which the card does not carry.

### `NR` is a deliberate divergence from ngspice

**ngspice ignores `NR`.** Its answer is bit-identical with and without `NR=2`, so it
hardcodes the default. `NR` is honoured here as a real parameter, which agrees with
ngspice on every card ngspice can represent — `NR` absent or 2 — and honours one it
cannot. A card that asks for `NR=1.5` gets `NR=1.5` rather than a silent 2.

This is the same shape as the diode's `AREA` in reverse breakdown, where ngspice
disagrees with its own parallel pair. Where the reference is self-inconsistent or
drops a parameter it accepts, this simulator follows the parameter.

### `CTA` and `VPT`, which ngspice ignores too

`CTA` would be a linear temperature coefficient on `CJO`. The capacitance at 1 V
reverse is bit-identical for `CTA` of 1e−3, 1e−2 and absent, at 27, 75 and 125 °C.
The junction capacitance *does* move with temperature here, by SPICE's `cjfact`
law — 6.546536e−12 at 27 °C against 6.866055e−12 at 125 °C — and that is the same
movement ngspice makes. So `CTA` is inert in the reference, not a missing
temperature dependence.

`VPT` would be a punch-through voltage in reverse breakdown. The current past a
`BV=50` knee is bit-identical for `VPT` of 10, 40, 49 and absent, so there is no
reference to match. Breakdown itself is modelled, and a card relying on
punch-through gets the plain `BV` knee.

Both stay on the unmodelled list with those measurements as their reasons, the same
call as the BJT's `FCS` and the MOSFET's LEVEL 2/3 mobility group.

### The base resistance, which is not constant

`RB` used to be a fixed resistance. A real transistor's base resistance collapses
under drive, and without that a hard-driven stage reads too little collector
current. Two laws, and `IRB` alone selects which:

```text
IRB == 0   rb_eff = RBM + (RB − RBM)/qb
IRB >  0   rb_eff = RBM + 3·(RB − RBM)·(tan z − z)/(z·tan²z)
           z = (sqrt(1 + 144/pi²·IB/IRB) − 1) / ((24/pi²)·sqrt(IB/IRB))
```

`RBM` defaults to `RB`, which makes both laws return `RB` exactly — measured, not
assumed: ngspice gives bit-identical currents for `RB=10k` and `RB=10k RBM=10k`.
Any other default would silently change every card that sets `RB`.

Extracted from ngspice by binary search on a *fixed* `RB` that reproduces the same
collector current, so the extraction assumes no law:

| `Vb` | `IRB` | `IB` | extracted | the `tan z` law |
|---|---|---|---|---|
| 0.9 | 1e−6 | 7.0213e−05 | 973.3028 | 973.7307 |
| 1.0 | 1e−5 | 6.5150e−05 | 2621.3721 | 2622.2215 |
| 1.5 | 1e−4 | 1.4017e−04 | 4590.6510 | 4591.7263 |

`rb_eff` is a function of the iterate's own node voltages, through `IB` and `qb`, so
it is **not** lagged internal state and the convergence test sees it stop moving when
`x` does. The Jacobian stamps `1/rb_eff` without differentiating it, which costs
Newton steps and not correctness: the residual defines the answer and both read the
same value.

Both limits of the `tan z` form are numerically hostile and each has a test. The
bracket is written as `1/(z·tan z) − 1/tan²z`, a difference of two terms of size
`z⁻²` whose difference is `1/3`, so the relative cancellation is `3z²` — eight
digits left at `z = 1e-4`, none at `1e-8`. Below `1e-4` the series limit `1/3` is
used instead. At the other end `tan z` overflows, and the same spelling makes both
terms underflow to zero rather than forming `∞/∞`.

### `PTF`, and why it is not modelled

Excess phase is real in ngspice rather than inert. At 1 GHz with `TF = 1n`, so
`1/(2·pi·TF)` is 159 MHz:

| card | `ph(v(c))` | `mag(v(c))` |
|---|---|---|
| absent | 1.678594e+00 | 1.582311e+00 |
| `PTF=0` | 1.678594e+00 | 1.582311e+00 |
| `PTF=30` | −1.611270e+00 | 1.582311e+00 |
| `PTF=90` | −1.907830e+00 | 1.582311e+00 |

Phase only, magnitude bit-identical, which is what excess phase should do. So this
is a gap and not an *ngspice ignores it too* case like `FCS`.

It is left out because it is different in kind from everything else here. A
capacitance is `jw·C` and fits the `G + jwC − L/w` assembly; excess phase is a
delay on the transport current, so its contribution is a frequency-dependent
*transconductance*. Expressing it needs a hook the small-signal assembly does not
have. A stage past `1/(2·pi·TF)` gets the right gain and the wrong phase, which the
warning says.

### The collector-substrate junction

Three terminals used to be all this model had. `terminals[3]` was dropped, so a
reverse-biased BJT read `1·gmin·V` against ngspice's `2·gmin·V` — the substrate
junction was the missing one. Confirmed by pinning the substrate at the collector
potential in ngspice, which removes exactly one `gmin·V` and nothing else. That is
closed: the junction is stamped between the substrate and the **internal**
collector, so a card with `RC` puts the series resistance between the two, which is
where SPICE puts it.

The junction exists whether or not the card gives it anything. ngspice's leakage is
`2·gmin·V` for a bare `IS`/`BF` card, so `gmin` crosses it with no `CJS` and no
`ISS`, and it does here too.

| | law | measured against ngspice |
|---|---|---|
| DC branch | `ISS·(exp(V/vt) − 1) + gmin·V` | six biases over five decades |
| reverse capacitance | `CJS·(1 − V/VJS)^−MJS` | 5e−8 at 0, −0.5, −1, −3 V |
| forward capacitance | `CJS·(1 + MJS·V/VJS)` | 1.7e−7 at 0…2 V |

Defaults measured, not read: `CJS = 0`, `MJS = 0`, `VJS = 0.75`, `ISS = 0`. So a
card that names none of them gets the junction, its `gmin`, and nothing else.

The DC branch is plain Shockley, **not** the flat reverse branch ngspice's MOS1
bulk diodes use. At −0.05 V with `ISS = 1e-15` ngspice reads 8.553040e-16 where
Shockley gives 8.553119e-16 and a flat `−ISS` would give 1e-15. The two junction
families in SPICE genuinely differ, and each was measured rather than assumed.

**`FCS` is inert in ngspice.** The forward capacitance is a straight line from the
*zero-bias* value, not from `FCS·VJS`, so `FCS` never enters. The capacitance at
0.5 V forward is bit-identical for `FCS` of 0.1, 0.5, 0.9 and absent. It stays on
the unmodelled list with that as its reason, because honouring it would move away
from the reference. The forward law holds out to 2 V, well past `VJS`, where the
depletion law is singular — which is what the linearisation is for.

`CJS` takes no temperature factor. `TNOM` moves `CJE` and `CJC` here, and the
substrate junction stays at its nominal value because nothing has measured the law
for it.

---

## 5. MOSFET (`M` / `.model … NMOS|PMOS`) — Level 1 Shichman-Hodges

### Model-card parameters

| Parameter | Parsed | Stamped | Validated |
|---|---|---|---|
| `VTO` / `VTH0` / `VTHO` | ✅ | ✅ | ✅ ngspice DC sweep |
| `KP` | ✅ | ✅ | ✅ ngspice DC sweep |
| `LAMBDA` | ✅ | ✅ | ⚠️ transitively |
| `GAMMA`, `PHI` | ✅ | ✅ | ⚠️ transitively (body effect) |
| `CGSO`, `CGDO`, `CGBO` | ✅ | ✅ | ✅ ngspice switching-time golden |
| `COX` | ✅ | ✅ | ⚠️ transitively (Meyer channel caps) |
| `TOX` | ✅ | ✅ | ⚠️ converted to `COX`, transitively |
| `CJ`, `CJSW` | ✅ | ✅ | ✅ ngspice switching-time golden |
| `PB`, `MJ` | ✅ | ✅ bottom of the junction | ⚠️ transitively |
| `FC` | ✅ | ✅ | ⚠️ transitively |
| `MJSW` | ✅ | ✅ sidewall, graded separately from `MJ` | ✅ closed form with `CJ=0`, where `MJ` cannot substitute |
| `KF`, `AF` | ✅ | ✅ flicker noise, `KF·|Id|^AF / (f·W·L·COX)`. A card with `KF` and no `TOX`/`COX` is refused by name — the density's denominator would be zero | ✅ closed form and structure; ngspice is **not** an anchor here, see below |
| `RD`, `RS` | ✅ | ✅ a real internal node each, `1/R` between the external terminal and it — not an analytic elimination, see below | ✅ ngspice DC, and equal to an external resistor of the same value |
| `TNOM` | ✅ | ✅ `KP(T)`, `PHI(T)`, `VTO(T)`, `PB(T)`, `CJ`/`CJSW`, and the bulk junctions' `Isat(T)` | ✅ ngspice at −40/27/75/125 °C, threshold and mobility separated by a two-point fit |
| `UO` | ✅ | ✅ derives `KP = UO·COX` when the card gives no `KP` | ✅ ngspice: `UO=300` gives exactly half the drain current of the 600 default |
| `IS`, `JS` | ✅ | ✅ bulk-source and bulk-drain diodes, `JS·AS` / `JS·AD` when the area is given, else `IS` | ✅ ngspice DC forward and reverse, plus a closed-form anchor on the reverse branch |
| `RSH`, `NSUB`, `NSS`, `TPG`, `XJ`, `LD`, `DELTA`, `PHP` | ⚠️ accepted, not modelled | ❌ | ✅ warning text pinned |
| `UCRIT`, `UEXP`, `UTRA`, `VMAX`, `NFS`, `THETA`, `ETA`, `KAPPA` | ⚠️ accepted, not modelled — **and not part of LEVEL 1**, see below | ❌ | ✅ warning text pinned |

### The body diodes, and where they differ from ngspice

The bulk-source and bulk-drain junctions are real pn junctions, so `gmin` crosses
them. That is what makes the three families consistent: before this a MOSFET had
no junction at all, and its `gmin` was a Jacobian-only channel floor while the
diode's and the BJT's were conductances. A reverse-biased MOSFET now leaks
`2·(Isat + gmin·V)`, which is ngspice's answer at every `gmin` tried.

The law is `Isat·(exp(V/vt) − 1) + gmin·V`, the same one
`ShockleyDiode::junction` uses. **ngspice's is not, and the difference is worth
naming.** ngspice's MOS1 reverse branch is flat at exactly `−Isat` from `−3·vt`
outward, and inside `±3·vt` its total over the two junctions measures as one
junction flat and one plain Shockley — matched to seven digits, and still there
with the bulk-drain junction held five volts reverse. That asymmetry is a
numerical convenience in the reference, not physics.

Both pure choices sit the same distance from it. At −0.01 V with `IS = 1e-14`,
ngspice reads 1.32e-14 A, Shockley-on-both 6.4e-15 and flat-on-both 2.0e-14 — a
6.8e-15 A difference either way. Outside `±3·vt` the question disappears: `exp(V/vt)`
underflows toward zero and Shockley *is* the flat answer, to 4.4e-4 relative at
−0.2 V and exactly by −0.5 V. So this takes the smooth branch. It is the junction
law, it is C¹ at zero where the flat form has a kink, and one law lives in one
place.

One consequence to know about: the junction current appears in a drain
measurement. A reverse-biased bulk-drain junction adds `Isat + gmin·V` to
`I(vd)`, about 3 pA at the default `gmin`, which is 9e-9 of a 0.3 mA drain
current and is why the Level-1 channel tests carry a 1e-6 tolerance rather than
1e-9.

`Isat` also moves with temperature, by a **third** law — neither the diode's nor
the BJT's:

```text
Isat(T) = Isat · exp(Eg(TNOM)/vt(TNOM) − Eg(T)/vt(T))
```

The other two families use a constant `EG` from the card in `exp(EG·(T/TNOM −
1)/vt(T))`, times `(T/TNOM)^XTI`. A MOSFET card carries neither `EG` nor `XTI`,
so SPICE puts the temperature-dependent bandgap in the exponent instead. Reusing
the diode's law here would be out by up to 2.4× over −40 to 125 °C. Measured
against ngspice at five temperatures spanning five decades of `Isat`, worst
residual 3.7e-4.

### `RD`/`RS`, and why they are rows rather than an elimination

Each non-zero series resistance allocates an **internal node**, exactly as the
BJT's `RB`/`RC`/`RE` do, and stamps `1/R` between the external terminal and it.
With no resistance the internal node aliases the external one, so a card without
them allocates no extra rows and stamps no extra conductances.

Not an analytic elimination, deliberately, and the diode's `RS` is why: it
eliminated the series drop by iterating on a junction voltage the outer Newton
could not see, read 2.7% low against ngspice, and the convergence test had no way
to notice — see *`gmin`, `RS`, and step limiting* under the diode. A row costs one
unknown. A hidden state costs a silent wrong answer.

`RSH` stays unmodelled: it is a *sheet* resistance and needs `NRD`/`NRS` (numbers
of squares) to become a resistance at all, and those are instance parameters this
model does not take. A card giving `RSH` without `RD`/`RS` is still told.

Only Level 1 exists. There is no `LEVEL` parameter and no BSIM — for foundry
PDKs the answer is the OSDI/Verilog-A path (see user guide §14).

### The mobility-degradation group belongs to LEVEL 2/3, not here

`UCRIT`, `UEXP`, `UTRA`, `VMAX`, `NFS`, `THETA`, `ETA` and `KAPPA` are on the
accepted-not-modelled list above, and that is not a gap in this Level 1
implementation — **they are not Level 1 parameters.** `UCRIT`/`UEXP`/`UTRA`/`NFS`
belong to SPICE's LEVEL 2 and `THETA`/`ETA`/`KAPPA` to LEVEL 3; `VMAX` to both.
Shichman-Hodges has no field-dependent mobility and no subthreshold region.

Measured, because the list read like a to-do: at LEVEL 1 ngspice's drain current
is **bit-identical** with and without every one of them —

| variation | ngspice LEVEL 1 | ngspice LEVEL 3 |
|---|---|---|
| `THETA=0.1` | 1.0000× | 0.9259× |
| `ETA=0.1` | 1.0000× | 1.3854× |
| `VMAX=1e5` | 1.0000× | 0.6917× |
| `NFS`, `UCRIT`, `UEXP`, `UTRA`, `KAPPA` | 1.0000× | — |
| `LAMBDA=0.05` *(control)* | **1.1500×** | — |

`LAMBDA` is the control: a genuine Level 1 parameter, and it moves the current, so
the probe works. At LEVEL 3 three of the group move it, which is where they live.

So implementing them here would move *away* from the reference. fairchild's
behaviour is already ngspice's, with the one difference that fairchild **says so**
and ngspice is silent. A deck that needs them wants LEVEL 2/3 — which is a separate
model, not a parameter — or the OSDI path (§9b), where BSIM and PSP carry the real
short-channel physics.

`UO` is the exception and is modelled: SPICE derives `KP = UO·COX` when the card
gives no `KP`, which is real Level 1 behaviour. Measured — with `TOX=20n` and no
`KP`, `UO=300` gives exactly half the current of the 600 default, and `UO·COX`
reproduces ngspice's 3.315020e-4 A exactly.

### Flicker noise, and why ngspice is not the anchor for it

`KF`/`AF` give `KF·|Id|^AF / (f·W·L·COX)`, the documented SPICE3 form, with `Id`
stored at eval time rather than reconstructed from the Norton offset (`bjt.rs`
carried exactly that bug once). `LD` is unmodelled, so `Leff` is the drawn `L`.

A card with `KF` and no `TOX` or `COX` is a **hard error naming both**: the
density's denominator is `W·L·COX`, and `COX` is zero unless the card gives one,
so the alternative is a non-finite noise density reaching the matrix.

ngspice's MOS1 flicker density could not be used as the reference. Measured at
`KF=1e-24, W=10u, L=1u, TOX=20n`, it returns **3.706770e-11 V²/Hz at every one
of `AF` = 0.5, 1.0, 1.2 and 2.0** — bit-identical, so its `AF` does nothing. Over
the same sweep the *diode's* `AF` moves as a clean power law, so this is a
property of ngspice's MOS1 rather than of the deck or the card syntax. It also
scales as `W¹·L⁻³` where the documented form gives `W⁰·L⁻²` (`W⁰` because
`Id ∝ W/L` cancels the `W` in the denominator).

Asserting ngspice's number would mean asserting a density that ignores a
parameter the card sets. So `mosfet_flicker_matches_the_closed_form` checks the
closed form, and `the_mosfet_normalisation_is_read` checks that `W`, `L` and `TOX`
each move the answer. The diode and the BJT *are* anchored on ngspice, which
agrees with both to 5e-3.

### `gmin`

`.options gmin` reaches this device as a floor under the **channel**
conductance (`gds`), and it is deliberately Jacobian-only: `jeq` subtracts
`gds_total·vds`, so the term cancels exactly out of the terminal current and the
operating point does not depend on it. That is what a conditioning floor is for,
and it is not what `gmin` does on the diode and BJT, where it crosses a pn
junction and carries current.

The reason for the difference is `IS`/`JS` on the unmodelled list above: **there
are no body diodes here**, so this device has no pn junction to put `gmin`
across. A reverse-biased drain-bulk carries no `gmin` leakage where ngspice's
does. What changed is only that the value now comes from the solve rather than
from a `const` in `mosfet1.rs`, so a deck raising `gmin` to get a stubborn
circuit through gets help from its MOSFETs as well as its junctions.

### Instance parameters

| Parameter | Parsed | Stamped | Validated |
|---|---|---|---|
| `W`, `L` | ✅ | ✅ | ✅ ngspice DC |
| `AS`, `AD`, `PS`, `PD` | ✅ | ✅ | ✅ via `CJ`/`CJSW` golden |

Any other instance parameter is named on stderr per instance. The return of
`set_instance_params` was discarded before #26, so `M1 … banana=3` was accepted
in silence on the one device family whose instance parameters did work.

Channel thermal noise (`8kT·gm/3`) is stamped for `.noise`.

---

## 6. Switches (`S`, `W` / `.model … SW|CSW`)

| Parameter | Parsed | Stamped | Validated |
|---|---|---|---|
| `RON` | ✅ | ✅ | ✅ ngspice DC (both `S` and `W`) |
| `ROFF` | ✅ | ✅ | ✅ ngspice DC |
| `VT` (`SW`) / `IT` (`CSW`) | ✅ | ✅ | ✅ ngspice DC, on both sides of the boundary |
| `VH` (`SW`) / `IH` (`CSW`) | ✅ | ✅ | ⚠️ in-tree (unit tests pin the band against measured ngspice behaviour) |
| `ON` / `OFF` instance keyword | ✅ | ✅ | ✅ in-tree |

Non-positive `RON`/`ROFF` is a hard error, not a clamp. Switching is a hard step
with no smoothing, matching ngspice. **Known divergence:** `.dc` sweep points run
in parallel and do not inherit each other's hysteresis state; transient does.
With the default `VH = 0` there is no difference.

---

## 7. Behavioral source (`B`)

`V=<expr>` and `I=<expr>` with the full expression grammar: `+ - * / ^`,
`sin cos tan asin acos atan sinh cosh tanh exp ln log10 sqrt abs sgn ceil floor
min max pow atan2 if`, SPICE suffixes (`meg k m u n p f g t`), `V(node)`,
`V(a,b)`, `I(vsource)`, and `.param` references. Parsed, stamped, and exercised
transitively by the integration tests. No ngspice comparison.

---

## 8. Transmission line (`T`) — lossless

| Parameter | Parsed | Stamped | Validated |
|---|---|---|---|
| `Z0` | ✅ | ✅ | ✅ ngspice (matched, open, shorted far end) |
| `TD` | ✅ | ✅ | ✅ ngspice |
| `F`, `NL` | ✅ | ✅ | ⚠️ desugars to `TD` |

Lossy lines (LTRA-style loss and dispersion) are **not** implemented.

---

## 9. Photonic devices

All native Rust, all using the slowly-varying-envelope `(re, im, λ)`
representation. λ is a **label**, not an unknown: every device declares where a
wavelength enters it, where it goes, and which of its terminals carry one, and
the wavelength on every λ net is resolved once before the solve. A device reads
its channel wavelength from that resolution, so nothing differentiates against
λ and no device chooses between two λ wires at run time. **None of these is
validated against an external simulator** —
the tests are analytic closed forms, bit-for-bit characterisation pins, and
equivalence tests between two spellings of the same circuit. That is a real
gap. The cheapest way to close it is a cross-check of the linear passives
against a frequency-domain S-matrix tool such as SAX, which shares no code and
no author with this one; the active devices have no obvious external reference,
since a time-domain electro-optic simulator is what fairchild exists to be.

Every parameter listed here is both parsed and stamped unless the note says
otherwise. What each parameter *means*, and which tier to pick, is in
[Photonic models](photonic-models.md).

| Device | Parameters | Validated |
|---|---|---|
| `fc_cw_laser` | `power_mW`, `power_W`, `wavelength_nm`, `wavelength_m`, `phi_0_deg`, `re_amp`, `im_amp`, `rin_db_hz` | ⚠️ analytic |
| `fc_driven_laser` | `slope_w_v`/`slope_mw_v`, `v_th`, `p_floor_w`, `r_in`, `phi_0_deg`, `wavelength_nm`, `rin_db_hz` | ⚠️ analytic L–V + a transient through a PD; the opto-electronic-loop test is what pins the stamped `dA/dV` |
| `fc_facet` | `reflectance`/`r`, `transmittance`/`t`, `loss`, `phase_deg` | ⚠️ analytic (`native_facet`): power fraction, phase rotation, round-trip loss, budget rejection |
| `fc_waveguide` | `l_um`, `l_m`/`length`, `n_g`, `n_eff`, `alpha_db_cm`, `wl_ref_nm`/`wl_ref_m`, `pin_at_ref` | ⚠️ analytic closed form (loss, phase, group delay) |
| `fc_dcoupler` | `kappa_per_m`, `l_um`, `l_m`, `kappa_l` | ⚠️ analytic |
| `fc_splitter` | `alpha`, `alpha_db`/`il_db`, `r`/`split_ratio` | ⚠️ analytic (`native_splitter_asymmetric`) |
| `fc_grating_coupler` | `alpha_db`/`il_db`, `alpha` | ⚠️ analytic |
| `fc_circulator` | (topology only) | ⚠️ analytic (`native_circulator`): each of the three routes, plus a power-conserving loop in `bidirectional_composition` |
| `fc_mux` / `fc_demux` | `il_db`, `lambda0_nm`, `df_ghz`, `fwhm_ghz`, `shape_p`, `dlambda_dt_pm_per_k`, `t_nom_k` | ⚠️ analytic; the filter is opt-in and the default is bit-for-bit identity. The backward route is pinned by a mirror-terminated power budget (`bidirectional_composition`) |
| `fc_awgr` | `lambda0_nm`/`lambda0_m`, `df_ghz`, `fsr_ghz`, `fwhm_ghz`, `shape_p`, `il_db`, `il_tilt_db`, `xt_adj_db`, `xt_bg_db`, `dlambda_dt_pm_per_k`, `t_nom_k`, `sfile` | ⚠️ analytic + a permutation-equivalence test vs demux/mux |
| `fc_mzm` | `v_pi`, `alpha`, `alpha_db`/`il_db`, `e_r`, `e_r_db`, `f_c` | ⚠️ analytic |
| `fc_optical_2x2` | `s11`, `s12`, `s21`, `s22`, `il_db`, `tau_s`, `allow_gain`, per-channel `w`/`dw_dv_<k>` | ⚠️ analytic + power conservation |
| `fc_photodetector` | `responsivity`, `i_dark`, `r_shunt`, `r_series`, `c_par`/`c_j0` | ⚠️ analytic (`native_pd_r_series`); owns no optical wire, and sums both propagation directions into one photocurrent |
| Phase shifters (`fc_pn_ps`, `fc_thermal_ps`, `fc_pn_th_ps`, … + `LEVEL`) | `dn_dv`, `da_dv`, `g_pn`, `v_pi_l`, `c_j0`, `v_bi`, `m_j`, `i_sat`, `n_diode`, `tau_carrier`, `dn_dv_inj`, `da_dv_inj`, `dn_di`, `da_di`, `beta_tpa`, `a_eff_m2`, `r_th`, `dn_dt`, `r_series`, `r_heater`, `p_pi`, `tau_th`, plus the `fc_waveguide` geometry set | ⚠️ characterisation pins + an equivalence test vs a discrete `C`; the LEVEL-4 optical back-action (`beta_tpa`, `r_th`, `dn_dt`) is pinned on a WDM bus against hand-computed absorbed-power and per-channel cross-TPA budgets, forward and backward (`native_pn_ps_full.rs`, `bidirectional_composition.rs`); fitted against measured chip data in `experiments/giona` but that comparison is not a committed regression |
| `fc_phase_shifter_expr` | `dneff`, `dalpha` (expression strings over `V`, `T`, `lambda`) | ⚠️ in-tree |

### Photonic capability gaps

| Gap | Status |
|---|---|
| Optical noise — laser RIN, PD shot noise | ✅ `2q(I_ph+I_dark)` from `fc_photodetector`, `RIN·P²` from the lasers (`rin_db_hz`, off unless set), in **both** `.noise` and `.tran` (`.options trannoise=1`). Both flat with frequency; no relaxation-oscillation peak, no APD excess-noise factor |
| Bidirectional propagation | ⚠️ infrastructure present (`enable_bidirectional`); `fc_facet` is the only source of backward light. No distributed backscatter, and no laser sensitivity to back-reflection — lasers absorb what returns, they do not respond to it |
| Waveguide group delay | ⚠️ opt-in via `.options waveguide_delay=1`, **off by default** |
| Reflections at facets | ✅ `fc_facet` (one port, flat with wavelength, needs `enable_bidirectional`). Grating-coupler and interface reflections are still ❌ |

---

## 9a. Verilog-A disciplines

What a `.va` model may declare a node to be, and what fairchild does with it.

| Discipline | Implemented | Validated |
|---|---|---|
| `electrical` | ✅ potential in V, flow in A; `vntol`/`abstol` | ✅ throughout |
| `thermal` | ✅ potential in K (a rise above `$temperature`), flow in W; the row is found from the OSDI descriptor's `units` and bounded by `temptol`, for a port and an internal node alike | ✅ `thermal_discipline.rs` — a self-heated resistor's fixed point against its closed form, and the row's tolerance read off the topology (the answer alone cannot show it: `vntol` on kelvin is *tighter* than needed, so a misclassified row still converges) |
| optical (`optical_bundle`, the bundle-port dialect) | ✅ `E_RE`/`E_IM` per channel on real MNA rows; λ is a label resolved before the solve, not a row | ✅ `bundle_dialect_e2e.rs`, `mrm_wdm_example.rs` |

A thermal network is written with ordinary `R`/`C`/`I`/`V` — on a thermal node
`R` is K/W, `C` is J/K, `I` is watts. There is deliberately no discipline check
refusing them, unlike on optical wires: the electrical primitives *are* the
thermal primitives, and refusing them would ban device-to-device thermal
crosstalk.

Gap: no native (Rust) device exposes a thermal port yet, so a thermal network
can only be entered through a Verilog-A model or the deck's own R/C. `fc_pn_ps`
and friends keep `r_th` as a parameter and their temperature is internal to the
device rather than a shared unknown.

---

## 9b. The OSDI v0.4 ABI surface

A compiled Verilog-A model is only as correct as the parts of the ABI the
simulator honours, and a part left unread fails silently — the model does its
job and the answer is wrong anyway. This is the per-field audit.

| Descriptor surface | Honoured | Validated |
|---|---|---|
| `setup_model` / `setup_instance` | ✅ parameters written first, then setup, which is OSDI's order | ✅ every OSDI test |
| `real` parameters | ✅ | ✅ `osdi_device.rs`, `osdi_model_card.rs` |
| `integer` parameters | ✅ written as `i32`, at the width the model declared | ✅ `osdi_abi_contract.rs` — including a negative value and a neighbouring field used as a clobber witness |
| `string` parameters | ⚠️ cannot be set from a deck; a numeric value is refused with a warning naming the parameter, and the model keeps its default |  |
| instance vs model parameter split (`num_instance_params`) | ✅ | ✅ `inst_param.va` fixture |
| a module name that is not all lower case | ✅ the registry folds case on the way in as well as on the way out — Verilog-A preserves the case an author wrote, SPICE does not distinguish | ✅ `osdi_abi_contract.rs`; `PSP102VA`, `PSP103VA`, `DIODE_CMC`, `JUNCAP200` and `hicumL2va` all used to load and then fail as `unknown model` |
| `$mfactor` from the deck's `m=` | ✅ | ✅ `osdi_device.rs` |
| **node collapsing** (`collapsible` / `collapsed_offset`) | ✅ a collapsed node shares its neighbour's MNA row, or is dropped when collapsed to ground. Its own row stays allocated and unstamped rather than being renumbered — `stamp_gmin` pins it, so the cost is one wasted row per collapsed node | ✅ `osdi_abi_contract.rs` closed form; BSIM4 against ngspice |
| resistive Jacobian (`write_jacobian_array_resist`) | ✅ read as the **packed** array it is: one value per entry carrying `JACOBIAN_ENTRY_RESIST`, not the first `num_resistive` entries | ✅ `osdi_abi_contract.rs` — a fixture whose first entry is reactive-only |
| reactive Jacobian (`write_jacobian_array_react`) | ✅ keyed off `react_ptr_off` | ✅ `osdi_reactive.rs` |
| `load_spice_rhs_dc` / `_tran` | ✅ | ✅ `osdi_dc_op.rs`, `osdi_tran.rs` |
| `load_residual_react` | ✅ charge read back per node for the integrator | ✅ `osdi_reactive.rs` |
| `$limit` (`num_states`, `OSDI_LIM_TABLE`) | ✅ `pnjlim` installed; state buffers swapped per eval | ✅ `osdi_limiting.rs` |
| noise sources (`load_noise`) | ✅ `white_noise` / `flicker_noise` | ⚠️ see §10 `.noise` |
| `thermal` nodes (via `units == "K"`) | ✅ | ✅ `thermal_discipline.rs` |
| operating-point variables (`num_opvars`) | ⚠️ readable through `OsdiDevice::read_opvar`, and nothing surfaces them to a deck or a probe |  |
| `OsdiSimParas` | ✅ `gmin`, `scale`, `shrink`, `simulatorSubversion`. A name not in the table falls back to the model's own default, as every name did before, so the list is monotone. `iteration`/`sourceScaleFactor` are left out deliberately — they change per Newton iteration and would cost an allocation on the hottest path for something no model in the sweep reads; `imax`/`imelt` are left out because fairchild does not enforce them and answering would be a claim it does not keep | ✅ `abi_simparam.va` makes two lookups the *only* conductance on separate branches, with absurd fallbacks, so a missing table is orders of magnitude and not a tolerance |
| `eval`'s return flags | ✅ `EVAL_RET_FLAG_FATAL` reaches the Newton loop through `Device::eval_status`, and joins a clamped step and a stalled line search as a reason this iterate cannot be the converged one. If a device is still saying so when the budget runs out, its message is the error rather than a bare non-convergence. `EVAL_RET_FLAG_LIM` is routine and ignored |  |
| `ANALYSIS_IC` / `ANALYSIS_STATIC` / `ANALYSIS_NODESET` | ❌ never set; a model that branches on them behaves as if none applied |  |
| `given_flag_model` / `given_flag_instance` | ❌ unused; "was this parameter given" is inferred from the deck instead |  |
| init errors (`OsdiInitInfo`) | ⚠️ surfaced, but only `INIT_ERR_OUT_OF_BOUNDS` is decoded; any other code is reported as a bare number |  |
| `bound_step_offset` | ✅ read on the variable-step transient path, where `h` is bounded by the smallest request across devices. Not on the fixed-step path: there the step is what the deck asked for. The sentinel is `u32::MAX`, **not** zero — measured across every fixture, which all report 0xffffffff against a valid 8-aligned 104 for one that calls `$bound_step`. Guarding on zero dereferences 0xffffffff and aborts | ✅ `abi_bound_step.va`, binding and not binding |
| `load_jacobian_resist` (the aliasing path) | ❌ the copy path is used instead, deliberately |  |

### What has been run through it

Every model in OpenVAF-Reloaded's `integration_tests/` was compiled through
fairchild's own `.va` path and then loaded by **ngspice-46 through `pre_osdi`**,
so both simulators evaluate one identical `.osdi`. Any disagreement is a
difference in how the two of them stamp it, with nothing else in the way. DC
sweep, comparison per point, mixed tolerance (a relative metric is meaningless
at the bottom of a sweep where both simulators sit on their own `gmin`).

| Family | Models | Result |
|---|---|---|
| BSIM | BSIM3, BSIM4, BSIMBULK, BSIMCMG, BSIMIMG, BSIMSOI | agree; worst difference 2 pA |
| PSP | PSP102, PSP103, JUNCAP200 | agree; worst 4 pA |
| HiSIM | HiSIM2, HiSIMSOTB | agree; worst 6 pA |
| EKV | EKV, EKV long-channel | agree; worst 1 pA |
| BJT | HICUM/L2, MEXTRAM 505 (3- and 4-terminal) | agree; worst 5 pA |
| Diode | DIODE, DIODE_CMC | agree; worst 1 pA |
| Other | ASM-HEMT, MVSG_CMC, resistor, strings, amplifier, current source | agree; worst 3 pA |

The per-family numbers above were measured **before** `gmin` moved from every
node's diagonal to across each pn junction, and that offset is exactly what they
were: a fixed sub-picoamp difference attributed at the time to "the two
simulators put `gmin` on slightly different sets of nodes". Re-measured for BSIM4
after the change, the disagreement is 9e-14 to 6e-13 relative — femtoamps against
milliamps, i.e. round-off at ngspice's twelve printed digits, and the in-tree
`bsim4_acceptance.rs` tolerance came down from 1e-7 to 1e-11 with it.

The rest of the table has not been re-run (it needs an OpenVAF-Reloaded
checkout), so read those picoamps as upper bounds rather than as current
measurements. The conclusion — every family agrees — is unaffected either way.

Not covered, and why:

- **HiSIM-HV** and the Verilog-A **resistor** fixture: ngspice fails on them
  itself (`"impossible error" in OSDI setup_instance`, and a transient-op
  failure), so there is no reference to compare against. fairchild runs both.
- **VCCS / CCCS**: these expose a DC-solver fault that is not about models at
  all — the `vmax` trust region interacting with the relative convergence test.
  An ideal VCCS into a 1 kΩ load is linear and one solve from its answer, and
  fairchild reports 0.0502 V on a node the deck pins at 0.1 V, as a converged
  operating point. `.options vmax=1e5` gives the right answer. Tracked in #90;
  see the note in `newton.rs` beside the convergence test.
- **MVSG_CMC** looked like a 0.07% disagreement and is not one: ngspice's *own*
  DC sweep drifts from its own operating point at its default `reltol`, and at
  `reltol=1e-10` it lands where fairchild already was.

The comparison needs a checkout of OpenVAF-Reloaded and is not in CI. The ABI
corners it found are covered in-tree by
`crates/fairchild-osdi/tests/models/abi_*.va` and run always.

On macOS/aarch64 every formatted-output and severity task (`$strobe`, `$error`,
`$fatal`, …) is stripped before compiling, loudly, because the compiler
miscompiles the call and the process dies. See `crates/fairchild-osdi/src/portability.rs`.
BSIM4 is 429 such calls, so on that platform the model cannot report anything
about its own parameters — which is why the ABI table above is the only channel
left.

---

## 10. Analyses

| Analysis | Implemented | Validated |
|---|---|---|
| `.op` | ✅ | ✅ ngspice |
| `.dc` | ✅ | ✅ ngspice |
| `.tran` (BE / TR / GEAR, fixed and variable step) | ✅ | ✅ ngspice (RC, RL, RLC, diode, BJT, CMOS, ring oscillator, switch) |
| `.ac` | ✅ magnitude **and** phase (see §2) | ✅ ngspice magnitude + phase on an RC corner; analytic RC Bode and RLC resonance |
| `.noise` | ✅ incl. Verilog-A `white_noise` / `flicker_noise` via OSDI | ⚠️ RC thermal vs analytic + an ngspice spot check; the Verilog-A path is pinned against `(i_n·z_t)²` in `examples/verilog_a/check.py`; device noise unvalidated externally |
| `.temp`, `.alter` | ✅ | ⚠️ transitively |
| `.tf`, `.pz`, `.disto`, `.sens` | ❌ | — |
| `.mc` Monte Carlo | ❌ (use `Circuit.sweep()` from Python) | — |

---

## How to update this document

It is an audit, not a design doc, so it goes stale silently — the worst
property a contract can have. When you add or change a model:

1. Add the row. A parameter with no row is invisible.
2. Be honest in the **Stamped** column. If a parameter is accepted so a foundry
   card loads but changes nothing, say so — that is precisely the information
   this file exists to carry.
3. Name the test in **Validated**, or write ⚠️/❌. "It probably works" is ❌.
