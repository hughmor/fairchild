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
| `RS` | ✅ | ✅ | ✅ ngspice DC (`diode_series_rd`) |
| `CJO` / `CJ0` | ✅ | ✅ | ✅ ngspice transient; equivalence test vs a discrete `C` |
| `VJ` | ✅ | ✅ | ⚠️ transitively via `CJO` |
| `M` / `MJ` | ✅ | ✅ | ⚠️ transitively; `M=0` is exercised directly by the integrator equivalence test |
| `FC` | ✅ | ✅ | ⚠️ transitively |
| `TT` | ✅ | ✅ | ⚠️ transitively (transit-time charge) |
| `BV`, `IBV`, `EG`, `XTI`, `KF`, `AF`, `ISR`, `NR`, `IKF`, `TNOM`, `TRS1`, `TRS2`, `CTA`, `VPT` | ⚠️ accepted, not modelled | ❌ | ✅ warning text pinned |

Anything not on either list is an unknown parameter and is warned about as one.

**Reverse breakdown is not modelled.** A Zener or an ESD clamp will simulate as
an ordinary diode with no breakdown knee. The parameters are accepted so a
foundry card loads, and a warning names them.

### Instance parameters

| Parameter | Parsed | Stamped | Validated |
|---|---|---|---|
| `AREA` | ✅ | ✅ scales IS and CJO, divides RS | ✅ exact agreement with N parallel diodes, incl. with RS |

Any other instance parameter is named on stderr per instance rather than
dropped. A value that cannot be read (`area=2x`) is a parse error — it used to
be discarded, which read as "this simulator ignores AREA".

Shot noise (`2q·|Id|`) is stamped for `.noise`.

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
| `CJS`, `VJS`, `MJS`, `FCS`, `XCJC`, `RBM`, `IRB`, `XTF`, `VTF`, `ITF`, `PTF`, `XTB`, `EG`, `XTI`, `KF`, `AF`, `TNOM` | ⚠️ accepted, not modelled | ❌ | ✅ warning text pinned |

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
| `IS`, `JS`, `RD`, `RS`, `RSH`, `NSUB`, `NSS`, `NFS`, `TPG`, `UO`, `UCRIT`, `UEXP`, `UTRA`, `VMAX`, `XJ`, `LD`, `DELTA`, `THETA`, `ETA`, `KAPPA`, `KF`, `AF`, `TNOM`, `PHP` | ⚠️ accepted, not modelled | ❌ | ✅ warning text pinned |

Note `RD`/`RS`/`RSH` in that list: a MOSFET card's series resistances are
**not** stamped, unlike the BJT's. A card that models its access resistance
there gets none of it, and now says so.

Only Level 1 exists. There is no `LEVEL` parameter and no BSIM — for foundry
PDKs the answer is the OSDI/Verilog-A path (see user guide §14).

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
| `OsdiSimParas` | ❌ a null table is passed, so a model's `$simparam("gmin")` and friends take the default written into the call rather than this simulator's value |  |
| `eval`'s return flags | ❌ discarded — `EVAL_RET_FLAG_FATAL` means the model is telling us this evaluation is invalid, and it is not heard |  |
| `ANALYSIS_IC` / `ANALYSIS_STATIC` / `ANALYSIS_NODESET` | ❌ never set; a model that branches on them behaves as if none applied |  |
| `given_flag_model` / `given_flag_instance` | ❌ unused; "was this parameter given" is inferred from the deck instead |  |
| init errors (`OsdiInitInfo`) | ⚠️ surfaced, but only `INIT_ERR_OUT_OF_BOUNDS` is decoded; any other code is reported as a bare number |  |
| `bound_step_offset` | ❌ unused — a model cannot limit the timestep |  |
| `load_jacobian_resist` (the aliasing path) | ❌ the copy path is used instead, deliberately |  |

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
