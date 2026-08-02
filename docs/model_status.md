# Model status — what is parsed, what is stamped, what is validated

*Audited against the source on 2026-08-01. If you find a disagreement between
this table and the simulator, the simulator is the bug.*

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

**The short version of what to watch out for:**

- BJT high-injection and leakage parameters (`IKF`, `IKR`, `ISE`, `ISC`, …) are
  accepted **silently** and do nothing. The diode at least warns about `BV`.
- MOSFET `MJSW` is parsed, stored, and never read — sidewall junction caps use
  `MJ`.
- `.noise` has no ngspice comparison. Optical noise (PD shot, laser RIN) exists
  and is checked against the analytic receiver budget, but only in `.noise` —
  nothing injects noise into `.tran`.
- The photonic models are validated against analytic forms and against
  themselves, never against an external simulator. That is the biggest gap in
  this document; see `_notes/sotu.md` §I.
- `.ac` ignores source phase, and there is no `AC <mag>` spec on source lines.

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

⚠️ **There is no `AC <mag> [phase]` spec on the source line.** The AC excitation
magnitude is an argument of the analysis, not a property of the source, and by
default *every* independent source is driven at that magnitude; pass a source
name to select one. A netlist written the SPICE way — `V1 in 0 DC 0 AC 1` —
parses, but the `AC 1` is not what sets the excitation.

⚠️ **The AC phase argument is accepted and ignored** (`build_ac_rhs` takes
`_phase_rad`). Every AC source is driven at zero phase, so a multi-source AC
analysis with intended relative phase is wrong.

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
| `BV`, `IBV` | ⚠️ warned | ❌ | — |
| `EG`, `XTI`, `KF`, `AF`, anything else | ⚠️ warned | ❌ | — |

**Reverse breakdown is not modelled.** A Zener or an ESD clamp will simulate as
an ordinary diode with no breakdown knee. The parameters are accepted so a
foundry card loads, and a warning names them.

⚠️ **Diode instance parameters are parsed and then dropped.** `D1 a k dm area=2`
parses, `Element::Diode` carries the list, and nothing reads it —
`ShockleyDiode` does not implement `Device::set_real_param`, so the default
(`false`) applies and `ParamSet::apply` consumes nothing. `AREA` in particular
is not modelled. Nothing warns.

Shot noise (`2q·|Id|`) is stamped for `.noise`.

---

## 4. BJT (`Q` / `.model … NPN|PNP`) — Gummel-Poon Level 1

| Parameter | Parsed | Stamped | Validated |
|---|---|---|---|
| `IS` | ✅ | ✅ | ✅ ngspice DC (3 op-point tests) |
| `BF`, `BR` | ✅ | ✅ | ✅ ngspice DC (forward-active + saturation) |
| `NF`, `NR` | ✅ | ✅ | ⚠️ transitively |
| `VAF` / `VA`, `VAR` / `VB` | ✅ | ✅ | ⚠️ transitively (Early effect) |
| `TF`, `TR` | ✅ | ✅ | ⚠️ transitively (transit-time charge) |
| `RB`, `RC`, `RE` | ✅ | ✅ | ✅ ngspice op-point to 6 significant figures |
| `CJE`, `VJE`, `MJE` | ✅ | ✅ | ✅ ngspice transient (CE stage, 5 %) |
| `CJC`, `VJC`, `MJC` | ✅ | ✅ | ✅ ngspice transient |
| `FC` | ✅ | ✅ | ⚠️ transitively |
| `IKF`, `IKR` | ✅ **silently** | ❌ | — |
| `ISE`, `ISC`, `NE`, `NC` | ✅ **silently** | ❌ | — |
| `CJS`, `VJS`, `MJS`, `XCJC` | ✅ **silently** | ❌ | — |
| `XTB`, `EG`, `XTI`, `PTF`, `TNOM` | ✅ **silently** | ❌ | — |
| `KF`, `AF` (flicker noise) | ✅ **silently** | ❌ | — |

⚠️ **These are accepted with no warning at all.** Unlike the diode, the BJT
matches its unmodelled parameters explicitly and discards them, so nothing is
printed. A foundry card with `IKF` loads and simulates as though high-injection
roll-off did not exist. Verified 2026-08-01: a deck with `IKF=10m ISE=1e-14
KF=1e-15` produces no output on stderr.

⚠️ **BJT instance parameters are parsed and then dropped**, and more thoroughly
than the diode's: `build_devices`' BJT arm does not even destructure `params`.

`RB` is a constant resistance. The current-dependent base resistance (`RBM`,
`IRB`) is not modelled.

Shot noise on both junctions is stamped for `.noise`.

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
| `PB`, `MJ` | ✅ | ✅ | ⚠️ transitively |
| `FC` | ✅ | ✅ | ⚠️ transitively |
| `MJSW` | ✅ | ❌ **stored, never read** | — |

⚠️ **`MJSW` does nothing.** The sidewall junction capacitance is computed with
`MJ`, not `MJSW`; the field is `#[allow(dead_code)]` in `mosfet1.rs`. Cards that
set them differently will not get what they asked for, and nothing warns.

Only Level 1 exists. There is no `LEVEL` parameter and no BSIM — for foundry
PDKs the answer is the OSDI/Verilog-A path (see user guide §14).

### Instance parameters

| Parameter | Parsed | Stamped | Validated |
|---|---|---|---|
| `W`, `L` | ✅ | ✅ | ✅ ngspice DC |
| `AS`, `AD`, `PS`, `PD` | ✅ | ✅ | ✅ via `CJ`/`CJSW` golden |

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
representation. **None of these is validated against an external simulator** —
the tests are analytic closed forms, bit-for-bit characterisation pins, and
equivalence tests between two spellings of the same circuit. That is a real gap
and it is tracked in `_notes/sotu.md` §I; a SAX cross-check for the linear
passives is the cheapest way to close it.

Every parameter listed here is both parsed and stamped unless the note says
otherwise.

| Device | Parameters | Validated |
|---|---|---|
| `fc_cw_laser` | `power_mW`, `power_W`, `wavelength_nm`, `wavelength_m`, `phi_0_deg`, `re_amp`, `im_amp`, `rin_db_hz` | ⚠️ analytic |
| `fc_driven_laser` | `slope_w_v`/`slope_mw_v`, `v_th`, `p_floor_w`, `r_in`, `phi_0_deg`, `wavelength_nm`, `rin_db_hz` | ⚠️ analytic L–V + a transient through a PD; the opto-electronic-loop test is what pins the stamped `dA/dV` |
| `fc_facet` | `reflectance`/`r`, `transmittance`/`t`, `loss`, `phase_deg` | ⚠️ analytic (`native_facet`): power fraction, phase rotation, round-trip loss, budget rejection |
| `fc_waveguide` | `l_um`, `l_m`/`length`, `n_g`, `n_eff`, `alpha_db_cm`, `wl_ref_nm`/`wl_ref_m`, `pin_at_ref` | ⚠️ analytic closed form (loss, phase, group delay) |
| `fc_dcoupler` | `kappa_per_m`, `l_um`, `l_m`, `kappa_l` | ⚠️ analytic |
| `fc_splitter` | `alpha`, `alpha_db`/`il_db`, `r`/`split_ratio` | ⚠️ analytic (`native_splitter_asymmetric`) |
| `fc_grating_coupler` | `alpha_db`/`il_db`, `alpha` | ⚠️ analytic |
| `fc_circulator` | (topology only) | ⚠️ in-tree, bidirectional |
| `fc_mux` / `fc_demux` | `il_db`, `lambda0_nm`, `df_ghz`, `fwhm_ghz`, `shape_p`, `dlambda_dt_pm_per_k`, `t_nom_k` | ⚠️ analytic; the filter is opt-in and the default is bit-for-bit identity |
| `fc_awgr` | `lambda0_nm`/`lambda0_m`, `df_ghz`, `fsr_ghz`, `fwhm_ghz`, `shape_p`, `il_db`, `il_tilt_db`, `xt_adj_db`, `xt_bg_db`, `dlambda_dt_pm_per_k`, `t_nom_k`, `sfile` | ⚠️ analytic + a permutation-equivalence test vs demux/mux |
| `fc_mzm` | `v_pi`, `alpha`, `alpha_db`/`il_db`, `e_r`, `e_r_db`, `f_c` | ⚠️ analytic |
| `fc_optical_2x2` | `s11`, `s12`, `s21`, `s22`, `il_db`, `tau_s`, `allow_gain`, per-channel `w`/`dw_dv_<k>` | ⚠️ analytic + power conservation |
| `fc_photodetector` | `responsivity`, `i_dark`, `r_shunt`, `r_series`, `c_par`/`c_j0` | ⚠️ analytic (`native_pd_r_series`) |
| Phase shifters (`fc_pn_ps`, `fc_thermal_ps`, `fc_pn_th_ps`, … + `LEVEL`) | `dn_dv`, `da_dv`, `g_pn`, `v_pi_l`, `c_j0`, `v_bi`, `m_j`, `i_sat`, `n_diode`, `tau_carrier`, `dn_dv_inj`, `da_dv_inj`, `dn_di`, `da_di`, `beta_tpa`, `a_eff_m2`, `r_th`, `dn_dt`, `r_series`, `r_heater`, `p_pi`, `tau_th`, plus the `fc_waveguide` geometry set | ⚠️ characterisation pins + an equivalence test vs a discrete `C`; fitted against measured chip data in `experiments/giona` but that comparison is not a committed regression |
| `fc_phase_shifter_expr` | `dneff`, `dalpha` (expression strings over `V`, `T`, `lambda`) | ⚠️ in-tree |

### Photonic capability gaps

| Gap | Status |
|---|---|
| Optical noise — laser RIN, PD shot noise | ✅ in `.noise`: `2q(I_ph+I_dark)` from `fc_photodetector`, `RIN·P²` from `fc_cw_laser` (`rin_db_hz`, off unless set). Both flat with frequency; no relaxation-oscillation peak, no APD excess-noise factor, no time-domain noise |
| Bidirectional propagation | ⚠️ infrastructure present (`enable_bidirectional`); `fc_facet` is the only source of backward light. No distributed backscatter, and `fc_cw_laser` *drives* its backward wires to zero rather than absorbing them — putting one at the far end of a reflecting chain over-determines that node |
| Waveguide group delay | ⚠️ opt-in via `.options waveguide_delay=1`, **off by default** |
| Reflections at facets | ✅ `fc_facet` (one port, flat with wavelength, needs `enable_bidirectional`). Grating-coupler and interface reflections are still ❌ |

---

## 10. Analyses

| Analysis | Implemented | Validated |
|---|---|---|
| `.op` | ✅ | ✅ ngspice |
| `.dc` | ✅ | ✅ ngspice |
| `.tran` (BE / TR / GEAR, fixed and variable step) | ✅ | ✅ ngspice (RC, RL, RLC, diode, BJT, CMOS, ring oscillator, switch) |
| `.ac` | ⚠️ magnitude only — the phase argument is ignored (see §2) | ⚠️ analytic only (RC Bode, RLC resonance) — no ngspice waveform comparison |
| `.noise` | ✅ | ⚠️ RC thermal vs analytic + an ngspice spot check; device noise unvalidated externally |
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
