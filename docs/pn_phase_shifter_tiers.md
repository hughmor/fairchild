# PN Phase Shifter — Tiered Model Reference

Fairchild ships several phase-shifter device classes that share a common
optical interface (bundle in/out + 2 electrical terminals) but model
increasingly complete PN-junction physics.  All tiers use the same
defaults for the **passive** waveguide properties (bent-rib at R = 8 µm,
n_eff = 2.7654, n_g = 4.02 at 1550 nm) and the same Soref-Bennett
free-carrier coefficients; they differ in **which voltage-dependent
effects they activate**.

> **One-line rule:** pick the lowest tier that captures the physics your
> circuit is sensitive to.  Use the same `pn_modulator.py` extraction for
> all tiers; the higher tiers just expose more of the same underlying
> physics.

---

## How to populate the parameters

For any tier, run

```bash
cd scripts/waveguide_simulations/pn_modulator
python pn_modulator.py
```

after editing the `Geometry` / `Doping` dataclasses at the top of
`pn_modulator.py` to match your device.  Output `pn_extracted.json` lists
every derived quantity; pick the ones relevant to the tier you're using.

The defaults already shipped in the Rust models come from this script run
with **5e17 N_A / 5e17 N_D, 100 nm junction offset toward the N side, T =
300 K, 1550 nm**.  Change the dataclass values and re-run to get
parameters for any other geometry / doping profile.

---

## L1 — `fc_pn_ps` (testbench, AC small-signal)

**Use when:** wiring up a circuit topology, AC small-signal response,
quick functional check at a fixed bias.

**Physics:**
- Linear electro-optic: φ_eo = 2π · L · (dn/dV) · V_pn / λ
- Single small-signal Δn/ΔV — `dn_dv` (sign and magnitude chosen by user)
- Linear junction conductance `g_pn` — no Shockley I-V
- Constant propagation loss `α` (= waveguide loss + a representative FCA
  value at the chosen bias)
- **No** C_j, **no** bias-dependent loss, **no** TPA, **no** carrier
  dynamics

**Geometry assumption:** symmetric (junction centred on the rib).

**Parameters mapped from `pn_extracted.json` (V_pn ≤ 0 design point):**

| device parameter | source | typical value |
|---|---|---|
| `n_eff`        | `mode.n_eff_rib_straight` | 2.7654 |
| `n_g`          | (constant 4.02 from bent-rib sweep) | 4.02 |
| `wl_ref_nm`    | (designer's choice) | 1550 |
| `alpha_dB_cm`  | `linearised.alpha_at_0V_dB_cm`  | ~3 dB/cm (FCA only); add bend / scattering on top |
| `dn_dv`        | `linearised.dn_dv_reverse_per_V` | 5 × 10⁻⁵ /V |
| `g_pn`         | `1/(R_s)`, where R_s ≈ 30–100 Ω·µm / L | 1 × 10⁻³ S |

---

## L2a — `fc_pn_ps_cap` (depletion-mode, reverse-bias only)

**Use when:** designing a reverse-biased depletion-mode modulator that
needs realistic C_j(V) and FCA(V) for DC + transient simulation in the
**reverse-bias regime only** (V_pn ≤ 0).

**Physics in addition to L1:**
- C_j(V) = C_j0 / (1 − V_pn/V_bi)^m_j (depletion approximation, integrator-
  managed companion model — already implemented).
- Bias-dependent loss α(V) (linear extrapolation via `da_dv` coefficient,
  already in the model).
- Linear Δn_eff(V) (single `dn_dv` is still small-signal).
- Shockley reverse I-V (just leakage; replaces linear `g_pn`).

**Geometry assumption:** 100 nm junction offset toward N.

**Parameters added beyond L1:**

| device parameter | source | typical value |
|---|---|---|
| `c_j0`        | `depletion.C_j_per_um_at_0V_F` × `L_um`  | 7 fF at L = 50 µm |
| `v_bi`        | `depletion.V_bi`                          | 0.92 V |
| `m_j`         | (abrupt: 0.5; linearly graded: 0.33)      | 0.5 |
| `da_dv`       | slope of `delta_alpha_vs_v` for V ≤ 0     | ~50 Np/m per V (negative-going) |
| `i_sat`       | Shockley reverse saturation current; from doping + lifetime | 1 × 10⁻¹² A |

**Out of scope:** forward bias.  Above V_pn ≈ +0.3 V the model is wrong
(it under-counts loss and wildly under-counts dn/dV, because injection
isn't modelled).

---

## L2b — `fc_pn_ps_inj` (carrier-injection, forward-bias only) [NEW]

**Use when:** designing a forward-biased carrier-injection modulator
(used in low-V, lower-speed applications: variable optical attenuators,
slow ring tuners, biomedical).  Operates in V_pn ∈ [0, ~0.8 V].

**Physics in addition to L1:**
- Shockley forward I-V (proper diode + diffusion current).
- Diffusion capacitance C_d = τ · g_d (replaces depletion C_j as V → V_bi).
- Carrier-injection Δn(V): much larger and exponential in V_pn.
  Two-coefficient fit:  Δn_eff(V) = a · (exp(V/V_t) − 1) + b · V .
- High injected-carrier α_FCA(V), also exponential.

**Geometry assumption:** 100 nm junction offset toward N (same as L2a).

**Parameters added beyond L1:**

| device parameter | source | typical value |
|---|---|---|
| `i_sat`       | Shockley saturation current      | 1 × 10⁻¹² A |
| `n_diode`     | ideality factor (1.0–2.0)        | 1.05 |
| `tau_carrier` | minority carrier lifetime        | 10 ns |
| `dn_dv_inj`   | exponential prefactor from `delta_neff_vs_v` fit on V ≥ 0 | 1 × 10⁻⁴ /V |
| `da_dv_inj`   | corresponding loss prefactor     | ~150 Np/m per V (positive-going) |

**Out of scope:** reverse bias.

> *Note*: L2b is a separate device class because injection vs depletion are
> physically distinct enough that one mediocre fit for both regimes is
> worse than two good fits for each.  Schematic designers pick whichever
> matches the bias they actually use.

---

## L3 — `fc_pn_ps_full` (combined depletion + injection + TPA, steady-state) [NEW]

**Use when:** the modulator may visit either bias regime in normal
operation, OR the optical power is high enough that TPA matters.  Most
"design-grade" simulations of high-performance modulators belong here.

**Physics in addition to L2a + L2b:**
- Smooth transition through V_pn = 0 (one model, both regimes).
- Two-photon absorption: α_TPA = β_TPA · (|A|² / A_eff), where |A|² is the
  per-channel SVEA intensity in the device.  Loss term scales with optical
  power on top of the FCA loss.
- TPA-generated free-carrier steady-state:
  N_TPA = (β_TPA · |A|²/A_eff · |A|²) · τ / (ℏω · A_eff) → adds to FCA.
- Static thermal Δn from absorbed power: ΔT_ss = R_th · P_abs;
  Δn_thermal = (dn/dT) · ΔT_ss ≈ 1.86 × 10⁻⁴ · ΔT/K.

**Geometry assumption:** 100 nm junction offset toward N.

**Parameters added beyond L2a + L2b:**

| device parameter | source | typical value |
|---|---|---|
| `beta_TPA`        | `tpa.beta_TPA_m_per_W`           | 7.9 × 10⁻¹² m/W |
| `a_eff_um2`       | `mode.A_eff_um2`                 | 0.126 µm² |
| `dn_dt`           | (constant for Si)                | 1.86 × 10⁻⁴ /K |
| `r_th_K_per_W`    | thermal resistance to substrate (from substrate-PDK if available) | ~10⁴ K/W per µm |

**Out of scope:** dynamic carrier transients (τ-decay of TPA carriers in
transient), pulsed self-heating.  Those belong in L4.

---

## L4 (planned) — `fc_pn_ps_carrier` (dynamic carriers + thermal RC)

**Use when:** doing transient analysis at high power where carrier
relaxation (τ ≈ 1 ns) and thermal time constants (τ_th ≈ 1 µs) matter.

**Physics in addition to L3:**
- dN/dt = G − N/τ_carrier rate equation, integrator-managed state row.
- Thermal RC (separately modelled via `fc_thermal_ps_rc`'s pattern).
- Self-pulsing / Q-switching regime accessible.

**Out of scope (everything is in L4):** quantum noise, RIN, ASE.  Those
belong to laser models, not modulators.

---

## Tier summary table

|                                  | L1 `fc_pn_ps` | L2a `_cap` | L2b `_inj` | L3 `_full` | L4 `_carrier` |
|----------------------------------|:-:|:-:|:-:|:-:|:-:|
| Linear EO Δn ∝ V                  | ✓ | ✓ | ✓ | ✓ | ✓ |
| Constant α                       | ✓ |   |   |   |   |
| C_j(V) (depletion)                |   | ✓ |   | ✓ | ✓ |
| C_d(V) (diffusion)                |   |   | ✓ | ✓ | ✓ |
| Bias-dependent α(V) — reverse     |   | ✓ |   | ✓ | ✓ |
| Bias-dependent α(V) — forward     |   |   | ✓ | ✓ | ✓ |
| Shockley I-V                     |   | ✓ (reverse) | ✓ (forward) | ✓ | ✓ |
| Nonlinear Δn_eff(V) — reverse     |   | ✓ |   | ✓ | ✓ |
| Nonlinear Δn_eff(V) — forward     |   |   | ✓ | ✓ | ✓ |
| Bidirectional bias regime         | ✓ |   |   | ✓ | ✓ |
| TPA + TPA-induced FCA             |   |   |   | ✓ | ✓ |
| Self-heating Δn                   |   |   |   | ✓ (static) | ✓ (dynamic) |
| Dynamic carrier rate equation     |   |   |   |   | ✓ |
| Junction-offset assumed            | symmetric | 100 nm | 100 nm | 100 nm | 100 nm |
| Implemented today                 | ✓ | ✓ | ✗ | ✗ | ✗ |

---

## What the script actually computes vs what's left in the model

`pn_modulator.py` already does:
- Femwell mode solve → exact n_eff
- Gaussian approximation for |E|² → A_eff
- 1-D abrupt-junction depletion → W(V), C_j(V)
- 2-D depletion mask on a fine x,y grid → carrier-density perturbation
- Soref-Bennett at 1550 nm → Δn(x,y,V), Δα(x,y,V)
- First-order perturbation theory overlap → Δn_eff(V), Δα(V)
- Low-injection minority-carrier injection at V > 0
- TPA loss coefficient from textbook β_TPA / extracted A_eff

What the script doesn't do (and the models will need to in L3/L4):
- High-injection regime — model breaks down for V_pn near V_bi
- 2-D Poisson coupled to drift-diffusion — would give more accurate
  forward-bias carrier profile (femwell can drive a real TCAD solver but
  it's overkill for this tier of model)
- Thermal solver — R_th and the static ΔT response need a separate
  thermo-physical sim (or a measurement)
- Carrier lifetime measurement — currently a hand-tuned parameter

For most designs, the script's first-order numbers are within ~30 % of
TCAD, which is good enough to get the simulator's behaviour qualitatively
right while keeping setup time to seconds rather than days.
