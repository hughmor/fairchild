# AWG router and WDM filter response

`fc_awgr` is an N×N cyclic arrayed-waveguide grating router. `fc_mux` /
`fc_demux` gained an optional per-channel spectral response at the same time.
Both share `models/photonic/spectrum.rs`.

---

## 1. What the model is, and what it is not

Each optical channel wire triple `(re, im, λ)` carries a **slowly-varying
complex envelope** at a carrier, so the physical field is
`E(t) = Σ_k A_k(t)·exp(jω_k t)`. A linear device with impulse response `h`
acts on each channel as a baseband convolution:

```
B_k(t) = ∫ h(τ)·A_k(t−τ)·e^{−jω_k τ} dτ = (h_k ⊛ A_k)(t),   H_k(Ω) = H(ω_k + Ω)
```

These devices keep the **zeroth-order** term: `B_k = H(ω_k)·A_k`, a static
complex coefficient evaluated at the carrier. That is not an approximation of
convenience — it is the exact narrowband limit — and the error is bounded:

| situation | relative field error |
|---|---|
| carrier at band centre, modulation bandwidth `B` | `≈ 2·ln2·(B/FWHM)²` → **1.4 %** at `B = FWHM/10` |
| carrier detuned by `Δ` | first-order tilt `≈ 4·ln2·Δ·B/FWHM²` → **7 %** at `Δ = FWHM/4`, `B = FWHM/10` |

**Reproduced exactly:** insertion loss, port non-uniformity, passband detuning
penalty, channel crosstalk (including its coherent accumulation), thermal grid
drift, FSR periodicity.

**Not reproduced, by construction:** sideband shaping / ISI from passband
narrowing, PM→AM conversion off a detuned passband slope, differential group
delay between channels. All three need a rational `H_k(Ω)` integrated as ODE
states in the MNA — see §7.

Rule of thumb: trust it while `B < FWHM/10`. A 100 GHz-grid AWG with a 40 GHz
passband against a ~1 GHz neural-network modulation rate is far inside that. A
25 Gbaud datacom link is not, and that is exactly the regime where the
literature quotes an "AWG bandwidth-narrowing penalty".

---

## 2. The rule the port shapes follow

> **Fields may only be summed inside one channel slot. Never across slots.**

Two carriers do not add as complex envelopes: `|a_j + a_k|²` would contribute a
static `2·Re(a_j a_k*)` term where the physical beat sits at `|f_j − f_k|`
(≈ 100 GHz) and is filtered out by any real photodiode. Summing them injects a
spurious DC offset into every downstream detector.

Adopt **channel slot index ≡ wavelength index** and the router becomes one
complex matrix per slot:

```
out_j[k] = Σ_i t_ij(λ_{i,k}) · in_i[k]
```

Both crosstalk mechanisms fit inside that single form:

- **Wrong-port crosstalk** arrives at the *same* wavelength → same slot →
  coherently sums with the intended signal. Correct, and it dominates AWG
  penalty budgets.
- **Wrong-wavelength crosstalk** lands in *its own* slot and stays separate.

Neither ever mixes carriers. This is why the router is exactly representable —
and why a demux with single-channel outputs is not (§6).

---

## 3. `fc_awgr`

Input port `i` channel `k` → output port `(i + k) mod N`, still in slot `k`.
Terminals are `2·wpc·N²`: all N input ports, then all N output ports, each
carrying its N channels in order.

```spice
.optical_port in0 8
* … in1 … in7, out0 … out7 …
Xr in0 in1 in2 in3 in4 in5 in6 in7  out0 out1 out2 out3 out4 out5 out6 out7
+   fc_awgr df_ghz=100 fwhm_ghz=40 il_db=3 xt_adj_db=-30
```

### Modes

Chosen by which parameters are present, not by a mode string.

- **ideal** (nothing set) — the exact cyclic permutation, lossless. A pure
  slot permutation that does *not* consult the grid, so an off-grid comb still
  routes rather than going silently dark.
- **gauss** (`fwhm_ghz > 0`) — super-Gaussian passbands on a periodic
  frequency grid, floored by the crosstalk spec.
- **table** (`.model … sfile="…"`) — measured N×N spectra, §5.

### Parameters

| name | default | meaning |
|---|---|---|
| `lambda0_nm` | `.options lambda_center_nm` | grid anchor: passband centre of the `(j−i) mod N == 0` pairs |
| `df_ghz` | 100 | channel spacing — a **frequency** grid |
| `fsr_ghz` | `N·df_ghz` | free spectral range; the default *is* the cyclic condition |
| `fwhm_ghz` | — | power FWHM; positive selects gauss mode, 0 stays ideal |
| `shape_p` | 1 | super-Gaussian order: 1 = Gaussian, 2–4 = flat-top |
| `il_db` | 3 (gauss) / 0 (ideal) | peak insertion loss |
| `il_tilt_db` | 0 | extra loss at the outermost channel vs the centre one |
| `xt_adj_db` | −30 | adjacent-channel crosstalk floor |
| `xt_bg_db` | −40 | non-adjacent crosstalk floor |
| `dlambda_dt_pm_per_k` | 0 | thermal grid drift (silica ≈ 11, SOI ≈ 80) |
| `t_nom_k` | 300.15 | reference temperature for the drift |
| `lambda_src` | 0 | which input port's λ tags the outputs mirror |

The grid is in **frequency** because an AWG is: its FSR is set by the arm
path-length increment, `FSR = c/(n_g·ΔL)`, a constant in Hz, and ITU spacings
are defined in GHz.

The **crosstalk floors are not optional realism**. A pure Gaussian tail three
channels out is below −1000 dB; a fabricated AWG floors at −25…−35 dB from
phase errors in the array. `max(gaussian, floor)` reproduces the datasheet
shape with the two numbers datasheets actually quote.

### Output λ tags

Output `(j, k)` mirrors input port `lambda_src`'s slot-`k` λ wire, falling back
to the device's own grid wavelength if that node does not exist. The test is
*structural*, not value-based: the first coefficient build happens on the zero
initial guess where nothing is driven yet, so a value-based test would freeze
every output onto the grid and discard the comb's detuning. Mirroring is what
keeps a detuned laser detuned for whatever resonant device sits downstream;
transmission itself is always evaluated at the actual per-input λ either way.

### Band limit

The FSR periodicity is bounded: the response is zero more than one FSR beyond
either end of the channel grid. Real star couplers roll off after a few FSRs,
and an unbounded wrap is not merely optimistic — it folds Newton's early
iterates (λ ≈ 1e-8 m on the way up from zero) straight back onto a passband,
so coefficients thrash between iterations until the line search collapses to
minimum steps and the solve never converges.

---

## 4. Solver tolerance — fixed, and why the old workaround is gone

**Nothing to do. Run an `fc_awgr` on default options.** Earlier revisions of
this document told you to put `.options vntol=1e-14 reltol=1e-12` in any deck
containing a router. That is obsolete; the deck no longer has to know about the
solver.

The problem was real. λ wires carry ~1.55e-6, and SPICE's default absolute
*voltage* tolerance `vntol = 1e-6` is the same order as the entire quantity, so
Newton's step test was satisfied while λ was still ~10 pm off — a real detuning
for a 40 GHz passband. The router then reported the transmission for the wrong
wavelength, with no error. Measured on the 8×8 example: `V(in0_0_wl)` settled
at 1.53783 µm instead of 1.55 µm and the routed output read 0 instead of 1.109.
At N ≤ 5 the first Newton step lands accurately enough that it never showed,
which is precisely what made it a trap.

This was never specific to this device — every λ-reading photonic model has the
same exposure, ring resonances included (12 pm is comparable to the 13.3 pm/V
depletion tuning measured on giona).

The fix is per-row tolerances in the solver: λ rows carry their own absolute
tolerance `lambdatol` (default 1e-13 m = 0.1 pm) and **no relative term**, since
`reltol·|λ|` is scale-invariant and would still permit 1.55 nm. See
`crates/fairchild-core/src/tolerance.rs` for the full argument, including why
respelling λ in µm does not solve it. `an_eight_by_eight_router_is_exact_on_default_tolerances`
is the regression.

---

## 5. Measured data

Real data is an N×N grid of transmission spectra on a common λ grid — from a
tunable laser + power meter sweep (magnitude only) or a swept interferometric
instrument (magnitude and phase). The reader takes CSV, because that is what
dumps out of numpy in one line:

```
wavelength_nm,t_0_0_db,t_0_1_db,t_1_0_db,t_1_1_db,t_0_0_deg
1549.0,-3.1,-31.4,-30.9,-3.2,0.0
1551.0,-3.4,-31.1,-31.2,-3.3,0.0
```

`t_<in>_<out>_db` is required per pair; `t_<in>_<out>_deg` is optional and
defaults to zero phase. Rows may be unordered. Missing pairs read as dark, so a
partially measured device still runs. Interpolation is **linear in dB and in
unwrapped degrees** — a spline rings on the steep skirts and can overshoot into
negative power, which surfaces as a NaN three devices downstream. Values are
held for 10 % of the span past each end and dark beyond that.

The path is a string, and X-line instance params are numeric only, so a
measured router must come in through a `.model` card:

```spice
.model awg8 fc_awgr sfile="awgr8.csv"
Xr in0 … in7 out0 … out7 awg8
```

---

## 6. `fc_mux` / `fc_demux`

Both still default to a **lossless identity route** — existing netlists are
bit-for-bit unaffected. Setting any of `il_db`, `lambda0_nm`, `df_ghz`,
`fwhm_ghz`, `shape_p`, `dlambda_dt_pm_per_k`, `t_nom_k` turns on a **diagonal**
filter: channel `k` scaled by its own passband, evaluated at the wavelength
that channel actually carries. λ labels are never attenuated.

```spice
Xmux bus c0 c1 c2 c3 fc_mux il_db=3 fwhm_ghz=40 df_ghz=100
```

Diagonal is the whole story for a mux — N inputs land in N distinct slots, so
cross-channel leakage has nowhere to go. It is deliberately **incomplete for a
demux**: a real demux does leak channel `k` into output port `j ≠ k`, but an
`fc_demux` output port is a single channel, so representing that leak would
mean summing two carriers into one envelope — forbidden by §2.

For a demux *with* crosstalk, use `fc_awgr` with `N−1` input ports left dark.
That is not a workaround: a demux **is** an AWG with one input used, a mux is
one with one output used, and a router is one with all of them used — the same
star-coupler-plus-array silicon. Dark ports contribute nothing regardless of
their λ tags. The router's output ports are N-channel buses, which is exactly
the somewhere the leakage needs to live.

---

## 7. Not implemented

Deliberate omissions, so nobody has to rediscover them:

- **Synthetic crosstalk phase, and the crosstalk-penalty distribution.**
  Analytic modes produce **purely real** transmission, so every crosstalk term
  adds in phase. That is the *pessimistic* bound — a conservative
  over-estimate, which is the right default — but it is not the physical
  spread. Real crosstalk arrives with essentially random phase, and
  magnitude-only measurements leave that phase undetermined, so the honest use
  of such data is to sample it: a `phase_mode` / `xt_phase_seed` pair
  generating deterministic pseudo-random per-pair phases, swept over ~20 seeds
  to get a penalty distribution rather than one number. **Not built.** Table
  mode does honour *measured* phase from `t_<i>_<j>_deg` columns.
  Also unbuilt, and the principled alternative for magnitude-only data:
  minimum-phase reconstruction via a Hilbert transform of `ln|H|`.
- **Sideband fidelity (tier 2).** Per §1. The shape is known: a one-pole
  baseband filter per coupling, `dB/dt = (−γ + jδ)·B + γ·t·A`, two extra MNA
  rows each, state owned by the device via `commit_timestep` (the `DelayLine`
  precedent, not `ReactiveBranchSpec`). Full fidelity is `N³` states; the
  pragmatic version puts the pole only on the intended coupling and leaves
  crosstalk static, for `N²`.
- **Bidirectional propagation** (`wpc = 5`). Backward fields need the
  transposed routing; `setup_instance` rejects it rather than leaving the
  backward wires undriven.
- **Latency** (`tau_s`). An AWG's few-ps array transit is far below any
  timestep this simulator runs at.

---

## 8. Cost

An N×N router owns `9N²` MNA rows (`2·wpc·N²` port wires + `3N²` branch rows).
That is inherent to the wire-per-channel representation; prefer KLU for N ≥ 8.

The coefficient table is `N³` transcendentals, rebuilt only when an input λ
changes — with CW lasers that is once. Stamping replays a precompiled list.

`Pattern::build` otherwise takes a **clique** over every row a device owns,
which for this device is `O(N⁴)` against a true `O(N³)` footprint. `fc_awgr`
implements `Device::stamp_pairs` to declare the real one:

| N | rows | declared nnz | clique nnz | ratio |
|---|---|---|---|---|
| 2 | 37 | 103 | 1 299 | 13× |
| 4 | 145 | 563 | 20 739 | 37× |
| 8 | 577 | 3 523 | 331 779 | **94×** |
| 16 | 2 305 | 24 323 | 5 308 419 | **218×** |

`stamp_pairs` returns `None` by default, so every other device keeps the
conservative clique. Any device implementing it must keep it in sync with its
stamping — cells outside the declared footprint are dropped by the sparse
solve (debug builds catch this via `MnaMatrix::debug_assert_covers`).

---

## 9. See also

- `crates/fairchild-core/src/models/photonic/awgr.rs` — device
- `crates/fairchild-core/src/models/photonic/spectrum.rs` — shared response model
- `crates/fairchild-core/tests/native_awgr.rs` — including the equivalence test
  against an `fc_demux` + `fc_mux` permutation built from primitives
- `examples/photonic/native_awgr_router.py` — routing map, crosstalk matrix,
  passband sweep (`--selftest` asserts all three)
- `_notes/awgr_design.md` — the design pass this came from
