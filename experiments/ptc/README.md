# ptc — hypermultiplexed photonic tensor core

Circuit-level transient simulation of the cascaded-MZM + AWGR architecture
analysed in [`ptc-scaling-analysis/manuscript.md`](https://github.com/Shastri-Lab/hypermultiplexed-ptc-scaling-analysis),
built to test that manuscript's closed-form models against a solved circuit.
Scaling parameters follow the manuscript: `N` wavelengths, `S` spatial fan-out,
`K = N`, `L` accumulated symbols, symbol rate `f_B`. Digital interfaces and the
analog readout mux are out of scope; the integrator is treated as ideal.

```bash
python3 experiments/ptc/ptc.py --selftest              # correctness
python3 experiments/ptc/ptc.py -N 4 -o core4.sp        # just the netlist
python3 experiments/ptc/enob_sweep.py --png sweep.png  # ENOB vs P_det
```

## What the netlist is

`ptc.py::build` emits, for a given `(N, S, L, f_B)` and a block of operand data:

| stage | devices | count |
|---|---|---|
| comb source | `fc_cw_laser` → `fc_mux` | `N` + 1 |
| broadcast | `fc_splitter` chain | `N-1` |
| weight bank | `fc_mzm` on an `N`-channel bundle | `N` |
| router | `fc_awgr` | 1 |
| spatial fan-out | `fc_splitter` chain | `N(S-1)` |
| input bank | `fc_mzm` | `SN` |
| receivers | `fc_demux` + `fc_photodetector` + `C` | `SN` + `2N²S` |

Receiver `(j, s, k)` accumulates weight row `i = (j - k) mod N` against input
column `(j, s)`, so an `L`-symbol block yields `Y[i,s,j] = Σ_l W[i,l] X[l,s,j]`.
`readout()` applies that index map; `decode()` subtracts the deterministic
linear MZM terms of eq:cascaded_mzm using a zero-drive calibration run.

Three things that are not obvious and will bite:

- **`uic=True` is mandatory.** The integrator node sees only `C_int` and
  `r_shunt`; a DC operating point puts `I_ph × 1e15 Ω` on it and Newton never
  converges. The physical initial state is a discharged integrator.
- **`method="tr"`.** The observable is the integral of a piecewise-linear
  photocurrent. Trapezoidal is exact on it; BE and GEAR are first-order and
  leak charge at every symbol edge.
- **`tr_frac × oversample ≥ 4`**, asserted in `build`. An unresolved symbol edge
  is indistinguishable from a slower one: the solver smears the transition over
  a timestep and the MAC error stops tracking `tr` and starts tracking `dt`.

## Verified so far

- Ideal router + predistorted drive reproduces the exact MAC. The only residual
  is the symbol edge: two operands ramping together across a boundary integrate
  to `(tr/6)·ΔW·ΔX` less than the rectangular product, **linear in `tr/f_B`**
  and independent of the timestep. That is ISI, which the manuscript defers.
- Noise (TIA + shot + RIN) matches eq:master_snr / eq:coefficients to **0.13 b**
  from the TIA floor, through the shot knee at `I* = A_0/A_1 ≈ 1 mA`, up to the
  RIN ceiling (`enob_sweep.py`, `results/enob_sweep.png`).
- AWGR port isolation degrades the MAC as expected when tightened to −15 dB.
- **Experiment E** (`accumulation.py`): per-MAC ENOB is flat in `L` to 0.08 b
  over `L = 1…32`, output-referred rises at **0.485 b per octave** against the
  predicted 0.5 (rerun on native noise sources). Confirms the `1/√L` converter amortization of eq:adc_power.

### A correction to A₂ that the sim keeps insisting on

eq:coefficients writes `⟨i_RIN²⟩ = r I_ph² B` against the **average**
photocurrent. RIN is multiplicative on the **instantaneous** power, and two
cascaded modulators raise `⟨P²⟩/⟨P⟩²` by `(1 + m²⟨u²⟩)²` — 1.47 for uniform
operands at `m = 0.8`, i.e. 1.7 dB or 0.28 b off the ceiling.

Measured three ways, and it survived all of them:

1. Unmodulated, RIN-only: the injected noise is `√(r I_ph² f_B / 2)` to **0.5 %**,
   so the generator itself is exactly the paper's expression.
2. Modulators running, `m` swept 0 → 1 at RIN-dominant power: the noise grows as
   `√(⟨P²⟩/⟨P⟩²)` to ~1 %, matching `(1 + m²/3)` per stage.
3. Full deck at `I_ph = 4.6 mA`, 1920 MACs over 3 seeds: **8.03 b** measured
   against **8.37 b** for the paper as written and **8.14 b** with the factor in.

So `A₂ = r f_B (1 + m²σ_op²)² + k_xt χ` — a small correction, in the direction of
*less* headroom. Shot noise is unaffected (its mean survives zero-mean
modulation untouched). It also matters more than 0.28 b suggests, because `A₂`
sets `SNR_∞`, which no amount of power moves.

Residual after the correction: **0.11 b**, still slightly systematic. Not yet
explained; candidates are the boxcar-bandwidth approximation on a
signal-correlated noise term, and shot noise still being 15 % of the total at
that current. Small enough not to block anything, big enough not to call closed.

## What fairchild is missing

### 1. Optical amplifier — the one real gap

No SOA/EDFA. Needed for the manuscript's central link result
(fig:link_limited: signal–spontaneous beat noise dominating the whole feasible
interior, and the saturation-forced amplification law). A `fc_soa` needs:

- gain `G`, noise figure `F`;
- **aggregate** output saturation `n_λ P_ch G ≤ P_sat` — a coupling *across*
  channels of one bundle, which no photonic device in the tree currently has;
- ASE injection as an additive complex field per channel slot. Signal–spontaneous
  beat then falls out of the photodiode for free, because ASE inside a slot
  legitimately sums with that slot's carrier. Spontaneous–spontaneous does not:
  it needs the optical bandwidth `B_o`, which the envelope does not carry, so it
  has to stay a receiver-side term.

Without it, the reachable `I_ph` stops at what a 13 dBm source can push through
the passive link — well below the shot-noise knee `I* = A_0/A_1 ≈ 1 mA`. The
sweep figure marks that point.

Interim: gain as a negative-loss attenuator plus an equivalent receiver-side ASE
current noise. That gets the ASE *noise* right at a fixed operating point and the
*placement* physics wrong, so it cannot reproduce fig:link_limited.

### 2. Transient noise — solved upstream

`.options trannoise=1 noiseseed=<n>` (master, `da8ac20`) injects every generator
the `.noise` table lists as a per-timestep random current. `fc_photodetector`
shots `2q(I_ph + I_dark)` off its *current* photocurrent, so shot noise tracks
the modulation; `fc_cw_laser rin_db_hz=` injects RIN at the source, so it
propagates through both modulator banks and arrives at every receiver of
wavelength `k` perfectly correlated — the physical behaviour of a shared laser,
for free. This replaced a hand-rolled PWL construction here and the deck got
smaller, faster, and more correct.

Two constraints it imposes, both fine:

- **Fixed step required.** `variable_step=1` is refused, not approximated.
  This deck already ran fixed-step for the symbol edges.
- **Noise bandwidth is the timestep's Nyquist.** Harmless here: the integrator
  band-limits, so the block charge is step-independent for any `h ≪ 1/f_B`.

The one source with no device is the **TIA's input-referred current noise**.
There is no TIA model and no element spelling "a current source of PSD `S`". A
resistor of `4kT/S_i²` = 51 Ω has the right spectrum but a 25 fs time constant
across the integrator. `_tia_noise` therefore picks a convenient `R`, and
buffers its `4kT·R` volts through a `B`-element transconductance sized to land
the wanted current PSD: the resistor sets the spectrum, the `B`-element sets the
amplitude, nothing loads the integrating node. It works, but a generic
`noise PSD=<S>` two-terminal element would be the honest fix, and it is the
natural home for ASE too.

**Size `C_int` to the photocurrent.** At 500 mW the default 500 fF ends the run
on a kilovolt node, where `reltol` alone exceeds the noise being measured and
ENOB collapses from 8.3 b to 2.1 b. Nothing physical — but it looks like a
physics result, so it is worth stating.

### 3. Spectral crosstalk — structurally out of scope

**This is the one to be clear about.** The manuscript's headline result
(eq:wall, the spectral packing wall) is a beat between a victim's carrier and
the *modulation sidebands* of its spectral neighbours. fairchild carries one
complex envelope per channel slot and forbids summing fields across slots — the
reason `fc_demux` has no crosstalk parameter at all. There are no sidebands, and
router transmission is a static coefficient at each carrier (the exact
narrowband limit). So `χ` and `k_xt` cannot *emerge* from a simulation here;
they can only be injected as an equivalent receiver-side noise, which restates
the theory rather than testing it.

Reproducing it honestly would mean representing one WDM channel as `M`
sub-carrier bins and giving the photodiode the full `|Σ E|²` including cross-bin
beats inside the electrical band. That is a different detector model and a
different governing rule — a project, not a device.

### 4. AWGR *spatial* crosstalk — in scope, and the opportunity

Leaked port images arrive at the victim's **own** wavelength, so they land in the
same slot and coherently sum. `fc_awgr` already models this
(`xt_adj_db`/`xt_bg_db`). The manuscript *defers* this term — it appears only as
the bound `ε_agg ≈ 2ε_adj + (N−3)ε_bg`, with the note that it "should be promoted
from bound to model term once device data are available". The simulator supplies
exactly that promotion, with the operands' real statistics instead of a
worst-case in-phase sum.

Caveat: `fc_awgr`'s analytic modes are purely real, so every leak still adds in
phase — the sim currently reproduces the pessimistic bound, just with correct
data statistics. The deferred `phase_mode` feature (random per-pair array phases)
is what turns it into a distribution. **That is the one fairchild feature worth
building for this paper**, and it is much smaller than an SOA.

### 5. Smaller gaps

`fc_mzm` accepts `f_c` but ignores it (no EO bandwidth limit); `fc_photodetector`
has no junction capacitance or transit time; there is no driver-nonlinearity
model; the laser has no linewidth or phase noise. None block the experiments
below. The **integrator is not missing** — it is a capacitor on the PD anode.

## Recommended experiments

Ordered by value to the manuscript. A–E need no new devices.

**A. AWGR spatial crosstalk → ENOB(N, ε).** Sweep `N ∈ {4,8,16}` and
`ε_adj ∈ [−40,−20] dB` with every other noise source off; measure ENOB from the
MAC residual. Deliverable: measured `ε_agg` against the `2ε_adj + (N−3)ε_bg`
bound, and how much the bound overstates for random operands. Promotes a
deferred effect to a model term. **Start here.**

**B. Modulator nonlinearity → what `m = 0.8` is worth.** Run `predistort=False`
so the MZM's real `cos` transfer compresses, and sweep `m`. Deliverable: ENOB vs
`m` showing the interior optimum between more signal (`m⁴` in the numerator) and
more distortion — i.e. whether 0.8 is the right baseline. Second deferred effect
priced.

**C. ISI vs rise time.** Nearly free: it falls out of the runs above. The
residual is already characterised as `(tr/6T)·ΔW·ΔX` per edge, linear in `tr/T`.
Third deferred effect priced.

**D. Validation of eq:master_snr.** Done — `enob_sweep.py`, 0.13 b agreement.
Supplementary-figure material: the receiver model is measured, not asserted.

**E. Accumulation bookkeeping.** Sweep `L ∈ {2,4,8,16}` and confirm per-MAC ENOB
is `L`-independent while output-referred gains `½log₂L`. That is the claim behind
eq:adc_power's `1/√L` amortization, which the manuscript flags for triple
checking. Cheap, and it settles a flagged equation.

**F. Link budget and ASE** (fig:link_limited: achievable ENOB vs scale under
saturation-forced placement). Needs `fc_soa` — the only experiment here that
requires building a device.

**Not recommended: spectral crosstalk.** See gap 3.

## Cost

`N=4, S=1`, 10 blocks, noise on: ~8 s. `N=8, S=1`: ~90 s. `N=8, S=2`: ~4 min.
Runtime is dominated by `oversample`, which the `tr ≥ 4 dt` rule ties to
`tr_frac`. When the ISI term is not the object of study, difference against the
same deck with noise off (as `enob_sweep.py` does) — it cancels exactly and lets
`tr_frac` be relaxed.
