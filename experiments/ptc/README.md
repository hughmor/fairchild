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
- **Experiment E** (`accumulation.py`): per-MAC ENOB is flat in `L` to 0.17 b
  over `L = 1…32`, output-referred rises at **0.478 b per octave** against the
  predicted 0.5. Confirms the `1/√L` converter amortization of eq:adc_power.

### One discrepancy worth chasing

In the RIN-dominated corner the sim lands **0.3–0.4 b below** the closed form.
About 0.28 b of that is accounted for by a definition: eq:coefficients writes
`⟨i_RIN²⟩ = r I_ph² B` against the *average* photocurrent, but RIN is
multiplicative on the *instantaneous* power, and two cascaded modulators inflate
`⟨P²⟩/⟨P⟩²` by `(1 + m²⟨u²⟩)² = 1.47` for uniform operands at `m = 0.8` — 1.7 dB,
i.e. 0.28 b off the ceiling. Substituting that factor into `A_2` drops the worst
disagreement from 0.375 b to 0.164 b, the remainder being sample statistics
(640 MACs) and the boxcar-bandwidth approximation.

If it holds up, `A_2 = r f_B (1 + m²σ_op²)² + k_xt χ` — a small correction to
`SNR_∞`, in the direction of *less* headroom. Shot noise is unaffected (its mean
is untouched by zero-mean modulation). Worth confirming with a dedicated run at
high `I_ph` and swept `m` before it goes anywhere near the manuscript.

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

### 2. Transient random-noise source — worked around

There is no `TRNOISE`-style primitive, so noise here is pre-sampled and injected
as `PWL` sources: one independent unit-variance node per receiver carrying TIA +
shot through a `B`-element product with the instantaneous optical power, and one
node per *wavelength* for RIN, shared by every receiver reading that laser
because RIN is a property of the source. This is exact but the netlist grows as
`(#noise nodes × duration × 2 f_B)` — already 0.5 MB at `N=8`, 10 blocks. A
`TRNOISE(rms, dt)` waveform would be ~60 lines in the waveform module and make
noise decks O(1) in size. Worth building before going past `N = 8`.

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
