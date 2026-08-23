# mpc — model-predictive control on the network

In progress. The method is de Lima et al., *Machine Learning With Neuromorphic
Photonics*, JLT 37(5) 2019, Appendix A: map MPC to a QP, solve the QP with a
continuous-time RNN whose fixed point is the QP solution, and run that RNN on
the network.

The mapping, in this chip's variables:

| CT-RNN | here |
|---|---|
| weight matrix `(Id − P)/α` | `W = (Id − P/s)/G`, `G` measured per neuron |
| external input `−q` | **optical**, on dedicated wavelengths, individually weighted |
| DC bias `P·y₀` | the row-sum rule, `../rnn_characterization/rnn_math.md` §3 |
| neuron output `y` | `I_D`, the modulator junction current |
| `τ` | `τ_c ≈ 10 ns` |

**Planned design point.** 1-D double integrator, regulation, `Hu = 6`:
6 neurons + 2 weighted optical input channels on the unchanged 8×8 hardware
(48 of 64 rings). `q = 2ΘᵀΨ x_k` has rank 2 for a 2-state plant, which is what
forces two individually-weighted inputs — an unweighted 9th wavelength can only
deliver a term common to every neuron, and neither column of `ΘᵀΨ` is uniform.

At `dt = 0.2 s`, `Hp = 16`, `λ_u = 0.3`: `cond(P) = 28`, settling **0.28 µs**
against a 0.2 s plant step, recurrent `max|GW| = 1.18` so `G ≥ 1.18` against
1.8–3.7 measured trimmed. Input weight columns scale independently, so each
normalises to ±1 on-chip with the scale on the off-chip modulator drive.

Weight precision is the binding constraint, not gain: 0.5 % of full scale gives
a 6 % solution error, 1 % gives 12 %, 5 % gives 68 %.

Two decks, as everywhere else here: ideal weights first to prove the loop
closes, then the same topology on the real rings.

## Why the state variable is the junction current

Worth writing down before the next attempt, because it is the one part of the
mapping that is settled. The hardware's fixed point (`rnn_math.md` eq 1) is

    I_i = b_i + eta*kappa * sum_j W_ij * T_j(I_j)

which puts the nonlinearity BEFORE the weighting, not after — so at first glance
it is not the paper's Eq. 8. Substituting `z_i = T_i(I_i)` turns it into
`z_i = T_i(b_i + sum_j W'_ij z_j)`, which is Eq. 8 exactly. So the network is a
CT-RNN, and linearising about the rest point,

    (Id - G*W) dy = db,   dy = dI_D,   G = eta*kappa*T'(I*)

so `Id - G*W = P/s` is the mapping and `dy` — the modulator junction current,
small-signal about wherever the network rests — is the QP variable.

Two consequences for the driver, once the gain measurement is trustworthy:

- **The input transfer should be measured, not derived.** Drive one input
  channel with the recurrent weights off and read `dy`; that column of B is
  known without assembling `eta`, `kappa`, `a` and `R`, each of which is its own
  calibration with its own error and only the product matters.
- **The rest point can be taken as found.** The paper's `y0` offset and the
  row-sum bias rule exist to PUT the rest point somewhere convenient; measuring
  where it lands and working relative to it is the same thing with one fewer
  stage that can be wrong.

## Status

| step | state |
|---|---|
| `build_mpc_deck.py` — both decks | done; both converge |
| `characterize.py` — operating point map | done; see below |
| `gain_matrix.py` — open-loop neuron-to-neuron gain | done; structure is clear, magnitudes are not reproducible |
| per-ring trim to equalise gain | **attempted and NOT achieved — see "the measurement problem"** |
| bias solve (row-sum rule) | not started |
| program `W`, check the fixed point against a digital QP | not started |
| close the plant loop | not started |

### Operating point

30 mW/channel, `PDB = −8 V`, `Iht ≈ 3.12 mA`, where `G = 2.49` on neuron 1
against the 1.18 the weights need.

**More laser power is not simply more gain.** The corrected modulator model has
real self-heating (`r_th = 3139.9`), and the absorbed light detunes the ring off
the resonance the modulation depends on:

| P/channel | thru/P on its own channel | ring ΔT | resonance shift |
|---|---|---|---|
| 0.03 mW | 1.8 % — the full notch | 0.03 K | 1.7 pm |
| 1 mW | 21 % | 0.65 K | 45 pm |
| 30 mW | 73 % — notch washed out | 3.4 K | **236 pm** |

The linewidth is ~167 pm, so 30 mW sits 1.4 linewidths off. `G` still rises with
power — 3× the power from 10 to 30 mW buys 2.8× the gain — but sublinearly, and
10 mW gives only `G = 0.89`, short of the 1.18 needed. The old `fc_pn_th_ps`
card had `r_th = 0` and was structurally blind to this, so it is a prediction of
the corrected model rather than a measurement. It is worth testing against the
chip, because it says the WTA's answer to low gain — more power — was pushing
against a wall it could not see.

### The measurement problem

Measuring the loop gain on this deck has defeated three estimators, and nothing
downstream should be built until one of them is pinned against a known answer.

1. **Closed loop** (`rnn_math.md` (11), `G = (1 - S0/S(w))/w`). Fails two ways
   at once: the probe weight has to keep `w*G` under 1 or the self-coupled
   neuron is at its own bistability threshold and `S(w)` diverges and flips
   sign; and `reltol = 1e-3` on a 977 mV node is a millivolt of slop against a
   1.4 mV response. Raising the bias step to 200 mV fixed the worst of it —
   neuron 5 went from -4255 to a consistent 5.6 — but four of six neurons still
   disagreed between probe weights by more than 15 %. Tightening `reltol` is not
   the answer: at 1e-10 the solver chases the gmin-floored dark optical nets
   (`ain1_re_7` asks to move 5e47 V) and stops converging.
2. **Open loop, by column** (`gain_matrix.py`). Well conditioned in principle —
   drive neuron j, read everyone's response, no small differences — and it gave
   the one structurally believable result here: every ROW of `G_ij` is identical,
   so the six neuron circuits are uniform and all the spread is in how hard each
   RING drives the bus (0.003, 0.046, -0.424, 1.425, 1.954, 0.018).
3. **The same measurement inside a heater sweep** (`trim.py`). Returns ~0.01 for
   every ring at every heater current, including the configuration where (2)
   read 1.4 and 2.0, with occasional spikes to 11. Two measurements of one
   quantity disagreeing by two orders of magnitude means the measurement is
   wrong.

Suspects, none yet eliminated: `source_gain` puts the probe on the driven neuron
as well (`W_jj`), so its reference response carries self-feedback; the divided
responses are again ~1 mV on 977 mV; and the rings couple hard enough through
the shared bus that a one-at-a-time sweep may not isolate what it intends.

**The next step is not another sweep.** It is to build a case with an arithmetic
answer — an ideal `fc_optical_2x2` at a programmed weight, where the expected
photocurrent is known — and make the estimator reproduce it before it is trusted
to set anything.

### The thing that motivated the trim

At a single shared trim, per-neuron `G` spans **1.58 to 5.60**. Every neuron
clears the threshold, so the loop will close, but the spread is the problem:
row `i` of `W` scales by `1/G_i`, so the strong neurons use a fifth of the
weight range and their weights land proportionally coarser. Against the 0.5 %
precision budget above, that is a 3.5× penalty on the best-off row.

Every ring got ring 1's `Iht`. Trimming each ring individually to equalise `G`
is what makes the precision budget reachable on every row rather than one.
