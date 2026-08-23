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

## Status

| step | state |
|---|---|
| `build_mpc_deck.py` — both decks | done; both converge |
| `characterize.py` — operating point, per-neuron `G` | done; see below |
| per-ring trim to equalise `G` | **next** |
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

### The thing to fix next

At a single shared trim, per-neuron `G` spans **1.58 to 5.60**. Every neuron
clears the threshold, so the loop will close, but the spread is the problem:
row `i` of `W` scales by `1/G_i`, so the strong neurons use a fifth of the
weight range and their weights land proportionally coarser. Against the 0.5 %
precision budget above, that is a 3.5× penalty on the best-off row.

Every ring got ring 1's `Iht`. Trimming each ring individually to equalise `G`
is what makes the precision budget reachable on every row rather than one.
