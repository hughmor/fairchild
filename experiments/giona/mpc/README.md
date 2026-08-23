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
