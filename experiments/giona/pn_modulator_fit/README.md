# pn_modulator_fit — the `fc_pn_th_ps` card

Produces [`../common/giona_pn_th_ps.inc`](../common/giona_pn_th_ps.inc), whose
defaults reproduce the measured device. Everything downstream reads that card.

| script | dataset | fits | result |
|---|---|---|---|
| `../common/ringfit.py` | May sparse joint IV+spectra (neuron2 capture → neuron7 params) | staged: passive → thermal → EO → injection | `results/giona_neuron7_pn_th_ps_full_fit.json` |
| `fit_may_injection.py` | ↑ same, forward-bias slice | `dn_di`, `da_di` with the diode pinned from the IV alone | `results/giona_dn_di_fit.json` |
| `fit_jul_neuron3.py` | Jul 2026 dense 50×100 HC×JV sweep (3.6 GB) | passive + thermal + linear EO from tracked notch positions | `results/giona_neuron3_pn_th_ps_fit.json` |
| `fit_transient.py` | (machinery + synthetic recovery selftest) | time-domain fitting vs an AWG-drive/scope-PD capture | — |
| `expt_forward_mismatch.py` | synthetic | model-**form** adequacy: linear vs full PN | 4.6× residual floor gap |

```bash
.venv/bin/python experiments/giona/common/ringfit.py --model fc_pn_th_ps_full
.venv/bin/python experiments/giona/pn_modulator_fit/fit_may_injection.py
.venv/bin/python experiments/giona/pn_modulator_fit/fit_jul_neuron3.py --cache
```

The card is a consolidation across datasets, not the output of any single fit —
`p_pi_th` is from July neuron 3, the EO and passive terms from May neuron 7, the
junction capacitance from FEM. Its header says which number came from where.

`results/mod_iv_*.png` are orphans: the scripts that drew them were deleted
before this directory was organised.
