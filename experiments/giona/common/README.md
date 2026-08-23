# common — what more than one experiment needs

Nothing here is an experiment. If a file in this directory is only used by one
experiment, it is in the wrong place.

| file | what |
|---|---|
| `ringfit.py` | dataset → observables (`load_sweep`, `extract_data`), netlist assembly, the ring wavelength sweep, the staged fitter, and the fit plotting |
| `rnn_drive.py` | the network driver: bias solve, `.op`, transient, `.param` substitution over the RNN deck |
| `giona_pn_th_ps.inc` | the device card — **the** deliverable of `pn_modulator_fit`, and an input to everything downstream |
| `netlists/` | decks no experiment owns: legacy hand-written captures of the chip kept for reference |

`ringfit.py` keeps its own CLI, which runs the staged fit. That run belongs to
`pn_modulator_fit`, and its `RESULTS` points there. A caller wanting figures
somewhere else assigns `ringfit.RESULTS` before calling.

The parser has no `.include` for the card: concatenate `giona_pn_th_ps.inc`
ahead of a netlist (the fit scripts do card + netlist via `load_str`).
