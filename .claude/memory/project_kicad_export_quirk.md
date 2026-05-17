---
name: project-kicad-export-quirk
description: KiCad's SPICE exporter emits X-element lines in TWO different shapes depending on whether the component has Sim.Params beyond model+type — both formats can appear in a single netlist
metadata:
  type: project
---

KiCad's SPICE exporter (KiCad 8+, Sim.* property namespace) emits X-element lines in two different shapes depending on whether the component has any `Sim.Params` value beyond the bare `type=X model=NAME` pair:

- **With extra params** (e.g. `Sim.Params = V_pi=1.0 f_c=10e9`):
  `MZM1 net0 net1 ... type=X model=fc_mzm V_pi=1.0 f_c=10e9`
  (refdes verbatim, kwargs incl. `type=X` and `model=NAME`)

- **Without extra params** (MUX, DEMUX, splitter — anything with no tunable parameters):
  `XSPL4 net0 net1 net2 fc_splitter`
  (X-prefix added to refdes; model is the rightmost positional token)

**Why:** confirmed by user observation in `examples/kicad_photonics/mrm_kicad_test_wdm.cir` (2026-05-17). The fc_mzm lines have `type=X model=fc_mzm V_pi=…` while the fc_mux line is just `XMUX1 IN_OPT ch0 ch1 fc_mux`. The first transpiler version only handled format (a) and silently fell through on format (b), which broke bundle-width inference and made WDM netlists panic with "got 5 terminals."

**How to apply:** any code that ingests KiCad SPICE exports (today: `scripts/kicad_to_fairchild.py`) must handle both formats. The robust approach is to pre-scan `.model NAME BASE_KIND` directives to build a known-models lookup, then dispatch on either:
- `type=X` kwarg present, OR
- refdes starts with `X` (case-insensitive)

and find the model by `model=NAME` kwarg in format (a) or as the unique positional token matching a known-models entry in format (b). Multi-match positional should warn (ambiguous), zero-match should pass through unchanged (probably a subcircuit / custom OSDI device — fairchild's parser will then either handle it or error with a clear line reference).

Related: [[project-fairchild]] has the KiCad integration status at row 3.7.
