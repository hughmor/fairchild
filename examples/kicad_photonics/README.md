# fairchild × KiCad

Draw a photonic circuit in KiCad's schematic editor, export it as a SPICE
netlist, simulate it in fairchild. See
[`kicad_integration.md`](../../kicad_integration.md) at the repo root for the
full setup guide.

## The symbol library

**`fairchild_photonics.kicad_sym`** — one symbol per native `fc_*` photonic
device. Optical (bundle) ports are a single pin per port; the post-processor
expands each into its `(re, im, λ)` wires at simulate time.
`sym-lib-table` registers it for the example project via `${KIPRJMOD}`.

Most symbols were drawn by hand. Newer ones are generated:

```bash
python3 scripts/gen_kicad_symbols.py          # write the generated symbols
python3 scripts/gen_kicad_symbols.py --check  # exit 1 if the library is stale
```

The generator only touches names listed in its `SYMBOLS` table and replaces
those wholesale, so it is safe to re-run — but do not add a spec for a symbol
whose art someone drew by hand unless you mean to lose it. What the table
really encodes is **pin order**, which has to match the positional net order on
the `X…` line. Keep it in step with `PORT_SCHEMA` in
`scripts/kicad_to_fairchild.py` and `register_native_photonics` in
`crates/fairchild-core/src/device_registry.rs`.

## Examples

**`mrm_single_channel.kicad_sch`** (+ `.kicad_pro`) — the worked example.
A single-channel add-drop micro-ring modulator: CW laser → grating coupler →
bus → ring → grating couplers → through/drop photodetectors with loads. The
ring itself is the hierarchical sub-sheet **`mrm_ring.kicad_sch`**: two
directional couplers closed by two `fc_pn_th_ps` half-rings (PN junction +
heater), parametrised so the same block can be instanced per wavelength
channel. `mrm_single_channel.cir` is its KiCad SPICE export, committed so you
can run the pipeline without opening KiCad.

Two flat reference exports, used to regression-test the post-processor:

- **`mrm_kicad_test.cir`** — single-channel PN-modulated ring
  (laser → waveguide → coupler ⇄ PN phase shifter → waveguide → PD).
- **`mrm_kicad_test_wdm.cir`** — the two-wavelength version; both channels
  share one bus and one ring, and the parser's per-channel replication does
  the routing.

## Running one

```bash
python3 scripts/kicad_to_fairchild.py examples/kicad_photonics/mrm_single_channel.cir \
    --tran "10n 5u" --method gear --output /tmp/run_mrm.sp
target/release/fairchild -f /tmp/run_mrm.sp --probe "V(V_THRU),V(V_DROP)"
```

Or from the GUI, which does export → transpile → simulate → plot in one shot:

```bash
python3 scripts/fairchild_gui.py examples/kicad_photonics/mrm_single_channel.kicad_sch
```

## Not committed

Transpiler wrapper netlists (`run_*.sp`), simulation outputs, KiCad `*.bak` /
`.history` / `*.kicad_prl`. Chip-specific fitting work lives in
[`experiments/`](../../experiments), not here.
