# fairchild × KiCad

The official symbol library + reference SPICE exports for the KiCad
integration. See [`kicad_integration.md`](../../kicad_integration.md) at
the repo root for the full setup guide.

## Files

- **`fairchild_photonics.kicad_sym`** — KiCad symbol library. One symbol
  per native `fc_*` photonic device. Bundle-port pins are single-pin per
  optical port; the post-processor expands them into `(re, im, λ)` wires
  at simulate time.

- **`mrm_kicad_test.cir`** — reference KiCad SPICE export of a single-
  channel PN-modulated micro-ring (laser → waveguide → directional coupler
  ⇄ PN phase shifter (ring) → waveguide → photodetector). Used as a
  regression test for the post-processor.

## Reproducing the reference simulation

```bash
python3 scripts/kicad_to_fairchild.py examples/kicad_photonics/mrm_kicad_test.cir \
    --tran "10n 5u" --method gear --output /tmp/run_mrm.sp
target/release/fairchild -f /tmp/run_mrm.sp --probe "V(Net-_PD1-cathode_),V(/Vmod)"
```

## Files we don't commit

- Wrapper netlists produced by `kicad_to_fairchild.py` (`run_*.sp`).
- Simulation outputs (`*.raw`, `*.csv`, transient PNG plots).
- `*.bak` files from KiCad.
