# fairchild ↔ KiCad — Setup Guide

This guide is the **native-device flow**. It assumes you draw a schematic in
KiCad using symbols from a `fairchild_photonics.kicad_sym` library whose models
map onto fairchild's native `fc_*` photonic devices and onto built-in
electrical primitives (R / L / C / V / I / D / M). The pre-Phase-B path that
used compiled `.osdi` Verilog-A models is no longer recommended; that file
lives in git history if you need to refer back to it.

The end-to-end flow:

1. Draw the schematic in KiCad, using fairchild-native symbols.
2. Export to SPICE: `File → Export → Netlist → SPICE`.
3. Run `scripts/kicad_to_fairchild.py` to produce a wrapper netlist with
   `.optical_port` declarations and the analysis directive.
4. Run `fairchild -f run_my_circuit.sp`.

There is no OSDI preamble step, no `.osdi` directives, no 12-pin coupler
symbol. A directional coupler is a 4-pin symbol, a waveguide is a 2-pin
symbol, and so on — one KiCad pin per fairchild optical port.

---

## 1. Prerequisites

- **KiCad 10** (uses the `Sim.*` simulation-model dialog)
- **fairchild** built: `cargo build --release` produces `target/release/fairchild`
- **Python 3.9+** to run the post-processor

No OpenVAF, no `va-models/build/`, no `.osdi` artefacts are required for the
native flow.

---

## 2. Project layout

The fairchild repo ships the symbol library and reference KiCad examples
under `examples/kicad/`:

```
fairchild/
└── examples/
    └── kicad/
        ├── fairchild_photonics.kicad_sym   ← native fc_* symbol library
        └── mrm_kicad_test.cir              ← reference MRR-modulator export
```

Your own project lives outside the fairchild repo and looks like:

```
my-photonic-project/
├── my_circuit.kicad_sch
├── my_circuit.kicad_pro
├── my_circuit.cir                  ← KiCad SPICE export (regenerated each save)
└── run_my_circuit.sp               ← wrapper from kicad_to_fairchild.py
```

Add `fairchild_photonics.kicad_sym` via **Preferences → Manage Symbol Libraries
→ Project Libraries → Add**, pointing at
`<fairchild>/examples/kicad/fairchild_photonics.kicad_sym`. Update the
library path whenever you pull fairchild updates that touch the library.

---

## 3. Symbol library conventions

### One KiCad pin per fairchild port

Native fairchild photonic devices treat each optical port as a single bundle
that the parser expands into `(re, im, λ)` wires under the hood. The KiCad
symbol therefore has **one pin per optical port**, not three. Electrical pins
(`anode`, `cathode`, `heat_p`, `heat_n`) are also one pin each.

| Card | Optical ports (KiCad pins) | Electrical pins |
|---|---|---|
| `fc_cw_laser` | `out` | — |
| `fc_waveguide` | `in`, `out` | — |
| `fc_dcoupler` | `a1`, `a2`, `b1`, `b2` | — |
| `fc_splitter` | `in`, `out_a`, `out_b` | — |
| `fc_pn_ps` | `in`, `out` | `anode`, `cathode` |
| `fc_thermal_ps` | `in`, `out` | `heat_p`, `heat_n` |
| `fc_photodetector` | `in` | `anode`, `cathode` |

Pin **number** order matters — it is the order positional arguments appear
in the exported SPICE line. The table order above is also the device's
positional order; pin 1 = first listed, pin 2 = second, etc.

Use distinct KiCad pin types so optical and electrical pins are visually
distinguishable on the schematic:

- Optical pins: `Bidirectional` (signals propagate both ways through optical
  components in the SVEA model).
- Electrical drive pins: `Input` (anode, heat_p) and `Output` for
  photodetector outputs.

### SPICE model in the symbol editor

In the Symbol Editor, **press E** on the symbol → **Simulation Model**:

1. **Model tab**: choose **Raw SPICE Element**.
2. **Spice element type** = `X`.
3. **Model name** = the device card (`fc_dcoupler`, `fc_pn_ps`, etc.).
4. **Library path**: leave empty. Native devices are built into the
   `fairchild` binary; no external library file.
5. **Parameter list**: add each device parameter as `name` / `value` pair
   (e.g. `kappa_L = 0.336`).
6. **Pin Assignments tab**: bind each KiCad pin to its positional slot
   following the order in the table above.

KiCad will write `Sim.Device = SPICE`, `Sim.Type = X`, and the params to the
symbol's `Sim.*` fields. You do not need to edit those fields manually.

### Parameter reference

Each device's parameter set is documented in
`crates/fairchild-core/src/models/photonic.rs`. The most common knobs:

| Device | Common parameters |
|---|---|
| `fc_cw_laser` | `power_mW`, `wavelength_nm` |
| `fc_waveguide` | `L_um`, `n_g`, `alpha_dB_cm`, `wavelength_nm` |
| `fc_dcoupler` | `kappa_L` |
| `fc_pn_ps` | `L_um`, `V_pi_L`, `g_pn`, `dn_dv`, `alpha_dB_cm`, `n_g`, `wavelength_nm` |
| `fc_thermal_ps` | `L_um`, `R_heater`, `R_thermal`, `dn_dT`, `n_g`, `wavelength_nm` |
| `fc_photodetector` | `responsivity`, `i_dark_a`, `r_shunt` |

Parameter names match case-insensitively but are case-preserved through the
parser. Use the canonical case shown in the docs.

### Electrical primitives

Use KiCad's stock R / L / C / V / I / D / M / B symbols — their SPICE export
already produces fairchild-compatible cards. Set `Sim.Device` to `R`, `L`,
etc. and use the standard SPICE syntax for waveforms (`PULSE`, `SIN`, `PWL`,
`EXP`, `SFFM`, `AM`).

---

## 4. SPICE export configuration

**File → Export → Netlist → SPICE tab**:

1. **Output file** → `my_circuit.net` (alongside the schematic).
2. **Adjust passive symbol values**: ON, so KiCad emits canonical values.
3. **Save all voltages**: ON (cheap; the wrapper script doesn't filter).
4. **Save all currents**: ON.
5. **Prolog text**: **leave empty**. The wrapper script (§5) supplies the
   analysis and `.optical_port` declarations. KiCad would otherwise inject
   prolog text into the netlist itself, which is then overwritten on every
   export.

KiCad emits a file shaped like:

```
* my_circuit/my_circuit.kicad_sch
Xlaser laser_out fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xwg1 laser_out wg1_out fc_waveguide L_um=50 ...
Xdc wg1_out dc_b dc_c pn_in fc_dcoupler kappa_L=0.336
Xpn pn_in dc_b vmod 0 fc_pn_ps ...
Xwg2 dc_c pd_in fc_waveguide L_um=50 ...
Xpd pd_in pd_anode 0 fc_photodetector responsivity=0.8
Vbias bias 0 DC 1.0
Rload pd_anode bias 1k
Vmod vmod 0 PULSE(0 4 100n 100n 100n 800n 2u)
.end
```

That file alone is **not yet runnable** — the optical bundle nets need to be
declared, and an analysis directive is needed. That's the wrapper script's
job.

---

## 5. The wrapper script

`scripts/kicad_to_fairchild.py` reads the KiCad export, scans for native
`fc_*` X-elements, looks up each device's port schema, and identifies which
positional nets are optical bundles. It emits a wrapper file that:

- Declares every bundle net via `.optical_port`.
- Adds any `.options` you supply (`--opt`, `--method`).
- Appends the analysis directive (`.op` / `.tran` / `.ac`) **before** the
  `.include`, because the KiCad netlist conventionally ends with `.end`
  which terminates parsing.
- `.include`s the unmodified KiCad netlist.

### Common invocations

```bash
# DC operating point
python3 scripts/kicad_to_fairchild.py my_circuit.net -o run_my_circuit.sp
fairchild -f run_my_circuit.sp

# Transient with GEAR integrator
python3 scripts/kicad_to_fairchild.py my_circuit.net -o run_my_circuit.sp \
    --tran "5n 2u" --method gear

# AC sweep with tightened tolerance
python3 scripts/kicad_to_fairchild.py my_circuit.net -o run_my_circuit.sp \
    --ac "dec 20 1 1G" --opt reltol=1e-5

# Stream wrapper to stdout (no file)
python3 scripts/kicad_to_fairchild.py my_circuit.net --op
```

### What the script detects

For a KiCad export containing `Xdc wg1_out dc_b dc_c pn_in fc_dcoupler ...`,
the script knows `fc_dcoupler` has four bundle ports — so `wg1_out`, `dc_b`,
`dc_c`, and `pn_in` each get a `.optical_port` declaration. For
`Xpn pn_in dc_b vmod 0 fc_pn_ps ...`, only `pn_in` and `dc_b` are bundles;
`vmod` and `0` (ground) are scalar electrical nets and are left alone.

The script warns if:

- An `fc_*` instance has the wrong number of positional nets for its model.
- A bundle pin is wired to ground (likely a wiring mistake — `0` and `gnd`
  are reserved for electrical references, not optical signals).

Non-`fc_*` X-elements (e.g. `.subckt`-based hierarchical photonic blocks,
third-party OSDI models) are passed through untouched and reported with an
info message under `-v`.

### WDM via `fc_mux` / `fc_demux`

KiCad doesn't allow connecting a bus directly to a single symbol pin —
buses are explicitly for grouping multiple individual wires. The fairchild
WDM convention sidesteps this entirely: two new native devices, `fc_mux`
and `fc_demux`, bridge between N single-channel optical bundles and one
N-channel bundle. Every wire on the schematic stays a normal single-line
wire; the "bus-ness" lives in fairchild's `.optical_port NAME N` declaration
which the post-processor auto-generates.

Symbol convention:

- **`fc_mux`** (combiner): pin 1 is the bus output; pins 2..N+1 are N
  single-channel inputs. One symbol per supported N (typically 2, 4, 8).
- **`fc_demux`** (splitter): pin 1 is the bus input; pins 2..N+1 are N
  single-channel outputs.

The post-processor infers each bundle's channel count from MUX/DEMUX
instance arity (N = pin count − 1) and propagates the width along the
bus to every device between the MUX and DEMUX. Single-channel bundles
on the periphery (where lasers and detectors connect) stay at width 1.

Example transpiler output for a 2-channel WDM circuit:

```spice
.optical_port wdm_bus 2          ← bus, auto-detected from MUX1's 3-pin arity
.optical_port ch0                ← laser-side single channel
.optical_port ch1                ← laser-side single channel
.optical_port wg_out 2           ← propagated from MUX → WG1
.optical_port d0                 ← detector-side single channel
.optical_port d1
```

The MRR-modulator + MUX/DEMUX combo gives you the "single physical ring,
multiple wavelengths simultaneously" topology described in
`examples/photonic/native_wdm_mrr_modulator.sp`, but now drawable as a
KiCad schematic.

### Limits today (first PR scope)

- **Mixed widths on one device** — if two `fc_mux` instances feed into
  different parts of the same bundle name with different N, the
  post-processor uses the first authoritative width it sees and warns
  on conflicts. Don't reuse bundle names across mismatched widths.
- **No KiCad-bus expansion** — `net[0..3]` syntax is unused; the MUX/DEMUX
  approach makes it unnecessary. If a future user wants buses on the
  schematic for visual grouping, that's a separate post-processor feature.

---

## 6. First circuits to build (in order)

These build progressively, each one exercising one more piece of the flow.

### Circuit 1 — laser → photodetector (smoke test)

```
[fc_cw_laser] ──opt──▶ [fc_photodetector] ──▶ Rload(1k) ──▶ GND
                                           │
                                           └─ Vbias = 1.0 V
```

Validates: symbol library, SPICE export, bundle-net detection, wrapper
generation, fairchild end-to-end.

Expected V(pd_anode): 1 mW laser × 0.8 A/W responsivity → 0.8 mA into 1 kΩ
biased at 1 V → **V(pd_anode) ≈ 1.8 V** (off-resonance, full transmission).

### Circuit 2 — laser → waveguide → photodetector

Validates: passive propagation loss + parameter passing.

With `fc_waveguide L_um=1000 alpha_dB_cm=3.0`:
loss = 10^(−3.0 × 0.1 / 10) ≈ 0.93 → V(pd_anode) ≈ 1.75 V.

### Circuit 3 — micro-ring modulator

Build the topology in `examples/photonic/native_mrr_modulator.sp`:
laser → wg → directional coupler ⇄ PN phase shifter (ring) → wg → PD →
Rload + Vbias. Sweep `Vmod` from 0 to 4 V (use `PULSE` or a `.dc` sweep)
and verify the transmission swing from notch (~1.1 V) to full
transmission (~1.8 V).

### Circuit 4 — add-drop MRR

Same as Circuit 3 but with the ring's "drop" port routed to a second
photodetector. Verify that through-port and drop-port outputs sum to the
input intensity (modulo loss).

### Circuit 5 — Mach-Zehnder modulator

Two `fc_dcoupler` instances with `fc_pn_ps` in each arm. The bar/cross
outputs swap at `V = V_pi`.

---

## 7. Things to watch for

- **Pin number order**: the SPICE export uses positional pin numbers, not
  pin names. Renumbering a pin in the Symbol Editor silently changes the
  emitted netlist. Always verify by exporting and running through the
  wrapper script.
- **Reserved nets**: `0`, `gnd`, `GND` always mean ground. Don't use these
  as optical bundle names.
- **One bundle per optical net**: a single KiCad net carries one optical
  bundle. If you connect three fc_dcoupler outputs to the same KiCad net,
  you're connecting three optical bundles together — that may be what you
  want (an N-way merge) but it isn't always.
- **The `.end` placement**: KiCad emits `.end` at the bottom of the export.
  The wrapper script puts analysis + `.optical_port` directives **before**
  the `.include`, so this just works. If you ever hand-edit the wrapper,
  preserve that ordering.

---

## 8. Roadmap

The next steps for this integration, in rough order:

1. **`fairchild_photonics.kicad_sym` symbol library** — committed alongside a
   regression-test schematic.
2. **Bus-syntax auto-expansion in the wrapper** — translate KiCad buses
   `opt[0..3]` into `.optical_port opt 4` plus per-channel wiring.
3. **One-button driver** — ✅ `scripts/kicad_fairchild.py` chains
   schematic → SPICE export → transpile → simulate in one command, and can be
   registered as a KiCad schematic-editor "button" (see §9). Remaining: an
   in-editor results viewer.
4. **`fairchild.viewer`** Python helper — wrap matplotlib for common signal
   shapes (V vs t, optical power vs t for `(re, im)` pairs, AC magnitude /
   phase, noise PSD).

---

## 9. One-button driver — `scripts/kicad_fairchild.py`

`kicad_fairchild.py` collapses the whole §4–§5 flow into a single command:

```
schematic (.kicad_sch)  ──kicad-cli──▶  SPICE export (.cir)
                        ──kicad_to_fairchild.py──▶  fairchild netlist (run_*.sp)
                        ──fairchild -f──▶  results (CSV / rawfile)
```

It locates `kicad-cli` automatically (PATH, then the standard KiCad install
locations on macOS / Linux / Windows; override with `--kicad-cli`). `kicad-cli`
ships with KiCad 7 and later.

### Command-line use

```bash
# From a schematic (runs the whole pipeline and simulates):
python3 scripts/kicad_fairchild.py my.kicad_sch --tran "5n 2u" --run \
        --probe "V(pd_anode)" --sim-output results.csv

# From an already-exported .cir (skips kicad-cli — useful on a sim-only box):
python3 scripts/kicad_fairchild.py my_export.cir --tran "5n 2u" --run

# Enable the optical group-delay model for high-speed links:
python3 scripts/kicad_fairchild.py link.kicad_sch --tran "1p 2n" \
        --opt waveguide_delay=1 --run
```

### The schematic-editor "button"

KiCad's stable, version-portable hook for a schematic-side button is the
**BOM / netlist-generator** mechanism (eeschema → *Tools → Generate Bill of
Materials…* on KiCad 7/8; the *BOM* / netlist generator dialog on 9/10), **not**
the `pcbnew` action-plugin API (that only adds buttons to the PCB editor). Add a
generator with the command line:

```
python3 "/path/to/fairchild/scripts/kicad_fairchild.py" "%I" -o "%O" --tran "5n 2u"
```

eeschema passes the intermediate netlist as `%I`; the driver detects the XML,
recovers the source schematic from its `<design source=…>` field, exports SPICE
via `kicad-cli`, and transpiles. Press *Generate* to produce the fairchild
netlist; add `--run` to simulate in the same click.

> The command-line and `.cir` paths are exercised by the repo's example
> circuits. The KiCad-side button (the `kicad-cli` export step and the eeschema
> generator registration) depends on your local KiCad install and should be
> validated there once — it is not covered by CI.

### KiCad 9+ IPC API (future)

KiCad 9 introduced an out-of-process **IPC plugin API** that can add a true
toolbar button and read the live schematic without a netlist round-trip. That is
the eventual home for an integrated "Run fairchild + plot" action; the
BOM-generator hook above is the portable interim that works on every version
back to KiCad 7.
