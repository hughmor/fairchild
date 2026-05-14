# fairchild KiCAD Integration — Setup Guide

## Overview

**Split of work:**

- **You (KiCAD side):** symbol library, schematics, SPICE export config
- **Me (fairchild side):** WDM bus vector expansion in the parser, post-processor script for KiCAD → fairchild netlist conversion

These can proceed in parallel. The validation step at the end brings them together.

---

## Prerequisites

- **KiCAD 10** (these instructions reflect KiCAD 10's `Sim.*` simulation model dialog)
- All fairchild VA models compiled to `.osdi` in `va-models/build/`
- `fairchild` CLI binary at `target/release/fairchild`

---

## 1. Project structure

Create a KiCAD project alongside your fairchild netlists:

```
photonic-circuits/
├── fairchild_photonics.kicad_sym   ← your symbol library
├── fairchild_preamble.sp           ← .osdi + .optical declarations (see §4)
├── my_circuit.kicad_sch            ← schematics
└── my_circuit.kicad_pro
```

---

## 2. Symbol library — conventions

### Create the library

**Preferences → Manage Symbol Libraries → Project Libraries tab → Add ("+").**  
Name it `fairchild_photonics`, point it to `fairchild_photonics.kicad_sym`.  
Then **open the Symbol Editor**, select that library, and create symbols.

### Optical pin convention

Each optical port is represented as **three adjacent KiCAD pins**:

| Suffix | Role | KiCAD Pin Type |
|--------|------|----------------|
| `_re` | Real amplitude (√W) | Bidirectional |
| `_im` | Imaginary amplitude (√W) | Bidirectional |
| `_wl` | Wavelength (µm) | Bidirectional |

Group these visually on the symbol body — e.g., draw a small `[ ]` bracket on the optical side and label it with the port name. Electrical pins use the standard Input/Output/BiDir types as appropriate.

**Pin naming is critical** — the pin name becomes the net name in the KiCAD schematic, which becomes the node name in the exported SPICE netlist. Use exactly these suffixes (`_re`, `_im`, `_wl`) because the `.optical` directive in the preamble (§4) uses them by pattern.

### SPICE simulation model (KiCAD 10)

On each symbol, press **E** (or right-click → Edit Simulation Model). In the dialog:

1. **Model tab**: select **Raw SPICE Element**
2. Set **Spice element type** = `X`
3. Set **Model name** = the OSDI model name (e.g. `cw_laser`)
4. Leave **Library path** empty — `.osdi` loading is handled by the fairchild wrapper, not KiCAD
5. In the **Parameter list** below: add each default parameter as `name` / `value` pairs
6. Switch to the **Pin Assignments tab**: assign node numbers to each KiCAD pin following the VA terminal order tables below

KiCAD auto-creates `Sim.Device = SPICE` and `Sim.Type = X` in the symbol's Fields. `Sim.Params` and `Sim.Pins` are populated from the dialog. You don't need to edit these Fields directly.

**Note on parameter case**: fairchild matches parameters case-insensitively (`Vpi_L` = `vpi_l`). Enter them with the canonical case shown below for readability.

---

## 3. Symbol catalog

For each symbol, number pins **strictly in the order the VA module declares them**. The pin number determines its position in the exported SPICE X-element line — position 1 is VA terminal 1, etc.

### `cw_laser`

VA terminals: `out_re, out_im, out_lambda`

| Pin # | Name | Type | Side |
|-------|------|------|------|
| 1 | `out_re` | BiDir | Right |
| 2 | `out_im` | BiDir | Right |
| 3 | `out_wl` | BiDir | Right |

Default parameters: `power_mW=1.0 wavelength_nm=1550.0`

Shape suggestion: rectangle with a laser diode chevron on the left, 3 optical pins on the right.

---

### `waveguide`

VA terminals: `in_re, in_im, in_lambda, out_re, out_im, out_lambda`

| Pin # | Name | Type | Side |
|-------|------|------|------|
| 1 | `in_re` | BiDir | Left |
| 2 | `in_im` | BiDir | Left |
| 3 | `in_wl` | BiDir | Left |
| 4 | `out_re` | BiDir | Right |
| 5 | `out_im` | BiDir | Right |
| 6 | `out_wl` | BiDir | Right |

Default parameters: `L_um=100 n_g=4.2 alpha_dB_cm=2.0 wavelength_nm=1550.0`

Shape suggestion: a horizontal rectangle (bus/wire symbol), inputs left, outputs right.

---

### `directional_coupler`

VA terminals: `a1_re, a1_im, a1_lambda, a2_re, a2_im, a2_lambda, b1_re, b1_im, b1_lambda, b2_re, b2_im, b2_lambda`

Port labeling: `a1`=input bus in, `a2`=input bus feedback/second port, `b1`=through, `b2`=cross/coupled

| Pin # | Name | Side |
|-------|------|------|
| 1–3 | `a1_re`, `a1_im`, `a1_wl` | Left-top |
| 4–6 | `a2_re`, `a2_im`, `a2_wl` | Left-bottom |
| 7–9 | `b1_re`, `b1_im`, `b1_wl` | Right-top |
| 10–12 | `b2_re`, `b2_im`, `b2_wl` | Right-bottom |

Default parameters: `kappa_0=0.5 wavelength_nm=1550.0`

Shape suggestion: two diagonal crossing lines (classic DC symbol) inside a rectangle. Top-left = a1 port group, bottom-left = a2, top-right = b1, bottom-right = b2.

---

### `pn_phase_shifter_l1`

VA terminals: `in_re, in_im, in_lambda, out_re, out_im, out_lambda, anode, cathode`

| Pin # | Name | Type | Side |
|-------|------|------|------|
| 1–3 | `in_re`, `in_im`, `in_wl` | BiDir | Left |
| 4–6 | `out_re`, `out_im`, `out_wl` | BiDir | Right |
| 7 | `anode` | Input | Bottom |
| 8 | `cathode` | Input | Bottom |

Default parameters: `L_um=500 n_g=4.2 alpha_dB_cm=3.0 Vpi_L=2.0 V_ref=0.0 wavelength_nm=1550.0`

---

### `thermo_phase_shifter_l1`

VA terminals: `in_re, in_im, in_lambda, out_re, out_im, out_lambda, heat_p, heat_n`

Same shape as PN PS, but `heat_p`/`heat_n` instead of `anode`/`cathode`. Use a resistor symbol on the bottom to indicate the heater.

Default parameters: `L_um=500 n_g=4.2 alpha_dB_cm=2.5 R_heater=1000 R_thermal=50000 dn_dT=1.86e-4 wavelength_nm=1550.0`

---

### `photodetector`

VA terminals: `in_re, in_im, in_lambda, anode, cathode`

| Pin # | Name | Type | Side |
|-------|------|------|------|
| 1–3 | `in_re`, `in_im`, `in_wl` | BiDir | Left |
| 4 | `anode` | Output | Right |
| 5 | `cathode` | Output | Right |

Default parameters: `responsivity=1.0`

Shape suggestion: photodiode triangle with an arrow for optical input on the left.

---

### Monolithic compound models

These are the models to use for MRR and MZI (the compound `.spc` files have a known convergence issue for MZI).

#### `mrr_modulator_l1` (all-pass)

VA terminals: `in_re, in_im, in_lambda, out_re, out_im, out_lambda, anode, cathode`

Same pin layout as `pn_phase_shifter_l1`. Label ports `IN` and `THRU`. Draw a ring resonator glyph inside the symbol body.

Default parameters: `kappa_0=0.1 L_ring_um=100 n_g=4.2 alpha_dB_cm=2.0 Vpi_L=2.0 V_ref=0.0 wavelength_nm=1550.0`

#### `mrr_modulator_l1_adddrop` (add-drop)

VA terminals: `in_re, in_im, in_lambda, th_re, th_im, th_lambda, dp_re, dp_im, dp_lambda, ad_re, ad_im, ad_lambda, anode, cathode`

| Pin group | Nets | Side |
|-----------|------|------|
| IN port | `in_re/im/wl` | Left-top |
| THRU port | `th_re/im/wl` | Right-top |
| DROP port | `dp_re/im/wl` | Right-bottom |
| ADD port | `ad_re/im/wl` | Left-bottom |
| `anode`, `cathode` | — | Bottom |

Default parameters: `kappa_0=0.1 L_ring_um=100 L_coup_um=10 n_g=4.2 alpha_dB_cm=2.0 Vpi_L=2.0 V_ref=0.0 wavelength_nm=1550.0`

#### `mzi_modulator_pn_l1`

VA terminals: `in_re, in_im, in_lambda, bar_re, bar_im, bar_lambda, cross_re, cross_im, cross_lambda, anode, cathode`

| Pin group | Nets | Side |
|-----------|------|------|
| IN port | `in_re/im/wl` | Left |
| BAR port | `bar_re/im/wl` | Right-top |
| CROSS port | `cross_re/im/wl` | Right-bottom |
| `anode`, `cathode` | — | Bottom |

Default parameters: `kappa_0=0.5 L_arm_um=500 n_g=4.2 alpha_dB_cm=3.0 Vpi_L=2.0 V_ref=0.0 wavelength_nm=1550.0`

#### `mzi_modulator_thermo_l1`

Same layout as `mzi_modulator_pn_l1` but terminals end with `heat_p, heat_n` instead of `anode, cathode`.

---

## 4. Preamble file

KiCAD cannot generate `.osdi` or `.optical` directives. Create `fairchild_preamble.sp` in your project directory:

```spice
* fairchild_preamble.sp — auto-included by fairchild wrapper netlist
* Update paths to match your fairchild repo location.

.osdi /path/to/fairchild/va-models/build/cw_laser.osdi
.osdi /path/to/fairchild/va-models/build/waveguide.osdi
.osdi /path/to/fairchild/va-models/build/directional_coupler.osdi
.osdi /path/to/fairchild/va-models/build/pn_phase_shifter_l1.osdi
.osdi /path/to/fairchild/va-models/build/thermo_phase_shifter_l1.osdi
.osdi /path/to/fairchild/va-models/build/photodetector.osdi
.osdi /path/to/fairchild/va-models/build/mrr_modulator_l1.osdi
.osdi /path/to/fairchild/va-models/build/mrr_modulator_l1_adddrop.osdi
.osdi /path/to/fairchild/va-models/build/mzi_modulator_pn_l1.osdi
.osdi /path/to/fairchild/va-models/build/mzi_modulator_thermo_l1.osdi

* Declare optical discipline for all nets that appear in your circuit.
* Add a .optical line for every optical net group (re/im/wl triples).
* Example for a single laser → PD circuit:
* .optical lre lim wl  ore oim

* NOTE: you will need one .optical line listing ALL optical nodes in your schematic.
* A post-processor script will generate this automatically from the netlist.
```

Create a wrapper netlist `run_circuit.sp` that includes both the preamble and the KiCAD export:

```spice
* run_circuit.sp — wrapper that adds fairchild directives around the KiCAD export
.include "fairchild_preamble.sp"
.include "my_circuit.net"     * ← KiCAD-exported SPICE netlist filename

.op
.end
```

Then run with:
```bash
DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib \
  /path/to/fairchild/target/release/fairchild -f run_circuit.sp --probe "V(ph_a)"
```

---

## 5. KiCAD SPICE export configuration

1. **File → Export → Netlist → SPICE tab**
2. Set **Output file** to `my_circuit.net`
3. In the **"Custom net format"** section:
   - Leave Format = `SPICE`
   - **Uncheck** "Use default reference designators" so `Xlaser` comes out as `Xlaser` not `XU1`
4. **Do NOT** add the prolog text in KiCAD itself — use the wrapper `.sp` approach from §4 instead. This keeps the preamble separate from the auto-generated file (which gets overwritten every export).
5. Click **Export Netlist**.

The exported `.net` file will look like:
```spice
* Created by KiCAD ...
Xlaser lre lim wl cw_laser power_mW=1.0 wavelength_nm=1550.0
Xpd    ore oim wl ph_a 0 photodetector responsivity=1.0
Rload  ph_a 0 1k
```

---

## 6. WDM bus naming convention (Tier 2)

Bus vector expansion (`net[0..3]` → `net_0, net_1, net_2, net_3`) is being added to the fairchild parser. **Use this naming convention in your schematics now** so everything works automatically once it lands.

For a 4-channel WDM bus:
- Use **net labels** in KiCAD: `opt_re_0`, `opt_re_1`, `opt_re_2`, `opt_re_3` (and `_im_*`, `_wl_*`)
- Or draw a **bus** in KiCAD with member nets `opt_re[0..3]`, which KiCAD expands to those individual net names in the SPICE export

For the `.optical` directive in the preamble, the pattern will be:
```spice
.optical_bus 4  opt_re opt_im opt_wl
* expands to: .optical opt_re_0 opt_im_0 opt_wl_0 opt_re_1 opt_im_1 opt_wl_1 ...
```

---

## 7. First circuits to build (in order)

Build and test these in sequence to validate the full flow:

### Circuit 1 — laser → PD (smoke test)

Validates: symbol library, SPICE export, preamble, `.optical` declaration.

```
[cw_laser] --lre,lim,wl--> [photodetector] --> Rload(1k) --> GND
```

Expected: `V(ph_a) ≈ 1.0` (1 mW × 1 A/W × 1 kΩ = 1 V)

### Circuit 2 — laser → waveguide → PD

Validates: passive propagation loss. With `L_um=1000, alpha_dB_cm=3.0`:
loss = 10^(−3.0 × 0.01 / 10) ≈ −0.069 dB → V(ph_a) ≈ 0.985 V

### Circuit 3 — all-pass MRR modulator

Use `mrr_modulator_l1` symbol. Sweep wavelength to verify resonance dip at the expected wavelength:
```
λ_res = L_ring × n_g / m
```
For L=100 µm, n_g=4.2: resonances every ~5.7 nm, nearest to 1550 nm ≈ 1544 nm.

### Circuit 4 — add-drop MRR filter

Use `mrr_modulator_l1_adddrop`. At resonance, power should split between THRU and DROP with the ratio determined by kappa and ring loss. Use Rload on both output PDs and check thru + drop ≈ input.

### Circuit 5 — MZI modulator

Use `mzi_modulator_pn_l1`. With `kappa_0=0.5` (balanced), at `V=0`: bar ≈ 1.0, cross ≈ 0. At `V = Vpi_L / (2 × L_arm_cm)` = 2.0 V·cm / (2 × 500e-4 cm) = 20 V: bar ≈ 0, cross ≈ 1.0.

---

## 8. Things to watch for

**Net name collisions**: KiCAD uses the reference designator prefix in subcircuit expansion. Make sure optical net names in your schematic are unique and don't accidentally match electrical net names.

**`wl` net sharing**: All devices on the same optical bus MUST share the same `wl` (wavelength) net. In a single-wavelength circuit, use one global `wl` net label and connect it to every optical device's `*_wl` pins.

**`.optical` declaration must list every optical net**: If you forget a net, discipline checking will either error or silently treat it as electrical. The post-processor script will generate this line automatically from the exported netlist.

**ADD port floating**: For the add-drop MRR, the ADD port is normally unused. Float it by connecting `add_re`, `add_im`, `add_wl` to unique undriven nets (not GND — they're optical nodes). In fairchild, undriven optical nets default to zero amplitude, which is correct for an empty input port.
