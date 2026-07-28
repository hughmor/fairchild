# Verilog-A in fairchild

Five worked examples and eight models, plus what the support amounts to.
The authoring guide is `docs/user-guide.md` §14; this is the worked set.

```sh
cargo build --release -p fairchild-cli          # from the repo root
OPENVAF=/path/to/openvaf-r ./build.sh           # compile models/*.va
./check.py                                      # run everything, assert the physics
```

---

## What the support is

fairchild does not parse Verilog-A. You compile it with
[OpenVAF-Reloaded](https://codeberg.org/arpadbuermen/OpenVAF-Reloaded) to a
`.osdi` shared library, and `crates/fairchild-osdi` `dlopen`s it and drives it
through the OSDI v0.4 ABI. This is the intended route to foundry models —
BSIM, PSP, HiCUM — which fairchild will never hand-write in Rust.

**Optical Verilog-A works too, and needs no compiler fork.** fairchild carries
an optical signal on ordinary real-valued MNA unknowns — three per channel
(`re`, `im`, `wl`) — so a custom `optical_field` / `optical_lambda` discipline
is just metadata that OSDI passes through untouched. `models/optical.vams` is a
self-contained copy of that discipline, with the rules that matter.

The two worlds interoperate **exactly**, not approximately: a native
`.optical_port p` expands to the wires `p_re_0 p_im_0 p_wl_0`, carrying field
amplitude in sqrt(W) and wavelength in metres, which is precisely what a
Verilog-A model on this discipline reads and writes. `wg_compare.sp` puts the
two waveguides side by side in one deck and `check.py` pins them to 1e-9.

What Verilog-A cannot reach is the rest of fairchild's optical abstraction
layer: WDM bundle-awareness, bidirectional propagation, `DelayLine` group
delay, and `PhotonicActiveModel` composition. A Verilog-A optical model is
single-channel and forward-only. Photonics needing those stay native Rust.

---

## Getting a model into a netlist

```spice
.osdi  build/va_diode.osdi                 ; path relative to the netlist
Xd1  a  out  va_diode  Is=1e-14 Rs=0.5     ; model name == Verilog-A module name
```

or, the foundry-PDK form — process parameters on a card, geometry per instance:

```spice
.osdi  build/va_nmos.osdi
.model nch  va_nmos (KP=120u VTH0=0.7 LAMBDA=0.02)
Mn  out in 0 0  nch  W=10u L=1u
```

`.osdi` registers every descriptor in the library under its module name;
`.model` then binds a card name to one, with the card's parameters as
defaults. Instance parameters win over the card. `X`, `M`, `Q` and `D` all
work.

---

## The examples

| netlist | what it shows |
|---|---|
| `rectifier.sp` | Verilog-A diode + native R, C, V — the simplest mixed deck |
| `cmos_inverter.sp` | Verilog-A transistors via `.model` cards — the PDK idiom |
| `eam_link.sp` | a Verilog-A modulator inside an otherwise native photonic link |
| `va_link.sp` | a Mach-Zehnder where *every* optical device is Verilog-A |
| `wg_compare.sp` | Verilog-A vs native waveguide, same deck, must agree |

### `rectifier.sp`

Half-wave rectifier. The Verilog-A diode (`models/va_diode.va`) has series
resistance, an internal node and a `ddt` junction charge, so it exercises the
three things the runtime has to get right: a nonlinear resistive branch,
`num_nodes > num_terminals`, and a reactive branch. Charges to Vpk − Vf ≈
4.17 V and droops ~0.35 V between peaks.

### `cmos_inverter.sp`

`.osdi` + `.model` card + `M` element, which is the shape every real PDK deck
has. Swap `va_nmos`/`va_pmos` for a `bsim4.osdi` and the deck does not change
shape. The edges overshoot a little — that is real Miller feedthrough through
the models' `CGDO`.

### `eam_link.sp`

An electro-absorption modulator (fairchild has no native EAM) between two
native `fc_waveguide` sections, fed by a native `fc_cw_laser` and read by a
native `fc_photodetector`. The coupling runs both ways: a native RC drive sets
the modulator bias, and the modulator's photocurrent — from the
electro-absorbed light only, so exactly zero unbiased — loads that same node.
9.0 dB extinction at the detector, matching `10^(−er_dB·(vr/v_full)²/10)` to
0.2 %.

### `va_link.sp`

A Mach-Zehnder with no native photonic device at all: `va_laser`, two
`va_coupler`s, two `va_waveguide`s of unequal length, two `va_photodetector`s.
Sweeping 1540–1575 nm walks one full FSR, and bar + cross holds at 0.9813 to
5e-4 — the whole-chain power budget, which is the real statement that the
coupler is unitary.

---

## Writing an optical model: two rules

**Take the wavelength from a parameter, never off the λ wire.** The wire
exists so a chain stays consistent and native devices can read it; do not use
`OWL(...)` inside an expression you contribute from. Propagation phase is
thousands of radians (400 µm at n_g = 4.2 is ~6800 rad), so OpenVAF
differentiating it against the λ unknown puts ∂φ/∂λ = φ/λ ≈ 1e9 per metre into
the Jacobian and Newton does not converge — at some wavelengths and not
others, which reads like a physics bug and is not one. Native devices freeze λ
at the previous iterate; Verilog-A has no way to say that. Sweep wavelength
with `--param` / `set_param`.

**dB is power dB, so amplitude divides by 20.** `10^(−dB/20)`. Everything
under `legacy/` predates the `0f689cb` fix and is a factor of two out.

---

## Limitations

All measured, not inferred:

- **OSDI reactive branches always integrate with Backward Euler**, whatever
  `.options method` says. Under `--method be` a Verilog-A `ddt(C*V)` matches a
  native `C` bit-for-bit; under `tr` it lags ~0.6 % on a 0.45 τ RC step. Pin
  `method=be` when a model carries charge and you care. (This is not
  OSDI-specific — every device-internal reactance in fairchild is BE.)
- **No limiting.** fairchild never calls `load_limit_rhs_resist` /
  `load_limit_rhs_react`, and OpenVAF 23.5 rejects `$limexp` outright — it is a
  compile error, not a silent no-op. In practice the Newton loop's Armijo line
  search carries a bare `exp()`; `va_diode.va` was checked to a 500 V drive.
  Clamp the exponent in the model if you do hit an overflow.
- **`prev_state` / `next_state` are null**, so state-carrying Verilog-A
  constructs are unavailable.
- **`--param` only reaches `X`, `R`, `C`, `L` elements**, so a Verilog-A
  transistor's parameters cannot be swept from the CLI. Use the Python
  bindings' `set_param`, or a `.param` in the netlist.

Fixed while writing these, so no longer limitations: `.model` cards resolving
to OSDI descriptors, `D`-element instance parameters, `$abstime` reading the
transient clock, and `ddt` reaching `.ac` / `.noise`.
