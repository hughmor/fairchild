# Photonic examples

## Featured — native Rust devices (Phase B)

These examples use only the B-phase native photonic primitives
(`fc_cw_laser`, `fc_waveguide`, `fc_dcoupler`, `fc_splitter`, `fc_pn_ps`,
`fc_thermal_ps`, `fc_photodetector`).  No `.osdi` import, no Verilog-A,
no Norton hack.  These are the recommended starting points.

- **`native_mrr_modulator.sp`** — single-channel electro-optic micro-ring
  modulator: CW laser → waveguide → directional coupler → PN-loaded ring →
  waveguide → photodetector + load.  Transient analysis with a single
  0→Vπ→0 voltage pulse on the PN phase shifter shows the through-port
  transmission swing from a deep notch (V=0, on-resonance) to near-unity
  (V=Vπ, off-resonance).

- **`native_mrr_modulator.py`** — same circuit driven through the Python
  bindings; produces a three-panel plot (PN voltage, optical power at
  PD input, PD anode voltage).

- **`native_mrr_wavelength_sweep.py`** — sweeps laser wavelength across a
  full FSR at V_pn = 0 and V_pn = Vπ, plotting both transmission spectra
  to visualise the resonance and how the PN bias shifts it.

- **`native_mrr_bias_heater_sweep.py`** — the two static sweeps you run to
  characterise a fabricated add-drop MRM, on a `LEVEL=4` card whose defaults
  came from fitting a real silicon device (`experiments/giona/`). PN bias
  −1 → +1 V and heater current 0 → 1 mA, with through- and drop-port spectra
  plus Δλ_res versus each knob. Shows the three distinct tuning mechanisms
  side by side: depletion (+13 pm at −1 V, linear in V), carrier injection
  (−184 pm at +1 V, tracking the diode *current* — and visibly shallowing the
  notch through free-carrier absorption), and thermo-optic (+78 pm at 1 mA,
  linear in I²R). `--selftest` asserts all three, no plotting.

- **`native_wdm_mrr_modulator.{sp,py}`** — two-wavelength WDM extension
  of the modulator: two lasers (±50 pm around the ring resonance) share
  one bus through one ring driven by one V_pn.  The two photodetectors
  show very different transmission profiles — the same modulator
  produces a sharp notch on the red-side channel while the blue-side
  channel monotonically rises.  No multiplexer / demultiplexer device
  is needed: each wavelength is its own bundle channel, and the
  parser's per-channel replication does the routing automatically.

- **`native_weight_bank.py`** — WDM weight bank built on `fc_optical_2x2`, the
  behavioural 2×2 transfer block: one instance stands in for a cascade of ring
  modulators sharing a through bus and a drop bus, with an independent bipolar
  weight per wavelength feeding a balanced photodetector pair. Shows the two
  things a static ring netlist can't do — one control wire per channel via
  `.electrical_port wctl N` (so the netlist scales by changing one number), and
  weights swept by those control voltages *inside* a single `.tran`. Reach for
  this when the interesting physics is downstream of the weights (the O/E/O
  nonlinearity, the receiver) and the rings are just free parameters and a
  timestep constraint. `--selftest` asserts the weights, passivity, and that the
  balanced output tracks Σ w_k(t).

## PCell library (`pcells/`)

Parameterized `.subckt` cells, one file per component, meant to be `.include`d.
Every sub-device parameter is a cell parameter, and each instance gets its own
`.model` card built from them — so two instances of the same cell can differ.

- **`mrm.sp`** — add-drop micro-ring modulator. Four optical ports
  (in/thru/add/drop) and four electrical (PN anode/cathode, heater ±). `radius`
  (default 8 µm) drives the arc length as `{pi*radius}`; the full LEVEL=4 EO
  model (depletion, current-driven injection, TPA, thermo-optic) is exposed
  parameter by parameter, defaulting to the fitted giona values. Single-channel
  by design — read the header for what that means on a WDM bus.

- **`source_bank.sp`** — 8-channel WDM signal generator: eight lasers, each
  through its own ideal MZM, multiplexed onto one 8-channel optical bus. One
  optical output, eight electrical drives. Wavelengths and per-laser powers are
  parameters, so `p3=0` turns channel 3 off.

- **`native_pcell_link.py`** (in this directory) — wires both together: the
  source bank feeding a cascade of eight `mrm` instances, each trimmed to its
  own channel, each with independent PN bias and heater. Shows the drop matrix
  is diagonal, that a PN bias detunes exactly one ring, and that the MZM
  transfer matches `(1+cos(πV/Vπ))/2`. `--selftest` asserts all of it.

## Legacy (`legacy/`)

The older examples use the pre-B1 OSDI-based photonic models with the
Norton-hack discipline.  These still work — the OSDI compatibility shim
keeps loading them — but new users should prefer the native examples
above.  See `crates/fairchild-osdi/src/lib.rs` for the deprecation
rationale.

The legacy `mrr_compound.sp` / `mrr_adddrop_compound.sp` /
`mzi_pn_compound.sp` files demonstrate `.subckt`-based compound photonic
devices; the same compositions can be built today from native primitives
directly in the netlist.
