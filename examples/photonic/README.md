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
