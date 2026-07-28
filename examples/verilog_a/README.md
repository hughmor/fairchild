# Verilog-A in fairchild

Two worked examples, plus what the support actually amounts to as of 2026-07-28.

```sh
cargo build --release -p fairchild-cli          # from the repo root
OPENVAF=/path/to/openvaf-r ./build.sh           # compile the .va models
./check.py                                      # run both examples, assert the physics
```

---

## Short answer

**Electrical Verilog-A: fully supported and load-bearing.** fairchild does not
parse Verilog-A itself. You compile it with
[OpenVAF-Reloaded](https://github.com/openvaf/OpenVAF-Reloaded) to a `.osdi`
shared library, and `crates/fairchild-osdi` `dlopen`s it and drives it through
the OSDI v0.4 ABI. This is the intended route to foundry models — BSIM, PSP,
HiCUM — which fairchild will never hand-write in Rust.

**Optical Verilog-A: also supported, and no compiler fork is needed.** This was
the open question. It works because fairchild carries an optical signal on
ordinary real-valued MNA unknowns — three per channel (`re`, `im`, `wl`) — so a
custom `optical_field` / `optical_lambda` discipline is just metadata that OSDI
passes through untouched. `models/optical.vams` is a self-contained copy of
that discipline, and `models/va_eam.va` is a working modulator built on it.

The two worlds interoperate **exactly**, not approximately: a native
`.optical_port p` expands to the wires `p_re_0 p_im_0 p_wl_0`, carrying field
amplitude in sqrt(W) and wavelength in metres, which is precisely what a
Verilog-A model on this discipline reads and writes. `eam_link.sp` puts a
Verilog-A modulator between a native `fc_waveguide` on each side and it all
solves in one Newton loop.

What is *not* reachable from Verilog-A is the rest of fairchild's optical
abstraction layer: WDM bundle-awareness, bidirectional propagation,
`crate::delay::DelayLine` group delay, and the `PhotonicActiveModel` drive
composition. A Verilog-A optical model is single-channel and forward-only. The
division of labour in `crates/fairchild-osdi/src/lib.rs` still holds —
photonics that need those features stay native Rust — but "can we write optical
Verilog-A at all" is a clear yes.

---

## How a netlist reaches a Verilog-A model

```spice
.osdi build/va_diode.osdi                       ; path relative to the netlist
Xd1  a  out  va_diode  Is=1e-14 Rs=0.5          ; model name == Verilog-A module name
```

`.osdi` loads the library and registers every descriptor in it into the same
`DeviceRegistry` the native models live in; from there the model name resolves
like any other. Both the CLI and the Python extension do this
(`fairchild-cli/src/main.rs:368`, `fairchild-py/src/lib.rs:757`); both build
with the `osdi` feature on by default.

**Use the `X` prefix.** An OSDI model with two terminals will also instantiate
as `D1 a b va_diode`, but the `D`/`M`/`Q` element parsers stop reading at the
model name (`spice/element.rs:75`), so instance parameters on those lines are
silently discarded — and a `.model` card does not reach an OSDI device either
(verified: `rs=500` on a card changed nothing). `X` is the only form that
parameterises one.

---

## The examples

### `rectifier.sp` — electrical

Half-wave rectifier: a Verilog-A diode (`models/va_diode.va`) with series
resistance, an internal node, and a `ddt` junction charge, feeding a native
`Rload`/`Cload` from a native `Vin`/`Rsrc`. Charges to Vpk − Vf ≈ 4.17 V and
droops ~0.35 V between peaks.

The model exercises the three things the OSDI runtime has to get right: a
nonlinear resistive branch, an internal node (`num_nodes > num_terminals`), and
a reactive branch.

### `eam_link.sp` — optical + electrical

An electro-absorption modulator — which fairchild has no native equivalent for
— dropped into a link that is otherwise entirely native:

```
fc_cw_laser ─► fc_waveguide ─► [ va_eam ] ─► fc_waveguide ─► fc_photodetector ─► Rload
                                   ▲
                         Rdrv + Cpar (native RC, 200 ps)
                                   ▲
                               Vdrv pulse
```

The coupling runs both ways. The native RC drive sets the modulator bias, and
the modulator's photocurrent — generated only by the electro-absorbed light, so
it is exactly zero unbiased — loads that same drive node. Measured 9.0 dB
extinction at the detector, matching `10^(-er_dB·(vr/v_full)²/10)` to 0.2 %.

---

## Limitations found while building these

All measured on this branch, not inferred:

- **OSDI reactive branches always integrate with Backward Euler**, whatever
  `.options method` says — `load_jacobian_tran` gets `alpha = 1/step`
  unconditionally. Under `--method be` a Verilog-A `ddt(C*V)` matches a native
  `C` bit-for-bit; under `tr` it lags by ~0.6 % on a 0.45 τ RC step. Pin
  `method=be` when a Verilog-A model carries charge and you care.
- **`ddt` is invisible to `.ac`.** `OsdiDevice` does not implement
  `Device::small_signal_reactances`, so only the resistive Jacobian is stamped.
  Changing `Cj0` from 2 pF to 2 µF leaves the AC response bit-identical.
  DC, `.dc` and `.tran` are all fine.
- **No limiting.** fairchild never calls `load_limit_rhs_resist` /
  `load_limit_rhs_react`, so there is no `pnjlim` equivalent for a Verilog-A
  junction. `$limexp` is not an escape hatch — OpenVAF 23.5 rejects it outright
  ("`$limexp` was not found in the current scope"). In practice the Newton
  loop's Armijo line search carries a bare `exp()` fine; `va_diode.va` was
  checked to a 500 V drive. Clamp the exponent in the model if you do hit an
  overflow.
- **`$abstime` always reads 0** — `OsdiSimInfo.abstime` is hardcoded
  (`fairchild-osdi/src/device.rs`). `prev_state` / `next_state` are null, so
  state-carrying Verilog-A constructs are not available either.
- **Instance parameters only via `X`** — see above.

## A bug this work fixed

`OsdiDevice::load_residual_tran` passed the *previous timestep's* solution as
OSDI's `prev_solve`. OpenVAF's `load_spice_rhs_tran` linearises about whatever
vector it is handed (`J·prev_solve − f`), so the residual and the Jacobian were
taken about different points, and any nonlinear Verilog-A model was
Newton-inconsistent the moment its operating point moved. Symptom: `.op` and
`.dc` were correct, a constant-source `.tran` was correct, and the first moving
source failed with "Newton-Raphson did not converge after 150 iterations". A
model with an internal node failed even on a constant source.

The fix takes the resistive part about the current iterate and recovers the
reactive history as the difference of the two OSDI entry points evaluated at
the previous timestep. Regression test:
`crates/fairchild-osdi/tests/osdi_tran.rs`.

## `legacy/va-models/` — read before reusing

That tree has 28 Verilog-A models, including a full photonic set
(waveguide, coupler, MRRs, MZIs, phase shifters, PD, laser), and they still
load and run. They are also unmaintained, and **their loss convention is a
factor of two out**: they divide dB by 8.6859 and then halve again in the
amplitude exponent. A 1 mm waveguide at 3 dB/cm passes 0.966 of its power
through the legacy model and 0.933 through native `fc_waveguide`; 0.933 is
correct (commit `0f689cb`). The models here use `10^(-dB/20)`, which agrees
with native.
