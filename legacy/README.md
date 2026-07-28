# legacy/

Pre-Phase-B Verilog-A model library. **Historical reference — do not start
from these.** The maintained Verilog-A models, and the authoring guide, live in
`examples/verilog_a/`.

## What was promoted out of here

| here | promoted to | changed |
|---|---|---|
| `electronic/nmos_l1.va` | `va_nmos.va` | module renamed |
| `electronic/pmos_l1.va` | `va_pmos.va` | module renamed |
| `photonic/cw_laser.va` | `va_laser.va` | module renamed |
| `photonic/photodetector.va` | `va_photodetector.va` | module renamed |
| `photonic/waveguide.va` | `va_waveguide.va` | **loss convention fixed**, λ from parameter |
| `photonic/directional_coupler.va` | `va_coupler.va` | **ported off the Norton scheme**, λ unit fixed |

`electronic/diode_shockley.va` was not promoted — `examples/verilog_a/models/va_diode.va`
supersedes it (series resistance, an internal node, a `ddt` junction charge).
It stays because the OSDI loader tests reference it. `photonic/optical_source.va`
was not promoted either: it is a strictly weaker `cw_laser`.

## Two reasons not to reuse what is left

**The loss convention is a factor of two out.** Every optical model here
converts dB with `α·1e2/8.6859` — the *amplitude* neper constant — and then
halves again in the exponent. A 1 mm waveguide at 3 dB/cm passes 0.9661 of its
power here; the right answer is 0.9333. Native fairchild had the same bug until
commit `0f689cb`; these predate that fix and never got it.

**Twenty of the twenty-four photonic models are on the superseded `optical`
discipline** — the Norton scheme, where an output is driven as
`Oamp(out) <+ G·(target − Ophase(out))` with a 1e6 conductance sitting in the
matrix and 1e-12 shunts on the inputs. `disciplines/optical.vams` still defines
it for back-compatibility. The replacement is potential-only (`optical_field` /
`optical_lambda`): the model contributes its output directly and puts nothing
fictitious in the matrix. `examples/verilog_a/models/optical.vams` is the
maintained copy, and carries the wavelength-handling rules as well.

Those twenty are the MRR, MZI and phase-shifter families. They were **not**
ported, deliberately: each duplicates a native Rust device that is faster,
WDM-aware, and — for the ring modulators — fitted against real chip data (see
`experiments/giona/`). Port one only if you specifically need it in Verilog-A.
`va_waveguide.va` and `va_coupler.va` are enough to build a ring or an MZI from
scratch; `examples/verilog_a/va_link.sp` does exactly that.

## Building the .osdi binaries

Not tracked in git (`*.osdi` is gitignored). To rebuild:

```sh
cd legacy/va-models
bash build.sh          # requires openvaf-r on PATH
```

Or individually:

```sh
openvaf-r legacy/va-models/photonic/waveguide.va \
    --output legacy/va-models/build/waveguide.osdi
```

The OSDI tests under `crates/fairchild-osdi/tests/` skip with a helpful message
when these are absent, so a clean checkout still passes `cargo test`.
