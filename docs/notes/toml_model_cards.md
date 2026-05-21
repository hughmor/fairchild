# Notes — TOML-driven model cards (parked)

Date: 2026-05-21.  Status: parked on branch `feature/toml-model-cards`
after working through a proof-of-concept; main branch reverted.

The goal was to let a physicist add a new device by writing one TOML
file + a KiCad symbol, no Rust boilerplate.  We got 7 of 19 devices
ported and the framework working end-to-end (build.rs codegen +
`Physics` trait split).  The user feedback was that the resulting TOML
ended up more complicated than the design sketch had suggested.

This note captures what we built, *why* the TOMLs grew, and the
specific design changes that would make them shorter when we revisit.

---

## What's on the parked branch

- `crates/fairchild-core/build.rs` (~410 lines) — TOML parser via
  serde+toml, code generator for struct + `new()` + full `Device` trait
  impl + registry factory.
- `crates/fairchild-core/src/models/physics.rs` — `Physics` trait
  with empty defaults for everything except `physics_eval` and
  `physics_load_jacobian` (those are required).
- `crates/fairchild-core/devices/*.toml` — 7 device cards:
  fc_grating_coupler, fc_splitter, fc_thermal_ps, fc_waveguide,
  fc_dcoupler, fc_pn_ps, fc_pn_th_ps.
- Generated code lives in `OUT_DIR/native_devices.rs`, included via
  `pub mod generated` in `models/mod.rs`.
- Build + all tests pass (except the pre-existing
  `wdm_transient_two_channels_diverge`).

Branch HEAD: `faaca1d refactor(devices): port fc_dcoupler, fc_pn_ps,
fc_pn_th_ps to TOML cards`.

---

## Where TOML complexity actually came from

The mental model was: "TOML lists ports + parameters; the rest is just
Rust".  Reality of the 7 ports we did:

| Pain point | Frequency | Cause |
|---|---|---|
| Multiple SPICE-card aliases per field (e.g. `length_m` / `L_m` / `length` / `l_um`) | every device | SPICE convention has wide-tolerance naming |
| Unit conversion on alias (`l_um` → `length_m` via ×1e-6) | every device | nm/µm/mm/dB-based aliases are universal |
| `default = "expr"` strings (e.g. `dB_per_cm_to_neper_per_m(2.0)`) | ≈half | derived constants that aren't trivially literal |
| `on_set = "method"` hooks | 1 device | derived fields that need recompute (tau_g_s on L/n_g change) |
| Setters that depend on other fields (`V_pi_L` → `dn_dv = wl_ref_m/(2·V_pi_L)`, `kappa_L` → `kappa_per_m = kappa_L/length_m`) | 2 devices | escape hatch — `physics_set_real_param` |
| Non-standard branch count (6/8 per channel, not wpc) | 2 devices | escape hatch — manual resize in `physics_setup_instance` |
| Caches section listing `Vec<f64>` fields | every device | should be auto from physics file usage |
| Per-field `docs` strings | every device | nice-to-have, accumulates lines |

The cards themselves grew to 50–90 lines.  The ratio of "this is the
unique physics of this device" to "this is parameter plumbing" was about
1 : 4.  Worse, the TOML schema **already had three escape hatches**
(`on_set`, `physics_set_real_param`, `physics_setup_instance` resize)
and we hit them on 3 of 7 simple devices.

---

## Specific design moves that would help if we revisit

### 1. Convention-over-config for unit aliases

Instead of declaring every variant:

```toml
[params.length_m]
default = 100e-6
aliases = ["L_m", "length"]
[params.l_um]
alias_of = "length_m"
scale = 1e-6
```

Make the canonical field `length_m` and have the codegen automatically
generate aliases for any name ending in `_m`:

- `length_um` / `L_um` / `l_um`  →  ×1e-6
- `length_nm`                    →  ×1e-9
- `length_mm`                    →  ×1e-3

Same for `_w`/_W, `_v`/_V, `_s`/_seconds.  Saves ~4 lines per param.

Similarly: for any `*_neper_m` field, auto-generate a `*_dB_cm` alias
that goes through `dB_per_cm_to_neper_per_m`.  Same for `*_amp` / `*_dB`.

### 2. Default = expr string is fine, just less of it

The "Rust expression as default" pattern is unavoidable for things like
`dB_per_cm_to_neper_per_m(2.0)`.  Just lean into it — every default is
a Rust expression; literal numbers are also valid expressions.  No need
for the `default = number` / `default = "string"` split serde currently
has.

### 3. Drop the `[caches]` section

The physics impl ALWAYS adds `c_cached`/`s_cached` for any rotating
device.  The codegen could instead provide a "physics-managed scratch
area" via a single `physics_scratch: Vec<f64>` field plus indexed
accessors — or just let the physics file add its own `impl Native<X> {
... }` block with extra inherent fields it manages.

Wait — that doesn't work, because Rust struct fields all live in one
place.  The closest reasonable thing: the TOML codegen always emits
the `Vec<f64>` caches as a standard pair `c_cached`, `s_cached` for
phase-rotation devices (declared by a `[caches] kind = "phase"`).
Devices that need something else (the cap device's `c_j_cached: f64`)
add a single line.

Actually I think the cleanest fix is: codegen generates the struct
fields from the TOML, but the physics file may declare *additional*
fields via a sibling `impl` and a `#[fairchild_aux_fields]` macro that
adds them post-codegen.  Too clever — drop this for v2.

For v2: just list cache field names but skip the type when it's
`Vec<f64>` (the default).  Two lines per cache becomes one.

### 4. Remove `on_set` — use a setter trait method instead

`on_set = "refresh_tau"` is used by ONE device.  Rather than a TOML
field, the codegen could always call `self.after_param_set(name)` (with
a default no-op).  Devices that need it override the method.  Saves a
TOML knob, moves one indirection to Rust where it's more visible.

### 5. Eliminate the `[branches]` section

Branch count modes (per_channel_wpc / per_channel / fixed / none) is a
quirk of the implementation, not a property of the device.  Better
default: codegen never allocates branches; the physics file's
`physics_setup_instance` is the canonical place to size them, with one-
or-two convenience helpers like `self.alloc_wpc_branches()`.  Most
devices write `self.alloc_wpc_branches();` in setup_instance and are
done.

### 6. Decide the SPICE-name policy

A real bikeshed but worth picking once:
- Lower-case in netlist, mixed-case in TOML field names?
- Case-insensitive matching always?
- Underscore vs no-underscore in aliases (`L_um` vs `Lum`)?

Right now the codegen lowercases everything before matching, which is
fine.  Document the rule and stop double-listing.

### 7. The big simplification — minimal TOML

After all the above, a fc_pn_ps card would look like:

```toml
[device]
name = "fc_pn_ps"
struct = "NativePnPhaseShifter"
docs = "Linear small-signal PN-junction phase shifter."

[ports]
bundles = ["in", "out"]
scalars = ["anode", "cathode"]

[params]
length_m       = 1e-3                                # auto-aliases: l_um, l_mm, l_nm
n_eff          = 2.7654
n_g            = 4.02                                # informational dispersion slope
wl_ref_m       = 1.55e-6                             # band-centre default applies
dn_dv          = 5.166666666666667e-5
g_pn           = 1e-3
alpha_neper_m  = "dB_per_cm_to_neper_per_m(20.0)"   # auto-aliases: alpha_dB_cm
pin_at_ref     = false
```

~20 lines vs. ~90.  The user wanted that and was right.

---

## Why we should still come back to this

The PoC isn't wrong, it's verbose.  The pain points above are all
mechanical and the framework's load-bearing parts (build.rs codegen,
Physics trait split, registry registrar) work end-to-end.  Specifically,
the auto-alias convention (item 1) accounts for ~60 % of every card's
length and is a pure codegen improvement — no schema change visible to
the device author.

For when we revisit:

1. Start by sketching what 3 representative cards (waveguide, pn_ps,
   photodetector) would look like *after* all the simplifications.  If
   they read clean, port-cost is real.
2. The hard cases that broke the schema cleanly are `kappa_L` /
   `V_pi_L` (length-aware setters).  Either accept the
   `physics_set_real_param` escape hatch or extend the `fn` form to
   accept `self`-aware expressions (`fn = "self.wl_ref_m / (2.0 *
   value)"`).
3. Photodetector + circulator + mux/demux were not even attempted.
   They're each ~1 special-case escape hatch.  If those each cost an
   `impl Physics` override that takes 10 lines, that's fine.  If they
   each break the schema, the schema needs a rework.

---

## Next steps if we ever resume

1. `git checkout feature/toml-model-cards` to get back to the working
   PoC state.
2. Apply schema simplifications 1–6 above (~1 day).
3. Re-port the 7 done devices on the new schema; if cards drop from
   ~70 lines to ~20 average, keep going.
4. Port `fc_cw_laser` (medium), `fc_thermal_ps_rc` (state variable),
   then the PN tier.  Stop after each batch for review.
5. The hard cases (photodetector, mux/demux, circulator) come last;
   if they each require a special-case TOML field, that's a signal the
   schema isn't right yet.

---

## What the user's reaction is telling us

User: "I’m not loving it so far.  The TOMLs are much more complicated
than we planned originally."

The original sketch was a 20-line card.  The actual ones are 50-90.
The framework works, but it broke the implicit promise of "you read
the TOML and *immediately* understand the device."  Until the cards
get back to "scannable in 30 seconds," the value isn't there.

Don't reach for proc-macros or Verilog-A.  The right answer is still
TOML-with-codegen; we just under-resourced the convention-over-config
pass before declaring victory.

---

*— filed by Claude with the parked PoC branch as supplementary evidence*
