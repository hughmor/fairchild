# fairchild — instructions for AI coding agents

Human contributors: see [`CONTRIBUTING.md`](CONTRIBUTING.md). This file is the
same information in the form an agent needs it.

## Build and test

```bash
cargo build --release
cargo test --workspace                 # ~2 min from clean; ngspice goldens skip
                                       # without ngspice, OSDI suites without
                                       # openvaf-r on PATH (FAIRCHILD_OPENVAF=<path>
                                       # also works). Run the whole thing — it is
                                       # cheaper than deciding what to skip.
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all
```

Integration tests are grouped one binary per subject
(`crates/*/tests/<subject>/main.rs`, with the tests themselves in files beside
it). Add a test to the file it belongs in, or add a file and one `mod` line —
**not** a new `tests/*.rs`, which Cargo compiles and links as its own binary.
That is what made the suite take 45 minutes to run 3 minutes of tests. Filters
read the same as ever: `cargo test --test native mzm`.

Enable the hooks once per clone — `pre-commit` runs the same `fmt` and `clippy`
gates CI fails on, over the whole working tree rather than the index, and
`pre-push` runs `cargo test --workspace`:

```bash
git config core.hooksPath .githooks
```

Python bindings need `maturin develop --release` (maturin ≥ 1.8). `cargo build`
does not rebuild them, so anything Python-facing — an example, a notebook, a fit
script — runs whatever was last built. `import fairchild` refuses a `.so` older
than the sources rather than answering from stale code; the fix is that same
command. The Rust toolchain is pinned by `rust-toolchain.toml`; do not pin a
version anywhere else.

## Where things live

| Path | What |
|---|---|
| `crates/fairchild-core/` | Solver: MNA, Newton, transient, AC, noise, native photonic devices |
| `crates/fairchild-parser/` | SPICE parser, expression grammar, bundle ports |
| `crates/fairchild-{cli,py,c,osdi,klu}/` | Binary, PyO3 bindings, C ABI, Verilog-A runtime, KLU backend |
| `docs/user-guide.md` | Everything a user needs to write and run a deck |
| `docs/photonic-models.md` | The photonic discipline and every `fc_*` device |
| `docs/model_status.md` | Per-parameter contract: parsed vs stamped vs validated |
| `docs/spice_support.md` | Which ngspice constructs work, and how the rest fail |

## Rules that matter here

**A silent wrong answer is the worst outcome.** This codebase has shipped
twelve of them and each cost more to find than to fix. Prefer a hard error
naming the fix. When you must choose between refusing to answer and answering
approximately, refuse — a solver may fail, it must not invent.

**Verify a test by sabotage.** After writing a test for non-trivial logic,
break the thing it covers and confirm it fails. Several tests in this tree
passed against a bug they were written to catch. An agreement invariant between
two subsystems cannot detect a fault common to both — it needs an absolute
anchor on one side.

**Tests are judged by what they would catch, not by how many there are.**
A test earns its place by failing on a plausible bug. Three shapes do not:

- *It ran.* `expect("solve")` with no value compared. `X_keyword_is_accepted`
  passes whether the parameter is implemented, dropped, or deleted.
- *It agrees with itself.* Two subsystems compared where a shared fault is
  invisible. One side must be an absolute anchor — a closed form, an analytic
  limit, another simulator.
- *It only covers the on state.* A feature test that always runs with the
  feature enabled cannot tell you the switch is wired up. `.options method=gear`
  parsed, stored, and then ran Backward Euler; every golden passed because none
  of them asked for anything but the default (#93).

The last one is the expensive one, because it is an absence rather than a
weakness — no existing test looks wrong. Prefer a **table with a completeness
gate**: enumerate the knobs, assert each does something, and make adding a knob
without an entry a failure. `tests/circuit/options_take_effect.rs` is the
pattern — it derives the field list from `SimOptions`'s own `Debug` output, so a
new option cannot be added silently, and it carries an explicit `NOT_OBSERVABLE`
list with a reason per entry rather than a silent gap.

Prefer one generalising test over five specific ones, and delete a test that
cannot fail rather than leaving it to be counted. A `#[test]` that only prints
is a debugging harness: mark it `#[ignore]`.

**One place interprets each concept.** `crate::reactive` owns "what is
reactive" and "what does `method` mean"; `crate::noise::NoiseSources` owns "what
is a noise source"; `crate::tolerance` owns "what unit is this row". Two lists
are two chances to disagree silently. When you add a source or a method, add it
there, not at each consumer.

**Adding a photonic device touches three places** and each will let you forget
the others: `register_native_photonics` in `device_registry.rs`, `PORT_SCHEMA`
in `scripts/kicad_to_fairchild.py`, and the `SYMBOLS` table in
`scripts/gen_kicad_symbols.py`. Pin order must match the positional net order on
the `X…` line.

**Update `docs/model_status.md` in the same commit** as any model change. It is
an audit, so it goes stale silently — the worst property a contract can have.

## Conventions

- Commit messages explain *why*, and name what was ruled out. Say what a fix
  does not fix.
- Never edit `LICENSE`; it is the canonical Apache-2.0 text and is checksum-verifiable.
- Goldens may change when the change increases correctness — say so in the message.
- Plots: always `MPLBACKEND=Agg`.
