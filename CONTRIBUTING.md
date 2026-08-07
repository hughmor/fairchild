# Contributing to fairchild

Bug reports, failing netlists, and models are all welcome. A netlist that
produces the wrong number is the single most valuable thing you can send —
attach the deck, what you expected, and what you got.

## Getting set up

```bash
git clone https://github.com/hughmor/fairchild
cd fairchild
cargo build --release
cargo test --workspace
git config core.hooksPath .githooks     # once per clone
```

The pre-commit hook runs the same `cargo fmt --check` and
`cargo clippy -D warnings` that CI fails on. It checks the whole working tree,
not just what you staged — so a partially staged commit still has to be clean
overall.

The Rust toolchain is pinned by `rust-toolchain.toml` and rustup installs it on
the first `cargo` call. Don't pin a version anywhere else; two sources of truth
is how they drift apart.

**ngspice** is needed for the comparison suites. Without it those tests *skip*
rather than fail, so you can develop without it — but CI installs it and asserts
it is on `PATH`, because nine suites silently comparing nothing is worse than a
missing dependency.

```bash
brew install ngspice suite-sparse        # macOS
sudo apt-get install ngspice libsuitesparse-dev
```

**Python bindings**: `maturin develop --release` (maturin ≥ 1.8 — 1.7 cannot
parse PEP 639 metadata and fails outright).

## The two rules that matter most

**A silent wrong answer is worse than a crash.** This project has shipped a
dozen of them and every one cost more to find than to fix: an inductor that was
an open circuit at DC, a line search that reported a 56 %-wrong operating point
as converged, a `.model` card that lost its first parameter and simulated the
defaults. They survive because the netlist runs and you get a plausible number.

So: prefer a hard error that names the fix. If a parameter cannot be honoured,
say so rather than ignoring it. If the solver cannot reach an answer, fail —
never stop somewhere and report it as converged. When you fix one of these,
write down in the commit message what it does *not* fix.

**A test is not finished until you have watched it fail.** Break the code the
test covers, run it, confirm it goes red, then put the code back. Several tests
in this tree passed against the exact bug they were written to catch. One
specific trap: a test that checks two subsystems agree with each other cannot
detect a fault they share — deleting a noise source from the shared list leaves
the frequency-domain and time-domain answers agreeing perfectly about a circuit
that is missing a generator. Agreement invariants need an absolute anchor on one
side.

## Conventions

**Commit messages explain why.** State the symptom, the cause, what you ruled
out, and what remains broken. `git log` is the design record here; commits are
long on purpose.

**Goldens may change** when the change increases correctness — say so explicitly
in the message and show the before and after.

**Update [`docs/model_status.md`](docs/model_status.md) in the same commit** as
any model change. It is an audit, so it goes stale silently, which is the worst
property a contract can have.

**One place interprets each concept.** `crate::reactive` owns "what is reactive"
and "what does `method` mean"; `crate::noise::NoiseSources` owns "what is a
noise source"; `crate::tolerance` owns "what unit is this row". Add to those,
don't reimplement at each consumer.

**Adding a photonic device touches three files**, and each will let you forget
the others:

| File | What |
|---|---|
| `crates/fairchild-core/src/device_registry.rs` | `register_native_photonics` — the card name |
| `scripts/kicad_to_fairchild.py` | `PORT_SCHEMA` — the port layout |
| `scripts/gen_kicad_symbols.py` | `SYMBOLS` — the schematic symbol |

Pin order must match the positional net order on the `X…` line. CI runs
`gen_kicad_symbols.py --check` and fails on a stale library.

The user guide's *Writing custom devices* section walks through the whole
process with a template.

## Pull requests

CI must be green: `fmt`, `clippy -D warnings`, and the full test suite on both
Linux and macOS, with and without the `klu` feature, plus a Python wheel build
and a C-API compile against the hand-written header.

Describe what you verified and how. "Tests pass" says less than "measured
−0.05 % against the closed form on three receivers, each dominated by a
different noise term".

## Licence

Contributions are accepted under the [Apache License 2.0](LICENSE).
