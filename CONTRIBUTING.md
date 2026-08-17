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

## Releasing

Three channels ship from one tag, because they serve different people:
**PyPI** (`fairchild-sim` — wheels, both the module and the `fairchild`
command, no toolchain needed), **crates.io** (source, for Rust consumers),
**GitHub Releases** (prebuilt CLI binaries and the C library, which has nowhere
else to live).

The version lives in exactly one editable place, `[workspace.package]` in the
root `Cargo.toml`. Everything else derives from it: crates via
`version.workspace = true`, the wheel via `dynamic = ["version"]`. The internal
dependency versions in `[workspace.dependencies]` are unavoidable copies —
Cargo has no inheritance there and crates.io will not publish without them —
so `scripts/check_versions.sh` exists to make a stale copy fail CI rather than
ship. Run it any time; it takes no arguments.

To cut a release:

1. Bump `[workspace.package] version`; open it as a PR.
2. **Merge it, then tag the merged commit.** Not before. A squash-merge
   rewrites the commit, so a tag pushed from the branch is left pointing at
   something that never reaches `master` — this has already happened once.
3. `git push origin vX.Y.Z`. `release.yml` does the rest, and refuses to
   publish anything if the tag and the workspace version disagree.

Exercise the pipeline without publishing via **Actions → Release → Run
workflow** with `dry_run` checked. Everything builds and the wheel is
install-tested; nothing is uploaded. Worth doing before any release you care
about — a wheel job that only breaks on macOS is otherwise discovered on
announcement day.

Versions are permanent. crates.io cannot delete a published version, only yank
it, and PyPI will not accept a re-upload of the same number. There is no
untagging a mistake, only a follow-up release.

## Licence

Contributions are accepted under the [Apache License 2.0](LICENSE).
