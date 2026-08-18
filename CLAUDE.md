# fairchild — instructions for AI coding agents

Human contributors: see [`CONTRIBUTING.md`](CONTRIBUTING.md). This file is the
same information in the form an agent needs it.

## Build and test

```bash
cargo build --release
cargo test --workspace                 # 579 tests; ngspice goldens skip if ngspice is absent
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all
```

Enable the pre-commit hook once per clone — it runs the same `fmt` and `clippy`
gates CI fails on, over the whole working tree rather than the index:

```bash
git config core.hooksPath .githooks
```

Python bindings need `maturin develop --release` (maturin ≥ 1.8). The Rust
toolchain is pinned by `rust-toolchain.toml`; do not pin a version anywhere else.

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
