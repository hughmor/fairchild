#!/usr/bin/env python3
"""What does the OSDI path cost against a native Rust model?

The OSDI/Verilog-A path is how every foundry model reaches fairchild — BSIM,
PSP, HiSIM, HICUM, MEXTRAM — so its overhead sets the cost of the whole
"bring your own PDK" story. This measures it against native models doing
identical arithmetic.

Three model pairs, chosen so the result can be attributed:

  resistive      `osdi_g_shunt.va` against a discrete R. Almost no arithmetic and
                 no `ddt`, so this is close to pure per-eval ABI overhead: one
                 indirect call, the residual and Jacobian copies, the sim-info
                 struct.
  reactive       `rc_shunt.va` against a discrete R and C. The same model plus one
                 `ddt()` term, so the gap between these two ratios is the cost of
                 OSDI's reactive path with the arithmetic held constant.
  nonlinear_exp  `osdi_shockley.va` against the native `D`. One `exp` per eval,
                 which is what a real model's inner loop looks like, and no
                 reactive term.

Correctness first: each pair is checked to agree before any timing is reported.
A timing comparison between two models that do not compute the same thing is
not a measurement of overhead. All three pairs currently agree *exactly* — the
Verilog-A Shockley diode and the native one are bit-identical.

Read the bottom panel, not the ratios. A ratio is relative to the native side, so
`nonlinear_exp` reads better than `resistive` (6.1x vs 7.4x) only because a native
diode is slower than a native resistor.

Follows benchmarks/METHODOLOGY.md — default options only, same netlist shape on
both sides, and failures are reported rather than dropped. Timings are
**interleaved** A/B rather than batched: this machine drifts up to ~45% run to
run, which is larger than the effect being measured.

Usage:
    cargo build --release
    python benchmarks/osdi_vs_native.py [--out FIG.png] [--reps N] [--json F]
"""

import argparse
import json
import os
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
FC = REPO / "target" / "release" / "fairchild"
MODELS = Path(__file__).resolve().parent / "models"
OSDI_MODELS = REPO / "crates" / "fairchild-osdi" / "tests" / "models"

# Device counts. A ladder of N devices, so the per-device cost is the slope and
# process startup is the intercept — startup dominates below ~20 nodes and the
# fit is what removes it.
COUNTS = [8, 32, 128, 512, 2048]


def die(msg):
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(1)


def compile_va(src: Path, cache: Path) -> Path:
    """Compile a `.va` to `.osdi`, returning the artefact path."""
    compiler = os.environ.get("FAIRCHILD_OPENVAF") or "openvaf-r"
    out = cache / (src.stem + ".osdi")
    r = subprocess.run(
        [compiler, str(src), "-o", str(out)], capture_output=True, text=True
    )
    if r.returncode != 0:
        die(f"compiling {src.name} failed:\n{r.stdout}\n{r.stderr}")
    return out


# ---------------------------------------------------------------------------
# Deck construction
# ---------------------------------------------------------------------------
#
# Every deck is a ladder: a source, then N copies of the device under test each
# with a series resistor to the next node. Same topology on both sides, so the
# matrix is the same size and sparsity and only the device differs.


def ladder(n, device_line, extra_cards="", analysis=".op\n"):
    lines = ["* osdi vs native ladder", extra_cards.rstrip("\n")]
    lines = [l for l in lines if l]
    lines.append("V1 n0 0 DC 0.7")
    for i in range(n):
        lines.append(f"R{i} n{i} n{i + 1} 1k")
        lines.append(device_line.format(i=i, node=f"n{i + 1}"))
    lines.append(analysis.rstrip("\n"))
    return "\n".join(lines) + "\n"


DECKS = {
    # The control. `resistive` and `reactive` differ only in a `ddt()` term, so
    # the gap between their ratios is the cost of OSDI's reactive path with the
    # model's arithmetic held constant. That contrast exists because the first run
    # of this benchmark showed the *linear* model at 15x and the nonlinear one —
    # which does strictly more work — at 6x, which cannot be the arithmetic.
    "resistive": {
        "native": lambda n, _: ladder(n, "R{i}g {node} 0 1k"),
        "osdi": lambda n, osdi: ladder(
            n,
            "X{i} {node} 0 osdi_g_shunt gd=1m",
            extra_cards=f".osdi {osdi['osdi_g_shunt']}",
        ),
        "va": "osdi_g_shunt",
        "probe": "n1",
    },
    "reactive": {
        # 1 mS in parallel with 1 nF, both ways.
        "native": lambda n, _: ladder(
            n, "R{i}g {node} 0 1k\nC{i}g {node} 0 1n"
        ),
        "osdi": lambda n, osdi: ladder(
            n,
            "X{i} {node} 0 rc_shunt gd=1m c=1n",
            extra_cards=f".osdi {osdi['rc_shunt']}",
        ),
        "va": "rc_shunt",
        "probe": "n1",
    },
    "nonlinear_exp": {
        "native": lambda n, _: ladder(
            n,
            "D{i} {node} 0 dm",
            extra_cards=".model dm D (IS=1e-14 N=1)",
        ),
        "osdi": lambda n, osdi: ladder(
            n,
            "X{i} {node} 0 osdi_shockley is=1e-14 n=1",
            extra_cards=f".osdi {osdi['osdi_shockley']}",
        ),
        "va": "osdi_shockley",
        "probe": "n1",
    },
}

TRAN = ".tran 1u 200u\n"


# ---------------------------------------------------------------------------
# Running
# ---------------------------------------------------------------------------


def run(deck_text, tmp, tag, want_csv=False):
    """Run one deck. Returns (seconds, stdout) or (None, error)."""
    path = tmp / f"{tag}.sp"
    path.write_text(deck_text)
    argv = [str(FC), "-f", str(path)]
    if not want_csv:
        argv += ["-o", os.devnull]
    t0 = time.perf_counter()
    r = subprocess.run(argv, capture_output=True, text=True)
    dt = time.perf_counter() - t0
    if r.returncode != 0:
        first = (r.stderr.strip().splitlines() or ["?"])[0]
        return None, first
    return dt, r.stdout


def probe_value(stdout, column):
    """Last row's value for `column`, for the agreement check."""
    lines = stdout.strip().splitlines()
    if len(lines) < 2:
        return None
    hdr = [h.strip() for h in lines[0].split(",")]
    want = f"V({column})"
    if want not in hdr:
        return None
    return float(lines[-1].split(",")[hdr.index(want)])


def check_agreement(tmp, osdi, results):
    """Both sides must compute the same circuit before any timing is quoted."""
    print("Agreement check (correctness gates the measurement)", file=sys.stderr)
    ok = True
    for name, spec in DECKS.items():
        vals = {}
        for side in ("native", "osdi"):
            deck = spec[side](8, osdi)
            dt, out = run(deck, tmp, f"agree_{name}_{side}", want_csv=True)
            if dt is None:
                print(f"  {name}/{side}: FAILED — {out}", file=sys.stderr)
                ok = False
                vals[side] = None
                continue
            vals[side] = probe_value(out, spec["probe"])
        a, b = vals.get("native"), vals.get("osdi")
        if a is None or b is None:
            ok = False
            continue
        rel = abs(a - b) / max(abs(a), abs(b), 1e-30)
        results["agreement"][name] = {"native": a, "osdi": b, "rel": rel}
        verdict = "ok" if rel < 1e-6 else "DISAGREE"
        print(
            f"  {name:10} native={a:.9e}  osdi={b:.9e}  rel={rel:.2e}  {verdict}",
            file=sys.stderr,
        )
        if rel >= 1e-6:
            ok = False
    return ok


def measure(tmp, osdi, reps, results):
    print(f"\nTiming, {reps} interleaved reps", file=sys.stderr)
    for analysis, label in ((".op\n", "dc"), (TRAN, "tran")):
        for name, spec in DECKS.items():
            for n in COUNTS:
                samples = {"native": [], "osdi": []}
                fails = {}
                for _ in range(reps):
                    for side in ("native", "osdi"):
                        deck = spec[side](n, osdi)
                        # Swap in the analysis card.
                        deck = deck.replace(".op\n", analysis)
                        dt, err = run(deck, tmp, f"t_{name}_{side}_{n}")
                        if dt is None:
                            fails[side] = err
                        else:
                            samples[side].append(dt)
                key = f"{label}/{name}/{n}"
                if fails or not samples["native"] or not samples["osdi"]:
                    results["timing"][key] = {"failed": fails or "no samples"}
                    print(f"  {key:24} FAILED {fails}", file=sys.stderr)
                    continue
                nat = statistics.median(samples["native"])
                osd = statistics.median(samples["osdi"])
                results["timing"][key] = {
                    "native_s": nat,
                    "osdi_s": osd,
                    "ratio": osd / nat,
                }
                print(
                    f"  {key:24} native={nat * 1e3:8.2f} ms  "
                    f"osdi={osd * 1e3:8.2f} ms  ratio={osd / nat:5.2f}x",
                    file=sys.stderr,
                )


def plot(results, out_path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1.15], hspace=0.42, wspace=0.28)
    fig.suptitle(
        "OSDI (compiled Verilog-A) against native Rust models, identical arithmetic",
        fontsize=14,
    )

    for row, label in enumerate(("dc", "tran")):
        for col, name in enumerate(DECKS):
            ax = fig.add_subplot(gs[row, col])
            xs, nat, osd = [], [], []
            for n in COUNTS:
                e = results["timing"].get(f"{label}/{name}/{n}", {})
                if "ratio" not in e:
                    continue
                xs.append(n)
                nat.append(e["native_s"] * 1e3)
                osd.append(e["osdi_s"] * 1e3)
            if not xs:
                ax.text(0.5, 0.5, "no data", ha="center", transform=ax.transAxes)
                ax.set_title(f"{label} — {name}")
                continue
            ax.plot(xs, nat, "o-", label="native", color="#1f77b4")
            ax.plot(xs, osd, "s-", label="OSDI", color="#d62728")
            ax.set_xscale("log", base=2)
            ax.set_yscale("log")
            ax.set_xlabel("devices in the ladder")
            ax.set_ylabel("wall clock (ms)")
            ratios = [o / m for o, m in zip(osd, nat)]
            ax.set_title(f"{label} — {name}   {min(ratios):.2f}–{max(ratios):.2f}×")
            ax.grid(True, which="both", alpha=0.3)
            ax.legend(fontsize=8)
            for x, o, m in zip(xs, osd, nat):
                ax.annotate(
                    f"{o / m:.2f}×",
                    (x, o),
                    textcoords="offset points",
                    xytext=(0, 7),
                    ha="center",
                    fontsize=8,
                    color="#d62728",
                )

    # The panel the ratios cannot show. A ratio is relative to the native side, so
    # `nonlinear_exp` looks *better* than `resistive` only because a native diode
    # is slower than a native resistor. Absolute extra time per device removes
    # that, and the three curves then say what the overhead is made of.
    ax = fig.add_subplot(gs[2, :])
    styles = {
        "resistive": ("#1f77b4", "o"),
        "reactive": ("#d62728", "s"),
        "nonlinear_exp": ("#2ca02c", "^"),
    }
    for name in DECKS:
        for label, dash in (("dc", ":"), ("tran", "-")):
            xs, per = [], []
            for n in COUNTS:
                e = results["timing"].get(f"{label}/{name}/{n}", {})
                if "ratio" not in e:
                    continue
                extra = (e["osdi_s"] - e["native_s"]) * 1e6  # microseconds
                xs.append(n)
                per.append(extra / n)
            if not xs:
                continue
            colour, marker = styles.get(name, ("#888", "x"))
            ax.plot(
                xs,
                per,
                dash,
                marker=marker,
                color=colour,
                label=f"{name} ({label})",
                alpha=1.0 if label == "tran" else 0.55,
            )
    ax.set_xscale("log", base=2)
    ax.set_xlabel("devices in the ladder")
    ax.set_ylabel("OSDI overhead (µs per device)")
    ax.set_title(
        "What the overhead is made of — absolute extra time per device, not a ratio"
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, ncol=3)
    ax.axhline(0.0, color="#444", lw=0.8)
    ax.text(
        0.012,
        0.94,
        "Below ~128 devices process startup dominates and these numbers are noise.\n"
        "At 2048, transient: resistive 232 µs/device and nonlinear_exp 231 — identical, so\n"
        "the overhead is the per-eval ABI call, not the model's arithmetic. reactive 594:\n"
        "one ddt() term more than doubles it (a second pair of ABI calls every step).\n"
        "But it RISES with N rather than flattening, so it is not a constant per device:\n"
        "512→2048 the OSDI side grows as N^1.5–1.7 where the native side grows as N^0.7–0.9.\n"
        "That gap is unexplained and is a lead, not a conclusion.",
        transform=ax.transAxes,
        fontsize=9,
        va="top",
        bbox={"boxstyle": "round", "fc": "#fffbe6", "ec": "#ccc"},
    )

    agree = results.get("agreement", {})
    note = "  ".join(f"{k}: rel {v['rel']:.1e}" for k, v in agree.items() if "rel" in v)
    fig.text(
        0.5,
        0.005,
        f"Agreement checked before any timing — {note}   |   interleaved A/B, "
        f"median of {results['reps']} reps, default options, benchmarks/METHODOLOGY.md",
        ha="center",
        fontsize=8,
    )
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    print(f"\nwrote {out_path}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="osdi_vs_native.png")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--json")
    args = ap.parse_args()

    if not FC.exists():
        die(f"{FC} missing — run `cargo build --release`")

    results = {"reps": args.reps, "counts": COUNTS, "agreement": {}, "timing": {}}
    with tempfile.TemporaryDirectory(prefix="fc_osdi_bench_") as td:
        tmp = Path(td)
        cache = tmp / "va"
        cache.mkdir()
        osdi = {}
        for spec in DECKS.values():
            stem = spec["va"]
            src = MODELS / f"{stem}.va"
            if not src.exists():
                src = OSDI_MODELS / f"{stem}.va"
            if not src.exists():
                die(f"cannot find {stem}.va")
            osdi[stem] = str(compile_va(src, cache))

        if not check_agreement(tmp, osdi, results):
            print(
                "\nrefusing to report timings: the two sides do not agree, so any "
                "ratio would be a comparison of different circuits",
                file=sys.stderr,
            )
            sys.exit(1)
        measure(tmp, osdi, args.reps, results)

    if args.json:
        Path(args.json).write_text(json.dumps(results, indent=2))
    plot(results, args.out)


if __name__ == "__main__":
    main()
