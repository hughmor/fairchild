#!/usr/bin/env python3
"""
Benchmark fairchild vs ngspice: wall-clock time and peak RSS memory.

Requirements:
    pip install tabulate
    ngspice must be on PATH
    cargo build --release must have been run

Usage:
    python scripts/benchmark.py [--runs N]

Output: Markdown table suitable for pasting into README.md
"""

import argparse
import os
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent.parent  # electronic/ -> examples/ -> repo root
EXAMPLES_DIR = Path(__file__).parent             # netlists live alongside this script

NETLISTS = [
    ("RC step (1k/1µF, 5ms tran)",   "rc_step.sp"),
    ("RLC resonator (1ms tran)",       "rlc_resonator.sp"),
    ("Diode rectifier (3µs tran)",     "diode_rectifier.sp"),
    ("CMOS inverter (120ns tran)",     "cmos_inverter.sp"),
    ("NMOS DC op",                     "nmos_dc_sweep.sp"),
]


def find_binary(name: str, release: bool = True) -> str:
    if name == "fairchild":
        profile = "release" if release else "debug"
        p = REPO_ROOT / "target" / profile / "fairchild"
        if not p.exists():
            sys.exit(f"fairchild not found at {p}. Run: cargo build --release")
        return str(p)
    candidate = shutil.which(name)
    if not candidate:
        for path in [f"/opt/homebrew/bin/{name}", f"/usr/local/bin/{name}"]:
            if Path(path).exists():
                return path
        sys.exit(f"{name} not found on PATH")
    return candidate


def time_and_rss(cmd: list[str], stdin_text: str | None = None) -> tuple[float, int]:
    """Return (wall_seconds, peak_rss_kb)."""
    start = time.perf_counter()
    proc = subprocess.run(
        cmd,
        input=stdin_text,
        capture_output=True,
        text=True,
    )
    elapsed = time.perf_counter() - start
    # macOS: getrusage MAXRSS is in bytes; Linux: kilobytes
    rss_bytes = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if sys.platform == "darwin":
        rss_kb = rss_bytes // 1024
    else:
        rss_kb = rss_bytes
    if proc.returncode != 0:
        print(f"  warning: exit code {proc.returncode}", file=sys.stderr)
        print(proc.stderr[:200], file=sys.stderr)
    return elapsed, rss_kb


def ngspice_cmd_for(netlist: Path) -> list[str]:
    src = netlist.read_text()
    stripped = "\n".join(l for l in src.splitlines() if l.strip().lower() not in (".end",))
    has_tran = any(".tran" in l.lower() for l in src.splitlines())
    has_op = any(l.strip().lower() == ".op" for l in src.splitlines())
    if has_tran:
        ctrl = ".control\ntran\nquit\n.endc\n.end\n"
    else:
        ctrl = ".control\nop\nquit\n.endc\n.end\n"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".sp", delete=False) as f:
        f.write(stripped + "\n" + ctrl)
        return ["ngspice", "-b", f.name]


def fmt_time(s: float) -> str:
    if s < 1e-3:
        return f"{s*1e6:.0f} µs"
    if s < 1:
        return f"{s*1e3:.1f} ms"
    return f"{s:.2f} s"


def fmt_rss(kb: int) -> str:
    if kb < 1024:
        return f"{kb} KB"
    return f"{kb/1024:.1f} MB"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs", type=int, default=5, help="Repetitions per benchmark")
    ap.add_argument("--no-release", action="store_true", help="Use debug build of fairchild")
    args = ap.parse_args()

    fairchild = find_binary("fairchild", release=not args.no_release)
    ngspice = find_binary("ngspice")

    print(f"fairchild: {fairchild}")
    print(f"ngspice:   {ngspice}")
    print(f"runs per benchmark: {args.runs}\n")

    rows = []
    for label, filename in NETLISTS:
        netlist = EXAMPLES_DIR / filename
        if not netlist.exists():
            print(f"skip {filename}: not found")
            continue

        fc_times, fc_rsses = [], []
        ng_times, ng_rsses = [], []

        ng_cmd = ngspice_cmd_for(netlist)

        for _ in range(args.runs):
            t, r = time_and_rss([fairchild, "-f", str(netlist)])
            fc_times.append(t)
            fc_rsses.append(r)

        for _ in range(args.runs):
            t, r = time_and_rss(ng_cmd)
            ng_times.append(t)
            ng_rsses.append(r)

        fc_t = min(fc_times)
        ng_t = min(ng_times)
        fc_r = min(fc_rsses)
        ng_r = min(ng_rsses)
        speedup = ng_t / fc_t if fc_t > 0 else float("inf")

        rows.append((label, fmt_time(fc_t), fmt_time(ng_t), f"{speedup:.1f}×",
                     fmt_rss(fc_r), fmt_rss(ng_r)))
        print(f"  {label}: fairchild {fmt_time(fc_t)}, ngspice {fmt_time(ng_t)}, speedup {speedup:.1f}×")

    # Print markdown table
    print("\n## Benchmark Results\n")
    header = ("Circuit", "fairchild", "ngspice", "Speedup", "fairchild RSS", "ngspice RSS")
    sep = ("---",) * len(header)
    print("| " + " | ".join(header) + " |")
    print("| " + " | ".join(sep) + " |")
    for row in rows:
        print("| " + " | ".join(row) + " |")


if __name__ == "__main__":
    main()
