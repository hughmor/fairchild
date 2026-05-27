#!/usr/bin/env python3
"""
Benchmark fairchild vs ngspice: accuracy and wall-clock time.

Outputs a JSON blob to stdout (or --output FILE) and a summary table to stderr.

Requirements:
    ngspice on PATH (optional — omitted metrics are null)
    fairchild built: cargo build --release

Usage:
    python benchmarks/run_all.py [--release] [--output results.json]
"""

import argparse
import csv
import io
import json
import os
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
CIRCUITS   = Path(__file__).parent / "circuits"
FC_BIN_RELEASE = REPO_ROOT / "target" / "release" / "fairchild"
FC_BIN_DEBUG   = REPO_ROOT / "target" / "debug"  / "fairchild"

# Circuits to benchmark.  Each entry: (label, filename, probe_node, t_sample)
# t_sample: time (s) at which to compare fairchild vs ngspice voltage, or None
# for DC-only circuits.
BENCHMARKS = [
    ("RC step response",       "rc_step.sp",       "out",  2e-3),
    ("RLC resonator",          "rlc_resonator.sp", "n2",   0.5e-3),
    ("Diode rectifier",        "diode_rectifier.sp","out", 2e-6),
    ("CMOS inverter",          "cmos_inverter.sp", "out",  60e-9),
    ("BJT CE amplifier",       "bjt_ce_amp.sp",    "c",    100e-9),
    ("Ring osc 3-stage",       "ring_osc_3.sp",    "n1",   None),
    ("Ring osc 11-stage",      "ring_osc_11.sp",   "n1",   None),
]


def find_binary(name: str, release: bool) -> str:
    if name == "fairchild":
        p = FC_BIN_RELEASE if release else FC_BIN_DEBUG
        if not p.exists():
            sys.exit(f"fairchild not found at {p}. Run: cargo build {'--release' if release else ''}")
        return str(p)
    candidate = shutil.which(name)
    if candidate:
        return candidate
    for path in [f"/opt/homebrew/bin/{name}", f"/usr/local/bin/{name}", f"/usr/bin/{name}"]:
        if Path(path).exists():
            return path
    return None  # optional tools


def wall_ms(cmd: list[str], **kwargs) -> tuple[float, subprocess.CompletedProcess]:
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    return (time.perf_counter() - t0) * 1000, proc


def parse_csv_node(csv_text: str, node: str, t_sample: float) -> float | None:
    """Find voltage of `node` at the time-row closest to t_sample."""
    reader = csv.DictReader(io.StringIO(csv_text))
    col = f"V({node})" if f"V({node})" in (reader.fieldnames or []) else None
    if col is None:
        # Try lowercase
        for fn in (reader.fieldnames or []):
            if fn.lower() == f"v({node.lower()})":
                col = fn
                break
    if col is None:
        return None
    best_t, best_v = None, None
    for row in reader:
        try:
            t = float(list(row.values())[0])
            v = float(row[col])
        except (ValueError, KeyError):
            continue
        if best_t is None or abs(t - t_sample) < abs(best_t - t_sample):
            best_t, best_v = t, v
    return best_v


def ngspice_sample(netlist_path: Path, node: str, t_sample: float, ng_bin: str) -> float | None:
    """Run ngspice and extract V(node) at t_sample via .meas."""
    netlist = netlist_path.read_text()
    # Strip .end, inject .meas and .control
    lines = [l for l in netlist.splitlines() if l.strip().lower() not in (".end",)]
    meas = f".meas tran vsample FIND V({node}) AT={t_sample:.6e}"
    ctrl = f".control\nrun\nprint vsample\n.endc\n.end"
    full = "\n".join(lines) + f"\n{meas}\n{ctrl}\n"

    with tempfile.NamedTemporaryFile(suffix=".sp", mode="w", delete=False) as f:
        f.write(full)
        tmp = f.name
    try:
        proc = subprocess.run([ng_bin, "-b", tmp], capture_output=True, text=True)
        combined = proc.stdout + proc.stderr
        for line in combined.splitlines():
            if "vsample" in line.lower() and "=" in line:
                try:
                    return float(line.split("=")[-1].strip().split()[0])
                except ValueError:
                    pass
    finally:
        os.unlink(tmp)
    return None


def run_benchmark(label: str, spice_file: str, node: str | None, t_sample: float | None,
                  fc_bin: str, ng_bin: str | None) -> dict:
    circuit = CIRCUITS / spice_file
    result = {"label": label, "circuit": spice_file}

    # fairchild
    fc_ms, fc_proc = wall_ms([fc_bin, "-f", str(circuit)])
    result["fairchild_ms"] = round(fc_ms, 1)
    fc_ok = fc_proc.returncode == 0
    result["fairchild_ok"] = fc_ok
    if not fc_ok:
        result["error"] = fc_proc.stderr.strip().splitlines()[-1] if fc_proc.stderr else "unknown"

    # Sample fairchild output at t_sample for accuracy comparison
    fc_v = None
    if fc_ok and node and t_sample is not None:
        fc_v = parse_csv_node(fc_proc.stdout, node, t_sample)
    result["fairchild_v"] = fc_v

    # ngspice
    if ng_bin:
        ng_ms, ng_proc = wall_ms([ng_bin, "-b", str(circuit)])
        result["ngspice_ms"] = round(ng_ms, 1)
        result["ngspice_ok"] = ng_proc.returncode == 0

        if fc_ok and ng_proc.returncode == 0 and node and t_sample is not None:
            ng_v = ngspice_sample(circuit, node, t_sample, ng_bin)
            result["ngspice_v"] = ng_v
            if fc_v is not None and ng_v is not None and ng_v != 0:
                result["rel_error"] = abs(fc_v - ng_v) / abs(ng_v)

    return result


def print_table(results: list[dict]) -> None:
    hdr = f"{'Circuit':<25} {'FC (ms)':>9} {'NG (ms)':>9} {'FC V':>12} {'NG V':>12} {'rel err':>9}"
    print(hdr, file=sys.stderr)
    print("-" * len(hdr), file=sys.stderr)
    for r in results:
        fc_v = f"{r['fairchild_v']:.4e}" if r.get("fairchild_v") is not None else "—"
        ng_v = f"{r.get('ngspice_v', None):.4e}" if r.get("ngspice_v") is not None else "—"
        rel  = f"{r['rel_error']:.2e}" if r.get("rel_error") is not None else "—"
        ng_ms = f"{r['ngspice_ms']:.0f}" if r.get("ngspice_ms") is not None else "—"
        ok = "✓" if r.get("fairchild_ok") else "✗"
        print(
            f"{ok} {r['label']:<23} {r['fairchild_ms']:>9.1f} {ng_ms:>9} {fc_v:>12} {ng_v:>12} {rel:>9}",
            file=sys.stderr
        )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", action="store_true", default=True)
    ap.add_argument("--debug",   action="store_true")
    ap.add_argument("--output",  default=None, help="JSON output file (default: stdout)")
    args = ap.parse_args()

    release = not args.debug
    fc_bin  = find_binary("fairchild", release)
    ng_bin  = find_binary("ngspice",  release)
    if ng_bin is None:
        print("ngspice not found — accuracy comparison will be skipped", file=sys.stderr)

    results = []
    for label, spice, node, t_sample in BENCHMARKS:
        print(f"  {label}...", end=" ", flush=True, file=sys.stderr)
        r = run_benchmark(label, spice, node, t_sample, fc_bin, ng_bin)
        results.append(r)
        status = "ok" if r.get("fairchild_ok") else f"FAILED: {r.get('error','?')}"
        print(status, file=sys.stderr)

    print_table(results)

    blob = json.dumps({"results": results}, indent=2)
    if args.output:
        Path(args.output).write_text(blob)
        print(f"\nResults written to {args.output}", file=sys.stderr)
    else:
        print(blob)


if __name__ == "__main__":
    main()
