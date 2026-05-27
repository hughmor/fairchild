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

# Circuits to benchmark.  Each entry: (label, filename, probe_node, t_sample, n_stages)
# t_sample: time (s) at which to compare fairchild vs ngspice voltage, or None
# n_stages: number of ring osc stages (for period scaling check), or None
BENCHMARKS = [
    ("RC step response",       "rc_step.sp",         "out",  2e-3,   None),
    ("RLC resonator",          "rlc_resonator.sp",   "n2",   0.5e-3, None),
    ("Diode rectifier",        "diode_rectifier.sp", "out",  2e-6,   None),
    ("CMOS inverter",          "cmos_inverter.sp",   "out",  60e-9,  None),
    ("BJT CE amplifier",       "bjt_ce_amp.sp",      "c",    100e-9, None),
    ("Schmitt trigger",        "schmitt_trigger.sp", "out",  0.2e-6, None),
    ("Ring osc 3-stage",       "ring_osc_3.sp",      "n1",   None,   3),
    ("Ring osc 5-stage",       "ring_osc_5.sp",      "n1",   None,   5),
    ("Ring osc 11-stage",      "ring_osc_11.sp",     "n1",   None,   11),
    ("Ring osc 21-stage",      "ring_osc_21.sp",     "n1",   None,   21),
    ("Ring osc 51-stage",      "ring_osc_51.sp",     "n1",   None,   51),
    ("Ring osc 101-stage",     "ring_osc_101.sp",    "n1",   None,   101),
    ("Ring osc 201-stage",     "ring_osc_201.sp",    "n1",   None,   201),
    ("Ring osc 499-stage",     "ring_osc_499.sp",    "n1",   None,   499),
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


def measure_period(csv_text: str, node: str) -> float | None:
    """Measure the oscillation period of `node` from zero-crossing detection.

    Uses mid-swing threshold (VDD/2 = 0.9V) with hysteresis. Returns period
    in seconds, or None if fewer than 2 full cycles are found.
    """
    reader = csv.DictReader(io.StringIO(csv_text))
    col = None
    for fn in (reader.fieldnames or []):
        if fn.lower() == f"v({node.lower()})":
            col = fn
            break
    if col is None:
        return None

    times, voltages = [], []
    for row in reader:
        try:
            times.append(float(list(row.values())[0]))
            voltages.append(float(row[col]))
        except (ValueError, KeyError):
            continue

    if len(times) < 10:
        return None

    vmax = max(voltages)
    vmin = min(voltages)
    swing = vmax - vmin
    if swing < 0.5:  # not oscillating
        return None
    threshold = vmin + 0.5 * swing

    # Find rising zero-crossings
    crossings = []
    for i in range(1, len(times)):
        if voltages[i - 1] < threshold <= voltages[i]:
            # Linear interpolation
            frac = (threshold - voltages[i - 1]) / (voltages[i] - voltages[i - 1])
            crossings.append(times[i - 1] + frac * (times[i] - times[i - 1]))

    if len(crossings) < 2:
        return None
    # Average period over all consecutive pairs
    periods = [crossings[i + 1] - crossings[i] for i in range(len(crossings) - 1)]
    return sum(periods) / len(periods)


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


def _ngspice_tmp_with_control(circuit: Path, node: str) -> str:
    """Return path to a temp .sp file with a .control block so ngspice runs in batch mode."""
    src = circuit.read_text()
    lines = [l for l in src.splitlines() if l.strip().lower() != ".end"]
    ctrl = f".control\nrun\nprint {node}\n.endc\n.end\n"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".sp", delete=False) as f:
        f.write("\n".join(lines) + "\n" + ctrl)
        return f.name


def run_benchmark(label: str, spice_file: str, node: str | None, t_sample: float | None,
                  n_stages: int | None, fc_bin: str, ng_bin: str | None,
                  n_timing_runs: int = 3) -> dict:
    circuit = CIRCUITS / spice_file
    result = {"label": label, "circuit": spice_file}
    if n_stages is not None:
        result["n_stages"] = n_stages

    # fairchild — median of n_timing_runs; capture stdout on first run for accuracy
    fc_times = []
    fc_stdout_first = None
    for i in range(n_timing_runs):
        fc_ms, fc_proc = wall_ms([fc_bin, "-f", str(circuit)])
        fc_times.append(fc_ms)
        if i == 0:
            fc_stdout_first = fc_proc.stdout
            fc_ok = fc_proc.returncode == 0
            if not fc_ok:
                result["error"] = fc_proc.stderr.strip().splitlines()[-1] if fc_proc.stderr else "unknown"

    fc_times.sort()
    result["fairchild_ms"] = round(fc_times[n_timing_runs // 2], 1)
    result["fairchild_ok"] = fc_ok

    # Sample fairchild output at t_sample for accuracy comparison
    fc_v = None
    if fc_ok and node and t_sample is not None:
        fc_v = parse_csv_node(fc_stdout_first, node, t_sample)
    result["fairchild_v"] = fc_v

    # Ring osc accuracy: measure period from waveform and verify oscillation
    if fc_ok and n_stages is not None and node and fc_stdout_first:
        period = measure_period(fc_stdout_first, node)
        result["fc_period_s"] = period
        if period is not None:
            result["oscillating"] = True
        else:
            result["oscillating"] = False
            result["accuracy_warning"] = "no oscillation detected — results may be incorrect"

    # ngspice — must inject a .control block so it actually runs in batch mode
    if ng_bin:
        ng_tmp = _ngspice_tmp_with_control(circuit, node or "v(0)")
        try:
            ng_times = []
            ng_ok = False
            for i in range(n_timing_runs):
                ng_ms, ng_proc = wall_ms([ng_bin, "-b", ng_tmp])
                ng_times.append(ng_ms)
                if i == 0:
                    ng_ok = ng_proc.returncode in (0, 1)  # ngspice exits 1 on some warnings
            ng_times.sort()
            result["ngspice_ms"] = round(ng_times[n_timing_runs // 2], 1)
            result["ngspice_ok"] = ng_ok

            if fc_ok and ng_ok and node and t_sample is not None:
                ng_v = ngspice_sample(circuit, node, t_sample, ng_bin)
                result["ngspice_v"] = ng_v
                if fc_v is not None and ng_v is not None and ng_v != 0:
                    result["rel_error"] = abs(fc_v - ng_v) / abs(ng_v)
        finally:
            try:
                os.unlink(ng_tmp)
            except OSError:
                pass

    return result


def print_table(results: list[dict]) -> None:
    hdr = f"{'Circuit':<25} {'FC (ms)':>9} {'NG (ms)':>9} {'FC V':>12} {'NG V':>12} {'rel err':>9} {'osc?':>6}"
    print(hdr, file=sys.stderr)
    print("-" * len(hdr), file=sys.stderr)
    for r in results:
        fc_v = f"{r['fairchild_v']:.4e}" if r.get("fairchild_v") is not None else "—"
        ng_v = f"{r.get('ngspice_v', None):.4e}" if r.get("ngspice_v") is not None else "—"
        rel  = f"{r['rel_error']:.2e}" if r.get("rel_error") is not None else "—"
        ng_ms = f"{r['ngspice_ms']:.0f}" if r.get("ngspice_ms") is not None else "—"
        ok = "✓" if r.get("fairchild_ok") else "✗"
        if "oscillating" in r:
            osc = "yes" if r["oscillating"] else "NO"
        else:
            osc = "—"
        print(
            f"{ok} {r['label']:<23} {r['fairchild_ms']:>9.1f} {ng_ms:>9} {fc_v:>12} {ng_v:>12} {rel:>9} {osc:>6}",
            file=sys.stderr
        )

    # Period scaling check: period must scale linearly with N (period ∝ N).
    # We compute ns/stage for each ring, then warn if any individual value
    # deviates by more than 50% from the median (indicating non-oscillation
    # or solver divergence, not a specific expected absolute value).
    ring_results = [r for r in results if r.get("n_stages") and r.get("fc_period_s") is not None]
    if len(ring_results) >= 2:
        ratios = [r["fc_period_s"] * 1e9 / r["n_stages"] for r in ring_results]
        ratios_sorted = sorted(ratios)
        median_ratio = ratios_sorted[len(ratios_sorted) // 2]
        print(f"\nPeriod scaling check (median {median_ratio:.3f} ns/stage; warn if >50% deviation):", file=sys.stderr)
        any_warn = False
        for r, ratio in zip(ring_results, ratios):
            n = r["n_stages"]
            p_ns = r["fc_period_s"] * 1e9
            dev = abs(ratio - median_ratio) / median_ratio if median_ratio > 0 else float("inf")
            flag = f"  ← WARN: {dev:.0%} deviation from median" if dev > 0.5 else ""
            if flag:
                any_warn = True
            print(f"  {r['label']:<25} N={n:>4}  period={p_ns:6.1f}ns  ns/stage={ratio:.3f}{flag}", file=sys.stderr)
        if any_warn:
            print("  SCALING LINEARITY FAILURE — check oscillation and solver correctness", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", action="store_true", default=True)
    ap.add_argument("--debug",   action="store_true")
    ap.add_argument("--output",  default=None, help="JSON output file (default: stdout)")
    ap.add_argument("--runs",    type=int, default=3, help="Timing runs per circuit (default: 3)")
    args = ap.parse_args()

    release = not args.debug
    fc_bin  = find_binary("fairchild", release)
    ng_bin  = find_binary("ngspice",  release)
    if ng_bin is None:
        print("ngspice not found — accuracy comparison will be skipped", file=sys.stderr)

    results = []
    for label, spice, node, t_sample, n_stages in BENCHMARKS:
        print(f"  {label}...", end=" ", flush=True, file=sys.stderr)
        r = run_benchmark(label, spice, node, t_sample, n_stages, fc_bin, ng_bin, args.runs)
        results.append(r)
        status = "ok" if r.get("fairchild_ok") else f"FAILED: {r.get('error','?')}"
        if r.get("accuracy_warning"):
            status += f"  WARNING: {r['accuracy_warning']}"
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
