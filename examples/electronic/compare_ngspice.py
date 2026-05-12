#!/usr/bin/env python3
"""
Compare fairchild vs ngspice on several example circuits.

Requirements:
    pip install matplotlib numpy
    ngspice must be on PATH
    cargo build (or --release) must have been run

Usage:
    python examples/compare_ngspice.py [--release]
"""

import argparse
import csv
import io
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless — save to file
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).parent.parent.parent  # electronic/ -> examples/ -> repo root
EXAMPLES_DIR = Path(__file__).parent             # netlists live alongside this script
PLOTS_DIR = REPO_ROOT / "docs" / "plots"

# SPICE SI suffix → multiplier
_SUFFIX = {"f": 1e-15, "p": 1e-12, "n": 1e-9, "u": 1e-6, "m": 1e-3,
           "k": 1e3, "meg": 1e6, "g": 1e9, "t": 1e12}


def spice_float(s: str) -> float:
    """Convert a SPICE number string like '50u' or '5m' to a float."""
    s = s.lower().strip()
    for suffix, mult in sorted(_SUFFIX.items(), key=lambda x: -len(x[0])):
        if s.endswith(suffix):
            return float(s[: -len(suffix)]) * mult
    return float(s)


def find_fairchild(release: bool) -> Path:
    profile = "release" if release else "debug"
    p = REPO_ROOT / "target" / profile / "fairchild"
    if not p.exists():
        sys.exit(f"fairchild binary not found at {p}.\nRun: cargo build {'--release' if release else ''}")
    return p


def find_ngspice() -> str:
    for candidate in ["ngspice", "/opt/homebrew/bin/ngspice", "/usr/local/bin/ngspice", "/usr/bin/ngspice"]:
        if shutil.which(candidate) or Path(candidate).exists():
            return candidate
    sys.exit("ngspice not found on PATH")


# ---------------------------------------------------------------------------
# fairchild runner
# ---------------------------------------------------------------------------

def run_fairchild(binary: Path, netlist: Path) -> dict[str, list[float]]:
    """Run fairchild and return dict of lowercase column_name → [values]."""
    result = subprocess.run(
        [str(binary), "-f", str(netlist)],
        capture_output=True, text=True, check=True,
    )
    reader = csv.DictReader(io.StringIO(result.stdout))
    data: dict[str, list[float]] = {}
    for row in reader:
        for k, v in row.items():
            try:
                data.setdefault(k.lower(), []).append(float(v))
            except ValueError:
                pass  # skip label columns like "analysis"
    return data


# ---------------------------------------------------------------------------
# ngspice runners
# ---------------------------------------------------------------------------

def _parse_ngspice_print_table(output: str) -> tuple[list[float], list[float]]:
    """Parse ngspice `print v(node)` column output.

    Expected format (after a header line containing "Index"):
        Index   time            v(out)
        ------  ----            ------
        0       0.000000e+00    0.000000e+00
        1       1.000000e-11    1.000000e-10
    Returns (time_list, value_list).
    """
    time_vals: list[float] = []
    node_vals: list[float] = []
    in_table = False
    for line in output.splitlines():
        stripped = line.strip()
        if not in_table:
            if "index" in stripped.lower() and "time" in stripped.lower():
                in_table = True
            continue
        if stripped.startswith("-"):  # separator line
            continue
        parts = stripped.split()
        if len(parts) >= 3 and parts[0].lstrip("-").isdigit():
            try:
                time_vals.append(float(parts[1]))
                node_vals.append(float(parts[2]))
            except ValueError:
                continue
    return time_vals, node_vals


def run_ngspice_tran(netlist: Path, node: str, step: float, stop: float) -> dict[str, list[float]]:
    """Run ngspice transient and return {time: [...], node: [...]}.

    We strip .tran and .end from the netlist and inject a .control block
    so we can pass step/stop explicitly and use `print` for clean output.
    """
    src = netlist.read_text()
    lines = [
        l for l in src.splitlines()
        if not l.strip().lower().startswith(".tran")
        and l.strip().lower() != ".end"
    ]
    body = "\n".join(lines)
    # Format step/stop as scientific notation to avoid unit ambiguity
    control = (
        f".control\n"
        f"tran {step:.6e} {stop:.6e}\n"
        f"print {node}\n"
        f".endc\n"
        f".end\n"
    )

    with tempfile.NamedTemporaryFile(mode="w", suffix=".sp", delete=False) as f:
        f.write(body + "\n" + control)
        tmp = f.name

    try:
        result = subprocess.run(
            [find_ngspice(), "-b", tmp],
            capture_output=True, text=True,
        )
        combined = result.stdout + result.stderr
        if result.returncode != 0 and "aborted" in combined:
            print(f"  [ngspice error] {combined[:200]}", file=sys.stderr)
        t, v = _parse_ngspice_print_table(combined)
        return {"time": t, node: v}
    finally:
        os.unlink(tmp)


def run_ngspice_op(netlist: Path, nodes: list[str]) -> dict[str, float]:
    """Run ngspice .op and return {node: voltage}."""
    src = netlist.read_text()
    lines = [
        l for l in src.splitlines()
        if not l.strip().lower().startswith(".tran")
        and l.strip().lower() != ".end"
    ]
    body = "\n".join(lines)
    node_list = " ".join(nodes)
    control = f".control\nop\nprint {node_list}\n.endc\n.end\n"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".sp", delete=False) as f:
        f.write(body + "\n" + control)
        tmp = f.name

    try:
        result = subprocess.run(
            [find_ngspice(), "-b", tmp],
            capture_output=True, text=True,
        )
        values: dict[str, float] = {}
        for line in (result.stdout + result.stderr).splitlines():
            line = line.strip()
            if "=" in line:
                lhs, rhs = line.split("=", 1)
                key = lhs.strip().lower()
                if any(key == n.lower() for n in nodes):
                    try:
                        values[key] = float(rhs.strip().split()[0])
                    except (ValueError, IndexError):
                        pass
        return values
    finally:
        os.unlink(tmp)


# ---------------------------------------------------------------------------
# Example definitions
# ---------------------------------------------------------------------------

TRAN_EXAMPLES = [
    {
        "netlist": "rc_step.sp",
        "node": "v(out)",
        "step": 50e-6,
        "stop": 5e-3,
        "title": "RC Step Response (τ = 1 ms)",
        "xlabel": "Time (ms)",
        "ylabel": "Voltage (V)",
        "xscale": 1e3,
        "plot": "rc_step_comparison.png",
    },
    {
        "netlist": "rlc_resonator.sp",
        "node": "v(n2)",
        "step": 10e-6,
        "stop": 1e-3,
        "title": "RLC Resonator Step Response (f₀ ≈ 5 kHz)",
        "xlabel": "Time (µs)",
        "ylabel": "Voltage (V)",
        "xscale": 1e6,
        "plot": "rlc_resonator_comparison.png",
    },
    {
        "netlist": "diode_rectifier.sp",
        "node": "v(out)",
        "step": 10e-9,
        "stop": 3e-6,
        "title": "Half-Wave Diode Rectifier (3 µs)",
        "xlabel": "Time (µs)",
        "ylabel": "Voltage (V)",
        "xscale": 1e6,
        "plot": "diode_rectifier_comparison.png",
    },
    {
        "netlist": "cmos_inverter.sp",
        "node": "v(out)",
        "step": 1e-9,
        "stop": 120e-9,
        "title": "CMOS Inverter Transient (Level 1 MOSFET)",
        "xlabel": "Time (ns)",
        "ylabel": "Voltage (V)",
        "xscale": 1e9,
        "plot": "cmos_inverter_comparison.png",
    },
]

DC_EXAMPLES = [
    {
        "netlist": "nmos_dc_sweep.sp",
        "nodes": ["v(d)", "v(g)", "v(vdd)"],
        "labels": ["V(d)", "V(g)", "V(vdd)"],
        "title": "NMOS DC Operating Point",
        "plot": "nmos_dc_comparison.png",
    },
]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--release", action="store_true", help="Use release build of fairchild")
    args = ap.parse_args()

    binary = find_fairchild(args.release)
    print(f"fairchild: {binary}")
    print(f"ngspice:   {find_ngspice()}")
    print()

    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    for ex in TRAN_EXAMPLES:
        netlist = EXAMPLES_DIR / ex["netlist"]
        node = ex["node"]
        print(f"Running {ex['netlist']} ({node}) ...")

        fc_data = run_fairchild(binary, netlist)
        ng_data = run_ngspice_tran(netlist, node, ex["step"], ex["stop"])

        fc_t = np.array(fc_data.get("time", []))
        fc_v = np.array(fc_data.get(node, []))
        ng_t = np.array(ng_data.get("time", []))
        ng_v = np.array(ng_data.get(node, []))

        print(f"  fairchild: {len(fc_t)} timepoints")
        print(f"  ngspice:   {len(ng_t)} timepoints")

        if fc_t.size == 0:
            print(f"  [warning] fairchild returned no data for {node}")
            print(f"  available keys: {list(fc_data.keys())}")
            continue
        if ng_t.size == 0:
            print(f"  [warning] ngspice returned no data for {node}")

        scale = ex["xscale"]
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(fc_t * scale, fc_v, label="fairchild", linewidth=2, color="#2196F3")
        if ng_t.size:
            ax.plot(ng_t * scale, ng_v, "--", label="ngspice", linewidth=1.5,
                    color="#F44336", alpha=0.85)

        ax.set_xlabel(ex["xlabel"])
        ax.set_ylabel(ex["ylabel"])
        ax.set_title(ex["title"])
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()

        out_path = PLOTS_DIR / ex["plot"]
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"  saved → {out_path.relative_to(REPO_ROOT)}")

    for ex in DC_EXAMPLES:
        netlist = EXAMPLES_DIR / ex["netlist"]
        nodes = ex["nodes"]
        labels = ex["labels"]
        print(f"Running {ex['netlist']} (DC op, nodes: {nodes}) ...")

        fc_data = run_fairchild(binary, netlist)
        ng_vals = run_ngspice_op(netlist, nodes)

        fc_vals = [fc_data.get(n, [np.nan])[0] for n in nodes]
        ng_vals_list = [ng_vals.get(n, np.nan) for n in nodes]

        x = np.arange(len(labels))
        width = 0.35
        fig, ax = plt.subplots(figsize=(7, 4))
        bars_fc = ax.bar(x - width / 2, fc_vals, width, label="fairchild", color="#2196F3")
        bars_ng = ax.bar(x + width / 2, ng_vals_list, width, label="ngspice",
                         color="#F44336", alpha=0.85)

        ax.set_ylabel("Voltage (V)")
        ax.set_title(ex["title"])
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend()
        ax.grid(True, axis="y", alpha=0.3)

        for bar in list(bars_fc) + list(bars_ng):
            h = bar.get_height()
            if not np.isnan(h):
                ax.annotate(f"{h:.3f}", xy=(bar.get_x() + bar.get_width() / 2, h),
                            xytext=(0, 3), textcoords="offset points",
                            ha="center", va="bottom", fontsize=8)

        fig.tight_layout()
        out_path = PLOTS_DIR / ex["plot"]
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"  saved → {out_path.relative_to(REPO_ROOT)}")

    print("\nDone. Plots saved to docs/plots/")


if __name__ == "__main__":
    main()
