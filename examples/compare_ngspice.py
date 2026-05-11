#!/usr/bin/env python3
"""
Compare fairchild vs ngspice on several example circuits.

Requirements:
    pip install matplotlib numpy
    ngspice must be on PATH
    cargo build --release (or debug) must have been run

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

REPO_ROOT = Path(__file__).parent.parent
EXAMPLES_DIR = REPO_ROOT / "examples"
PLOTS_DIR = REPO_ROOT / "docs" / "plots"


def find_fairchild(release: bool) -> Path:
    profile = "release" if release else "debug"
    p = REPO_ROOT / "target" / profile / "fairchild"
    if not p.exists():
        sys.exit(f"fairchild binary not found at {p}. Run: cargo build {'--release' if release else ''}")
    return p


def find_ngspice() -> str:
    for candidate in ["ngspice", "/opt/homebrew/bin/ngspice", "/usr/local/bin/ngspice", "/usr/bin/ngspice"]:
        if shutil.which(candidate) or Path(candidate).exists():
            return candidate
    sys.exit("ngspice not found on PATH")


def run_fairchild(binary: Path, netlist: Path) -> dict[str, list[float]]:
    """Run fairchild and return dict of column_name → [values]."""
    result = subprocess.run(
        [str(binary), "-f", str(netlist)],
        capture_output=True, text=True, check=True,
    )
    reader = csv.DictReader(io.StringIO(result.stdout))
    data: dict[str, list[float]] = {}
    for row in reader:
        for k, v in row.items():
            data.setdefault(k, []).append(float(v))
    return data


def run_ngspice_tran(netlist: Path, time_col: str = "time") -> dict[str, list[float]]:
    """Run ngspice batch on a .tran netlist with all-node print and return data."""
    src = netlist.read_text()
    # Strip .end; inject .control block
    stripped = "\n".join(
        l for l in src.splitlines()
        if l.strip().lower() not in (".end",)
    )
    control = ".control\ntran\nwrite /dev/stdout\n.endc\n.end\n"
    # Use print instead of write for text output
    control = ".control\ntran\nprint all\n.endc\n.end\n"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".sp", delete=False) as f:
        f.write(stripped + "\n" + control)
        tmp = f.name

    try:
        result = subprocess.run(
            ["ngspice", "-b", tmp],
            capture_output=True, text=True,
        )
        return _parse_ngspice_print(result.stdout + result.stderr)
    finally:
        os.unlink(tmp)


def run_ngspice_op(netlist: Path) -> dict[str, float]:
    """Run ngspice .op and return node voltages."""
    src = netlist.read_text()
    stripped = "\n".join(l for l in src.splitlines() if l.strip().lower() not in (".end",))
    control = ".control\nop\nprint all\n.endc\n.end\n"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".sp", delete=False) as f:
        f.write(stripped + "\n" + control)
        tmp = f.name

    try:
        result = subprocess.run(["ngspice", "-b", tmp], capture_output=True, text=True)
        parsed = _parse_ngspice_print(result.stdout + result.stderr)
        return {k: v[0] for k, v in parsed.items() if v}
    finally:
        os.unlink(tmp)


def _parse_ngspice_print(output: str) -> dict[str, list[float]]:
    """Parse ngspice `print all` output into column → [values]."""
    data: dict[str, list[float]] = {}
    for line in output.splitlines():
        line = line.strip()
        if "=" in line:
            lhs, rhs = line.split("=", 1)
            key = lhs.strip().lower()
            if key.startswith(("v(", "i(", "time")):
                try:
                    val = float(rhs.strip().split()[0])
                    data.setdefault(key, []).append(val)
                except (ValueError, IndexError):
                    pass
    return data


def plot_tran_comparison(fc_data: dict, ng_data: dict, node: str,
                          title: str, ylabel: str, out_path: Path):
    fig, ax = plt.subplots(figsize=(8, 4))

    fc_t = np.array(fc_data.get("time", []))
    fc_v = np.array(fc_data.get(node, []))

    ng_t = np.array(ng_data.get("time", []))
    ng_v = np.array(ng_data.get(node, []))

    if fc_t.size:
        ax.plot(fc_t * 1e3, fc_v, label="fairchild", linewidth=2)
    if ng_t.size:
        ax.plot(ng_t * 1e3, ng_v, "--", label="ngspice", linewidth=1.5, alpha=0.8)

    ax.set_xlabel("Time (ms)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  saved {out_path.relative_to(REPO_ROOT)}")


EXAMPLES = [
    {
        "netlist": "rc_step.sp",
        "type": "tran",
        "node": "v(out)",
        "title": "RC Step Response (τ = 1 ms)",
        "ylabel": "Voltage (V)",
        "plot": "rc_step_comparison.png",
        "time_scale": 1e3,
        "time_unit": "ms",
    },
    {
        "netlist": "rlc_resonator.sp",
        "type": "tran",
        "node": "v(n2)",
        "title": "RLC Resonator Step Response (f₀ ≈ 5 kHz)",
        "ylabel": "Voltage (V)",
        "plot": "rlc_resonator_comparison.png",
        "time_scale": 1e6,
        "time_unit": "µs",
    },
    {
        "netlist": "cmos_inverter.sp",
        "type": "tran",
        "node": "v(out)",
        "title": "CMOS Inverter Transient",
        "ylabel": "Voltage (V)",
        "plot": "cmos_inverter_comparison.png",
        "time_scale": 1e9,
        "time_unit": "ns",
    },
]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--release", action="store_true", help="Use release build of fairchild")
    args = ap.parse_args()

    binary = find_fairchild(args.release)
    ngspice = find_ngspice()
    print(f"fairchild: {binary}")
    print(f"ngspice:   {ngspice}")
    print()

    for ex in EXAMPLES:
        netlist = EXAMPLES_DIR / ex["netlist"]
        print(f"Running {ex['netlist']} ...")

        fc_data = run_fairchild(binary, netlist)
        ng_data = run_ngspice_tran(netlist)

        node = ex["node"]
        scale = ex["time_scale"]
        unit = ex["time_unit"]

        fig, ax = plt.subplots(figsize=(8, 4))
        fc_t = np.array(fc_data.get("time", []))
        fc_v = np.array(fc_data.get(node, []))
        ng_t = np.array(ng_data.get("time", []))
        ng_v = np.array(ng_data.get(node, []))

        if fc_t.size:
            ax.plot(fc_t * scale, fc_v, label="fairchild", linewidth=2)
        if ng_t.size:
            ax.plot(ng_t * scale, ng_v, "--", label="ngspice", linewidth=1.5, alpha=0.8)
        else:
            print(f"  warning: ngspice returned no data for {node}")

        ax.set_xlabel(f"Time ({unit})")
        ax.set_ylabel(ex["ylabel"])
        ax.set_title(ex["title"])
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()

        out_path = PLOTS_DIR / ex["plot"]
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"  saved {out_path.relative_to(REPO_ROOT)}")

    print("\nDone. Plots saved to docs/plots/")


if __name__ == "__main__":
    main()
