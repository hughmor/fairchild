#!/usr/bin/env python3
"""
Generate benchmark plots comparing fairchild vs ngspice.

Outputs (to docs/plots/):
  accuracy_analog.png    — 6-panel waveform overlay with RMS error
  scaling_wall_time.png  — log-log wall-clock vs circuit size

Usage:
    python benchmarks/plot.py [--release] [--no-ngspice]

Methodology: the accuracy panels run fairchild and ngspice with default
options on identical netlists. The scaling plot forces each fairchild solver
backend in turn (--solver dense|sparse|klu) to compare them head-to-head, with
ngspice on its default solver. Timing is median of 3 end-to-end runs (including
binary load + parse). KLU requires a build with --features klu; it is probed at
runtime and skipped if absent. Dense LU is capped at DENSE_MAX_NODES nodes
(O(N^3)). Failed circuits appear as labelled gray panels rather than dropped.
"""

import argparse
import csv
import io
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from statistics import median

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).parent.parent
CIRCUITS   = Path(__file__).parent / "circuits"
PLOTS_DIR  = REPO_ROOT / "docs" / "plots"

FC_BIN_RELEASE = REPO_ROOT / "target" / "release" / "fairchild"
FC_BIN_DEBUG   = REPO_ROOT / "target" / "debug"   / "fairchild"

# Accuracy panel configuration: (title, netlist, fc_col, ng_node, x_label, x_scale)
ACCURACY_PANELS = [
    ("RC Step Response",     "rc_step.sp",        "V(out)", "v(out)", "Time (ms)",  1e3,  "Voltage (V)"),
    ("RLC Resonator",        "rlc_resonator.sp",  "V(n2)",  "v(n2)", "Time (µs)",  1e6,  "Voltage (V)"),
    ("Diode Rectifier",      "diode_rectifier.sp","V(out)", "v(out)", "Time (µs)",  1e6,  "Voltage (V)"),
    ("CMOS Inverter",        "cmos_inverter.sp",  "V(out)", "v(out)", "Time (ns)",  1e9,  "Voltage (V)"),
    ("BJT CE Amplifier",     "bjt_ce_amp.sp",     "V(c)",   "v(c)",  "Time (ns)",  1e9,  "Voltage (V)"),
    ("Schmitt Trigger",      "schmitt_trigger.sp","V(out)", "v(out)", "Time (µs)",  1e6,  "Voltage (V)"),
]

# Scaling circuits: (netlist, stages)
SCALING_CIRCUITS = [
    ("ring_osc_3.sp",    3),
    ("ring_osc_5.sp",    5),
    ("ring_osc_11.sp",   11),
    ("ring_osc_21.sp",   21),
    ("ring_osc_51.sp",   51),
    ("ring_osc_101.sp",  101),
    ("ring_osc_201.sp",  201),
    ("ring_osc_499.sp",  499),
]

N_TIMING_RUNS = 3   # median of this many runs for wall-clock timing

# Solver backends compared in the scaling plot. KLU requires the binary to be
# built with `--features klu`; it is probed at runtime and skipped if absent.
SCALING_BACKENDS = ["dense", "sparse", "klu"]
# Dense LU is O(N^3) per factorisation; cap it so the suite stays quick. Beyond
# this node count dense is skipped (and the cap is noted on the plot) — its
# steep curve up to the cap already shows the asymptotic disadvantage.
DENSE_MAX_NODES = 250
# (colour, marker, legend label) per backend.
BACKEND_STYLE = {
    "dense":  ("#9C27B0", "o", "fairchild — dense LU (faer)"),
    "sparse": ("#2196F3", "s", "fairchild — sparse LU (faer)"),
    "klu":    ("#FF9800", "^", "fairchild — KLU (SuiteSparse)"),
}

# Matplotlib style
BLUE   = "#2196F3"
RED    = "#E53935"
GRAY   = "#9E9E9E"
plt.rcParams.update({
    "font.size":        9,
    "axes.linewidth":   0.8,
    "grid.alpha":       0.25,
    "figure.dpi":       150,
})


# ---------------------------------------------------------------------------
# Simulator runners
# ---------------------------------------------------------------------------

def find_binary(name: str, release: bool) -> str | None:
    if name == "fairchild":
        p = FC_BIN_RELEASE if release else FC_BIN_DEBUG
        return str(p) if p.exists() else None
    for cand in [name, f"/opt/homebrew/bin/{name}", f"/usr/local/bin/{name}", f"/usr/bin/{name}"]:
        if shutil.which(cand) or Path(cand).exists():
            return cand
    return None


def run_fairchild(fc_bin: str, netlist: Path) -> tuple[dict[str, np.ndarray] | None, str]:
    try:
        proc = subprocess.run(
            [fc_bin, "-f", str(netlist)],
            capture_output=True, text=True, timeout=120,
        )
        if proc.returncode != 0:
            return None, proc.stderr.strip().splitlines()[-1] if proc.stderr else "non-zero exit"
        reader = csv.DictReader(io.StringIO(proc.stdout))
        data: dict[str, list[float]] = {}
        for row in reader:
            for k, v in row.items():
                try:
                    data.setdefault(k, []).append(float(v))
                except (ValueError, TypeError):
                    pass
        return {k: np.array(v) for k, v in data.items()}, ""
    except Exception as e:
        return None, str(e)


def _parse_ngspice_print(output: str) -> tuple[list[float], list[float]]:
    """Parse 'print v(node)' tabular output from ngspice batch mode."""
    times: list[float] = []
    vals:  list[float] = []
    in_table = False
    for line in output.splitlines():
        s = line.strip()
        if not in_table:
            if "index" in s.lower() and "time" in s.lower():
                in_table = True
            continue
        if s.startswith("-"):
            continue
        parts = s.split()
        if len(parts) >= 3 and parts[0].lstrip("-").isdigit():
            try:
                times.append(float(parts[1]))
                vals.append(float(parts[2]))
            except ValueError:
                pass
    return times, vals


def run_ngspice(ng_bin: str, netlist: Path, node: str) -> tuple[tuple[np.ndarray, np.ndarray] | None, str]:
    src = netlist.read_text()
    lines = [l for l in src.splitlines() if l.strip().lower() != ".end"]
    body = "\n".join(lines)
    control = f".control\nrun\nprint {node}\n.endc\n.end\n"
    tmp = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".sp", delete=False) as f:
            f.write(body + "\n" + control)
            tmp = f.name
        proc = subprocess.run(
            [ng_bin, "-b", tmp],
            capture_output=True, text=True, timeout=120,
        )
        if proc.returncode not in (0, 1):  # ngspice exits 1 on some warnings
            combined = proc.stdout + proc.stderr
            if "aborted" in combined.lower() or "error" in combined.lower():
                msg = next((l for l in combined.splitlines() if "error" in l.lower()), "failed")
                return None, msg
        t, v = _parse_ngspice_print(proc.stdout + proc.stderr)
        if not t:
            return None, "no waveform data in ngspice output"
        return (np.array(t), np.array(v)), ""
    except Exception as e:
        return None, str(e)
    finally:
        if tmp:
            try:
                os.unlink(tmp)
            except OSError:
                pass


def wall_ms(cmd: list[str]) -> float:
    t0 = time.perf_counter()
    subprocess.run(cmd, capture_output=True)
    return (time.perf_counter() - t0) * 1000


def time_median_ms(cmd: list[str], n: int = N_TIMING_RUNS) -> float:
    return median(wall_ms(cmd) for _ in range(n))


# ---------------------------------------------------------------------------
# RMS accuracy metric
# ---------------------------------------------------------------------------

def rms_error(t_ref: np.ndarray, v_ref: np.ndarray,
              t_test: np.ndarray, v_test: np.ndarray) -> float | None:
    if t_ref.size < 2 or t_test.size < 2:
        return None
    # Interpolate test onto reference timepoints within the overlapping window
    t_lo = max(t_ref[0], t_test[0])
    t_hi = min(t_ref[-1], t_test[-1])
    mask = (t_ref >= t_lo) & (t_ref <= t_hi)
    if mask.sum() < 2:
        return None
    t_common = t_ref[mask]
    v_ref_c  = v_ref[mask]
    v_test_i = np.interp(t_common, t_test, v_test)
    rms = float(np.sqrt(np.mean((v_ref_c - v_test_i) ** 2)))
    return rms


# ---------------------------------------------------------------------------
# Plot 1: accuracy panels
# ---------------------------------------------------------------------------

def plot_accuracy(fc_bin: str, ng_bin: str | None) -> Path:
    n_panels = len(ACCURACY_PANELS)
    ncols = 3
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(13, 3.5 * nrows))
    axes = axes.flatten()

    summary_rows = []

    for idx, (title, fname, fc_col, ng_node, xlabel, xscale, ylabel) in enumerate(ACCURACY_PANELS):
        ax = axes[idx]
        netlist = CIRCUITS / fname
        print(f"  {title} ...", end=" ", flush=True)

        # --- fairchild ---
        fc_data, fc_err = run_fairchild(fc_bin, netlist)
        fc_t = fc_data.get("time") if fc_data else None
        fc_v = fc_data.get(fc_col) if fc_data else None

        # --- ngspice ---
        ng_result, ng_err = (None, "ngspice not available")
        if ng_bin:
            ng_result, ng_err = run_ngspice(ng_bin, netlist, ng_node)

        if fc_t is None:
            # Fairchild failed — gray panel
            ax.set_facecolor("#f5f5f5")
            ax.text(0.5, 0.5, f"FAILED\n{fc_err[:60]}", transform=ax.transAxes,
                    ha="center", va="center", fontsize=7, color="#c62828",
                    bbox=dict(fc="white", ec="#c62828", pad=3))
            ax.set_title(title, fontsize=9, pad=3)
            print(f"fc FAILED: {fc_err[:40]}")
            summary_rows.append((title, "FAILED", "—", "—"))
            continue

        ax.plot(fc_t * xscale, fc_v, color=BLUE, lw=1.5, label="fairchild", zorder=3)

        rms_mv = None
        if ng_result is not None:
            ng_t, ng_v = ng_result
            ax.plot(ng_t * xscale, ng_v, "--", color=RED, lw=1.2,
                    label="ngspice", alpha=0.85, zorder=2)
            rms = rms_error(ng_t, ng_v, fc_t, fc_v)
            if rms is not None:
                rms_mv = rms * 1000
                ax.annotate(f"RMS {rms_mv:.2f} mV", xy=(0.97, 0.05),
                            xycoords="axes fraction", ha="right", va="bottom",
                            fontsize=7.5, color="#555")
        elif ng_bin:
            ax.annotate(f"ngspice error", xy=(0.97, 0.05),
                        xycoords="axes fraction", ha="right", va="bottom",
                        fontsize=7.5, color=RED)

        ax.set_title(title, fontsize=9, pad=3)
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.legend(fontsize=7, loc="upper right", framealpha=0.8)
        ax.grid(True, lw=0.4)
        ax.tick_params(labelsize=7)

        rms_str = f"{rms_mv:.2f} mV" if rms_mv is not None else ("—" if ng_bin else "no ngspice")
        summary_rows.append((title, "ok", rms_str, f"{len(fc_t)} pts"))
        print(f"ok  RMS={rms_str}")

    # Hide unused axes
    for ax in axes[n_panels:]:
        ax.set_visible(False)

    fig.suptitle("fairchild vs ngspice — waveform accuracy", fontsize=11, y=1.01)
    fig.tight_layout()

    # Print summary table
    print("\n  {:<25} {:>8} {:>12}".format("Circuit", "Status", "RMS error"))
    print("  " + "-" * 48)
    for row in summary_rows:
        print("  {:<25} {:>8} {:>12}".format(*row[:3]))

    out = PLOTS_DIR / "accuracy_analog.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Plot 2: scaling wall-clock
# ---------------------------------------------------------------------------

def backend_available(fc_bin: str, backend: str) -> bool:
    """Probe whether the fairchild binary supports a solver backend (KLU needs
    the `klu` cargo feature)."""
    probe = subprocess.run(
        [fc_bin, "-f", str(CIRCUITS / "ring_osc_3.sp"), "--solver", backend],
        capture_output=True, text=True, timeout=30,
    )
    return probe.returncode == 0


def plot_scaling(fc_bin: str, ng_bin: str | None) -> Path:
    print("\n  Timing ring oscillator scaling per solver backend "
          f"(median of {N_TIMING_RUNS} runs each) ...")

    # Which backends does this binary actually support?
    backends = []
    for b in SCALING_BACKENDS:
        if backend_available(fc_bin, b):
            backends.append(b)
        else:
            print(f"    backend '{b}' unavailable (rebuild with --features klu?) — skipping")

    # data[backend] = (nodes, times); ngspice tracked separately
    data: dict[str, tuple[list[int], list[float]]] = {b: ([], []) for b in backends}
    ng_nodes, ng_times = [], []
    dense_capped = False

    for fname, stages in SCALING_CIRCUITS:
        netlist = CIRCUITS / fname
        n_nodes = 2 * stages + 1  # NMOS + PMOS per stage + Vdd node

        # Build ngspice temp netlist (with control block)
        tmp_ng = None
        if ng_bin:
            src = netlist.read_text()
            lines = [l for l in src.splitlines() if l.strip().lower() != ".end"]
            body = "\n".join(lines)
            control = ".control\nrun\nprint v(n1)\n.endc\n.end\n"
            with tempfile.NamedTemporaryFile(mode="w", suffix=".sp", delete=False) as f:
                f.write(body + "\n" + control)
                tmp_ng = f.name

        print(f"    ring_osc_{stages} (n≈{n_nodes}): ", end="", flush=True)

        for b in backends:
            if b == "dense" and n_nodes > DENSE_MAX_NODES:
                dense_capped = True
                print(f"{b}=skip ", end="", flush=True)
                continue
            ms = time_median_ms([fc_bin, "-f", str(netlist), "--solver", b])
            data[b][0].append(n_nodes)
            data[b][1].append(ms)
            print(f"{b}={ms:.0f} ", end="", flush=True)

        if ng_bin and tmp_ng:
            ng_ms = time_median_ms([ng_bin, "-b", tmp_ng])
            ng_nodes.append(n_nodes)
            ng_times.append(ng_ms)
            print(f"ngspice={ng_ms:.0f}", end="")
        print()

        if tmp_ng:
            try:
                os.unlink(tmp_ng)
            except OSError:
                pass

    # ---- plot ----
    fig, ax = plt.subplots(figsize=(7.5, 5))

    for b in backends:
        nodes, times = data[b]
        if not nodes:
            continue
        color, marker, label = BACKEND_STYLE[b]
        ax.loglog(nodes, times, marker=marker, ls="-", color=color, lw=1.8, ms=6,
                  label=label, zorder=3)
    if ng_nodes:
        ax.loglog(ng_nodes, ng_times, marker="x", ls="--", color=RED, lw=1.6, ms=7,
                  label="ngspice (default solver)", alpha=0.9, zorder=2)

    if dense_capped:
        ax.text(0.02, 0.98,
                f"dense LU capped at n ≤ {DENSE_MAX_NODES} (O(N³); larger N omitted)",
                transform=ax.transAxes, ha="left", va="top", fontsize=7.5,
                color=BACKEND_STYLE["dense"][0])

    ax.set_xlabel("Approximate node count (2N+1 for N-stage ring oscillator)", fontsize=9)
    ax.set_ylabel(f"Wall-clock time (ms, median of {N_TIMING_RUNS} runs)", fontsize=9)
    ax.set_title("Transient simulation time vs circuit size, by solver backend\n"
                 "CMOS ring oscillator family (Level-1 MOSFET), end-to-end wall-clock",
                 fontsize=10)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(True, which="both", lw=0.4)
    ax.tick_params(labelsize=8)
    fig.tight_layout()

    out = PLOTS_DIR / "scaling_wall_time.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--release",    action="store_true", default=True,
                    help="use release build of fairchild (default)")
    ap.add_argument("--debug",      action="store_true",
                    help="use debug build instead")
    ap.add_argument("--no-ngspice", action="store_true",
                    help="skip ngspice; plot fairchild only")
    args = ap.parse_args()

    release = not args.debug
    fc_bin  = find_binary("fairchild", release)
    if not fc_bin:
        sys.exit(f"fairchild binary not found. Run: cargo build {'--release' if release else ''}")

    ng_bin = None if args.no_ngspice else find_binary("ngspice", release)
    if not args.no_ngspice and ng_bin is None:
        print("ngspice not found — plotting fairchild only", file=sys.stderr)

    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    print(f"\nfairchild: {fc_bin}")
    print(f"ngspice:   {ng_bin or '(not available)'}")
    print(f"plots →    {PLOTS_DIR.relative_to(REPO_ROOT)}\n")

    print("=== Accuracy panels ===")
    acc_out = plot_accuracy(fc_bin, ng_bin)
    print(f"\n  → {acc_out.relative_to(REPO_ROOT)}")

    print("\n=== Scaling plot ===")
    sc_out = plot_scaling(fc_bin, ng_bin)
    print(f"\n  → {sc_out.relative_to(REPO_ROOT)}")

    print("\nDone.")


if __name__ == "__main__":
    main()
