#!/usr/bin/env python3
"""
MRR modulator voltage sweep — parametric Python driver.

Sweeps the PN junction reverse bias on an MRR modulator from 0 V to -5 V
and plots the through-port optical power vs bias voltage.

Requires: fairchild Python package (maturin build/install), numpy, matplotlib.

Usage:
    cd examples/photonic
    python3 mrr_voltage_sweep.py
"""

import os
import sys
import subprocess
import numpy as np
import csv

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
FAIRCHILD_BIN = os.path.join(ROOT, "target", "release", "fairchild")
NETLIST = os.path.join(HERE, "mrr_modulator_dc.sp")
DYLD = f"DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib"

# ── Voltage sweep using the CLI --param flag ──────────────────────────────────

def run_cli_sweep(v_values):
    """Run one DC op per voltage bias via the CLI, collect V(ph_a)."""
    results = []
    for v in v_values:
        cmd = (
            f'{DYLD} {FAIRCHILD_BIN} -f {NETLIST} '
            f'--param "Vbias.dc={v}" --probe "V(ph_a)" -q'
        )
        try:
            out = subprocess.check_output(cmd, shell=True, text=True, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            print(f"  warning: simulation failed at V={v}", file=sys.stderr)
            results.append(None)
            continue

        # Parse CSV: header row, then data row
        lines = [l for l in out.strip().splitlines() if l]
        if len(lines) < 2:
            results.append(None)
            continue
        header = lines[0].split(",")
        data   = lines[1].split(",")
        row = dict(zip(header, data))
        v_pha = float(row.get("V(ph_a)", "nan"))
        # V(ph_a) is voltage across R_load=1kΩ; P_opt = V²/(R_load × responsivity²)
        i_ph  = v_pha / 1e3        # photocurrent (A) with R_load=1kΩ
        p_opt = i_ph / 1.0         # responsivity=1 A/W → optical power (W)
        results.append(p_opt)
    return results


def main():
    # Bias sweep: 0 V to -5 V in 0.25 V steps
    v_sweep = np.arange(0.0, -5.25, -0.25)

    print(f"Sweeping {len(v_sweep)} bias points: {v_sweep[0]:.2f} V → {v_sweep[-1]:.2f} V")
    print("Running fairchild CLI simulations...")

    p_through = run_cli_sweep(v_sweep)

    # Print results table
    print(f"\n{'Vbias (V)':>12}  {'P_through (mW)':>16}")
    print("-" * 32)
    for v, p in zip(v_sweep, p_through):
        p_mw = p * 1e3 if p is not None else float("nan")
        print(f"{v:>12.3f}  {p_mw:>16.4f}")

    # Find approximate resonance (minimum through-port power)
    p_arr = np.array([p if p is not None else np.nan for p in p_through])
    idx_min = int(np.nanargmin(p_arr))
    print(f"\nMinimum through-port power at Vbias = {v_sweep[idx_min]:.2f} V "
          f"(P = {p_arr[idx_min]*1e3:.4f} mW)")

    # Save to CSV
    out_csv = os.path.join(HERE, "mrr_voltage_sweep_result.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Vbias_V", "P_through_mW"])
        for v, p in zip(v_sweep, p_through):
            p_mw = p * 1e3 if p is not None else float("nan")
            w.writerow([f"{v:.3f}", f"{p_mw:.6f}"])
    print(f"\nResults saved to: {out_csv}")

    # Plot if matplotlib is available
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(v_sweep, p_arr * 1e3, "o-", ms=4, lw=1.5)
        ax.set_xlabel("Bias voltage (V)")
        ax.set_ylabel("Through-port power (mW)")
        ax.set_title("MRR Modulator — Voltage Sweep\n"
                     "L=100 µm, n_g=4.2, α=2 dB/cm, κ=0.1, Vpi_rt=10 V")
        ax.grid(True, alpha=0.4)
        ax.axvline(v_sweep[idx_min], color="r", linestyle="--", alpha=0.6,
                   label=f"resonance @ {v_sweep[idx_min]:.2f} V")
        ax.legend()
        plt.tight_layout()
        plot_path = os.path.join(HERE, "mrr_voltage_sweep.png")
        plt.savefig(plot_path, dpi=150)
        print(f"Plot saved to: {plot_path}")
        plt.show()
    except ImportError:
        print("(matplotlib not available; skipping plot)")


if __name__ == "__main__":
    main()
