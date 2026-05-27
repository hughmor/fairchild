#!/usr/bin/env python3
"""
Generate a CMOS ring oscillator SPICE netlist for arbitrary N stages.

Usage:
    python benchmarks/gen_ring_osc.py <N> [--output FILE]

Writes to FILE if given, else to benchmarks/circuits/ring_osc_<N>.sp.
N must be odd and >= 3.

Parameters match the existing ring_osc_*.sp benchmark circuits:
    NMOS: vto=0.5  kp=200u  lambda=0.05  W=10u  L=1u
    PMOS: vto=-0.5 kp=80u   lambda=0.05  W=20u  L=1u
    Load: C=100f per node
    VDD = 1.8V
    tran: step=50p, stop chosen so that ~3 cycles are simulated
          (t_pd ≈ 0.5ns per stage → period ≈ N ns, stop ≈ 3N ns)
"""

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
CIRCUITS = Path(__file__).parent / "circuits"


def t_stop_ns(n: int) -> float:
    """Simulate ~1.5 full oscillation periods. t_pd ≈ 0.5 ns/stage empirically."""
    period_ns = n * 1.0  # conservative: 1 ns/stage
    return max(1.5 * period_ns, 10.0)


def t_step_ps(n: int) -> float:
    """Scale timestep with circuit size so each run takes ~3000 timesteps.

    For N≤51: 50ps (matches existing benchmark circuits, captures edge transitions).
    For N>51: step = period_ns/3000 * 1000ps ≈ n/3 ps, rounded to nice values.
    """
    if n <= 51:
        return 50.0
    # period ≈ n ns; step = n_ns / 3000 * 1000 ps/ns = n/3 ps
    step_ps = n / 3.0
    # Round to nearest 50ps
    return max(50.0, round(step_ps / 50) * 50)


def gen_netlist(n: int) -> str:
    if n < 3 or n % 2 == 0:
        raise ValueError(f"N must be odd and >= 3, got {n}")

    stop_ns = t_stop_ns(n)
    stop_str = f"{stop_ns:.0f}n"

    lines = [
        f"* {n}-stage CMOS ring oscillator",
        f"* f ≈ 1/(2·N·t_pd); t_pd set by C_load/I_drive.",
        ".model nm NMOS (vto=0.5 kp=200u lambda=0.05)",
        ".model pm PMOS (vto=-0.5 kp=80u lambda=0.05)",
        "Vdd  vdd 0   DC 1.8",
        "",
    ]

    for i in range(1, n + 1):
        prev = n if i == 1 else i - 1
        lines += [
            f"Mn{i}  n{i} n{prev} 0   0   nm w=10u l=1u",
            f"Mp{i}  n{i} n{prev} vdd vdd pm w=20u l=1u",
            f"C{i}   n{i} 0  100f",
            "",
        ]

    # Alternating initial conditions: odd nodes HIGH (1.6V), even nodes LOW (0.1V)
    ic_pairs = [f"V(n{i})={'1.6' if i % 2 == 1 else '0.1'}" for i in range(1, n + 1)]
    step_ps = t_step_ps(n)
    step_str = f"{step_ps:.0f}p"

    lines.append(".ic " + " ".join(ic_pairs))
    lines.append(".options method=gear")
    lines.append(f".tran {step_str} {stop_str} UIC")
    lines.append(".end")

    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("n", type=int, help="Number of stages (odd, >= 3)")
    ap.add_argument("--output", "-o", default=None, help="Output file path")
    args = ap.parse_args()

    n = args.n
    netlist = gen_netlist(n)

    if args.output:
        out = Path(args.output)
    else:
        out = CIRCUITS / f"ring_osc_{n}.sp"

    out.write_text(netlist)
    print(f"Wrote {n}-stage ring oscillator to {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
