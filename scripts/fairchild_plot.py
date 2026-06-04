#!/usr/bin/env python3
"""fairchild_plot.py — matplotlib results viewer for fairchild CSV output.

Reads a fairchild CSV (from `fairchild -f … --format csv`, file or stdin),
auto-detects the analysis from the first column, and renders sensible plots:

    time  → transient: optical power (|A|²), node voltages, branch currents
    freq  → AC: Bode (magnitude dB + phase)
    analysis (one row) → DC operating point: optical-power + voltage bars

Photonic signals exported as `V(<name>_re_<ch>)` / `_im_<ch>` / `_wl_<ch>`
triplets are recognised and collapsed into optical power
`P = re² + im²` (W → mW) per (name, channel) — usually what you want to see,
not the raw quadratures.

Usage:
    fairchild -f ckt.sp --format csv | fairchild_plot.py -            # stdin
    fairchild_plot.py results.csv --save out.png
    fairchild_plot.py results.csv --probe "V(out),P(drop)"           # filter
    fairchild_plot.py results.csv --raw                              # don't derive power

This is the visualization half of the KiCad "run + view" flow: the driver
(`kicad_fairchild.py --plot`) pipes its CSV here.
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict

import matplotlib

# Use a non-interactive backend automatically when saving / headless.
if "--save" in sys.argv or not sys.stdout.isatty():
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


# ── CSV parsing ────────────────────────────────────────────────────────────
def read_csv(stream):
    """Return (x_name, x_vals, {col_name: [floats]}) from a fairchild CSV."""
    header = stream.readline().strip()
    if not header:
        raise SystemExit("fairchild_plot: empty input")
    cols = header.split(",")
    data: list[list[float]] = [[] for _ in cols]
    for line in stream:
        line = line.strip()
        if not line:
            continue
        parts = line.split(",")
        if len(parts) != len(cols):
            continue
        for i, p in enumerate(parts):
            try:
                data[i].append(float(p))
            except ValueError:
                data[i].append(float("nan"))
    x_name = cols[0]
    series = {cols[i]: data[i] for i in range(1, len(cols))}
    return x_name, data[0], series


# ── Optical-power derivation ───────────────────────────────────────────────
# Matches V(<base>_re_<ch>) and the _im_ partner.
_OPT = re.compile(r"^V\((?P<base>.+)_(?P<qi>re|im)_(?P<ch>\d+)\)$")


def derive_optical_power(series: dict[str, list[float]]):
    """Collapse re/im quadrature pairs into power traces P(base.ch) in mW.

    Returns (powers, consumed) where `consumed` is the set of raw column names
    folded into a power trace (re/im/wl), so the caller can drop them.
    """
    re_im: dict[tuple[str, str], dict[str, list[float]]] = defaultdict(dict)
    consumed: set[str] = set()
    for name, vals in series.items():
        m = _OPT.match(name)
        if m:
            re_im[(m["base"], m["ch"])][m["qi"]] = vals
            consumed.add(name)
        elif name.endswith(")") and "_wl_" in name:
            consumed.add(name)  # hide the λ pass-through wire from plots
    powers: dict[str, list[float]] = {}
    for (base, ch), qi in re_im.items():
        if "re" in qi and "im" in qi:
            label = f"P({base}.{ch})" if ch != "0" else f"P({base})"
            powers[label] = [
                1e3 * (r * r + i * i) for r, i in zip(qi["re"], qi["im"])
            ]  # W → mW
        else:
            consumed.discard(name)  # incomplete pair: keep raw
    return powers, consumed


# ── Plotters ───────────────────────────────────────────────────────────────
def _split_electrical(series):
    """Partition non-optical signals into voltages and currents."""
    volts = {k: v for k, v in series.items() if k.upper().startswith("V(")}
    currs = {k: v for k, v in series.items() if k.upper().startswith("I(")}
    other = {
        k: v for k, v in series.items() if k not in volts and k not in currs
    }
    return volts, currs, other


def plot_transient(x, series, raw, title):
    powers, consumed = ({}, set()) if raw else derive_optical_power(series)
    rest = {k: v for k, v in series.items() if k not in consumed}
    volts, currs, other = _split_electrical(rest)
    panels = [p for p in (("Optical power (mW)", powers),
                          ("Node voltage (V)", {**volts, **other}),
                          ("Branch current (A)", currs)) if p[1]]
    if not panels:
        panels = [("Signals", series)]
    fig, axes = plt.subplots(len(panels), 1, sharex=True, figsize=(9, 2.6 * len(panels)))
    if len(panels) == 1:
        axes = [axes]
    x_us = [t * 1e6 for t in x]  # µs is a friendly transient unit
    for ax, (ylabel, sigs) in zip(axes, panels):
        for name, vals in sigs.items():
            ax.plot(x_us, vals, label=name, lw=1.3)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        if len(sigs) <= 12:
            ax.legend(fontsize=8, ncol=2, loc="best")
    axes[-1].set_xlabel("time (µs)")
    fig.suptitle(title or "fairchild — transient")
    fig.tight_layout()
    return fig


def plot_ac(x, series, title):
    import math

    mags = {k: v for k, v in series.items() if k.lower().startswith("mag")}
    phases = {k: v for k, v in series.items() if k.lower().startswith("phase")}
    if not mags:  # fall back: treat all as magnitude
        mags = series
    fig, (a_mag, a_ph) = plt.subplots(2, 1, sharex=True, figsize=(9, 6))
    for name, vals in mags.items():
        db = [20.0 * math.log10(v) if v > 0 else float("nan") for v in vals]
        a_mag.semilogx(x, db, label=name, lw=1.3)
    a_mag.set_ylabel("magnitude (dB)")
    a_mag.grid(True, which="both", alpha=0.3)
    a_mag.legend(fontsize=8)
    for name, vals in phases.items():
        a_ph.semilogx(x, vals, label=name, lw=1.3)
    a_ph.set_ylabel("phase (deg)")
    a_ph.set_xlabel("frequency (Hz)")
    a_ph.grid(True, which="both", alpha=0.3)
    fig.suptitle(title or "fairchild — AC")
    fig.tight_layout()
    return fig


def plot_op(series, raw, title):
    powers, consumed = ({}, set()) if raw else derive_optical_power(series)
    rest = {k: v[0] for k, v in series.items() if k not in consumed and v}
    fig, axes = plt.subplots(1, 2 if powers else 1, figsize=(10, 4.5), squeeze=False)
    axes = axes[0]
    if powers:
        names = list(powers)
        axes[0].barh(names, [p[0] for p in powers.values()], color="tab:orange")
        axes[0].set_xlabel("optical power (mW)")
        axes[0].set_title("Optical")
    ax_e = axes[-1]
    names = list(rest)
    ax_e.barh(names, [rest[n] for n in names], color="tab:blue")
    ax_e.set_xlabel("value (V / A)")
    ax_e.set_title("Electrical")
    fig.suptitle(title or "fairchild — operating point")
    fig.tight_layout()
    return fig


def main():
    ap = argparse.ArgumentParser(description="Plot fairchild CSV results.")
    ap.add_argument("csv", help="CSV file, or '-' for stdin")
    ap.add_argument("--save", metavar="PNG", help="save to PNG instead of showing")
    ap.add_argument("--probe", help="comma-separated signals to keep")
    ap.add_argument("--raw", action="store_true", help="don't derive optical power from re/im")
    ap.add_argument("--title", help="figure title")
    args = ap.parse_args()

    stream = sys.stdin if args.csv == "-" else open(args.csv)
    x_name, x, series = read_csv(stream)
    if args.probe:
        keep = {s.strip() for s in args.probe.split(",")}
        series = {k: v for k, v in series.items() if k in keep}

    xl = x_name.lower()
    if xl.startswith("time"):
        fig = plot_transient(x, series, args.raw, args.title)
    elif xl.startswith("freq"):
        fig = plot_ac(x, series, args.title)
    else:  # "analysis" → operating point (single row)
        fig = plot_op(series, args.raw, args.title)

    if args.save:
        fig.savefig(args.save, dpi=130)
        print(f"[fairchild_plot] wrote {args.save}", file=sys.stderr)
    else:
        plt.show()


if __name__ == "__main__":
    main()
