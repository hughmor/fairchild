#!/usr/bin/env python3
"""
kicad_to_fairchild.py — KiCad SPICE export → fairchild wrapper netlist.

Reads a KiCad-exported SPICE netlist that uses native fairchild photonic
devices (`fc_cw_laser`, `fc_waveguide`, `fc_dcoupler`, `fc_splitter`,
`fc_pn_ps`, `fc_thermal_ps`, `fc_photodetector`) and emits a wrapper netlist
that:

  • Declares every detected optical bundle net via `.optical_port NAME`.
  • Includes the KiCad-exported netlist.
  • Appends the requested analysis (`.op` or `.tran ...`).

The convention assumed throughout: in the KiCad symbol library, each native
device's "optical port" is a SINGLE KiCad pin, and the connected net name is
a bundle name that the fairchild parser expands into (re, im, wl) wires.
Electrical pins (anode, cathode, heat_p, heat_n) are ordinary scalar nets.

Usage
-----
    python3 scripts/kicad_to_fairchild.py my_circuit.net \\
        --tran "5n 2u" \\
        --output run_my_circuit.sp

    fairchild -f run_my_circuit.sp

Legacy OSDI-based photonic models are not handled. Pin a pre-Phase-B commit
of this script if you need that path; new work should be on native devices.
"""
import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ── Port schema for every native photonic device ────────────────────────────
#
# Each entry maps model name → list of positional ports. Each port is one of:
#   "bundle" — connected net is an optical bundle name. The fairchild parser
#              expands it into three underlying wires (re, im, wl) per channel
#              when a matching `.optical_port` directive is declared.
#   "scalar" — connected net is an ordinary electrical / ground net. No
#              `.optical_port` declaration needed.
#
# The list length equals the number of positional nets in a bundle-style
# X-element instance. Keep in sync with `register_native_photonics` in
# crates/fairchild-core/src/device_registry.rs.

PORT_SCHEMA: Dict[str, List[str]] = {
    "fc_cw_laser":      ["bundle"],                                   # out
    "fc_waveguide":     ["bundle", "bundle"],                         # in, out
    "fc_dcoupler":      ["bundle", "bundle", "bundle", "bundle"],     # a1, a2, b1, b2
    "fc_splitter":      ["bundle", "bundle", "bundle"],               # in, out_a, out_b
    "fc_photodetector": ["bundle", "scalar", "scalar"],               # in, anode, cathode
    "fc_pn_ps":         ["bundle", "bundle", "scalar", "scalar"],     # in, out, anode, cathode
    "fc_thermal_ps":    ["bundle", "bundle", "scalar", "scalar"],     # in, out, heat_p, heat_n
}

# Nets that must never be declared as bundles even if a bundle-position pin
# is wired to them (degenerate user error or wiring an unused port to ground).
GROUND_NETS = {"0", "gnd"}


# ── KiCad SPICE export parsing ──────────────────────────────────────────────

def join_continuations(raw_lines: List[str]) -> List[str]:
    """Merge SPICE continuation lines (`+`) into logical lines."""
    logical: List[str] = []
    for raw in raw_lines:
        if raw.startswith("+") and logical:
            logical[-1] = logical[-1].rstrip() + " " + raw[1:].strip()
        else:
            logical.append(raw)
    return logical


def parse_x_elements(lines: List[str]) -> List[Tuple[str, str, List[str]]]:
    """
    Extract every X-element instance as (instance_name, model_name, [nets]).
    Skips comments and blank lines. `nets` excludes the model name and any
    `key=value` parameter tokens.
    """
    out: List[Tuple[str, str, List[str]]] = []
    for line in lines:
        s = line.strip()
        if not s or s.startswith("*") or s.startswith(";"):
            continue
        if s[0].lower() != "x":
            continue
        tokens = s.split()
        if len(tokens) < 3:
            continue
        positional = [t for t in tokens[1:] if "=" not in t]
        if len(positional) < 2:
            continue
        instance = tokens[0]
        model = positional[-1].lower()
        nets = positional[:-1]
        out.append((instance, model, nets))
    return out


# ── Bundle-net detection ────────────────────────────────────────────────────

def collect_bundle_nets(
    x_elements: List[Tuple[str, str, List[str]]],
    warn,
) -> Tuple[List[str], List[str]]:
    """
    Walk every native-device X-element, look up its port schema, and collect
    the set of nets that appear in a bundle position. Returns (bundle_nets,
    unknown_models) where bundle_nets preserves first-seen insertion order.
    """
    seen = set()
    bundle_nets: List[str] = []
    unknown: List[str] = []

    for inst, model, nets in x_elements:
        if model not in PORT_SCHEMA:
            unknown.append(model)
            continue
        schema = PORT_SCHEMA[model]
        if len(nets) != len(schema):
            warn(
                f"{inst} ({model}): expected {len(schema)} positional nets, "
                f"got {len(nets)}; skipping bundle detection"
            )
            continue
        for kind, net in zip(schema, nets):
            if kind != "bundle":
                continue
            nlc = net.lower()
            if nlc in GROUND_NETS:
                warn(f"{inst} ({model}): bundle pin tied to ground ('{net}') — likely wiring error")
                continue
            if nlc in seen:
                continue
            seen.add(nlc)
            bundle_nets.append(net)

    return bundle_nets, sorted(set(unknown))


# ── Wrapper emission ────────────────────────────────────────────────────────

def emit_wrapper(
    kicad_path: Path,
    output_path: Path,
    analysis: str,
    tran_args: str,
    options_lines: List[str],
    probe: Optional[str],
    verbose: bool,
) -> Tuple[List[str], List[str]]:
    raw = kicad_path.read_text().splitlines()
    logical = join_continuations(raw)
    x_elements = parse_x_elements(logical)

    warnings: List[str] = []
    warn = warnings.append

    bundle_nets, unknown_models = collect_bundle_nets(x_elements, warn)

    out_dir = output_path.parent if str(output_path) != "-" else Path.cwd()
    try:
        kicad_rel = str(kicad_path.relative_to(out_dir.resolve()))
    except ValueError:
        kicad_rel = str(kicad_path.resolve())

    lines: List[str] = []
    lines.append(f"* Auto-generated by kicad_to_fairchild.py")
    lines.append(f"* Source: {kicad_path.name}")
    lines.append(f"* Detected {len(bundle_nets)} optical bundle net(s) across "
                 f"{len({m for _, m, _ in x_elements if m in PORT_SCHEMA})} "
                 f"native photonic instance type(s).")
    lines.append("")

    if bundle_nets:
        lines.append("* ── Optical bundle declarations ──────────────────────────────")
        for net in bundle_nets:
            lines.append(f".optical_port {net}")
        lines.append("")
    else:
        lines.append("* No native photonic devices detected — circuit is pure electrical.")
        lines.append("")

    for opt in options_lines:
        lines.append(opt)
    if options_lines:
        lines.append("")

    # Analysis directives go BEFORE the .include, because the KiCad-exported
    # file conventionally ends with `.end`, which the fairchild parser treats
    # as a hard EOF — anything after the include is dropped on the floor.
    lines.append("* ── Analysis ────────────────────────────────────────────────")
    if analysis == "op":
        lines.append(".op")
    elif analysis == "tran":
        lines.append(f".tran {tran_args}")
    elif analysis == "ac":
        # Caller passes the .ac arguments through tran_args slot for simplicity.
        lines.append(f".ac {tran_args}")
    lines.append("")

    lines.append("* ── KiCad netlist ────────────────────────────────────────────")
    lines.append(f'.include "{kicad_rel}"')
    lines.append("")

    text = "\n".join(lines)
    if str(output_path) == "-":
        print(text)
    else:
        output_path.write_text(text)
        if verbose:
            print(f"wrote {output_path}", file=sys.stderr)

    # Side-channel warnings & info.
    if unknown_models and verbose:
        print(
            f"info: {len(unknown_models)} non-native X-element model(s) in netlist; "
            f"passing through untouched: {', '.join(unknown_models)}",
            file=sys.stderr,
        )
    for w in warnings:
        print(f"warning: {w}", file=sys.stderr)
    if verbose and probe:
        print(f"\nTo run:\n  fairchild -f {output_path} --probe \"{probe}\"", file=sys.stderr)

    return bundle_nets, unknown_models


# ── CLI ─────────────────────────────────────────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("netlist", help="KiCad-exported SPICE netlist (.net or .sp)")
    ap.add_argument("-o", "--output", default="-",
                    help="Wrapper netlist path (default: stdout)")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--op", action="store_true", help="Emit a .op analysis (default if --tran/--ac omitted)")
    g.add_argument("--tran", metavar="ARGS",
                   help="Emit a .tran with these arguments, e.g. '5n 2u'")
    g.add_argument("--ac", metavar="ARGS",
                   help="Emit a .ac with these arguments, e.g. 'dec 20 1 1G'")
    ap.add_argument("--opt", action="append", default=[], metavar="KEY=VALUE",
                    help="Pass a SimOptions override into the wrapper as .options KEY=VALUE")
    ap.add_argument("--method", metavar="be|tr|gear",
                    help="Convenience: equivalent to --opt method=...")
    ap.add_argument("--probe",
                    help="Probe expression(s) to surface in the run-helper hint (informational)")
    ap.add_argument("-v", "--verbose", action="store_true")

    args = ap.parse_args()

    if args.tran:
        analysis, tran_args = "tran", args.tran
    elif args.ac:
        analysis, tran_args = "ac", args.ac
    else:
        analysis, tran_args = "op", ""

    options_lines: List[str] = []
    flat_opts = [f"method={args.method}"] if args.method else []
    flat_opts.extend(args.opt)
    if flat_opts:
        options_lines.append(".options " + " ".join(flat_opts))

    kicad_path = Path(args.netlist).resolve()
    output_path = Path(args.output) if args.output != "-" else Path("-")

    bundle_nets, unknown = emit_wrapper(
        kicad_path=kicad_path,
        output_path=output_path,
        analysis=analysis,
        tran_args=tran_args,
        options_lines=options_lines,
        probe=args.probe,
        verbose=args.verbose,
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
