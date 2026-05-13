#!/usr/bin/env python3
"""
kicad_to_fairchild.py — KiCAD SPICE netlist → fairchild wrapper generator

Reads a KiCAD-exported SPICE netlist (.net or .sp), auto-detects optical nets
by scanning for _re/_im/_wl triplets, and emits a fairchild-ready wrapper that
includes:
  • .osdi directives for every model referenced in the netlist
  • .optical declaration covering all detected optical nets
  • .include of the KiCAD netlist
  • .op (or .tran if --tran is given)

Usage:
    python3 scripts/kicad_to_fairchild.py my_circuit.net \\
        --osdi-dir va-models/build \\
        --output run_my_circuit.sp

The wrapper can be run directly:
    DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib \\
        ./target/release/fairchild -f run_my_circuit.sp --probe "V(ph_a)"
"""

import argparse
import os
import re
import sys
from pathlib import Path

# ──────────────────────────────────────────────────────────────────────────────
# All fairchild photonic model names → their .osdi filename (sans extension).
# Models that ship as compound .spc subckts are listed under SUBCKT_MODELS.
# ──────────────────────────────────────────────────────────────────────────────

OSDI_MODELS = {
    "cw_laser",
    "waveguide",
    "directional_coupler",
    "pn_phase_shifter_l1",
    "pn_phase_shifter_l2",
    "thermo_phase_shifter_l1",
    "thermo_phase_shifter_l2",
    "photodetector",
    "photodetector_l2",
    "mrr_modulator_l1",
    "mrr_modulator_l2",
    "mrr_modulator_l3",
    "mrr_modulator_l1_adddrop",
    "mrr_modulator_l2_adddrop",
    "mrr_modulator_l3_adddrop",
    "mrr_heater_l1",
    "mrr_heater_l2",
    "mrr_heater_l1_adddrop",
    "mrr_heater_l2_adddrop",
    "mzi_modulator_pn_l1",
    "mzi_modulator_pn_l2",
    "mzi_modulator_thermo_l1",
    "mzi_modulator_thermo_l2",
}

SUBCKT_MODELS = {
    "mrr_allpass_pn_l1":    "va-models/photonic/subckts/mrr_allpass_pn_l1.spc",
    "mrr_allpass_thermo_l1":"va-models/photonic/subckts/mrr_allpass_thermo_l1.spc",
    "mrr_adddrop_pn_l1":    "va-models/photonic/subckts/mrr_adddrop_pn_l1.spc",
}


def parse_x_lines(lines):
    """
    Collect all X-element model names and net lists from logical SPICE lines
    (continuation lines already joined).
    Returns list of (model_name, [nets]).
    """
    results = []
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("*"):
            continue
        if stripped[0].lower() != "x":
            continue
        tokens = stripped.split()
        positional = [t for t in tokens[1:] if "=" not in t]
        if len(positional) < 2:
            continue
        model  = positional[-1].lower()
        nets   = positional[:-1]
        results.append((model, nets))
    return results


def join_continuation(raw_lines):
    """Join SPICE continuation lines (lines starting with '+') into logical lines."""
    logical = []
    for raw in raw_lines:
        if raw.startswith("+") and logical:
            logical[-1] = logical[-1].rstrip() + " " + raw[1:].strip()
        else:
            logical.append(raw)
    return logical


def detect_optical_nets(all_nets_by_model):
    """
    From a list of (model_name, [nets]) pairs, find all nets that belong to
    optical triplets: any net whose name ends in _re has a matching _im and _wl
    sibling (same base name).

    Also handles WDM bus nets: foo_re_0, foo_im_0, foo_wl_0, etc.

    Returns a list of optical net names in the order they are first seen.
    """
    all_nets = []
    for _, nets in all_nets_by_model:
        all_nets.extend(nets)

    # Build sets for fast lookup
    net_set = set(n.lower() for n in all_nets)

    seen = set()
    optical = []

    for net in all_nets:
        nl = net.lower()
        if nl in seen:
            continue
        # Check if this net ends in _re (with optional _N suffix for WDM)
        m = re.match(r"^(.+)_re(_\d+)?$", nl)
        if m:
            base  = m.group(1)
            idx   = m.group(2) or ""
            im_n  = f"{base}_im{idx}"
            wl_n  = f"{base}_wl{idx}"
            if im_n in net_set and wl_n in net_set:
                # Found a complete triplet
                for n in [nl, im_n, wl_n]:
                    if n not in seen:
                        optical.append(n)
                        seen.add(n)
    return optical


def find_osdi(model_name, osdi_dir):
    """Return path to .osdi file, or None if not found."""
    candidate = Path(osdi_dir) / f"{model_name}.osdi"
    if candidate.exists():
        return str(candidate)
    return None


def generate_wrapper(
    kicad_netlist_path,
    osdi_dir,
    fairchild_repo_root,
    analysis,
    tran_args,
    output_path,
    probe,
):
    kicad_path = Path(kicad_netlist_path)
    raw_lines  = kicad_path.read_text().splitlines()
    logical    = join_continuation(raw_lines)
    model_nets = parse_x_lines(logical)

    models_used = sorted({m for m, _ in model_nets})
    optical_nets = detect_optical_nets(model_nets)

    # Compute relative path from output dir to KiCAD netlist
    out_dir = Path(output_path).parent if output_path != "-" else Path.cwd()
    try:
        kicad_rel = os.path.relpath(kicad_path, out_dir)
    except ValueError:
        kicad_rel = str(kicad_path.resolve())

    lines = []
    lines.append(f"* Auto-generated by kicad_to_fairchild.py")
    lines.append(f"* Source: {kicad_path.name}")
    lines.append("")

    # .osdi directives
    lines.append("* ── OSDI models ──────────────────────────────────────────────")
    missing_osdi = []
    for model in models_used:
        if model in OSDI_MODELS:
            p = find_osdi(model, osdi_dir)
            if p:
                lines.append(f".osdi {p}")
            else:
                lines.append(f"* WARNING: {model}.osdi not found in {osdi_dir}")
                missing_osdi.append(model)
        elif model in SUBCKT_MODELS:
            spc_rel = os.path.join(fairchild_repo_root, SUBCKT_MODELS[model])
            lines.append(f".include {spc_rel}  * subckt model")
        # else: unknown model — leave for user to sort out

    lines.append("")

    # .optical declaration
    if optical_nets:
        lines.append("* ── Optical discipline ───────────────────────────────────────")
        # Break into groups of 12 nets per line for readability
        chunk_size = 12
        for i in range(0, len(optical_nets), chunk_size):
            chunk = " ".join(optical_nets[i:i + chunk_size])
            prefix = ".optical" if i == 0 else "+"
            lines.append(f"{prefix} {chunk}")
        lines.append("")
    else:
        lines.append("* No optical net triplets (_re/_im/_wl) detected.")
        lines.append("* Add a .optical directive manually if needed.")
        lines.append("")

    # Include the KiCAD netlist
    lines.append("* ── KiCAD netlist ────────────────────────────────────────────")
    lines.append(f".include \"{kicad_rel}\"")
    lines.append("")

    # Analysis
    if analysis == "op":
        lines.append(".op")
    elif analysis == "tran":
        lines.append(f".tran {tran_args}")

    lines.append(".end")
    lines.append("")

    output = "\n".join(lines)

    if output_path == "-":
        print(output)
    else:
        Path(output_path).write_text(output)
        print(f"Wrote {output_path}", file=sys.stderr)

    if missing_osdi:
        print(f"\nWARNING: {len(missing_osdi)} model(s) not found in {osdi_dir}:", file=sys.stderr)
        for m in missing_osdi:
            print(f"  {m}.osdi", file=sys.stderr)
        print("Run 'make' in va-models/ to compile all models.", file=sys.stderr)

    return optical_nets, models_used


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("netlist", help="KiCAD-exported SPICE netlist (.net or .sp)")
    ap.add_argument("--osdi-dir", default="va-models/build",
                    help="Directory containing compiled .osdi files (default: va-models/build)")
    ap.add_argument("--repo-root", default=".",
                    help="Path to fairchild repo root for resolving .spc subckt files (default: .)")
    ap.add_argument("--output", "-o", default="-",
                    help="Output wrapper netlist path (default: stdout)")
    ap.add_argument("--tran",
                    help="Add a .tran directive with these arguments, e.g. '1p 10n'")
    ap.add_argument("--probe", help="Probe expression(s) to display (informational only)")
    args = ap.parse_args()

    analysis  = "tran" if args.tran else "op"
    tran_args = args.tran or ""

    osdi_dir = os.path.join(args.repo_root, args.osdi_dir) \
               if not os.path.isabs(args.osdi_dir) else args.osdi_dir

    optical, models = generate_wrapper(
        kicad_netlist_path = args.netlist,
        osdi_dir           = osdi_dir,
        fairchild_repo_root= args.repo_root,
        analysis           = analysis,
        tran_args          = tran_args,
        output_path        = args.output,
        probe              = args.probe,
    )

    if args.output != "-":
        print(f"Detected {len(optical)} optical net(s) across {len(models)} model(s).", file=sys.stderr)
        if optical:
            print(f"Optical nets: {' '.join(optical[:9])}{'...' if len(optical) > 9 else ''}", file=sys.stderr)
        if args.probe:
            print(f"\nTo run:", file=sys.stderr)
            print(f"  DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib \\", file=sys.stderr)
            print(f"    ./target/release/fairchild -f {args.output} --probe \"{args.probe}\"", file=sys.stderr)


if __name__ == "__main__":
    main()
