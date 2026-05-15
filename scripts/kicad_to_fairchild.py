#!/usr/bin/env python3
"""
kicad_to_fairchild.py — KiCad SPICE export → fairchild wrapper netlist.

Reads a KiCad-exported SPICE netlist that uses native fairchild photonic
devices (`fc_cw_laser`, `fc_waveguide`, `fc_dcoupler`, `fc_splitter`,
`fc_pn_ps`, `fc_thermal_ps`, `fc_photodetector`) and emits a self-contained
wrapper netlist that:

  • Transpiles KiCad's X-element lines (`REF net1 ... type=X model=NAME …`)
    into fairchild form (`XREF net1 ... NAME …`).
  • Declares every detected optical bundle net via `.optical_port NAME`.
  • Comments out KiCad-specific directives fairchild doesn't recognise
    (`.save`, `.probe`, `.title`).
  • Inlines the transpiled netlist body and appends the analysis directive.

Why transpile instead of `.include`-ing the raw export? fairchild's SPICE
parser dispatches by the first letter of an element line (`R*` → resistor,
`V*` → voltage source, etc.). KiCad reference designators like `CWL1` or
`WG1` would be misclassified as capacitors / inductors without the leading
`X` prefix. Rather than teach the parser KiCad's quirks, the post-processor
hides them.

Convention: in the KiCad symbol library, each native device's "optical port"
is a SINGLE KiCad pin; the connected net name is a bundle name that the
fairchild parser expands into (re, im, wl) wires. Electrical pins (anode,
cathode, heat_p, heat_n) are ordinary scalar nets.

Usage
-----
    python3 scripts/kicad_to_fairchild.py my_circuit.cir \\
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

GROUND_NETS = {"0", "gnd"}

# KiCad directives we silently comment out — fairchild's parser doesn't
# recognise them and would otherwise warn or error.
KICAD_PROLOG_DIRECTIVES = {".save", ".probe", ".title"}

# KiCad emits voltage / current sources in the canonical SPICE combined form
#     V<name> pos neg DC <op_value> SIN(...) AC <ac_mag>
# so the same source can drive a DC operating point, a transient, and an AC
# sweep. fairchild's parser doesn't yet handle this combined form — it picks
# up the `DC <op_value>` and silently drops the transient waveform. Strip the
# `DC <value>` prefix when a non-DC waveform follows; the OP will then read
# the waveform's first sample (zero crossing or PWL's t=0 point) which is the
# right behaviour. Long-term fix lives in fairchild's parser.
DC_BEFORE_WAVEFORM_RE = re.compile(
    r"\s+DC\s+\S+\s+(?=(?:SIN|PULSE|EXP|PWL|SFFM|AM)\s*\()",
    re.IGNORECASE,
)


# ── Logical-line collection (handles SPICE continuation lines) ──────────────

def strip_dc_prefix_before_waveform(line: str) -> str:
    """Drop `DC <value>` from `V*` / `I*` lines that also carry a transient
    waveform. See DC_BEFORE_WAVEFORM_RE for context."""
    s = line.lstrip()
    if not s or s[0].lower() not in ("v", "i"):
        return line
    return DC_BEFORE_WAVEFORM_RE.sub(" ", line)


def join_continuations(raw_lines: List[str]) -> List[str]:
    """Merge SPICE continuation lines (`+`) into logical lines."""
    logical: List[str] = []
    for raw in raw_lines:
        if raw.startswith("+") and logical:
            logical[-1] = logical[-1].rstrip() + " " + raw[1:].strip()
        else:
            logical.append(raw)
    return logical


# ── KiCad X-element detection + transpilation ───────────────────────────────

def try_transpile_x_element(
    s: str,
) -> Optional[Tuple[str, str, List[str]]]:
    """
    Try to interpret `s` as a KiCad-style X-element. KiCad emits a line like:

        REFDES net1 net2 ... type=X model=NAME k1=v1 k2=v2 ...

    On success, returns (transpiled_line, model_lowercase, [nets]).
    On failure (not an X-element), returns None.

    The transpiled line has the form:

        XREFDES net1 net2 ... NAME k1=v1 k2=v2 ...

    so fairchild's parser dispatches correctly on the leading `X`.
    """
    tokens = s.split()
    if len(tokens) < 2:
        return None

    # Scan kwargs for `type=X` and `model=NAME`.
    has_type_x = False
    model: Optional[str] = None
    for t in tokens[1:]:
        if "=" not in t:
            continue
        k, v = t.split("=", 1)
        kl = k.lower()
        if kl == "type" and v.strip().upper() == "X":
            has_type_x = True
        elif kl == "model":
            model = v.strip()

    if not has_type_x or model is None:
        return None

    # Split positional nets from "other" kwargs (drop type= and model=).
    refdes = tokens[0]
    nets: List[str] = []
    other_kwargs: List[str] = []
    for t in tokens[1:]:
        if "=" not in t:
            nets.append(t)
        else:
            k, _ = t.split("=", 1)
            if k.lower() in ("type", "model"):
                continue
            other_kwargs.append(t)

    if not refdes.lower().startswith("x"):
        refdes = "X" + refdes

    new_line = " ".join([refdes] + nets + [model] + other_kwargs)
    return new_line, model.lower(), nets


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

    for refdes, model, nets in x_elements:
        if model not in PORT_SCHEMA:
            unknown.append(model)
            continue
        schema = PORT_SCHEMA[model]
        if len(nets) != len(schema):
            warn(
                f"{refdes} ({model}): expected {len(schema)} positional nets, "
                f"got {len(nets)}; skipping bundle detection for this instance"
            )
            continue
        for kind, net in zip(schema, nets):
            if kind != "bundle":
                continue
            nlc = net.lower()
            if nlc in GROUND_NETS:
                warn(f"{refdes} ({model}): bundle pin tied to ground ('{net}') — likely wiring error")
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
    analysis_args: str,
    options_lines: List[str],
    probe: Optional[str],
    verbose: bool,
) -> Tuple[List[str], List[str]]:
    raw = kicad_path.read_text().splitlines()
    logical = join_continuations(raw)

    warnings: List[str] = []
    warn = warnings.append

    transpiled_body: List[str] = []
    x_elements: List[Tuple[str, str, List[str]]] = []
    n_transpiled = 0

    for line in logical:
        s = line.strip()

        # Pass-through: blank lines and comments.
        if not s or s.startswith("*") or s.startswith(";"):
            transpiled_body.append(line)
            continue

        lc = s.lower()

        # Hard EOF — strip; we'll add our own .end at the bottom.
        if lc == ".end":
            continue

        # KiCad prolog directives fairchild doesn't recognise.
        first_token = lc.split()[0]
        if first_token in KICAD_PROLOG_DIRECTIVES:
            transpiled_body.append("* [stripped by kicad_to_fairchild] " + line)
            continue

        # Try KiCad-style X-element transpile.
        result = try_transpile_x_element(s)
        if result is not None:
            new_line, model, nets = result
            transpiled_body.append(new_line)
            refdes = new_line.split(maxsplit=1)[0]
            x_elements.append((refdes, model, nets))
            n_transpiled += 1
            continue

        # Anything else: pass through, but rewrite V/I sources to strip the
        # combined `DC X SIN(...)` form that fairchild's parser can't handle.
        transpiled_body.append(strip_dc_prefix_before_waveform(line))

    bundle_nets, unknown_models = collect_bundle_nets(x_elements, warn)

    lines: List[str] = []
    lines.append("* Auto-generated by kicad_to_fairchild.py — DO NOT EDIT BY HAND.")
    lines.append(f"* Source: {kicad_path.name}")
    lines.append(f"* Transpiled {n_transpiled} KiCad X-element(s); "
                 f"detected {len(bundle_nets)} optical bundle net(s).")
    lines.append("")

    if bundle_nets:
        lines.append("* ── Optical bundle declarations ──────────────────────────────")
        for net in bundle_nets:
            lines.append(f".optical_port {net}")
        lines.append("")

    if options_lines:
        lines.extend(options_lines)
        lines.append("")

    lines.append("* ── Analysis ────────────────────────────────────────────────")
    if analysis == "op":
        lines.append(".op")
    elif analysis == "tran":
        lines.append(f".tran {analysis_args}")
    elif analysis == "ac":
        lines.append(f".ac {analysis_args}")
    lines.append("")

    lines.append("* ── Transpiled KiCad netlist body ────────────────────────────")
    lines.extend(transpiled_body)
    if not lines[-1].strip() == "":
        lines.append("")
    lines.append(".end")
    lines.append("")

    text = "\n".join(lines)
    if str(output_path) == "-":
        print(text)
    else:
        output_path.write_text(text)
        if verbose:
            print(f"wrote {output_path}", file=sys.stderr)

    if unknown_models and verbose:
        print(
            f"info: {len(unknown_models)} non-native X-element model(s) detected; "
            f"passed through untouched: {', '.join(unknown_models)}",
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
    ap.add_argument("netlist", help="KiCad-exported SPICE netlist (.cir / .net / .sp)")
    ap.add_argument("-o", "--output", default="-",
                    help="Wrapper netlist path (default: stdout)")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--op", action="store_true",
                   help="Emit a .op analysis (default if --tran/--ac omitted)")
    g.add_argument("--tran", metavar="ARGS",
                   help="Emit .tran with these arguments, e.g. '5n 2u'")
    g.add_argument("--ac", metavar="ARGS",
                   help="Emit .ac with these arguments, e.g. 'dec 20 1 1G'")
    ap.add_argument("--opt", action="append", default=[], metavar="KEY=VALUE",
                    help="SimOptions override emitted as .options KEY=VALUE (repeatable)")
    ap.add_argument("--method", metavar="be|tr|gear",
                    help="Convenience: equivalent to --opt method=...")
    ap.add_argument("--probe",
                    help="Probe expression(s) to surface in run-helper hint (informational)")
    ap.add_argument("-v", "--verbose", action="store_true")

    args = ap.parse_args()

    if args.tran:
        analysis, analysis_args = "tran", args.tran
    elif args.ac:
        analysis, analysis_args = "ac", args.ac
    else:
        analysis, analysis_args = "op", ""

    options_lines: List[str] = []
    flat_opts: List[str] = []
    if args.method:
        flat_opts.append(f"method={args.method}")
    flat_opts.extend(args.opt)
    if flat_opts:
        options_lines.append(".options " + " ".join(flat_opts))

    kicad_path = Path(args.netlist).resolve()
    output_path = Path(args.output) if args.output != "-" else Path("-")

    emit_wrapper(
        kicad_path=kicad_path,
        output_path=output_path,
        analysis=analysis,
        analysis_args=analysis_args,
        options_lines=options_lines,
        probe=args.probe,
        verbose=args.verbose,
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
