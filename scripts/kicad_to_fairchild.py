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

# `fc_mux` / `fc_demux` are variable-arity bundle bridges. Pin 1 is the
# multi-channel bus side; pins 2..N+1 are the single-channel side. The
# channel count N is inferred from instance pin count: `1 + N` positional
# nets. The bus net needs a multi-channel `.optical_port NAME N` declaration
# and every device on the bus inherits that width through per-channel
# replication.
BUNDLE_BRIDGE_MODELS = {"fc_mux", "fc_demux"}

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
) -> Tuple[List[Tuple[str, int]], List[str]]:
    """
    Walk every native-device X-element, look up its port schema, and collect
    the set of nets that appear in a bundle position. Returns
    (bundle_nets_with_width, unknown_models) where bundle_nets_with_width is
    an ordered list of (net_name, channel_count) pairs.

    Channel-count inference:
      - `fc_mux` / `fc_demux` with M positional nets implies N = M − 1
        channels on the bus side (pin 1). All single-channel pins on the
        bridge stay at N = 1.
      - After MUX/DEMUX widths are known, propagate the width along bus
        wires by walking the X-element graph: any non-bridge device that
        connects bundle-side to a known multi-channel net inherits that
        width on its other bundle pins.
    """
    bundle_widths: Dict[str, int] = {}
    bundle_first_seen_order: List[str] = []
    unknown: List[str] = []

    def record(net: str, width: int, *, authoritative: bool = False):
        """Register a bundle net at the given channel width.  Non-authoritative
        records (regular devices defaulting to 1) lose to authoritative ones
        (MUX/DEMUX or width-propagation passes)."""
        nlc = net.lower()
        if nlc in GROUND_NETS:
            return
        existing = bundle_widths.get(nlc)
        if existing is None:
            bundle_widths[nlc] = width
            bundle_first_seen_order.append(net)
        elif existing != width:
            if authoritative:
                bundle_widths[nlc] = width
            else:
                # Quietly defer to the existing (likely authoritative) width.
                pass

    # ── Pass 1: fc_mux / fc_demux first — bus widths are authoritative.
    for refdes, model, nets in x_elements:
        if model not in BUNDLE_BRIDGE_MODELS:
            continue
        if len(nets) < 2:
            warn(f"{refdes} ({model}): need ≥ 2 nets (1 bus + ≥ 1 channels); got {len(nets)}")
            continue
        n_channels = len(nets) - 1
        bus_net = nets[0]
        record(bus_net, n_channels, authoritative=True)
        for ch_net in nets[1:]:
            record(ch_net, 1, authoritative=True)

    # ── Pass 2: regular devices, default to width 1.
    for refdes, model, nets in x_elements:
        if model in BUNDLE_BRIDGE_MODELS:
            continue
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
            if kind == "bundle":
                if net.lower() in GROUND_NETS:
                    warn(f"{refdes} ({model}): bundle pin tied to ground ('{net}') — likely wiring error")
                    continue
                record(net, 1)

    # ── Pass 3: propagate widths through non-bridge devices ─────────────
    # Any non-bridge device with multiple bundle pins should have matching
    # widths on every bundle pin (the parser enforces this). If one bundle
    # pin is already known to be N-channel and others are still at default
    # width 1, upgrade those to N.
    changed = True
    while changed:
        changed = False
        for refdes, model, nets in x_elements:
            if model in BUNDLE_BRIDGE_MODELS or model not in PORT_SCHEMA:
                continue
            schema = PORT_SCHEMA[model]
            if len(nets) != len(schema):
                continue
            bundle_pins = [n for kind, n in zip(schema, nets) if kind == "bundle"]
            if len(bundle_pins) < 2:
                continue
            widths_here = [bundle_widths.get(p.lower(), 1) for p in bundle_pins]
            max_w = max(widths_here)
            if max_w == 1:
                continue
            for pin in bundle_pins:
                pl = pin.lower()
                if bundle_widths.get(pl, 1) != max_w:
                    bundle_widths[pl] = max_w
                    changed = True

    ordered: List[Tuple[str, int]] = []
    for net in bundle_first_seen_order:
        w = bundle_widths.get(net.lower(), 1)
        ordered.append((net, w))
    # Also include any nets that appeared only through MUX/DEMUX processing
    # but weren't in first-seen order yet (rare; defensive).
    for net_lc, w in bundle_widths.items():
        if not any(n.lower() == net_lc for n, _ in ordered):
            ordered.append((net_lc, w))

    return ordered, sorted(set(unknown))


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

    n_wdm = sum(1 for _, w in bundle_nets if w > 1)
    lines: List[str] = []
    lines.append("* Auto-generated by kicad_to_fairchild.py — DO NOT EDIT BY HAND.")
    lines.append(f"* Source: {kicad_path.name}")
    lines.append(f"* Transpiled {n_transpiled} KiCad X-element(s); "
                 f"detected {len(bundle_nets)} optical bundle net(s)"
                 + (f", {n_wdm} of which are WDM (multi-channel)." if n_wdm else "."))
    lines.append("")

    if bundle_nets:
        lines.append("* ── Optical bundle declarations ──────────────────────────────")
        for net, width in bundle_nets:
            if width > 1:
                lines.append(f".optical_port {net} {width}")
            else:
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
