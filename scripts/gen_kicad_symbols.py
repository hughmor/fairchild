#!/usr/bin/env python3
"""Generate KiCad symbols for fairchild photonic models.

The `.kicad_sym` s-expression is regular enough that a symbol is worth
generating rather than drawing: pins on a grid, one body rectangle, a
`Sim.Device` / `Sim.Params` pair, and some line art. What is NOT regular is the
pin *order* — it has to match the positional net order on the `X…` line, and
KiCad orders by pin number — so that is what this file actually encodes.

    python3 scripts/gen_kicad_symbols.py            # write into the library
    python3 scripts/gen_kicad_symbols.py --check    # exit 1 if it would change
    python3 scripts/gen_kicad_symbols.py --print fc_facet

**Only symbols listed in `SYMBOLS` are touched.** A symbol already in the
library is replaced wholesale, so do not add a spec for one whose art someone
drew by hand unless you mean to lose it. Everything else is left alone, which
is why this can be re-run safely on a library that is mostly hand-drawn.

Keep `SYMBOLS` in sync with `PORT_SCHEMA` in `scripts/kicad_to_fairchild.py`
and `register_native_photonics` in `crates/fairchild-core/src/device_registry.rs`.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Tuple

LIB = Path(__file__).resolve().parents[1] / "examples/kicad_photonics/fairchild_photonics.kicad_sym"

# Pin sides. `L` connects on the left and points right; `R` is the mirror.
LEFT, RIGHT = "L", "R"

# Pin electrical type. Optical bundles are `bidirectional` throughout the
# existing library — a bundle carries light both ways under
# `enable_bidirectional`, and KiCad's ERC has no optical discipline anyway.
OPTICAL, ELECTRICAL = "bidirectional", "passive"


def pin(name: str, side: str, y: float, kind: str = OPTICAL) -> dict:
    return {"name": name, "side": side, "y": y, "kind": kind}


SYMBOLS = [
    {
        "name": "fc_driven_laser",
        "ref": "LD",
        "params": "type=X model=fc_driven_laser slope_w_v=4e-3 v_th=0.9 wavelength_nm=1550.0",
        # The body is wider than the art needs: KiCad draws pin NAMES inside
        # the box, so leaving room for them is what stops "out" landing on the
        # glyph. Same reason the emission arrows sit off the y = 0 centreline.
        "box": (-5.08, 2.54, 5.08, -2.54),
        # Laser-diode glyph: a filled triangle into a bar, then emission.
        "art": [
            ("poly", [(-2.54, -1.27), (-2.54, 1.27), (-0.635, 0.0), (-2.54, -1.27)], "outline"),
            ("poly", [(-0.635, 1.27), (-0.635, -1.27)], "none"),
            ("poly", [(0.0, 1.905), (1.524, 1.905), (1.143, 2.159)], "none"),
            ("poly", [(1.524, 1.905), (1.143, 1.651)], "none"),
            ("poly", [(0.0, -1.905), (1.524, -1.905), (1.143, -1.651)], "none"),
            ("poly", [(1.524, -1.905), (1.143, -2.159)], "none"),
        ],
        # Positional order on the X line: out, p, n.
        "pins": [
            pin("out", RIGHT, 0.0),
            pin("p", LEFT, 1.27, ELECTRICAL),
            pin("n", LEFT, -1.27, ELECTRICAL),
        ],
    },
    {
        "name": "fc_facet",
        "ref": "FCT",
        "params": "type=X model=fc_facet reflectance=0.0",
        "box": (-5.08, 2.54, 2.54, -2.54),
        # A mirror: hatched face on the right, light in above and back out below.
        "art": [
            ("poly", [(1.27, 2.54), (1.27, -2.54)], "none"),
            ("poly", [(1.27, -2.54), (2.54, -1.27)], "none"),
            ("poly", [(1.27, -1.27), (2.54, 0.0)], "none"),
            ("poly", [(1.27, 0.0), (2.54, 1.27)], "none"),
            ("poly", [(1.27, 1.27), (2.54, 2.54)], "none"),
            ("poly", [(-1.524, 1.905), (0.635, 1.905), (0.254, 2.159)], "none"),
            ("poly", [(0.635, 1.905), (0.254, 1.651)], "none"),
            ("poly", [(0.635, -1.905), (-1.524, -1.905), (-1.143, -1.651)], "none"),
            ("poly", [(-1.524, -1.905), (-1.143, -2.159)], "none"),
        ],
        "pins": [pin("port", LEFT, 0.0)],
    },
    {
        "name": "fc_tw_ps",
        "ref": "TWPS",
        "params": "type=X model=fc_tw_ps l_um=3000 v_pi_l=0.012 n_m=4.2 z0=35 f_max=50G",
        "box": (-6.35, 3.81, 6.35, -3.81),
        # A travelling-wave electrode over a guide: the RF crosses the top and
        # the light crosses the bottom, with the taps that couple them between.
        # The picture is the topology, which is the thing a reader needs.
        "art": [
            ("poly", [(-6.35, 2.54), (6.35, 2.54)], "none"),
            ("poly", [(-6.35, -2.54), (6.35, -2.54)], "none"),
            ("poly", [(-3.81, 2.54), (-3.81, -2.54)], "none"),
            ("poly", [(-1.27, 2.54), (-1.27, -2.54)], "none"),
            ("poly", [(1.27, 2.54), (1.27, -2.54)], "none"),
            ("poly", [(3.81, 2.54), (3.81, -2.54)], "none"),
        ],
        # Positional order on the X line: in, out, rf_in, rf_out.
        "pins": [
            pin("in", LEFT, -2.54),
            pin("out", RIGHT, -2.54),
            pin("rf_in", LEFT, 2.54, ELECTRICAL),
            pin("rf_out", RIGHT, 2.54, ELECTRICAL),
        ],
    },
]

PIN_LENGTH = 2.54
FONT = "\t\t\t\t\t(font\n\t\t\t\t\t\t(size 1.27 1.27)\n\t\t\t\t\t)"


def num(v: float) -> str:
    """KiCad writes trailing-zero-free decimals; match it so --check is stable."""
    s = f"{v:.6f}".rstrip("0").rstrip(".")
    return s if s not in ("", "-0") else "0"


def render_property(name: str, value: str, y: float, hide: bool, justify_left: bool) -> str:
    out = [f'\t\t(property "{name}" "{value}"', f"\t\t\t(at -6.35 {num(y)} 0)",
           "\t\t\t(show_name no)", "\t\t\t(do_not_autoplace no)"]
    if hide:
        out.append("\t\t\t(hide yes)")
    out += ["\t\t\t(effects", "\t\t\t\t(font", "\t\t\t\t\t(size 1.27 1.27)", "\t\t\t\t)"]
    if justify_left:
        out.append("\t\t\t\t(justify left)")
    out += ["\t\t\t)", "\t\t)"]
    return "\n".join(out)


def render_poly(points: List[Tuple[float, float]], fill: str) -> str:
    pts = " ".join(f"(xy {num(x)} {num(y)})" for x, y in points)
    return "\n".join([
        "\t\t\t(polyline", "\t\t\t\t(pts", f"\t\t\t\t\t{pts}", "\t\t\t\t)",
        "\t\t\t\t(stroke", "\t\t\t\t\t(width 0)", "\t\t\t\t\t(type default)", "\t\t\t\t)",
        "\t\t\t\t(fill", f"\t\t\t\t\t(type {fill})", "\t\t\t\t)", "\t\t\t)",
    ])


def render_pin(spec: dict, number: int, box: Tuple[float, float, float, float]) -> str:
    x1, _, x2, _ = box
    if spec["side"] == LEFT:
        x, rot = x1 - PIN_LENGTH, 0
    else:
        x, rot = x2 + PIN_LENGTH, 180
    return "\n".join([
        f"\t\t\t(pin {spec['kind']} line",
        f"\t\t\t\t(at {num(x)} {num(spec['y'])} {rot})",
        f"\t\t\t\t(length {num(PIN_LENGTH)})",
        f'\t\t\t\t(name "{spec["name"]}"', "\t\t\t\t\t(effects", FONT, "\t\t\t\t\t)", "\t\t\t\t)",
        f'\t\t\t\t(number "{number}"', "\t\t\t\t\t(effects", FONT, "\t\t\t\t\t)", "\t\t\t\t)",
        "\t\t\t)",
    ])


def render_symbol(spec: dict) -> str:
    name, box = spec["name"], spec["box"]
    x1, y1, x2, y2 = box
    lines = [
        f'\t(symbol "{name}"',
        "\t\t(exclude_from_sim no)", "\t\t(in_bom yes)", "\t\t(on_board yes)",
        "\t\t(in_pos_files yes)",
        render_property("Reference", spec["ref"], 3.81, False, False),
        render_property("Value", "", 0.508, False, False),
        render_property("Footprint", "", 0.508, True, False),
        render_property("Datasheet", "", 0.508, True, False),
        render_property("Description", "", 0.508, True, False),
        render_property("Sim.Device", "SPICE", -3.81, True, True),
        render_property("Sim.Params", spec["params"], -6.35, True, True),
        f'\t\t(symbol "{name}_1_1"',
        "\t\t\t(rectangle",
        f"\t\t\t\t(start {num(x1)} {num(y1)})",
        f"\t\t\t\t(end {num(x2)} {num(y2)})",
        "\t\t\t\t(stroke", "\t\t\t\t\t(width 0)", "\t\t\t\t\t(type default)", "\t\t\t\t)",
        "\t\t\t\t(fill", "\t\t\t\t\t(type none)", "\t\t\t\t)", "\t\t\t)",
    ]
    for kind, points, fill in spec["art"]:
        assert kind == "poly", kind
        lines.append(render_poly(points, fill))
    for i, p in enumerate(spec["pins"], start=1):
        lines.append(render_pin(p, i, box))
    lines += ["\t\t)", "\t\t(embedded_fonts no)", "\t)"]
    return "\n".join(lines)


def split_top_level(text: str) -> Tuple[str, List[Tuple[str, str]], str]:
    """(header, [(symbol_name, block)], footer) — one tab of indent is top level."""
    lines = text.split("\n")
    header: List[str] = []
    blocks: List[Tuple[str, str]] = []
    cur: List[str] | None = None
    cur_name = ""
    footer: List[str] = []
    for line in lines:
        if line.startswith('\t(symbol "'):
            if cur is not None:
                blocks.append((cur_name, "\n".join(cur)))
            cur_name = line.split('"')[1]
            cur = [line]
        elif cur is not None:
            cur.append(line)
            if line == "\t)":
                blocks.append((cur_name, "\n".join(cur)))
                cur = None
        elif blocks:
            footer.append(line)
        else:
            header.append(line)
    if cur is not None:
        raise SystemExit("unbalanced (symbol …) block — refusing to rewrite the library")
    return "\n".join(header), blocks, "\n".join(footer)


def rebuild(text: str) -> str:
    """Replace generated symbols in place; insert new ones alphabetically.

    Existing order is otherwise preserved — sorting the whole file would work
    too, but it would bury two new symbols in a diff that moved everything.
    """
    header, blocks, footer = split_top_level(text)
    generated = {s["name"]: render_symbol(s) for s in SYMBOLS}
    pending = sorted(n for n in generated if n not in {b[0] for b in blocks})
    out: List[str] = []
    for name, block in blocks:
        while pending and pending[0] < name:
            out.append(generated[pending.pop(0)])
        out.append(generated.get(name, block))
    out.extend(generated[n] for n in pending)
    return "\n".join([header] + out + [footer])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="exit 1 if the library is stale")
    ap.add_argument("--print", metavar="NAME", help="print one symbol and exit")
    ap.add_argument("--lib", type=Path, default=LIB)
    args = ap.parse_args()

    if args.print:
        for spec in SYMBOLS:
            if spec["name"] == args.print:
                print(render_symbol(spec))
                return 0
        print(f"no spec for {args.print}; have "
              f"{', '.join(s['name'] for s in SYMBOLS)}", file=sys.stderr)
        return 1

    old = args.lib.read_text()
    new = rebuild(old)
    if old == new:
        print(f"{args.lib.name}: up to date")
        return 0
    if args.check:
        print(f"{args.lib.name}: STALE — run scripts/gen_kicad_symbols.py", file=sys.stderr)
        return 1
    args.lib.write_text(new)
    print(f"{args.lib.name}: wrote {', '.join(s['name'] for s in SYMBOLS)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
