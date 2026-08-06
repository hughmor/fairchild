"""Entry point for the eeschema toolbar buttons.

Installed by `install_plugin.sh` into KiCad's user plugin path
(`~/Documents/KiCad/<ver>/plugins/fairchild/`). KiCad launches the action with
`KICAD_API_SOCKET` and `KICAD_API_TOKEN` already in the environment, so no
discovery is needed here.

Two actions:

  annotate-op   scope an .op to the selection's sheet, label its node voltages
  clear         remove every annotation fairchild has written

`annotate-op` deliberately works off the *selection*: labelling 2896 pins would
be unreadable, and scoping to one sheet is also the only way a design of this
size solves at all right now.
"""
from __future__ import annotations

import sys

from fairchild.kicad import connect
from fairchild.kicad.annotate import (PREFIX, Annotation, op_labels,
                                      selected_symbols)

TAGS = ("op", "params", "probe", "field", "run")


def annotate_op(scope: str | None = None) -> int:
    """Simulate one sheet and label its node voltages.

    Scope, in order of preference: an explicit name, then the selection (once
    KiCad implements GetSelection for the schematic), then the sheet the editor
    is currently showing. The last is the intended day-to-day path — descend
    into a block, press the button.
    """
    sch = connect()
    if scope:
        syms = sch.scoped(scope)
    else:
        syms = selected_symbols(sch)          # [] until KiCad implements it
        if not syms:
            sheet = sch.current_sheet
            if sheet is None:
                print("You are on the root sheet, which for a design this size "
                      "will not converge. Descend into a block (or pass a scope "
                      "name) and press the button again.\n"
                      f"available: {', '.join(sorted(sch.scopes()))}")
                return 1
            syms = sch.scoped(sheet.name)
            print(f"current sheet: {sheet.name} ({sheet.filename})")
    if not syms:
        print("Nothing to annotate here.")
        return 1

    sheets = {s.path for s in syms}
    if len(sheets) > 1:
        print(f"Scope spans {len(sheets)} sheets; annotating each in turn.")

    total = 0
    for path in sheets:
        sheet = sch.sheet_params.get(path)
        scope = sheet.name if sheet else None
        group = [s for s in syms if s.path == path]
        try:
            ckt = sch.circuit(scope=scope) if scope else sch.circuit()
        except Exception as e:
            print(f"{scope or 'root'}: could not build a deck — {e}")
            continue
        try:
            result = ckt.run("op")
        except Exception as e:
            print(f"{scope or 'root'}: .op did not converge — {e}")
            continue
        n = op_labels(sch, result, tag="op", symbols=group)
        print(f"{scope or 'root'}: labelled {n} node(s) on {len(group)} symbol(s)")
        total += n

    if not total:
        print("No node voltages matched — the selection may be optical-only, "
              "or its nets are outside the scoped deck.")
    return 0


def clear() -> int:
    sch = connect()
    removed = 0
    for tag in TAGS:
        removed += Annotation(sch, tag=tag).clear()
    print(f"removed {removed} item(s) written by {PREFIX}")
    return 0


def main(argv: list[str] | None = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    action = argv[0] if argv else "annotate-op"
    try:
        if action == "annotate-op":
            return annotate_op(argv[1] if len(argv) > 1 else None)
        if action == "clear":
            return clear()
    except Exception as e:  # a traceback in KiCad's log window helps nobody
        print(f"fairchild: {type(e).__name__}: {e}")
        return 1
    print(f"unknown action {action!r}; expected annotate-op or clear")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
