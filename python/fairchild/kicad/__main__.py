"""REPL entrypoint: `python -i -m fairchild.kicad` drops you in with `sch` live.

    $ .venv/bin/python -i -m fairchild.kicad
    >>> print(sch.report())
    >>> ckt = sch.circuit(); r = ckt.run("op")
    >>> sch.refresh()          # after editing in KiCad

`--demo` runs the self-check instead, `--deck` prints the netlist and exits.
"""
import sys

from fairchild.kicad import connect, demo

if "--demo" in sys.argv:
    demo()
    sys.exit(0)

sch = connect()
if "--deck" in sys.argv:
    print(sch.deck())
    sys.exit(0)

print(sch.report())
for w in sch.warnings:
    print(f"  ! {w}", file=sys.stderr)
print("\nsch = live schematic   sch.refresh()   sch.deck()   sch.circuit()")
