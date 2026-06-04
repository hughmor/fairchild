#!/usr/bin/env python3
"""
kicad_fairchild.py — one-button KiCad → fairchild driver.

Wraps the whole schematic→simulation pipeline behind a single command:

    KiCad schematic (.kicad_sch)
        │  kicad-cli sch export netlist --format spice
        ▼
    KiCad SPICE export (.cir)
        │  kicad_to_fairchild.py   (transpile + .optical_port + analysis)
        ▼
    fairchild netlist (run_*.sp)
        │  fairchild -f            (optional, with --run)
        ▼
    results (CSV / rawfile)

This is the engine behind the "Run fairchild" button described in
`kicad_integration.md`. It runs three ways:

1. Standalone on a schematic (what a developer types):
       python3 scripts/kicad_fairchild.py my.kicad_sch --tran "5n 2u" --run

2. As a KiCad *BOM / netlist-exporter* generator (the version-stable schematic
   "button" hook — KiCad 7/8/9/10). eeschema calls the generator with the
   intermediate netlist path; this script detects the XML, recovers the source
   schematic from <design source=…>, and runs the pipeline:
       python3 scripts/kicad_fairchild.py "%I" -o "%O" --tran "5n 2u"
   See `kicad_integration.md` for the exact dialog setup.

3. On an already-exported SPICE netlist (skips kicad-cli entirely — handy when
   KiCad isn't installed on the machine doing the simulation):
       python3 scripts/kicad_fairchild.py my_export.cir --tran "5n 2u" --run

kicad-cli is located automatically (PATH first, then the standard KiCad install
locations on macOS / Linux / Windows); override with --kicad-cli.
"""
import argparse
import os
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONVERTER = HERE / "kicad_to_fairchild.py"


# ---------------------------------------------------------------------------
# Locating kicad-cli
# ---------------------------------------------------------------------------

def _candidate_cli_paths():
    """Standard kicad-cli locations across platforms / installed versions."""
    yield from (shutil.which(n) for n in ("kicad-cli", "kicad-cli.exe"))
    # macOS nightly / dev bundles (named KiCad-dev.app, KiCad-nightly.app).
    for app in ("KiCad-dev", "KiCad-nightly"):
        yield f"/Applications/{app}.app/Contents/MacOS/kicad-cli"
    # macOS application bundles (newest versions first).
    for ver in ("", "_11", "_10", "_9", "_8", "_7"):
        yield f"/Applications/KiCad{ver}/KiCad.app/Contents/MacOS/kicad-cli"
    # Linux package installs.
    yield from ("/usr/bin/kicad-cli", "/usr/local/bin/kicad-cli")
    # Windows default install.
    for pf in (os.environ.get("ProgramFiles", r"C:\Program Files"),):
        for ver in ("9.0", "8.0", "7.0"):
            yield rf"{pf}\KiCad\{ver}\bin\kicad-cli.exe"


def find_kicad_cli(override: str | None) -> str | None:
    if override:
        return override if Path(override).exists() else None
    for cand in _candidate_cli_paths():
        if cand and Path(cand).exists():
            return cand
    return None


# ---------------------------------------------------------------------------
# Pipeline stages
# ---------------------------------------------------------------------------

def schematic_from_intermediate_netlist(xml_path: Path) -> Path | None:
    """KiCad's BOM/intermediate netlist is XML with <design source="…sch">.
    Recover the source schematic so we can ask kicad-cli for a SPICE export."""
    try:
        root = ET.parse(xml_path).getroot()
    except ET.ParseError:
        return None
    design = root.find("design")
    if design is None:
        return None
    src = design.findtext("source")
    return Path(src) if src else None


def export_spice(sch: Path, cli: str, verbose: bool) -> Path:
    """kicad-cli sch export netlist --format spice → .cir next to the schematic."""
    out = sch.with_suffix(".cir")
    cmd = [cli, "sch", "export", "netlist", "--format", "spice", "-o", str(out), str(sch)]
    if verbose:
        print(f"[kicad_fairchild] {' '.join(cmd)}", file=sys.stderr)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"kicad-cli netlist export failed:\n{proc.stderr or proc.stdout}")
    return out


def transpile(cir: Path, out: Path, passthru: list[str], verbose: bool) -> Path:
    """Run kicad_to_fairchild.py on the SPICE export → fairchild netlist."""
    cmd = [sys.executable, str(CONVERTER), str(cir), "-o", str(out), *passthru]
    if verbose:
        cmd.append("-v")
        print(f"[kicad_fairchild] {' '.join(cmd)}", file=sys.stderr)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.stderr:
        print(proc.stderr, file=sys.stderr, end="")
    if proc.returncode != 0:
        raise RuntimeError(f"kicad_to_fairchild.py failed:\n{proc.stderr or proc.stdout}")
    return out


def find_fairchild() -> str | None:
    for cand in (
        HERE.parent / "target" / "release" / "fairchild",
        HERE.parent / "target" / "debug" / "fairchild",
    ):
        if cand.exists():
            return str(cand)
    return shutil.which("fairchild")


def looks_like_spice(path: Path) -> bool:
    return path.suffix.lower() in (".cir", ".sp", ".net", ".spice")


def looks_like_intermediate_xml(path: Path) -> bool:
    if path.suffix.lower() != ".xml":
        return False
    try:
        head = path.read_text(errors="ignore")[:400]
    except OSError:
        return False
    return "<export" in head or "<design" in head


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_pipeline(args) -> int:
    inp = Path(args.input)
    if not inp.exists():
        print(f"error: input not found: {inp}", file=sys.stderr)
        return 2

    # Resolve the input down to a SPICE export (.cir), invoking kicad-cli if the
    # input is a schematic or an intermediate netlist.
    sch: Path | None = None
    if inp.suffix.lower() == ".kicad_sch":
        sch = inp
    elif looks_like_intermediate_xml(inp):
        sch = schematic_from_intermediate_netlist(inp)
        if sch is None or not sch.exists():
            print(f"error: could not locate source schematic from {inp}", file=sys.stderr)
            return 2

    if sch is not None:
        cli = find_kicad_cli(args.kicad_cli)
        if cli is None:
            print("error: kicad-cli not found. Install KiCad 7+ or pass --kicad-cli PATH.\n"
                  "       (Or export the SPICE netlist yourself and pass the .cir.)",
                  file=sys.stderr)
            return 3
        cir = export_spice(sch, cli, args.verbose)
    elif looks_like_spice(inp):
        cir = inp
    else:
        print(f"error: unrecognised input type {inp.suffix!r}; expected .kicad_sch, "
              f".xml (intermediate netlist), or a SPICE .cir/.sp", file=sys.stderr)
        return 2

    # Output fairchild netlist path.
    if args.output:
        out = Path(args.output)
        if out.suffix == "":  # eeschema passes %O without an extension
            out = out.with_suffix(".fairchild.sp")
    else:
        out = cir.with_name(f"run_{cir.stem}.sp")

    # Pass analysis/options flags straight through to the converter.
    passthru: list[str] = []
    if args.op:
        passthru.append("--op")
    if args.tran:
        passthru += ["--tran", args.tran]
    if args.ac:
        passthru += ["--ac", args.ac]
    if args.method:
        passthru += ["--method", args.method]
    for kv in args.opt or []:
        passthru += ["--opt", kv]

    transpile(cir, out, passthru, args.verbose)
    print(f"fairchild netlist → {out}")

    if args.run:
        fc = find_fairchild()
        if fc is None:
            print("warning: --run requested but fairchild binary not found "
                  "(build with `cargo build --release`).", file=sys.stderr)
            return 0
        cmd = [fc, "-f", str(out)]
        if args.probe:
            cmd += ["--probe", args.probe]
        if args.sim_output:
            cmd += ["--output", args.sim_output]
        if args.verbose:
            print(f"[kicad_fairchild] {' '.join(cmd)}", file=sys.stderr)
        proc = subprocess.run(cmd)
        return proc.returncode
    return 0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="KiCad schematic (.kicad_sch), intermediate netlist "
                                   "(.xml), or SPICE export (.cir/.sp)")
    ap.add_argument("-o", "--output", help="fairchild netlist output path")
    ap.add_argument("--op", action="store_true", help="emit a .op analysis")
    ap.add_argument("--tran", metavar="ARGS", help='emit .tran, e.g. --tran "5n 2u"')
    ap.add_argument("--ac", metavar="ARGS", help='emit .ac, e.g. --ac "dec 20 1 1G"')
    ap.add_argument("--method", choices=["be", "tr", "gear"], help="integration method")
    ap.add_argument("--opt", action="append", metavar="KEY=VAL",
                    help="extra .options token (repeatable), e.g. --opt waveguide_delay=1")
    ap.add_argument("--run", action="store_true", help="run fairchild on the result")
    ap.add_argument("--probe", help="probe expression passed to fairchild --probe")
    ap.add_argument("--sim-output", help="fairchild --output path (CSV/rawfile)")
    ap.add_argument("--kicad-cli", help="path to kicad-cli (auto-detected if omitted)")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()
    sys.exit(run_pipeline(args))


if __name__ == "__main__":
    main()
