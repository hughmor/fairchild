"""Style matplotlib with the schematic editor's own colour theme.

The IPC API exposes no colour settings, but the active theme is named in
`eeschema.json` (`appearance.color_theme`) and the theme itself is a plain JSON
file, so a plot written back into the schematic can match the canvas it lands
on. `op_voltages` / `op_currents` are in there too — the exact colours KiCad
uses for its own `${OP}` annotations — so our node labels can agree with them.

    from fairchild.kicad.theme import styled, palette
    with styled():
        fig, ax = plt.subplots()      # dark canvas, theme line cycle
"""
from __future__ import annotations

import json
import re
from contextlib import contextmanager
from pathlib import Path

_RGB = re.compile(r"rgba?\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*(?:,\s*([\d.]+)\s*)?\)")

#: Where KiCad keeps per-version config on each platform. macOS first because
#: that is what is tested; the others are the documented locations.
_CONFIG_DIRS = (
    Path.home() / "Library/Preferences/kicad",           # macOS
    Path.home() / ".config/kicad",                       # Linux
    Path.home() / "AppData/Roaming/kicad",               # Windows
)


def _config_dir(version: str | None = None) -> Path | None:
    for base in _CONFIG_DIRS:
        if not base.is_dir():
            continue
        if version and (base / version).is_dir():
            return base / version
        # Highest version number wins, so a dev build beats an old stable one.
        subs = sorted((d for d in base.iterdir() if d.is_dir()
                       and re.match(r"^\d+\.\d+$", d.name)),
                      key=lambda d: [int(x) for x in d.name.split(".")])
        if subs:
            return subs[-1]
    return None


def active_theme_path(version: str | None = None) -> Path | None:
    """The JSON file for the theme eeschema is currently drawing with."""
    cfg = _config_dir(version)
    if cfg is None:
        return None
    try:
        appearance = json.loads((cfg / "eeschema.json").read_text())
        name = appearance.get("appearance", {}).get("color_theme", "")
    except (OSError, ValueError):
        return None
    if not name:
        return None
    p = Path(name)
    if p.is_absolute():          # 3rd-party themes are stored by full path
        return p if p.exists() else None
    for cand in (cfg / "colors" / f"{name}.json", cfg / "colors" / name):
        if cand.exists():
            return cand
    return None


def _hex(value: str) -> str:
    m = _RGB.match(value.strip())
    if not m:
        return value
    r, g, b = (int(m.group(i)) for i in (1, 2, 3))
    return f"#{r:02x}{g:02x}{b:02x}"


def palette(version: str | None = None) -> dict[str, str]:
    """{colour role: #rrggbb} for the schematic, plus `cycle` (a list).

    Empty if no theme can be found, so callers can fall back to matplotlib
    defaults rather than crash on a machine with no KiCad config.
    """
    path = active_theme_path(version)
    if path is None:
        return {}
    try:
        theme = json.loads(path.read_text())
    except (OSError, ValueError):
        return {}
    out = {k: _hex(v) for k, v in theme.get("schematic", {}).items()
           if isinstance(v, str)}
    out["cycle"] = [_hex(c) for c in theme.get("palette", []) if isinstance(c, str)]
    out["_name"] = theme.get("meta", {}).get("name", path.stem)
    return out


def mpl_style(version: str | None = None) -> dict:
    """rcParams matching the schematic canvas. Empty dict if no theme found."""
    p = palette(version)
    if not p:
        return {}
    bg = p.get("background", "#131218")
    fg = p.get("note", "#f8f8f0")
    grid = p.get("grid", "#716799")
    rc = {
        "figure.facecolor": bg,
        "axes.facecolor": p.get("sheet_background", bg),
        "savefig.facecolor": bg,
        "axes.edgecolor": p.get("component_outline", fg),
        "axes.labelcolor": fg,
        "text.color": fg,
        "xtick.color": fg,
        "ytick.color": fg,
        "grid.color": grid,
        "grid.alpha": 0.4,
        "legend.facecolor": p.get("component_body", bg),
        "legend.edgecolor": p.get("component_outline", fg),
        "axes.titlecolor": p.get("label_global", fg),
    }
    if p.get("cycle"):
        from matplotlib import cycler
        rc["axes.prop_cycle"] = cycler(color=p["cycle"])
    return rc


@contextmanager
def styled(version: str | None = None):
    """Temporarily draw with the schematic's theme."""
    import matplotlib.pyplot as plt
    rc = mpl_style(version)
    if not rc:
        yield {}
        return
    with plt.rc_context(rc):
        yield palette(version)


def demo() -> None:
    p = palette()
    if not p:
        print("no KiCad theme found")
        return
    print(f"theme: {p['_name']}  ({active_theme_path()})")
    for key in ("background", "sheet_background", "note", "wire", "grid",
                "op_voltages", "op_currents", "label_global", "pin"):
        print(f"  {key:<18} {p.get(key)}")
    print(f"  cycle ({len(p['cycle'])}): {p['cycle'][:6]}")


if __name__ == "__main__":
    demo()
