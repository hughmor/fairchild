"""Write simulation results back onto the live KiCad schematic.

The interesting one: `SCH_BITMAP_T` is in `API_HANDLER_SCH::s_allowedTypes` and
`SchematicImage` carries raw `bytes image_data`, so a rendered matplotlib figure
can be dropped straight onto the drawing. A resonance sweep can live next to the
ring it belongs to.

Everything a run writes goes into one `SCH_GROUP` named `fairchild:<tag>`, inside
a single `BeginCommit`/`EndCommit`. That makes the whole annotation one Ctrl-Z,
and `clear()` can find and remove it by name without touching anything you drew.

    from fairchild.kicad import connect
    from fairchild.kicad.annotate import Annotation

    sch = connect()
    a = Annotation(sch, tag="op")
    a.text((120, 40), "V(w11) = 1.83 V")
    a.figure((150, 60), fig, scale=0.4)
    a.commit()          # one undo step
    a.clear()           # or take it back off
"""
from __future__ import annotations

import io
import json
from dataclasses import dataclass, field
from pathlib import Path

from google.protobuf.any_pb2 import Any as PbAny

from fairchild.kicad._proto.common.commands.editor_commands_pb2 import (
    BeginCommit, BeginCommitResponse, CommitAction, CreateItems,
    CreateItemsResponse, DeleteItems, DeleteItemsResponse, EndCommit,
    EndCommitResponse, GetItems, GetItemsResponse, UpdateItems,
    UpdateItemsResponse)
from fairchild.kicad._proto.common.types.base_types_pb2 import KIID, Vector2
from fairchild.kicad._proto.common.types.enums_pb2 import KiCadObjectType
from fairchild.kicad._proto.schematic.schematic_types_pb2 import (
    Group, SchematicImage, SchematicSymbolInstance, SchematicText)

MM = 1_000_000  # KiCad API distances are nanometres

#: Group-name prefix. Anything under it is ours and safe to delete; anything
#: else is the user's and must never be touched.
PREFIX = "fairchild"


def _pack(msg) -> PbAny:
    any_ = PbAny()
    any_.Pack(msg)
    return any_


@dataclass
class Annotation:
    """A batch of schematic items to be committed as one undoable group."""

    sch: object
    tag: str = "run"
    #: Sheet instance name to annotate; None targets whatever sheet the
    #: document specifier points at (the root).
    sheet: str | None = None
    items: list[PbAny] = field(default_factory=list)

    @property
    def group_name(self) -> str:
        return f"{PREFIX}:{self.tag}"

    # ── building ──────────────────────────────────────────────────────────
    def text(self, at: tuple[float, float], content: str, size_mm: float = 1.27,
             bold: bool = False) -> "Annotation":
        t = SchematicText()
        t.text.text = content
        t.text.position.CopyFrom(Vector2(x_nm=int(at[0] * MM), y_nm=int(at[1] * MM)))
        t.text.attributes.size.CopyFrom(
            Vector2(x_nm=int(size_mm * MM), y_nm=int(size_mm * MM)))
        t.text.attributes.bold = bold
        # Annotations are outputs, not stimulus — never let one feed the netlist.
        t.exclude_from_sim = True
        self.items.append(_pack(t))
        return self

    def image(self, at: tuple[float, float], png: bytes,
              scale: float = 1.0) -> "Annotation":
        img = SchematicImage()
        img.position.CopyFrom(Vector2(x_nm=int(at[0] * MM), y_nm=int(at[1] * MM)))
        img.image_scale.value = scale
        img.image_data = png
        self.items.append(_pack(img))
        return self

    def figure(self, at: tuple[float, float], fig, scale: float = 1.0,
               dpi: int = 150) -> "Annotation":
        """Embed a matplotlib figure. This is the fun one."""
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
        return self.image(at, buf.getvalue(), scale)

    # ── committing ────────────────────────────────────────────────────────
    def _header(self, req) -> None:
        req.header.document.CopyFrom(self.sch.doc)
        if self.sheet is not None:
            paths = [sh.path for sh in self.sch.sheet_params.values()
                     if sh.name == self.sheet]
            if not paths:
                raise KeyError(f"no sheet named {self.sheet!r}")
            del req.header.document.sheet_path.path[:]
            for kiid in paths[0]:
                req.header.document.sheet_path.path.add().value = kiid

    def commit(self, label: str | None = None) -> list[str]:
        """Create every queued item plus their group, as one undo step."""
        if not self.items:
            return []
        k = self.sch.k
        token = self._begin()
        ok = False
        try:
            req = CreateItems()
            self._header(req)
            req.items.extend(self.items)
            resp = k.send(req, CreateItemsResponse)
            made = [r.item.type_url for r in resp.created_items]
            ids = [_created_id(r) for r in resp.created_items]
            ids = [i for i in ids if i]

            # Group them so the annotation is one object on the canvas. Membership
            # cannot be set at creation: SCH_GROUP::Deserialize clears m_items and
            # then returns false on `if( !schematic )`, and a group being created
            # isn't attached to the schematic yet — so the name survives and the
            # members are silently dropped. Setting them in a second Update, once
            # both group and items exist, gives ResolveItem something to find.
            if ids:
                g = Group()
                g.name = self.group_name
                greq = CreateItems()
                self._header(greq)
                greq.items.append(_pack(g))
                gresp = k.send(greq, CreateItemsResponse)
                gid = next((_created_id(r) for r in gresp.created_items), "")
                if gid:
                    g.id.value = gid
                    for i in ids:
                        g.items.add().value = i
                    ureq = UpdateItems()
                    self._header(ureq)
                    ureq.items.append(_pack(g))
                    k.send(ureq, UpdateItemsResponse)
            self._remember(ids + ([gid] if ids and gid else []))
            ok = True
        finally:
            # Drop rather than commit on the way out of an exception, so a failed
            # annotation never leaves half a group on someone's schematic.
            self._end(token, label or f"fairchild: annotate ({self.tag})", ok)

        self.items.clear()
        return made

    def _begin(self) -> KIID:
        req = BeginCommit()
        self._header(req)
        return self.sch.k.send(req, BeginCommitResponse).id

    def _end(self, token: KIID, message: str, commit: bool = True) -> None:
        end = EndCommit()
        end.id.CopyFrom(token)
        end.action = CommitAction.CMA_COMMIT if commit else CommitAction.CMA_DROP
        end.message = message
        self._header(end)
        self.sch.k.send(end, EndCommitResponse)

    # ── cleanup ───────────────────────────────────────────────────────────
    def _ledger(self) -> Path:
        """Sidecar recording what we created, so cleanup never has to guess."""
        return Path(self.sch.doc.project.path) / ".fairchild-annotations.json"

    def _remember(self, ids: list[str]) -> None:
        path, data = self._ledger(), {}
        if path.exists():
            try:
                data = json.loads(path.read_text())
            except ValueError:
                data = {}
        data.setdefault(self.tag, [])
        data[self.tag] = sorted(set(data[self.tag]) | set(ids))
        try:
            path.write_text(json.dumps(data, indent=1))
        except OSError:
            pass  # cleanup still works off the text prefix; don't fail a commit

    def clear(self) -> int:
        """Remove annotations carrying this tag.

        Three sources, because none alone is sufficient: the ledger (catches
        bitmaps, which carry no identifying text), the `fairchild:` text prefix
        (catches anything written before a ledger existed), and the group name.
        Never deletes an item that fails all three — the user's own work is
        untouchable.
        """
        k = self.sch.k
        wanted = set()
        if self._ledger().exists():
            try:
                wanted = set(json.loads(self._ledger().read_text())
                             .get(self.tag, []))
            except ValueError:
                pass

        req = GetItems(types=[KiCadObjectType.KOT_SCH_TEXT,
                              KiCadObjectType.KOT_SCH_BITMAP,
                              KiCadObjectType.KOT_SCH_GROUP])
        req.header.document.CopyFrom(self.sch.doc)
        req.header.document.ClearField("sheet_path")
        victims: list[str] = []
        seen: set[str] = set()

        def take(kiid: str) -> None:
            # A shared sheet file reports its items once per instance; dedupe or
            # the same KIID is queued for deletion eight times.
            if kiid and kiid not in seen:
                seen.add(kiid)
                victims.append(kiid)

        for any_it in k.send(req, GetItemsResponse).items:
            t = SchematicText()
            if any_it.Unpack(t) and t.text.text:
                if t.id.value in wanted or t.text.text.startswith(PREFIX + ":"):
                    take(t.id.value)
                continue
            img = SchematicImage()
            if any_it.Unpack(img) and img.image_data:
                if img.id.value in wanted:
                    take(img.id.value)
                continue
            g = Group()
            if any_it.Unpack(g) and (g.name == self.group_name
                                     or g.id.value in wanted):
                take(g.id.value)

        if not victims:
            return 0
        token = self._begin()
        ok = False
        try:
            dreq = DeleteItems()
            self._header(dreq)
            for v in victims:
                dreq.item_ids.add().value = v
            k.send(dreq, DeleteItemsResponse)
            ok = True
        finally:
            self._end(token, f"fairchild: clear annotation ({self.tag})", ok)
        if ok and self._ledger().exists():
            try:
                data = json.loads(self._ledger().read_text())
                data.pop(self.tag, None)
                self._ledger().write_text(json.dumps(data, indent=1))
            except (OSError, ValueError):
                pass
        return len(victims)


def _created_id(result) -> str:
    """Pull the new item's KIID out of an ItemCreationResult, shape-agnostically.

    The result wraps the created item in an Any; every schematic item puts its
    KIID in field 1, but the concrete type varies, so decode generically rather
    than branching on all 16 possibilities.
    """
    for msg in (SchematicText, SchematicImage, Group):
        probe = msg()
        if result.item.Unpack(probe) and probe.id.value:
            return probe.id.value
    return ""


def demo(tag: str = "params") -> None:
    """Annotate each PN modulator sheet with its resolved ring length, and paste
    a plot of the ladder across all of them onto the parent sheet.

    Deliberately needs no simulation: it annotates data already extracted from
    the schematic, so it exercises the whole write-back path (text + bitmap +
    group + commit) without depending on a circuit that converges.
    """
    import matplotlib
    matplotlib.use("Agg")  # never pop a window; this runs headless by design
    import matplotlib.pyplot as plt

    from fairchild.kicad import connect

    sch = connect()
    rings = {}
    for s in sch.symbols:
        if s.model == "fc_pn_th_ps" and s.path in sch.sheet_params:
            sheet = sch.sheet_params[s.path]
            lm = s.param_dict.get("l_m")
            if lm:
                rings[sheet.name] = lm
    if not rings:
        print("no fc_pn_th_ps with an l_m found; nothing to annotate")
        return

    def um(v: str) -> float:
        return float(v.rstrip("u")) if v.endswith("u") else float(v) * 1e6

    order = sorted(rings, key=lambda n: (len(n), n))
    vals = [um(rings[n]) for n in order]

    # Per-sheet text, one commit each so a single sheet can be undone alone.
    for name, v in zip(order, vals):
        a = Annotation(sch, tag=tag, sheet=name)
        a.clear()
        a.text((5, 5), f"fairchild: l_m = {v:g} um", size_mm=1.5, bold=True)
        a.commit()
        print(f"  annotated {name}: l_m = {v:g} um")

    fig, ax = plt.subplots(figsize=(4.2, 2.4))
    ax.plot(range(1, len(vals) + 1), vals, "o-", lw=1.4, ms=5)
    ideal = [vals[0] - 0.03 * i for i in range(len(vals))]
    ax.plot(range(1, len(vals) + 1), ideal, "k--", lw=1,
            label="uniform -0.03 um ladder")
    ax.set_xlabel("PN modulator", fontsize=8)
    ax.set_ylabel("ring half length (um)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=7)
    ax.grid(alpha=.3)
    fig.tight_layout()

    parent = Annotation(sch, tag=tag, sheet="Giona Chip")
    parent.clear()
    parent.text((5, 5), "fairchild: ring-length ladder", size_mm=2.0, bold=True)
    parent.figure((5, 12), fig, scale=0.5)
    made = parent.commit()
    print(f"pasted ladder plot onto 'Giona Chip' ({len(made)} items)")
    print(f"undo with Ctrl-Z, or Annotation(sch, tag={tag!r}, sheet=...).clear()")


if __name__ == "__main__":
    demo()


# ── selection ─────────────────────────────────────────────────────────────────
def selection(sch, types=None) -> list:
    """Currently selected schematic items, as unpacked protos.

    Annotating everything on a chip this size is unreadable, so the selection is
    the natural scope: mark what you are looking at.

    Not available yet: API_HANDLER_SCH registers no GetSelection handler in
    10.99 (only the PCB side does), so this returns [] rather than throwing.
    Use `sch.current_sheet` for scoping in the meantime — it tracks which sheet
    the editor is showing, which is a better signal for a toolbar action anyway.
    """
    from fairchild.kicad import ApiError
    from fairchild.kicad._proto.common.commands.editor_commands_pb2 import (
        GetSelection, SelectionResponse)
    req = GetSelection()
    req.header.document.CopyFrom(sch.doc)
    req.header.document.ClearField("sheet_path")
    for t in types or ():
        req.types.append(t)
    try:
        resp = sch.k.send(req, SelectionResponse)
    except ApiError as e:
        if "UNHANDLED" not in str(e):
            raise
        return []
    out = []
    for any_it in resp.items:
        si = SchematicSymbolInstance()
        if any_it.Unpack(si) and si.id.value:
            out.append(si)
    return out


def selected_symbols(sch) -> list:
    """Selection intersected with the symbols we already know about.

    Matching on (path, kiid) rather than kiid alone, for the usual reason: a
    reused sheet repeats its symbol KIIDs across instances.
    """
    want = {(tuple(k.value for k in si.path.path), si.id.value)
            for si in selection(sch)}
    return [s for s in sch.symbols if (s.path, s.kiid) in want]


def op_labels(sch, result, tag: str = "op", size_mm: float = 1.0,
              offset_mm: tuple[float, float] = (1.5, -1.5),
              symbols=None) -> int:
    """Label each selected symbol's pins with its node voltage from `result`.

    Uses the theme's own `op_voltages` colour so the labels agree with KiCad's
    native ${OP} annotations. One text per pin whose net actually appears in the
    result, so unsolved or optical-only nets are quietly skipped rather than
    littering the sheet with blanks.
    """
    from fairchild.kicad import _sanitise
    syms = symbols if symbols is not None else selected_symbols(sch)
    if not syms:
        return 0
    sigs = set(result.signals())
    # Group by sheet: items live in a sheet's file, so each sheet is its own call.
    by_sheet: dict[tuple, list] = {}
    for s in syms:
        by_sheet.setdefault(s.path, []).append(s)

    written = 0
    for path, group in by_sheet.items():
        sheet = sch.sheet_params.get(path)
        a = Annotation(sch, tag=tag, sheet=sheet.name if sheet else None)
        for s in group:
            for p in s.sorted_pins():
                net = _sanitise(p.net)
                key = f"V({net.lower()})"
                if key not in sigs:
                    continue
                v = float(result[key][0])
                at = (p.position_mm[0] + offset_mm[0],
                      p.position_mm[1] + offset_mm[1])
                a.text(at, f"{v:.4g} V", size_mm=size_mm)
                written += 1
        if a.items:
            a.commit(f"fairchild: .op labels ({len(a.items)} nodes)")
    return written


# ── per-instance values on a reused sheet ────────────────────────────────────
def set_sheet_field(sch, sheet_name: str, field: str, value: str) -> None:
    """Set one `user_field` on ONE sheet instance, preserving the others.

    This is how a reused sheet gets per-instance annotation. Items belong to a
    sheet's *file*, so a literal text placed inside `pn_mrm_mod.kicad_sch` shows
    up in all 8 instances. But a text reading `${FC_RESULT}` resolves through
    SCH_SHEET::ResolveTextVar per sheet path, so writing a different
    `FC_RESULT` on each of the 8 sheet symbols — which ARE distinct items, on
    the parent's screen — renders a different value in each instance. Same trick
    KiCad uses for ${OP}, with us supplying the value.

    Read-modify-write is mandatory: SCH_SHEET::Deserialize clears m_fields and
    rebuilds it from the proto, so sending only the new field would delete
    RING_HALF_LENGTH and friends.
    """
    from fairchild.kicad._proto.common.commands.editor_commands_pb2 import (
        UpdateItems, UpdateItemsResponse)
    from fairchild.kicad._proto.schematic.schematic_types_pb2 import SheetSymbol

    req = GetItems(types=[KiCadObjectType.KOT_SCH_SHEET])
    req.header.document.CopyFrom(sch.doc)
    req.header.document.ClearField("sheet_path")
    target = None
    for any_it in sch.k.send(req, GetItemsResponse).items:
        sh = SheetSymbol()
        if any_it.Unpack(sh) and sh.name_field.text.text == sheet_name:
            target = sh
            break
    if target is None:
        raise KeyError(f"no sheet instance named {sheet_name!r}")

    existing = {f.name: f for f in target.user_fields}
    if field in existing:
        existing[field].text.text = value
    else:
        f = target.user_fields.add()
        f.name = field
        f.text.text = value
        f.visible = False          # the ${} text shows it; the field itself needn't
        # Park it on the sheet origin; invisible, but a wild position is still rude.
        f.text.position.CopyFrom(target.name_field.text.position)

    a = Annotation(sch, tag="field")
    token = a._begin()
    ok = False
    try:
        ureq = UpdateItems()
        ureq.header.document.CopyFrom(sch.doc)
        # A sheet symbol lives on its PARENT's screen, so the update has to be
        # scoped there or the handler answers ISC_NONEXISTENT for a UUID it just
        # gave us. SheetSymbol.path is the container's path, which is exactly
        # what belongs here — the one place not to clear sheet_path.
        del ureq.header.document.sheet_path.path[:]
        for kiid in target.path.path:
            ureq.header.document.sheet_path.path.add().CopyFrom(kiid)
        ureq.items.append(_pack(target))
        resp = sch.k.send(ureq, UpdateItemsResponse)
        bad = [r.status.error_message for r in resp.updated_items
               if r.status.code != 1]  # 1 == ISC_OK
        if bad:
            raise RuntimeError(f"sheet field update rejected: {bad[0]}")
        ok = True
    finally:
        a._end(token, f"fairchild: set {sheet_name}.{field}", ok)
