#!/usr/bin/env python3
"""fairchild_gui.py — a standalone desktop launcher for the KiCad → fairchild flow.

This is "Option C" from the KiCad-integration plan: a small control-panel app
that drives the existing pipeline and shows the results, with zero coupling to
KiCad's (PCB-only) IPC API and no need to fork eeschema. Pick a schematic or a
netlist, choose an analysis, hit Run — it does

    schematic (.kicad_sch / .xml)         native fairchild netlist (.sp)
        │  kicad-cli export                   │
        │  kicad_to_fairchild.py (transpile)  │  (run as-is)
        ▼                                     ▼
                    fairchild -f … --format csv
                              │
                              ▼
                 embedded matplotlib figure + log

Built on PyQt6 (already present) and the two existing scripts: `kicad_fairchild`
supplies the pipeline plumbing (kicad-cli location, schematic resolution,
transpile) and `fairchild_plot` supplies the plotters (which now accept a
`fig=` to draw into the embedded canvas).

Architecture: the pipeline + render logic lives in Qt-free module functions
(`run_to_csv`, `render`) so it is testable headless; the Qt classes are a thin
shell. `--selftest [INPUT]` exercises the whole wiring under the offscreen Qt
platform (no display needed).

Run it:
    python3 scripts/fairchild_gui.py
    python3 scripts/fairchild_gui.py examples/kicad_photonics/mrm_single_channel.kicad_sch
"""
from __future__ import annotations

import io
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

# The two sibling scripts. Import side-effect-free (fairchild_plot no longer
# locks a matplotlib backend at import; kicad_fairchild is pure functions).
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fairchild_plot as fp          # noqa: E402
import kicad_fairchild as kf         # noqa: E402


# ── Qt-free core ─────────────────────────────────────────────────────────────
# These functions do the real work and never touch Qt, so they can be unit
# tested without a display or an event loop.

#: Extensions that are KiCad exports / schematics → default to transpiling.
_TRANSPILE_EXTS = {".kicad_sch", ".xml", ".cir", ".net", ".spice"}
ANALYSES = ("op", "tran", "ac")


@dataclass
class RunOptions:
    input_path: str
    transpile: bool = True       # KiCad path (transpile + inject analysis)
    analysis: str = "op"         # one of ANALYSES (ignored when transpile is off)
    tran_args: str = "5n 2u"
    ac_args: str = "dec 20 1 1G"
    method: str = "default"      # default | be | tr | gear
    waveguide_delay: bool = False
    probe: str = ""              # display filter, comma-separated exact column names
    raw_quadratures: bool = False
    extra_opts: str = ""         # space-separated KEY=VAL tokens for fairchild --opt
    kicad_cli: str | None = None


def default_transpile_for(path: str) -> bool:
    """A native fairchild `.sp` carries its own analysis → run as-is; a KiCad
    export / schematic needs transpiling."""
    return Path(path).suffix.lower() in _TRANSPILE_EXTS


def _analysis_passthru(opts: RunOptions) -> list[str]:
    """Analysis flags for kicad_to_fairchild.py (transpile path only)."""
    if opts.analysis == "tran":
        return ["--tran", opts.tran_args]
    if opts.analysis == "ac":
        return ["--ac", opts.ac_args]
    return ["--op"]


def _fairchild_opt_tokens(opts: RunOptions) -> list[str]:
    """Solver-option overrides passed straight to `fairchild --opt`. Only the
    controls the user actually set are sent, so they override and an untouched
    control leaves the netlist's own `.options` alone."""
    toks: list[str] = []
    if opts.method and opts.method != "default":
        toks.append(f"method={opts.method}")
    if opts.waveguide_delay:
        toks.append("waveguide_delay=1")
    toks += opts.extra_opts.split()
    return toks


def run_to_csv(opts: RunOptions, log=lambda _m: None) -> str:
    """Drive the pipeline to a fairchild CSV string. `log(msg)` receives one
    progress line per stage. Raises kf.PipelineError / RuntimeError with a
    user-facing message (no stack trace) on any failure."""
    inp = Path(opts.input_path)
    if not opts.input_path:
        raise RuntimeError("Choose a schematic or netlist first.")

    if opts.transpile:
        cir = kf.resolve_to_spice(inp, opts.kicad_cli, verbose=False)
        if cir != inp:
            log(f"kicad-cli export → {cir.name}")
        netlist = Path(tempfile.gettempdir()) / f"{inp.stem}.fairchild.sp"
        kf.transpile(cir, netlist, _analysis_passthru(opts), verbose=False)
        log(f"transpiled ({opts.analysis}) → {netlist.name}")
    else:
        if not inp.exists():
            raise RuntimeError(f"input not found: {inp}")
        netlist = inp
        log(f"native netlist {inp.name} (using its own analysis directive)")

    fc = kf.find_fairchild()
    if fc is None:
        raise RuntimeError("fairchild binary not found — build it with "
                           "`cargo build --release`.")
    cmd = [fc, "-f", str(netlist), "--format", "csv"]
    for kv in _fairchild_opt_tokens(opts):
        cmd += ["--opt", kv]
    log("fairchild " + " ".join(cmd[1:]))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or "fairchild exited non-zero")
    if proc.stderr.strip():
        log(proc.stderr.strip())
    return proc.stdout


def render(csv_text: str, fig, opts: RunOptions, title: str | None = None) -> str:
    """Parse a fairchild CSV and draw it into `fig` (cleared first). Returns the
    detected analysis kind ('transient' | 'AC' | 'operating point')."""
    x_name, x, series = fp.read_csv(io.StringIO(csv_text))
    if opts.probe:
        keep = {s.strip() for s in opts.probe.split(",") if s.strip()}
        series = {k: v for k, v in series.items() if k in keep}
    xl = x_name.lower()
    if xl.startswith("time"):
        fp.plot_transient(x, series, opts.raw_quadratures, title, fig=fig)
        return "transient"
    if xl.startswith("freq"):
        fp.plot_ac(x, series, title, fig=fig)
        return "AC"
    fp.plot_op(series, opts.raw_quadratures, title, fig=fig)
    return "operating point"


# ── Qt shell ─────────────────────────────────────────────────────────────────
def _build_qt():
    """Import Qt + the matplotlib Qt canvas lazily, so the Qt-free core above
    can be imported and tested without Qt being importable at module load."""
    from PyQt6 import QtCore, QtWidgets
    from matplotlib.backends.backend_qtagg import (
        FigureCanvasQTAgg, NavigationToolbar2QT)
    from matplotlib.figure import Figure
    return QtCore, QtWidgets, FigureCanvasQTAgg, NavigationToolbar2QT, Figure


def make_app_and_window(initial_input: str | None = None):
    QtCore, QtWidgets, FigureCanvasQTAgg, NavigationToolbar2QT, Figure = _build_qt()

    class Worker(QtCore.QThread):
        done = QtCore.pyqtSignal(str)      # csv text
        failed = QtCore.pyqtSignal(str)    # user-facing message
        logged = QtCore.pyqtSignal(str)

        def __init__(self, opts: RunOptions):
            super().__init__()
            self._opts = opts

        def run(self):  # executes off the GUI thread
            try:
                csv = run_to_csv(self._opts, self.logged.emit)
                self.done.emit(csv)
            except Exception as exc:  # PipelineError / RuntimeError / anything
                self.failed.emit(str(exc))

    class MainWindow(QtWidgets.QMainWindow):
        def __init__(self):
            super().__init__()
            self.setWindowTitle("fairchild")
            self.resize(1180, 760)
            self._worker: Worker | None = None

            split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
            split.addWidget(self._build_controls(QtWidgets))
            split.addWidget(self._build_results(QtWidgets, QtCore,
                                                FigureCanvasQTAgg,
                                                NavigationToolbar2QT, Figure))
            split.setStretchFactor(1, 1)
            split.setSizes([330, 850])
            self.setCentralWidget(split)
            self.statusBar().showMessage("Ready.")

            if initial_input:
                self._set_input(initial_input)
            self._sync_enabled()

        # -- control panel -------------------------------------------------
        def _build_controls(self, W):
            panel = W.QWidget()
            v = W.QVBoxLayout(panel)

            # Input file
            box = W.QGroupBox("Input")
            form = W.QFormLayout(box)
            self.input_edit = W.QLineEdit()
            self.input_edit.setPlaceholderText("schematic or netlist…")
            browse = W.QPushButton("Browse…")
            browse.clicked.connect(self._on_browse)
            row = W.QHBoxLayout()
            row.addWidget(self.input_edit, 1)
            row.addWidget(browse)
            rw = W.QWidget(); rw.setLayout(row)
            form.addRow(rw)
            self.transpile_cb = W.QCheckBox("Transpile (KiCad export / schematic)")
            self.transpile_cb.setChecked(True)
            self.transpile_cb.toggled.connect(self._sync_enabled)
            form.addRow(self.transpile_cb)
            v.addWidget(box)

            # Analysis
            self.analysis_box = W.QGroupBox("Analysis")
            af = W.QFormLayout(self.analysis_box)
            self.analysis_combo = W.QComboBox()
            self.analysis_combo.addItems(["Operating point", "Transient", "AC sweep"])
            self.analysis_combo.currentIndexChanged.connect(self._sync_enabled)
            af.addRow("Type", self.analysis_combo)
            self.tran_edit = W.QLineEdit("5n 2u")
            af.addRow("Transient (step stop)", self.tran_edit)
            self.ac_edit = W.QLineEdit("dec 20 1 1G")
            af.addRow("AC (type n f0 f1)", self.ac_edit)
            v.addWidget(self.analysis_box)

            # Options
            opt = W.QGroupBox("Options")
            of = W.QFormLayout(opt)
            self.method_combo = W.QComboBox()
            self.method_combo.addItems(["default", "be", "tr", "gear"])
            of.addRow("Integration", self.method_combo)
            self.wgdelay_cb = W.QCheckBox("Waveguide group delay")
            of.addRow(self.wgdelay_cb)
            self.extra_edit = W.QLineEdit()
            self.extra_edit.setPlaceholderText("reltol=1e-5 gmin=1e-11")
            of.addRow("Extra .opt", self.extra_edit)
            self.probe_edit = W.QLineEdit()
            self.probe_edit.setPlaceholderText("V(out),V(in) — blank = all")
            of.addRow("Probe (plot filter)", self.probe_edit)
            self.raw_cb = W.QCheckBox("Raw quadratures (don't derive optical power)")
            of.addRow(self.raw_cb)
            v.addWidget(opt)

            self.run_btn = W.QPushButton("▶  Run")
            self.run_btn.setStyleSheet("font-weight: bold; padding: 8px;")
            self.run_btn.clicked.connect(self._on_run)
            v.addWidget(self.run_btn)
            v.addStretch(1)
            panel.setMaximumWidth(380)
            return panel

        # -- results pane --------------------------------------------------
        def _build_results(self, W, QtCore, Canvas, NavBar, Figure):
            wrap = W.QWidget()
            v = W.QVBoxLayout(wrap)
            self.figure = Figure(figsize=(7, 5), layout="tight")
            self.canvas = Canvas(self.figure)
            v.addWidget(NavBar(self.canvas, wrap))
            rsplit = W.QSplitter(QtCore.Qt.Orientation.Vertical)
            rsplit.addWidget(self.canvas)
            self.log = W.QPlainTextEdit()
            self.log.setReadOnly(True)
            self.log.setMaximumBlockCount(2000)
            self.log.setPlaceholderText("Pipeline log…")
            rsplit.addWidget(self.log)
            rsplit.setStretchFactor(0, 4)
            rsplit.setStretchFactor(1, 1)
            v.addWidget(rsplit, 1)
            return wrap

        # -- behaviour -----------------------------------------------------
        def _set_input(self, path: str):
            self.input_edit.setText(path)
            self.transpile_cb.setChecked(default_transpile_for(path))

        def _on_browse(self):
            path, _ = QtWidgets.QFileDialog.getOpenFileName(
                self, "Open schematic or netlist", "",
                "KiCad / SPICE (*.kicad_sch *.cir *.sp *.net *.xml *.spice);;All files (*)")
            if path:
                self._set_input(path)
                self._sync_enabled()

        def _sync_enabled(self, *_):
            transpiling = self.transpile_cb.isChecked()
            self.analysis_box.setEnabled(transpiling)
            idx = self.analysis_combo.currentIndex()
            self.tran_edit.setEnabled(transpiling and idx == 1)
            self.ac_edit.setEnabled(transpiling and idx == 2)

        def _gather(self) -> RunOptions:
            return RunOptions(
                input_path=self.input_edit.text().strip(),
                transpile=self.transpile_cb.isChecked(),
                analysis=ANALYSES[self.analysis_combo.currentIndex()],
                tran_args=self.tran_edit.text().strip(),
                ac_args=self.ac_edit.text().strip(),
                method=self.method_combo.currentText(),
                waveguide_delay=self.wgdelay_cb.isChecked(),
                probe=self.probe_edit.text().strip(),
                raw_quadratures=self.raw_cb.isChecked(),
                extra_opts=self.extra_edit.text().strip(),
            )

        def _append_log(self, msg: str):
            self.log.appendPlainText(msg)

        def _on_run(self):
            opts = self._gather()
            if not opts.input_path:
                self.statusBar().showMessage("Choose a schematic or netlist first.")
                return
            self.run_btn.setEnabled(False)
            self.log.clear()
            self.statusBar().showMessage("Running…")
            self._opts_in_flight = opts
            self._worker = Worker(opts)
            self._worker.logged.connect(self._append_log)
            self._worker.done.connect(self._on_done)
            self._worker.failed.connect(self._on_failed)
            self._worker.start()

        def _on_done(self, csv: str):
            # Runs on the GUI thread: now safe to touch matplotlib / the canvas.
            try:
                title = Path(self._opts_in_flight.input_path).stem
                kind = render(csv, self.figure, self._opts_in_flight, title)
                self.canvas.draw_idle()
                self._append_log(f"rendered {kind}.")
                self.statusBar().showMessage(f"Done — {kind}.")
            except Exception as exc:
                self._on_failed(f"plotting failed: {exc}")
            finally:
                self.run_btn.setEnabled(True)

        def _on_failed(self, msg: str):
            self._append_log("ERROR: " + msg)
            self.statusBar().showMessage("Failed — see log.")
            self.run_btn.setEnabled(True)

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv[:1])
    win = MainWindow()
    return app, win


# ── entry points ─────────────────────────────────────────────────────────────
def _selftest(initial_input: str | None) -> int:
    """Headless wiring smoke test: build the window under the offscreen Qt
    platform, fire a run, pump the event loop until it finishes, save the
    canvas. Catches signal/slot + thread-marshaling breakage with no display."""
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PyQt6 import QtCore
    sample = initial_input or str(
        HERE.parent / "examples" / "kicad_photonics" / "mrm_single_channel.kicad_sch")
    app, win = make_app_and_window(sample)
    win.show()

    loop = QtCore.QEventLoop()
    result = {"ok": False, "msg": ""}
    # _on_run builds + starts the worker; fairchild's process-spawn latency is
    # orders of magnitude longer than the signal hookup below, so the terminal
    # signals can't be missed in practice (the timer is the belt-and-braces).
    win._on_run()
    if win._worker is None:
        print("selftest: worker never started", file=sys.stderr)
        return 1
    win._worker.done.connect(lambda _c: (result.__setitem__("ok", True), loop.quit()))
    win._worker.failed.connect(lambda m: (result.__setitem__("msg", m), loop.quit()))
    QtCore.QTimer.singleShot(60_000, loop.quit)  # backstop
    loop.exec()

    if not result["ok"]:
        print(f"selftest FAILED: {result['msg'] or 'timed out'}", file=sys.stderr)
        return 1
    out = Path(tempfile.gettempdir()) / "fairchild_gui_selftest.png"
    win.figure.savefig(out, dpi=110)
    print(f"selftest OK → {out}")
    return 0


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser(description="fairchild desktop launcher")
    ap.add_argument("input", nargs="?", help="schematic/netlist to pre-load")
    ap.add_argument("--selftest", action="store_true",
                    help="headless wiring test (offscreen Qt); no display needed")
    args = ap.parse_args()
    if args.selftest:
        return _selftest(args.input)
    app, win = make_app_and_window(args.input)
    win.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
