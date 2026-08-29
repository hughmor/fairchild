"""link.py — drive the CMIM analogue deck from Python.

The deck at `link/netlists/cmim_link.sp` is the physical layer: laser, two TFLN
modulators, the noise arm, the quantum dot amplifier and the receiver, solved in
one Newton loop.  This module runs it.

Why a subprocess and not `import fairchild`.  The compiled Python extension in
the shared virtual environment belongs to whichever checkout last ran
`maturin develop`, and rebuilding it from a worktree would silently replace the
one the main checkout uses.  The command-line binary carries no such coupling.
The simulator is called a few dozen times per characterisation run and never
inside the iteration loop, so the process cost does not matter.

    from link import Link
    lk = Link()
    lk.op()                          # operating point, volts at every probe
    f, h = lk.ac(1e8, 2e11, 40)      # small-signal response, spins to receiver
    t, v = lk.tran(t_s, vx, vw)      # a waveform through the whole chain

Every method takes `**params` that override the deck's `.param` values.
"""
from __future__ import annotations

import csv
import pathlib
import re
import subprocess
import tempfile

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
DECK = HERE.parent / "link" / "netlists" / "cmim_link.sp"
# The worktree's own binary, then the repository's, then whatever is on PATH.
_CANDIDATES = [
    HERE.parents[2] / "target" / "release" / "fairchild",
    HERE.parents[2] / "target" / "debug" / "fairchild",
]


def _binary() -> str:
    for c in _CANDIDATES:
        if c.exists():
            return str(c)
    return "fairchild"


class LinkError(RuntimeError):
    """The simulator refused the deck, or the deck did not converge.

    Raised rather than returning a partial result.  A link that did not solve
    has no operating point, and a caller that gets numbers back cannot tell the
    difference.
    """


class Link:
    def __init__(self, deck: pathlib.Path = DECK, binary: str | None = None):
        self.text = deck.read_text()
        self.dir = deck.parent
        self.binary = binary or _binary()

    # ── deck assembly ────────────────────────────────────────────────────
    def _src(self, params: dict, analysis: str, stim: str | None = None) -> str:
        """Substitute `.param` values, replace the analysis card, drive it.

        `.param` values are substituted at parse time, so they cannot be reached
        with a command-line override.  Rewriting the line is the honest way.
        """
        s = self.text
        for k, v in params.items():
            s, n = re.subn(
                rf"^(\.param\s+{re.escape(k)}\s*=\s*)\S+",
                lambda m: f"{m.group(1)}{float(v):.10g}",
                s,
                flags=re.M,
            )
            if not n:
                raise LinkError(f"no .param named {k!r} in the deck")

        if stim is not None:
            # The two drive lines are single-line by contract; see the deck.
            s = re.sub(r"^Vx\s+xin\s+0\s+.*$", stim.split("\n")[0], s, flags=re.M)
            s = re.sub(r"^Vw\s+win\s+0\s+.*$", stim.split("\n")[1], s, flags=re.M)

        s = re.sub(r"^\.op\s*$", analysis, s, flags=re.M)
        return s

    def _run(self, src: str, probes: list[str], opts: list[str] | None = None):
        with tempfile.TemporaryDirectory() as td:
            # Written into the deck's own directory so that the relative `.va`
            # paths still resolve.
            path = self.dir / f".cmim_run_{pathlib.Path(td).name}.sp"
            path.write_text(src)
            try:
                cmd = [self.binary, "-f", str(path), "--probe", ",".join(probes)]
                cmd += opts or []
                p = subprocess.run(cmd, capture_output=True, text=True)
                if p.returncode != 0:
                    raise LinkError(p.stderr.strip() or "simulator failed")
                rows = list(csv.DictReader(p.stdout.splitlines()))
                if not rows:
                    raise LinkError(f"no output rows; stderr was:\n{p.stderr}")
                return rows
            finally:
                path.unlink(missing_ok=True)

    # ── analyses ─────────────────────────────────────────────────────────
    def op(self, probes: list[str] | None = None, **params) -> dict[str, float]:
        probes = probes or ["v(rxout)"]
        rows = self._run(self._src(params, ".op"), probes)
        return {k: float(v) for k, v in rows[0].items() if k != "analysis"}

    def noise_params(self, **params) -> dict[str, float]:
        """Measure the bias the receiver's noise densities must be evaluated at.

        fairchild evaluates a Verilog-A `white_noise()` power argument at a
        fixed point and never updates it, so a density written against the
        photocurrent is frozen near the dark-current value.  `va_receiver`
        therefore takes the bias as parameters, and this measures them from a
        first pass over the same deck rather than asserting them.

        Returns `{i_ph_n, p_sig_n, rho_ase_n}`, ready to splat into any of the
        analyses below.  Call it at the operating point you care about: a link
        carrying a modulated signal has no single bias, so characterise it at
        each rail, exactly as `.noise` demands.
        """
        op = self.op(["v(soao_re_0)", "v(soao_im_0)", "v(asesoa)"], **params)
        lo = {k.lower(): v for k, v in op.items()}
        p_sig = lo["v(soao_re_0)"] ** 2 + lo["v(soao_im_0)"] ** 2
        rho = max(lo["v(asesoa)"], 0.0)
        # Same optical noise bandwidth the receiver uses; both come from the
        # amplifier's 27.3 nm window at 1310 nm.
        bw_opt = 4.77e12
        i_ph = 0.8 * (p_sig + rho * bw_opt) + 10e-9
        return {"i_ph_n": i_ph, "p_sig_n": p_sig, "rho_ase_n": rho}

    def ac(self, f0: float, f1: float, n_per_decade: int = 20,
           node: str = "rxout", **params):
        """Small-signal response from the spin channel input to `node`.

        Returns `(f, H)` with `H` complex, normalised to nothing: take the
        magnitude at the lowest point yourself if you want a relative curve.

        The linearisation point matters.  With `v_x = 0` the modulator sits at
        quadrature, where the slope of its transfer function is largest and the
        second derivative is zero.  That is the operating point the machine
        actually uses, and it is the only one where a small-signal number means
        what a data sheet means by it.
        """
        src = self._src(params, f".ac dec {int(n_per_decade)} {f0:.6g} {f1:.6g}")
        src = re.sub(r"^Vx\s+xin\s+0\s+DC\s+(\S+)\s*$",
                     r"Vx xin 0 DC \1 AC 1", src, flags=re.M)
        rows = self._run(src, [f"v({node})"])
        # The CSV header preserves the netlist's case, which is not the case
        # the probe was written in.  Match on the lowered name.
        col = {k.lower(): k for k in rows[0]}
        f = np.array([float(r[col["freq_hz"]]) for r in rows])
        mag = np.array([float(r[col[f"mag_v({node})"]]) for r in rows])
        ph = np.array([float(r[col[f"phase_deg_v({node})"]]) for r in rows])
        return f, mag * np.exp(1j * np.deg2rad(ph))

    def tran(self, t: np.ndarray, vx: np.ndarray, vw: np.ndarray,
             tstep: float | None = None, node: str = "rxout",
             method: str = "be", trannoise: bool = False, **params):
        """Push a waveform through the chain and read the receiver.

        `t` is in seconds, `vx` and `vw` are the two arbitrary-waveform-generator
        channels ahead of their drivers.  Returns `(t_out, v_out)` on the
        simulator's own time grid, which is not `t`: resample before comparing.

        Backward Euler by default.  Trapezoidal rings on the steep edges of a
        pulse-shaped waveform, and that ringing is not in the hardware.
        """
        if not (len(t) == len(vx) == len(vw)):
            raise LinkError("t, vx and vw must be the same length")
        tstep = tstep or float(np.min(np.diff(t)))
        stim = (
            "Vx xin 0 PWL(" + " ".join(f"{a:.6e} {b:.6e}" for a, b in zip(t, vx)) + ")\n"
            "Vw win 0 PWL(" + " ".join(f"{a:.6e} {b:.6e}" for a, b in zip(t, vw)) + ")"
        )
        card = f".options method={method}\n.tran {tstep:.6e} {t[-1]:.6e}"
        if trannoise:
            card = ".options trannoise=1\n" + card
        rows = self._run(self._src(params, card, stim=stim), [f"v({node})"])
        col = {k.lower(): k for k in rows[0]}
        tcol = col.get("time_s") or col.get("time")
        return (np.array([float(r[tcol]) for r in rows]),
                np.array([float(r[col[f"v({node})"]]) for r in rows]))


if __name__ == "__main__":
    lk = Link()
    print("operating point:")
    for k, v in lk.op(["v(rxout)", "v(soao_re_0)", "v(soao_im_0)"]).items():
        print(f"  {k:20s} {v: .6e}")
