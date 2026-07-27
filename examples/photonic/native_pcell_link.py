#!/usr/bin/env python3
"""
Hierarchical photonic link from two PCells — `source_bank` + a bank of `mrm`.

Demonstrates fairchild's subcircuit support used the way a PDK would: each
component is a parameterized `.subckt` in its own file, included by the deck,
instantiated with per-instance parameters.

```
 source_bank (8 lasers → 8 ideal MZMs → mux)
        │  one 8-channel optical bus
        ▼
   ring 1 ─ ring 2 ─ … ─ ring 8      each MRM PCell trimmed to its own channel,
        │ (shared thru bus)          each with its own PN bias + heater
        ▼
   thru / drop → photodetectors
```

What this exercises that a flat netlist can't:
  * `.include` of a PCell library, one file per component
  * per-instance `.model` cards — the LEVEL=4 EO model is built from each
    instance's own parameters, so every ring can differ
  * `{…}` parameter arithmetic — the ring's arc length is `{pi*radius}`, so
    radius is the knob a designer actually wants
  * two subcircuits that each take a whole WDM bus — source_bank (8 lasers,
    8 MZMs, a mux) and mrm_wdm8 (one ring, one junction, one heater, all 8
    wavelengths) — matched to the bus by their declared port count

Run:      .venv/bin/python examples/photonic/native_pcell_link.py
Selftest: same, with --selftest (asserts the physics, no plotting)
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np

try:
    import fairchild as fc
except ImportError as e:
    sys.exit(f"fairchild Python package not installed: {e}\n"
             "Build with: maturin develop --release -m crates/fairchild-py/Cargo.toml")

HERE = pathlib.Path(__file__).resolve().parent
PCELLS = HERE / "pcells"

N_CH = 8
# 100 GHz grid, matching source_bank's defaults.
LAMBDAS_NM = [1546.12, 1546.92, 1547.72, 1548.51,
              1549.32, 1550.12, 1550.92, 1551.72]
RADIUS_M = 8e-6
N_G = 4.2
N_EFF_NOM = 2.2810


def _single_ring_deck(n_eff: float, lambda_nm: float) -> str:
    """One MRM alone, one wavelength — the probe used to trim n_eff."""
    return (
        f"* trim probe\n"
        f".include {PCELLS / 'mrm.sp'}\n"
        ".optical_port pin\n.optical_port pth\n"
        ".optical_port pad\n.optical_port pdr\n"
        f"Xl pin fc_cw_laser power_mW=1.0 wavelength_nm={lambda_nm:.6f}\n"
        f"Var pad_re 0 DC 0\nVai pad_im 0 DC 0\n"
        f"Vaw pad_wl 0 DC {lambda_nm * 1e-9:.6e}\n"
        f"Xr pin pth pad pdr vpn 0 hc 0 mrm"
        f" radius={RADIUS_M:.6g} n_eff={n_eff:.9f}\n"
        "Vpn vpn 0 DC 0\nIhc 0 hc DC 0\n.op\n"
    )


def trim_n_eff(lambda_nm: float, n_pts: int = 121) -> float:
    """Effective index that puts a ring's resonance on `lambda_nm`.

    Found numerically rather than from m·λ = n_eff·L, because the segment adds
    first-order dispersion (n_eff walks toward n_g away from wl_ref), so the
    naive comb formula misses. One FSR of n_eff at fixed L is λ/L, so scanning
    that interval is guaranteed to bracket exactly one resonance; a parabolic
    refinement on the drop-port peak then lands it. This mirrors what a real
    ring bank needs anyway — absolute n_eff is a per-ring trim, since fab
    variation moves it far more than a shared model card can carry.
    """
    L = 2.0 * np.pi * RADIUS_M
    span = lambda_nm * 1e-9 / L          # one free spectral range, in n_eff
    grid = N_EFF_NOM + np.linspace(0.0, span, n_pts)
    drop = np.empty(n_pts)
    c = fc.Circuit()
    for i, ne in enumerate(grid):
        c.load_str(_single_ring_deck(float(ne), lambda_nm))
        r = c.run("op")
        drop[i] = (float(r["V(pdr_re_0)"][0]) ** 2
                   + float(r["V(pdr_im_0)"][0]) ** 2)
    i = int(np.argmax(drop))
    if 0 < i < n_pts - 1:                # parabolic sub-step on the peak
        y0, y1, y2 = drop[i - 1], drop[i], drop[i + 1]
        den = y0 - 2 * y1 + y2
        shift = 0.0 if den == 0 else 0.5 * (y0 - y2) / den
        return float(grid[i] + shift * (grid[1] - grid[0]))
    return float(grid[i])


TRIMMED = None


def trims() -> list[float]:
    """n_eff per channel, computed once."""
    global TRIMMED
    if TRIMMED is None:
        TRIMMED = [trim_n_eff(wl) for wl in LAMBDAS_NM]
    return TRIMMED


def deck(*, ring_bias=None, heater_mA=None, drive=None, analysis=".op") -> str:
    """8-channel source bank into a cascade of 8 trimmed MRM PCells."""
    ring_bias = ring_bias if ring_bias is not None else [0.0] * N_CH
    heater_mA = heater_mA if heater_mA is not None else [0.0] * N_CH
    drive = drive if drive is not None else [0.0] * N_CH

    L = [f"* hierarchical PCell link, {N_CH} channels",
         f".include {PCELLS / 'source_bank.sp'}",
         f".include {PCELLS / 'mrm_wdm8.sp'}",
         f".optical_port bus {N_CH}"]
    # One bus segment per ring, plus the add/drop ports of each ring.
    for i in range(N_CH + 1):
        L.append(f".optical_port seg{i} {N_CH}")
    for i in range(N_CH):
        L.append(f".optical_port add{i} {N_CH}")
        L.append(f".optical_port drp{i} {N_CH}")

    # ── stimulus: one source_bank instance carrying the whole bus ───────────
    L.append("Xsrc seg0 " + " ".join(f"d{k + 1}" for k in range(N_CH))
             + " 0 source_bank")
    for k in range(N_CH):
        L.append(f"Vd{k + 1} d{k + 1} 0 DC {drive[k]:.6g}")

    # ── the ring bank: one WDM MRM PCell per ring, each trimmed to its λ ────
    # Separate instances because each ring carries different parameters — that
    # is what a PCell is for. Each instance is the 8-channel cell: one junction
    # and one heater serving the whole bus, which is what a ring on a WDM bus
    # physically is. (The single-channel mrm.sp on an 8-channel bundle would
    # describe eight rings sharing two electrical nodes; the parser rejects it.)
    for k in range(N_CH):
        L.append(
            f"Xr{k} seg{k} seg{k + 1} add{k} drp{k} vpn{k} 0 hc{k} 0 mrm_wdm8"
            f" radius={RADIUS_M:.6g} n_eff={trims()[k]:.9f}"
        )
        L.append(f"Vpn{k} vpn{k} 0 DC {ring_bias[k]:.6g}")
        L.append(f"Ihc{k} 0 hc{k} DC {heater_mA[k] * 1e-3:.6g}")
        # Each ring's add port is dark; its wires still need drivers.
        for w, val in (("re", 0.0), ("im", 0.0)):
            for ch in range(N_CH):
                L.append(f"Va{k}{w}{ch} add{k}_{w}_{ch} 0 DC {val}")
        for ch in range(N_CH):
            L.append(f"Va{k}w{ch} add{k}_wl_{ch} 0 DC {LAMBDAS_NM[ch] * 1e-9:.6e}")

    L.append(analysis)
    return "\n".join(L) + "\n"


def channel_powers(res, port: str, n=N_CH) -> np.ndarray:
    """Per-channel optical power (mW) on a bundle port."""
    return np.array([
        (float(res[f"V({port}_re_{k})"][0]) ** 2
         + float(res[f"V({port}_im_{k})"][0]) ** 2) / 1e-3
        for k in range(n)
    ])


def run(**kw) -> dict:
    c = fc.Circuit()
    c.load_str(deck(**kw))
    r = c.run("op")
    return {
        "out": channel_powers(r, f"seg{N_CH}"),
        "drops": np.array([channel_powers(r, f"drp{k}") for k in range(N_CH)]),
    }


def selftest() -> int:
    # 1. Every ring on resonance with its own channel: each ring should pull its
    #    OWN wavelength into its drop port far more than any other channel.
    base = run()
    drops = base["drops"]
    for k in range(N_CH):
        own = drops[k, k]
        others = np.delete(drops[k], k)
        assert own > 5 * others.max(), (
            f"ring {k} drops {own:.4f} mW of its own channel vs "
            f"{others.max():.4f} mW leakage — trim failed")
    # 2. Through-bus output is therefore suppressed on every channel.
    assert np.all(base["out"] < 0.35), f"thru not suppressed: {base['out']}"

    # 3. The MZM drives gate their channels independently. v_pi = 1 V ⇒ 1 V off.
    off = [0.0] * N_CH
    off[3] = 1.0
    r_off = run(drive=off)
    assert r_off["drops"][3, 3] < 1e-9, \
        f"channel 3 should be dark, got {r_off['drops'][3, 3]:.3e} mW"
    untouched = np.array([r_off["drops"][k, k] for k in range(N_CH) if k != 3])
    ref = np.array([drops[k, k] for k in range(N_CH) if k != 3])
    assert np.allclose(untouched, ref, rtol=1e-8), "driving ch3 disturbed others"

    # 4. Heater on ring 5 detunes ring 5 — per-instance electrical control
    #    through the PCell's own heater terminals. r_heater is per arc and the
    #    arcs are in series, so I²·2R is the whole-ring power against a 26.4 mW
    #    p_pi: 4 mA gives ~5.9 mW, comfortably past the half-FSR point where the
    #    ring has walked off its channel. Do NOT reach for a bigger current —
    #    12 mA is ~2π and brings the ring back ON resonance.
    r_ht = run(heater_mA=[0, 0, 0, 0, 0, 4.0, 0, 0])
    assert r_ht["drops"][5, 5] < 0.2 * drops[5, 5], (
        f"heater did not detune ring 5: {r_ht['drops'][5, 5]:.4f} vs "
        f"{drops[5, 5]:.4f} mW")
    # The ring UPSTREAM of it sees an identical input, so it must not budge at
    # all. The one downstream does shift slightly — a detuned ring 5 passes more
    # of its channel along the bus, which is real cascade crosstalk, not a bug.
    # (relative, not absolute: re-solving the whole circuit moves the answer by
    # Newton's convergence tolerance, ~1e-11 of the value.)
    assert np.isclose(r_ht["drops"][4, 4], drops[4, 4], rtol=1e-9), \
        "upstream ring 4 moved when only ring 5 was heated"
    downstream = abs(r_ht["drops"][6, 6] - drops[6, 6]) / drops[6, 6]
    assert downstream < 0.02, \
        f"downstream ring 6 shifted {downstream:.1%}, expected < 2%"

    print(f"selftest OK — {N_CH} PCell rings trimmed to their own channels "
          f"(drop {drops.diagonal().min():.3f}–{drops.diagonal().max():.3f} mW, "
          f"crosstalk < {max(np.delete(drops[k], k).max() for k in range(N_CH)):.4f} mW); "
          f"MZM gating and per-ring heaters work per instance")
    return 0


def plot(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    base = run()
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))

    im = ax[0].imshow(base["drops"], origin="lower", cmap="magma",
                      aspect="auto")
    ax[0].set_xlabel("wavelength channel")
    ax[0].set_ylabel("ring index")
    ax[0].set_title("Drop power (mW): each ring takes its own channel")
    fig.colorbar(im, ax=ax[0])

    # Bias sweep on one ring: its own channel's drop response.
    v = np.linspace(-2.0, 0.9, 30)
    own, neighbour = [], []
    for vv in v:
        bias = [0.0] * N_CH
        bias[4] = float(vv)
        d = run(ring_bias=bias)["drops"]
        own.append(d[4, 4])
        neighbour.append(d[3, 3])
    ax[1].plot(v, own, "o-", label="ring 4, its channel")
    ax[1].plot(v, neighbour, "s--", label="ring 3 (untouched)")
    ax[1].set_xlabel("PN bias on ring 4 (V)")
    ax[1].set_ylabel("drop power (mW)")
    ax[1].set_title("Per-instance PN bias detunes one ring")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)

    # MZM drive sweep on one channel.
    vd = np.linspace(0, 1.0, 21)
    got = []
    for vv in vd:
        drive = [0.0] * N_CH
        drive[2] = float(vv)
        got.append(run(drive=drive)["drops"][2, 2])
    got = np.array(got)
    ax[2].plot(vd, got / got[0], "o-", label="simulated")
    ax[2].plot(vd, (1 + np.cos(np.pi * vd)) / 2, "k--", lw=0.8,
               label="(1+cos(πV/Vπ))/2")
    ax[2].set_xlabel("MZM drive on channel 2 (V)")
    ax[2].set_ylabel("normalised drop power")
    ax[2].set_title("source_bank MZM transfer")
    ax[2].legend(fontsize=8)
    ax[2].grid(alpha=0.3)

    fig.suptitle("Hierarchical PCell link — source_bank + 8× mrm subcircuits",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    print(f"wrote {out}")


def main() -> int:
    if "--selftest" in sys.argv:
        return selftest()
    selftest()
    plot(HERE / "native_pcell_link.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
