#!/usr/bin/env python3
"""
8x8 cyclic AWG router — routing map, crosstalk matrix, and passband sweep.

Three figures from one device:

  1. Routing map.  Input port i channel k leaves on output port (i+k) mod 8.
     Every output receives exactly one wavelength from every input, which is
     what makes an AWGR an all-to-all interconnect with no switching.

  2. Crosstalk matrix.  Power at each output port with a single input lit,
     showing the intended channel against the -30/-40 dB floors that array
     phase errors impose on a real device.  Note where the leakage lands: in
     the *slots the intended channel does not occupy*, never summed into it.

  3. Passband sweep.  One laser tuned across a channel, showing insertion
     loss, the -3 dB width, and the penalty a detuned laser pays.

Run:  python3 examples/photonic/native_awgr_router.py [--png out.png]

These decks set no solver options.  They used to need a tight `vntol`, because
the lambda wires carry ~1.55e-6 and a *voltage* tolerance of 1e-6 let Newton
stop with lambda ~10 pm out -- a real detuning for a 40 GHz passband.  Lambda
rows now carry their own `lambdatol` in the solver, so the deck no longer has to
know about it.  See docs/photonic_awgr.md, "Solver tolerance".
"""
import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np

try:
    import fairchild
except ImportError as e:
    sys.exit(f"fairchild Python package not installed: {e}\n"
             "Build with: cd crates/fairchild-py && maturin develop --release")

C0 = 299_792_458.0
N = 8
LAMBDA0_NM = 1550.0
DF_GHZ = 100.0


def grid_nm(n=N, lambda0_nm=LAMBDA0_NM, df_ghz=DF_GHZ):
    """The device's channel grid, in nm. Uniform in frequency, not wavelength."""
    f0 = C0 / (lambda0_nm * 1e-9)
    return np.array([C0 / (f0 + k * df_ghz * 1e9) * 1e9 for k in range(n)])


def deck(lit_ports, lam_nm, extra="", n=N):
    """An n x n router with `lit_ports` driven at unit amplitude on every channel.

    `lam_nm[k]` is the wavelength put on channel k of every lit port.
    """
    lines = ["* AWGR"]
    for i in lit_ports:
        for k in range(n):
            lines += [f"V{i}_{k}r in{i}_{k}_re 0 DC 1.0",
                      f"V{i}_{k}i in{i}_{k}_im 0 DC 0.0",
                      f"V{i}_{k}w in{i}_{k}_wl 0 DC {lam_nm[k] * 1e-9:.12e}"]
    wires = "".join(
        f" {side}{p}_{k}_re {side}{p}_{k}_im {side}{p}_{k}_wl"
        for side in ("in", "out") for p in range(n) for k in range(n)
    )
    # No .options: λ wires carry their own `lambdatol` in the solver now.
    lines += [f"Xr{wires} fc_awgr{extra}", ".op", ".end"]
    return "\n".join(lines) + "\n"


def solve(text):
    c = fairchild.Circuit()
    c.load_str(text)
    return c.run("op")


def power(res, node_re, node_im):
    return res[f"V({node_re})"][0] ** 2 + res[f"V({node_im})"][0] ** 2


def routing_map(lam):
    """Ideal mode: which input reaches which output, per channel."""
    res = solve(deck(range(N), lam))
    m = np.zeros((N, N))  # [output port, channel] -> power
    for j in range(N):
        for k in range(N):
            m[j, k] = power(res, f"out{j}_{k}_re", f"out{j}_{k}_im")
    return m


def crosstalk_matrix(lam, lit=0):
    """Realistic mode, one input lit: power at every (output port, channel)."""
    res = solve(deck([lit], lam,
                     extra=f" df_ghz={DF_GHZ} fwhm_ghz=40 il_db=3"
                           " xt_adj_db=-30 xt_bg_db=-40"))
    m = np.zeros((N, N))
    for j in range(N):
        for k in range(N):
            m[j, k] = power(res, f"out{j}_{k}_re", f"out{j}_{k}_im")
    return m


def passband_sweep(lam, npts=161):
    """Tune channel 0's laser across its passband; read the routed output."""
    f0 = C0 / (LAMBDA0_NM * 1e-9)
    offsets = np.linspace(-1.5 * DF_GHZ, 1.5 * DF_GHZ, npts)  # GHz
    out = {}
    for shape_p, label in [(1.0, "Gaussian (p=1)"), (3.0, "flat-top (p=3)")]:
        p = []
        for d in offsets:
            probe = lam.copy()
            probe[0] = C0 / (f0 + d * 1e9) * 1e9
            res = solve(deck([0], probe,
                             extra=f" df_ghz={DF_GHZ} fwhm_ghz=40 il_db=3"
                                   f" shape_p={shape_p}"
                                   " xt_adj_db=-30 xt_bg_db=-40"))
            # Channel 0 of input 0 routes to output 0.
            p.append(power(res, "out0_0_re", "out0_0_im"))
        out[label] = np.array(p)
    return offsets, out


def fwhm_ghz(offsets, p):
    """Half-power width, with the crossings interpolated off the sweep grid.

    Reading the width off the sample points alone quantises it to the sweep
    step (1.9 GHz here), which is coarser than the number being checked.
    """
    half = p.max() / 2.0
    above = np.where(p >= half)[0]
    lo, hi = above[0], above[-1]
    left = np.interp(half, [p[lo - 1], p[lo]], [offsets[lo - 1], offsets[lo]])
    right = np.interp(half, [p[hi + 1], p[hi]], [offsets[hi + 1], offsets[hi]])
    return right - left


def db(x, floor=-70.0):
    with np.errstate(divide="ignore"):
        return np.maximum(10.0 * np.log10(np.maximum(x, 1e-30)), floor)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--png", nargs="?", const="native_awgr_router.png")
    ap.add_argument("--selftest", action="store_true",
                    help="assert the routing and crosstalk numbers, no plot")
    args = ap.parse_args()

    lam = grid_nm()
    routes = routing_map(lam)
    xt = crosstalk_matrix(lam)
    offsets, sweep = passband_sweep(lam)

    if args.selftest:
        # Ideal mode is a lossless permutation: exactly one unit-power entry
        # per (output, channel), sourced from input (j - k) mod N.
        assert np.allclose(routes, 1.0), routes
        # Realistic mode, input 0 lit: output j carries channel j at the full
        # -3 dB, everything else on a floor.
        il = 10 ** (-3.0 / 10.0)
        for j in range(N):
            assert abs(xt[j, j] / il - 1.0) < 1e-6, (j, xt[j, j])
            adj = (j + 1) % N
            want = il * 10 ** (-30.0 / 10.0)
            assert abs(xt[j, adj] / want - 1.0) < 1e-6, (j, xt[j, adj], want)
        # The passband is 40 GHz wide at -3 dB below its own peak.
        g = sweep["Gaussian (p=1)"]
        width = fwhm_ghz(offsets, g)
        assert abs(width - 40.0) < 0.2, width
        # A flat-top passes more at the band edge and less far out.
        f = sweep["flat-top (p=3)"]
        near = np.argmin(np.abs(offsets - 15.0))
        far = np.argmin(np.abs(offsets - 45.0))
        assert f[near] > g[near] and f[far] < g[far]
        print("selftest OK: routing, crosstalk floors, 40 GHz FWHM, flat-top shape")
        return

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

    ax = axes[0]
    src = np.array([[(j - k) % N for k in range(N)] for j in range(N)])
    im = ax.imshow(src, cmap="twilight", origin="lower")
    for j in range(N):
        for k in range(N):
            ax.text(k, j, str(src[j, k]), ha="center", va="center", fontsize=8)
    ax.set(xlabel="channel k", ylabel="output port j",
           title="Ideal routing: source input port\n$i = (j-k)\\ \\mathrm{mod}\\ N$")
    fig.colorbar(im, ax=ax, label="input port i")

    ax = axes[1]
    im = ax.imshow(db(xt), cmap="magma", origin="lower", vmin=-45, vmax=0)
    ax.set(xlabel="channel k", ylabel="output port j",
           title="Input 0 lit: output power (dB)\n"
                 "3 dB IL, $-30$/$-40$ dB floors")
    fig.colorbar(im, ax=ax, label="dB")

    ax = axes[2]
    for label, p in sweep.items():
        ax.plot(offsets, db(p), label=label)
    ax.axhline(-3.0, ls=":", c="0.5", lw=1)
    ax.axvline(-20.0, ls=":", c="0.5", lw=1)
    ax.axvline(20.0, ls=":", c="0.5", lw=1)
    ax.set(xlabel="laser detuning from channel 0 centre (GHz)",
           ylabel="routed output power (dB)",
           title="Passband: 40 GHz FWHM, 3 dB IL", ylim=(-45, 1))
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    if args.png:
        fig.savefig(args.png, dpi=140)
        print(f"wrote {args.png}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
