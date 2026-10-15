#!/usr/bin/env python3
"""Measure laser coherence in a delayed Mach-Zehnder interferometer.

    cargo build --release -p fairchild-cli
    python3 examples/photonic/laser_linewidth.py --selftest
    python3 examples/photonic/laser_linewidth.py --png laser_linewidth.png

Two balanced couplers and a delayed waveguide convert phase diffusion into
complementary output powers. Their mean contrast is exp(-pi linewidth delay).
This uses the CLI so no Python extension rebuild is needed. --fairchild selects
another binary. The delays are integer timesteps to isolate coherence from the
waveguide history interpolator's approximation to a stochastic field.
"""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "fairchild-mplconfig"))

import numpy as np

C0 = 299_792_458.0
STEP = 1e-11
LINEWIDTH = 1e9
STOP = 600e-9
DELAYS = np.array([0, 40e-12, 160e-12, 320e-12, 640e-12])


def deck(delay, linewidth):
    return f"""* delayed laser interferometer
.options waveguide_delay=1 trannoise=1 noiseseed=23
.optical_port src
.optical_port unused
.optical_port short
.optical_port long
.optical_port delayed
.optical_port dark
.optical_port bright
Xlaser src fc_cw_laser power_mW=1 linewidth_hz={linewidth:.12e}
Xsplit src unused short long fc_dcoupler kappa_L=0.7853981633974483
Xdelay long delayed fc_waveguide l_m={C0*delay/4.2:.12e} n_g=4.2 alpha_dB_cm=0 pin_at_ref=1
Xjoin short delayed dark bright fc_dcoupler kappa_L=0.7853981633974483
.tran {STEP:.12e} {STOP:.12e}
"""


def measure(binary, delay, linewidth):
    signals = [f"V({p}_{q}_0)" for p in ("bright", "dark") for q in ("re", "im")]
    with tempfile.TemporaryDirectory(prefix="fairchild-linewidth-") as tmp:
        src = Path(tmp) / "interferometer.sp"
        out = Path(tmp) / "interferometer.csv"
        src.write_text(deck(delay, linewidth))
        result = subprocess.run([str(binary), "-f", str(src), "-o", str(out), "--quiet",
                                 "--probe", ",".join(signals)], capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(result.stderr or result.stdout)
        with out.open() as stream:
            header = stream.readline().strip().split(",")
        data = np.loadtxt(out, delimiter=",", skiprows=1)
    columns = {name: data[:, i] for i, name in enumerate(header)}
    keep = columns["time"] >= max(10 * delay, 1e-9)
    powers = [columns[f"V({p}_re_0)"][keep] ** 2 + columns[f"V({p}_im_0)"][keep] ** 2
              for p in ("bright", "dark")]
    bright, dark = powers
    contrast = np.mean(bright - dark) / np.mean(bright + dark)
    return contrast, np.mean(bright), np.mean(dark)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    executable = "fairchild.exe" if os.name == "nt" else "fairchild"
    ap.add_argument("--fairchild", type=Path,
                    default=Path(__file__).resolve().parents[2] / "target/release" / executable)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--png", type=Path, default=Path("laser_linewidth.png"))
    args = ap.parse_args()
    binary = args.fairchild.resolve()
    if not binary.is_file():
        ap.error("build the CLI first: cargo build --release -p fairchild-cli")
    measured, bright, dark = [], [], []
    for delay in DELAYS:
        visibility, p_bright, p_dark = measure(binary, delay, LINEWIDTH)
        measured.append(visibility)
        bright.append(p_bright)
        dark.append(p_dark)
        expected = np.exp(-np.pi * LINEWIDTH * delay)
        print(f"delay {delay*1e12:6.0f} ps: contrast {visibility:.5f}, analytic {expected:.5f}")
    measured = np.asarray(measured)
    expected = np.exp(-np.pi * LINEWIDTH * DELAYS)
    assert np.max(np.abs(measured - expected)) < 0.04, (measured, expected)
    assert np.max(np.abs(np.asarray(bright) + dark - 1e-3)) < 1e-9
    coherent, _, _ = measure(binary, DELAYS[-1], 0.0)
    assert abs(coherent - 1) < 1e-10, "zero linewidth must retain full contrast"
    if args.selftest:
        print("selftest OK — coherence decays exponentially; the zero-linewidth control stays coherent")
        return
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    tau = np.linspace(0, DELAYS[-1], 300)
    axes[0].plot(tau * 1e12, np.exp(-np.pi * LINEWIDTH * tau), label="analytic coherence")
    axes[0].plot(DELAYS * 1e12, measured, "o", label="simulated interferometer")
    axes[0].axhline(coherent, color="0.6", linestyle="--", label="zero-linewidth control")
    axes[0].set(xlabel="Arm delay (ps)", ylabel="Fringe visibility", ylim=(0, 1.06))
    axes[0].legend(fontsize=8)
    axes[1].plot(DELAYS * 1e12, np.asarray(bright) * 1e3, "o-", label="bright port")
    axes[1].plot(DELAYS * 1e12, np.asarray(dark) * 1e3, "o-", label="dark port")
    axes[1].set(xlabel="Arm delay (ps)", ylabel="Mean optical power (mW)")
    axes[1].legend()
    fig.suptitle("1 GHz laser linewidth: phase diffusion closes the interference fringe")
    fig.savefig(args.png, dpi=180)
    print(f"wrote {args.png}")


if __name__ == "__main__":
    main()
