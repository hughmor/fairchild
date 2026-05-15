#!/usr/bin/env python3
"""
WDM micro-ring modulator transient — two wavelength channels, one ring.

Companion to native_wdm_mrr_modulator.sp.  The two lasers are placed
symmetrically (±50 pm) around the ring's V_pn = 0 resonance.  As the
PN-junction voltage ramps 0 → 4 V → 0 the ring resonance walks ~570 pm
red across the spectrum, producing two very different transmission
traces at the two photodetectors:

  - Channel 0 (laser at 1549.95 nm, blue-side):
       resonance moves AWAY → monotonic increase to high transmission.
  - Channel 1 (laser at 1550.05 nm, red-side):
       resonance walks RIGHT ACROSS the wavelength → sharp dip event
       when ring λ_res matches the laser, then recovery.

Same modulator, same drive, different responses.  No MUX/DEMUX device
needed — the bundle-port abstraction makes each wavelength its own
independent SVEA channel.
"""
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np

try:
    import fairchild
except ImportError as e:
    sys.exit(f"fairchild Python package not installed: {e}\n"
             "Build with: cd crates/fairchild-py && maturin develop --release")

HERE = pathlib.Path(__file__).resolve().parent
NETLIST = HERE / "native_wdm_mrr_modulator.sp"


def main():
    c = fairchild.Circuit()
    c.load(str(NETLIST))
    result = c.run("tran", step=5e-9, stop=2e-6, variable_step=True)

    t      = result.time()
    v_mod  = result["V(vmod)"]
    v_pd1  = result["V(pd1_anode)"]
    v_pd2  = result["V(pd2_anode)"]
    # Optical power at each photodetector input (channel-by-channel).
    p1     = result["V(pd_in_re_0)"] ** 2 + result["V(pd_in_im_0)"] ** 2
    p2     = result["V(pd_in_re_1)"] ** 2 + result["V(pd_in_im_1)"] ** 2

    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8.5, 7.5))

    axes[0].plot(t * 1e9, v_mod, color="tab:gray", label="V(pn)")
    axes[0].set_ylabel("PN voltage (V)")
    axes[0].set_title("WDM micro-ring modulator — two channels, one ring")
    axes[0].grid(alpha=0.3)
    axes[0].legend(loc="upper right")

    axes[1].plot(t * 1e9, p1 * 1e3, color="tab:blue",
                 label="ch0 — λ=1549.95 nm (blue-side)")
    axes[1].plot(t * 1e9, p2 * 1e3, color="tab:orange",
                 label="ch1 — λ=1550.05 nm (red-side)")
    axes[1].set_ylabel("Optical power at PD (mW)")
    axes[1].grid(alpha=0.3)
    axes[1].legend(loc="lower right")

    axes[2].plot(t * 1e9, v_pd1, color="tab:blue",  label="V(pd1)")
    axes[2].plot(t * 1e9, v_pd2, color="tab:orange", label="V(pd2)")
    axes[2].set_xlabel("time (ns)")
    axes[2].set_ylabel("PD anode voltage (V)")
    axes[2].grid(alpha=0.3)
    axes[2].legend(loc="lower right")

    plt.tight_layout()
    out = HERE / "native_wdm_mrr_modulator.png"
    fig.savefig(out, dpi=120)
    print(f"wrote {out}")
    plt.show()


if __name__ == "__main__":
    main()
