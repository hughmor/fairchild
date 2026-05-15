#!/usr/bin/env python3
"""
Native micro-ring modulator transient — plot the optical-power modulation
at the photodetector as the PN-junction phase-shifter voltage swings
through 0 → Vπ → 0.

Uses only B-phase native Rust devices via the fairchild Python bindings:
no .osdi import, no Verilog-A.  Run after `maturin develop --release`.

Topology:

    laser ──► wg1 ──► port_a ┌────────────┐ port_c ──► wg2 ──► PD
                             │  DC coupler│
                  port_b ◄───┤   κL≈0.336 ├──► port_d
                       │     └────────────┘    │
                       └───── PN phase shifter ─┘ (the ring)

Companion netlist: native_mrr_modulator.sp  (drives the same circuit
through the CLI; this Python script just wraps it with plotting).
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
NETLIST_PATH = HERE / "native_mrr_modulator.sp"


def main():
    c = fairchild.Circuit()
    c.load(str(NETLIST_PATH))
    # Use GEAR (configured in the .sp via `.options method=gear`) — robust
    # against the sharp PWL slopes of the modulation signal.
    result = c.run("tran", step=5e-9, stop=2e-6, variable_step=True)

    t        = result.time()
    v_pn     = result["V(vmod)"]
    v_pd     = result["V(pd_anode)"]
    pd_re    = result["V(pd_in_re_0)"]
    pd_im    = result["V(pd_in_im_0)"]
    pd_power = pd_re ** 2 + pd_im ** 2  # |A|² in W (V/m units → W with our unit choice)

    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 7))

    axes[0].plot(t * 1e9, v_pn, label="V(pn)")
    axes[0].set_ylabel("PN voltage (V)")
    axes[0].set_title("Native MRR modulator transient")
    axes[0].grid(alpha=0.3)
    axes[0].legend(loc="upper right")

    axes[1].plot(t * 1e9, pd_power * 1e3, color="tab:orange",
                 label="|A|² at PD input")
    axes[1].set_ylabel("Optical power (mW)")
    axes[1].grid(alpha=0.3)
    axes[1].legend(loc="upper right")

    axes[2].plot(t * 1e9, v_pd, color="tab:green", label="V(pd_anode)")
    axes[2].set_ylabel("PD anode (V)")
    axes[2].set_xlabel("time (ns)")
    axes[2].grid(alpha=0.3)
    axes[2].legend(loc="upper right")

    plt.tight_layout()
    out_path = HERE / "native_mrr_modulator.png"
    fig.savefig(out_path, dpi=120)
    print(f"wrote {out_path}")
    plt.show()


if __name__ == "__main__":
    main()
