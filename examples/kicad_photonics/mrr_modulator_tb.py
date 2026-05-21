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

#%%
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
NETLIST_PATH = HERE / "wdm_mrm_transpiled_netlist.sp"

#%%

c = fairchild.Circuit()
c.load(str(NETLIST_PATH))

#%%

result = c.run("tran", step=1e-9, stop=2e-6, variable_step=True, method="gear")

t = result.time()
v_pn = result["V(/Vmod)"]
v_pd = result["V(V_THRU)"]

pd_re = result["V(THRU_OPT_re_0)"]
pd_im = result["V(THRU_OPT_im_0)"]
pd_power = pd_re ** 2 + pd_im ** 2  # |A|² in W (V/m units → W with our unit choice)

mrr_re = result["V(Net-_CPL1-b1__re_0)"]
mrr_im = result["V(Net-_CPL1-b1__im_0)"]
mrr_power = mrr_re ** 2 + mrr_im ** 2


#%%

fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 7))

axes[0].plot(t * 1e9, v_pn, label="V(pn)")
axes[0].set_ylabel("PN voltage (V)")
axes[0].set_title("Native MRR modulator transient")
axes[0].grid(alpha=0.3)
axes[0].legend(loc="upper right")

axes[1].plot(t * 1e9, pd_power * 1e3, color="tab:orange",
                label="|A|² at PD input")
# axes[1].plot(t * 1e9, mrr_power * 1e3, color="tab:red",
#                 label="|A|² in MRR")
axes[1].set_ylabel("Optical power (mW)")
axes[1].grid(alpha=0.3)
axes[1].legend(loc="upper right")

axes[2].plot(t * 1e9, v_pd, color="tab:green", label="V(pd_anode)")
axes[2].set_ylabel("PD anode (V)")
axes[2].set_xlabel("time (ns)")
axes[2].grid(alpha=0.3)
axes[2].legend(loc="upper right")

plt.tight_layout()


# %%

results = c.sweep("XCWL1.wavelength_nm", np.linspace(1549, 1551, 5),
                  "tran", step=5e-9, stop=2e-6, variable_step=True)

fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 7))

for result in results:
    t = result.time()
    v_pn = result["V(/Vmod)"]
    v_pd = result["V(V_THRU)"]

    pd_re = result["V(THRU_OPT_re_0)"]
    pd_im = result["V(THRU_OPT_im_0)"]
    pd_power = pd_re ** 2 + pd_im ** 2  # |A|² in W (V/m units → W with our unit choice)

    mrr_re = result["V(Net-_CPL1-b1__re_0)"]
    mrr_im = result["V(Net-_CPL1-b1__im_0)"]
    mrr_power = mrr_re ** 2 + mrr_im ** 2

    axes[0].plot(t * 1e9, v_pn, label="V(pn)")

    axes[1].plot(t * 1e9, pd_power * 1e3,
                    label="|A|² at PD input")
    # axes[1].plot(t * 1e9, mrr_power * 1e3, color="tab:red",
    #                 label="|A|² in MRR")

    axes[2].plot(t * 1e9, v_pd, label="V(pd_anode)")


axes[0].set_ylabel("PN voltage (V)")
axes[0].set_title("Native MRR modulator transient")
axes[0].grid(alpha=0.3)
axes[0].legend(loc="upper right")
axes[1].set_ylabel("Optical power (mW)")
axes[1].grid(alpha=0.3)
axes[1].legend(loc="upper right")
axes[2].set_ylabel("PD anode (V)")
axes[2].set_xlabel("time (ns)")
axes[2].grid(alpha=0.3)
axes[2].legend(loc="upper right")
plt.tight_layout()

# %%


wavelenghts = np.linspace(1549, 1551, 1000)
results = c.sweep("XCWL1.wavelength_nm", wavelenghts, "op")

inp = np.zeros_like(wavelenghts)
thru = np.zeros_like(wavelenghts)
drop = np.zeros_like(wavelenghts)
add = np.zeros_like(wavelenghts)

fig, ax = plt.subplots(1, 1, sharex=True, figsize=(8, 7))

for i, (wl, result) in enumerate(zip(wavelenghts, results)):
    v_pn = result["V(/Vmod)"]
    v_pd = result["V(V_THRU)"]

    in_re = result["V(IN_OPT_re_0)"]
    in_im = result["V(IN_OPT_im_0)"]
    in_power = in_re ** 2 + in_im ** 2  # |A|² in W (V/m units → W with our unit choice)

    thru_re = result["V(THRU_OPT_re_0)"]
    thru_im = result["V(THRU_OPT_im_0)"]
    thru_power = thru_re ** 2 + thru_im ** 2  # |A|² in W (V/m units → W with our unit choice)
    
    # add_re = result["V(ADD_OPT_re_0)"]
    # add_im = result["V(ADD_OPT_im_0)"]
    # add_power = add_re ** 2 + add_im ** 2  # |A|² in W (V/m units → W with our unit choice)
    
    drop_re = result["V(DROP_OPT_re_0)"]
    drop_im = result["V(DROP_OPT_im_0)"]
    drop_power = drop_re ** 2 + drop_im ** 2  # |A|² in W (V/m units → W with our unit choice)
    
    thru[i] = thru_power[0]
    drop[i] = drop_power[0]
    # add[i] = add_power[0]
    inp[i] = in_power[0]


ax.plot(wavelenghts, thru/inp, label='THRU')
ax.plot(wavelenghts, drop/inp, label='DROP')
# ax.plot(wavelenghts, add/inp, label='ADD')
ax.set_ylabel("Transmission")
ax.set_xlabel("Wavelength (nm)")
ax.grid(alpha=0.3)
ax.legend(loc="upper right")
plt.tight_layout()


# %%


### GIONA

GIONA_NETLIST_PATH = HERE / "giona_ff_netlist.sp"
giona = fairchild.Circuit()
giona.load(str(GIONA_NETLIST_PATH))
giona_results = giona.run(
    # "tran", step=1e-12, stop=1e-6,
    "op",
    method="gear",
    reltol=1e-2,
    abstol=1e-10,
    variable_step=True,
    itl1=200,
    itl4=100,
    max_rejections=50,
    srcsteps=21,
    solver="sparse",
    verbose=True,
)




# %%
