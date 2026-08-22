* Native Rust micro-ring modulator (Phase B reference example).
*
* This is what a user building an MRR on top of fairchild's B-phase
* primitives would write today — entirely native devices, no .osdi import,
* no Verilog-A authoring.  Six device classes from `fairchild_core::models::
* photonic` compose into a complete electro-optic transmitter:
*
*   laser ──► wg1 ──► port_a ┌────────────┐ port_c ──► wg2 ──► photodetector
*                            │  DC coupler│
*                  port_b ◄──┤   κL≈0.5   ├──► port_d
*                       │    └────────────┘    │
*                       └─── PN phase shifter ─┘   (the ring loop)
*
* Bundle-port syntax (B2) keeps the netlist readable: every optical port is
* a single named token that the parser expands to its three underlying
* (re, im, λ) wires.
*
* Run:
*   fairchild -f examples/photonic/native_mrr_modulator.sp \
*             --probe "v(pd_anode),v(pd_in_re_0),v(pd_in_im_0)" \
*             --format csv -o /tmp/mrr.csv

.optical_port laser_out
.optical_port wg1_out
.optical_port dc_b
.optical_port dc_c
.optical_port pn_in
.optical_port pd_in

* CW laser: 1 mW @ 1550 nm
Xlaser laser_out fc_cw_laser power_mW=1.0 wavelength_nm=1550

* Coupling waveguide into the ring
Xwg1 laser_out wg1_out fc_waveguide
+    L_um=50 n_g=4.2 alpha_dB_cm=2.0 wavelength_nm=1550

* 2×2 directional coupler tuned for near-critical coupling with the ring
* loss below (κL chosen so t = cos(κL) ≈ ring round-trip amplitude α).
* At t ≈ α the on-resonance transmission |c/a|² → 0.
Xdc wg1_out dc_b dc_c pn_in fc_dcoupler kappa_L=0.336

* The ring loop: PN-junction phase shifter, 500 µm long, with built-in
* propagation loss so the resonance has finite extinction.  Without loss
* the ring is all-pass (|T| = 1 for every φ) — adding 10 dB/cm × 500 µm
* gives α ≈ 0.944, near-critical with t = 0.944 above.
*   V_pi_L = 2 V·mm → 4 V across the PN junction adds π to the
*   round-trip phase, sweeping from on-resonance (V=0, deep notch) to
*   off-resonance (V=Vπ, ~unity transmission).
Xpn pn_in dc_b vmod 0 fc_pn_ps
+   L_um=500 V_pi_L=2e-3 g_pn=1e-3 alpha_dB_cm=10

* Output waveguide carries the through-port field to the detector
Xwg2 dc_c pd_in fc_waveguide
+    L_um=50 n_g=4.2 alpha_dB_cm=2.0 wavelength_nm=1550

* Reverse-biased photodetector + transimpedance load
Xpd pd_in pd_anode 0 fc_photodetector
+   responsivity=0.8 i_dark_a=1e-9 r_shunt=1Meg

Vbias bias 0 DC 1.0
Rload pd_anode bias 1k

* Modulation signal: a single 0→4 V→0 pulse over 2 µs (slow vs the
* photonic settling time — the modulator response is quasi-static).
Vmod vmod 0 PULSE(0 4 100n 100n 100n 800n 2u)

.options method=gear
.tran 5n 2u
