* mrr_allpass_pn_l1 — All-pass ring resonator modulator (PN, Level 1)
*
* Composition: 1× directional_coupler + 1× pn_phase_shifter_l1
*
* Topology (feedback loop):
*
*          ┌─── pn_phase_shifter ───┐
*          │    (full ring, L_ring) │
*     in ──┤DC a1              a2 ──┘  (ring_out feeds back to DC.a2)
*          │   b1 ──→ out (through)
*          └── b2 ──→ ring_in ──→ PS ──→ ring_out
*
* The NR solver closes the feedback loop algebraically; because the optical
* physics is linear in the field amplitudes, it converges in one Newton step.
*
* Resonance condition: phi_rt = 2·π·n_g·L_ring / λ + Δφ_PN = 2πm
*   → through-port power minimum at resonance
*
* Ports:
*   in_re in_im in_wl   — bus input  (3-wire)
*   out_re out_im out_wl — bus output / through port (3-wire)
*   anode cathode        — PN junction (electrical)
*
* Parameters (defaults match mrr_modulator_l1):
*   kappa_0=0.1        bus–ring power coupling fraction
*   L_ring_um=100      ring circumference (µm)
*   n_g=4.2            group index
*   alpha_dB_cm=2.0    ring loss (dB/cm)
*   Vpi_L=2.0          PN Vπ·L product (V·cm)
*   V_ref=0.0          zero-phase bias point (V)
*   wavelength_nm=1550 design wavelength (nm)

.subckt mrr_allpass_pn_l1
+ in_re in_im in_wl
+ out_re out_im out_wl
+ anode cathode
+ kappa_0=0.1 L_ring_um=100 n_g=4.2 alpha_dB_cm=2.0
+ Vpi_L=2.0 V_ref=0.0 wavelength_nm=1550.0

* ——— Bus–ring directional coupler ———
* a1=(bus input), a2=(ring output, feedback), b1=(through), b2=(ring input)
Xdc  in_re in_im in_wl
+    ring_out_re ring_out_im in_wl
+    out_re out_im out_wl
+    ring_in_re ring_in_im in_wl
+    directional_coupler kappa_0={kappa_0} wavelength_nm={wavelength_nm}

* ——— Ring phase shifter (full round trip) ———
* Feedback: input comes from DC.b2, output drives DC.a2
Xps  ring_in_re ring_in_im in_wl
+    ring_out_re ring_out_im in_wl
+    anode cathode
+    pn_phase_shifter_l1
+    L_um={L_ring_um} n_g={n_g} alpha_dB_cm={alpha_dB_cm}
+    Vpi_L={Vpi_L} V_ref={V_ref} wavelength_nm={wavelength_nm}

.ends mrr_allpass_pn_l1
