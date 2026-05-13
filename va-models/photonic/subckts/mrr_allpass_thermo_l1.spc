* mrr_allpass_thermo_l1 — All-pass ring resonator with thermal tuning (Level 1)
*
* Composition: 1× directional_coupler + 1× thermo_phase_shifter_l1
*
* Same feedback topology as mrr_allpass_pn_l1 but uses Joule heating
* to shift the ring resonance.  The heater acts as a resistor; the
* thermal resistance converts dissipated power to a temperature rise
* and hence a phase shift via the thermo-optic effect.
*
* Ports:
*   in_re in_im in_wl    — bus input  (3-wire)
*   out_re out_im out_wl  — through port (3-wire)
*   heat_p heat_n         — heater electrical terminals
*
* Parameters (defaults):
*   kappa_0=0.1        bus–ring coupling fraction
*   L_ring_um=100      ring circumference (µm)
*   n_g=4.2            group index
*   alpha_dB_cm=2.0    ring loss (dB/cm)
*   R_heater=1000      heater resistance (Ω)
*   R_thermal=50000    thermal resistance (K/W)
*   dn_dT=1.86e-4      thermo-optic coefficient (K⁻¹; Si)
*   wavelength_nm=1550 design wavelength (nm)

.subckt mrr_allpass_thermo_l1
+ in_re in_im in_wl
+ out_re out_im out_wl
+ heat_p heat_n
+ kappa_0=0.1 L_ring_um=100 n_g=4.2 alpha_dB_cm=2.0
+ R_heater=1000.0 R_thermal=50000.0 dn_dT=1.86e-4 wavelength_nm=1550.0

* ——— Bus–ring directional coupler ———
Xdc  in_re in_im in_wl
+    ring_out_re ring_out_im in_wl
+    out_re out_im out_wl
+    ring_in_re ring_in_im in_wl
+    directional_coupler kappa_0={kappa_0} wavelength_nm={wavelength_nm}

* ——— Ring thermo-optic phase shifter (full round trip) ———
Xps  ring_in_re ring_in_im in_wl
+    ring_out_re ring_out_im in_wl
+    heat_p heat_n
+    thermo_phase_shifter_l1
+    L_um={L_ring_um} n_g={n_g} alpha_dB_cm={alpha_dB_cm}
+    R_heater={R_heater} R_thermal={R_thermal} dn_dT={dn_dT}
+    wavelength_nm={wavelength_nm}

.ends mrr_allpass_thermo_l1
