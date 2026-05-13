* All-pass MRR modulator — compound subckt demo
*
* The ring resonator is built from a directional coupler and a PN phase
* shifter connected in a feedback loop.  The NR solver closes the loop
* algebraically — identical physics to the monolithic mrr_modulator_l1.
*
* Ring parameters:
*   L_ring = 100 µm, n_g = 4.2 → FSR ≈ 5.72 nm
*   Nearest resonance to 1550 nm: ≈ 1544.1 nm (same as mrr_modulator_dc.sp)
*
* Run:
*   fairchild -f mrr_compound.sp --probe "V(ph_a)"
* Sweep:
*   for V in 0 1 2 3 4 5; do
*     fairchild -f mrr_compound.sp --param "Vbias.value=$V" --probe "V(ph_a)"
*   done

* ——— Load compiled OSDI models ———
.osdi ../../va-models/build/cw_laser.osdi
.osdi ../../va-models/build/directional_coupler.osdi
.osdi ../../va-models/build/pn_phase_shifter_l1.osdi
.osdi ../../va-models/build/photodetector.osdi

* ——— Load compound MRR subckt definition ———
.include ../../va-models/photonic/subckts/mrr_allpass_pn_l1.spc

* ——— CW laser at 1544.1 nm (near ring resonance) ———
Xlaser  lre lim wl  cw_laser  power_mW=1.0 wavelength_nm=1544.1

* ——— All-pass ring resonator modulator (compound) ———
* Ports: in  out  anode cathode
Xmrr  lre lim wl
+     ore oim wl
+     vbias 0
+     mrr_allpass_pn_l1
+     kappa_0=0.1 L_ring_um=100.0 n_g=4.2 alpha_dB_cm=2.0
+     Vpi_L=10.0 V_ref=0.0 wavelength_nm=1544.1

* ——— Photodetector on through port ———
Xpd  ore oim wl  ph_a 0  photodetector  responsivity=1.0
Rload  ph_a 0  1k

* ——— PN reverse bias ———
Vbias  vbias 0  DC 0.0

* ——— Optical net declarations ———
.optical  lre lim wl  ore oim

.op
.end
