* Add-drop MRR — compound subckt demo
*
* The ring is built from two directional couplers + a PN phase shifter
* (active section) + a passive waveguide (return path), all connected in a
* feedback loop solved by the NR engine.
*
* Ring parameters:
*   L_active = 50 µm  (PN-modulated section)
*   L_passive = 50 µm  (passive return waveguide)
*   Total L_ring = 100 µm → FSR ≈ 5.72 nm at n_g=4.2
*
* Run:
*   fairchild -f mrr_adddrop_compound.sp --probe "V(thru_p),V(drop_p)"

.osdi ../../legacy/va-models/build/cw_laser.osdi
.osdi ../../legacy/va-models/build/directional_coupler.osdi
.osdi ../../legacy/va-models/build/pn_phase_shifter_l1.osdi
.osdi ../../legacy/va-models/build/waveguide.osdi
.osdi ../../legacy/va-models/build/photodetector.osdi

.include ../../legacy/va-models/photonic/subckts/mrr_adddrop_pn_l1.spc

Xlaser lre lim wl cw_laser power_mW=1.0 wavelength_nm=1544.1

* Add-drop ring: in=laser, add=dark(float), thru+drop=output ports
* Ports: in  thru  drop  add  anode cathode
Xmrr  lre lim wl
+     thru_re thru_im wl
+     drop_re drop_im wl
+     add_re add_im wl
+     vbias 0
+     mrr_adddrop_pn_l1
+     kappa_0=0.1 L_active_um=50 L_passive_um=50
+     n_g=4.2 alpha_dB_cm=2.0
+     Vpi_L=10.0 V_ref=0.0 wavelength_nm=1544.1

Xpd_thru thru_re thru_im wl  thru_p 0  photodetector responsivity=1.0
Xpd_drop drop_re drop_im wl  drop_p 0  photodetector responsivity=1.0
Rthru thru_p 0 1k
Rdrop drop_p 0 1k

Vbias vbias 0 DC 0.0

.optical lre lim wl thru_re thru_im drop_re drop_im add_re add_im

.op
