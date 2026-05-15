* MZI push-pull modulator — compound subckt demo
*
* NOTE: This example currently DOES NOT CONVERGE.  See the known-bug note in
* va-models/photonic/subckts/mzi_pn_l1.spc.  Use the monolithic model instead:
*   examples/photonic/sweep_mzi_pn.py  (uses mzi_modulator_pn_l1.osdi directly)
*
* Sweeps the differential PN bias to trace the MZI transfer function.
* Both arms are driven push-pull: V1=+Vbias, V2=−Vbias.
*
* Expected result: P_cross(V) = 0.5·P_in·sin²(π·Vbias·L / Vpi_L)
*   → null at V=0, maximum near V = Vpi_L/(2L) ≈ 2.0 V·cm / (2 × 500e-4 cm) = 20 V
*
* Run:
*   fairchild -f mzi_pn_compound.sp --probe "V(cross_p)"
* Sweep (shell loop):
*   for V in 0 5 10 15 20; do
*     fairchild -f mzi_pn_compound.sp --param "Vbias.value=$V" --probe "V(cross_p)"
*   done

* ——— Load compiled OSDI models ———
.osdi ../../va-models/build/cw_laser.osdi
.osdi ../../va-models/build/directional_coupler.osdi
.osdi ../../va-models/build/pn_phase_shifter_l1.osdi
.osdi ../../va-models/build/photodetector.osdi

* ——— Load compound MZI subckt definition ———
.include ../../va-models/photonic/subckts/mzi_pn_l1.spc

* ——— CW laser at 1550 nm, 1 mW ———
Xlaser  lre lim wl  cw_laser  power_mW=1.0 wavelength_nm=1550.0

* ——— Push-pull MZI modulator ———
* Ports: in  bar  cross  anode1 cathode1  anode2 cathode2
Xmzi  lre lim wl
+     bar_re bar_im wl
+     cross_re cross_im wl
+     va1 0  va2 0
+     mzi_pn_l1
+     kappa_0=0.5 L_arm_um=500 n_g=4.2 alpha_dB_cm=3.0
+     Vpi_L=2.0 V_ref=0.0 wavelength_nm=1550.0

* ——— Photodetectors on bar and cross ports ———
Xpd_bar    bar_re bar_im wl   bar_p 0   photodetector  responsivity=1.0
Xpd_cross  cross_re cross_im wl  cross_p 0  photodetector  responsivity=1.0

Rbar    bar_p 0  1k
Rcross  cross_p 0  1k

* ——— Push-pull drive: arm1 = +V, arm2 = −V ———
Vbias1  va1 0  DC 0.0
Vbias2  va2 0  DC 0.0

* ——— Optical net declarations ———
.optical  lre lim wl  bar_re bar_im  cross_re cross_im

.op
.end
