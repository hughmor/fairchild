* PN junction phase shifter — DC operating point sweep
*
* Example: measure the optical phase shift vs applied reverse bias
* on a 500 µm PN-doped Si waveguide at 1550 nm.
*
* Circuit:
*   CW laser → PN phase shifter → photodetector
*   (phase shift doesn't change detected power for direct detection,
*    but demonstrates the discipline co-simulation)
*
* Parameters (L1 model):
*   Vpi_L = 2.0 V·cm, L = 500 µm → Vpi ≈ 4.0 V for π phase shift
*   Reverse bias: 0 V → -3 V gives ~3π/4 phase shift

.osdi ../../legacy/va-models/build/cw_laser.osdi
.osdi ../../legacy/va-models/build/pn_phase_shifter_l1.osdi
.osdi ../../legacy/va-models/build/photodetector.osdi

* CW laser: 1 mW at 1550 nm
Xlaser  lre lim wl  cw_laser  power_mW=1.0 wavelength_nm=1550.0

* PN phase shifter L1: 500 µm segment
* (with a shared lambda wire from the laser)
Xpnps   lre lim wl  ore oim wl  vbias 0  pn_phase_shifter_l1 \
        L_um=500.0 n_g=4.2 alpha_dB_cm=3.0 Vpi_L=2.0 V_ref=0.0 \
        wavelength_nm=1550.0

* Photodetector
Xpd     ore oim wl  ph_a 0  photodetector  responsivity=1.0

* Load resistor
Rload   ph_a 0  1k

* PN junction reverse bias (typically 0 to -5 V for Si carrier depletion)
Vbias   vbias 0  DC -2.0

.optical  lre lim wl ore oim

.op
