* Waveguide group-delay demonstration.
*
* Shows the `waveguide_delay` SimOption: when enabled, an optical waveguide
* delays its output envelope by the group delay τ_g = L·n_g/c instead of
* transmitting it instantaneously.  This matters whenever the optical
* modulation bandwidth approaches 1/τ_g — e.g. high-speed links, or probing
* the transit time of a long delay line.
*
*   laser ──► MZM (fast electrical drive) ──► long waveguide ──► photodetector
*
* The waveguide is 1 cm long (n_g = 4.2), so τ_g = 1e-2 · 4.2 / 3e8 ≈ 140 ps.
* The MZM is driven by a 20 ps-edge step at t = 500 ps.
*
* Run it both ways and compare V(pd_anode):
*   fairchild -f waveguide_delay_demo.sp                       # instantaneous
*   fairchild -f waveguide_delay_demo.sp --opt waveguide_delay=1   # delayed
*
* Without the option the detector edge tracks the MZM at ~510 ps; with it the
* edge appears ~140 ps later (≈ τ_g), the finite optical transit time.

.optical_port laser_out
.optical_port mzm_out
.optical_port pd_in

* CW pump: 1 mW @ 1550 nm
Xlaser laser_out fc_cw_laser power_mW=1.0 wavelength_nm=1550

* Mach-Zehnder modulator driven by the fast electrical step below
Xmzm laser_out mzm_out vmod 0 fc_mzm

* 1 cm strip waveguide ⇒ τ_g ≈ 140 ps (the delay under test)
Xwg mzm_out pd_in fc_waveguide L_um=10000 n_g=4.2 alpha_dB_cm=0.5 wavelength_nm=1550

* Reverse-biased photodetector + transimpedance load
Xpd pd_in pd_anode 0 fc_photodetector responsivity=0.8 r_shunt=1Meg
Vbias bias 0 DC 1.0
Rload pd_anode bias 1k

* Fast modulation: 0→3 V step at 500 ps, 20 ps edges
Vmod vmod 0 PULSE(0 3 500p 20p 20p 5n 10n)

* Variable-step gear integration resolves the sub-100 ps dynamics.
.options method=gear variable_step=1
.tran 5p 1.5n
