* Ring resonator wavelength sweep — Phase 2 photonic co-simulation example
*
* Models: cw_laser, directional_coupler, waveguide, photodetector (OSDI v0.4)
* Compile VA models before running:
*   cd va-models && make
*
* Circuit:
*   CW laser → directional coupler (kappa=0.1) → through port → photodetector
*                     ↑                                  ↓
*               ring waveguide  ←←←←←←←←←←←←←←←←←←←←←←
*
* Physics: L_ring=100 µm, n_g=4.2, alpha=2 dB/cm
*   FSR   ≈ 5.72 nm
*   lambda_res ≈ 1544.12 nm, 1549.82 nm, 1555.56 nm  (within sweep range)
*   T_min ≈ 0.916  (8.4% dip, kappa=0.1)
*   FWHM  ≈ 0.098 nm
*
* Run a parametric sweep with the fairchild CLI (Python driver, see below)
* or invoke directly with a single wavelength:
*   fairchild -f ring_resonator_sweep.sp

.osdi ../../va-models/build/cw_laser.osdi
.osdi ../../va-models/build/directional_coupler.osdi
.osdi ../../va-models/build/waveguide.osdi
.osdi ../../va-models/build/photodetector.osdi

Xlaser     laser_re laser_im                                           cw_laser \
           power_mW=1.0 wavelength_nm=1550.0

Xcoupler   laser_re laser_im  ring_fb_re ring_fb_im  \
           through_re through_im  ring_in_re ring_in_im  directional_coupler \
           kappa_0=0.1 wavelength_nm=1550.0

Xring      ring_in_re ring_in_im  ring_fb_re ring_fb_im  waveguide \
           L_um=100.0 n_g=4.2 alpha_dB_cm=2.0 wavelength_nm=1550.0

Xpd        through_re through_im  ph_a 0  photodetector

Rload      ph_a 0  1k

.optical   laser_re laser_im  ring_in_re ring_in_im  ring_fb_re ring_fb_im \
           through_re through_im

.op
.end
