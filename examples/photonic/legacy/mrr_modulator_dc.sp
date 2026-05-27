* MRR modulator — DC operating point
*
* Demonstrates the MRR modulator L1 model:
*   CW laser → all-pass ring resonator modulator → photodetector
*
* The ring is biased at Vbias.  At resonance, the through-port power
* drops; applying reverse bias shifts the resonance.
*
* Ring parameters:
*   L_ring = 100 µm, n_g = 4.2, alpha = 2 dB/cm, kappa = 0.1
*   FSR ≈ 5.72 nm,  λ_res ≈ 1544 nm (nearest to 1550 nm)
*   Vpi_rt = 10 V  (voltage for π round-trip phase shift)
*
* Try: fairchild -f mrr_modulator_dc.sp --probe "V(ph_a)"
* Or with bias sweep:
*   for V in 0 1 2 3 4 5; do
*     fairchild -f mrr_modulator_dc.sp --param "Vbias.value=$V" --probe "V(ph_a)"
*   done

.osdi ../../legacy/va-models/build/cw_laser.osdi
.osdi ../../legacy/va-models/build/mrr_modulator_l1.osdi
.osdi ../../legacy/va-models/build/photodetector.osdi

* CW laser at 1544.12 nm (near ring resonance)
Xlaser  lre lim wl  cw_laser  power_mW=1.0 wavelength_nm=1544.12

* MRR modulator L1
* Ports: in_re in_im in_lambda  out_re out_im out_lambda  anode cathode
Xmod    lre lim wl  ore oim wl  vbias 0  mrr_modulator_l1 \
        kappa_0=0.1 L_ring_um=100.0 n_g=4.2 alpha_dB_cm=2.0 \
        Vpi_rt=10.0 V_ref=0.0 wavelength_nm=1544.12

* Photodetector on the through port
Xpd     ore oim wl  ph_a 0  photodetector  responsivity=1.0

* 1 kΩ transimpedance load
Rload   ph_a 0  1k

* PN junction reverse bias (0 V = unbiased; −V tunes resonance blue)
Vbias   vbias 0  DC 0.0

.optical  lre lim wl ore oim

.op
.end
