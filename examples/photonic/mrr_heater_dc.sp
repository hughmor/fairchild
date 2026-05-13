* N-doped heater MRR — DC operating point
*
* Demonstrates thermal tuning of a ring resonator using an N-doped
* waveguide heater.  The heater resistance R_heater = 500 Ω; applying
* 3 V gives P = 18 mW, ΔT = 540 K (with R_thermal = 30 kK/W).
*
* Note: ΔT = 540 K is unphysically large — it demonstrates the model
* sensitivity.  In practice, thermal runaway limits usable ΔT to ~100 K.
* For realistic heater voltages (~1 V) with a less thermally isolated ring,
* reduce R_thermal to ~5000 K/W.
*
* Ring parameters (L1 model, steady-state):
*   L_ring = 100 µm, n_g = 4.2, alpha = 10 dB/cm (N-doped FCA)
*   R_heater = 500 Ω, R_thermal = 30 kK/W, dn/dT = 1.86e-4 K⁻¹
*
* At λ = 1544 nm, ΔT = 10 K:
*   Δφ = 2π × 1.86e-4 × 10 × 100e-6 / 1544e-9 = 0.76 rad ≈ λ/8 shift

.osdi ../../va-models/build/cw_laser.osdi
.osdi ../../va-models/build/mrr_heater_l1.osdi
.osdi ../../va-models/build/photodetector.osdi

* CW laser near ring resonance
Xlaser  lre lim wl  cw_laser  power_mW=1.0 wavelength_nm=1544.12

* N-doped heater-integrated MRR L1
Xring   lre lim wl  ore oim wl  vh 0  mrr_heater_l1 \
        kappa_0=0.1 L_ring_um=100.0 n_g=4.2 alpha_dB_cm=10.0 \
        R_heater=500.0 R_thermal=30000.0 dn_dT=1.86e-4 \
        wavelength_nm=1544.12

* Photodetector
Xpd     ore oim wl  ph_a 0  photodetector  responsivity=1.0
Rload   ph_a 0  1k

* Heater bias: 1.0 V → P = 2 mW, ΔT = 60 K
Vheater vh 0  DC 1.0

.optical  lre lim wl ore oim

.op
.end
