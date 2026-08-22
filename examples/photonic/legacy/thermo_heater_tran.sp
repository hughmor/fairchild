* Thermo-optic phase shifter L2 — transient heating step response
*
* Demonstrates the L2 thermo_phase_shifter_l2 model with an external
* thermal RC network.  At t=0, the heater is stepped on; the optical
* phase shifts with the thermal time constant τ = R_thermal × C_thermal.
*
* Thermal parameters:
*   R_heater  = 1000 Ω  (TiN heater)
*   R_thermal = 50 kΩ → 50 kK/W thermal resistance
*   C_thermal = 1 µF  → 1 µJ/K heat capacity
*   τ_thermal = R_thermal × C_thermal = 50 ms
*   At Vheat = 2 V: P = 4 mW, ΔT_ss = 200 K
*
* Run: fairchild -f thermo_heater_tran.sp --probe "V(ph_a),V(T_dev)"
* Plot the optical output V(ph_a) and temperature V(T_dev) over time.

.osdi ../../legacy/va-models/build/cw_laser.osdi
.osdi ../../legacy/va-models/build/thermo_phase_shifter_l2.osdi
.osdi ../../legacy/va-models/build/photodetector.osdi

* CW laser
Xlaser  lre lim wl  cw_laser  power_mW=1.0 wavelength_nm=1550.0

* Thermo-optic phase shifter L2 with external T_node
Xheater  lre lim wl  ore oim wl  heat_p 0  T_dev  thermo_phase_shifter_l2 \
         L_um=200.0 n_g=4.2 alpha_dB_cm=2.5 \
         R_heater=1000.0 dn_dT=1.86e-4 wavelength_nm=1550.0

* External thermal circuit: parallel Rth and Cth to ground.
* V(T_dev) = ΔT (K).  τ = Rth × Cth = 50 ms.
Rth  T_dev  0      50000   ; 50 kΩ → 50 kK/W thermal resistance to ambient
Cth  T_dev  0      1e-6    ; 1 µF → 1 µJ/K heat capacity

* Photodetector
Xpd     ore oim wl  ph_a 0  photodetector  responsivity=1.0
Rload   ph_a 0  1k

* Heater: step from 0 to 2 V at t = 10 ms
Vheater heat_p 0  PULSE(0 2 10m 1u 1u 1 1)

.optical  lre lim wl ore oim

* Simulate for 200 ms; step = 1 ms
.tran 1m 200m
