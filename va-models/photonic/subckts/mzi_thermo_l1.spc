* mzi_thermo_l1 — Thermo-optic MZI modulator (Level 1)
*
* KNOWN BUG: This compound subckt does NOT converge in fairchild's NR solver.
* Two directional couplers with both outputs of DC1 feeding both inputs of DC2
* causes NR divergence regardless of kappa value (including kappa=0).
* Root cause: competing cross-Jacobian entries from the two OSDI instances on
* the same shared nodes.  Use the monolithic mzi_modulator_thermo_l1.osdi instead.
*
* Composition: 2× directional_coupler + 2× thermo_phase_shifter_l1
*
* Topology:
*   Input DC splits to arm 1 and arm 2 (dark port floats near zero)
*   Independent thermo-optic phase shifters on each arm
*   Output DC recombines arms into bar and cross ports
*
* Single-arm operation is typical: arm 1 heated, arm 2 unbiased.
* Differential operation (one arm hotter, one cooler) doubles efficiency.
*
* Ports:
*   in_re in_im in_wl          — optical input  (3-wire)
*   bar_re bar_im bar_wl       — bar-port output (3-wire)
*   cross_re cross_im cross_wl — cross-port output (3-wire)
*   heat1_p heat1_n            — arm-1 heater electrical terminals
*   heat2_p heat2_n            — arm-2 heater electrical terminals
*
* Parameters (defaults):
*   kappa_0=0.5        50:50 coupler
*   L_arm_um=500       arm (heater) length (µm)
*   n_g=4.2            group index
*   alpha_dB_cm=2.5    arm loss (dB/cm)
*   R_heater=1000      heater resistance (Ω)
*   R_thermal=50000    thermal resistance (K/W)
*   dn_dT=1.86e-4      thermo-optic coefficient (K⁻¹; Si)
*   wavelength_nm=1550 design wavelength (nm)

.subckt mzi_thermo_l1
+ in_re in_im in_wl
+ bar_re bar_im bar_wl
+ cross_re cross_im cross_wl
+ heat1_p heat1_n  heat2_p heat2_n
+ kappa_0=0.5 L_arm_um=500 n_g=4.2 alpha_dB_cm=2.5
+ R_heater=1000.0 R_thermal=50000.0 dn_dT=1.86e-4 wavelength_nm=1550.0

* ——— Input 50:50 directional coupler ———
Xdc1  in_re in_im in_wl
+     dark_re dark_im in_wl
+     arm1_re arm1_im in_wl
+     arm2_re arm2_im in_wl
+     directional_coupler kappa_0={kappa_0} wavelength_nm={wavelength_nm}

* ——— Arm-1 thermo-optic phase shifter ———
Xps1  arm1_re arm1_im in_wl
+     ps1_re ps1_im in_wl
+     heat1_p heat1_n
+     thermo_phase_shifter_l1
+     L_um={L_arm_um} n_g={n_g} alpha_dB_cm={alpha_dB_cm}
+     R_heater={R_heater} R_thermal={R_thermal} dn_dT={dn_dT}
+     wavelength_nm={wavelength_nm}

* ——— Arm-2 thermo-optic phase shifter ———
Xps2  arm2_re arm2_im in_wl
+     ps2_re ps2_im in_wl
+     heat2_p heat2_n
+     thermo_phase_shifter_l1
+     L_um={L_arm_um} n_g={n_g} alpha_dB_cm={alpha_dB_cm}
+     R_heater={R_heater} R_thermal={R_thermal} dn_dT={dn_dT}
+     wavelength_nm={wavelength_nm}

* ——— Output 50:50 directional coupler ———
Xdc2  ps1_re ps1_im in_wl
+     ps2_re ps2_im in_wl
+     bar_re bar_im bar_wl
+     cross_re cross_im cross_wl
+     directional_coupler kappa_0={kappa_0} wavelength_nm={wavelength_nm}

.ends mzi_thermo_l1
