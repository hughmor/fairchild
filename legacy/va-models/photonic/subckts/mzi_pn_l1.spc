* mzi_pn_l1 — Push-pull MZI modulator (PN junctions, Level 1)
*
* KNOWN BUG: This compound subckt does NOT converge in fairchild's NR solver.
* Two directional couplers with both outputs of DC1 feeding both inputs of DC2
* causes NR divergence regardless of kappa value (including kappa=0).
* Root cause: competing cross-Jacobian entries from the two OSDI instances on
* the same shared nodes.  Use the monolithic mzi_modulator_pn_l1.osdi instead.
*
* Composition: 2× directional_coupler + 2× pn_phase_shifter_l1
*
* Topology:
*   Input DC splits to arm 1 and arm 2 (dark port floats near zero)
*   Independent PN phase shifters on each arm
*   Output DC recombines arms into bar and cross ports
*
* Signal flow (push-pull: V1 = +V, V2 = −V):
*   phi1 = phi_base + π·V·L / Vpi_L
*   phi2 = phi_base − π·V·L / Vpi_L
*   P_bar   ∝ cos²(phi1 − phi2) = cos²(2·π·V·L / Vpi_L)
*   P_cross ∝ sin²(phi1 − phi2)
*
* Ports:
*   in_re in_im in_wl          — optical input  (3-wire)
*   bar_re bar_im bar_wl       — bar-port output (3-wire)
*   cross_re cross_im cross_wl — cross-port output (3-wire)
*   anode1 cathode1            — arm-1 PN junction (electrical)
*   anode2 cathode2            — arm-2 PN junction (electrical)
*
* Parameters (defaults):
*   kappa_0=0.5        50:50 coupler (nominal)
*   L_arm_um=500       arm length (µm)
*   n_g=4.2            group index
*   alpha_dB_cm=3.0    arm loss (dB/cm)
*   Vpi_L=2.0          PN Vπ·L product (V·cm)
*   V_ref=0.0          zero-phase bias point (V)
*   wavelength_nm=1550 design wavelength (nm)
*
* Usage: place this file alongside your top-level netlist and add
*   .include "mzi_pn_l1.spc"
* then instantiate with:
*   Xmzi  in_re in_im wl  bar_re bar_im wl  cross_re cross_im wl
*   +      va1 vc1  va2 vc2
*   +      mzi_pn_l1  L_arm_um=500 Vpi_L=2.0 wavelength_nm=1550

.subckt mzi_pn_l1
+ in_re in_im in_wl
+ bar_re bar_im bar_wl
+ cross_re cross_im cross_wl
+ anode1 cathode1  anode2 cathode2
+ kappa_0=0.5 L_arm_um=500 n_g=4.2 alpha_dB_cm=3.0
+ Vpi_L=2.0 V_ref=0.0 wavelength_nm=1550.0

* ——— Input 50:50 directional coupler ———
* a1=(signal in), a2=(dark, floating → 0), b1→arm1, b2→arm2
Xdc1  in_re in_im in_wl
+     dark_re dark_im in_wl
+     arm1_re arm1_im in_wl
+     arm2_re arm2_im in_wl
+     directional_coupler kappa_0={kappa_0} wavelength_nm={wavelength_nm}

* ——— Arm-1 PN phase shifter ———
Xps1  arm1_re arm1_im in_wl
+     ps1_re ps1_im in_wl
+     anode1 cathode1
+     pn_phase_shifter_l1
+     L_um={L_arm_um} n_g={n_g} alpha_dB_cm={alpha_dB_cm}
+     Vpi_L={Vpi_L} V_ref={V_ref} wavelength_nm={wavelength_nm}

* ——— Arm-2 PN phase shifter ———
Xps2  arm2_re arm2_im in_wl
+     ps2_re ps2_im in_wl
+     anode2 cathode2
+     pn_phase_shifter_l1
+     L_um={L_arm_um} n_g={n_g} alpha_dB_cm={alpha_dB_cm}
+     Vpi_L={Vpi_L} V_ref={V_ref} wavelength_nm={wavelength_nm}

* ——— Output 50:50 directional coupler ———
* a1=(arm1 out), a2=(arm2 out), b1→bar, b2→cross
Xdc2  ps1_re ps1_im in_wl
+     ps2_re ps2_im in_wl
+     bar_re bar_im bar_wl
+     cross_re cross_im cross_wl
+     directional_coupler kappa_0={kappa_0} wavelength_nm={wavelength_nm}

.ends mzi_pn_l1
