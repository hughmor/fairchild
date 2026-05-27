* mrr_adddrop_pn_l1 — Add-drop ring resonator modulator (PN, Level 1)
*
* Composition: 2× directional_coupler + 1× pn_phase_shifter_l1 + 1× waveguide
*
* Topology (two-bus ring with feedback):
*
*    in ─── DC1 ─── thru              add ─── DC2 ─── drop
*             │                                 │
*           ring1_in → PS(active) → ring1_out ──┘
*             │                                 │
*           ring2_out ←── waveguide ←── ring2_in
*
* DC1 couples the input bus to the ring; DC2 couples the ring to the drop bus.
* The ring splits into an active section (PN phase shifter) and a passive
* return path (straight waveguide).  Total ring circumference = L_active + L_passive.
*
* Feedback loop: ring2_out → DC1.a2; DC2.b2 → ring2_in → waveguide → ring2_out.
*
* Ports:
*   in_re in_im in_wl      — input bus, input port  (3-wire)
*   thru_re thru_im thru_wl — input bus, through port (3-wire)
*   drop_re drop_im drop_wl — drop bus, drop port    (3-wire)
*   add_re add_im add_wl   — drop bus, add port (float for standard filter)
*   anode cathode           — PN junction (electrical)
*
* Parameters (defaults):
*   kappa_0=0.1         bus–ring coupling fraction (both couplers)
*   L_active_um=50      active (PN) ring section length (µm)
*   L_passive_um=50     passive ring section length (µm)
*   n_g=4.2             group index
*   alpha_dB_cm=2.0     ring loss (dB/cm)
*   Vpi_L=2.0           PN Vπ·L product (V·cm)
*   V_ref=0.0           zero-phase bias (V)
*   wavelength_nm=1550  design wavelength (nm)

.subckt mrr_adddrop_pn_l1
+ in_re in_im in_wl
+ thru_re thru_im thru_wl
+ drop_re drop_im drop_wl
+ add_re add_im add_wl
+ anode cathode
+ kappa_0=0.1 L_active_um=50 L_passive_um=50
+ n_g=4.2 alpha_dB_cm=2.0
+ Vpi_L=2.0 V_ref=0.0 wavelength_nm=1550.0

* ——— Input-bus coupler (DC1) ———
* a1=(bus in), a2=(ring2 output, feedback), b1=(thru), b2=(ring1 input)
Xdc1  in_re in_im in_wl
+     ring2_out_re ring2_out_im in_wl
+     thru_re thru_im thru_wl
+     ring1_in_re ring1_in_im in_wl
+     directional_coupler kappa_0={kappa_0} wavelength_nm={wavelength_nm}

* ——— Active ring section: PN phase shifter ———
Xps   ring1_in_re ring1_in_im in_wl
+     ring1_out_re ring1_out_im in_wl
+     anode cathode
+     pn_phase_shifter_l1
+     L_um={L_active_um} n_g={n_g} alpha_dB_cm={alpha_dB_cm}
+     Vpi_L={Vpi_L} V_ref={V_ref} wavelength_nm={wavelength_nm}

* ——— Drop-bus coupler (DC2) ———
* a1=(ring1 output), a2=(add port, float for filter), b1=(drop), b2=(ring2 input)
Xdc2  ring1_out_re ring1_out_im in_wl
+     add_re add_im add_wl
+     drop_re drop_im drop_wl
+     ring2_in_re ring2_in_im in_wl
+     directional_coupler kappa_0={kappa_0} wavelength_nm={wavelength_nm}

* ——— Passive ring return path ———
Xwg   ring2_in_re ring2_in_im in_wl
+     ring2_out_re ring2_out_im in_wl
+     waveguide
+     L_um={L_passive_um} n_g={n_g} alpha_dB_cm={alpha_dB_cm}
+     wavelength_nm={wavelength_nm}

.ends mrr_adddrop_pn_l1
