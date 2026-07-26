* mrm.sp — add-drop micro-ring modulator PCell.
*
* Four optical ports (in / thru / add / drop) and four electrical terminals
* (PN anode/cathode, heater +/-). Every sub-device parameter is exposed, so an
* instance is a fully-specified device; defaults are the fitted giona neuron
* values (see experiments/giona/giona_pn_th_ps.inc for provenance).
*
*   .include examples/photonic/pcells/mrm.sp
*   Xr1 in_re in_im in_wl  th_re th_im th_wl  ad_re ad_im ad_wl  dr_re dr_im dr_wl
*   +   vpn 0  heat 0  mrm  radius=8e-6 kappa_l=0.183 dn_dv=-3.62e-5
*
* SINGLE CHANNEL by design. Hand the instance a multi-channel .optical_port and
* the parser replicates it once per wavelength (12 optical ports = one channel's
* worth), which is the right thing for a ring bank: pair it with an
* .electrical_port so each ring also gets its own drive wire. What you must NOT
* do is expect one instance to serve N wavelengths off one junction — that is
* what fc_optical_2x2 is for.
*
* Topology — add-drop, two couplers with the ring split into two arcs, exactly
* the form the parameters were fitted against (experiments/giona/ringfit.py):
*
*     IN  --> CPLB ---------------------> THRU
*              |  ^
*            PS2  PS1        each arc = half the ring: PN junction + heater
*              v  |
*     DROP <-- CPLD <-- ADD
*
* Both arcs are driven from the same anode/cathode and their heaters sit in
* SERIES, so one current drives both. Consequences worth knowing:
*   - l_m per arc is pi*radius (half the circumference), derived below.
*   - p_pi_th is a whole-ring number: heaters in series put half the total power
*     in each arc, and two arcs each contributing pi*P_arc/p_pi sums to
*     pi*P_total/p_pi. So the fitted 26.4 mW carries over unchanged.
*   - r_heater is PER ARC; the terminal resistance is 2*r_heater.
*   - i_sat / dn_di / da_di are PER ARC, matching the fit. The two junctions are
*     electrically in parallel, so the terminal current is 2x a single arc's.
*     That is inherited from the fit, not introduced here — the pending refit
*     from a clean on-die IV (see the card header) is what fixes it. Do not
*     "correct" i_sat alone: it is only meaningful paired with dn_di/da_di.
*
* LEVEL=4 (depletion + current-driven injection + TPA) needs a .model card, and
* the card is declared INSIDE the subckt so every instance gets its own copy
* built from its own parameters. That is what makes this a PCell rather than a
* fixed cell.

.subckt mrm
+ in_re in_im in_wl  th_re th_im th_wl  ad_re ad_im ad_wl  dr_re dr_im dr_wl
+ pn_a pn_c  ht_p ht_n
* ── geometry ────────────────────────────────────────────────────────────────
+ radius=8e-6 n_g=4.2 n_eff=2.2810 alpha_db_cm=10.7 kappa_l=0.183
+ wl_ref_nm=1550
* ── heater ──────────────────────────────────────────────────────────────────
+ r_heater=184.4 p_pi_th=26.4e-3 dn_dt=1.86e-4 r_th=0
* ── PN depletion ────────────────────────────────────────────────────────────
+ dn_dv=-3.62e-5 da_dv=3.29e-4 c_j0=1.375e-13 v_bi=0.917 m_j=0.5
* ── PN diode + carrier injection ────────────────────────────────────────────
+ i_sat=5.099e-8 n_diode=5.0 r_series=0 tau_carrier=10e-9
+ dn_di=3.99 da_di=4.63e6 dn_dv_inj=0 da_dv_inj=0
* ── nonlinear optics ────────────────────────────────────────────────────────
+ beta_tpa=7.9e-12 a_eff_m2=1.257e-13 pin_at_ref=0

* Per-instance phase-shifter card. Every parameter above that the device
* understands is forwarded; {pi*radius} is the arc length.
.model arc_ps fc_pn_th_ps LEVEL=4
+ l_m={pi*radius} n_g={n_g} n_eff={n_eff} alpha_db_cm={alpha_db_cm}
+ wl_ref_nm={wl_ref_nm}
+ pin_at_ref={pin_at_ref}
+ r_heater={r_heater} p_pi_th={p_pi_th} dn_dt={dn_dt} r_th={r_th}
+ dn_dv={dn_dv} da_dv={da_dv} c_j0={c_j0} v_bi={v_bi} m_j={m_j}
+ i_sat={i_sat} n_diode={n_diode} r_series={r_series} tau_carrier={tau_carrier}
+ dn_di={dn_di} da_di={da_di} dn_dv_inj={dn_dv_inj} da_dv_inj={da_dv_inj}
+ beta_tpa={beta_tpa} a_eff_m2={a_eff_m2}

* Bus coupler (in/thru side) and drop coupler (add/drop side).
Xcplb in_re in_im in_wl  ra_re ra_im ra_wl  th_re th_im th_wl  rb_re rb_im rb_wl
+     fc_dcoupler kappa_L={kappa_l}
Xcpld rc_re rc_im rc_wl  ad_re ad_im ad_wl  rd_re rd_im rd_wl  dr_re dr_im dr_wl
+     fc_dcoupler kappa_L={kappa_l}

* Arc 1: drop coupler -> bus coupler.  Arc 2: bus coupler -> drop coupler.
* Heaters in series through the internal node ht_mid.
Xps1 rd_re rd_im rd_wl  ra_re ra_im ra_wl  pn_a pn_c  ht_p  ht_mid  arc_ps
Xps2 rb_re rb_im rb_wl  rc_re rc_im rc_wl  pn_a pn_c  ht_mid ht_n   arc_ps
.ends
