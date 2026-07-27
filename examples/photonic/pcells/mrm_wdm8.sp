* mrm_wdm8.sp — 8-channel add-drop micro-ring modulator PCell.
*
* The WDM sibling of mrm.sp: same topology, same parameters, same fitted
* defaults — but ONE instance serving all eight wavelengths off ONE PN junction
* and ONE heater, which is what a real ring on a WDM bus is.
*
*   .include examples/photonic/pcells/mrm_wdm8.sp
*   .optical_port bus 8
*   Xr1 bus_in bus_th bus_ad bus_dr  vpn 0  heat 0  mrm_wdm8 radius=8e-6
*
* Why a separate file rather than instancing mrm.sp on an 8-channel bundle:
* mrm.sp declares one channel's worth of ports, so that connection describes
* eight independent rings, each with its own junction and heater, all wired to
* the same two electrical nodes — eight times the terminal current. The parser
* rejects it (see BundleArity::Scalar); it used to replicate silently.
*
* Ports: in / th / ad / dr, each 8 channels x (re, im, wl), then
*        pn_a pn_c ht_p ht_n.
*
* See mrm.sp for the topology diagram, the per-arc vs whole-ring parameter
* conventions (r_heater and the junctions are per arc, p_pi_th is whole-ring),
* and the fit provenance. Everything there applies unchanged.

.subckt mrm_wdm8
+ in_re_0 in_im_0 in_wl_0 in_re_1 in_im_1 in_wl_1 in_re_2 in_im_2 in_wl_2
+ in_re_3 in_im_3 in_wl_3 in_re_4 in_im_4 in_wl_4 in_re_5 in_im_5 in_wl_5
+ in_re_6 in_im_6 in_wl_6 in_re_7 in_im_7 in_wl_7 th_re_0 th_im_0 th_wl_0
+ th_re_1 th_im_1 th_wl_1 th_re_2 th_im_2 th_wl_2 th_re_3 th_im_3 th_wl_3
+ th_re_4 th_im_4 th_wl_4 th_re_5 th_im_5 th_wl_5 th_re_6 th_im_6 th_wl_6
+ th_re_7 th_im_7 th_wl_7 ad_re_0 ad_im_0 ad_wl_0 ad_re_1 ad_im_1 ad_wl_1
+ ad_re_2 ad_im_2 ad_wl_2 ad_re_3 ad_im_3 ad_wl_3 ad_re_4 ad_im_4 ad_wl_4
+ ad_re_5 ad_im_5 ad_wl_5 ad_re_6 ad_im_6 ad_wl_6 ad_re_7 ad_im_7 ad_wl_7
+ dr_re_0 dr_im_0 dr_wl_0 dr_re_1 dr_im_1 dr_wl_1 dr_re_2 dr_im_2 dr_wl_2
+ dr_re_3 dr_im_3 dr_wl_3 dr_re_4 dr_im_4 dr_wl_4 dr_re_5 dr_im_5 dr_wl_5
+ dr_re_6 dr_im_6 dr_wl_6 dr_re_7 dr_im_7 dr_wl_7 pn_a pn_c ht_p ht_n
+ radius=8e-6 n_g=4.2 n_eff=2.2810 alpha_db_cm=10.7 kappa_l=0.183
+ wl_ref_nm=1550
+ r_heater=184.4 p_pi_th=26.4e-3 dn_dt=1.86e-4 r_th=0
+ dn_dv=-3.62e-5 da_dv=3.29e-4 c_j0=1.375e-13 v_bi=0.917 m_j=0.5
+ i_sat=5.099e-8 n_diode=5.0 r_series=0 tau_carrier=10e-9
+ dn_di=3.99 da_di=4.63e6 dn_dv_inj=0 da_dv_inj=0
+ beta_tpa=7.9e-12 a_eff_m2=1.257e-13 pin_at_ref=0

* Per-instance phase-shifter card, built from this instance's parameters —
* {pi*radius} is the arc length, half the circumference.
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
Xcplb
+ in_re_0 in_im_0 in_wl_0 in_re_1 in_im_1 in_wl_1 in_re_2 in_im_2 in_wl_2
+ in_re_3 in_im_3 in_wl_3 in_re_4 in_im_4 in_wl_4 in_re_5 in_im_5 in_wl_5
+ in_re_6 in_im_6 in_wl_6 in_re_7 in_im_7 in_wl_7 ra_re_0 ra_im_0 ra_wl_0
+ ra_re_1 ra_im_1 ra_wl_1 ra_re_2 ra_im_2 ra_wl_2 ra_re_3 ra_im_3 ra_wl_3
+ ra_re_4 ra_im_4 ra_wl_4 ra_re_5 ra_im_5 ra_wl_5 ra_re_6 ra_im_6 ra_wl_6
+ ra_re_7 ra_im_7 ra_wl_7 th_re_0 th_im_0 th_wl_0 th_re_1 th_im_1 th_wl_1
+ th_re_2 th_im_2 th_wl_2 th_re_3 th_im_3 th_wl_3 th_re_4 th_im_4 th_wl_4
+ th_re_5 th_im_5 th_wl_5 th_re_6 th_im_6 th_wl_6 th_re_7 th_im_7 th_wl_7
+ rb_re_0 rb_im_0 rb_wl_0 rb_re_1 rb_im_1 rb_wl_1 rb_re_2 rb_im_2 rb_wl_2
+ rb_re_3 rb_im_3 rb_wl_3 rb_re_4 rb_im_4 rb_wl_4 rb_re_5 rb_im_5 rb_wl_5
+ rb_re_6 rb_im_6 rb_wl_6 rb_re_7 rb_im_7 rb_wl_7 fc_dcoupler
+ kappa_L={kappa_l}
Xcpld
+ rc_re_0 rc_im_0 rc_wl_0 rc_re_1 rc_im_1 rc_wl_1 rc_re_2 rc_im_2 rc_wl_2
+ rc_re_3 rc_im_3 rc_wl_3 rc_re_4 rc_im_4 rc_wl_4 rc_re_5 rc_im_5 rc_wl_5
+ rc_re_6 rc_im_6 rc_wl_6 rc_re_7 rc_im_7 rc_wl_7 ad_re_0 ad_im_0 ad_wl_0
+ ad_re_1 ad_im_1 ad_wl_1 ad_re_2 ad_im_2 ad_wl_2 ad_re_3 ad_im_3 ad_wl_3
+ ad_re_4 ad_im_4 ad_wl_4 ad_re_5 ad_im_5 ad_wl_5 ad_re_6 ad_im_6 ad_wl_6
+ ad_re_7 ad_im_7 ad_wl_7 rd_re_0 rd_im_0 rd_wl_0 rd_re_1 rd_im_1 rd_wl_1
+ rd_re_2 rd_im_2 rd_wl_2 rd_re_3 rd_im_3 rd_wl_3 rd_re_4 rd_im_4 rd_wl_4
+ rd_re_5 rd_im_5 rd_wl_5 rd_re_6 rd_im_6 rd_wl_6 rd_re_7 rd_im_7 rd_wl_7
+ dr_re_0 dr_im_0 dr_wl_0 dr_re_1 dr_im_1 dr_wl_1 dr_re_2 dr_im_2 dr_wl_2
+ dr_re_3 dr_im_3 dr_wl_3 dr_re_4 dr_im_4 dr_wl_4 dr_re_5 dr_im_5 dr_wl_5
+ dr_re_6 dr_im_6 dr_wl_6 dr_re_7 dr_im_7 dr_wl_7 fc_dcoupler
+ kappa_L={kappa_l}

* Arc 1: drop coupler -> bus coupler.  Arc 2: bus coupler -> drop coupler.
* Heaters in series through the internal node ht_mid.
Xps1
+ rd_re_0 rd_im_0 rd_wl_0 rd_re_1 rd_im_1 rd_wl_1 rd_re_2 rd_im_2 rd_wl_2
+ rd_re_3 rd_im_3 rd_wl_3 rd_re_4 rd_im_4 rd_wl_4 rd_re_5 rd_im_5 rd_wl_5
+ rd_re_6 rd_im_6 rd_wl_6 rd_re_7 rd_im_7 rd_wl_7 ra_re_0 ra_im_0 ra_wl_0
+ ra_re_1 ra_im_1 ra_wl_1 ra_re_2 ra_im_2 ra_wl_2 ra_re_3 ra_im_3 ra_wl_3
+ ra_re_4 ra_im_4 ra_wl_4 ra_re_5 ra_im_5 ra_wl_5 ra_re_6 ra_im_6 ra_wl_6
+ ra_re_7 ra_im_7 ra_wl_7 pn_a pn_c ht_p ht_mid arc_ps
Xps2
+ rb_re_0 rb_im_0 rb_wl_0 rb_re_1 rb_im_1 rb_wl_1 rb_re_2 rb_im_2 rb_wl_2
+ rb_re_3 rb_im_3 rb_wl_3 rb_re_4 rb_im_4 rb_wl_4 rb_re_5 rb_im_5 rb_wl_5
+ rb_re_6 rb_im_6 rb_wl_6 rb_re_7 rb_im_7 rb_wl_7 rc_re_0 rc_im_0 rc_wl_0
+ rc_re_1 rc_im_1 rc_wl_1 rc_re_2 rc_im_2 rc_wl_2 rc_re_3 rc_im_3 rc_wl_3
+ rc_re_4 rc_im_4 rc_wl_4 rc_re_5 rc_im_5 rc_wl_5 rc_re_6 rc_im_6 rc_wl_6
+ rc_re_7 rc_im_7 rc_wl_7 pn_a pn_c ht_mid ht_n arc_ps
.ends
