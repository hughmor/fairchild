* giona_pn_th_ps.inc — best-known giona neuron PN+heater phase-shifter card.
*
* Combines the best-measured piece from each dataset:
*   [May26 n7] staged CW fit vs the live-junction sparse sweep (fit_photonic.py,
*              giona_neuron7_pn_th_ps_full_fit.json):
*              l_m, n_eff, alpha_db_cm=10.7 (post-fix units; fitted 21.4 pre-fix), r_heater=184.4,
*              dn_dv=-3.62e-5, da_dv=3.29e-4
*   [Jul26 n3] dense 50x100 HC x JV sweep, fit_neuron3.py (junction path broken
*              that day - only passive/thermal trusted):
*              p_pi_th=26.4mW (clean I^2 law, 6.6pm rms over 1.6nm; neuron
*              spread ~2x - n7 May gave 53mW)
*   [FEM]      pn_modulator/pn_extracted.json: c_j0, v_bi, m_j, beta_tpa, a_eff
*
* Injection: CALIBRATED from the May live-junction data (fit_may_dn_di.py,
* giona_dn_di_fit.json): diode trio pinned from the measured IV alone (2k
* shunt subtracted), then dn_di/da_di (delta_n = -dn_di*I_fwd, delta_alpha =
* da_di*I_fwd) fit as two linear coefficients vs the forward spectra
* (0<JV<=0.9, hc=0). Residual 0.20 vs 0.98 with injection off; captures the
* forward blue-shift + notch shallowing at 0.71/0.84 V.
* CAVEATS: (1) the diode trio is ASSEMBLY-effective (junction + PCB network;
* n_diode sits at the 5.0 bound) - refit from a clean on-die IV when one
* exists, and refit dn_di/da_di with it (they are only meaningful paired with
* their I(V)). (2) Do NOT fit dn_dv_inj/da_dv_inj jointly with i_sat - they
* are degenerate through delta_n = (dn_dv_inj/i_sat)*I; that degeneracy is
* what produced the May pathology this parametrization replaced.
*
* Absolute n_eff must be trimmed per ring (fab variation exceeds anything a
* shared card can carry). Testbench (not device) notes: ring coupler
* kappa_l ~ 0.18 [May26 n7]; the PCB puts 10k in series and 2k in parallel
* with the diode - model that network in the testbench netlist, never in
* r_series.
*
* The parser has no .include yet: concatenate this file ahead of the netlist
* (fit scripts do card + netlist via load_str).
*
* Units: the historic factor-2 loss bug (sotu.md section F) is FIXED —
* alpha_db_cm is physical power dB/cm. Alphas fitted pre-fix were halved to
* convert (21.4 pre-fix == 10.7 physical).

.model giona_ps fc_pn_th_ps LEVEL=4
+ l_m=2.5378e-5 n_g=4.2 n_eff=2.2810 pin_at_ref=0 alpha_db_cm=10.7
+ r_heater=184.4 p_pi_th=26.4e-3
+ dn_dv=-3.62e-5 da_dv=3.29e-4 c_j0=1.375e-13 v_bi=0.917 m_j=0.5
+ i_sat=5.099e-8 n_diode=5.0 r_series=0 tau_carrier=10e-9
+ dn_dv_inj=0 da_dv_inj=0 dn_di=3.99 da_di=4.63e6
+ beta_tpa=7.9e-12 a_eff_m2=1.257e-13 dn_dt=1.86e-4 r_th=0
