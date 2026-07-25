* giona_pn_th_ps.inc — best-known giona neuron PN+heater phase-shifter card.
*
* Combines the best-measured piece from each dataset:
*   [May26 n7] staged CW fit vs the live-junction sparse sweep (fit_photonic.py,
*              giona_neuron7_pn_th_ps_full_fit.json):
*              l_m, n_eff, alpha_db_cm=21.4, r_heater=184.4,
*              dn_dv=-3.62e-5, da_dv=3.29e-4
*   [Jul26 n3] dense 50x100 HC x JV sweep, fit_neuron3.py (junction path broken
*              that day - only passive/thermal trusted):
*              p_pi_th=26.4mW (clean I^2 law, 6.6pm rms over 1.6nm; neuron
*              spread ~2x - n7 May gave 53mW)
*   [FEM]      pn_modulator/pn_extracted.json: c_j0, v_bi, m_j, beta_tpa, a_eff
*
* NOT calibrated (needs a live-junction dense forward sweep): the injection
* terms and the diode electricals. They are set to physical placeholders
* (i_sat=1e-12, n_diode=1.05, r_series=60) with dn_di=da_di=0. Fit workflow
* when the data lands: (1) pin i_sat/n_diode/r_series from the measured IV
* alone, then (2) fit dn_di/da_di (delta_n = -dn_di*I_fwd, delta_alpha =
* da_di*I_fwd) as two linear coefficients vs the optical response. Do NOT fit
* dn_dv_inj/da_dv_inj jointly with i_sat - they are degenerate through
* delta_n = (dn_dv_inj/i_sat)*I, which is what produced n_diode=5-at-bound and
* i_sat=1e-4 in the May fit.
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
* UNIT WARNING (known factor-2 issue, _notes/sotu.md): fairchild currently
* attenuates power by HALF the specified alpha_db_cm (dB->field-neper converter
* and the segment's exp(-aL/2) both halve). All alphas in this card were FITTED
* through fairchild, so in-sim behaviour is correct; but alpha_db_cm=21.4 here
* corresponds to ~10.7 physical dB/cm. Do not paste datasheet dB/cm values
* until the converter is fixed (then halve this card's alphas).

.model giona_ps fc_pn_th_ps LEVEL=4
+ l_m=2.5378e-5 n_g=4.2 n_eff=2.2810 pin_at_ref=0 alpha_db_cm=21.4
+ r_heater=184.4 p_pi_th=26.4e-3
+ dn_dv=-3.62e-5 da_dv=3.29e-4 c_j0=1.375e-13 v_bi=0.917 m_j=0.5
+ i_sat=1e-12 n_diode=1.05 r_series=60 tau_carrier=10e-9
+ dn_dv_inj=0 da_dv_inj=0 dn_di=0 da_di=0
+ beta_tpa=7.9e-12 a_eff_m2=1.257e-13 dn_dt=1.86e-4 r_th=0
