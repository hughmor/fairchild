* CMIM analogue link — the feedforward path of the cascaded-modulator Ising machine
*
* Al-Kayed et al., Nature 648, 576 (2025), Figure 1 and Extended Data Figure 1.
* Every device value is one of that paper's numbers or a stated assumption.  See
* ../../SPECS.md for the citation of each one.
*
*   DFB ─► [MZM1: spins x] ─► [MZM2: weights w] ─► 50:50 ─► [QD SOA] ─► [PD+TIA]
*                                                    ▲
*                                              [bulk SOA ASE]
*
* The first line of a SPICE deck is a title, so this comment is load bearing.
*
* Drive.  Vx and Vw are the two AWG channels ahead of the drivers.  They are DC
* sources here so that .op and .ac work out of the box.  common/link.py rewrites
* those two lines into PWL when it needs a waveform, which is why each is
* written on one line with no continuation.
*
* Run:
*   fairchild -f experiments/cmim/link/netlists/cmim_link.sp --probe "v(rxout)"

* Sources, not artefacts.  Four of these five declare `optical_bundle` ports,
* whose width the deck decides, so fairchild expands the source for the width in
* use and compiles it.  There is nothing to build by hand and no build/ to keep
* in step.  `optical.vams` is found beside the models through their `include`.
.va ../../models/va_tfln_mzm.va
.va ../../models/va_qd_soa.va
.va ../../models/va_ase_source.va
.va ../../models/va_mod_driver.va
.va ../../models/va_receiver.va

* ── operating conditions ─────────────────────────────────────────────────
* The machine runs at room temperature with the laser and the QD SOA held at
* 20 C (S1.2).  SPICE defaults to 27 C, which moves the amplifier's gain peak
* by 4 nm and costs 0.8 dB of gain, so this line is not cosmetic.
* The band centre is O-band here, and it decides the wavelength of any port no
* source reaches — the noise arm is one.
.options temp=20 lambda_center_nm=1310

* ── knobs ────────────────────────────────────────────────────────────────
.param p_laser_mW = 4.47      $ 6.5 dBm, S1.2
.param lambda_nm  = 1310
.param v_x        = 0         $ spin channel drive into the driver (V)
.param v_w        = 0         $ weight channel drive into the driver (V)
.param v_th_x     = 0         $ MZM1 heater bias (V), 0 keeps the quadrature default
.param v_th_w     = 0         $ MZM2 heater bias (V)
.param i_noise_mA = 0         $ bulk SOA bias, the annealing knob
.param p_ase_ref  = 0         $ bulk SOA ASE at i_ref, dBm — a fitted constant
* Driver gains, referred to the AWG's full scale.  The driver's 50 ohm source
* into the modulator's 50 ohm termination halves the swing, so a_v_x = 1.5 puts
* full scale at 0.75 V, which is pi/2 of phase: the full nonlinear half-wave,
* reached and not passed.  a_v_w = 0.4 puts the weight channel at 0.42 rad,
* within 3 % of linear.  The swings sit well above both, so the soft clip
* limits an overdrive rather than costing gain at the nominal drive.
.param a_v_x      = 1.5       $ spin driver gain: spans the nonlinear half-wave
.param a_v_w      = 0.4       $ weight driver gain: stays linear
.param sw_x       = 4.0       $ spin driver output swing, half of peak-to-peak
.param sw_w       = 4.0       $ weight driver output swing

* The bias the receiver's noise densities are evaluated at.  Zero means no
* noise.  `Link.noise_params()` measures them from a first pass over this same
* deck and hands them back, because fairchild freezes a Verilog-A noise density
* that is written against a solved quantity.  See models/va_receiver.va.
.param i_ph_n     = 0         $ photocurrent for shot noise (A)
.param p_sig_n    = 0         $ signal power for the beat terms (W)
.param rho_ase_n  = 0         $ ASE density for the beat terms (W/Hz)

.optical_port las
.optical_port m1o
.optical_port m2o
.optical_port ase
.optical_port comb
.optical_port dark

* ── source ───────────────────────────────────────────────────────────────
Xlaser las fc_cw_laser power_mW={p_laser_mW} wavelength_nm={lambda_nm}

* ── spin channel: AWG ch1 ► driver ► MZM1 ────────────────────────────────
Vx    xin 0 DC {v_x}
Xdrvx xin xdrv 0 va_mod_driver a_v={a_v_x} v_swing={sw_x}
Vthx  thx 0 DC {v_th_x}
Xmzm1 las m1o xdrv 0 thx 0 va_tfln_mzm

* ── weight channel: AWG ch2 ► driver ► MZM2 ──────────────────────────────
Vw    win 0 DC {v_w}
Xdrvw win wdrv 0 va_mod_driver a_v={a_v_w} v_swing={sw_w}
Vthw  thw 0 DC {v_th_w}
Xmzm2 m1o m2o wdrv 0 thw 0 va_tfln_mzm

* ── controllable optical noise, combined 50:50 ───────────────────────────
* The signal pays 3 dB at the combiner to let the noise arm in.  kappa_L = pi/4
* is the 3 dB point.  The noise arm's own optical field is identically zero —
* spontaneous emission has no mean amplitude — so the dark port carries it and
* the arm itself contributes only a spectral density, on the `asein` net.
Xcomb m2o dark comb ase fc_dcoupler kappa_L=0.785398163397

* asein and asesoa carry watts per hertz, not volts.  See the model headers.
Xase  asein va_ase_source i_bias_mA={i_noise_mA} p_ase_ref_dBm={p_ase_ref}
+     wavelength_nm={lambda_nm}

* ── quantum dot amplifier ────────────────────────────────────────────────
.optical_port soao
Vsoa  soa_a 0 DC 1.7
Xsoa  comb soao asein asesoa soa_a 0 va_qd_soa

* ── receiver ─────────────────────────────────────────────────────────────
Xrx   soao asesoa rxout 0 va_receiver
+     i_ph_n={i_ph_n} p_sig_n={p_sig_n} rho_ase_n={rho_ase_n}
Rload rxout 0 1Meg

.op
