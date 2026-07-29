* Optical link with a Verilog-A modulator dropped into a native photonic chain.
*
* Native fairchild devices provide the laser, the routing waveguides, the
* photodetector and the whole electrical side.  The modulator is Verilog-A
* (models/va_eam.va) — fairchild has no native EAM.  All of it solves in one
* Newton-Raphson loop.
*
*   fc_cw_laser ─► fc_waveguide ─► [ va_eam ] ─► fc_waveguide ─► fc_photodetector ─► Rload
*                                      ▲
*                            Rdrv + Cpar (native RC, tau = 200 ps)
*                                      ▲
*                                  Vdrv pulse
*
* Bridging the two worlds is just naming: `.optical_port p` declares a bundle
* the native devices take as one token `p`, and expands it to the three wires
* `p_re_0 p_im_0 p_wl_0` that the Verilog-A module's terminals bind to.  Same
* units on both sides — sqrt(W) fields, wavelength in metres.
*
* The EAM's photocurrent flows back into the native drive network, so V(eam_a)
* does not sit exactly at the driven level: the link is coupled both ways.
*
* Run:
*   fairchild -f examples/verilog_a/eam_link.sp \
*             --probe "v(drv),v(eam_a),v(pd_out)" --format csv -o /tmp/eam.csv

.osdi build/va_eam.osdi

.optical_port las
.optical_port mod_in
.optical_port mod_out
.optical_port pd_in

* ── optical chain, native either side of the Verilog-A block ──────────────
Xlaser las fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xwg1   las mod_in  fc_waveguide L_um=200 n_g=4.2 alpha_dB_cm=2.0

Xeam   mod_in_re_0 mod_in_im_0 mod_in_wl_0
+      mod_out_re_0 mod_out_im_0 mod_out_wl_0
+      eam_a 0
+      va_eam L_um=100 il_dB=3.0 er_dB=10.0 v_full=2.0 responsivity=0.8

Xwg2   mod_out pd_in fc_waveguide L_um=200 n_g=4.2 alpha_dB_cm=2.0
Xpd    pd_in pd_out 0 fc_photodetector responsivity=0.8 i_dark_a=1e-9 r_shunt=1Meg
Rload  pd_out 0 1k

* ── native drive network: 200 ps RC in front of the modulator ─────────────
Vdrv  drv   0     PULSE(0 -2 200p 50p 50p 800p 2n)
Rdrv  drv   eam_a 200
Cpar  eam_a 0     1p

* Backward Euler here is a choice, not a workaround: the modulator's response
* is quasi-static over a 200 ps RC, so first order is plenty and it keeps the
* waveform free of trapezoidal ringing at the pulse edges.  Verilog-A `ddt`
* honours whichever method you pick.
.options method=be
.tran 5p 4n
.end
