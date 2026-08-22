* An optical link built entirely in Verilog-A — a Mach-Zehnder interferometer.
*
* eam_link.sp shows Verilog-A interoperating with native photonics.  This one
* uses no native photonic device at all: source, propagation, splitting,
* recombining and detection are all compiled Verilog-A.  Only the load
* resistors are built in.
*
*                     ┌── wg_top (L = 400 µm) ──┐
*   va_laser ── split ┤                         ├ combine ── va_photodetector
*                     └── wg_bot (L = 400+ΔL) ──┘     │
*                                                     └──── va_photodetector
*
* Both couplers are 3 dB (κ=0.5).  The arm imbalance ΔL sets the relative
* phase, so sweeping wavelength walks the two outputs through the MZI fringe:
*
*   P_bar   = P_in · cos²(Δφ/2)
*   P_cross = P_in · sin²(Δφ/2)      Δφ = 2π·n_g·ΔL/λ
*
* FSR = λ²/(n_g·ΔL) = 1550² nm / (4.2 · 20 µm) ≈ 28.6 nm, so a 1540–1570 nm
* sweep covers just over one period.  check.py asserts the two ports stay
* complementary — that is the statement that the coupler is unitary and the
* whole chain conserves power.
*
* Every optical model takes its wavelength as a PARAMETER, never off the λ
* wire — see the long note in models/optical.vams.  `.param wl` keeps the one
* number in one place.  There is no CLI override for a `.param`, so a sweep
* overrides the three element params instead:
*
*   for wl in 1540 1545 1550 1555 1560 1565 1570; do
*     fairchild -f examples/verilog_a/va_link.sp \
*       --param "Xlaser.wavelength_nm=$wl" \
*       --param "Xwgtop.wavelength_nm=$wl" --param "Xwgbot.wavelength_nm=$wl" \
*       --probe "v(bar_i),v(cross_i)" --format csv
*   done

.param wl=1550

.osdi build/va_laser.osdi
.osdi build/va_waveguide.osdi
.osdi build/va_coupler.osdi
.osdi build/va_photodetector.osdi

* Named as bundles purely so the wires get optical-discipline names; the
* Verilog-A elements address the underlying <port>_re_0 / _im_0 / _wl_0 wires.
.optical_port las
.optical_port top_in
.optical_port bot_in
.optical_port top_out
.optical_port bot_out
.optical_port bar
.optical_port cross
.optical_port dark

Xlaser  las_re_0 las_im_0 las_wl_0
+       va_laser power_mW=1.0 wavelength_nm={wl}

* Split: port 2 of the input coupler is unconnected (dark), which is what an
* undriven optical bundle means — no source, field 0.
Xsplit  las_re_0 las_im_0 las_wl_0
+       dark_re_0 dark_im_0 dark_wl_0
+       top_in_re_0 top_in_im_0 top_in_wl_0
+       bot_in_re_0 bot_in_im_0 bot_in_wl_0
+       va_coupler kappa_0=0.5 wavelength_nm={wl}

Xwgtop  top_in_re_0 top_in_im_0 top_in_wl_0
+       top_out_re_0 top_out_im_0 top_out_wl_0
+       va_waveguide L_um=400 n_g=4.2 alpha_dB_cm=2.0 wavelength_nm={wl}

* +20 µm of path length: the whole interferometer is this one number.
Xwgbot  bot_in_re_0 bot_in_im_0 bot_in_wl_0
+       bot_out_re_0 bot_out_im_0 bot_out_wl_0
+       va_waveguide L_um=420 n_g=4.2 alpha_dB_cm=2.0 wavelength_nm={wl}

Xcomb   top_out_re_0 top_out_im_0 top_out_wl_0
+       bot_out_re_0 bot_out_im_0 bot_out_wl_0
+       bar_re_0 bar_im_0 bar_wl_0
+       cross_re_0 cross_im_0 cross_wl_0
+       va_coupler kappa_0=0.5 wavelength_nm={wl}

Xpd1    bar_re_0 bar_im_0 bar_wl_0     bar_i 0    va_photodetector
+       responsivity=1.0 I_dark_A=0 R_shunt=1e9
Xpd2    cross_re_0 cross_im_0 cross_wl_0 cross_i 0 va_photodetector
+       responsivity=1.0 I_dark_A=0 R_shunt=1e9

* 1 kΩ transimpedance: V = 1 V per mW at R = 1 A/W.
Rbar    bar_i   0  1k
Rcross  cross_i 0  1k

.op
