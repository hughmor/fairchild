* Verilog-A waveguide against the native one, same parameters, one deck.
*
* Two independent chains from identical native lasers into identical native
* detectors; only the middle section differs.  They must agree exactly — the
* Verilog-A model implements the same physics, and after the loss-convention
* fix (see models/va_waveguide.va) the same convention too.
*
* 1 mm at 3 dB/cm is 0.3 dB, so 10^(-0.3/10) = 0.9333 of the power arrives.  The
* unmaintained legacy/va-models/photonic/waveguide.va gives 0.9661 here: it
* converts dB with the amplitude constant 8.6859 and then halves again in the
* exponent, counting the factor of two twice.  Native fairchild had the same
* bug until 0f689cb.  check.py pins both the agreement and the absolute value.
*
* Run:
*   fairchild -f examples/verilog_a/wg_compare.sp --probe "v(va_i),v(native_i)"

.osdi build/va_waveguide.osdi

.optical_port va_src
.optical_port va_out
.optical_port nat_src
.optical_port nat_out

* ── Verilog-A path ─────────────────────────────────────────────────────────
Xvalas  va_src fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xvawg   va_src_re_0 va_src_im_0 va_src_wl_0
+       va_out_re_0 va_out_im_0 va_out_wl_0
+       va_waveguide L_um=1000 n_g=4.2 alpha_dB_cm=3.0 wavelength_nm=1550
Xvapd   va_out va_i 0 fc_photodetector responsivity=0.8 i_dark_a=0 r_shunt=1G
Rva     va_i 0 1k

* ── native path ────────────────────────────────────────────────────────────
Xnatlas nat_src fc_cw_laser power_mW=1.0 wavelength_nm=1550
Xnatwg  nat_src nat_out fc_waveguide L_um=1000 n_g=4.2 alpha_dB_cm=3.0
Xnatpd  nat_out native_i 0 fc_photodetector responsivity=0.8 i_dark_a=0 r_shunt=1G
Rnat    native_i 0 1k

.op
