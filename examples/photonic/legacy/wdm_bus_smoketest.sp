* WDM bus notation smoke test
* Tests .optical_bus and bus vector expansion in X element nets.
* Uses a 2-channel optical bus, each channel driven by its own laser
* and detected by its own photodetector.
*
* Expected: V(ph_a_0) ≈ 1.0, V(ph_a_1) ≈ 1.0
*
* Run:
*   fairchild -f wdm_bus_smoketest.sp --probe "V(ph_a_0),V(ph_a_1)"

.osdi ../../va-models/build/cw_laser.osdi
.osdi ../../va-models/build/photodetector.osdi

* Declare 2-channel optical bus (expands to ch_re_0 ch_im_0 ch_wl_0 ch_re_1 ch_im_1 ch_wl_1)
.optical_bus 2 ch_re ch_im ch_wl

* Channel 0: laser at 1550 nm
Xlaser0  ch_re_0 ch_im_0 ch_wl_0  cw_laser  power_mW=1.0 wavelength_nm=1550.0
Xpd0     ch_re_0 ch_im_0 ch_wl_0  ph_a_0 0  photodetector  responsivity=1.0
Rload0   ph_a_0 0  1k

* Channel 1: laser at 1544 nm (different wavelength, same channel structure)
Xlaser1  ch_re_1 ch_im_1 ch_wl_1  cw_laser  power_mW=1.0 wavelength_nm=1544.0
Xpd1     ch_re_1 ch_im_1 ch_wl_1  ph_a_1 0  photodetector  responsivity=1.0
Rload1   ph_a_1 0  1k

.op
.end
