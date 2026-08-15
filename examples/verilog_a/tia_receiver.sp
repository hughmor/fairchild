* tia_receiver — a native photodiode read by a Verilog-A transimpedance amp.
*
* The receiver front end a real link uses, and the reason it exists: a load
* resistor big enough to give a readable voltage is slow and noisy, because
* 1/(2*pi*R*C) and 4kT/R pull in opposite directions. A TIA breaks that trade
* by holding the summing node at a virtual ground.
*
* Everything optical here is native; everything electrical after the diode is
* Verilog-A. They share one Newton iteration — the photocurrent reaches the
* amplifier and the amplifier's input impedance loads the diode, in the same
* solve.
*
* The interesting number is `i_n_in`: 15 pA/sqrt(Hz) of input-referred current
* noise, which `.noise` sees through OSDI's load_noise and which sets this
* receiver's floor. `check.py` pins the output PSD against i_n_in * z_t.

.osdi build/va_tia.osdi

.optical_port beam

* 20 uW landing on the diode — a real receiver's input, not a demo's.
Xlas  beam fc_cw_laser power_mW=0.02
Xpd   beam sum 0 fc_photodetector responsivity=0.9 r_shunt=1Meg i_dark_a=10n

* z_t = 2 kohm into 10 GHz, 50 ohm in and out.
Xtia  sum tout 0 va_tia z_t=2000 r_in=50 f_3db=10e9 i_n_in=15e-12
+     v_out_dc=0 v_swing=1.0 r_out=50

Rload tout 0 50

* An AC source that injects at the summing node, so `.noise` has an input to
* refer to. 1 Mohm makes it a current source of 1 uA/V and loads nothing.
Vac   acs 0 DC 0 AC 1
Racs  acs sum 1meg

.op
.noise V(tout) Vac dec 10 1e6 1e11
.end
