* WDM micro-ring modulator: two wavelength channels through one ring.
*
* Same MRR topology as native_mrr_modulator.sp, but the input bus carries
* TWO independent wavelength channels.  Both channels share the physical
* ring (and therefore the same V_pn modulation), but each laser sees a
* different transmission profile because the channels' detunings from the
* ring resonance are different — and they evolve differently as the ring
* resonance walks across the spectrum with V_pn.
*
* No multiplexer or demultiplexer device is needed.  Each wavelength is its
* own 3-wire bundle (re, im, λ); the bundle-port directive declares an N=2
* bus and the parser replicates every photonic device along the bus into
* two parallel single-channel instances.  Two lasers wire to channel 0 and
* channel 1 by explicit underlying-wire names ("MUX"); two photodetectors
* read each channel's output the same way ("DEMUX").
*
* Laser placement: symmetric ±50 pm around the ring's reference wavelength
* (1550 nm).  As V_pn ramps 0→4 V the ring resonance walks blue by ≈ 570 pm
* (≈ half-FSR), so:
*   - Channel 0 (1549.95 nm) is initially −50 pm from resonance; resonance
*     reaches the laser around V_pn ≈ 0.35 V, producing a sharp notch.
*   - Channel 1 (1550.05 nm) is initially +50 pm; resonance moves AWAY as
*     V_pn rises, so the channel monotonically rises to high transmission
*     with no notch event.

.optical_port bus_in 2
.optical_port wg1_out 2
.optical_port dc_b 2
.optical_port dc_c 2
.optical_port pn_in 2
.optical_port pd_in 2

* Two lasers, each driving one channel of the input bus via explicit wires.
* This is the "MUX" — a real photonic chip would use an AWG here, but in
* the SVEA model each channel is structurally independent so a naming
* convention suffices.
Xlaser1 bus_in_re_0 bus_in_im_0 bus_in_wl_0 fc_cw_laser
+       power_mW=1.0 wavelength_nm=1549.95
Xlaser2 bus_in_re_1 bus_in_im_1 bus_in_wl_1 fc_cw_laser
+       power_mW=1.0 wavelength_nm=1550.05

* Coupling waveguide: bundle ports → parser replicates per channel.
Xwg1 bus_in wg1_out fc_waveguide
+    L_um=50 n_g=4.2 alpha_dB_cm=2.0 wavelength_nm=1550

* 2×2 directional coupler — replicates per channel (so each wavelength
* sees its own dedicated coupler with identical κL).
Xdc wg1_out dc_b dc_c pn_in fc_dcoupler kappa_L=0.336

* The ring (PN phase shifter) — replicates per channel, but vmod and the
* ground reference are plain nets so both rings share the same drive.
* That's the whole point: one physical ring, multiple wavelengths.
Xpn pn_in dc_b vmod 0 fc_pn_ps
+   L_um=500 V_pi_L=2e-3 g_pn=1e-3 alpha_dB_cm=10
+   n_g=4.2 wavelength_nm=1550

* Output waveguide
Xwg2 dc_c pd_in fc_waveguide
+    L_um=50 n_g=4.2 alpha_dB_cm=2.0 wavelength_nm=1550

* Two photodetectors, one per channel — the "DEMUX".  Explicit wire names
* on each instance line so each PD reads exactly one wavelength channel.
Xpd1 pd_in_re_0 pd_in_im_0 pd_in_wl_0 pd1_anode 0 fc_photodetector
+    responsivity=0.8 i_dark_a=1e-9 r_shunt=1Meg
Xpd2 pd_in_re_1 pd_in_im_1 pd_in_wl_1 pd2_anode 0 fc_photodetector
+    responsivity=0.8 i_dark_a=1e-9 r_shunt=1Meg

Vbias  bias 0 DC 1.0
Rload1 pd1_anode bias 1k
Rload2 pd2_anode bias 1k

* Modulation: 0 → 4 V pulse sweeps the ring resonance across the bus.
Vmod vmod 0 PULSE(0 4 100n 100n 100n 800n 2u)

.options method=gear
.tran 5n 2u
.end
