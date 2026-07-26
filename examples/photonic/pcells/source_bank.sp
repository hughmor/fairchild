* source_bank.sp — 8-channel WDM signal-generation PCell.
*
* Eight CW lasers, each through its own ideal Mach-Zehnder modulator, all
* multiplexed onto one 8-channel optical bus. Use it as the stimulus end of a
* photonic link: apply a waveform to each electrical drive and you get eight
* independently modulated wavelength channels on one output.
*
*   .include examples/photonic/pcells/source_bank.sp
*   .optical_port src 8
*   Xsrc src  d1 d2 d3 d4 d5 d6 d7 d8  0  source_bank
*   +      wl1=1546.0 wl2=1547.0 ... p1=1.0 p2=1.0 ...
*
* Ports — 24 optical (the 8-channel bus, 3 wires per channel in channel order)
* then 8 electrical drives then one shared drive return:
*
*   o0r o0i o0w  o1r o1i o1w  …  o7r o7i o7w   d1 … d8   gnd
*
* Pass the bus as a single 8-channel `.optical_port` token and the parser
* flattens it into those 24 wires for ONE instance (the port count decides:
* 24+9 = 33 matches this subckt, so it flattens rather than replicating).
*
*                    ┌──────┐   ┌──────┐
*   laser λ1 (p1) ──►│ MZM1 │──►│      │
*   laser λ2 (p2) ──►│ MZM2 │──►│ mux  │──► out (8-channel bus)
*        …           │  …   │   │      │
*   laser λ8 (p8) ──►│ MZM8 │──►│      │
*                    └───▲──┘   └──────┘
*                      d1…d8
*
* DRIVE POLARITY. fc_mzm intensity transmission is
*   T(V) = alpha · [ (1 − 1/E_r)·(1 + cos(pi·V/v_pi))/2 + 1/E_r ]
* so V = 0 is FULLY ON and V = v_pi is fully off. Drive 0 → v_pi for an
* inverted-NRZ eye, or bias at v_pi/2 and swing ±v_pi/2 for a linear-ish drive.
*
* Turning a channel off: set its power to zero (`p3=0`). The laser, its
* modulator and its mux slot all remain, so the channel count and every wire
* name stay put — the channel simply carries no light.
*
* The mux is a topology marker, not a physical device: each wavelength is its
* own bundle channel, so muxing is identity routing per channel and costs no
* loss. Add `fc_grating_coupler` or a lossy waveguide if you want the real
* insertion loss of an AWG.

.subckt source_bank
+ o0r o0i o0w  o1r o1i o1w  o2r o2i o2w  o3r o3i o3w
+ o4r o4i o4w  o5r o5i o5w  o6r o6i o6w  o7r o7i o7w
+ d1 d2 d3 d4 d5 d6 d7 d8  gnd
* ── laser wavelengths (nm) — 100 GHz grid around 1550 nm by default ─────────
+ wl1=1546.12 wl2=1546.92 wl3=1547.72 wl4=1548.51
+ wl5=1549.32 wl6=1550.12 wl7=1550.92 wl8=1551.72
* ── laser powers (mW) — set one to 0 to turn that channel off ──────────────
+ p1=1.0 p2=1.0 p3=1.0 p4=1.0 p5=1.0 p6=1.0 p7=1.0 p8=1.0
* ── modulator: ideal by default (no loss, ~90 dB extinction) ────────────────
+ v_pi=1.0 il_db=0 e_r=1e9 f_c=1e11

Xl1 a1r a1i a1w fc_cw_laser power_mW={p1} wavelength_nm={wl1}
Xl2 a2r a2i a2w fc_cw_laser power_mW={p2} wavelength_nm={wl2}
Xl3 a3r a3i a3w fc_cw_laser power_mW={p3} wavelength_nm={wl3}
Xl4 a4r a4i a4w fc_cw_laser power_mW={p4} wavelength_nm={wl4}
Xl5 a5r a5i a5w fc_cw_laser power_mW={p5} wavelength_nm={wl5}
Xl6 a6r a6i a6w fc_cw_laser power_mW={p6} wavelength_nm={wl6}
Xl7 a7r a7i a7w fc_cw_laser power_mW={p7} wavelength_nm={wl7}
Xl8 a8r a8i a8w fc_cw_laser power_mW={p8} wavelength_nm={wl8}

Xm1 a1r a1i a1w b1r b1i b1w d1 gnd fc_mzm v_pi={v_pi} il_db={il_db} e_r={e_r} f_c={f_c}
Xm2 a2r a2i a2w b2r b2i b2w d2 gnd fc_mzm v_pi={v_pi} il_db={il_db} e_r={e_r} f_c={f_c}
Xm3 a3r a3i a3w b3r b3i b3w d3 gnd fc_mzm v_pi={v_pi} il_db={il_db} e_r={e_r} f_c={f_c}
Xm4 a4r a4i a4w b4r b4i b4w d4 gnd fc_mzm v_pi={v_pi} il_db={il_db} e_r={e_r} f_c={f_c}
Xm5 a5r a5i a5w b5r b5i b5w d5 gnd fc_mzm v_pi={v_pi} il_db={il_db} e_r={e_r} f_c={f_c}
Xm6 a6r a6i a6w b6r b6i b6w d6 gnd fc_mzm v_pi={v_pi} il_db={il_db} e_r={e_r} f_c={f_c}
Xm7 a7r a7i a7w b7r b7i b7w d7 gnd fc_mzm v_pi={v_pi} il_db={il_db} e_r={e_r} f_c={f_c}
Xm8 a8r a8i a8w b8r b8i b8w d8 gnd fc_mzm v_pi={v_pi} il_db={il_db} e_r={e_r} f_c={f_c}

* fc_mux: 3N bus wires (output) first, then the N single-channel inputs.
Xmux o0r o0i o0w  o1r o1i o1w  o2r o2i o2w  o3r o3i o3w
+    o4r o4i o4w  o5r o5i o5w  o6r o6i o6w  o7r o7i o7w
+    b1r b1i b1w  b2r b2i b2w  b3r b3i b3w  b4r b4i b4w
+    b5r b5i b5w  b6r b6i b6w  b7r b7i b7w  b8r b8i b8w
+    fc_mux
.ends
