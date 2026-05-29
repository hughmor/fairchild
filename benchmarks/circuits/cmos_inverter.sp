* CMOS inverter switching transient — Level-1 MOSFET with Meyer+junction caps
.model nm NMOS (VTO=0.7 KP=100u CGSO=2.5e-10 CGDO=2.5e-10 CJ=2e-4 CJSW=5e-10)
.model pm PMOS (VTO=-0.7 KP=100u CGSO=2.5e-10 CGDO=2.5e-10 CJ=2e-4 CJSW=5e-10)
VDD  vdd 0   DC 3.3
VIN  in  0   PULSE(0 3.3 10n 1n 1n 40n 100n)
MN   out in  0   0   nm  W=10u L=1u AS=50p AD=50p PS=20u PD=20u
MP   out in  vdd vdd pm  W=10u L=1u AS=50p AD=50p PS=20u PD=20u
* Step is 1/5 the 1 ns input edge so both simulators resolve the transition.
* Residual RMS after this is a ~40 ps Meyer-cap edge-shape difference (inherent
* to the approximate Meyer gate-charge model), not a DC/level error — see
* docs/benchmarks.md.
.tran 0.2n 120n
.end
