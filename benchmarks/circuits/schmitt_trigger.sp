* CMOS Schmitt trigger — 4-transistor inverting topology
* Demonstrates hysteretic switching: V_T+(high→low) ≠ V_T-(low→high) relative to VDD/2
* Tests convergence on positive-feedback nonlinear circuit.
* Realistic parasitic caps (overlap + junction) give the internal nodes real
* capacitance. Without them, the held-LOW output is a zero-capacitance floating
* island (MN2 feedback gate cuts off near VTO, MP1 off) — an ill-posed node that
* fairchild's diagonal gmin discharges to 0 V while ngspice holds it at ~VTO.
* With caps the node holds charge and both settle at ~0.5 V (= VTO). See
* docs/benchmarks.md "Schmitt trigger held-low state".
.model nm NMOS (VTO=0.5 KP=200u lambda=0.05 CGSO=2.5e-10 CGDO=2.5e-10 CJ=2e-4 CJSW=5e-10)
.model pm PMOS (VTO=-0.5 KP=80u lambda=0.05 CGSO=2.5e-10 CGDO=2.5e-10 CJ=2e-4 CJSW=5e-10)
VDD  vdd 0  DC 3.3

* Slow ramp: 0→3.3V then 3.3→0V; passes through both switching thresholds
VIN  in  0  PULSE(0 3.3 100n 600n 600n 1000n 2800n)

* Pull-up stack: MP2 (feedback gate=out) in series with MP1 (main gate=in)
* When out=LOW → MP2 ON → extra pull-up strength → V_T-(0→1) shifted up
MP2  pint out  vdd vdd pm  W=10u L=1u AS=20p AD=20p PS=14u PD=14u
MP1  out  in   pint vdd pm  W=20u L=1u AS=40p AD=40p PS=24u PD=24u

* Pull-down stack: MN1 (main gate=in) in series with MN2 (feedback gate=out)
* When out=HIGH → MN2 ON → extra pull-down strength → V_T+(1→0) shifted down
MN1  out  in   nint 0   nm  W=10u L=1u AS=20p AD=20p PS=14u PD=14u
MN2  nint out  0    0   nm  W=5u  L=1u AS=10p AD=10p PS=9u  PD=9u

.tran 5n 2800n
