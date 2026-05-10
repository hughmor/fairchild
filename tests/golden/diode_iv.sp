* Diode DC I-V: current-source bias
* 1 mA forced through D1; V(b) = Vt*ln(Ib/Is) ≈ 0.655 V
Ib 0 b 1m
D1 b 0 myd
.model myd D (Is=1e-14 N=1)
.op
.end
