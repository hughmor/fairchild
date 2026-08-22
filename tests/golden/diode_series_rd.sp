* Diode DC: R-D series with Vdd
* Vdd=5V, R1=10k, D1 from b to ground
* Requires NR iteration; V(b) ≈ 0.63 V
Vdd a 0 DC 5
R1 a b 10k
D1 b 0 myd
.model myd D (Is=1e-14 N=1)
.op
