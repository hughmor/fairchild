* Half-wave rectifier: ideal Shockley diode + 1kΩ load
* V(out) ≈ max(0, V(in) - Vf) where Vf ≈ 0.6 V
.model D1N4148 D (IS=2.52n N=1.752 RS=0.568 CJO=4p)
V1  in  0  PULSE(-2 2 0 100n 100n 400n 1u)
D1  in  out  D1N4148
R1  out 0    1k
.tran 10n 3u
.end
