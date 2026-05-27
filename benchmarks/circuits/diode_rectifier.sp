* Half-wave diode rectifier — nonlinear transient baseline
* f = 1 MHz; RC load τ = 10µs; peak output ≈ Vpk - 0.7V.
.model dmod D (IS=1e-14 N=1)
V1   in  0   SIN(0 2 1MEG)
D1   in  out dmod
R1   out 0   1k
C1   out 0   10n
.ic V(out)=0
.tran 3n 3u
.end
