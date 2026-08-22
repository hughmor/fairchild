* RLC series resonator: R=10Ω, L=1mH, C=1µF
* Resonant frequency f0 = 1/(2π√LC) ≈ 5.03 kHz
* Step response shows damped oscillation.
V1  in  0  PULSE(0 1 0 1n 1n 1 2)
R1  in  n1  10
L1  n1  n2  1m
C1  n2  0   1u
.tran 10u 1m
