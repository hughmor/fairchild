* RLC series resonator — tests underdamped transient
* f0 = 1/(2π√(LC)) = 1/(2π√(1m·1µ)) ≈ 5.03 kHz; Q = (1/R)√(L/C) ≈ 10.
V1   in  0   PULSE(0 1 0 1n 1n 0.5m 1m)
R1   in  n1  10
L1   n1  n2  1m
C1   n2  0   1u
.ic V(n2)=0
.tran 5u 1m
.end
