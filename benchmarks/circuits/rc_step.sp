* RC step response — simple linear transient baseline
* τ = R*C = 1kΩ * 1µF = 1ms; step to 1V, simulate 5τ.
V1  in  0  PULSE(0 1 0 1n 1n 10m 20m)
R1  in  out 1k
C1  out 0   1u
.ic V(out)=0
.tran 10u 5m
