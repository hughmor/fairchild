* RC step response: 1kΩ / 1µF, τ = 1 ms
* V(out) charges to 1 V with time constant 1 ms.
V1  in  0  PULSE(0 1 0 1n 1n 10m 20m)
R1  in  out  1k
C1  out 0    1u
.tran 50u 5m
