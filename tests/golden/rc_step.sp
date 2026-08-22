* RC step response: τ = RC = 1kΩ × 1µF = 1ms
* V(out) = 1 - e^{-t/τ}  →  0.6321 at t=τ, 0.8647 at t=2τ, 0.9933 at t=5τ
V1 in 0 PULSE(0 1 0 1n 1n 10m 20m)
R1 in out 1k
C1 out 0 1u
.tran 1u 5m
.meas tran v_1tau FIND v(out) AT=1e-3
.meas tran v_2tau FIND v(out) AT=2e-3
.meas tran v_5tau FIND v(out) AT=5e-3
