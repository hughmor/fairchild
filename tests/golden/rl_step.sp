* RL step response: τ = L/R = 1H / 1kΩ = 1ms
* V(out) = V across L = e^{-t/τ}  →  0.3679 at t=τ, 0.1353 at t=2τ, 0.0067 at t=5τ
V1 in 0 PULSE(0 1 0 1n 1n 10m 20m)
R1 in out 1k
L1 out 0 1
.tran 1u 5m
.meas tran v_1tau FIND v(out) AT=1e-3
.meas tran v_2tau FIND v(out) AT=2e-3
.meas tran v_5tau FIND v(out) AT=5e-3
.end
