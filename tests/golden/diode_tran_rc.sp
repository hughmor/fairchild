* Diode RC transient: half-wave rectifier
* V1 square wave → R1 → D1 → C1 to ground
* C1 charges to ~0.63V (Vfwd) on positive half; holds charge on negative half.
* Expected: V(cap) > 0.55 after first positive half-cycle (t ≈ 0.5 ms)
V1 in 0 PULSE(0 5 0 1n 1n 500u 1m)
R1 in a 1k
D1 a cap myd
C1 cap 0 100n
.model myd D (Is=1e-14 N=1)
.tran 10n 600u
