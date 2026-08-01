* Switched-capacitor sample-and-hold — voltage-controlled switch (S)
*
* A clock gates a ramp onto a hold capacitor: while V(clk) > VT the switch is
* RON = 10 Ω and C1 tracks the input; while it is below, the switch is
* ROFF = 1 GΩ and C1 holds the last sampled value. The staircase on V(out) is
* the point — each plateau is the input at the previous falling clock edge.
*
* `.model <name> SW (VT= VH= RON= ROFF=)`
*   VT   threshold on V(NC+,NC-)          VH   hysteresis half-width (default 0)
*   RON  on-resistance (Ω)                ROFF off-resistance (Ω)
*
* Switching is a hard step, not a smooth transition — the same model ngspice
* implements. Two practical consequences:
*   - Pick a timestep fine enough that the hold cap's companion conductance
*     (2C/h under TR) dominates 1/RON, or the node can move far enough in one
*     step to re-cross the threshold and stall Newton.
*   - Give VH a non-zero value on any switch whose own output can reach its
*     control input; that feedback path is what a hard switch chatters on.
*
* The current-controlled twin is `W<name> N+ N- <vsource> <model>` with
* `.model <name> CSW (IT= IH= RON= ROFF=)`, watching a source's branch current.

.model swmod SW (VT=2.5 VH=0 RON=10 ROFF=1e9)

Vin  in  0   PWL(0 0 100u 5)
Vclk clk 0   PULSE(0 5 0 1n 1n 10u 20u)
S1   in  out clk 0 swmod OFF
C1   out 0   1n

.tran 0.2u 100u
.end
