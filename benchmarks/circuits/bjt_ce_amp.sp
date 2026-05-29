* NPN common-emitter amplifier — switching transient
* Demonstrates BJT Gummel-Poon Level 1 on a realistic circuit.
* VCC=5V, RB=10k, RC=3.3k; PULSE drives base from cutoff to active.
.model npn1 NPN (IS=1e-15 BF=100 BR=1)
VCC  cc  0   DC 5
VIN  in  0   PULSE(0 0.8 10n 1n 1n 40n 100n)
RB   in  b   10k
RC   cc  c   3.3k
Q1   c b 0 0 npn1
* Step is 1/10 the 1 ns input edge so both simulators resolve the transition.
* A 1 ns step on a 1 ns edge under-resolves it (fairchild fixed-step vs ngspice
* adaptive), which inflates pointwise RMS even though DC levels match exactly.
.tran 0.1n 200n
.end
