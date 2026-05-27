* NPN common-emitter amplifier — switching transient
* Demonstrates BJT Gummel-Poon Level 1 on a realistic circuit.
* VCC=5V, RB=10k, RC=3.3k; PULSE drives base from cutoff to active.
.model npn1 NPN (IS=1e-15 BF=100 BR=1)
VCC  cc  0   DC 5
VIN  in  0   PULSE(0 0.8 10n 1n 1n 40n 100n)
RB   in  b   10k
RC   cc  c   3.3k
Q1   c b 0 0 npn1
.tran 1n 200n
.end
