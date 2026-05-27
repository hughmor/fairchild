* NPN common-emitter amplifier — BJT Gummel-Poon Level 1
*
* VCC=5V, RB=10kΩ, RC=3.3kΩ.  A PULSE on the base drives the transistor
* from cutoff (Vin=0) through active region to saturation (Vin=0.8V).
*
* Key nodes:
*   b — base (input through RB)
*   c — collector (output; swings 0V..VCC)
*
* Try:
*   fairchild -f examples/electronic/bjt_ce_amplifier.sp \
*             --probe "V(c),V(b)" --opt method=gear

.model npn1 NPN (IS=1e-15 BF=100 BR=1)

VCC  cc  0   DC 5
VIN  in  0   PULSE(0 0.8 10n 1n 1n 40n 100n)
RB   in  b   10k
RC   cc  c   3.3k
Q1   c b 0 0 npn1

.tran 1n 200n
.end
