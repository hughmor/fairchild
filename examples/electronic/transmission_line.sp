* Lossless transmission line (T element) — reflection demo.
*
* A 1 ns line (Z0 = 50 Ω) driven through a matched source resistance, with an
* OPEN far end so the wave reflects (+1) and doubles. Watch:
*   V(a): launched half-step (0.5 V) at the 0.5 ns source edge, then steps to
*         1.0 V when the reflection returns at 2·TD = 2.5 ns.
*   V(b): nothing until the wave arrives at TD = 1.5 ns, then jumps to 1.0 V
*         (open-end doubling), settling at the source value.
*
* Syntax:  T<name> A+ A- B+ B- Z0=<ohms> TD=<seconds>
*          (or F=<Hz> [NL=<wavelengths>], TD = NL/F, default NL = 0.25)
*
* The delay is intrinsic — always modelled in transient (Branin's method).
* Compare against ngspice with the same netlist; they agree at every plateau.
*
* Run:
*   fairchild -f transmission_line.sp --probe "V(a),V(b)"

Vs s 0 PULSE(0 1 0.5n 10p 10p 100n 200n)
Rs s a 50
T1 a 0 b 0 Z0=50 TD=1n
* Open far end (large resistor); change to "RL b 0 50" for a matched, clean
* delayed step, or "RL b 0 1m" for an inverted (-1) reflection.
RL b 0 1e9

.tran 20p 4n
.end
