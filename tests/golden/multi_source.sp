* Two voltage sources with shared ground, current source load
* V(a) = 3V (from V1), V(b) = 7V (from V2)
* I1 draws 1mA from node a; R1 connects a to b
* I(R1) = (V(b) - V(a)) / R1 = (7-3)/2k = 2mA (from a to b)
V1 a 0 DC 3.0
V2 b 0 DC 7.0
R1 a b 2k
I1 a 0 DC 1m

.op
.end
