* 3-section resistor ladder: V1=5V, all R=1k
* V(n1) = 5 * (R2||...) / (R1 + R2||...)  -- solved numerically by both sims
V1 in 0 DC 5.0
R1 in n1 1k
R2 n1 0 1k
R3 n1 n2 1k
R4 n2 0 1k
R5 n2 n3 1k
R6 n3 0 1k

.op
.end
