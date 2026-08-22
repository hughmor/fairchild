* NMOS resistive-load DC operating point
* Used for DC validation against ngspice
.model nm1 NMOS (VTO=0.7 KP=100u)
VDD vdd 0 DC 3.3
VG  g   0 DC 2.5
R1  vdd d 10k
M1  d g 0 0 nm1 W=10u L=1u
.op
