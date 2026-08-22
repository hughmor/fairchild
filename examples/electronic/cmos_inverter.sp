* CMOS inverter: NMOS + PMOS Level 1 (Shichman-Hodges)
* VIN sweeps from 0 → VDD; VOUT inverts from VDD → 0
.model nm NMOS (VTO=0.7  KP=100u)
.model pm PMOS (VTO=-0.7 KP=100u)
VDD vdd 0  DC 3.3
VIN in  0  PULSE(0 3.3 10n 1n 1n 40n 100n)
MN  out in 0   0   nm  W=10u L=1u
MP  out in vdd vdd pm  W=10u L=1u
.tran 1n 120n
