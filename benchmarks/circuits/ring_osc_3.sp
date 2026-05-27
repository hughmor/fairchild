* 3-stage CMOS ring oscillator — small scaling point
* f ≈ 1/(2·N·t_pd) = 1/(6·t_pd); t_pd set by C_load/I_drive.
.model nm NMOS (vto=0.5 kp=200u lambda=0.05)
.model pm PMOS (vto=-0.5 kp=80u lambda=0.05)
Vdd  vdd 0   DC 1.8

Mn1  n1 n3 0   0   nm w=10u l=1u
Mp1  n1 n3 vdd vdd pm w=20u l=1u
C1   n1 0  100f

Mn2  n2 n1 0   0   nm w=10u l=1u
Mp2  n2 n1 vdd vdd pm w=20u l=1u
C2   n2 0  100f

Mn3  n3 n2 0   0   nm w=10u l=1u
Mp3  n3 n2 vdd vdd pm w=20u l=1u
C3   n3 0  100f

.ic V(n1)=1.6 V(n2)=0.1 V(n3)=1.6
.options method=gear
.tran 50p 10n UIC
.end
