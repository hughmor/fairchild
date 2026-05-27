* 5-stage CMOS ring oscillator — scaling benchmark
.model nm NMOS (vto=0.5 kp=200u lambda=0.05)
.model pm PMOS (vto=-0.5 kp=80u lambda=0.05)
Vdd  vdd 0   DC 1.8

Mn1  n1 n5 0   0   nm w=10u l=1u
Mp1  n1 n5 vdd vdd pm w=20u l=1u
C1   n1 0  100f

Mn2  n2 n1 0   0   nm w=10u l=1u
Mp2  n2 n1 vdd vdd pm w=20u l=1u
C2   n2 0  100f

Mn3  n3 n2 0   0   nm w=10u l=1u
Mp3  n3 n2 vdd vdd pm w=20u l=1u
C3   n3 0  100f

Mn4  n4 n3 0   0   nm w=10u l=1u
Mp4  n4 n3 vdd vdd pm w=20u l=1u
C4   n4 0  100f

Mn5  n5 n4 0   0   nm w=10u l=1u
Mp5  n5 n4 vdd vdd pm w=20u l=1u
C5   n5 0  100f

.ic V(n1)=1.6 V(n2)=0.1 V(n3)=1.6 V(n4)=0.1 V(n5)=1.6
.options method=gear
.tran 50p 30n UIC
.end
