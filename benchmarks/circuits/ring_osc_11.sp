* 11-stage CMOS ring oscillator — larger scaling point
* Tests solver performance on a 24-node nonlinear transient.
.model nm NMOS (vto=0.5 kp=200u lambda=0.05)
.model pm PMOS (vto=-0.5 kp=80u lambda=0.05)
Vdd  vdd 0   DC 1.8

Mn1  n1  n11 0   0   nm w=10u l=1u
Mp1  n1  n11 vdd vdd pm w=20u l=1u
C1   n1  0   100f

Mn2  n2  n1  0   0   nm w=10u l=1u
Mp2  n2  n1  vdd vdd pm w=20u l=1u
C2   n2  0   100f

Mn3  n3  n2  0   0   nm w=10u l=1u
Mp3  n3  n2  vdd vdd pm w=20u l=1u
C3   n3  0   100f

Mn4  n4  n3  0   0   nm w=10u l=1u
Mp4  n4  n3  vdd vdd pm w=20u l=1u
C4   n4  0   100f

Mn5  n5  n4  0   0   nm w=10u l=1u
Mp5  n5  n4  vdd vdd pm w=20u l=1u
C5   n5  0   100f

Mn6  n6  n5  0   0   nm w=10u l=1u
Mp6  n6  n5  vdd vdd pm w=20u l=1u
C6   n6  0   100f

Mn7  n7  n6  0   0   nm w=10u l=1u
Mp7  n7  n6  vdd vdd pm w=20u l=1u
C7   n7  0   100f

Mn8  n8  n7  0   0   nm w=10u l=1u
Mp8  n8  n7  vdd vdd pm w=20u l=1u
C8   n8  0   100f

Mn9  n9  n8  0   0   nm w=10u l=1u
Mp9  n9  n8  vdd vdd pm w=20u l=1u
C9   n9  0   100f

Mn10 n10 n9  0   0   nm w=10u l=1u
Mp10 n10 n9  vdd vdd pm w=20u l=1u
C10  n10 0   100f

Mn11 n11 n10 0   0   nm w=10u l=1u
Mp11 n11 n10 vdd vdd pm w=20u l=1u
C11  n11 0   100f

.ic V(n1)=1.6 V(n2)=0.1 V(n3)=1.6 V(n4)=0.1 V(n5)=1.6
.ic V(n6)=0.1 V(n7)=1.6 V(n8)=0.1 V(n9)=1.6 V(n10)=0.1 V(n11)=1.6
.options method=gear
.tran 50p 25n UIC
.end
