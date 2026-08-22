* 51-stage CMOS ring oscillator — scaling benchmark
.model nm NMOS (vto=0.5 kp=200u lambda=0.05)
.model pm PMOS (vto=-0.5 kp=80u lambda=0.05)
Vdd  vdd 0   DC 1.8

Mn1  n1 n51 0   0   nm w=10u l=1u
Mp1  n1 n51 vdd vdd pm w=20u l=1u
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

Mn6  n6 n5 0   0   nm w=10u l=1u
Mp6  n6 n5 vdd vdd pm w=20u l=1u
C6   n6 0  100f

Mn7  n7 n6 0   0   nm w=10u l=1u
Mp7  n7 n6 vdd vdd pm w=20u l=1u
C7   n7 0  100f

Mn8  n8 n7 0   0   nm w=10u l=1u
Mp8  n8 n7 vdd vdd pm w=20u l=1u
C8   n8 0  100f

Mn9  n9 n8 0   0   nm w=10u l=1u
Mp9  n9 n8 vdd vdd pm w=20u l=1u
C9   n9 0  100f

Mn10 n10 n9 0   0   nm w=10u l=1u
Mp10 n10 n9 vdd vdd pm w=20u l=1u
C10  n10 0  100f

Mn11 n11 n10 0   0   nm w=10u l=1u
Mp11 n11 n10 vdd vdd pm w=20u l=1u
C11  n11 0  100f

Mn12 n12 n11 0   0   nm w=10u l=1u
Mp12 n12 n11 vdd vdd pm w=20u l=1u
C12  n12 0  100f

Mn13 n13 n12 0   0   nm w=10u l=1u
Mp13 n13 n12 vdd vdd pm w=20u l=1u
C13  n13 0  100f

Mn14 n14 n13 0   0   nm w=10u l=1u
Mp14 n14 n13 vdd vdd pm w=20u l=1u
C14  n14 0  100f

Mn15 n15 n14 0   0   nm w=10u l=1u
Mp15 n15 n14 vdd vdd pm w=20u l=1u
C15  n15 0  100f

Mn16 n16 n15 0   0   nm w=10u l=1u
Mp16 n16 n15 vdd vdd pm w=20u l=1u
C16  n16 0  100f

Mn17 n17 n16 0   0   nm w=10u l=1u
Mp17 n17 n16 vdd vdd pm w=20u l=1u
C17  n17 0  100f

Mn18 n18 n17 0   0   nm w=10u l=1u
Mp18 n18 n17 vdd vdd pm w=20u l=1u
C18  n18 0  100f

Mn19 n19 n18 0   0   nm w=10u l=1u
Mp19 n19 n18 vdd vdd pm w=20u l=1u
C19  n19 0  100f

Mn20 n20 n19 0   0   nm w=10u l=1u
Mp20 n20 n19 vdd vdd pm w=20u l=1u
C20  n20 0  100f

Mn21 n21 n20 0   0   nm w=10u l=1u
Mp21 n21 n20 vdd vdd pm w=20u l=1u
C21  n21 0  100f

Mn22 n22 n21 0   0   nm w=10u l=1u
Mp22 n22 n21 vdd vdd pm w=20u l=1u
C22  n22 0  100f

Mn23 n23 n22 0   0   nm w=10u l=1u
Mp23 n23 n22 vdd vdd pm w=20u l=1u
C23  n23 0  100f

Mn24 n24 n23 0   0   nm w=10u l=1u
Mp24 n24 n23 vdd vdd pm w=20u l=1u
C24  n24 0  100f

Mn25 n25 n24 0   0   nm w=10u l=1u
Mp25 n25 n24 vdd vdd pm w=20u l=1u
C25  n25 0  100f

Mn26 n26 n25 0   0   nm w=10u l=1u
Mp26 n26 n25 vdd vdd pm w=20u l=1u
C26  n26 0  100f

Mn27 n27 n26 0   0   nm w=10u l=1u
Mp27 n27 n26 vdd vdd pm w=20u l=1u
C27  n27 0  100f

Mn28 n28 n27 0   0   nm w=10u l=1u
Mp28 n28 n27 vdd vdd pm w=20u l=1u
C28  n28 0  100f

Mn29 n29 n28 0   0   nm w=10u l=1u
Mp29 n29 n28 vdd vdd pm w=20u l=1u
C29  n29 0  100f

Mn30 n30 n29 0   0   nm w=10u l=1u
Mp30 n30 n29 vdd vdd pm w=20u l=1u
C30  n30 0  100f

Mn31 n31 n30 0   0   nm w=10u l=1u
Mp31 n31 n30 vdd vdd pm w=20u l=1u
C31  n31 0  100f

Mn32 n32 n31 0   0   nm w=10u l=1u
Mp32 n32 n31 vdd vdd pm w=20u l=1u
C32  n32 0  100f

Mn33 n33 n32 0   0   nm w=10u l=1u
Mp33 n33 n32 vdd vdd pm w=20u l=1u
C33  n33 0  100f

Mn34 n34 n33 0   0   nm w=10u l=1u
Mp34 n34 n33 vdd vdd pm w=20u l=1u
C34  n34 0  100f

Mn35 n35 n34 0   0   nm w=10u l=1u
Mp35 n35 n34 vdd vdd pm w=20u l=1u
C35  n35 0  100f

Mn36 n36 n35 0   0   nm w=10u l=1u
Mp36 n36 n35 vdd vdd pm w=20u l=1u
C36  n36 0  100f

Mn37 n37 n36 0   0   nm w=10u l=1u
Mp37 n37 n36 vdd vdd pm w=20u l=1u
C37  n37 0  100f

Mn38 n38 n37 0   0   nm w=10u l=1u
Mp38 n38 n37 vdd vdd pm w=20u l=1u
C38  n38 0  100f

Mn39 n39 n38 0   0   nm w=10u l=1u
Mp39 n39 n38 vdd vdd pm w=20u l=1u
C39  n39 0  100f

Mn40 n40 n39 0   0   nm w=10u l=1u
Mp40 n40 n39 vdd vdd pm w=20u l=1u
C40  n40 0  100f

Mn41 n41 n40 0   0   nm w=10u l=1u
Mp41 n41 n40 vdd vdd pm w=20u l=1u
C41  n41 0  100f

Mn42 n42 n41 0   0   nm w=10u l=1u
Mp42 n42 n41 vdd vdd pm w=20u l=1u
C42  n42 0  100f

Mn43 n43 n42 0   0   nm w=10u l=1u
Mp43 n43 n42 vdd vdd pm w=20u l=1u
C43  n43 0  100f

Mn44 n44 n43 0   0   nm w=10u l=1u
Mp44 n44 n43 vdd vdd pm w=20u l=1u
C44  n44 0  100f

Mn45 n45 n44 0   0   nm w=10u l=1u
Mp45 n45 n44 vdd vdd pm w=20u l=1u
C45  n45 0  100f

Mn46 n46 n45 0   0   nm w=10u l=1u
Mp46 n46 n45 vdd vdd pm w=20u l=1u
C46  n46 0  100f

Mn47 n47 n46 0   0   nm w=10u l=1u
Mp47 n47 n46 vdd vdd pm w=20u l=1u
C47  n47 0  100f

Mn48 n48 n47 0   0   nm w=10u l=1u
Mp48 n48 n47 vdd vdd pm w=20u l=1u
C48  n48 0  100f

Mn49 n49 n48 0   0   nm w=10u l=1u
Mp49 n49 n48 vdd vdd pm w=20u l=1u
C49  n49 0  100f

Mn50 n50 n49 0   0   nm w=10u l=1u
Mp50 n50 n49 vdd vdd pm w=20u l=1u
C50  n50 0  100f

Mn51 n51 n50 0   0   nm w=10u l=1u
Mp51 n51 n50 vdd vdd pm w=20u l=1u
C51  n51 0  100f

.ic V(n1)=1.6 V(n2)=0.1 V(n3)=1.6 V(n4)=0.1 V(n5)=1.6 V(n6)=0.1 V(n7)=1.6 V(n8)=0.1 V(n9)=1.6 V(n10)=0.1 V(n11)=1.6 V(n12)=0.1 V(n13)=1.6 V(n14)=0.1 V(n15)=1.6 V(n16)=0.1 V(n17)=1.6 V(n18)=0.1 V(n19)=1.6 V(n20)=0.1 V(n21)=1.6 V(n22)=0.1 V(n23)=1.6 V(n24)=0.1 V(n25)=1.6 V(n26)=0.1 V(n27)=1.6 V(n28)=0.1 V(n29)=1.6 V(n30)=0.1 V(n31)=1.6 V(n32)=0.1 V(n33)=1.6 V(n34)=0.1 V(n35)=1.6 V(n36)=0.1 V(n37)=1.6 V(n38)=0.1 V(n39)=1.6 V(n40)=0.1 V(n41)=1.6 V(n42)=0.1 V(n43)=1.6 V(n44)=0.1 V(n45)=1.6 V(n46)=0.1 V(n47)=1.6 V(n48)=0.1 V(n49)=1.6 V(n50)=0.1 V(n51)=1.6
.options method=gear
.tran 50p 200n UIC
