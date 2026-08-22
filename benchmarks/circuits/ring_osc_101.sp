* 101-stage CMOS ring oscillator
* f ≈ 1/(2·N·t_pd); t_pd set by C_load/I_drive.
.model nm NMOS (vto=0.5 kp=200u lambda=0.05)
.model pm PMOS (vto=-0.5 kp=80u lambda=0.05)
Vdd  vdd 0   DC 1.8

Mn1  n1 n101 0   0   nm w=10u l=1u
Mp1  n1 n101 vdd vdd pm w=20u l=1u
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

Mn10  n10 n9 0   0   nm w=10u l=1u
Mp10  n10 n9 vdd vdd pm w=20u l=1u
C10   n10 0  100f

Mn11  n11 n10 0   0   nm w=10u l=1u
Mp11  n11 n10 vdd vdd pm w=20u l=1u
C11   n11 0  100f

Mn12  n12 n11 0   0   nm w=10u l=1u
Mp12  n12 n11 vdd vdd pm w=20u l=1u
C12   n12 0  100f

Mn13  n13 n12 0   0   nm w=10u l=1u
Mp13  n13 n12 vdd vdd pm w=20u l=1u
C13   n13 0  100f

Mn14  n14 n13 0   0   nm w=10u l=1u
Mp14  n14 n13 vdd vdd pm w=20u l=1u
C14   n14 0  100f

Mn15  n15 n14 0   0   nm w=10u l=1u
Mp15  n15 n14 vdd vdd pm w=20u l=1u
C15   n15 0  100f

Mn16  n16 n15 0   0   nm w=10u l=1u
Mp16  n16 n15 vdd vdd pm w=20u l=1u
C16   n16 0  100f

Mn17  n17 n16 0   0   nm w=10u l=1u
Mp17  n17 n16 vdd vdd pm w=20u l=1u
C17   n17 0  100f

Mn18  n18 n17 0   0   nm w=10u l=1u
Mp18  n18 n17 vdd vdd pm w=20u l=1u
C18   n18 0  100f

Mn19  n19 n18 0   0   nm w=10u l=1u
Mp19  n19 n18 vdd vdd pm w=20u l=1u
C19   n19 0  100f

Mn20  n20 n19 0   0   nm w=10u l=1u
Mp20  n20 n19 vdd vdd pm w=20u l=1u
C20   n20 0  100f

Mn21  n21 n20 0   0   nm w=10u l=1u
Mp21  n21 n20 vdd vdd pm w=20u l=1u
C21   n21 0  100f

Mn22  n22 n21 0   0   nm w=10u l=1u
Mp22  n22 n21 vdd vdd pm w=20u l=1u
C22   n22 0  100f

Mn23  n23 n22 0   0   nm w=10u l=1u
Mp23  n23 n22 vdd vdd pm w=20u l=1u
C23   n23 0  100f

Mn24  n24 n23 0   0   nm w=10u l=1u
Mp24  n24 n23 vdd vdd pm w=20u l=1u
C24   n24 0  100f

Mn25  n25 n24 0   0   nm w=10u l=1u
Mp25  n25 n24 vdd vdd pm w=20u l=1u
C25   n25 0  100f

Mn26  n26 n25 0   0   nm w=10u l=1u
Mp26  n26 n25 vdd vdd pm w=20u l=1u
C26   n26 0  100f

Mn27  n27 n26 0   0   nm w=10u l=1u
Mp27  n27 n26 vdd vdd pm w=20u l=1u
C27   n27 0  100f

Mn28  n28 n27 0   0   nm w=10u l=1u
Mp28  n28 n27 vdd vdd pm w=20u l=1u
C28   n28 0  100f

Mn29  n29 n28 0   0   nm w=10u l=1u
Mp29  n29 n28 vdd vdd pm w=20u l=1u
C29   n29 0  100f

Mn30  n30 n29 0   0   nm w=10u l=1u
Mp30  n30 n29 vdd vdd pm w=20u l=1u
C30   n30 0  100f

Mn31  n31 n30 0   0   nm w=10u l=1u
Mp31  n31 n30 vdd vdd pm w=20u l=1u
C31   n31 0  100f

Mn32  n32 n31 0   0   nm w=10u l=1u
Mp32  n32 n31 vdd vdd pm w=20u l=1u
C32   n32 0  100f

Mn33  n33 n32 0   0   nm w=10u l=1u
Mp33  n33 n32 vdd vdd pm w=20u l=1u
C33   n33 0  100f

Mn34  n34 n33 0   0   nm w=10u l=1u
Mp34  n34 n33 vdd vdd pm w=20u l=1u
C34   n34 0  100f

Mn35  n35 n34 0   0   nm w=10u l=1u
Mp35  n35 n34 vdd vdd pm w=20u l=1u
C35   n35 0  100f

Mn36  n36 n35 0   0   nm w=10u l=1u
Mp36  n36 n35 vdd vdd pm w=20u l=1u
C36   n36 0  100f

Mn37  n37 n36 0   0   nm w=10u l=1u
Mp37  n37 n36 vdd vdd pm w=20u l=1u
C37   n37 0  100f

Mn38  n38 n37 0   0   nm w=10u l=1u
Mp38  n38 n37 vdd vdd pm w=20u l=1u
C38   n38 0  100f

Mn39  n39 n38 0   0   nm w=10u l=1u
Mp39  n39 n38 vdd vdd pm w=20u l=1u
C39   n39 0  100f

Mn40  n40 n39 0   0   nm w=10u l=1u
Mp40  n40 n39 vdd vdd pm w=20u l=1u
C40   n40 0  100f

Mn41  n41 n40 0   0   nm w=10u l=1u
Mp41  n41 n40 vdd vdd pm w=20u l=1u
C41   n41 0  100f

Mn42  n42 n41 0   0   nm w=10u l=1u
Mp42  n42 n41 vdd vdd pm w=20u l=1u
C42   n42 0  100f

Mn43  n43 n42 0   0   nm w=10u l=1u
Mp43  n43 n42 vdd vdd pm w=20u l=1u
C43   n43 0  100f

Mn44  n44 n43 0   0   nm w=10u l=1u
Mp44  n44 n43 vdd vdd pm w=20u l=1u
C44   n44 0  100f

Mn45  n45 n44 0   0   nm w=10u l=1u
Mp45  n45 n44 vdd vdd pm w=20u l=1u
C45   n45 0  100f

Mn46  n46 n45 0   0   nm w=10u l=1u
Mp46  n46 n45 vdd vdd pm w=20u l=1u
C46   n46 0  100f

Mn47  n47 n46 0   0   nm w=10u l=1u
Mp47  n47 n46 vdd vdd pm w=20u l=1u
C47   n47 0  100f

Mn48  n48 n47 0   0   nm w=10u l=1u
Mp48  n48 n47 vdd vdd pm w=20u l=1u
C48   n48 0  100f

Mn49  n49 n48 0   0   nm w=10u l=1u
Mp49  n49 n48 vdd vdd pm w=20u l=1u
C49   n49 0  100f

Mn50  n50 n49 0   0   nm w=10u l=1u
Mp50  n50 n49 vdd vdd pm w=20u l=1u
C50   n50 0  100f

Mn51  n51 n50 0   0   nm w=10u l=1u
Mp51  n51 n50 vdd vdd pm w=20u l=1u
C51   n51 0  100f

Mn52  n52 n51 0   0   nm w=10u l=1u
Mp52  n52 n51 vdd vdd pm w=20u l=1u
C52   n52 0  100f

Mn53  n53 n52 0   0   nm w=10u l=1u
Mp53  n53 n52 vdd vdd pm w=20u l=1u
C53   n53 0  100f

Mn54  n54 n53 0   0   nm w=10u l=1u
Mp54  n54 n53 vdd vdd pm w=20u l=1u
C54   n54 0  100f

Mn55  n55 n54 0   0   nm w=10u l=1u
Mp55  n55 n54 vdd vdd pm w=20u l=1u
C55   n55 0  100f

Mn56  n56 n55 0   0   nm w=10u l=1u
Mp56  n56 n55 vdd vdd pm w=20u l=1u
C56   n56 0  100f

Mn57  n57 n56 0   0   nm w=10u l=1u
Mp57  n57 n56 vdd vdd pm w=20u l=1u
C57   n57 0  100f

Mn58  n58 n57 0   0   nm w=10u l=1u
Mp58  n58 n57 vdd vdd pm w=20u l=1u
C58   n58 0  100f

Mn59  n59 n58 0   0   nm w=10u l=1u
Mp59  n59 n58 vdd vdd pm w=20u l=1u
C59   n59 0  100f

Mn60  n60 n59 0   0   nm w=10u l=1u
Mp60  n60 n59 vdd vdd pm w=20u l=1u
C60   n60 0  100f

Mn61  n61 n60 0   0   nm w=10u l=1u
Mp61  n61 n60 vdd vdd pm w=20u l=1u
C61   n61 0  100f

Mn62  n62 n61 0   0   nm w=10u l=1u
Mp62  n62 n61 vdd vdd pm w=20u l=1u
C62   n62 0  100f

Mn63  n63 n62 0   0   nm w=10u l=1u
Mp63  n63 n62 vdd vdd pm w=20u l=1u
C63   n63 0  100f

Mn64  n64 n63 0   0   nm w=10u l=1u
Mp64  n64 n63 vdd vdd pm w=20u l=1u
C64   n64 0  100f

Mn65  n65 n64 0   0   nm w=10u l=1u
Mp65  n65 n64 vdd vdd pm w=20u l=1u
C65   n65 0  100f

Mn66  n66 n65 0   0   nm w=10u l=1u
Mp66  n66 n65 vdd vdd pm w=20u l=1u
C66   n66 0  100f

Mn67  n67 n66 0   0   nm w=10u l=1u
Mp67  n67 n66 vdd vdd pm w=20u l=1u
C67   n67 0  100f

Mn68  n68 n67 0   0   nm w=10u l=1u
Mp68  n68 n67 vdd vdd pm w=20u l=1u
C68   n68 0  100f

Mn69  n69 n68 0   0   nm w=10u l=1u
Mp69  n69 n68 vdd vdd pm w=20u l=1u
C69   n69 0  100f

Mn70  n70 n69 0   0   nm w=10u l=1u
Mp70  n70 n69 vdd vdd pm w=20u l=1u
C70   n70 0  100f

Mn71  n71 n70 0   0   nm w=10u l=1u
Mp71  n71 n70 vdd vdd pm w=20u l=1u
C71   n71 0  100f

Mn72  n72 n71 0   0   nm w=10u l=1u
Mp72  n72 n71 vdd vdd pm w=20u l=1u
C72   n72 0  100f

Mn73  n73 n72 0   0   nm w=10u l=1u
Mp73  n73 n72 vdd vdd pm w=20u l=1u
C73   n73 0  100f

Mn74  n74 n73 0   0   nm w=10u l=1u
Mp74  n74 n73 vdd vdd pm w=20u l=1u
C74   n74 0  100f

Mn75  n75 n74 0   0   nm w=10u l=1u
Mp75  n75 n74 vdd vdd pm w=20u l=1u
C75   n75 0  100f

Mn76  n76 n75 0   0   nm w=10u l=1u
Mp76  n76 n75 vdd vdd pm w=20u l=1u
C76   n76 0  100f

Mn77  n77 n76 0   0   nm w=10u l=1u
Mp77  n77 n76 vdd vdd pm w=20u l=1u
C77   n77 0  100f

Mn78  n78 n77 0   0   nm w=10u l=1u
Mp78  n78 n77 vdd vdd pm w=20u l=1u
C78   n78 0  100f

Mn79  n79 n78 0   0   nm w=10u l=1u
Mp79  n79 n78 vdd vdd pm w=20u l=1u
C79   n79 0  100f

Mn80  n80 n79 0   0   nm w=10u l=1u
Mp80  n80 n79 vdd vdd pm w=20u l=1u
C80   n80 0  100f

Mn81  n81 n80 0   0   nm w=10u l=1u
Mp81  n81 n80 vdd vdd pm w=20u l=1u
C81   n81 0  100f

Mn82  n82 n81 0   0   nm w=10u l=1u
Mp82  n82 n81 vdd vdd pm w=20u l=1u
C82   n82 0  100f

Mn83  n83 n82 0   0   nm w=10u l=1u
Mp83  n83 n82 vdd vdd pm w=20u l=1u
C83   n83 0  100f

Mn84  n84 n83 0   0   nm w=10u l=1u
Mp84  n84 n83 vdd vdd pm w=20u l=1u
C84   n84 0  100f

Mn85  n85 n84 0   0   nm w=10u l=1u
Mp85  n85 n84 vdd vdd pm w=20u l=1u
C85   n85 0  100f

Mn86  n86 n85 0   0   nm w=10u l=1u
Mp86  n86 n85 vdd vdd pm w=20u l=1u
C86   n86 0  100f

Mn87  n87 n86 0   0   nm w=10u l=1u
Mp87  n87 n86 vdd vdd pm w=20u l=1u
C87   n87 0  100f

Mn88  n88 n87 0   0   nm w=10u l=1u
Mp88  n88 n87 vdd vdd pm w=20u l=1u
C88   n88 0  100f

Mn89  n89 n88 0   0   nm w=10u l=1u
Mp89  n89 n88 vdd vdd pm w=20u l=1u
C89   n89 0  100f

Mn90  n90 n89 0   0   nm w=10u l=1u
Mp90  n90 n89 vdd vdd pm w=20u l=1u
C90   n90 0  100f

Mn91  n91 n90 0   0   nm w=10u l=1u
Mp91  n91 n90 vdd vdd pm w=20u l=1u
C91   n91 0  100f

Mn92  n92 n91 0   0   nm w=10u l=1u
Mp92  n92 n91 vdd vdd pm w=20u l=1u
C92   n92 0  100f

Mn93  n93 n92 0   0   nm w=10u l=1u
Mp93  n93 n92 vdd vdd pm w=20u l=1u
C93   n93 0  100f

Mn94  n94 n93 0   0   nm w=10u l=1u
Mp94  n94 n93 vdd vdd pm w=20u l=1u
C94   n94 0  100f

Mn95  n95 n94 0   0   nm w=10u l=1u
Mp95  n95 n94 vdd vdd pm w=20u l=1u
C95   n95 0  100f

Mn96  n96 n95 0   0   nm w=10u l=1u
Mp96  n96 n95 vdd vdd pm w=20u l=1u
C96   n96 0  100f

Mn97  n97 n96 0   0   nm w=10u l=1u
Mp97  n97 n96 vdd vdd pm w=20u l=1u
C97   n97 0  100f

Mn98  n98 n97 0   0   nm w=10u l=1u
Mp98  n98 n97 vdd vdd pm w=20u l=1u
C98   n98 0  100f

Mn99  n99 n98 0   0   nm w=10u l=1u
Mp99  n99 n98 vdd vdd pm w=20u l=1u
C99   n99 0  100f

Mn100  n100 n99 0   0   nm w=10u l=1u
Mp100  n100 n99 vdd vdd pm w=20u l=1u
C100   n100 0  100f

Mn101  n101 n100 0   0   nm w=10u l=1u
Mp101  n101 n100 vdd vdd pm w=20u l=1u
C101   n101 0  100f

.ic V(n1)=1.6 V(n2)=0.1 V(n3)=1.6 V(n4)=0.1 V(n5)=1.6 V(n6)=0.1 V(n7)=1.6 V(n8)=0.1 V(n9)=1.6 V(n10)=0.1 V(n11)=1.6 V(n12)=0.1 V(n13)=1.6 V(n14)=0.1 V(n15)=1.6 V(n16)=0.1 V(n17)=1.6 V(n18)=0.1 V(n19)=1.6 V(n20)=0.1 V(n21)=1.6 V(n22)=0.1 V(n23)=1.6 V(n24)=0.1 V(n25)=1.6 V(n26)=0.1 V(n27)=1.6 V(n28)=0.1 V(n29)=1.6 V(n30)=0.1 V(n31)=1.6 V(n32)=0.1 V(n33)=1.6 V(n34)=0.1 V(n35)=1.6 V(n36)=0.1 V(n37)=1.6 V(n38)=0.1 V(n39)=1.6 V(n40)=0.1 V(n41)=1.6 V(n42)=0.1 V(n43)=1.6 V(n44)=0.1 V(n45)=1.6 V(n46)=0.1 V(n47)=1.6 V(n48)=0.1 V(n49)=1.6 V(n50)=0.1 V(n51)=1.6 V(n52)=0.1 V(n53)=1.6 V(n54)=0.1 V(n55)=1.6 V(n56)=0.1 V(n57)=1.6 V(n58)=0.1 V(n59)=1.6 V(n60)=0.1 V(n61)=1.6 V(n62)=0.1 V(n63)=1.6 V(n64)=0.1 V(n65)=1.6 V(n66)=0.1 V(n67)=1.6 V(n68)=0.1 V(n69)=1.6 V(n70)=0.1 V(n71)=1.6 V(n72)=0.1 V(n73)=1.6 V(n74)=0.1 V(n75)=1.6 V(n76)=0.1 V(n77)=1.6 V(n78)=0.1 V(n79)=1.6 V(n80)=0.1 V(n81)=1.6 V(n82)=0.1 V(n83)=1.6 V(n84)=0.1 V(n85)=1.6 V(n86)=0.1 V(n87)=1.6 V(n88)=0.1 V(n89)=1.6 V(n90)=0.1 V(n91)=1.6 V(n92)=0.1 V(n93)=1.6 V(n94)=0.1 V(n95)=1.6 V(n96)=0.1 V(n97)=1.6 V(n98)=0.1 V(n99)=1.6 V(n100)=0.1 V(n101)=1.6
.options method=gear
.tran 50p 152n UIC
