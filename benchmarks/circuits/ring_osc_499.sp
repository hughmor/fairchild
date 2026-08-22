* 499-stage CMOS ring oscillator
* f ≈ 1/(2·N·t_pd); t_pd set by C_load/I_drive.
.model nm NMOS (vto=0.5 kp=200u lambda=0.05)
.model pm PMOS (vto=-0.5 kp=80u lambda=0.05)
Vdd  vdd 0   DC 1.8

Mn1  n1 n499 0   0   nm w=10u l=1u
Mp1  n1 n499 vdd vdd pm w=20u l=1u
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

Mn102  n102 n101 0   0   nm w=10u l=1u
Mp102  n102 n101 vdd vdd pm w=20u l=1u
C102   n102 0  100f

Mn103  n103 n102 0   0   nm w=10u l=1u
Mp103  n103 n102 vdd vdd pm w=20u l=1u
C103   n103 0  100f

Mn104  n104 n103 0   0   nm w=10u l=1u
Mp104  n104 n103 vdd vdd pm w=20u l=1u
C104   n104 0  100f

Mn105  n105 n104 0   0   nm w=10u l=1u
Mp105  n105 n104 vdd vdd pm w=20u l=1u
C105   n105 0  100f

Mn106  n106 n105 0   0   nm w=10u l=1u
Mp106  n106 n105 vdd vdd pm w=20u l=1u
C106   n106 0  100f

Mn107  n107 n106 0   0   nm w=10u l=1u
Mp107  n107 n106 vdd vdd pm w=20u l=1u
C107   n107 0  100f

Mn108  n108 n107 0   0   nm w=10u l=1u
Mp108  n108 n107 vdd vdd pm w=20u l=1u
C108   n108 0  100f

Mn109  n109 n108 0   0   nm w=10u l=1u
Mp109  n109 n108 vdd vdd pm w=20u l=1u
C109   n109 0  100f

Mn110  n110 n109 0   0   nm w=10u l=1u
Mp110  n110 n109 vdd vdd pm w=20u l=1u
C110   n110 0  100f

Mn111  n111 n110 0   0   nm w=10u l=1u
Mp111  n111 n110 vdd vdd pm w=20u l=1u
C111   n111 0  100f

Mn112  n112 n111 0   0   nm w=10u l=1u
Mp112  n112 n111 vdd vdd pm w=20u l=1u
C112   n112 0  100f

Mn113  n113 n112 0   0   nm w=10u l=1u
Mp113  n113 n112 vdd vdd pm w=20u l=1u
C113   n113 0  100f

Mn114  n114 n113 0   0   nm w=10u l=1u
Mp114  n114 n113 vdd vdd pm w=20u l=1u
C114   n114 0  100f

Mn115  n115 n114 0   0   nm w=10u l=1u
Mp115  n115 n114 vdd vdd pm w=20u l=1u
C115   n115 0  100f

Mn116  n116 n115 0   0   nm w=10u l=1u
Mp116  n116 n115 vdd vdd pm w=20u l=1u
C116   n116 0  100f

Mn117  n117 n116 0   0   nm w=10u l=1u
Mp117  n117 n116 vdd vdd pm w=20u l=1u
C117   n117 0  100f

Mn118  n118 n117 0   0   nm w=10u l=1u
Mp118  n118 n117 vdd vdd pm w=20u l=1u
C118   n118 0  100f

Mn119  n119 n118 0   0   nm w=10u l=1u
Mp119  n119 n118 vdd vdd pm w=20u l=1u
C119   n119 0  100f

Mn120  n120 n119 0   0   nm w=10u l=1u
Mp120  n120 n119 vdd vdd pm w=20u l=1u
C120   n120 0  100f

Mn121  n121 n120 0   0   nm w=10u l=1u
Mp121  n121 n120 vdd vdd pm w=20u l=1u
C121   n121 0  100f

Mn122  n122 n121 0   0   nm w=10u l=1u
Mp122  n122 n121 vdd vdd pm w=20u l=1u
C122   n122 0  100f

Mn123  n123 n122 0   0   nm w=10u l=1u
Mp123  n123 n122 vdd vdd pm w=20u l=1u
C123   n123 0  100f

Mn124  n124 n123 0   0   nm w=10u l=1u
Mp124  n124 n123 vdd vdd pm w=20u l=1u
C124   n124 0  100f

Mn125  n125 n124 0   0   nm w=10u l=1u
Mp125  n125 n124 vdd vdd pm w=20u l=1u
C125   n125 0  100f

Mn126  n126 n125 0   0   nm w=10u l=1u
Mp126  n126 n125 vdd vdd pm w=20u l=1u
C126   n126 0  100f

Mn127  n127 n126 0   0   nm w=10u l=1u
Mp127  n127 n126 vdd vdd pm w=20u l=1u
C127   n127 0  100f

Mn128  n128 n127 0   0   nm w=10u l=1u
Mp128  n128 n127 vdd vdd pm w=20u l=1u
C128   n128 0  100f

Mn129  n129 n128 0   0   nm w=10u l=1u
Mp129  n129 n128 vdd vdd pm w=20u l=1u
C129   n129 0  100f

Mn130  n130 n129 0   0   nm w=10u l=1u
Mp130  n130 n129 vdd vdd pm w=20u l=1u
C130   n130 0  100f

Mn131  n131 n130 0   0   nm w=10u l=1u
Mp131  n131 n130 vdd vdd pm w=20u l=1u
C131   n131 0  100f

Mn132  n132 n131 0   0   nm w=10u l=1u
Mp132  n132 n131 vdd vdd pm w=20u l=1u
C132   n132 0  100f

Mn133  n133 n132 0   0   nm w=10u l=1u
Mp133  n133 n132 vdd vdd pm w=20u l=1u
C133   n133 0  100f

Mn134  n134 n133 0   0   nm w=10u l=1u
Mp134  n134 n133 vdd vdd pm w=20u l=1u
C134   n134 0  100f

Mn135  n135 n134 0   0   nm w=10u l=1u
Mp135  n135 n134 vdd vdd pm w=20u l=1u
C135   n135 0  100f

Mn136  n136 n135 0   0   nm w=10u l=1u
Mp136  n136 n135 vdd vdd pm w=20u l=1u
C136   n136 0  100f

Mn137  n137 n136 0   0   nm w=10u l=1u
Mp137  n137 n136 vdd vdd pm w=20u l=1u
C137   n137 0  100f

Mn138  n138 n137 0   0   nm w=10u l=1u
Mp138  n138 n137 vdd vdd pm w=20u l=1u
C138   n138 0  100f

Mn139  n139 n138 0   0   nm w=10u l=1u
Mp139  n139 n138 vdd vdd pm w=20u l=1u
C139   n139 0  100f

Mn140  n140 n139 0   0   nm w=10u l=1u
Mp140  n140 n139 vdd vdd pm w=20u l=1u
C140   n140 0  100f

Mn141  n141 n140 0   0   nm w=10u l=1u
Mp141  n141 n140 vdd vdd pm w=20u l=1u
C141   n141 0  100f

Mn142  n142 n141 0   0   nm w=10u l=1u
Mp142  n142 n141 vdd vdd pm w=20u l=1u
C142   n142 0  100f

Mn143  n143 n142 0   0   nm w=10u l=1u
Mp143  n143 n142 vdd vdd pm w=20u l=1u
C143   n143 0  100f

Mn144  n144 n143 0   0   nm w=10u l=1u
Mp144  n144 n143 vdd vdd pm w=20u l=1u
C144   n144 0  100f

Mn145  n145 n144 0   0   nm w=10u l=1u
Mp145  n145 n144 vdd vdd pm w=20u l=1u
C145   n145 0  100f

Mn146  n146 n145 0   0   nm w=10u l=1u
Mp146  n146 n145 vdd vdd pm w=20u l=1u
C146   n146 0  100f

Mn147  n147 n146 0   0   nm w=10u l=1u
Mp147  n147 n146 vdd vdd pm w=20u l=1u
C147   n147 0  100f

Mn148  n148 n147 0   0   nm w=10u l=1u
Mp148  n148 n147 vdd vdd pm w=20u l=1u
C148   n148 0  100f

Mn149  n149 n148 0   0   nm w=10u l=1u
Mp149  n149 n148 vdd vdd pm w=20u l=1u
C149   n149 0  100f

Mn150  n150 n149 0   0   nm w=10u l=1u
Mp150  n150 n149 vdd vdd pm w=20u l=1u
C150   n150 0  100f

Mn151  n151 n150 0   0   nm w=10u l=1u
Mp151  n151 n150 vdd vdd pm w=20u l=1u
C151   n151 0  100f

Mn152  n152 n151 0   0   nm w=10u l=1u
Mp152  n152 n151 vdd vdd pm w=20u l=1u
C152   n152 0  100f

Mn153  n153 n152 0   0   nm w=10u l=1u
Mp153  n153 n152 vdd vdd pm w=20u l=1u
C153   n153 0  100f

Mn154  n154 n153 0   0   nm w=10u l=1u
Mp154  n154 n153 vdd vdd pm w=20u l=1u
C154   n154 0  100f

Mn155  n155 n154 0   0   nm w=10u l=1u
Mp155  n155 n154 vdd vdd pm w=20u l=1u
C155   n155 0  100f

Mn156  n156 n155 0   0   nm w=10u l=1u
Mp156  n156 n155 vdd vdd pm w=20u l=1u
C156   n156 0  100f

Mn157  n157 n156 0   0   nm w=10u l=1u
Mp157  n157 n156 vdd vdd pm w=20u l=1u
C157   n157 0  100f

Mn158  n158 n157 0   0   nm w=10u l=1u
Mp158  n158 n157 vdd vdd pm w=20u l=1u
C158   n158 0  100f

Mn159  n159 n158 0   0   nm w=10u l=1u
Mp159  n159 n158 vdd vdd pm w=20u l=1u
C159   n159 0  100f

Mn160  n160 n159 0   0   nm w=10u l=1u
Mp160  n160 n159 vdd vdd pm w=20u l=1u
C160   n160 0  100f

Mn161  n161 n160 0   0   nm w=10u l=1u
Mp161  n161 n160 vdd vdd pm w=20u l=1u
C161   n161 0  100f

Mn162  n162 n161 0   0   nm w=10u l=1u
Mp162  n162 n161 vdd vdd pm w=20u l=1u
C162   n162 0  100f

Mn163  n163 n162 0   0   nm w=10u l=1u
Mp163  n163 n162 vdd vdd pm w=20u l=1u
C163   n163 0  100f

Mn164  n164 n163 0   0   nm w=10u l=1u
Mp164  n164 n163 vdd vdd pm w=20u l=1u
C164   n164 0  100f

Mn165  n165 n164 0   0   nm w=10u l=1u
Mp165  n165 n164 vdd vdd pm w=20u l=1u
C165   n165 0  100f

Mn166  n166 n165 0   0   nm w=10u l=1u
Mp166  n166 n165 vdd vdd pm w=20u l=1u
C166   n166 0  100f

Mn167  n167 n166 0   0   nm w=10u l=1u
Mp167  n167 n166 vdd vdd pm w=20u l=1u
C167   n167 0  100f

Mn168  n168 n167 0   0   nm w=10u l=1u
Mp168  n168 n167 vdd vdd pm w=20u l=1u
C168   n168 0  100f

Mn169  n169 n168 0   0   nm w=10u l=1u
Mp169  n169 n168 vdd vdd pm w=20u l=1u
C169   n169 0  100f

Mn170  n170 n169 0   0   nm w=10u l=1u
Mp170  n170 n169 vdd vdd pm w=20u l=1u
C170   n170 0  100f

Mn171  n171 n170 0   0   nm w=10u l=1u
Mp171  n171 n170 vdd vdd pm w=20u l=1u
C171   n171 0  100f

Mn172  n172 n171 0   0   nm w=10u l=1u
Mp172  n172 n171 vdd vdd pm w=20u l=1u
C172   n172 0  100f

Mn173  n173 n172 0   0   nm w=10u l=1u
Mp173  n173 n172 vdd vdd pm w=20u l=1u
C173   n173 0  100f

Mn174  n174 n173 0   0   nm w=10u l=1u
Mp174  n174 n173 vdd vdd pm w=20u l=1u
C174   n174 0  100f

Mn175  n175 n174 0   0   nm w=10u l=1u
Mp175  n175 n174 vdd vdd pm w=20u l=1u
C175   n175 0  100f

Mn176  n176 n175 0   0   nm w=10u l=1u
Mp176  n176 n175 vdd vdd pm w=20u l=1u
C176   n176 0  100f

Mn177  n177 n176 0   0   nm w=10u l=1u
Mp177  n177 n176 vdd vdd pm w=20u l=1u
C177   n177 0  100f

Mn178  n178 n177 0   0   nm w=10u l=1u
Mp178  n178 n177 vdd vdd pm w=20u l=1u
C178   n178 0  100f

Mn179  n179 n178 0   0   nm w=10u l=1u
Mp179  n179 n178 vdd vdd pm w=20u l=1u
C179   n179 0  100f

Mn180  n180 n179 0   0   nm w=10u l=1u
Mp180  n180 n179 vdd vdd pm w=20u l=1u
C180   n180 0  100f

Mn181  n181 n180 0   0   nm w=10u l=1u
Mp181  n181 n180 vdd vdd pm w=20u l=1u
C181   n181 0  100f

Mn182  n182 n181 0   0   nm w=10u l=1u
Mp182  n182 n181 vdd vdd pm w=20u l=1u
C182   n182 0  100f

Mn183  n183 n182 0   0   nm w=10u l=1u
Mp183  n183 n182 vdd vdd pm w=20u l=1u
C183   n183 0  100f

Mn184  n184 n183 0   0   nm w=10u l=1u
Mp184  n184 n183 vdd vdd pm w=20u l=1u
C184   n184 0  100f

Mn185  n185 n184 0   0   nm w=10u l=1u
Mp185  n185 n184 vdd vdd pm w=20u l=1u
C185   n185 0  100f

Mn186  n186 n185 0   0   nm w=10u l=1u
Mp186  n186 n185 vdd vdd pm w=20u l=1u
C186   n186 0  100f

Mn187  n187 n186 0   0   nm w=10u l=1u
Mp187  n187 n186 vdd vdd pm w=20u l=1u
C187   n187 0  100f

Mn188  n188 n187 0   0   nm w=10u l=1u
Mp188  n188 n187 vdd vdd pm w=20u l=1u
C188   n188 0  100f

Mn189  n189 n188 0   0   nm w=10u l=1u
Mp189  n189 n188 vdd vdd pm w=20u l=1u
C189   n189 0  100f

Mn190  n190 n189 0   0   nm w=10u l=1u
Mp190  n190 n189 vdd vdd pm w=20u l=1u
C190   n190 0  100f

Mn191  n191 n190 0   0   nm w=10u l=1u
Mp191  n191 n190 vdd vdd pm w=20u l=1u
C191   n191 0  100f

Mn192  n192 n191 0   0   nm w=10u l=1u
Mp192  n192 n191 vdd vdd pm w=20u l=1u
C192   n192 0  100f

Mn193  n193 n192 0   0   nm w=10u l=1u
Mp193  n193 n192 vdd vdd pm w=20u l=1u
C193   n193 0  100f

Mn194  n194 n193 0   0   nm w=10u l=1u
Mp194  n194 n193 vdd vdd pm w=20u l=1u
C194   n194 0  100f

Mn195  n195 n194 0   0   nm w=10u l=1u
Mp195  n195 n194 vdd vdd pm w=20u l=1u
C195   n195 0  100f

Mn196  n196 n195 0   0   nm w=10u l=1u
Mp196  n196 n195 vdd vdd pm w=20u l=1u
C196   n196 0  100f

Mn197  n197 n196 0   0   nm w=10u l=1u
Mp197  n197 n196 vdd vdd pm w=20u l=1u
C197   n197 0  100f

Mn198  n198 n197 0   0   nm w=10u l=1u
Mp198  n198 n197 vdd vdd pm w=20u l=1u
C198   n198 0  100f

Mn199  n199 n198 0   0   nm w=10u l=1u
Mp199  n199 n198 vdd vdd pm w=20u l=1u
C199   n199 0  100f

Mn200  n200 n199 0   0   nm w=10u l=1u
Mp200  n200 n199 vdd vdd pm w=20u l=1u
C200   n200 0  100f

Mn201  n201 n200 0   0   nm w=10u l=1u
Mp201  n201 n200 vdd vdd pm w=20u l=1u
C201   n201 0  100f

Mn202  n202 n201 0   0   nm w=10u l=1u
Mp202  n202 n201 vdd vdd pm w=20u l=1u
C202   n202 0  100f

Mn203  n203 n202 0   0   nm w=10u l=1u
Mp203  n203 n202 vdd vdd pm w=20u l=1u
C203   n203 0  100f

Mn204  n204 n203 0   0   nm w=10u l=1u
Mp204  n204 n203 vdd vdd pm w=20u l=1u
C204   n204 0  100f

Mn205  n205 n204 0   0   nm w=10u l=1u
Mp205  n205 n204 vdd vdd pm w=20u l=1u
C205   n205 0  100f

Mn206  n206 n205 0   0   nm w=10u l=1u
Mp206  n206 n205 vdd vdd pm w=20u l=1u
C206   n206 0  100f

Mn207  n207 n206 0   0   nm w=10u l=1u
Mp207  n207 n206 vdd vdd pm w=20u l=1u
C207   n207 0  100f

Mn208  n208 n207 0   0   nm w=10u l=1u
Mp208  n208 n207 vdd vdd pm w=20u l=1u
C208   n208 0  100f

Mn209  n209 n208 0   0   nm w=10u l=1u
Mp209  n209 n208 vdd vdd pm w=20u l=1u
C209   n209 0  100f

Mn210  n210 n209 0   0   nm w=10u l=1u
Mp210  n210 n209 vdd vdd pm w=20u l=1u
C210   n210 0  100f

Mn211  n211 n210 0   0   nm w=10u l=1u
Mp211  n211 n210 vdd vdd pm w=20u l=1u
C211   n211 0  100f

Mn212  n212 n211 0   0   nm w=10u l=1u
Mp212  n212 n211 vdd vdd pm w=20u l=1u
C212   n212 0  100f

Mn213  n213 n212 0   0   nm w=10u l=1u
Mp213  n213 n212 vdd vdd pm w=20u l=1u
C213   n213 0  100f

Mn214  n214 n213 0   0   nm w=10u l=1u
Mp214  n214 n213 vdd vdd pm w=20u l=1u
C214   n214 0  100f

Mn215  n215 n214 0   0   nm w=10u l=1u
Mp215  n215 n214 vdd vdd pm w=20u l=1u
C215   n215 0  100f

Mn216  n216 n215 0   0   nm w=10u l=1u
Mp216  n216 n215 vdd vdd pm w=20u l=1u
C216   n216 0  100f

Mn217  n217 n216 0   0   nm w=10u l=1u
Mp217  n217 n216 vdd vdd pm w=20u l=1u
C217   n217 0  100f

Mn218  n218 n217 0   0   nm w=10u l=1u
Mp218  n218 n217 vdd vdd pm w=20u l=1u
C218   n218 0  100f

Mn219  n219 n218 0   0   nm w=10u l=1u
Mp219  n219 n218 vdd vdd pm w=20u l=1u
C219   n219 0  100f

Mn220  n220 n219 0   0   nm w=10u l=1u
Mp220  n220 n219 vdd vdd pm w=20u l=1u
C220   n220 0  100f

Mn221  n221 n220 0   0   nm w=10u l=1u
Mp221  n221 n220 vdd vdd pm w=20u l=1u
C221   n221 0  100f

Mn222  n222 n221 0   0   nm w=10u l=1u
Mp222  n222 n221 vdd vdd pm w=20u l=1u
C222   n222 0  100f

Mn223  n223 n222 0   0   nm w=10u l=1u
Mp223  n223 n222 vdd vdd pm w=20u l=1u
C223   n223 0  100f

Mn224  n224 n223 0   0   nm w=10u l=1u
Mp224  n224 n223 vdd vdd pm w=20u l=1u
C224   n224 0  100f

Mn225  n225 n224 0   0   nm w=10u l=1u
Mp225  n225 n224 vdd vdd pm w=20u l=1u
C225   n225 0  100f

Mn226  n226 n225 0   0   nm w=10u l=1u
Mp226  n226 n225 vdd vdd pm w=20u l=1u
C226   n226 0  100f

Mn227  n227 n226 0   0   nm w=10u l=1u
Mp227  n227 n226 vdd vdd pm w=20u l=1u
C227   n227 0  100f

Mn228  n228 n227 0   0   nm w=10u l=1u
Mp228  n228 n227 vdd vdd pm w=20u l=1u
C228   n228 0  100f

Mn229  n229 n228 0   0   nm w=10u l=1u
Mp229  n229 n228 vdd vdd pm w=20u l=1u
C229   n229 0  100f

Mn230  n230 n229 0   0   nm w=10u l=1u
Mp230  n230 n229 vdd vdd pm w=20u l=1u
C230   n230 0  100f

Mn231  n231 n230 0   0   nm w=10u l=1u
Mp231  n231 n230 vdd vdd pm w=20u l=1u
C231   n231 0  100f

Mn232  n232 n231 0   0   nm w=10u l=1u
Mp232  n232 n231 vdd vdd pm w=20u l=1u
C232   n232 0  100f

Mn233  n233 n232 0   0   nm w=10u l=1u
Mp233  n233 n232 vdd vdd pm w=20u l=1u
C233   n233 0  100f

Mn234  n234 n233 0   0   nm w=10u l=1u
Mp234  n234 n233 vdd vdd pm w=20u l=1u
C234   n234 0  100f

Mn235  n235 n234 0   0   nm w=10u l=1u
Mp235  n235 n234 vdd vdd pm w=20u l=1u
C235   n235 0  100f

Mn236  n236 n235 0   0   nm w=10u l=1u
Mp236  n236 n235 vdd vdd pm w=20u l=1u
C236   n236 0  100f

Mn237  n237 n236 0   0   nm w=10u l=1u
Mp237  n237 n236 vdd vdd pm w=20u l=1u
C237   n237 0  100f

Mn238  n238 n237 0   0   nm w=10u l=1u
Mp238  n238 n237 vdd vdd pm w=20u l=1u
C238   n238 0  100f

Mn239  n239 n238 0   0   nm w=10u l=1u
Mp239  n239 n238 vdd vdd pm w=20u l=1u
C239   n239 0  100f

Mn240  n240 n239 0   0   nm w=10u l=1u
Mp240  n240 n239 vdd vdd pm w=20u l=1u
C240   n240 0  100f

Mn241  n241 n240 0   0   nm w=10u l=1u
Mp241  n241 n240 vdd vdd pm w=20u l=1u
C241   n241 0  100f

Mn242  n242 n241 0   0   nm w=10u l=1u
Mp242  n242 n241 vdd vdd pm w=20u l=1u
C242   n242 0  100f

Mn243  n243 n242 0   0   nm w=10u l=1u
Mp243  n243 n242 vdd vdd pm w=20u l=1u
C243   n243 0  100f

Mn244  n244 n243 0   0   nm w=10u l=1u
Mp244  n244 n243 vdd vdd pm w=20u l=1u
C244   n244 0  100f

Mn245  n245 n244 0   0   nm w=10u l=1u
Mp245  n245 n244 vdd vdd pm w=20u l=1u
C245   n245 0  100f

Mn246  n246 n245 0   0   nm w=10u l=1u
Mp246  n246 n245 vdd vdd pm w=20u l=1u
C246   n246 0  100f

Mn247  n247 n246 0   0   nm w=10u l=1u
Mp247  n247 n246 vdd vdd pm w=20u l=1u
C247   n247 0  100f

Mn248  n248 n247 0   0   nm w=10u l=1u
Mp248  n248 n247 vdd vdd pm w=20u l=1u
C248   n248 0  100f

Mn249  n249 n248 0   0   nm w=10u l=1u
Mp249  n249 n248 vdd vdd pm w=20u l=1u
C249   n249 0  100f

Mn250  n250 n249 0   0   nm w=10u l=1u
Mp250  n250 n249 vdd vdd pm w=20u l=1u
C250   n250 0  100f

Mn251  n251 n250 0   0   nm w=10u l=1u
Mp251  n251 n250 vdd vdd pm w=20u l=1u
C251   n251 0  100f

Mn252  n252 n251 0   0   nm w=10u l=1u
Mp252  n252 n251 vdd vdd pm w=20u l=1u
C252   n252 0  100f

Mn253  n253 n252 0   0   nm w=10u l=1u
Mp253  n253 n252 vdd vdd pm w=20u l=1u
C253   n253 0  100f

Mn254  n254 n253 0   0   nm w=10u l=1u
Mp254  n254 n253 vdd vdd pm w=20u l=1u
C254   n254 0  100f

Mn255  n255 n254 0   0   nm w=10u l=1u
Mp255  n255 n254 vdd vdd pm w=20u l=1u
C255   n255 0  100f

Mn256  n256 n255 0   0   nm w=10u l=1u
Mp256  n256 n255 vdd vdd pm w=20u l=1u
C256   n256 0  100f

Mn257  n257 n256 0   0   nm w=10u l=1u
Mp257  n257 n256 vdd vdd pm w=20u l=1u
C257   n257 0  100f

Mn258  n258 n257 0   0   nm w=10u l=1u
Mp258  n258 n257 vdd vdd pm w=20u l=1u
C258   n258 0  100f

Mn259  n259 n258 0   0   nm w=10u l=1u
Mp259  n259 n258 vdd vdd pm w=20u l=1u
C259   n259 0  100f

Mn260  n260 n259 0   0   nm w=10u l=1u
Mp260  n260 n259 vdd vdd pm w=20u l=1u
C260   n260 0  100f

Mn261  n261 n260 0   0   nm w=10u l=1u
Mp261  n261 n260 vdd vdd pm w=20u l=1u
C261   n261 0  100f

Mn262  n262 n261 0   0   nm w=10u l=1u
Mp262  n262 n261 vdd vdd pm w=20u l=1u
C262   n262 0  100f

Mn263  n263 n262 0   0   nm w=10u l=1u
Mp263  n263 n262 vdd vdd pm w=20u l=1u
C263   n263 0  100f

Mn264  n264 n263 0   0   nm w=10u l=1u
Mp264  n264 n263 vdd vdd pm w=20u l=1u
C264   n264 0  100f

Mn265  n265 n264 0   0   nm w=10u l=1u
Mp265  n265 n264 vdd vdd pm w=20u l=1u
C265   n265 0  100f

Mn266  n266 n265 0   0   nm w=10u l=1u
Mp266  n266 n265 vdd vdd pm w=20u l=1u
C266   n266 0  100f

Mn267  n267 n266 0   0   nm w=10u l=1u
Mp267  n267 n266 vdd vdd pm w=20u l=1u
C267   n267 0  100f

Mn268  n268 n267 0   0   nm w=10u l=1u
Mp268  n268 n267 vdd vdd pm w=20u l=1u
C268   n268 0  100f

Mn269  n269 n268 0   0   nm w=10u l=1u
Mp269  n269 n268 vdd vdd pm w=20u l=1u
C269   n269 0  100f

Mn270  n270 n269 0   0   nm w=10u l=1u
Mp270  n270 n269 vdd vdd pm w=20u l=1u
C270   n270 0  100f

Mn271  n271 n270 0   0   nm w=10u l=1u
Mp271  n271 n270 vdd vdd pm w=20u l=1u
C271   n271 0  100f

Mn272  n272 n271 0   0   nm w=10u l=1u
Mp272  n272 n271 vdd vdd pm w=20u l=1u
C272   n272 0  100f

Mn273  n273 n272 0   0   nm w=10u l=1u
Mp273  n273 n272 vdd vdd pm w=20u l=1u
C273   n273 0  100f

Mn274  n274 n273 0   0   nm w=10u l=1u
Mp274  n274 n273 vdd vdd pm w=20u l=1u
C274   n274 0  100f

Mn275  n275 n274 0   0   nm w=10u l=1u
Mp275  n275 n274 vdd vdd pm w=20u l=1u
C275   n275 0  100f

Mn276  n276 n275 0   0   nm w=10u l=1u
Mp276  n276 n275 vdd vdd pm w=20u l=1u
C276   n276 0  100f

Mn277  n277 n276 0   0   nm w=10u l=1u
Mp277  n277 n276 vdd vdd pm w=20u l=1u
C277   n277 0  100f

Mn278  n278 n277 0   0   nm w=10u l=1u
Mp278  n278 n277 vdd vdd pm w=20u l=1u
C278   n278 0  100f

Mn279  n279 n278 0   0   nm w=10u l=1u
Mp279  n279 n278 vdd vdd pm w=20u l=1u
C279   n279 0  100f

Mn280  n280 n279 0   0   nm w=10u l=1u
Mp280  n280 n279 vdd vdd pm w=20u l=1u
C280   n280 0  100f

Mn281  n281 n280 0   0   nm w=10u l=1u
Mp281  n281 n280 vdd vdd pm w=20u l=1u
C281   n281 0  100f

Mn282  n282 n281 0   0   nm w=10u l=1u
Mp282  n282 n281 vdd vdd pm w=20u l=1u
C282   n282 0  100f

Mn283  n283 n282 0   0   nm w=10u l=1u
Mp283  n283 n282 vdd vdd pm w=20u l=1u
C283   n283 0  100f

Mn284  n284 n283 0   0   nm w=10u l=1u
Mp284  n284 n283 vdd vdd pm w=20u l=1u
C284   n284 0  100f

Mn285  n285 n284 0   0   nm w=10u l=1u
Mp285  n285 n284 vdd vdd pm w=20u l=1u
C285   n285 0  100f

Mn286  n286 n285 0   0   nm w=10u l=1u
Mp286  n286 n285 vdd vdd pm w=20u l=1u
C286   n286 0  100f

Mn287  n287 n286 0   0   nm w=10u l=1u
Mp287  n287 n286 vdd vdd pm w=20u l=1u
C287   n287 0  100f

Mn288  n288 n287 0   0   nm w=10u l=1u
Mp288  n288 n287 vdd vdd pm w=20u l=1u
C288   n288 0  100f

Mn289  n289 n288 0   0   nm w=10u l=1u
Mp289  n289 n288 vdd vdd pm w=20u l=1u
C289   n289 0  100f

Mn290  n290 n289 0   0   nm w=10u l=1u
Mp290  n290 n289 vdd vdd pm w=20u l=1u
C290   n290 0  100f

Mn291  n291 n290 0   0   nm w=10u l=1u
Mp291  n291 n290 vdd vdd pm w=20u l=1u
C291   n291 0  100f

Mn292  n292 n291 0   0   nm w=10u l=1u
Mp292  n292 n291 vdd vdd pm w=20u l=1u
C292   n292 0  100f

Mn293  n293 n292 0   0   nm w=10u l=1u
Mp293  n293 n292 vdd vdd pm w=20u l=1u
C293   n293 0  100f

Mn294  n294 n293 0   0   nm w=10u l=1u
Mp294  n294 n293 vdd vdd pm w=20u l=1u
C294   n294 0  100f

Mn295  n295 n294 0   0   nm w=10u l=1u
Mp295  n295 n294 vdd vdd pm w=20u l=1u
C295   n295 0  100f

Mn296  n296 n295 0   0   nm w=10u l=1u
Mp296  n296 n295 vdd vdd pm w=20u l=1u
C296   n296 0  100f

Mn297  n297 n296 0   0   nm w=10u l=1u
Mp297  n297 n296 vdd vdd pm w=20u l=1u
C297   n297 0  100f

Mn298  n298 n297 0   0   nm w=10u l=1u
Mp298  n298 n297 vdd vdd pm w=20u l=1u
C298   n298 0  100f

Mn299  n299 n298 0   0   nm w=10u l=1u
Mp299  n299 n298 vdd vdd pm w=20u l=1u
C299   n299 0  100f

Mn300  n300 n299 0   0   nm w=10u l=1u
Mp300  n300 n299 vdd vdd pm w=20u l=1u
C300   n300 0  100f

Mn301  n301 n300 0   0   nm w=10u l=1u
Mp301  n301 n300 vdd vdd pm w=20u l=1u
C301   n301 0  100f

Mn302  n302 n301 0   0   nm w=10u l=1u
Mp302  n302 n301 vdd vdd pm w=20u l=1u
C302   n302 0  100f

Mn303  n303 n302 0   0   nm w=10u l=1u
Mp303  n303 n302 vdd vdd pm w=20u l=1u
C303   n303 0  100f

Mn304  n304 n303 0   0   nm w=10u l=1u
Mp304  n304 n303 vdd vdd pm w=20u l=1u
C304   n304 0  100f

Mn305  n305 n304 0   0   nm w=10u l=1u
Mp305  n305 n304 vdd vdd pm w=20u l=1u
C305   n305 0  100f

Mn306  n306 n305 0   0   nm w=10u l=1u
Mp306  n306 n305 vdd vdd pm w=20u l=1u
C306   n306 0  100f

Mn307  n307 n306 0   0   nm w=10u l=1u
Mp307  n307 n306 vdd vdd pm w=20u l=1u
C307   n307 0  100f

Mn308  n308 n307 0   0   nm w=10u l=1u
Mp308  n308 n307 vdd vdd pm w=20u l=1u
C308   n308 0  100f

Mn309  n309 n308 0   0   nm w=10u l=1u
Mp309  n309 n308 vdd vdd pm w=20u l=1u
C309   n309 0  100f

Mn310  n310 n309 0   0   nm w=10u l=1u
Mp310  n310 n309 vdd vdd pm w=20u l=1u
C310   n310 0  100f

Mn311  n311 n310 0   0   nm w=10u l=1u
Mp311  n311 n310 vdd vdd pm w=20u l=1u
C311   n311 0  100f

Mn312  n312 n311 0   0   nm w=10u l=1u
Mp312  n312 n311 vdd vdd pm w=20u l=1u
C312   n312 0  100f

Mn313  n313 n312 0   0   nm w=10u l=1u
Mp313  n313 n312 vdd vdd pm w=20u l=1u
C313   n313 0  100f

Mn314  n314 n313 0   0   nm w=10u l=1u
Mp314  n314 n313 vdd vdd pm w=20u l=1u
C314   n314 0  100f

Mn315  n315 n314 0   0   nm w=10u l=1u
Mp315  n315 n314 vdd vdd pm w=20u l=1u
C315   n315 0  100f

Mn316  n316 n315 0   0   nm w=10u l=1u
Mp316  n316 n315 vdd vdd pm w=20u l=1u
C316   n316 0  100f

Mn317  n317 n316 0   0   nm w=10u l=1u
Mp317  n317 n316 vdd vdd pm w=20u l=1u
C317   n317 0  100f

Mn318  n318 n317 0   0   nm w=10u l=1u
Mp318  n318 n317 vdd vdd pm w=20u l=1u
C318   n318 0  100f

Mn319  n319 n318 0   0   nm w=10u l=1u
Mp319  n319 n318 vdd vdd pm w=20u l=1u
C319   n319 0  100f

Mn320  n320 n319 0   0   nm w=10u l=1u
Mp320  n320 n319 vdd vdd pm w=20u l=1u
C320   n320 0  100f

Mn321  n321 n320 0   0   nm w=10u l=1u
Mp321  n321 n320 vdd vdd pm w=20u l=1u
C321   n321 0  100f

Mn322  n322 n321 0   0   nm w=10u l=1u
Mp322  n322 n321 vdd vdd pm w=20u l=1u
C322   n322 0  100f

Mn323  n323 n322 0   0   nm w=10u l=1u
Mp323  n323 n322 vdd vdd pm w=20u l=1u
C323   n323 0  100f

Mn324  n324 n323 0   0   nm w=10u l=1u
Mp324  n324 n323 vdd vdd pm w=20u l=1u
C324   n324 0  100f

Mn325  n325 n324 0   0   nm w=10u l=1u
Mp325  n325 n324 vdd vdd pm w=20u l=1u
C325   n325 0  100f

Mn326  n326 n325 0   0   nm w=10u l=1u
Mp326  n326 n325 vdd vdd pm w=20u l=1u
C326   n326 0  100f

Mn327  n327 n326 0   0   nm w=10u l=1u
Mp327  n327 n326 vdd vdd pm w=20u l=1u
C327   n327 0  100f

Mn328  n328 n327 0   0   nm w=10u l=1u
Mp328  n328 n327 vdd vdd pm w=20u l=1u
C328   n328 0  100f

Mn329  n329 n328 0   0   nm w=10u l=1u
Mp329  n329 n328 vdd vdd pm w=20u l=1u
C329   n329 0  100f

Mn330  n330 n329 0   0   nm w=10u l=1u
Mp330  n330 n329 vdd vdd pm w=20u l=1u
C330   n330 0  100f

Mn331  n331 n330 0   0   nm w=10u l=1u
Mp331  n331 n330 vdd vdd pm w=20u l=1u
C331   n331 0  100f

Mn332  n332 n331 0   0   nm w=10u l=1u
Mp332  n332 n331 vdd vdd pm w=20u l=1u
C332   n332 0  100f

Mn333  n333 n332 0   0   nm w=10u l=1u
Mp333  n333 n332 vdd vdd pm w=20u l=1u
C333   n333 0  100f

Mn334  n334 n333 0   0   nm w=10u l=1u
Mp334  n334 n333 vdd vdd pm w=20u l=1u
C334   n334 0  100f

Mn335  n335 n334 0   0   nm w=10u l=1u
Mp335  n335 n334 vdd vdd pm w=20u l=1u
C335   n335 0  100f

Mn336  n336 n335 0   0   nm w=10u l=1u
Mp336  n336 n335 vdd vdd pm w=20u l=1u
C336   n336 0  100f

Mn337  n337 n336 0   0   nm w=10u l=1u
Mp337  n337 n336 vdd vdd pm w=20u l=1u
C337   n337 0  100f

Mn338  n338 n337 0   0   nm w=10u l=1u
Mp338  n338 n337 vdd vdd pm w=20u l=1u
C338   n338 0  100f

Mn339  n339 n338 0   0   nm w=10u l=1u
Mp339  n339 n338 vdd vdd pm w=20u l=1u
C339   n339 0  100f

Mn340  n340 n339 0   0   nm w=10u l=1u
Mp340  n340 n339 vdd vdd pm w=20u l=1u
C340   n340 0  100f

Mn341  n341 n340 0   0   nm w=10u l=1u
Mp341  n341 n340 vdd vdd pm w=20u l=1u
C341   n341 0  100f

Mn342  n342 n341 0   0   nm w=10u l=1u
Mp342  n342 n341 vdd vdd pm w=20u l=1u
C342   n342 0  100f

Mn343  n343 n342 0   0   nm w=10u l=1u
Mp343  n343 n342 vdd vdd pm w=20u l=1u
C343   n343 0  100f

Mn344  n344 n343 0   0   nm w=10u l=1u
Mp344  n344 n343 vdd vdd pm w=20u l=1u
C344   n344 0  100f

Mn345  n345 n344 0   0   nm w=10u l=1u
Mp345  n345 n344 vdd vdd pm w=20u l=1u
C345   n345 0  100f

Mn346  n346 n345 0   0   nm w=10u l=1u
Mp346  n346 n345 vdd vdd pm w=20u l=1u
C346   n346 0  100f

Mn347  n347 n346 0   0   nm w=10u l=1u
Mp347  n347 n346 vdd vdd pm w=20u l=1u
C347   n347 0  100f

Mn348  n348 n347 0   0   nm w=10u l=1u
Mp348  n348 n347 vdd vdd pm w=20u l=1u
C348   n348 0  100f

Mn349  n349 n348 0   0   nm w=10u l=1u
Mp349  n349 n348 vdd vdd pm w=20u l=1u
C349   n349 0  100f

Mn350  n350 n349 0   0   nm w=10u l=1u
Mp350  n350 n349 vdd vdd pm w=20u l=1u
C350   n350 0  100f

Mn351  n351 n350 0   0   nm w=10u l=1u
Mp351  n351 n350 vdd vdd pm w=20u l=1u
C351   n351 0  100f

Mn352  n352 n351 0   0   nm w=10u l=1u
Mp352  n352 n351 vdd vdd pm w=20u l=1u
C352   n352 0  100f

Mn353  n353 n352 0   0   nm w=10u l=1u
Mp353  n353 n352 vdd vdd pm w=20u l=1u
C353   n353 0  100f

Mn354  n354 n353 0   0   nm w=10u l=1u
Mp354  n354 n353 vdd vdd pm w=20u l=1u
C354   n354 0  100f

Mn355  n355 n354 0   0   nm w=10u l=1u
Mp355  n355 n354 vdd vdd pm w=20u l=1u
C355   n355 0  100f

Mn356  n356 n355 0   0   nm w=10u l=1u
Mp356  n356 n355 vdd vdd pm w=20u l=1u
C356   n356 0  100f

Mn357  n357 n356 0   0   nm w=10u l=1u
Mp357  n357 n356 vdd vdd pm w=20u l=1u
C357   n357 0  100f

Mn358  n358 n357 0   0   nm w=10u l=1u
Mp358  n358 n357 vdd vdd pm w=20u l=1u
C358   n358 0  100f

Mn359  n359 n358 0   0   nm w=10u l=1u
Mp359  n359 n358 vdd vdd pm w=20u l=1u
C359   n359 0  100f

Mn360  n360 n359 0   0   nm w=10u l=1u
Mp360  n360 n359 vdd vdd pm w=20u l=1u
C360   n360 0  100f

Mn361  n361 n360 0   0   nm w=10u l=1u
Mp361  n361 n360 vdd vdd pm w=20u l=1u
C361   n361 0  100f

Mn362  n362 n361 0   0   nm w=10u l=1u
Mp362  n362 n361 vdd vdd pm w=20u l=1u
C362   n362 0  100f

Mn363  n363 n362 0   0   nm w=10u l=1u
Mp363  n363 n362 vdd vdd pm w=20u l=1u
C363   n363 0  100f

Mn364  n364 n363 0   0   nm w=10u l=1u
Mp364  n364 n363 vdd vdd pm w=20u l=1u
C364   n364 0  100f

Mn365  n365 n364 0   0   nm w=10u l=1u
Mp365  n365 n364 vdd vdd pm w=20u l=1u
C365   n365 0  100f

Mn366  n366 n365 0   0   nm w=10u l=1u
Mp366  n366 n365 vdd vdd pm w=20u l=1u
C366   n366 0  100f

Mn367  n367 n366 0   0   nm w=10u l=1u
Mp367  n367 n366 vdd vdd pm w=20u l=1u
C367   n367 0  100f

Mn368  n368 n367 0   0   nm w=10u l=1u
Mp368  n368 n367 vdd vdd pm w=20u l=1u
C368   n368 0  100f

Mn369  n369 n368 0   0   nm w=10u l=1u
Mp369  n369 n368 vdd vdd pm w=20u l=1u
C369   n369 0  100f

Mn370  n370 n369 0   0   nm w=10u l=1u
Mp370  n370 n369 vdd vdd pm w=20u l=1u
C370   n370 0  100f

Mn371  n371 n370 0   0   nm w=10u l=1u
Mp371  n371 n370 vdd vdd pm w=20u l=1u
C371   n371 0  100f

Mn372  n372 n371 0   0   nm w=10u l=1u
Mp372  n372 n371 vdd vdd pm w=20u l=1u
C372   n372 0  100f

Mn373  n373 n372 0   0   nm w=10u l=1u
Mp373  n373 n372 vdd vdd pm w=20u l=1u
C373   n373 0  100f

Mn374  n374 n373 0   0   nm w=10u l=1u
Mp374  n374 n373 vdd vdd pm w=20u l=1u
C374   n374 0  100f

Mn375  n375 n374 0   0   nm w=10u l=1u
Mp375  n375 n374 vdd vdd pm w=20u l=1u
C375   n375 0  100f

Mn376  n376 n375 0   0   nm w=10u l=1u
Mp376  n376 n375 vdd vdd pm w=20u l=1u
C376   n376 0  100f

Mn377  n377 n376 0   0   nm w=10u l=1u
Mp377  n377 n376 vdd vdd pm w=20u l=1u
C377   n377 0  100f

Mn378  n378 n377 0   0   nm w=10u l=1u
Mp378  n378 n377 vdd vdd pm w=20u l=1u
C378   n378 0  100f

Mn379  n379 n378 0   0   nm w=10u l=1u
Mp379  n379 n378 vdd vdd pm w=20u l=1u
C379   n379 0  100f

Mn380  n380 n379 0   0   nm w=10u l=1u
Mp380  n380 n379 vdd vdd pm w=20u l=1u
C380   n380 0  100f

Mn381  n381 n380 0   0   nm w=10u l=1u
Mp381  n381 n380 vdd vdd pm w=20u l=1u
C381   n381 0  100f

Mn382  n382 n381 0   0   nm w=10u l=1u
Mp382  n382 n381 vdd vdd pm w=20u l=1u
C382   n382 0  100f

Mn383  n383 n382 0   0   nm w=10u l=1u
Mp383  n383 n382 vdd vdd pm w=20u l=1u
C383   n383 0  100f

Mn384  n384 n383 0   0   nm w=10u l=1u
Mp384  n384 n383 vdd vdd pm w=20u l=1u
C384   n384 0  100f

Mn385  n385 n384 0   0   nm w=10u l=1u
Mp385  n385 n384 vdd vdd pm w=20u l=1u
C385   n385 0  100f

Mn386  n386 n385 0   0   nm w=10u l=1u
Mp386  n386 n385 vdd vdd pm w=20u l=1u
C386   n386 0  100f

Mn387  n387 n386 0   0   nm w=10u l=1u
Mp387  n387 n386 vdd vdd pm w=20u l=1u
C387   n387 0  100f

Mn388  n388 n387 0   0   nm w=10u l=1u
Mp388  n388 n387 vdd vdd pm w=20u l=1u
C388   n388 0  100f

Mn389  n389 n388 0   0   nm w=10u l=1u
Mp389  n389 n388 vdd vdd pm w=20u l=1u
C389   n389 0  100f

Mn390  n390 n389 0   0   nm w=10u l=1u
Mp390  n390 n389 vdd vdd pm w=20u l=1u
C390   n390 0  100f

Mn391  n391 n390 0   0   nm w=10u l=1u
Mp391  n391 n390 vdd vdd pm w=20u l=1u
C391   n391 0  100f

Mn392  n392 n391 0   0   nm w=10u l=1u
Mp392  n392 n391 vdd vdd pm w=20u l=1u
C392   n392 0  100f

Mn393  n393 n392 0   0   nm w=10u l=1u
Mp393  n393 n392 vdd vdd pm w=20u l=1u
C393   n393 0  100f

Mn394  n394 n393 0   0   nm w=10u l=1u
Mp394  n394 n393 vdd vdd pm w=20u l=1u
C394   n394 0  100f

Mn395  n395 n394 0   0   nm w=10u l=1u
Mp395  n395 n394 vdd vdd pm w=20u l=1u
C395   n395 0  100f

Mn396  n396 n395 0   0   nm w=10u l=1u
Mp396  n396 n395 vdd vdd pm w=20u l=1u
C396   n396 0  100f

Mn397  n397 n396 0   0   nm w=10u l=1u
Mp397  n397 n396 vdd vdd pm w=20u l=1u
C397   n397 0  100f

Mn398  n398 n397 0   0   nm w=10u l=1u
Mp398  n398 n397 vdd vdd pm w=20u l=1u
C398   n398 0  100f

Mn399  n399 n398 0   0   nm w=10u l=1u
Mp399  n399 n398 vdd vdd pm w=20u l=1u
C399   n399 0  100f

Mn400  n400 n399 0   0   nm w=10u l=1u
Mp400  n400 n399 vdd vdd pm w=20u l=1u
C400   n400 0  100f

Mn401  n401 n400 0   0   nm w=10u l=1u
Mp401  n401 n400 vdd vdd pm w=20u l=1u
C401   n401 0  100f

Mn402  n402 n401 0   0   nm w=10u l=1u
Mp402  n402 n401 vdd vdd pm w=20u l=1u
C402   n402 0  100f

Mn403  n403 n402 0   0   nm w=10u l=1u
Mp403  n403 n402 vdd vdd pm w=20u l=1u
C403   n403 0  100f

Mn404  n404 n403 0   0   nm w=10u l=1u
Mp404  n404 n403 vdd vdd pm w=20u l=1u
C404   n404 0  100f

Mn405  n405 n404 0   0   nm w=10u l=1u
Mp405  n405 n404 vdd vdd pm w=20u l=1u
C405   n405 0  100f

Mn406  n406 n405 0   0   nm w=10u l=1u
Mp406  n406 n405 vdd vdd pm w=20u l=1u
C406   n406 0  100f

Mn407  n407 n406 0   0   nm w=10u l=1u
Mp407  n407 n406 vdd vdd pm w=20u l=1u
C407   n407 0  100f

Mn408  n408 n407 0   0   nm w=10u l=1u
Mp408  n408 n407 vdd vdd pm w=20u l=1u
C408   n408 0  100f

Mn409  n409 n408 0   0   nm w=10u l=1u
Mp409  n409 n408 vdd vdd pm w=20u l=1u
C409   n409 0  100f

Mn410  n410 n409 0   0   nm w=10u l=1u
Mp410  n410 n409 vdd vdd pm w=20u l=1u
C410   n410 0  100f

Mn411  n411 n410 0   0   nm w=10u l=1u
Mp411  n411 n410 vdd vdd pm w=20u l=1u
C411   n411 0  100f

Mn412  n412 n411 0   0   nm w=10u l=1u
Mp412  n412 n411 vdd vdd pm w=20u l=1u
C412   n412 0  100f

Mn413  n413 n412 0   0   nm w=10u l=1u
Mp413  n413 n412 vdd vdd pm w=20u l=1u
C413   n413 0  100f

Mn414  n414 n413 0   0   nm w=10u l=1u
Mp414  n414 n413 vdd vdd pm w=20u l=1u
C414   n414 0  100f

Mn415  n415 n414 0   0   nm w=10u l=1u
Mp415  n415 n414 vdd vdd pm w=20u l=1u
C415   n415 0  100f

Mn416  n416 n415 0   0   nm w=10u l=1u
Mp416  n416 n415 vdd vdd pm w=20u l=1u
C416   n416 0  100f

Mn417  n417 n416 0   0   nm w=10u l=1u
Mp417  n417 n416 vdd vdd pm w=20u l=1u
C417   n417 0  100f

Mn418  n418 n417 0   0   nm w=10u l=1u
Mp418  n418 n417 vdd vdd pm w=20u l=1u
C418   n418 0  100f

Mn419  n419 n418 0   0   nm w=10u l=1u
Mp419  n419 n418 vdd vdd pm w=20u l=1u
C419   n419 0  100f

Mn420  n420 n419 0   0   nm w=10u l=1u
Mp420  n420 n419 vdd vdd pm w=20u l=1u
C420   n420 0  100f

Mn421  n421 n420 0   0   nm w=10u l=1u
Mp421  n421 n420 vdd vdd pm w=20u l=1u
C421   n421 0  100f

Mn422  n422 n421 0   0   nm w=10u l=1u
Mp422  n422 n421 vdd vdd pm w=20u l=1u
C422   n422 0  100f

Mn423  n423 n422 0   0   nm w=10u l=1u
Mp423  n423 n422 vdd vdd pm w=20u l=1u
C423   n423 0  100f

Mn424  n424 n423 0   0   nm w=10u l=1u
Mp424  n424 n423 vdd vdd pm w=20u l=1u
C424   n424 0  100f

Mn425  n425 n424 0   0   nm w=10u l=1u
Mp425  n425 n424 vdd vdd pm w=20u l=1u
C425   n425 0  100f

Mn426  n426 n425 0   0   nm w=10u l=1u
Mp426  n426 n425 vdd vdd pm w=20u l=1u
C426   n426 0  100f

Mn427  n427 n426 0   0   nm w=10u l=1u
Mp427  n427 n426 vdd vdd pm w=20u l=1u
C427   n427 0  100f

Mn428  n428 n427 0   0   nm w=10u l=1u
Mp428  n428 n427 vdd vdd pm w=20u l=1u
C428   n428 0  100f

Mn429  n429 n428 0   0   nm w=10u l=1u
Mp429  n429 n428 vdd vdd pm w=20u l=1u
C429   n429 0  100f

Mn430  n430 n429 0   0   nm w=10u l=1u
Mp430  n430 n429 vdd vdd pm w=20u l=1u
C430   n430 0  100f

Mn431  n431 n430 0   0   nm w=10u l=1u
Mp431  n431 n430 vdd vdd pm w=20u l=1u
C431   n431 0  100f

Mn432  n432 n431 0   0   nm w=10u l=1u
Mp432  n432 n431 vdd vdd pm w=20u l=1u
C432   n432 0  100f

Mn433  n433 n432 0   0   nm w=10u l=1u
Mp433  n433 n432 vdd vdd pm w=20u l=1u
C433   n433 0  100f

Mn434  n434 n433 0   0   nm w=10u l=1u
Mp434  n434 n433 vdd vdd pm w=20u l=1u
C434   n434 0  100f

Mn435  n435 n434 0   0   nm w=10u l=1u
Mp435  n435 n434 vdd vdd pm w=20u l=1u
C435   n435 0  100f

Mn436  n436 n435 0   0   nm w=10u l=1u
Mp436  n436 n435 vdd vdd pm w=20u l=1u
C436   n436 0  100f

Mn437  n437 n436 0   0   nm w=10u l=1u
Mp437  n437 n436 vdd vdd pm w=20u l=1u
C437   n437 0  100f

Mn438  n438 n437 0   0   nm w=10u l=1u
Mp438  n438 n437 vdd vdd pm w=20u l=1u
C438   n438 0  100f

Mn439  n439 n438 0   0   nm w=10u l=1u
Mp439  n439 n438 vdd vdd pm w=20u l=1u
C439   n439 0  100f

Mn440  n440 n439 0   0   nm w=10u l=1u
Mp440  n440 n439 vdd vdd pm w=20u l=1u
C440   n440 0  100f

Mn441  n441 n440 0   0   nm w=10u l=1u
Mp441  n441 n440 vdd vdd pm w=20u l=1u
C441   n441 0  100f

Mn442  n442 n441 0   0   nm w=10u l=1u
Mp442  n442 n441 vdd vdd pm w=20u l=1u
C442   n442 0  100f

Mn443  n443 n442 0   0   nm w=10u l=1u
Mp443  n443 n442 vdd vdd pm w=20u l=1u
C443   n443 0  100f

Mn444  n444 n443 0   0   nm w=10u l=1u
Mp444  n444 n443 vdd vdd pm w=20u l=1u
C444   n444 0  100f

Mn445  n445 n444 0   0   nm w=10u l=1u
Mp445  n445 n444 vdd vdd pm w=20u l=1u
C445   n445 0  100f

Mn446  n446 n445 0   0   nm w=10u l=1u
Mp446  n446 n445 vdd vdd pm w=20u l=1u
C446   n446 0  100f

Mn447  n447 n446 0   0   nm w=10u l=1u
Mp447  n447 n446 vdd vdd pm w=20u l=1u
C447   n447 0  100f

Mn448  n448 n447 0   0   nm w=10u l=1u
Mp448  n448 n447 vdd vdd pm w=20u l=1u
C448   n448 0  100f

Mn449  n449 n448 0   0   nm w=10u l=1u
Mp449  n449 n448 vdd vdd pm w=20u l=1u
C449   n449 0  100f

Mn450  n450 n449 0   0   nm w=10u l=1u
Mp450  n450 n449 vdd vdd pm w=20u l=1u
C450   n450 0  100f

Mn451  n451 n450 0   0   nm w=10u l=1u
Mp451  n451 n450 vdd vdd pm w=20u l=1u
C451   n451 0  100f

Mn452  n452 n451 0   0   nm w=10u l=1u
Mp452  n452 n451 vdd vdd pm w=20u l=1u
C452   n452 0  100f

Mn453  n453 n452 0   0   nm w=10u l=1u
Mp453  n453 n452 vdd vdd pm w=20u l=1u
C453   n453 0  100f

Mn454  n454 n453 0   0   nm w=10u l=1u
Mp454  n454 n453 vdd vdd pm w=20u l=1u
C454   n454 0  100f

Mn455  n455 n454 0   0   nm w=10u l=1u
Mp455  n455 n454 vdd vdd pm w=20u l=1u
C455   n455 0  100f

Mn456  n456 n455 0   0   nm w=10u l=1u
Mp456  n456 n455 vdd vdd pm w=20u l=1u
C456   n456 0  100f

Mn457  n457 n456 0   0   nm w=10u l=1u
Mp457  n457 n456 vdd vdd pm w=20u l=1u
C457   n457 0  100f

Mn458  n458 n457 0   0   nm w=10u l=1u
Mp458  n458 n457 vdd vdd pm w=20u l=1u
C458   n458 0  100f

Mn459  n459 n458 0   0   nm w=10u l=1u
Mp459  n459 n458 vdd vdd pm w=20u l=1u
C459   n459 0  100f

Mn460  n460 n459 0   0   nm w=10u l=1u
Mp460  n460 n459 vdd vdd pm w=20u l=1u
C460   n460 0  100f

Mn461  n461 n460 0   0   nm w=10u l=1u
Mp461  n461 n460 vdd vdd pm w=20u l=1u
C461   n461 0  100f

Mn462  n462 n461 0   0   nm w=10u l=1u
Mp462  n462 n461 vdd vdd pm w=20u l=1u
C462   n462 0  100f

Mn463  n463 n462 0   0   nm w=10u l=1u
Mp463  n463 n462 vdd vdd pm w=20u l=1u
C463   n463 0  100f

Mn464  n464 n463 0   0   nm w=10u l=1u
Mp464  n464 n463 vdd vdd pm w=20u l=1u
C464   n464 0  100f

Mn465  n465 n464 0   0   nm w=10u l=1u
Mp465  n465 n464 vdd vdd pm w=20u l=1u
C465   n465 0  100f

Mn466  n466 n465 0   0   nm w=10u l=1u
Mp466  n466 n465 vdd vdd pm w=20u l=1u
C466   n466 0  100f

Mn467  n467 n466 0   0   nm w=10u l=1u
Mp467  n467 n466 vdd vdd pm w=20u l=1u
C467   n467 0  100f

Mn468  n468 n467 0   0   nm w=10u l=1u
Mp468  n468 n467 vdd vdd pm w=20u l=1u
C468   n468 0  100f

Mn469  n469 n468 0   0   nm w=10u l=1u
Mp469  n469 n468 vdd vdd pm w=20u l=1u
C469   n469 0  100f

Mn470  n470 n469 0   0   nm w=10u l=1u
Mp470  n470 n469 vdd vdd pm w=20u l=1u
C470   n470 0  100f

Mn471  n471 n470 0   0   nm w=10u l=1u
Mp471  n471 n470 vdd vdd pm w=20u l=1u
C471   n471 0  100f

Mn472  n472 n471 0   0   nm w=10u l=1u
Mp472  n472 n471 vdd vdd pm w=20u l=1u
C472   n472 0  100f

Mn473  n473 n472 0   0   nm w=10u l=1u
Mp473  n473 n472 vdd vdd pm w=20u l=1u
C473   n473 0  100f

Mn474  n474 n473 0   0   nm w=10u l=1u
Mp474  n474 n473 vdd vdd pm w=20u l=1u
C474   n474 0  100f

Mn475  n475 n474 0   0   nm w=10u l=1u
Mp475  n475 n474 vdd vdd pm w=20u l=1u
C475   n475 0  100f

Mn476  n476 n475 0   0   nm w=10u l=1u
Mp476  n476 n475 vdd vdd pm w=20u l=1u
C476   n476 0  100f

Mn477  n477 n476 0   0   nm w=10u l=1u
Mp477  n477 n476 vdd vdd pm w=20u l=1u
C477   n477 0  100f

Mn478  n478 n477 0   0   nm w=10u l=1u
Mp478  n478 n477 vdd vdd pm w=20u l=1u
C478   n478 0  100f

Mn479  n479 n478 0   0   nm w=10u l=1u
Mp479  n479 n478 vdd vdd pm w=20u l=1u
C479   n479 0  100f

Mn480  n480 n479 0   0   nm w=10u l=1u
Mp480  n480 n479 vdd vdd pm w=20u l=1u
C480   n480 0  100f

Mn481  n481 n480 0   0   nm w=10u l=1u
Mp481  n481 n480 vdd vdd pm w=20u l=1u
C481   n481 0  100f

Mn482  n482 n481 0   0   nm w=10u l=1u
Mp482  n482 n481 vdd vdd pm w=20u l=1u
C482   n482 0  100f

Mn483  n483 n482 0   0   nm w=10u l=1u
Mp483  n483 n482 vdd vdd pm w=20u l=1u
C483   n483 0  100f

Mn484  n484 n483 0   0   nm w=10u l=1u
Mp484  n484 n483 vdd vdd pm w=20u l=1u
C484   n484 0  100f

Mn485  n485 n484 0   0   nm w=10u l=1u
Mp485  n485 n484 vdd vdd pm w=20u l=1u
C485   n485 0  100f

Mn486  n486 n485 0   0   nm w=10u l=1u
Mp486  n486 n485 vdd vdd pm w=20u l=1u
C486   n486 0  100f

Mn487  n487 n486 0   0   nm w=10u l=1u
Mp487  n487 n486 vdd vdd pm w=20u l=1u
C487   n487 0  100f

Mn488  n488 n487 0   0   nm w=10u l=1u
Mp488  n488 n487 vdd vdd pm w=20u l=1u
C488   n488 0  100f

Mn489  n489 n488 0   0   nm w=10u l=1u
Mp489  n489 n488 vdd vdd pm w=20u l=1u
C489   n489 0  100f

Mn490  n490 n489 0   0   nm w=10u l=1u
Mp490  n490 n489 vdd vdd pm w=20u l=1u
C490   n490 0  100f

Mn491  n491 n490 0   0   nm w=10u l=1u
Mp491  n491 n490 vdd vdd pm w=20u l=1u
C491   n491 0  100f

Mn492  n492 n491 0   0   nm w=10u l=1u
Mp492  n492 n491 vdd vdd pm w=20u l=1u
C492   n492 0  100f

Mn493  n493 n492 0   0   nm w=10u l=1u
Mp493  n493 n492 vdd vdd pm w=20u l=1u
C493   n493 0  100f

Mn494  n494 n493 0   0   nm w=10u l=1u
Mp494  n494 n493 vdd vdd pm w=20u l=1u
C494   n494 0  100f

Mn495  n495 n494 0   0   nm w=10u l=1u
Mp495  n495 n494 vdd vdd pm w=20u l=1u
C495   n495 0  100f

Mn496  n496 n495 0   0   nm w=10u l=1u
Mp496  n496 n495 vdd vdd pm w=20u l=1u
C496   n496 0  100f

Mn497  n497 n496 0   0   nm w=10u l=1u
Mp497  n497 n496 vdd vdd pm w=20u l=1u
C497   n497 0  100f

Mn498  n498 n497 0   0   nm w=10u l=1u
Mp498  n498 n497 vdd vdd pm w=20u l=1u
C498   n498 0  100f

Mn499  n499 n498 0   0   nm w=10u l=1u
Mp499  n499 n498 vdd vdd pm w=20u l=1u
C499   n499 0  100f

.ic V(n1)=1.6 V(n2)=0.1 V(n3)=1.6 V(n4)=0.1 V(n5)=1.6 V(n6)=0.1 V(n7)=1.6 V(n8)=0.1 V(n9)=1.6 V(n10)=0.1 V(n11)=1.6 V(n12)=0.1 V(n13)=1.6 V(n14)=0.1 V(n15)=1.6 V(n16)=0.1 V(n17)=1.6 V(n18)=0.1 V(n19)=1.6 V(n20)=0.1 V(n21)=1.6 V(n22)=0.1 V(n23)=1.6 V(n24)=0.1 V(n25)=1.6 V(n26)=0.1 V(n27)=1.6 V(n28)=0.1 V(n29)=1.6 V(n30)=0.1 V(n31)=1.6 V(n32)=0.1 V(n33)=1.6 V(n34)=0.1 V(n35)=1.6 V(n36)=0.1 V(n37)=1.6 V(n38)=0.1 V(n39)=1.6 V(n40)=0.1 V(n41)=1.6 V(n42)=0.1 V(n43)=1.6 V(n44)=0.1 V(n45)=1.6 V(n46)=0.1 V(n47)=1.6 V(n48)=0.1 V(n49)=1.6 V(n50)=0.1 V(n51)=1.6 V(n52)=0.1 V(n53)=1.6 V(n54)=0.1 V(n55)=1.6 V(n56)=0.1 V(n57)=1.6 V(n58)=0.1 V(n59)=1.6 V(n60)=0.1 V(n61)=1.6 V(n62)=0.1 V(n63)=1.6 V(n64)=0.1 V(n65)=1.6 V(n66)=0.1 V(n67)=1.6 V(n68)=0.1 V(n69)=1.6 V(n70)=0.1 V(n71)=1.6 V(n72)=0.1 V(n73)=1.6 V(n74)=0.1 V(n75)=1.6 V(n76)=0.1 V(n77)=1.6 V(n78)=0.1 V(n79)=1.6 V(n80)=0.1 V(n81)=1.6 V(n82)=0.1 V(n83)=1.6 V(n84)=0.1 V(n85)=1.6 V(n86)=0.1 V(n87)=1.6 V(n88)=0.1 V(n89)=1.6 V(n90)=0.1 V(n91)=1.6 V(n92)=0.1 V(n93)=1.6 V(n94)=0.1 V(n95)=1.6 V(n96)=0.1 V(n97)=1.6 V(n98)=0.1 V(n99)=1.6 V(n100)=0.1 V(n101)=1.6 V(n102)=0.1 V(n103)=1.6 V(n104)=0.1 V(n105)=1.6 V(n106)=0.1 V(n107)=1.6 V(n108)=0.1 V(n109)=1.6 V(n110)=0.1 V(n111)=1.6 V(n112)=0.1 V(n113)=1.6 V(n114)=0.1 V(n115)=1.6 V(n116)=0.1 V(n117)=1.6 V(n118)=0.1 V(n119)=1.6 V(n120)=0.1 V(n121)=1.6 V(n122)=0.1 V(n123)=1.6 V(n124)=0.1 V(n125)=1.6 V(n126)=0.1 V(n127)=1.6 V(n128)=0.1 V(n129)=1.6 V(n130)=0.1 V(n131)=1.6 V(n132)=0.1 V(n133)=1.6 V(n134)=0.1 V(n135)=1.6 V(n136)=0.1 V(n137)=1.6 V(n138)=0.1 V(n139)=1.6 V(n140)=0.1 V(n141)=1.6 V(n142)=0.1 V(n143)=1.6 V(n144)=0.1 V(n145)=1.6 V(n146)=0.1 V(n147)=1.6 V(n148)=0.1 V(n149)=1.6 V(n150)=0.1 V(n151)=1.6 V(n152)=0.1 V(n153)=1.6 V(n154)=0.1 V(n155)=1.6 V(n156)=0.1 V(n157)=1.6 V(n158)=0.1 V(n159)=1.6 V(n160)=0.1 V(n161)=1.6 V(n162)=0.1 V(n163)=1.6 V(n164)=0.1 V(n165)=1.6 V(n166)=0.1 V(n167)=1.6 V(n168)=0.1 V(n169)=1.6 V(n170)=0.1 V(n171)=1.6 V(n172)=0.1 V(n173)=1.6 V(n174)=0.1 V(n175)=1.6 V(n176)=0.1 V(n177)=1.6 V(n178)=0.1 V(n179)=1.6 V(n180)=0.1 V(n181)=1.6 V(n182)=0.1 V(n183)=1.6 V(n184)=0.1 V(n185)=1.6 V(n186)=0.1 V(n187)=1.6 V(n188)=0.1 V(n189)=1.6 V(n190)=0.1 V(n191)=1.6 V(n192)=0.1 V(n193)=1.6 V(n194)=0.1 V(n195)=1.6 V(n196)=0.1 V(n197)=1.6 V(n198)=0.1 V(n199)=1.6 V(n200)=0.1 V(n201)=1.6 V(n202)=0.1 V(n203)=1.6 V(n204)=0.1 V(n205)=1.6 V(n206)=0.1 V(n207)=1.6 V(n208)=0.1 V(n209)=1.6 V(n210)=0.1 V(n211)=1.6 V(n212)=0.1 V(n213)=1.6 V(n214)=0.1 V(n215)=1.6 V(n216)=0.1 V(n217)=1.6 V(n218)=0.1 V(n219)=1.6 V(n220)=0.1 V(n221)=1.6 V(n222)=0.1 V(n223)=1.6 V(n224)=0.1 V(n225)=1.6 V(n226)=0.1 V(n227)=1.6 V(n228)=0.1 V(n229)=1.6 V(n230)=0.1 V(n231)=1.6 V(n232)=0.1 V(n233)=1.6 V(n234)=0.1 V(n235)=1.6 V(n236)=0.1 V(n237)=1.6 V(n238)=0.1 V(n239)=1.6 V(n240)=0.1 V(n241)=1.6 V(n242)=0.1 V(n243)=1.6 V(n244)=0.1 V(n245)=1.6 V(n246)=0.1 V(n247)=1.6 V(n248)=0.1 V(n249)=1.6 V(n250)=0.1 V(n251)=1.6 V(n252)=0.1 V(n253)=1.6 V(n254)=0.1 V(n255)=1.6 V(n256)=0.1 V(n257)=1.6 V(n258)=0.1 V(n259)=1.6 V(n260)=0.1 V(n261)=1.6 V(n262)=0.1 V(n263)=1.6 V(n264)=0.1 V(n265)=1.6 V(n266)=0.1 V(n267)=1.6 V(n268)=0.1 V(n269)=1.6 V(n270)=0.1 V(n271)=1.6 V(n272)=0.1 V(n273)=1.6 V(n274)=0.1 V(n275)=1.6 V(n276)=0.1 V(n277)=1.6 V(n278)=0.1 V(n279)=1.6 V(n280)=0.1 V(n281)=1.6 V(n282)=0.1 V(n283)=1.6 V(n284)=0.1 V(n285)=1.6 V(n286)=0.1 V(n287)=1.6 V(n288)=0.1 V(n289)=1.6 V(n290)=0.1 V(n291)=1.6 V(n292)=0.1 V(n293)=1.6 V(n294)=0.1 V(n295)=1.6 V(n296)=0.1 V(n297)=1.6 V(n298)=0.1 V(n299)=1.6 V(n300)=0.1 V(n301)=1.6 V(n302)=0.1 V(n303)=1.6 V(n304)=0.1 V(n305)=1.6 V(n306)=0.1 V(n307)=1.6 V(n308)=0.1 V(n309)=1.6 V(n310)=0.1 V(n311)=1.6 V(n312)=0.1 V(n313)=1.6 V(n314)=0.1 V(n315)=1.6 V(n316)=0.1 V(n317)=1.6 V(n318)=0.1 V(n319)=1.6 V(n320)=0.1 V(n321)=1.6 V(n322)=0.1 V(n323)=1.6 V(n324)=0.1 V(n325)=1.6 V(n326)=0.1 V(n327)=1.6 V(n328)=0.1 V(n329)=1.6 V(n330)=0.1 V(n331)=1.6 V(n332)=0.1 V(n333)=1.6 V(n334)=0.1 V(n335)=1.6 V(n336)=0.1 V(n337)=1.6 V(n338)=0.1 V(n339)=1.6 V(n340)=0.1 V(n341)=1.6 V(n342)=0.1 V(n343)=1.6 V(n344)=0.1 V(n345)=1.6 V(n346)=0.1 V(n347)=1.6 V(n348)=0.1 V(n349)=1.6 V(n350)=0.1 V(n351)=1.6 V(n352)=0.1 V(n353)=1.6 V(n354)=0.1 V(n355)=1.6 V(n356)=0.1 V(n357)=1.6 V(n358)=0.1 V(n359)=1.6 V(n360)=0.1 V(n361)=1.6 V(n362)=0.1 V(n363)=1.6 V(n364)=0.1 V(n365)=1.6 V(n366)=0.1 V(n367)=1.6 V(n368)=0.1 V(n369)=1.6 V(n370)=0.1 V(n371)=1.6 V(n372)=0.1 V(n373)=1.6 V(n374)=0.1 V(n375)=1.6 V(n376)=0.1 V(n377)=1.6 V(n378)=0.1 V(n379)=1.6 V(n380)=0.1 V(n381)=1.6 V(n382)=0.1 V(n383)=1.6 V(n384)=0.1 V(n385)=1.6 V(n386)=0.1 V(n387)=1.6 V(n388)=0.1 V(n389)=1.6 V(n390)=0.1 V(n391)=1.6 V(n392)=0.1 V(n393)=1.6 V(n394)=0.1 V(n395)=1.6 V(n396)=0.1 V(n397)=1.6 V(n398)=0.1 V(n399)=1.6 V(n400)=0.1 V(n401)=1.6 V(n402)=0.1 V(n403)=1.6 V(n404)=0.1 V(n405)=1.6 V(n406)=0.1 V(n407)=1.6 V(n408)=0.1 V(n409)=1.6 V(n410)=0.1 V(n411)=1.6 V(n412)=0.1 V(n413)=1.6 V(n414)=0.1 V(n415)=1.6 V(n416)=0.1 V(n417)=1.6 V(n418)=0.1 V(n419)=1.6 V(n420)=0.1 V(n421)=1.6 V(n422)=0.1 V(n423)=1.6 V(n424)=0.1 V(n425)=1.6 V(n426)=0.1 V(n427)=1.6 V(n428)=0.1 V(n429)=1.6 V(n430)=0.1 V(n431)=1.6 V(n432)=0.1 V(n433)=1.6 V(n434)=0.1 V(n435)=1.6 V(n436)=0.1 V(n437)=1.6 V(n438)=0.1 V(n439)=1.6 V(n440)=0.1 V(n441)=1.6 V(n442)=0.1 V(n443)=1.6 V(n444)=0.1 V(n445)=1.6 V(n446)=0.1 V(n447)=1.6 V(n448)=0.1 V(n449)=1.6 V(n450)=0.1 V(n451)=1.6 V(n452)=0.1 V(n453)=1.6 V(n454)=0.1 V(n455)=1.6 V(n456)=0.1 V(n457)=1.6 V(n458)=0.1 V(n459)=1.6 V(n460)=0.1 V(n461)=1.6 V(n462)=0.1 V(n463)=1.6 V(n464)=0.1 V(n465)=1.6 V(n466)=0.1 V(n467)=1.6 V(n468)=0.1 V(n469)=1.6 V(n470)=0.1 V(n471)=1.6 V(n472)=0.1 V(n473)=1.6 V(n474)=0.1 V(n475)=1.6 V(n476)=0.1 V(n477)=1.6 V(n478)=0.1 V(n479)=1.6 V(n480)=0.1 V(n481)=1.6 V(n482)=0.1 V(n483)=1.6 V(n484)=0.1 V(n485)=1.6 V(n486)=0.1 V(n487)=1.6 V(n488)=0.1 V(n489)=1.6 V(n490)=0.1 V(n491)=1.6 V(n492)=0.1 V(n493)=1.6 V(n494)=0.1 V(n495)=1.6 V(n496)=0.1 V(n497)=1.6 V(n498)=0.1 V(n499)=1.6
.options method=gear
.tran 150p 748n UIC
