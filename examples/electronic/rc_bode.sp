* RC low-pass — Bode plot demo (AC sweep).
*
* Analytic: f_c = 1/(2π·R·C) = 159.155 Hz with R = 1 kΩ, C = 1 µF.
* Magnitude at f_c is 1/√2 ≈ 0.7071; phase is −45°.
*
* Try:
*   fairchild -f rc_bode.sp --format csv -o /tmp/bode.csv
*   python -c "import numpy as np, matplotlib.pyplot as plt; \
*              d = np.loadtxt('/tmp/bode.csv', delimiter=',', skiprows=1); \
*              plt.semilogx(d[:,0], 20*np.log10(d[:,1])); plt.show()"

V1 in 0 DC 1 AC 1
R1 in out 1k
C1 out 0 1u
.ac dec 50 1 1Meg
