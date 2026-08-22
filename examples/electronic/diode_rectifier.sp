* Half-wave rectifier — D1N4148-like diode + 1kΩ load
* Full Shockley model: IS, N, RS (series resistance), CJO (junction cap),
* VJ, MJ, FC (depletion cap shape), TT (transit time diffusion charge).
* Vf ≈ 0.6 V at ~1 mA matches D1N4148 datasheet.
* CJO+TT produce capacitive transients at switching edges (negative dip on
* falling edge, early rise on rising edge) matching ngspice/LTspice behaviour.
.model D1N4148 D (IS=2.52n N=1.752 RS=0.568 CJO=4p VJ=0.5 MJ=0.37 FC=0.5 TT=5n)
V1  in  0  PULSE(-2 2 0 100n 100n 400n 1u)
D1  in  out  D1N4148
R1  out 0    1k
.tran 10n 3u
