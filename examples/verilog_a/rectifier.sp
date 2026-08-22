* Half-wave rectifier: Verilog-A diode + native SPICE R, C, V.
*
* The only non-native device here is Xd1, whose model comes from
* models/va_diode.va compiled to build/va_diode.osdi by ./build.sh.
* Everything else is a fairchild built-in, and they all share one
* Newton-Raphson loop — there is no co-simulation boundary.
*
*   Vin ~ 5 V, 1 kHz ──[Rsrc 50]── Xd1 ─┬── Cload 1u ── 0
*                                       └── Rload 10k ── 0
*
* Note the `X` prefix.  An OSDI model can also be instantiated as `D1 a b
* va_diode`, but the D/M/Q element parsers stop reading at the model name, so
* instance parameters on those lines are silently discarded — and a `.model`
* card does not reach an OSDI device either.  `X` is the only form that
* parameterises one.
*
* Expect V(out) to charge to about (5 - Vf) on the first positive half
* cycle and then droop with tau = Rload*Cload = 10 ms between peaks.
*
* Run:
*   fairchild -f examples/verilog_a/rectifier.sp \
*             --probe "v(in),v(out)" --format csv -o /tmp/rect.csv

.osdi build/va_diode.osdi

Vin   in  0   SIN(0 5 1k)
Rsrc  in  a   50
Xd1   a   out va_diode Is=1e-14 N=1.0 Rs=0.5 Cj0=2p
Cload out 0   1u
Rload out 0   10k

.options method=gear
.tran 5u 5m
