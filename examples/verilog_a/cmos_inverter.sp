* CMOS inverter with Verilog-A transistors — the foundry-PDK idiom.
*
* This is the shape every real PDK deck has, and the reason `.model` support
* for OSDI matters: you get a compiled `.osdi` from the foundry, name a card
* after its module, put the process parameters on the card, and put the
* geometry on the instance line.
*
*   .osdi  <module>.osdi
*   .model <card> <module> (<process params>)
*   M1 d g s b <card> W=… L=…
*
* Swap va_nmos/va_pmos for a real bsim4.osdi and the deck does not change
* shape.  Card params are model-level; instance params are per-device and win
* over the card, exactly as SPICE expects.
*
* Run:
*   fairchild -f examples/verilog_a/cmos_inverter.sp \
*             --probe "v(in),v(out)" --format csv
*   fairchild -f examples/verilog_a/cmos_inverter.sp \
*             --probe "v(out)" --format csv --opt itl1=200   # (DC transfer)

.osdi build/va_nmos.osdi
.osdi build/va_pmos.osdi

* Process cards.  KP is silicon; W/L is layout, so it lives on the instances.
.model nch va_nmos (KP=120u VTH0=0.7  LAMBDA=0.02 CGSO=200p CGDO=200p)
.model pch va_pmos (KP=40u  VTH0=-0.8 LAMBDA=0.02 CGSO=200p CGDO=200p)

Vdd  vdd 0   DC 3.3
Vin  in  0   PULSE(0 3.3 1n 200p 200p 4n 10n)

* 3x wider PMOS to compensate its lower KP — the usual beta-matching.
Mp   out in vdd vdd pch W=30u L=1u
Mn   out in 0   0   nch W=10u L=1u

Cload out 0 50f

.options method=gear
.tran 20p 20n
