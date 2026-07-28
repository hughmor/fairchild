# giona RNN — physics to recurrence

Symbols and every quoted number are from `mrm.sp` / `neuron_junction_wdm8.sp` /
`giona_rnn_perfectW.sp`, verified against `fairchild` solves.

## 1. Notation

| symbol | meaning | value |
|---|---|---|
| $N$ | neurons = WDM channels | 8 |
| $I_D$ | modulator terminal (junction) current | state variable |
| $I_{\rm arc}$ | current per ring arc, $=I_D/2$ | two junctions in parallel |
| $V$ | $V(\texttt{mod\_cathode})$; junction forward voltage $=-V$ | |
| $I_s,\,n$ | diode saturation current, ideality (per arc) | $5.099\times10^{-8}$ A, 5.0 |
| $V_T$ | $kT/q$ at 300.15 K | $25.865$ mV |
| $R_1$ | on-chip shunt, in parallel with the modulator | 2 kΩ |
| $R_b$ | off-chip bias resistor | 10 kΩ |
| $R_p$ | each PD shunt path, $17\,{\rm k}+1+50$ | 17.051 kΩ |
| $P_{\pm}$ | PD supplies | $\pm3$ V |
| $\mathcal{R}$ | PD responsivity | 0.8 A/W |
| $a$ | tap + 1:8 tree attenuation | $\tfrac12\cdot\tfrac18=\tfrac1{16}$ |
| $L$ | ring circumference $2\pi r$, $r=8\,\mu$m | 50.27 µm |
| $n_g,\,n_{\rm eff}$ | group / effective index | 4.2, 2.302 |
| $\kappa L$ | coupler strength, both couplers | 0.183 |
| $\alpha_0$ | background loss, 10.7 dB/cm | 246.4 Np/m |
| $\partial n/\partial I$ | injection index coefficient (per arc) | 3.99 A$^{-1}$ |
| $\partial\alpha/\partial I$ | injection loss coefficient (per arc) | $4.63\times10^{6}$ Np m$^{-1}$A$^{-1}$ |
| $\partial n/\partial V$ | depletion index coefficient | $-3.62\times10^{-5}$ V$^{-1}$ |
| $P_\pi^{\rm th},\,R_h$ | heater $\pi$-power (whole ring), heater R per arc | 26.4 mW, 184.4 Ω |
| $\tau_c$ | carrier lifetime | 10 ns |

## 2. The recurrence

Target form:

$$\tau\,\dot{x}_i = -x_i + \varphi\!\Big(\sum_j W_{ij}x_j + b_i\Big),\qquad i=1\dots N .$$

Physical realisation: $x_i \equiv I_{D,i}$, $\varphi$ is the ring's
current-to-transmission response $T(\cdot)$, $W$ is the programmed weight matrix
(one $\texttt{fc\_optical\_2x2}$ per neuron, $w=0$, $\mathrm{d}w/\mathrm{d}V=1$,
so $V(W_{ij})$ *is* $W_{ij}$), and $b_i$ is a bias **current** set by $\mathrm{PD\_B}_i$.

The loop, per neuron:

$$
I_{D,i}\ \xrightarrow{\ \Delta n,\Delta\alpha\ }\ T_i
\ \xrightarrow{\ P_0 \ }\ P_i^{\rm bus}
\ \xrightarrow{\ a\,W_{ij}\ }\ P^{\rm blk}_{ij}
\ \xrightarrow{\ \mathcal{R}\ }\ I_{{\rm sig},i}
\ \xrightarrow{\ \eta\ }\ I_{D,i}
$$

giving the fixed-point equation

$$\boxed{\;I_{D,i} \;=\; I_{b,i} \;+\; \eta\,\kappa \sum_j W_{ij}\,T_j(I_{D,j}),
\qquad \kappa \equiv \mathcal{R}\,a\,P_0 \;}\tag{1}$$

$\kappa$ has units of amps: the photocurrent a fully-transmitting channel
delivers at unit weight.

## 3. Bias: the row-sum rule

Demand that **all** neurons rest at the same $I_D=I^\ast$. Then every
$T_j(I^\ast)=T^\ast$ is a known constant and (1) collapses:

$$I^\ast = I_{b,i} + \eta\,\kappa\,T^\ast \underbrace{\textstyle\sum_j W_{ij}}_{\displaystyle s_i}
\qquad\Longrightarrow\qquad
\boxed{\;I_{b,i} \;=\; I^\ast - \eta\,\kappa\,T^\ast s_i\;}\tag{2}$$

with $s_i$ the $i$-th **row sum** of $W$. No iteration on the recurrence is
needed — choosing the rest state for all neurons at once makes the feedback term
a constant.

Consequences:

- $s_i = 0\ \forall i$ $\Rightarrow$ $I_{b,i}=I^\ast$: one bias for every
  neuron, independent of $W$. Reprogramming the weights leaves the bias alone.
- $W\mathbf{1}=0 \Rightarrow \det W = 0$. For $N=2$ the spectrum is
  $\{0,\operatorname{tr}W\}$, both real, so **a zero-row-sum 2-neuron network
  cannot oscillate** (§7).

## 4. Bias network $\to$ $\mathrm{PD\_B}$

KCL at $\texttt{mod\_cathode}$, junction anode grounded, $I_{\rm sig}$ drawn out
of the node:

$$I_D(V) \;=\; \frac{V}{R_1} + \frac{2V}{R_p} - \frac{\mathrm{PD\_B}-V}{R_b} + I_{\rm sig}
\;=\; \frac{V}{R_{\rm sh}} - \frac{\mathrm{PD\_B}}{R_b} + I_{\rm sig}\tag{3}$$

$$\frac{1}{R_{\rm sh}} = \frac{1}{R_1} + \frac{1}{R_b} + \frac{2}{R_p}
\qquad R_{\rm sh} = 1394\ \Omega$$

The two $\pm3$ V PD paths appear as $2V/R_p$: a balanced pair cancelling only at
$V=0$, so at a forward-biased node they are a further $R_p/2 = 8.53$ kΩ to
ground. Omitting them understates the shunt by 20 %.

Two junctions in parallel, so the terminal diode law is

$$I_D = 2I_s\Big(e^{-V/nV_T}-1\Big)
\qquad\Longleftrightarrow\qquad
V = -\,nV_T\ln\!\Big(1+\frac{I_D}{2I_s}\Big)\tag{4}$$

Setting $I_{\rm sig}=0$ in (3) and substituting (4) gives the bias in closed
form:

$$\boxed{\;\mathrm{PD\_B} \;=\; R_b\left(\frac{V^\ast}{R_{\rm sh}} - I^\ast\right),
\qquad V^\ast = -nV_T\ln\!\left(1+\frac{I^\ast}{2I_s}\right)\;}\tag{5}$$

Verified — (4) and (5) reproduce the solver to the printed precision:

| $\mathrm{PD\_B}$ | $I^\ast$ (µA) | $V^\ast$ measured | $V^\ast$ from (4) | (5) recovers |
|---|---|---|---|---|
| −5.0 V | 18.3 | −671.61 mV | −671.61 mV | −5.00 V |
| −8.0 V | 133.9 | −928.65 mV | −928.65 mV | −8.00 V |
| −10.0 V | 269.1 | −1018.90 mV | −1018.90 mV | −10.00 V |
| −13.4 V | 543.9 | −1109.86 mV | −1109.86 mV | −13.40 V |
| −17.8 V | 933.8 | −1179.74 mV | −1179.74 mV | −17.80 V |
| −22.0 V | 1321.6 | −1224.66 mV | −1224.66 mV | −22.00 V |

Combining (2) and (5) is the full bias recipe: pick $I^\ast$, read $T^\ast$ off
the activation, compute $I_{b,i}$ from the row sums, convert each to
$\mathrm{PD\_B}_i$.

## 5. Small-signal division at the node ($\eta$)

Differentiate (3) at fixed $\mathrm{PD\_B}$, with
$\mathrm{d}I_D/\mathrm{d}V = -1/r_d$ and $r_d = nV_T/I_D$:

$$-\frac{\mathrm{d}V}{r_d} = \frac{\mathrm{d}V}{R_{\rm sh}} + \mathrm{d}I_{\rm sig}
\;\Longrightarrow\;
\mathrm{d}V = -\,\mathrm{d}I_{\rm sig}\,(r_d\!\parallel\!R_{\rm sh})$$

$$\boxed{\;\eta \equiv \frac{\mathrm{d}I_D}{\mathrm{d}I_{\rm sig}}
= \frac{R_{\rm sh}}{R_{\rm sh}+r_d},\qquad r_d = \frac{nV_T}{I_D}\;}\tag{6}$$

Measured against (6):

| $I^\ast$ (µA) | 18.3 | 133.9 | 269.1 | 543.9 | 933.8 | 1321.6 |
|---|---|---|---|---|---|---|
| $r_d$ (Ω) | 7065 | 966 | 481 | 238 | 138 | 98 |
| $\eta$ predicted | 16.4 % | 59.1 % | 74.4 % | 85.4 % | 91.0 % | 93.4 % |
| $\eta$ measured | 16.6 % | 59.1 % | 74.4 % | 85.4 % | 91.0 % | 93.4 % |

$\eta\to1$ as $I^\ast$ rises: the stiffer the diode, the less signal current the
shunts steal. This is the first half of the rest-point trade-off.

**On the $2$ kΩ.** $R_1$ is across the modulator, so current into it is not inert
— it develops the junction *field*. Driven far enough negative that is
depletion-mode plasma dispersion, a different modulation mechanism, opposite in
sign. At the operating point it is a small correction:

$$\frac{|\partial n/\partial V\cdot V^\ast|}{|\partial n/\partial I\cdot I_{\rm arc}^\ast|}
= \frac{4.02\times10^{-5}}{1.085\times10^{-3}} = 3.7\ \%
\qquad (\mathrm{PD\_B}=-13.4\ {\rm V})$$

so §6 keeps only the injection term.

## 6. Junction current $\to$ optical response

Each arc carries $I_{\rm arc}=I_D/2$ over length $L/2$, so the round trip sees
the arc coefficients once over the full $L$:

$$\Delta n_{\rm eff} = -\frac{\partial n}{\partial I}\,\frac{I_D}{2},
\qquad
\Delta\alpha = \frac{\partial\alpha}{\partial I}\,\frac{I_D}{2}$$

Phase and resonance shift:

$$\phi = \frac{2\pi n_{\rm eff}L}{\lambda},\qquad
\frac{\mathrm{d}\phi}{\mathrm{d}I_D} = -\frac{\pi L}{\lambda}\frac{\partial n}{\partial I}
= -407.5\ {\rm rad/A}$$

$$\frac{\delta\lambda_{\rm res}}{\lambda} = \frac{\Delta n_{\rm eff}}{n_g}
\;\Longrightarrow\;
\boxed{\;\frac{\mathrm{d}\lambda_{\rm res}}{\mathrm{d}I_D}
= -\frac{\lambda}{2n_g}\frac{\partial n}{\partial I} = -0.734\ {\rm pm/\mu A}\;}\tag{7}$$

Injection blue-shifts; the heater red-shifts, with whole-ring power
$P_h = I_h^2\,(2R_h)$ and $\Delta\phi_{\rm th} = \pi P_h/P_\pi^{\rm th}$.
Cross-check at $\mathrm{PD\_B}=-13.4$ V, where the measured re-centring current
is $I_h = 2.29$ mA:

$$\delta\lambda_{\rm inj} = -399\ {\rm pm},\qquad
\delta\lambda_{\rm th} = +415\ {\rm pm},\qquad
\text{sum } +15\ {\rm pm}\ (4\ \%)$$

Two independent coefficient chains cancelling to 4 % is the check that (7) and
the thermal model are consistent.

Add-drop response, $t=\cos\kappa L$, $k^2=1-t^2$, amplitude round trip
$A=e^{-\alpha L/2}$ with $\alpha = \alpha_0 + (\partial\alpha/\partial I)I_{\rm arc}$:

$$T_{\rm notch} = \frac{t^2A^2 - 2t^2A\cos\phi + t^2}{1-2t^2A\cos\phi + t^4A^2},
\qquad
T_{\rm peak} = \frac{k^4 A}{1-2t^2A\cos\phi+t^4A^2}\tag{8}$$

`in→thru` and `add→drop` are notches; `in→drop` and `add→thru` are peaks.

Finesse $\mathcal{F}\approx \pi/(1-t^2A)$ — injection loss destroys it:

| $I_D$ | $\alpha L$ | round-trip loss | $t^2A$ | $\mathcal{F}$ |
|---|---|---|---|---|
| 0 | 0.0124 Np | 1.23 % | 0.9609 | 80.4 |
| 133.9 µA | 0.0280 Np | 2.76 % | 0.9535 | 67.5 |
| 543.9 µA | 0.0757 Np | 7.29 % | 0.9310 | 45.5 |

This is the second half of the trade-off: free-carrier absorption at high $I^\ast$
washes out the very resonance the modulation relies on. Measured peak
transmission falls $0.71\to0.10$ mW over the bias range while $\eta$ climbs
$17\to93$ %; the product $|\mathrm{d}T/\mathrm{d}I_{\rm sig}|$ peaks at
$I^\ast\approx134\ \mu$A ($\mathrm{PD\_B}=-8$ V).

Natural input scale — one linewidth of junction current:

$$I_{\rm lw} = \frac{\delta\lambda_{\rm FWHM}}{|\mathrm{d}\lambda_{\rm res}/\mathrm{d}I_D|}
= \frac{\mathrm{FSR}/\mathcal{F}}{0.734\ {\rm pm/\mu A}}
\approx \frac{11.32\ {\rm nm}/80}{0.734\ {\rm pm/\mu A}} \approx 192\ \mu{\rm A}$$

One FSR is $\approx15.4$ mA, so the $\pm800\ \mu$A sweeps in
`rnn_explore.py --space` span $\sim4$ linewidths and $5\ \%$ of an FSR.

## 7. Loop gain and eigenvalues

Linearise (1) about the rest point, $u_i = I_{D,i}-I^\ast$:

$$u_i = \eta\,\kappa\,T'(I^\ast)\sum_j W_{ij}u_j
\qquad\Longrightarrow\qquad
u = G\,W u,\qquad
\boxed{\;G \equiv \eta\,\kappa\,T'(I^\ast) = \eta\,\mathcal{R}\,a\,P_0\,\frac{\mathrm{d}T}{\mathrm{d}I_D}\;}\tag{9}$$

$G$ is dimensionless: the round-trip small-signal gain at unit weight.
Dynamics, with the dominant pole $\tau$ (§8):

$$\tau\,\dot{u} = -u + G\,W u
\qquad\Longrightarrow\qquad
\mu_k = \frac{-1+G\lambda_k}{\tau},\quad \lambda_k \in \operatorname{spec}W$$

$$\boxed{\;\text{instability} \iff \operatorname{Re}(G\lambda_k) > 1 \ \text{for some }k\;}\tag{10}$$

(The discrete map $x_{t+1}=\varphi(Wx_t+b)$ instead requires $|G\lambda_k|>1$.)
A complex pair crossing with $\operatorname{Im}(G\lambda_k)\neq0$ is a Hopf
bifurcation — oscillation; a real crossing is a saddle-node — bistability.

| $W$ | $\lambda$ | threshold, (10) | behaviour |
|---|---|---|---|
| $\begin{pmatrix}1&-1\\1&1\end{pmatrix}$ | $1\pm i$ | $G>1$ | Hopf $\Rightarrow$ oscillator |
| $\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$ | $2,\ 0$ | $G>\tfrac12$ | real $\Rightarrow$ WTA / bistable |
| $\begin{pmatrix}0&-1\\1&0\end{pmatrix}$ | $\pm i$ | none | $\operatorname{Re}(G\lambda)=0$: single pole never destabilises |
| any $s_i\equiv0$, $N=2$ | $0,\ \operatorname{tr}W$ | $G\operatorname{tr}W>1$ | real only $\Rightarrow$ no oscillation |

Two structural facts worth keeping separate:

- **Row sums set the bias** (2) — a *DC* statement.
- **Eigenvalues set the required gain** (10) — an *AC* statement.

They collide in the zero-row-sum case, which is attractive for biasing
($I_{b,i}$ weight-independent) and useless for 2-neuron oscillation
(spectrum real). Three neurons remove the conflict: $W\mathbf{1}=0$ leaves a
$2\times2$ block on $\mathbf{1}^\perp$ that can carry a complex pair.

## 8. Time constants

| pole | expression | value |
|---|---|---|
| carrier lifetime | $\tau_c$ | 10 ns |
| node $RC$ | $R_{\rm sh}(2C_{j0}+C_j)$, $C_{j0}=0.1375$ pF | 0.40 ns |
| photon lifetime | $Q\lambda/2\pi c$, $Q=\lambda\mathcal{F}/\mathrm{FSR}\approx1.1\times10^{4}$ | 9 ps |

$\tau\approx\tau_c$ dominates. Hopf frequency at onset:

$$f = \frac{\operatorname{Im}(G\lambda_k)}{2\pi\tau}
\;\overset{W=[[1,-1],[1,1]],\,G\approx1}{\approx}\; \frac{1}{2\pi\cdot10\ {\rm ns}} \approx 16\ {\rm MHz}$$

An order below the 200–400 MHz weight drives used in `rnn_wta.py` — those
experiments probe the forced response, not the natural mode.

## 9. Measuring $G$

Two independent routes; both need the ring **trimmed**, because $G\propto \mathrm{d}T/\mathrm{d}I_D$
is zero at a resonance peak and changes sign across it.

**(a) Closed-loop bias sensitivity.** With self-coupling $W_{11}=w$ and
$S \equiv \partial V(\texttt{mod\_cathode}_1)/\partial \mathrm{PD\_B}_1$, one
round trip adds $wG$:

$$u = S_0\,\delta\mathrm{PD\_B} + wG\,u
\;\Longrightarrow\;
S(w) = \frac{S_0}{1-wG}
\;\Longrightarrow\;
\boxed{\;G = \frac{1}{w}\left(1-\frac{S_0}{S(w)}\right)\;}\tag{11}$$

Requirements, each of which produced a wrong answer when violated:

- Perturb **one** neuron's bias. A common-mode $\delta\mathrm{PD\_B}$ moves all
  eight rings on the shared bus and folds crosstalk into $G$.
- $w$ large enough that $wG$ exceeds the solver's resolution on $S$. At
  $w=0.05$, $S_0$ and $S(w)$ differed in the 6th decimal and (11) returned
  $\sim10^{-3}$ regardless of bias — noise, not gain.
- $w$ small enough to stay small-signal. At $w=0.35$, $w>0$ and $w<0$ gave
  $-0.95$ and $+2.55$.
- $w\in[0.1,0.3]$ works: $G$ varies $<10\%$ over that range, which is the test
  that it *is* a small-signal quantity.

**(b) Direct transconductance.** Sweep a real weight, measure output power
against input photocurrent, and assemble (9) from measured factors:

$$G = a\,\mathcal{R}\,\frac{\mathrm{d}P^{\rm bus}}{\mathrm{d}I_{\rm in}}$$

No division by a small number, and every factor is separately checkable.

**Measured** (full 8-neuron deck, $P_0=30$ mW/channel, heater trimmed per bias
for maximum $|G|$):

| $\mathrm{PD\_B}$ | heater | $G$ from (11) | $G$ from (b) | $2G$ vs WTA threshold 1 |
|---|---|---|---|---|
| −5.0 V | 0.81 mA | — | 1.16 | 2.3 |
| −8.0 V | 0.81 mA | 2.06 | 3.51 | 7.0 |
| −10.0 V | 1.31 mA | 2.30 | 3.71 | 7.4 |
| −13.4 V | 1.94 mA | 1.80 | 2.89 | 5.8 |
| −17.8 V | 2.69 mA | — | 1.89 | 3.8 |

Caveats, stated rather than smoothed over:

- The two methods agree in sign and within $\sim1.7\times$ in magnitude. They
  are not the same estimator — (11) is a derivative at one operating point,
  (b) a chord over a weight range — so a single $G$ is only meaningful with the
  trim and the probe width quoted.
- $|G|$ is $O(1)$ across the whole practical bias range once trimmed, so
  $\operatorname{Re}(G\lambda)>1$ is reachable for both $W$ above. Bias choice
  is an efficiency question, not a feasibility one.
- `giona_rnn_perfectW.sp` hardcodes `Iht<k> ... DC 0`. Every WTA and AC-weight
  run to date was therefore **untrimmed**, where $|G|\approx0.24$–$1.0$ and the
  sign is opposite to the trimmed case. Trimming the heaters is the first thing
  to change before reading anything else into those transients.

## 10. Provenance

| quantity | produced by |
|---|---|
| (4), (5), (6), §6 tables | `rnn_explore.py --space` |
| $\eta$, $\mathrm{d}T/\mathrm{d}I_{\rm sig}$, rest-point trade-off | `rnn_explore.py --space` |
| $G$ from (b) | `rnn_explore.py --insitu` |
| $G$ from (11) | single-neuron $\mathrm{PD\_B}$ perturbation, §9 |
| radii, $\mathrm{d}\lambda/\mathrm{d}r$ | `rnn_explore.py --radii` |
| bias solve for a given $I^\ast$ | `rnn_drive.solve_bias` (coupled Newton on $\partial V_i/\partial\mathrm{PD\_B}_j$) |

Device coefficients are the `mrm.sp` defaults, fitted in
`experiments/giona/ringfit.py`; the $i_s$/$\partial n/\partial I$/$\partial\alpha/\partial I$
trio is only meaningful as a set and is pending a refit from a clean on-die IV.
