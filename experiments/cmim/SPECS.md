# CMIM specification card

Every number the models in this directory use, with the place it came from.
The source is Al-Kayed et al., *Programmable 200 GOPS Hopfield-inspired photonic
Ising machine*, Nature **648**, 576 (2025), doi:10.1038/s41586-025-09838-7, plus
its Supplementary Information.

Citations use these tags:

| Tag | Means |
|---|---|
| `M` | Methods section of the main article |
| `EDF n` / `EDT n` | Extended Data Figure / Table n |
| `S n` | Supplementary Information section n |
| **ASSUMED** | Not in the paper. A stated default, with the reason. |
| **GAP** | Not in the paper and it matters. Listed again in [Open gaps](#open-gaps). |

A parameter marked ASSUMED or GAP must never be reported as a measured value.

---

## 1. Optical source

Quantum-well DFB laser, O-band.

| Quantity | Value | Source |
|---|---|---|
| Operating output power | 6.5 dBm (4.47 mW) | S1.2 |
| Rated output power | 16 dBm at 1300 nm, 25 °C | S1.2 |
| Operating wavelength | 1310 nm | M |
| Wavelength at 15 °C | 1297.3 nm | S1.2, Fig. S2c |
| Wavelength at 45 °C | 1301.1 nm | S1.2, Fig. S2c |
| Wavelength drift | 0.127 nm/K (derived from the two points above) | derived |
| Power drift | 1.4 dB over 15 °C to 45 °C | S1.2 |
| Operating temperature | 20 °C, thermoelectrically controlled | S1.2 |
| Threshold current | 8 mA | S4.2 |
| Slope efficiency | 0.352 W/A | S4.2 |
| Bias voltage | 1.6 V | S4.2 |
| Relative intensity noise | **GAP** | S1.6 names RIN as a noise source and gives no value |

The two wavelength points disagree with the "1310 nm operating" statement in
Methods by about 10 nm. The laser is tuned per experiment. Treat 1310 nm as the
operating point and the S2c points as the temperature slope only.

## 2. TFLN Mach-Zehnder modulators

Two cascaded devices from HyperLight on one photonic chip. Each has RF
electrodes and one integrated thermal phase shifter on a single arm.

| Quantity | Value | Source |
|---|---|---|
| DC half-wave voltage | 1.5 V at low MHz | S1.1, Table S1, EDF 1 |
| Electro-optic 3 dB bandwidth | > 110 GHz | S1.1, Table S1 |
| DC extinction ratio | 42 dB | Table S1 |
| Electrode length | 1.4 cm | Table S1 |
| Half-wave voltage-length product | 2.1 V·cm | Table S1 |
| RF half-wave voltage at 100 GHz | about 4.5 V | Fig. S1e, right axis |
| RF drive at 100 GHz, quasi-linear | 2.7 Vpp, 0.955 Vrms | S4.2 |
| RF termination | 50 Ω on chip | S4.2 |
| Insertion loss, device | 2 dB | S4.2 |
| Coupling loss, chip | 4 dB | S4.2 |
| Quadrature bias loss | 3 dB | S4.2 |
| Total loss per modulator | 9 dB | S4.2 |
| Thermal phase shifter power | 30 mW average, so about 60 mW for π | S4.2 |
| RF power dissipated | 18.25 mW at 2.7 Vpp | S4.2 |

`Vπ(f)` rises from 1.5 V at DC to about 4.5 V at 100 GHz (Fig. S1e). The S21
curve on the same figure runs from about +2 dB to about −4 dB over 0 to 100 GHz.
The model in `models/va_tfln_mzm.va` reproduces this with a one-pole roll-off
plus a square-root-of-frequency microwave loss term. See
[Open gaps](#open-gaps) item 2.

**Roles.** The first modulator carries the spin vector `x(t)` and is driven hard
enough to span the full nonlinear half-wave. The second modulator carries the
flattened weight vector `w(t)` and is driven small-signal, so it stays linear.
(M, "Experimental test bed".)

## 3. Quantum dot SOA

Innolume GmbH. Signal amplifier, temperature controlled at 20 °C.

| Quantity | Value | Source |
|---|---|---|
| Small-signal gain | 18.2 dB | S1.1 |
| Input saturation power | 4.1 dBm | S1.1 |
| Noise figure at −5 dBm input | 5.6 dB | S1.1 |
| Operating input power | −13 dBm | S1.1 |
| Optical 3 dB bandwidth | 27.3 nm | S1.1 |
| Gain peak wavelength | about 1300 nm at 20 °C | Fig. S1c |
| Linewidth enhancement factor | 0.86 to 1.2 | S1.1 |
| Bias | 700 mA at 1.7 V, 1.19 W | S4.2 |
| ASE power spectral density | 1.513 fW/Hz | S4.2 |
| Polarisation correction factor | 2 | S4.2 |
| Gain peak at 15 °C | 17.9 dB at 1293 nm | S1.2 |
| Gain peak at 45 °C | 14.6 dB at 1310 nm | S1.2 |
| Gain peak drift | −0.11 dB/K, +0.567 nm/K | derived |

Figure S1b gives gain compression against input power. Gain falls from about
18.2 dB at −25 dBm to about 14 dB at +5 dBm. Noise figure rises from about
4.5 dB to about 7 dB across the same range. `models/va_qd_soa.va` fits both with
a standard saturation form. Section S4.2 quotes 18 dB and 5 dB for the same
device, rounded. The model uses the S1.1 numbers.

## 4. Bulk SOA, used as a controlled noise source

The bulk SOA injects amplified spontaneous emission through a 50:50 fibre
coupler ahead of the QD SOA. It is the annealing knob.

| Quantity | Value | Source |
|---|---|---|
| Bias range | 0 mA to 110 mA | M, "Controllable optical noise source" |
| Noise distribution at the scope | Gaussian | EDF 2b |
| Best constant bias, 400-node lattice | 60 mA | EDF 2d |
| Annealing schedule | `η² = η₀² · exp(−γ·t)` | M |
| Initial noise variance `η₀²` | 18.2 mV² | M |
| Decay rate `γ` | 0.0235 per iteration | M |
| Bias at the start of the schedule | 100 mA | EDF 2e |
| Iterations to ground state, no noise | 145 ± 15 | M |
| Iterations to ground state, 60 mA | 120 ± 30 | M |
| Iterations to ground state, schedule | 100 ± 25 | M |
| Electrical power | about 100 mW | S4.2 |
| Noise variance against bias current | **GAP** | Only the plot in EDF 2a |

The units of `η₀²` are mV² at the scope, not spin units. The mapping from scope
millivolts to the normalised spin amplitude is not given. See
[Open gaps](#open-gaps) item 4.

## 5. Photodetector and receiver

| Quantity | Value | Source |
|---|---|---|
| Photodetector 3 dB bandwidth | 100 GHz | M, EDF 1 |
| Photodetector bias power | 1 mW | S4.2 |
| Responsivity | 0.8 A/W | **ASSUMED**, typical InGaAs PIN at 1310 nm |
| Dark current | 10 nA | **ASSUMED**, typical for the class |
| Receiver amplifier | **GAP** | Not named in the paper |

Section S4.2 quotes a 45.5 GHz, 22 nm CMOS, 11.2 mW transimpedance amplifier
with 2.7 µA input-referred noise. That is a projection for a scaled integrated
system and **not** the laboratory receiver. `models/va_rx_amp.va` uses it as a
stated placeholder. See [Open gaps](#open-gaps) item 3.

## 6. Full electro-optic link response

Measured end to end with transmitter pre-emphasis applied (Fig. S5a).

| Quantity | Value | Source |
|---|---|---|
| 3 dB point | about 55 GHz | S1.5.1 |
| 6 dB point | about 75 GHz | S1.5.1 |
| Attenuation at 64 GHz | about 5 dB | S1.5.1 |
| Effective number of bits at DC | 6 bits | S1.6 |
| Effective number of bits at 100 GHz | 5.2 bits | S1.6 |

Section S1.6 prints these two ENOB figures with the unit "dB". ENOB is measured
in bits. The model reads them as bits.

This response is the single most useful calibration target in the paper. It is
the product of the driver, the modulator, the photodetector, the receiver and
the two data converters. `link/link_response.py` compares the deck against it.

## 7. Instruments

### Arbitrary waveform generator, Keysight M8199B

| Quantity | Value | Source |
|---|---|---|
| Channels used | 2 | M |
| Sample rate | 256 GSa/s | M |
| Memory | 2²⁰ samples per channel | S1.3 |
| Maximum symbols at 106 GBaud | 434,176 at 2.415 samples per symbol | S1.3 |
| Upload time | `4.4981e-6 · N_samples + 0.1628` s | S1.3, eq. 3 |
| Analogue bandwidth | 65 GHz | **ASSUMED**, M8199B data sheet |
| Vertical resolution | 8 bits nominal | **ASSUMED**, M8199B data sheet |
| Effective number of bits | 6 at DC, 5.2 at 100 GHz | S1.6, whole system |

Channel 1 carries the spins. Channel 2 carries the flattened weight matrix. The
relative delay between the two channels is trimmed in software so that the
element-wise product lines up. The optical travel time between the two
modulators, `ΔT`, sets the nominal offset (M, eq. 5).

### Real-time oscilloscope, Keysight UXR

| Quantity | Value | Source |
|---|---|---|
| Sample rate | 256 GSa/s | M |
| Capture time | 0.81 ± 0.02 s for 2²² samples | S1.3 |
| Analogue bandwidth | 110 GHz | **ASSUMED**, UXR series data sheet |
| Vertical resolution | 10 bits nominal | **ASSUMED**, UXR series data sheet |

### Digital signal processing host

Intel i7-8700K. Offline, not pipelined.

| Quantity | Value | Source |
|---|---|---|
| Processing time, 65,536 symbols at 106 GBaud | 24.47 ms | S1.3 |
| Maximum measured processing time | 58.28 ms for 606,308 symbols at 148 GBaud | S1.3 |
| Projected pipelined latency | 230 ns, 115-deep pipeline at 500 MHz | S1.3 |

## 8. Digital signal processing chain

Both stacks come from Methods, "Experimental test bed". The tap counts come
from the pipeline-depth estimate in S1.3.

**Transmitter, per iteration**

1. Build the spin vector for the current iteration and the weight vector.
2. Prepend 8,192 alternating one and zero pilot symbols.
3. Upsample to 4 samples per symbol.
4. Shape with a root-raised-cosine filter. 51 taps.
5. Apply a pre-emphasis filter for the known transmission-path loss.
6. Resample to the 256 GSa/s AWG grid.

**Receiver, per iteration**

1. Remove the mean.
2. Resample to 2 samples per symbol. A 20-tap anti-aliasing filter runs here.
3. Synchronise against the known pilot sequence and align the spin indices.
4. Train a 51-tap feedforward equaliser on the pilot symbols only.
5. Convolve the trained equaliser with the whole received sequence.
6. Resample to 1 sample per symbol and discard the pilots.
7. Clip the signal to bound amplitude inhomogeneity.
8. Sum the samples belonging to each spin index to finish the multiply-accumulate.

Root-raised-cosine roll-off is per task. See the hyperparameter table below.

Step 7 is called out in the main text as one of the two reasons the solution
quality beats earlier work. The paper does not give the clipping threshold. See
[Open gaps](#open-gaps) item 5.

## 9. Signal encoding

The weight matrix `W` is flattened to a vector `w̃`. The spin vector `x'` is
repeated `N` times to make a vector `x̃` of the same length. The two drive the
two modulators, so the optical output is their element-wise product. Summation
over each index block completes the matrix-vector product (M, eqs. 4, 5, 9).

Detection uses time interleaving. One auxiliary symbol follows every data
symbol, which recovers the sign of a product that a square-law detector would
otherwise lose. This halves the usable AWG memory (S2).

The photodetected voltage before filtering is (M, eq. 6):

```
V_out(t) = I0·R · cos²( (π/4)·(V1/Vπ)·x̃(t) − 1 ) · cos²( (π/4)·(V2/Vπ)·w̃(t) − 1 )
```

With the second modulator in its linear regime this reduces to (M, eq. 8):

```
y(t) = (I0·R/4) · w(t) · sin(x(t))
```

so the spin nonlinearity `σ(·)` the algorithm sees is a sine, bounded to its
half-wave. That is the "half-wave-bounded sinusoidal function" of Methods.

## 10. Algorithm

State update (M, eq. 3):

```
W = [ αI − βJ | −βh ]      x' = [ x | 1 ]
x(t+1) = W · σ(x'(t)) + η(ζ)      η ~ N(0, ζ²)
```

Hyperparameters, all measured at 106 GBaud (EDT 1):

| Task | Spins `N` | Non-zero couplings | `α` | `β` | RRC roll-off |
|---|---|---|---|---|---|
| Square lattice 101×101 | 10,201 | 50,601 | 0.86 | 1 | 0.2 |
| Square lattice 203×203 | 41,209 | 205,233 | 0.86 | 1 | 0.2 |
| Max-cut G22 | 2,000 | 41,980 | 0.7 | 0.78 | 0.4 |
| Max-cut G81 | 20,000 | 100,000 | 0.7 | 0.78 | 0.4 |
| Number partitioning | 256 | 65,536 | 1 | 0.3 | 0.2 |
| Protein folding S30 | 630 | 31,920 | 1 | 0.04 | 0.2 |

Simulation-only values for the large lattices are `α = 0.78`, `β = 0.86`
(S1.7.1). With those, the 143×143 lattice reaches the ground state in 2,653
iterations and the 203×203 lattice in 5,366.

Bifurcation runs use `J = 0` and sweep `α` (Fig. 2a). The worked example uses
`α = 3.5`, 64 GBaud, 262,144 uncoupled spins and 50 iterations.

### Problem mappings

| Problem | Mapping | Source |
|---|---|---|
| Max-cut | `J_ij = k_ij / 2`, `h = 0` | M, eq. 12 |
| Number partitioning | `J_ij = s_i·s_j`, `h = 0` | M, eq. 18 |
| HP protein folding | `J_ij = Q_ij/4`, `h_i = Q_ii/2 + Σ_j Q_ij/4` | M, eqs. 14, 15 |
| Square lattice | `J_ij = −1` on every edge | M |

The QUBO matrix `Q` for protein folding follows the construction of Perdomo-Ortiz
et al. (reference 19 of the paper), which the paper cites and does not restate.

## 11. Published results to reproduce

| Figure | What it shows | Target |
|---|---|---|
| 2a | Bifurcation against `α`, 64 GBaud, `N = 262,144`, `J = 0` | Split above `α₀` |
| 2b | Spin evolution over 50 iterations at `α = 3.5` | Two fixed points |
| 2c | Analogue MVM, 64 GBaud, 500 random 128×128 matrices | 96.2 ± 0.4 % |
| 2d | Bit precision against matrix size, 32 to 200, 64 GBaud | Flat near 4.5 bits |
| 2e | Bit precision against baud rate, 4 to 148 GBaud, size 128 | 5.0342 down to 2.7885 |
| 3a | 101×101 lattice, 106 GBaud | Ground state at about 671 iterations |
| 3b | Solution quality against node count, 1,000 iterations | 97.2 % at 41,209 nodes |
| 3c | Max-cut G22, 800 iterations | 13,289 cuts, 99.48 % |
| 3d | Max-cut G81, 800 iterations | 13,516 cuts, 96.34 % |
| 4c, 4d | HP folding S4 to S30 | 100 % hit rate at S4, > 99 % quality at S30 |
| 5a, 5b | Number partitioning, `N` 16 to 256 | Ground state found, TTS against `N` |
| S7a | Iterations to 90, 93, 95 and 97 % against node count | Near linear in `N` |
| S7b | Iterations to ground state against baud, 400 nodes | 144 at 32 GBaud, 89 at 106, 240 at 128 |
| EDF 2d | Iterations against bulk SOA bias, 400 nodes, 64 GBaud | Minimum at 60 mA |
| EDF 2e | Constant noise against the exponential schedule | 145, 120, 100 iterations |

Accuracy at other baud rates (main text): 98.16 ± 0.31 % at 4 GBaud falling to
90.7 ± 0.94 % at 148 GBaud. Effective precision 4.5 bits at 32 GBaud and
3.3 bits at 106 GBaud.

**Two numbers in the paper disagree with themselves.** The main text calls the
G22 best-known cut 13,352 and Figure 3c calls it 13,359. Only 13,359 is
consistent with the stated 99.48 %, and 13,359 is the accepted literature value,
so this directory uses 13,359. The G81 quality appears as both 96.34 % and
96.36 %. Table S4 and the Figure 3d caption both say 96.34 %.

---

## Open gaps

These are the places where the paper does not say enough, ordered by how much
the answer changes the model.

1. **The time-interleaving scheme, and the rest of the DSP detail.** This is the
   one that currently blocks a result. Methods extracts the wanted product "by
   a combination of DC filtering and a time-interleaving encoding scheme" and
   cites reference 14 rather than stating it. `common/dsp.py` reconstructs a
   scheme that inverts Methods equation 7 exactly, and it is ill conditioned:
   see README, "What does not work yet".

   The paper also does not give the pre-emphasis response, the equaliser
   training length, the clipping threshold, the resampling filter or the
   synchronisation metric. Each has a stated default, marked in the code.

   Two of its stated details do not work as written. The pilot sequence,
   "8,192 alternating ones and zeros", is a single tone: its autocorrelation is
   periodic with a two-symbol period, so it fixes alignment only modulo two
   symbols, and a 51-tap least-squares fit against it is rank deficient.
   Separately, root-raised-cosine shaping is free of intersymbol interference
   only once a matched filter completes it, and this link multiplies the two
   waveforms optically before any matched filter could act.

2. **`Vπ(f)` and `S21(f)` are read off a plot.** Figure S1e is the only source.
   The model fits a smooth curve to the two endpoints and the shape. Any
   measured S-parameter file would replace this.

3. **Receiver amplifier and modulator drivers.** Neither is named. The paper
   gives one driver number, 100 mW, in a projection for a different system
   (S4.2), and one transimpedance amplifier that is also a projection. Both
   models here are placeholders that reproduce the *system* response of Fig. S5a
   rather than any real part.

4. **Bulk SOA noise against bias current.** Only the plot in EDF 2a. The
   annealing schedule is given in scope millivolts squared with no stated
   conversion to spin amplitude. `common/noise.py` calibrates the conversion
   against the three iteration counts of EDF 2e instead.

5. **Clipping threshold.** Named as a main reason for the solution quality and
   never given a value.

6. **Laser relative intensity noise.** Named as a noise source and never given a
   value.

7. **Photodetector responsivity.** Not given. Assumed.

8. **The protein-folding QUBO.** Methods defers the construction of `Q` to
   reference 19 and describes its variables as position-encoded. The paper
   reports 630 spins for a 30-residue sequence, and no position encoding of 30
   residues over a useful 2D lattice gives 630. `common/problems.hp_qubo` builds
   a working position-encoded QUBO that is not the paper's.

9. **Scaling of the couplings.** Nothing in the paper says how `J` is written to
   the converter's full scale, and the published hyperparameters only mean
   something once it is. Number partitioning shows this most sharply: with raw
   integers from {0..16} the couplings reach 256, the `alpha*I` term is two
   orders of magnitude too small to be heard, and every spin takes the same
   sign on iteration one. `problems.number_partition` normalises by the largest
   value. See its docstring.

The published data set at <https://github.com/Shastri-Lab/tfln-ising-nature-paper-2025>
may close several of these. It has not been read for this work.
