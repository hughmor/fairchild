# cmim — the cascaded-modulator photonic Ising machine

A testbench for Al-Kayed et al., *Programmable 200 GOPS Hopfield-inspired
photonic Ising machine*, Nature **648**, 576 (2025),
doi:10.1038/s41586-025-09838-7, and its Supplementary Information.

The machine is an optoelectronic oscillator. Two cascaded thin-film lithium
niobate Mach-Zehnder modulators multiply a time-multiplexed spin vector by a
flattened weight matrix in the optical domain. A quantum dot semiconductor
optical amplifier restores the level, a photodetector sums, and a digital
signal processing stack closes the loop. A second, bulk amplifier injects
controlled noise, which is how the machine anneals.

Two things are modelled here, and keeping them apart is the point.

| | what it is | what it is for |
|---|---|---|
| `models/` + `link/` | the physical layer, as a fairchild deck | frequency response, transfer functions, noise densities, power budget |
| `common/channel.py` | the same physics in numpy, one pass | the 800-iteration solver loop |

One iteration of the paper's largest problem is 434,176 symbols at 106 GBaud,
which is about a million timesteps of circuit simulation. A thousand of them
will not run, in this simulator or any other. So the fast model exists, and
`channel.validate()` pushes one waveform through both and reports the
difference rather than asserting the two agree.

**Read [`SPECS.md`](SPECS.md) first.** It carries every number the models use,
with the citation for each, and it marks what the paper does not say.

## Layout

| path | what |
|---|---|
| [`SPECS.md`](SPECS.md) | every parameter, its source, and the open gaps |
| `models/*.va` | five Verilog-A devices, bundle-port dialect |
| `link/netlists/cmim_link.sp` | the feedforward chain as one deck |
| `link/check.py` | 16 checks of that deck against closed form |
| `common/link.py` | drives the deck: `.op`, `.ac`, `.tran` |
| `common/instruments.py` | the AWG, the oscilloscope, effective resolution, latency |
| `common/dsp.py` | the transmit and receive stacks |
| `common/channel.py` | the fast link model, and its validation |
| `common/problems.py` | lattice, max-cut, number partitioning, HP folding |
| `common/ising.py` | the state update, three ways to evaluate it |
| `replicate.py` | works through the paper's figures |
| `results/` | JSON and PNG per figure |

## Running it

```bash
cargo build --release --bin fairchild          # the deck driver needs it
.venv/bin/python experiments/cmim/link/check.py            # the physical layer
.venv/bin/python experiments/cmim/common/channel.py        # fast model vs deck
MPLBACKEND=Agg .venv/bin/python experiments/cmim/replicate.py
```

Every module under `common/` runs as a script and checks itself:

```bash
cd experiments/cmim/common
.venv/bin/python instruments.py    # derived numbers against the paper's
.venv/bin/python dsp.py            # the interleaving algebra
.venv/bin/python problems.py       # the mappings against hand-computed energies
.venv/bin/python ising.py          # the bifurcation threshold, a 20x20 lattice
```

There is no build step for the models. Four of the five declare
`optical_bundle` ports, whose width the deck decides, so fairchild expands and
compiles the source itself. The deck names them with `.va`, not `.osdi`.

## What works

**The physical layer.** `link/check.py` runs 16 checks, all against an absolute
anchor rather than against another part of the deck, and each names the sabotage
that makes it fail. The power budget matches arithmetic on the paper's loss
figures to five significant figures. The modulator transfer is the analytic
bounded sine to 7e-8. The extinction ratio comes out at 42.00 dB against the
paper's 42. The small-signal response matches an analytic five-pole cascade to
4e-5 dB. Transient noise tracks its analytic discrete-time variance to 3 %.

**The instruments.** Every derived number lands on the paper's: 434,176 symbols
at 106 GBaud, 2.4151 samples per symbol, 9.76 s to fill both AWG channels,
4.11 us of feedforward latency.

**The fast model.** `channel.validate()` refines the deck's timestep and watches
it converge: at 16 samples per AWG sample the two agree on amplitude to 0.5 %
with an 8 % residual. At one sample per AWG sample they disagree by 41 %, and
that is the integrator, not the model. See "Traps" below.

**Bifurcation**, Figure 2(a) and 2(b). A clean pitchfork. The paper reports the
critical feedback strength as a measurement and gives no closed form; there is
one. With `J = 0` the update is `x <- alpha*sin(pi*x/2)`, so the fixed point at
zero loses stability at `alpha_0 = 2/pi = 0.6366`. Measured: 0.6254.

**Precision against solution quality**, Supplementary S3.2 and Figure S14. This
is the paper's headline comparison and it reproduces without forcing. On a 20x20
lattice, over eight seeds:

| per-symbol precision | undithered | with analogue dither |
|---|---|---|
| 2.5 bits | 82.2 % | 92.6 % |
| 3.3 bits | 87.5 % | **98.3 %** |
| 4.5 bits | 94.1 % | 99.3 % |
| 6.0 bits | **98.5 %** | 99.3 % |

Read the two bold entries against each other. That is S3.2's claim — "the
simulation required 6-bit precision to match the solution quality of our
analogue hardware, which achieved similar performance at 3.3-bit effective
resolution" — arriving on its own, from a dither term that is there because the
hardware has one.

**Square lattice scaling**, Figure 3. The right shape, at the paper's own
`alpha = 0.86`, `beta = 1`:

| nodes | our quality | ground state hit |
|---|---|---|
| 100 | 100.00 % | 100 % of runs, iteration 70 |
| 196 | 100.00 % | 100 %, iteration 93 |
| 400 | 98.95 % | 80 %, iteration 344 |
| 900 | 98.51 % | 60 %, iteration 574 |
| 1,600 | 96.97 % | 0 % |
| 2,500 | 96.67 % | 0 % |

The paper reaches the ground state for every lattice up to 10,201 nodes and
97.2 % at 41,209. We fall off around 1,600. The curve has the right shape and
sits about an order of magnitude to the left.

## What does not work yet

**Analogue matrix-vector multiplication, Figures 2(c) to 2(e).** The end-to-end
path through the real transmit DSP, the channel and the real receive DSP
recovers the wanted product with a correlation of about 0.9 at low baud rates
and worse above 64 GBaud, against the paper's 98.16 % accuracy at 4 GBaud.

The cause is isolated and it is not tuning. Turning the noise and the
quantisation off at 4 GBaud, where the link is flat to within 0.01 dB, makes it
*worse*, not better. The fault is the time-interleaving scheme, and that scheme
is a reconstruction: Methods says the wanted product is extracted "by a
combination of DC filtering and a time-interleaving encoding scheme" and cites
the group's earlier cascaded-modulator work rather than stating it.

The scheme in `dsp.interleave` sends `(x, w)` then `(-x, -w)` and takes the sum.
It is algebraically exact — `dsp.py`'s self-check inverts Methods equation 7 to
3e-16 — and it is the only scheme of that shape that works: differencing, or
inverting one channel alone, each leave a single-channel term behind. But it is
ill conditioned. It recovers a product term by cancelling two single-channel
terms that are two to three times larger and that alternate sign at the symbol
rate, which is the Nyquist frequency, which is exactly where a feedforward
equaliser is least accurate. A few per cent of residual intersymbol
interference on the larger terms swamps the smaller one.

Two further findings sit underneath it:

* Root-raised-cosine shaping is free of intersymbol interference only once a
  matched filter completes it into a raised cosine, and this link cannot apply
  one: the two waveforms are multiplied optically before anything matched could
  act. Measured at 8 GBaud, where bandwidth costs nothing, the raw product
  correlates 0.89 with the wanted value at roll-off 0.2 and 0.98 at roll-off
  1.0. Moving the equaliser to after the deinterleave lifts the low-baud
  correlation from about 0 to 0.77 and does not survive to 106 GBaud.
* The paper's stated pilot sequence, "8,192 alternating ones and zeros", cannot
  do either job the paper gives it. Its autocorrelation is periodic with a
  two-symbol period, so it fixes the alignment only modulo two symbols, and a
  51-tap least-squares fit against a single tone is rank deficient. `dsp.pilots`
  defaults to a maximal-length sequence and keeps the paper's version behind
  `kind="alternating"`.

Resolving this needs the interleaving scheme from reference 14, or the published
data set at <https://github.com/Shastri-Lab/tfln-ising-nature-paper-2025>.

**Number partitioning beyond about 32 numbers.** `N = 16` reaches the ground
state on every run. Above that the result is strongly seed dependent — at
`N = 256` the best of five seeds reaches a difference of 3 against a total of
2,131, and the median reaches 1,173. Annealing noise, an exponential schedule
and small initial amplitudes were each tried and none fixed it. The paper
reports ground-state solutions for every size to 256, in 572 iterations.

**Max-cut on G22 and G81.** `problems.read_gset` is written and untested: the
Gset files are not in the tree. Fetch them and it should run.

**HP lattice protein folding.** `problems.hp_qubo` builds a working
position-encoded QUBO, and it is not the paper's encoding. The paper reports 630
spins for a 30-residue sequence, which no position encoding of 30 residues over
a useful 2D lattice produces, and it defers the construction to its reference
19. Any Figure 4 comparison is qualitative until that is pinned down.

## Traps found along the way

Four of these cost real time and none of them is specific to this experiment.

1. **A Verilog-A `white_noise()` written against a solved quantity is silently
   frozen.** With `white_noise(2*q*i_ph)` against the real photocurrent, the
   output noise stayed at 3.2046e-6 V for 0, 0.01, 0.1, 1.0 and 2.0 mW of
   optical power, against 1.05e-3 V predicted at 1 mW. `.noise` returns
   bit-identical densities at 0 mW and 1 mW. A constant argument is exact: swept
   over four decades of `i_n_in` it tracks the analytic discrete-time prediction
   to 0.6 %. `va_receiver` therefore takes its noise operating point as
   parameters, and `Link.noise_params()` measures them from a first pass.
2. **A `white_noise()` inside a potential contribution contributes nothing.**
   fairchild realises noise as random currents, and an optical field wire only
   takes potential contributions, so 20 dBm of spontaneous emission on the field
   moved the receiver's variance by less than one part in ten thousand. The
   simulator has no field-domain optical noise at all — its native budget is
   `4kT/R + 2qI + RIN*I^2`, every term electrical, at the detector.
3. **`OWL(out) <+ OWL(in)` is a no-op**, and the example models under
   `examples/verilog_a/models/` still do it. Wavelength is resolved before the
   solve, so a fixed-port model that writes to a lambda wire compiles, runs and
   propagates nothing — which showed up here as a singular matrix the moment a
   native coupler needed a wavelength label. The `optical_bundle` dialect
   declares the routing and is the fix; it also made the models shorter.
4. **`examples/verilog_a/build.sh` exports `DYLD_LIBRARY_PATH=llvm@18`**, and
   the macOS `openvaf-r` carries its own LLVM 22 at
   `@executable_path/../lib/libLLVM.dylib`. The export shadows it and every
   model dies with "failed to parse bitcode", including ones that compiled the
   day before. Unset it.

And two that are physics rather than tooling:

5. **A transient at one step per AWG sample does not resolve this link's poles.**
   The fastest is at 127 GHz, whose time constant is 1.25 ps, against a 3.9 ps
   step. Backward Euler damps hard: the deck swings 35 % less than it should,
   and comparing a model against it at that step would "prove" the model wrong
   by exactly that much. Sixteen samples per AWG sample is the honest setting.
6. **The second modulator contributes no pole to the spin path.** With its drive
   at DC it is a constant attenuator, and its bandwidth applies to its own
   electrode, not to light passing through. Including it moves the small-signal
   fit from 0.03 dB to 7.6 dB.

## Where the paper disagrees with itself

Each of these is recorded at the point in the code that has to choose.

* **The lattice coupling sign.** Methods assigns "antiferromagnetic coupling
  (J_ij = -1)" and expects a checkerboard ground state. Under Methods
  equation 1, `H = sum J_ij s_i s_j` with no leading minus, `J = -1` is
  minimised by aligned neighbours: a ferromagnet. The paper's own max-cut
  mapping, `J_ij = k_ij / 2`, settles which convention is meant, so
  `problems.square_lattice` returns `+1`.
* **The G22 best-known cut** is 13,352 in the main text and 13,359 in the
  Figure 3c caption. Only 13,359 is consistent with the 99.48 % reported for
  13,289 cuts, and it is the accepted literature value.
* **The G81 solution quality** is 96.36 % in the main text and 96.34 % in both
  Table S4 and the Figure 3d caption.
* **The ASE spectral density** is printed as 1.513 fW/Hz in S4.2. The paper's
  own equation 20, with the paper's own noise figure and gain, gives 1.513e-17
  W/Hz. The mantissa is right and the exponent is out by 100.
* **The effective number of bits** is given in S1.6 as "6 dB" at DC and "5.2 dB"
  at 100 GHz. ENOB is measured in bits.
* **The modulator bandwidth.** Table S1 and the S1.1 text say the electro-optic
  3 dB bandwidth exceeds 110 GHz. The right-hand axis of Figure S1(e) reads
  about 4.5 V of half-wave voltage at 100 GHz against 1.5 V at DC, a factor of
  three, which needs a 3 dB point near 35 GHz. The models take 110 GHz and the
  system response of Figure S5(a) is the calibration target that matters.

## Next

1. Get the time-interleaving scheme from reference 14, or from the published
   data set, and redo Figures 2(c) to 2(e).
2. Fetch the Gset files and run G22 and G81.
3. Get the protein-folding QUBO construction from reference 19.
4. Fit `p_ase_ref_dBm` in `va_ase_source` against Extended Data Figure 2(e)'s
   three iteration counts, which is the only route to the annealing schedule's
   absolute scale.
5. Find why the lattice falls off an order of magnitude earlier in problem size
   than the paper's does.
