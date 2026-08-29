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

**The measurements run in the simulator.** `replicate_deck.py` drives
`link/netlists/cmim_link.sp` through `Link.tran` for every number it reports:
the analogue dot product against baud rate, and the bifurcation. A single-shot
measurement is a few thousand symbols and takes seconds, so there is no excuse
for doing it any other way.

`common/channel.py` holds a numpy model of the same physics, and it is for the
iteration loop and nothing else. One iteration of the paper's largest problem is
434,176 symbols at 106 GBaud, about a million timesteps; a thousand of them will
not run, here or anywhere. `channel.validate()` pushes one waveform through both
and reports the difference rather than asserting they agree, and
`DeckChannel` wears the same interface so anything written against one runs
against the other.

| | what it is | what it is for |
|---|---|---|
| `models/` + `link/` | the physical layer, as a fairchild deck | everything measurable |
| `channel.DeckChannel` | that deck, as a channel object | the reported measurements |
| `channel.Channel` | the same physics in numpy | the 800-iteration solver loop only |

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
| `replicate_deck.py` | the dot product and the bifurcation, **in the deck** |
| `replicate.py` | the solver benchmarks, on the numpy channel |
| `results/` | JSON and PNG per figure |

## Running it

```bash
cargo build --release --bin fairchild          # the deck driver needs it
.venv/bin/python experiments/cmim/link/check.py            # the physical layer
.venv/bin/python experiments/cmim/common/channel.py        # fast model vs deck

# the measurements, every sample out of a Newton solve
MPLBACKEND=Agg .venv/bin/python experiments/cmim/replicate_deck.py

# the solver benchmarks, which use the numpy channel because they cannot not
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

**Bifurcation**, Figure 2(a) and 2(b). The paper reports the critical feedback
strength as a measurement and gives no closed form; there is one. With `J = 0`
the update is `x <- alpha*sin(pi*x/2)`, so the fixed point at zero loses
stability at `alpha_0 = 2/pi = 0.6366`. `ising.bifurcation` measures 0.6254 —
but read the next paragraph before quoting it.

**That numpy pitchfork carries no noise of any kind.** `ising.bifurcation`
iterates the map directly: no channel, no shot noise, no quantisation, no
transient. It is a check on the nonlinearity and nothing more, which is why it
is suspiciously clean. `replicate_deck.bifurcation` is the real one: every
iteration is a transient through the modulators, the amplifier and the receiver
with `.options trannoise=1`, so the spins carry shot noise, amplifier noise and
signal-spontaneous beat noise.

It shows what the noiseless map cannot: spins clustered near zero below `2/pi`,
a broad partially-separated region through the transition, and full pinning at
+/-1 above about 1.6, with occasional noise-driven escapes at low feedback. The
"fully pinned" threshold measures 1.23 against the noiseless onset at 0.637 —
noise and finite precision push complete bifurcation above the analytic
instability, which is the behaviour Figure 2(a)'s heatmap shows.

Three things had to be right before it showed a threshold at all, and all three
are recorded in that function:

* **The loop gain is a hardware constant, calibrated once.** Refitting it each
  iteration by least squares against a noisy measurement is regression
  dilution: the gain comes out systematically small, the loop sags below
  threshold, and nothing separates at any alpha.
* **A row of `alpha*I` has one non-zero**, so pruning the zeros gives a
  single-element dot product with no block for the interleaving to cancel over.
  Sending each row as `block` symbols of `alpha/block` against the same spin
  accumulates the product coherently while the noise adds as its square root.
  The two encodings compute the same row; only one can be read back.
* **"Bifurcated" has to mean pinned, not displaced.** With real noise in the
  loop a sub-threshold spin wanders, and a loose threshold counts that as a
  bifurcation: at `|x| > 0.25` every feedback strength including zero came out
  bifurcated while the scatter plainly showed a transition at `2/pi`.

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

**Analogue dot product, Figures 2(c) to 2(e).** Runs in the deck. Accuracy falls
monotonically with baud rate, which is the behaviour the paper reports, and sits
below the paper's absolute numbers:

| baud | ours, spin channel nonlinear | ours, both linear | paper |
|---|---|---|---|
| 4 GBaud | 92.8 % (3.79 bits) | 89.7 % | 98.16 % (5.03 bits) |
| 32 GBaud | 82.1 % (2.48 bits) | 67.5 % | — (4.5 bits) |
| 64 GBaud | 67.2 % (1.61 bits) | 50.7 % | 96.2 % |
| 106 GBaud | 16.1 % | 5.8 % | — (3.3 bits) |
| 148 GBaud | 1.1 % | 1.7 % | 90.7 % (2.79 bits) |

Four things were wrong before this worked, and all four were mine:

1. **The extraction is the accumulation.** The interleaving sends `(x, w)` then
   `(-x, -w)`; summing a whole block cancels both single-channel terms and
   doubles the product. Averaging each pair first and accumulating afterwards is
   algebraically identical and numerically hopeless — it cancels terms several
   times larger than the wanted one, symbol by symbol, at the Nyquist frequency.
   See `dsp.accumulate_interleaved`.
2. **The receiver integrates.** Methods equation 9 offers the summation "in the
   analogue domain by means of an integrator before sampling" and Supplementary
   S4.2 builds the scaled system around "the integrating photoreceiver". Reading
   two samples per symbol and equalising instead gives 23.7 % where integration
   gives 99.1 %, on an ideal detector with no channel at all. The equaliser is
   actively harmful: trained on zero-mean pilots against a stream with a large
   DC, its DC gain came out at 0.20, and a block sum is a low-frequency
   quantity.
3. **The transmit pulse must be flat at the symbol centre.** The two channels
   are multiplied optically, so each channel's own intersymbol interference
   becomes a cross term that is second order in the data and cancels nowhere.
   Root-raised cosine gives 21 % with no channel impairment at all. The paper
   says NRZ where it compares itself with other machines (S3.3); `dsp.nrz` is
   the default.
4. **The integration boundaries are the converter's.** At 106 GBaud a symbol is
   2.4151 samples, so a DAC holds symbols for two or three samples and the
   integrator has to use the same edges. Truncating them to integers made every
   non-integer samples-per-symbol rate collapse while 4, 16, 32, 64 and
   128 GBaud looked fine. It reads exactly like a bandwidth limit.

What still separates us from the paper is smaller and named: the transmit
**pre-emphasis is written and not applied**. Figure S5(a) is measured with it, at
3 dB down at 55 GHz against 33 GHz for the raw link, so it is worth several bits
at the top of the range. `dsp.preemphasis` builds the filter from `Link.ac`.

Two notes on the paper's stated DSP, both recorded rather than worked around:

* The pilot sequence, "8,192 alternating ones and zeros", is a single tone. Its
  autocorrelation is periodic with a two-symbol period, so it fixes alignment
  only modulo two symbols, and a 51-tap least-squares fit against it is rank
  deficient. `dsp.pilots` defaults to a maximal-length sequence and keeps the
  paper's version behind `kind="alternating"`.
* The 51-tap feedforward equaliser is kept, behind `rx(mode="equalise")`, for
  symbol-level readout. It is not on the path that computes a dot product.

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
