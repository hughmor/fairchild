"""dsp.py — the transmitter and receiver signal processing of the CMIM.

Al-Kayed et al., Nature 648, 576 (2025) puts the digital signal processing at
the centre of the result: "embedding DSP ... within optical computation enhances
convergence and solution quality".  Methods, "Experimental test bed", gives the
stage list.  Supplementary S1.3 gives two tap counts, in passing, inside a
latency estimate.  Between them that is everything the paper says, so read the
open-gap notes below before trusting any default here.

Transmitter, per iteration (Methods)
    1. build the spin vector and the flattened weight vector
    2. prepend 8,192 alternating one and zero pilot symbols
    3. upsample to 4 samples per symbol
    4. root-raised-cosine pulse shaping
    5. pre-emphasis for the known transmission-path loss
    6. resample to the 256 GSa/s AWG grid

Receiver, per iteration (Methods)
    1. remove the mean
    2. resample to 2 samples per symbol
    3. synchronise against the pilot sequence
    4. train a 51-tap feedforward equaliser on the pilots ALONE
    5. convolve it with the whole received sequence
    6. resample to 1 sample per symbol and discard the pilots
    7. clip
    8. sum the samples belonging to each spin

WHAT THE PAPER DOES NOT SAY, and what this module assumes instead:

  * The pre-emphasis response.  Here it is the regularised inverse of the link's
    own measured small-signal response, capped at `PREEMPH_MAX_DB`.  That is
    what "compensate for known transmission-path losses" has to mean, and the
    measurement comes from the deck rather than from a guess.
  * The equaliser training length and the synchronisation metric.  Least
    squares on the pilots, and a cross-correlation peak.  Both are the obvious
    choice and neither is stated.
  * The clipping threshold, which the main text calls one of the two reasons
    the solution quality beats earlier work and never gives a value.
    `CLIP_SIGMA` is a free parameter.
  * The anti-aliasing filter, mentioned only as "one 20-tap" filter in a
    pipeline-depth estimate.  `resample_poly` designs its own.
  * The time-interleaving scheme.  See `interleave` below: the paper cites its
    earlier work for it and states only the consequence.
"""
from __future__ import annotations

from fractions import Fraction

import numpy as np
from scipy import signal

from instruments import AWG_FS, PILOT_SYMBOLS

SPS_SHAPE = 4          # samples per symbol for pulse shaping, Methods
RRC_TAPS = 51          # S1.3's pipeline estimate says 51 for pulse shaping
FFE_TAPS = 51          # and 51 for the feedforward equaliser
SPS_RX = 2             # the receiver resamples to two, Methods
PREEMPH_MAX_DB = 12.0  # ASSUMED: a real driver runs out of swing eventually
CLIP_SIGMA = 2.5       # ASSUMED: see the note above


# ── pulse shaping ───────────────────────────────────────────────────────────
def rrc(beta: float, sps: int = SPS_SHAPE, taps: int = RRC_TAPS) -> np.ndarray:
    """Root-raised-cosine impulse response, normalised to unit energy.

    `beta` is the roll-off.  Extended Data Table 1 uses 0.2 for the lattice,
    number-partitioning and folding tasks and 0.4 for both max-cut graphs.

    The two removable singularities are handled by substitution rather than by
    a small-denominator guard, because a guard leaves a step in the impulse
    response whose size depends on the tolerance you picked.
    """
    if not 0.0 <= beta <= 1.0:
        raise ValueError("roll-off must be in [0, 1]")
    n = np.arange(taps) - (taps - 1) / 2.0
    t = n / sps
    h = np.empty_like(t)

    zero = np.isclose(t, 0.0)
    h[zero] = 1.0 + beta * (4.0 / np.pi - 1.0)

    if beta > 0.0:
        sing = np.isclose(np.abs(t), 1.0 / (4.0 * beta))
        h[sing] = (beta / np.sqrt(2.0)) * (
            (1.0 + 2.0 / np.pi) * np.sin(np.pi / (4.0 * beta))
            + (1.0 - 2.0 / np.pi) * np.cos(np.pi / (4.0 * beta)))
    else:
        sing = np.zeros_like(zero)

    rest = ~(zero | sing)
    tr = t[rest]
    h[rest] = (np.sin(np.pi * tr * (1.0 - beta))
               + 4.0 * beta * tr * np.cos(np.pi * tr * (1.0 + beta))) \
        / (np.pi * tr * (1.0 - (4.0 * beta * tr) ** 2))
    return h / np.sqrt(np.sum(h ** 2))


def _rational(ratio: float, max_den: int = 4096) -> tuple[int, int]:
    f = Fraction(ratio).limit_denominator(max_den)
    return f.numerator, f.denominator


def to_awg_grid(x: np.ndarray, baud: float, sps_in: int = SPS_SHAPE
                ) -> np.ndarray:
    """Resample a pulse-shaped waveform onto the 256 GSa/s AWG grid.

    At 106 GBaud the target is 2.4151 samples per symbol, which is not an
    integer and never will be: the instrument has one clock and the baud rate is
    chosen freely against it.  A rational approximation to 1/4096 is exact to
    better than a part in ten thousand of a symbol.
    """
    up, down = _rational((AWG_FS / baud) / sps_in)
    return signal.resample_poly(x, up, down)


def from_awg_grid(x: np.ndarray, baud: float, sps_out: int = SPS_RX
                  ) -> np.ndarray:
    """The inverse: oscilloscope samples down to `sps_out` per symbol."""
    up, down = _rational(sps_out / (AWG_FS / baud))
    return signal.resample_poly(x, up, down)


# ── the pilot sequence ──────────────────────────────────────────────────────
def pilots(n: int = PILOT_SYMBOLS, kind: str = "prbs") -> np.ndarray:
    """The training sequence, as +/-1 symbols.

    Methods says "8,192 alternating ones and zeros in all experiments", and
    `kind="alternating"` is exactly that.  It is not the default, and the reason
    is worth stating rather than hiding.

    An alternating sequence is a single tone at the Nyquist frequency.  It
    cannot do either job the paper gives the pilots:

      * SYNCHRONISATION.  Its autocorrelation is periodic with a period of two
        symbols, so a correlation peak fixes the alignment only modulo two
        symbols.  Getting the spin indices wrong by one symbol swaps every spin
        with its neighbour, which the machine cannot detect and the answer does
        not survive.
      * EQUALISER TRAINING.  A 51-tap least-squares fit against a one-tone input
        is rank deficient by construction: the design matrix has two distinct
        rows repeated 4,000 times.  The equaliser it produces corrects one
        frequency and does what it likes at the rest.

    So the default is a maximal-length shift-register sequence, which is flat in
    frequency, has a single sharp autocorrelation peak, and trains a full-rank
    equaliser.  Any real high-speed link uses one.  The paper's phrase probably
    describes a preamble whose job is coarse framing, with the training done
    elsewhere; it does not say so, so this is recorded as an open gap rather
    than asserted.  ../SPECS.md open gap 1.
    """
    if kind == "alternating":
        return np.where(np.arange(n) % 2 == 0, 1.0, -1.0)
    if kind != "prbs":
        raise ValueError("pilot kind is 'prbs' or 'alternating'")
    # PRBS-15, x^15 + x^14 + 1.  Long enough that 8,192 symbols never repeat.
    reg, out = 0x7F5B, np.empty(n)
    for i in range(n):
        bit = ((reg >> 14) ^ (reg >> 13)) & 1
        reg = ((reg << 1) | bit) & 0x7FFF
        out[i] = 1.0 if bit else -1.0
    return out


# ── time interleaving ───────────────────────────────────────────────────────
def interleave(x: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """One auxiliary symbol after every data symbol, both channels inverted.

    RECONSTRUCTED, NOT QUOTED.  Methods says the wanted product is extracted "by
    a combination of DC filtering and a time-interleaving encoding scheme" and
    cites the group's earlier cascaded-modulator work for the scheme itself.
    Supplementary S2 says the scheme costs one auxiliary symbol per data symbol
    and exists "to correctly compute the results involving negative coupling
    matrix and spin values".  This is the scheme that satisfies both statements.

    The detected voltage is Methods equation 7, up to a constant:

        v = 1 + sin(x) + w + w*sin(x)

    Send (x, w), then (-x, -w):

        v_even = 1 + sin(x) + w + w*sin(x)
        v_odd  = 1 - sin(x) - w + w*sin(x)
        (v_even + v_odd)/2 - 1 = w*sin(x)

    The sum cancels both single-channel terms exactly, for any x and any w, and
    the leftover constant is the DC that Methods removes.  Neither channel alone
    would do: differencing leaves one of them behind.

    `deinterleave` undoes it.  Test them as a pair before trusting either.
    """
    out_x = np.empty(2 * len(x))
    out_w = np.empty(2 * len(w))
    out_x[0::2], out_x[1::2] = x, -x
    out_w[0::2], out_w[1::2] = w, -w
    return out_x, out_w


def accumulate_interleaved(y: np.ndarray, starts: np.ndarray,
                           parity: int = 0) -> np.ndarray:
    """Sum each interleaved block, which is where the linear terms cancel.

    THIS IS THE EXTRACTION, and pairing symbols off two at a time is not.

    After DC removal the received symbol stream is, up to constants,

        y_j = a*x_j + b*w_j + c*x_j*w_j

    and `interleave` has arranged that every element appears twice, once as
    `(x_i, w_i)` and once as `(-x_i, -w_i)`.  Summing a whole block therefore
    gives

        sum_j y_j = a*sum(x_i - x_i) + b*sum(w_i - w_i) + 2c*sum(x_i*w_i)

    so both single-channel terms cancel across the block and the wanted product
    comes out at twice its amplitude.  The multiply-accumulate the algorithm
    needs is the same operation that does the cancelling: there is no separate
    extraction step and nothing to condition.

    Intersymbol interference cancels with them, which is the part that matters.
    The electrical path is linear, so its contribution to the linear terms is a
    convolution, and the sum of a convolution over an antisymmetric sequence is
    itself antisymmetric — zero to the edges of the block.

    Averaging each pair first, then accumulating, is algebraically identical and
    numerically far worse: it cancels two terms several times larger than the
    wanted one, symbol by symbol, at the Nyquist frequency where an equaliser is
    least accurate, and the ISI no longer has a whole block to cancel over.
    Measured on this model at 106 GBaud, pairwise recovers the product with an
    accuracy of about 26 %; block accumulation recovers the dot product to a few
    per cent.

    `starts` are the block boundaries in DATA-element units, as
    `scipy.sparse.csr_matrix.indptr` gives them.  They are doubled here because
    the transmitted stream carries two symbols per element.
    """
    y = y[parity:]
    idx = 2 * np.asarray(starts, dtype=int)
    idx = np.clip(idx, 0, max(len(y) - 1, 0))
    if len(y) == 0:
        return np.zeros(len(starts))
    return 0.5 * np.add.reduceat(y, idx)


def deinterleave(y: np.ndarray, parity: int = 0) -> np.ndarray:
    """Average each data symbol with its auxiliary partner.

    Kept because it is the clearest way to SHOW that the interleaving inverts
    Methods equation 7, and the module self-check uses it for exactly that.  Do
    not use it to extract a matrix-vector product: see
    `accumulate_interleaved`, which is both better conditioned and the operation
    the algorithm wanted anyway.

    `parity` says which of the two symbols starts a pair.  Getting it wrong does
    not degrade the answer, it destroys it.  Use `find_parity`.
    """
    y = y[parity:]
    n = (len(y) // 2) * 2
    return 0.5 * (y[:n:2] + y[1:n:2])


def find_parity(y: np.ndarray) -> int:
    """Which offset pairs each data symbol with its own auxiliary.

    The receiver cannot know this from the synchronisation alone.  A correlation
    peak fixes the sample the pilots start on, and one symbol is two samples at
    two samples per symbol, so an error of a single symbol flips the pairing
    while leaving the pilot fit untouched.  It has to be resolved from the data.

    The discriminator comes from the scheme itself.  Under the right pairing the
    DIFFERENCE of a pair carries the two single-channel terms, which are several
    times larger than the product, and the SUM carries only the product.  Under
    the wrong pairing neither cancels and the two come out the same size.  So
    take the offset that maximises the ratio of the two variances.

    Measured on this model, over 4 to 148 GBaud: the right parity gives a ratio
    of 19 to 190, the wrong one gives 0.9 to 1.1.  There is no ambiguous case in
    the operational configuration.  With both modulators driven small-signal the
    margin narrows to about 2, which is a symptom of the same conditioning
    problem described in ../README.md.
    """
    best, arg = -1.0, 0
    for p in (0, 1):
        s = y[p:]
        n = (len(s) // 2) * 2
        if n < 4:
            continue
        d = 0.5 * (s[:n:2] - s[1:n:2])
        m = 0.5 * (s[:n:2] + s[1:n:2])
        r = float(np.var(d) / max(np.var(m), 1e-30))
        if r > best:
            best, arg = r, p
    return arg


# ── pre-emphasis ────────────────────────────────────────────────────────────
def preemphasis(f: np.ndarray, h: np.ndarray, fs: float = AWG_FS,
                taps: int = RRC_TAPS,
                max_db: float = PREEMPH_MAX_DB) -> np.ndarray:
    """An FIR that inverts the link's measured response, within a boost limit.

    `f` and `h` come from `Link.ac`.  The inverse is regularised by the boost
    cap rather than by a Tikhonov term, because a cap is what a real driver
    imposes: past its swing it stops boosting and starts clipping.

    `fs` is the grid the filter will run on, which for a zero-order-hold
    transmitter is the converter's own sample rate.  The result is windowed to
    `taps` and normalised to unit gain at DC, so it changes the shape of the
    spectrum and not the size of the signal.

    Methods puts a pre-emphasis filter after the pulse shaping, "to compensate
    for known transmission-path losses", and Figure S5(a) is measured with it
    applied: the loop is 3 dB down at 55 GHz there against 33 GHz for the raw
    link.  Leaving it out is the single largest reason a modelled accuracy sits
    below the paper's at high baud rates.
    """
    grid = np.fft.rfftfreq(4096, d=1.0 / fs)
    mag = np.interp(grid, f, np.abs(h), left=np.abs(h[0]), right=np.abs(h[-1]))
    mag = np.maximum(mag / mag[0], 1e-12)
    boost = np.minimum(1.0 / mag, 10.0 ** (max_db / 20.0))

    imp = np.fft.irfft(boost, n=4096)
    imp = np.roll(imp, taps // 2)[:taps] * np.hamming(taps)
    return imp / np.sum(imp)


# ── transmitter ─────────────────────────────────────────────────────────────
def nrz(sym: np.ndarray, baud: float) -> np.ndarray:
    """Zero-order hold onto the AWG grid: one held level per symbol.

    THE TRANSMIT PULSE MUST BE FREE OF INTERSYMBOL INTERFERENCE AT THE SAMPLING
    INSTANT, and that is a hard requirement of this architecture rather than a
    preference.  The two channels are multiplied OPTICALLY, before anything on
    the receive side can act.  Write the shaped waveforms as x + e and w + d,
    where e and d are each channel's own intersymbol interference:

        (x + e)(w + d) = x*w + x*d + w*e + e*d

    The two cross terms are products of the data with the interference.  They
    are second order in the data, so the interleaving does not cancel them —
    it cancels terms that are LINEAR in the transmitted sequence — and neither
    does the block accumulation.  Receive-side filtering is harmless for exactly
    the same reason: it acts after the product, it is linear in the product
    stream, and a linear filter of unit DC gain conserves a block sum.

    Measured at 4 GBaud with no channel impairment at all, recovering
    accumulated dot products: root-raised-cosine shaping gives 21 %, and this
    gives what the ideal symbols give.

    The paper says NRZ where it compares itself with other machines
    (Supplementary S3.3, "operating a non-return-to-zero (NRZ) signal for spin
    encoding at > 100 GBaud"), and root-raised cosine where it lists the DSP
    stages (Methods).  Both are true of a real transmitter — a DAC holds a level
    for a sample and the pre-emphasis shapes what follows — but the pulse the
    MODULATOR sees has to be flat at the symbol centre, so the hold is what
    matters here.

    Samples per symbol is not an integer in general: at 106 GBaud it is 2.4151,
    because the instrument has one clock and the baud rate is chosen against it.
    So the symbol index is computed per sample rather than by repeating.
    """
    sps = AWG_FS / baud
    n_out = int(np.floor(len(sym) * sps))
    idx = np.minimum((np.arange(n_out) / sps).astype(int), len(sym) - 1)
    return np.asarray(sym, dtype=float)[idx]


def tx(sym_x: np.ndarray, sym_w: np.ndarray, baud: float, rolloff: float,
       preemph: np.ndarray | None = None, n_pilot: int = PILOT_SYMBOLS,
       shape: str = "nrz") -> tuple[np.ndarray, np.ndarray]:
    """Both AWG channels, from spin symbols and weight symbols.

    `sym_x` and `sym_w` must already be interleaved and the same length.  The
    pilot sequence is prepended to both, as Methods describes, so the receiver's
    synchronisation and equaliser training see the same channel the data does.

    `shape` is `"nrz"` or `"rrc"`.  Read `nrz` before choosing `"rrc"`: root
    raised cosine is what Methods lists, and it costs most of the answer,
    because its intersymbol interference reaches the optical multiplication.
    """
    if len(sym_x) != len(sym_w):
        raise ValueError("the two channels must carry the same symbol count")
    p = pilots(n_pilot)
    h = rrc(rolloff)
    # The preamble goes on the SPIN channel, and the weight channel is held at
    # full scale beside it.  Putting the same pilots on both channels makes the
    # preamble useless: the detected value is w*sin(pi*x/2) with w = x = p, which
    # for a +/-1 alphabet is a constant, and an equaliser trained against a
    # constant learns nothing.  Held instead, the second modulator is a fixed
    # attenuator through the preamble, so what arrives is proportional to
    # sin(pi*p/2) — the pilot sequence itself, which is what the receiver's
    # training target assumes.
    out = []
    for s in (np.concatenate([p, sym_x]),
              np.concatenate([np.ones(n_pilot), sym_w])):
        if shape == "nrz":
            grid = nrz(s, baud)
            if preemph is not None:
                grid = np.convolve(grid, preemph, mode="same")
        elif shape == "rrc":
            up = np.zeros(len(s) * SPS_SHAPE)
            up[::SPS_SHAPE] = s
            shaped = np.convolve(up, h, mode="same")
            if preemph is not None:
                shaped = np.convolve(shaped, preemph, mode="same")
            grid = to_awg_grid(shaped, baud)
        else:
            raise ValueError("shape is 'nrz' or 'rrc'")
        # Scale to the converter's full scale, which is what a transmitter does
        # and what makes the drive levels downstream mean anything.  Pulse
        # shaping overshoots a symbol by 10 to 20 %, so without this the spin
        # channel is driven past the modulator's half-wave, where the sine turns
        # back on itself and two different spin amplitudes give one output.
        peak = float(np.max(np.abs(grid)))
        out.append(grid / peak if peak > 0 else grid)
    return out[0], out[1]


# ── receiver ────────────────────────────────────────────────────────────────
def _windows(rx2: np.ndarray, n_sym: int, taps: int, off: int) -> np.ndarray:
    """One row of `taps` samples per symbol, centred on that symbol.

    The equaliser is written as a correlation over an explicit window rather
    than as a convolution.  A convolution has to be handed its taps in the right
    order and read at the right delay, and getting either wrong produces a
    working-looking equaliser that is quietly misaligned.  A window has no
    orientation to get wrong.
    """
    pad = taps + abs(off)
    xp = np.pad(rx2, (pad, pad))
    starts = pad + off + SPS_RX * np.arange(n_sym) - taps // 2
    return np.lib.stride_tricks.sliding_window_view(xp, taps)[starts]


def _sync(rx2: np.ndarray, ref: np.ndarray, search: int = 4096,
          refine: int = 6) -> tuple[int, np.ndarray]:
    """Where the pilots start, in 2-samples-per-symbol units, and the equaliser.

    A correlation peak alone is not enough, and that cost a day.  It locates the
    preamble to within a symbol or two, which is fine for reading the pilots and
    useless for the data: the transmitted stream carries two symbols per matrix
    element, so an offset error of one symbol moves every block boundary and
    flips which auxiliary symbol pairs with which datum.  The accumulated
    products then come out uncorrelated with the answer while the pilot fit
    still looks perfect, which is the worst way for a bug to present.

    So the correlation gives a starting point and the offset is then chosen by
    the only criterion that reflects the data: train the equaliser at each
    candidate and keep the one whose PILOT RESIDUAL is smallest.  That pins the
    symbol grid, not just the preamble, because the equaliser can absorb a
    fractional delay and cannot absorb a whole-symbol one.

    Returns `(offset, coefficients)` so the caller does not train twice.
    """
    up = np.zeros(len(ref) * SPS_RX)
    up[::SPS_RX] = ref
    span = min(len(rx2), len(up) + search)
    c = signal.correlate(rx2[:span], up, mode="full")
    coarse = int(np.argmax(np.abs(c))) - (len(up) - 1)

    best = None
    for d in range(-refine, refine + 1):
        off = coarse + d * SPS_RX
        try:
            coef = _train_ffe(rx2, ref, off)
        except np.linalg.LinAlgError:
            continue
        resid = _windows(rx2, len(ref), FFE_TAPS, off) @ coef - ref
        score = float(np.mean(resid ** 2))
        if best is None or score < best[0]:
            best = (score, off, coef)
    if best is None:
        raise ValueError("could not train an equaliser at any candidate offset")
    return best[1], best[2]


def _train_ffe(rx2: np.ndarray, ref: np.ndarray, off: int,
               taps: int = FFE_TAPS, ridge: float = 1e-6) -> np.ndarray:
    """Least-squares feedforward equaliser, trained on the pilots alone.

    Methods is explicit that only the pilots train it, which matters: an
    equaliser trained on the data would adapt to the spin pattern and quietly
    become part of the algorithm instead of part of the channel.

    A small ridge term keeps the solve finite when the pilot is rank deficient,
    which the paper's own alternating sequence is.  It is not there to improve
    the fit and at 1e-6 it does not.
    """
    A = _windows(rx2, len(ref), taps, off)
    g = A.T @ A + ridge * np.trace(A.T @ A) / taps * np.eye(taps)
    return np.linalg.solve(g, A.T @ ref)


def integrate_symbols(samples: np.ndarray, baud: float, n_sym: int,
                      offset: float = 0.0) -> np.ndarray:
    """Average the samples inside each symbol.  The integrating photoreceiver.

    THIS IS THE RECEIVER, and a feedforward equaliser reading two samples per
    symbol is not.  Methods equation 9: the summation that finishes the
    multiply-accumulate "can be performed digitally after ADC sampling or in the
    analogue domain by means of an integrator before sampling".  Supplementary
    S4.2 builds the scaled architecture around "the integrating photoreceiver",
    whose "bandwidth requirements are relaxed" precisely because it integrates.

    Why it is not interchangeable with sampling and equalising, measured on an
    ideal square-law detector with no channel impairment at all, recovering
    accumulated dot products:

        integrate over each symbol                        99.1 %
        two samples per symbol, then a 51-tap equaliser   23.7 %
        two samples per symbol, no equaliser              51.3 %

    The equaliser is not merely unnecessary here, it is harmful, and the reason
    is worth keeping.  It is trained on zero-mean pilot symbols against a stream
    carrying a large DC term, so least squares learns to suppress low
    frequencies: its DC gain came out at 0.20.  The quantity being recovered is
    a block sum, which is a low-frequency quantity, so the equaliser attenuates
    exactly what the accumulation is trying to measure.  Symbol-level fidelity
    stayed at r = 0.9986 while the block sums fell to r = 0.65 — the residual
    was small but correlated, and summing 256 symbols amplifies a correlated
    error by 256 where it amplifies a white one by 16.

    Integration also rejects the transmit-side intersymbol interference that the
    optical multiplication would otherwise turn into a second-order error: with
    `nrz` shaping the detected waveform is flat across a symbol, and the mean of
    a flat thing is that thing.

    `offset` is the sub-symbol phase of the boundaries, in samples.

    THE BOUNDARIES MUST BE THE CONVERTER'S, NOT THE IDEAL SYMBOL PERIOD'S.  A
    digital-to-analogue converter holds a level for whole samples, so symbol `k`
    occupies samples `ceil(k*sps)` up to `ceil((k+1)*sps)`, and at 106 GBaud,
    where `sps` is 2.4151, consecutive symbols are therefore two or three
    samples wide.  Integrating instead over the exact interval `[k*sps,
    (k+1)*sps)` straddles the staircase and mixes up to 40 % of a neighbouring
    symbol into every result.

    The symptom was diagnostic once seen: every baud rate with an INTEGER number
    of samples per symbol worked — 4, 16, 32, 64 and 128 GBaud — and every one
    without collapsed.  It reads exactly like a bandwidth limit and is not one.
    Integrating on the converter's own boundaries is also what an integrating
    photoreceiver does in hardware, where the staircase is the thing arriving.

    `offset` is the sub-symbol phase of the boundaries, in samples.
    """
    sps = AWG_FS / baud
    x = np.asarray(samples, dtype=float)
    cum = np.concatenate([[0.0], np.cumsum(x)])
    edges = offset + np.arange(n_sym + 1) * sps
    e = np.clip(np.ceil(edges - 1e-9).astype(int), 0, len(x))
    width = np.maximum(e[1:] - e[:-1], 1)
    return (cum[e[1:]] - cum[e[:-1]]) / width


def _sync_integrated(x: np.ndarray, baud: float, ref: np.ndarray,
                     n_total: int) -> tuple[np.ndarray, float]:
    """Find the symbol phase, then integrate on it.

    The boundaries the integrator needs are a sub-sample phase, not a lag in
    whole samples, so this searches the phase directly and scores each candidate
    by how well the integrated preamble matches the known pilots.  Cheap: one
    pass per candidate over a preamble.
    """
    sps = AWG_FS / baud
    best = None
    # A FIXED number of phases across one symbol, never a number derived from
    # the sample count.  Deriving it gave 7 candidates at 106 GBaud, where a
    # symbol is 2.4 samples, and 65 at 4 GBaud, where it is 64 — so the phase
    # was resolved to a sixtieth of a symbol at the easy end and a fifth at the
    # hard one.  That alone took the accuracy at 106 GBaud from usable to noise,
    # and it looked exactly like a bandwidth limit.
    for frac in np.linspace(0.0, sps, 65)[:-1]:
        for lag in range(0, 8):
            sym = integrate_symbols(x, baud, len(ref), frac + lag * sps)
            a = sym - sym.mean()
            b = ref - ref.mean()
            d = float(a @ a)
            r = abs(float(a @ b)) / np.sqrt(d * float(b @ b)) if d > 0 else 0.0
            if best is None or r > best[0]:
                best = (r, frac + lag * sps)
    off = best[1]
    return integrate_symbols(x, baud, n_total, off), off


def rx(samples: np.ndarray, baud: float, n_data: int,
       n_pilot: int = PILOT_SYMBOLS, clip_sigma: float | None = CLIP_SIGMA,
       pilot_kind: str = "prbs", mode: str = "integrate") -> np.ndarray:
    """Recover `n_data` symbols from an oscilloscope capture.

    Returns the equalised, clipped symbol sequence — still interleaved.  Pass it
    to `deinterleave` and then sum over each spin's block to finish the
    multiply-accumulate, which is step 8 of the Methods list.

    The equaliser has real work to do here, and not only because of the
    channel's bandwidth.  Root-raised-cosine shaping is free of intersymbol
    interference only when the matched filter completes it into a raised
    cosine, and this link cannot: the two waveforms are multiplied optically
    before anything matched can be applied, so what arrives carries the ISI of
    one root filter.  Measured on this model at 8 GBaud, where bandwidth costs
    nothing, the raw product correlates 0.89 with the wanted value at roll-off
    0.2 and 0.98 at roll-off 1.0.  The equaliser closes that gap.
    """
    x = np.asarray(samples, dtype=float)
    x = x - np.mean(x)                                    # 1
    ref = pilots(n_pilot, pilot_kind)
    total = n_pilot + n_data
    if mode == "integrate":
        sym, _ = _sync_integrated(x, baud, ref, total)    # 2, 3
    elif mode == "equalise":
        x2 = from_awg_grid(x, baud, SPS_RX)               # 2
        off, coef = _sync(x2, ref)                        # 3, 4
        sym = _windows(x2, total, FFE_TAPS, off) @ coef   # 5, 6
    else:
        raise ValueError("mode is 'integrate' or 'equalise'")
    sym = sym[n_pilot:total]
    if clip_sigma is not None and len(sym):               # 7
        lim = clip_sigma * np.std(sym)
        sym = np.clip(sym, -lim, lim)
    return sym


def accumulate(sym: np.ndarray, block: int) -> np.ndarray:
    """Sum each block of `block` symbols.  Step 8, and the "accumulate" of MAC.

    `block` is N+1 for a dense N-spin problem: N couplings plus the local field
    that the appended 1 in the spin vector multiplies.
    """
    n = len(sym) // block
    return sym[:n * block].reshape(n, block).sum(axis=1)


if __name__ == "__main__":
    # The interleaving algebra, against the transfer function it claims to
    # invert.  If this fails, every multiply-accumulate downstream is wrong.
    rng = np.random.default_rng(0)
    x, w = rng.uniform(-1.2, 1.2, 5000), rng.uniform(-0.35, 0.35, 5000)
    ix, iw = interleave(x, w)
    v = 1.0 + np.sin(ix) + iw + iw * np.sin(ix)     # Methods equation 7
    got = deinterleave(v) - 1.0
    err = np.max(np.abs(got - w * np.sin(x)))
    print(f"interleaving residual        : {err:.3e}   (must be ~1e-16)")
    assert err < 1e-12, "the interleaving scheme does not invert equation 7"

    h = rrc(0.2)
    print(f"RRC energy                   : {np.sum(h**2):.6f}   (must be 1)")
    # Two cascaded root-raised-cosine filters make a raised cosine, which is
    # zero at every symbol instant but its own.  That is the property that makes
    # the shaping free of inter-symbol interference.
    rc = np.convolve(h, h)
    peak = len(rc) // 2
    isi = np.abs(rc[peak - SPS_SHAPE::-SPS_SHAPE][1:]).max()
    print(f"worst intersymbol term       : {isi/rc[peak]:.4f}   (want << 1)")
    assert isi / rc[peak] < 0.05

    from instruments import samples_per_symbol
    a, b = tx(*interleave(np.sign(rng.normal(size=64)),
                          rng.uniform(-0.3, 0.3, 64)), 106e9, 0.2, n_pilot=256)
    want = (256 + 128) * samples_per_symbol(106e9)
    print(f"AWG samples for 384 symbols  : {len(a)}   (want about {want:.0f})")
    assert abs(len(a) - want) < 8
    print("dsp self-check passed")
