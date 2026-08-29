"""instruments.py — the test and measurement set-up, as the paper describes it.

The CMIM's feedback path is a Keysight M8199B arbitrary waveform generator and a
Keysight UXR real-time oscilloscope, both at 256 GSa/s, with an Intel i7-8700K
doing the signal processing offline between them.  Al-Kayed et al., Nature 648,
576 (2025), Methods "Experimental test bed", and Supplementary S1.3 and S1.6.

The instruments matter for three reasons and only three:

1. **Memory bounds the problem.**  The AWG holds 2^20 samples per channel, and
   the time-interleaving scheme spends one auxiliary symbol per data symbol, so
   the largest Ising problem the machine can hold is set here and nowhere else.
2. **Effective resolution bounds the answer.**  Not the nominal 8 and 10 bits,
   but the 6 bits at DC falling to 5.2 bits at 100 GHz that Supplementary S1.6
   quotes for the whole system.  That is the dominant precision limit in the
   loop, and Figure 2(e) is a measurement of it.
3. **Upload and capture bound the wall clock.**  They do not bound the machine:
   the paper is careful to separate the 4.1 us feedforward latency it claims
   from the ~10 s of instrument overhead it measures.

What is NOT here: the analogue bandwidth of either instrument.  The transmit
path's bandwidth lives in `models/va_mod_driver.va` and the receive path's in
`models/va_receiver.va`, so that one place decides what the link's frequency
response is.  Filtering here as well would count it twice.

Supplementary S1.6 prints the two ENOB figures with the unit "dB".  Effective
number of bits is measured in bits.  They are read as bits.
"""
from __future__ import annotations

import numpy as np

# ── Keysight M8199B ─────────────────────────────────────────────────────────
AWG_FS = 256e9            # samples per second, Methods
AWG_MEMORY = 2 ** 20      # samples per channel, S1.3
AWG_BITS = 8              # nominal, data sheet — ASSUMED, not in the paper
AWG_CHANNELS = 2          # channel 1 spins, channel 2 weights

# ── Keysight UXR ────────────────────────────────────────────────────────────
RTO_FS = 256e9            # samples per second, Methods
RTO_BITS = 10             # nominal, data sheet — ASSUMED, not in the paper
RTO_CAPTURE_S = 0.81      # for 2^22 samples, S1.3

# ── system effective resolution, S1.6 ───────────────────────────────────────
ENOB_DC = 6.0             # bits
ENOB_100GHZ = 5.2         # bits

# ── the paper's own operating point ─────────────────────────────────────────
BAUD_NOMINAL = 106e9      # every benchmark in Extended Data Table 1
SPS_NOMINAL = AWG_FS / BAUD_NOMINAL       # 2.4151, matching S1.3's "2.415"
PILOT_SYMBOLS = 8192      # alternating ones and zeros, Methods


def samples_per_symbol(baud: float) -> float:
    """How many AWG samples one symbol occupies.  Not an integer in general."""
    return AWG_FS / baud


def max_symbols(baud: float, interleaved: bool = True) -> int:
    """The largest symbol count the AWG can hold at this baud rate.

    S1.3 quotes 434,176 symbols at 106 GBaud, which is 2^20 / 2.4151.  That
    figure is the raw count; `interleaved` halves it, because Methods sends one
    auxiliary symbol after every data symbol so that a square-law detector can
    still recover the sign of a product (S2).
    """
    n = int(AWG_MEMORY / samples_per_symbol(baud))
    return n // 2 if interleaved else n


def max_spins(baud: float, dense: bool = True) -> int:
    """The largest Ising problem this memory holds.

    A dense N-spin problem flattens an N x (N+1) weight matrix, so it needs
    N*(N+1) symbols.  A sparse one needs one symbol per non-zero element.
    """
    m = max_symbols(baud)
    return int((np.sqrt(1 + 4 * m) - 1) / 2) if dense else m


def enob(f_hz: float | np.ndarray) -> float | np.ndarray:
    """System effective number of bits against frequency, S1.6.

    Two points and a straight line between them.  The paper gives no third
    point, so nothing here claims a shape.  Held flat above 100 GHz rather than
    extrapolated to zero, because a linear fit run past its data is a fiction
    that keeps its units.
    """
    return ENOB_DC + (ENOB_100GHZ - ENOB_DC) * np.clip(
        np.asarray(f_hz, dtype=float) / 100e9, 0.0, 1.0)


def enob_for_baud(baud: float, rolloff: float = 0.2) -> float:
    """The effective resolution a signal at this baud rate sees.

    A root-raised-cosine signal occupies (1 + rolloff) * baud / 2 of bandwidth.
    Take the resolution at the mid-band of that, which is where most of the
    energy is.  S1.5.1 uses the same "half the baud rate" reasoning for the
    maximum frequency component.
    """
    return float(enob((1.0 + rolloff) * baud / 4.0))


def quantise(x: np.ndarray, bits: float, full_scale: float | None = None
             ) -> np.ndarray:
    """Uniform mid-tread quantisation to a possibly fractional number of bits.

    Fractional bits are the point.  An effective number of bits is not a
    converter setting, it is a signal-to-noise ratio written in bits, and the
    honest way to impose 3.3 of them is a step size of full_scale / 2^3.3.

    Values outside the full scale are clipped, not wrapped.  A converter clips.
    """
    fs = float(np.max(np.abs(x))) if full_scale is None else float(full_scale)
    if fs <= 0.0:
        return np.zeros_like(x)
    step = 2.0 * fs / (2.0 ** bits)
    return np.clip(np.round(x / step) * step, -fs, fs)


def awg_upload_time(n_samples: int, channels: int = AWG_CHANNELS) -> float:
    """Seconds to push a waveform from the host to the AWG.  S1.3, equation 3.

    Fitted by the authors to their own measurements.  It is the largest single
    term in the experimental iteration time, about 9.76 s for a full memory on
    two channels, and it is an artefact of the instrument rather than of the
    architecture.
    """
    return channels * (4.4981e-6 * n_samples + 0.1628)


def dsp_time(n_symbols: int) -> float:
    """Seconds of host signal processing per iteration.

    S1.3 measures 24.47 ms for 65,536 symbols at 106 GBaud and 58.28 ms for
    606,308 symbols at 148 GBaud.  Those two points are not proportional, so
    this scales the first one linearly and is therefore an estimate: the paper
    itself scales the same way when it computes time-to-solution.
    """
    return 24.47e-3 * n_symbols / 65536


def feedforward_latency(n_symbols: int, baud: float, path_m: float = 1.0,
                        n_eff: float = 1.5, t_conv: float = 5e-9) -> float:
    """The latency the paper actually claims, S1.3 equation 1.

    Modulating every symbol, plus the time of flight, plus one sample of data
    conversion.  At 106 GBaud with 434,176 symbols this gives 4.11 us, which is
    the number in the abstract.  It excludes the AWG and oscilloscope entirely,
    and the paper says so.
    """
    return n_symbols / baud + n_eff * path_m / 299792458.0 + t_conv


def time_to_solution(n_iter: int, n_symbols: int, baud: float,
                     p_success: float = 1.0, pipelined: bool = True) -> float:
    """Time to reach the ground state with 99 % probability.  S1.3, equation 6.

    With `pipelined`, the iteration is `max(feedforward, DSP)`, which is what
    the paper's starred figures assume.  Without it, the two add.
    """
    ff = feedforward_latency(n_symbols, baud)
    ds = dsp_time(n_symbols)
    per_iter = max(ff, ds) if pipelined else ff + ds
    t_ann = n_iter * per_iter
    if p_success >= 1.0:
        return t_ann
    return t_ann * np.log(1 - 0.99) / np.log(1 - p_success)


if __name__ == "__main__":
    print(f"samples per symbol at 106 GBaud : {samples_per_symbol(106e9):.4f}"
          f"   (S1.3 says 2.415)")
    print(f"max symbols at 106 GBaud        : {max_symbols(106e9, False):,}"
          f"   (S1.3 says 434,176)")
    print(f"max dense spins at 106 GBaud    : {max_spins(106e9):,}"
          f"   (the paper demonstrates 256)")
    print(f"max sparse spins at 106 GBaud   : {max_spins(106e9, False):,}"
          f"   (the paper demonstrates 205,233 couplings)")
    print(f"ENOB at 53 GHz                  : {enob(53e9):.3f} bits")
    print(f"AWG upload, full memory, 2 ch   : {awg_upload_time(2**20):.2f} s"
          f"   (S1.3 says 9.76 s)")
    print(f"feedforward latency, 106 GBaud  : "
          f"{feedforward_latency(434176, 106e9)*1e6:.2f} us"
          f"   (S1.3 says 4.11 us)")
