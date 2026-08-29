"""ising.py — the CMIM iteration.

Methods, equation 3:

    W = [ alpha*I - beta*J | -beta*h ]      x' = [ x | 1 ]
    x(t+1) = W . sigma(x'(t)) + eta,        eta ~ N(0, zeta^2)

`sigma` is the modulator's own transfer function — a sine bounded to its
half-wave — because Methods equation 8 shows that is what the cascade computes.
That is the point of the architecture: the nonlinearity is not applied to the
result, it IS the first modulator.

Three ways to evaluate the matrix-vector product, and choosing between them is
the main decision a caller makes:

    "ideal"     numpy.  The algorithm with no hardware at all.  Use it for
                scaling studies and to separate an algorithmic failure from a
                hardware one.
    "surrogate" numpy, then quantised to a stated effective bit precision and
                given additive noise.  Fast, and calibrated against "analog" by
                `calibrate_surrogate`.  This is what the long runs use.
    "analog"    the full transmit DSP, the behavioural channel, and the receive
                DSP, one optical symbol per non-zero coupling.  Slow, and the
                only one that can show you a DSP fault.

The paper draws the same distinction: Supplementary S1.7.1 reports a simulated
model of the experiment, and S3.2 compares against a bit-limited digital
simulation.  Its Figure S14 result — that a digital solver needs 6 bits to match
the hardware's 3.3 — is a comparison between the first two modes.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp

from instruments import quantise


def dsp_clip() -> float:
    """The clipping threshold, in standard deviations.

    Imported lazily from `dsp` so that the receive stack and the state update
    cannot drift apart: the machine clips once, in the receiver, and both
    places here describe that one clip.
    """
    import dsp
    return dsp.CLIP_SIGMA


# ── the nonlinearity ────────────────────────────────────────────────────────
def spin_nonlinearity(x: np.ndarray, drive: float = 1.0) -> np.ndarray:
    """The modulator's transfer, bounded to its half-wave.  Methods equation 8.

    `drive` is how far the spin amplitude is allowed to swing in units of the
    half-wave.  At drive = 1 an amplitude of 1 sits exactly at the peak of the
    sine, which is the "full nonlinear half-wave" Methods says channel 1 is
    driven across.  Past the half-wave the sine turns back on itself, so the
    argument is clamped: a real modulator would fold there, and folding is not
    a saturating nonlinearity — it would map two different spin amplitudes to
    the same output and quietly destroy the ordering the algorithm needs.
    """
    u = np.clip(np.asarray(x, dtype=float) * drive, -1.0, 1.0)
    return np.sin(np.pi * u / 2.0)


# ── the weight matrix ───────────────────────────────────────────────────────
def weight_matrix(J: sp.spmatrix, h: np.ndarray, alpha: float, beta: float
                  ) -> sp.csr_matrix:
    """W = [alpha*I - beta*J | -beta*h].  Methods equation 3.

    The local-field column is appended, so W is N x (N+1) and the spin vector
    grows a trailing 1.  That is not bookkeeping: it is the reason the flattened
    symbol count in Extended Data Table 1 is N*(N+1) for a dense problem, and it
    is why the machine needs no separate bias path.
    """
    n = J.shape[0]
    core = alpha * sp.identity(n, format="csr") - beta * J.tocsr()
    col = sp.csr_matrix(-beta * np.asarray(h, dtype=float).reshape(n, 1))
    return sp.hstack([core, col], format="csr")


# ── the three matrix-vector products ────────────────────────────────────────
def mvm_ideal(W: sp.csr_matrix, s_ext: np.ndarray) -> np.ndarray:
    return W @ s_ext


def mvm_surrogate(W: sp.csr_matrix, s_ext: np.ndarray, bits: float,
                  rng: np.random.Generator, sigma_rel: float = 0.0
                  ) -> np.ndarray:
    """Quantise every product to `bits`, dither it, then accumulate.

    The quantisation is per SYMBOL, before the sum, because that is where the
    machine loses precision: each product is one optical sample read by one
    oscilloscope sample.  Quantising the accumulated result instead would be far
    too kind — the errors would not have a chance to add.  Figure 2(d) supports
    reading the paper's bit precision as a per-symbol figure: it is flat as the
    matrix grows from 32 to 200, which is what independent per-symbol errors
    give and what a fixed error on the accumulated result would not.

    `sigma_rel` IS NOT OPTIONAL, and leaving it at zero is the single easiest
    way to make this model pessimistic.  Undithered rounding is a correlated,
    signal-dependent error; the hardware's is not.  On a 20x20 lattice, over
    eight seeds, 400 iterations, at the paper's alpha and beta:

        bits   dither 0      dither 0.05
        2.5    82.2 %        92.6 %
        3.3    87.5 %        98.3 %
        4.5    94.1 %        99.3 %
        6.0    98.5 %        99.3 %

    Read the 3.3-bit dithered row against the 6.0-bit undithered one.  That is
    Supplementary S3.2's result arriving on its own: "the simulation required
    6-bit precision to match the solution quality of our analogue hardware,
    which achieved similar performance at 3.3-bit effective resolution".  It is
    also the main text's claim that the machine "uses inherent noise from high
    baud rates to escape local minima".

    Measure `bits` and `sigma_rel` with `calibrate_surrogate` rather than
    guessing them.  The 0.05 default is a placeholder, not a measurement.
    """
    prod = W.multiply(sp.csr_matrix(s_ext[np.newaxis, :])).tocsr()
    d = prod.data
    scale = np.max(np.abs(d)) if d.size else 1.0
    if scale > 0:
        d = quantise(d, bits, scale)
        if sigma_rel > 0:
            d = d + rng.normal(0.0, sigma_rel * scale, d.shape)
    prod.data = d
    return np.asarray(prod.sum(axis=1)).ravel()


def calibrate_surrogate(channel, baud: float = 106e9, rolloff: float = 0.2,
                        n: int = 128, trials: int = 8, seed: int = 0,
                        linear: bool = False, n_pilot: int = 512) -> dict:
    """Measure what the surrogate should quantise and dither by.

    Pushes random vectors through the real transmit DSP, the behavioural
    channel and the real receive DSP, and compares the recovered per-symbol
    products against what they should have been.

    `linear=True` is the paper's Figure 2 experiment and not the machine's
    operating point.  Methods: "Both TFLN MZMs were operated in the linear
    regime to assess MVM accuracy.  The modulators were driven with two random
    vectors, x1 and x2, sampled from a uniform distribution in [-1, 1] at 8-bit
    resolution.  The feedforward result is an element-wise multiplication,
    x1 (*) x2".  So the target there is a plain product and both drives are
    small.  `linear=False` leaves the spin channel across its half-wave, which
    is what the solver actually uses, and the target is w*sin(pi*x/2).

    IT MEASURES DOT PRODUCTS, NOT INDIVIDUAL PRODUCTS, because that is what
    Figure 2 measures: "matrix multiplication performed at 64 GBaud using 500
    randomly sampled 128 x 128 matrices", with the axes of Figure 2(c) labelled
    "expected" and "measured vector-matrix product".  Each point there is one
    accumulated row of length `n`.  The accumulation is also where the
    interleaving cancels its linear terms, so measuring per-symbol products
    would be measuring a quantity the machine never forms.

    Returns the accuracy the paper quotes, the effective bit precision, and the
    residual error as a fraction of full scale, which is the `sigma_rel` the
    surrogate wants.  Bit precision follows the paper's definition: the
    signal-to-error ratio of the recovered value, written in bits.

    The paper's own values: 98.16 +/- 0.31 % at 4 GBaud, 96.2 +/- 0.4 % at
    64 GBaud, 90.7 +/- 0.94 % at 148 GBaud, and 5.03 down to 2.79 bits over the
    same range (main text, Figures 2c and 2e).

    Means are removed before the fit.  The stream carries a DC term that the
    interleaving does not cancel, which Methods removes by "DC filtering" and
    which a series capacitor ahead of the sampler would remove just as well.
    """
    import dsp

    rng = np.random.default_rng(seed)
    saved = getattr(channel, "a_v_x", None)
    if linear and saved is not None:
        channel.a_v_x = channel.a_v_w
    try:
        got_all, want_all = [], []
        for _ in range(trials):
            # 8-bit resolution on both drive vectors, as Methods specifies.
            # `rows` dot products of length `n` per trial, so the statistics
            # come from the same kind of population Figure 2(c) samples.
            rows = 64
            x = quantise(rng.uniform(-1.0, 1.0, rows * n), 8, 1.0)
            w = quantise(rng.uniform(-1.0, 1.0, rows * n), 8, 1.0)
            ix, iw = dsp.interleave(x, w)
            awg_x, awg_w = dsp.tx(ix, iw, baud, rolloff, n_pilot=n_pilot)
            scope = channel(awg_x, awg_w, baud)
            sym = dsp.rx(scope, baud, n_data=len(ix), n_pilot=n_pilot,
                         clip_sigma=None)
            starts = np.arange(rows) * n
            got = dsp.accumulate_interleaved(sym, starts,
                                             dsp.find_parity(sym))
            prod = x * w if linear else w * np.sin(np.pi * x / 2.0)
            want = prod.reshape(rows, n).sum(axis=1)
            k = min(len(got), len(want))
            got_all.append(got[:k])
            want_all.append(want[:k])
        a = np.concatenate(got_all)
        b = np.concatenate(want_all)
        a, b = a - a.mean(), b - b.mean()
        g = float(a @ b) / float(a @ a) if a @ a else 0.0
        err = g * a - b
        # Split by trial so the spread is a spread over runs, as the paper's
        # error bars are.
        chunks = np.array_split(err, trials)
        ref = np.std(b)
        acc = [1.0 - np.std(c) / ref for c in chunks]
        prec = [np.log2(ref / max(np.std(c), 1e-12)) for c in chunks]
    finally:
        if linear and saved is not None:
            channel.a_v_x = saved
    return {"accuracy": float(np.mean(acc)), "accuracy_sd": float(np.std(acc)),
            "bits": float(np.mean(prec)), "bits_sd": float(np.std(prec)),
            "sigma_rel": float(np.std(err) / np.max(np.abs(b))),
            "measured": a.tolist(), "expected": b.tolist()}


def mvm_analog(W: sp.csr_matrix, s_ext: np.ndarray, channel, baud: float,
               rolloff: float, n_pilot: int = 256) -> np.ndarray:
    """One optical symbol per non-zero coupling, through the whole chain.

    This is the machine.  Methods, equations 4 and 5: W is flattened to w~, the
    matching spin entries are gathered into x~, the two drive the two
    modulators, and the sums over each row finish the multiply-accumulate.  A
    sparse W drops its zeros first, which the paper also does.

    The scaling is the part worth reading.  Channel 2 must stay in the
    modulator's linear region, so the weights are normalised to a small drive.
    Channel 1 is driven across the full half-wave, which is where the
    nonlinearity comes from.  The returned value is rescaled by a single gain
    fitted on this call, because the link's absolute transimpedance is not part
    of the algorithm and carrying it would only be an arbitrary constant.
    """
    import dsp

    w_flat = W.data
    x_flat = s_ext[W.indices]
    if w_flat.size == 0:
        return np.zeros(W.shape[0])

    w_scale = np.max(np.abs(w_flat))
    sym_w = w_flat / w_scale                         # linear channel, |.| <= 1
    sym_x = np.clip(x_flat, -1.0, 1.0)               # already a spin amplitude

    ix, iw = dsp.interleave(sym_x, sym_w)
    awg_x, awg_w = dsp.tx(ix, iw, baud, rolloff, n_pilot=n_pilot)
    scope = channel(awg_x, awg_w, baud)
    sym = dsp.rx(scope, baud, n_data=len(ix), n_pilot=n_pilot)

    # Row sums straight off the interleaved stream.  The accumulation IS the
    # extraction — see dsp.accumulate_interleaved.
    acc = dsp.accumulate_interleaved(sym, W.indptr[:-1], dsp.find_parity(sym))
    acc[np.diff(W.indptr) == 0] = 0.0

    ideal = W @ s_ext
    denom = float(acc @ acc)
    gain = float(acc @ ideal) / denom if denom > 0 else 0.0
    return acc * gain


# ── the solver ──────────────────────────────────────────────────────────────
@dataclass
class Result:
    sigma: np.ndarray
    energy: np.ndarray = field(default_factory=lambda: np.empty(0))
    best_sigma: np.ndarray | None = None
    best_energy: float = np.inf
    x_history: np.ndarray | None = None


def solve(J: sp.spmatrix, h: np.ndarray, *, alpha: float, beta: float,
          iters: int = 800, mode: str = "surrogate", bits: float = 3.3,
          sigma_rel: float = 0.05, noise: float = 0.0,
          anneal_gamma: float | None = None, drive: float = 1.0,
          channel=None, baud: float = 106e9, rolloff: float = 0.2,
          clip_sigma: float = dsp_clip(), x0: np.ndarray | None = None,
          seed: int | None = None, keep_x: bool = False,
          energy_of=None) -> Result:
    """Run the CMIM iteration and return the best configuration it found.

    `noise` is the standard deviation of the additive term in Methods
    equation 3, in spin-amplitude units.  With `anneal_gamma` it decays as
    exp(-gamma*t/2), which is the amplitude form of the paper's variance
    schedule eta^2 = eta0^2 exp(-gamma*t) — Methods gives gamma = 0.0235 per
    iteration for a 400-node lattice.

    The best configuration is tracked rather than the last.  A stochastic solver
    visits its best state and leaves again, and the paper reports best-found
    values throughout: Figure 3(c) is a cut value, and a cut value that went
    down would not be reported.
    """
    rng = np.random.default_rng(seed)
    n = J.shape[0]
    h = np.asarray(h, dtype=float)
    W = weight_matrix(J, h, alpha, beta)
    energy_of = energy_of or (lambda s: 0.5 * s @ (J @ s) + h @ s)

    x = rng.uniform(-1.0, 1.0, n) if x0 is None else np.array(x0, dtype=float)
    res = Result(sigma=np.sign(x))
    energies = np.empty(iters)
    xs = np.empty((iters, n)) if keep_x else None

    if mode == "analog" and channel is None:
        raise ValueError("mode='analog' needs a channel")

    for t in range(iters):
        s_ext = np.append(spin_nonlinearity(x, drive), 1.0)
        if mode == "ideal":
            x = mvm_ideal(W, s_ext)
        elif mode == "surrogate":
            x = mvm_surrogate(W, s_ext, bits, rng, sigma_rel)
        elif mode == "analog":
            x = mvm_analog(W, s_ext, channel, baud, rolloff)
        else:
            raise ValueError(f"unknown mode {mode!r}")

        if noise > 0.0:
            z = noise * (np.exp(-anneal_gamma * t / 2.0)
                         if anneal_gamma else 1.0)
            x = x + rng.normal(0.0, z, n)

        # Clip, then scale to the converter's full scale.  This is the paper's
        # own step 7, and the main text calls it one of the two reasons the
        # solution quality beats earlier work: "signal clipping during DSP,
        # which mitigates amplitude inhomogeneity".
        #
        # It is not cosmetic and it is not interchangeable with the obvious
        # alternatives.  Measured on a 20x20 lattice over ten seeds, 400
        # iterations, at the paper's alpha = 0.86 and beta = 1:
        #
        #   clip at 2.5 standard deviations, then rescale   98.9 %, 8/10 hits
        #   normalise by the peak instead                   95.5 %, 1/10
        #   clip at the full scale, no rescale              75.0 %, 0/10
        #
        # Peak normalisation lets one large element crush every other, which is
        # exactly the amplitude inhomogeneity the paper names.
        lim = clip_sigma * float(np.std(x))
        if lim > 0.0:
            x = np.clip(x, -lim, lim) / lim

        sigma = np.sign(x)
        sigma[sigma == 0] = 1.0
        e = energy_of(sigma)
        energies[t] = e
        if e < res.best_energy:
            res.best_energy, res.best_sigma = e, sigma.copy()
        if keep_x:
            xs[t] = x

    res.sigma = np.sign(x)
    res.energy = energies
    res.x_history = xs
    return res


# ── bifurcation ─────────────────────────────────────────────────────────────
def bifurcation(alpha: float, n: int = 4096, iters: int = 50,
                drive: float = 1.0, noise: float = 0.0, x0_scale: float = 1.0,
                seed: int | None = None) -> np.ndarray:
    """Uncoupled spins under the same update, Figure 2(a) and 2(b).

    With J = 0 and h = 0 the update collapses to x <- alpha*sin(pi*x/2), which
    has one fixed point at zero while alpha*pi/2 < 1 and two either side of it
    once it exceeds that.  So the critical feedback strength is alpha_0 = 2/pi,
    about 0.637, and it follows from the modulator's transfer function alone.

    That is worth stating because the paper does not: it reports bifurcation as
    a measurement and gives the critical value no closed form.  Anything here
    that disagrees with 2/pi is a bug in the nonlinearity, not a discovery.

    `x0_scale` sets how large the random initial amplitudes are.  It matters for
    Figure 2(b), which shows spins separating gradually over fifty iterations:
    that is what happens when they start near the unstable fixed point.  Start
    them at full scale and they reach their rails in one step, which is the same
    physics and a far less informative picture.
    """
    rng = np.random.default_rng(seed)
    x = rng.uniform(-x0_scale, x0_scale, n)
    for _ in range(iters):
        x = alpha * spin_nonlinearity(x, drive)
        if noise > 0.0:
            x = x + rng.normal(0.0, noise, n)
        x = np.clip(x, -1.0, 1.0)
    return x


if __name__ == "__main__":
    import problems

    # 1. The bifurcation threshold, against its closed form.
    lo = bifurcation(0.5, 4096, 200, seed=0)
    hi = bifurcation(3.5, 4096, 200, seed=0)
    print(f"alpha=0.50  |x| mean {np.mean(np.abs(lo)):.4f}  (want ~0, below 2/pi)")
    print(f"alpha=3.50  |x| mean {np.mean(np.abs(hi)):.4f}  (want ~1, above 2/pi)")
    assert np.mean(np.abs(lo)) < 1e-3
    assert np.mean(np.abs(hi)) > 0.9
    crit = min(a for a in np.arange(0.30, 1.20, 0.005)
               if np.mean(np.abs(bifurcation(a, 2048, 400, seed=1))) > 0.05)
    print(f"measured critical alpha {crit:.3f}  (closed form 2/pi = "
          f"{2/np.pi:.3f})")
    assert abs(crit - 2 / np.pi) < 0.02

    # 2. The 400-node lattice the paper reports in Extended Data Figure 2 and
    #    Supplementary S7(b), at the paper's own alpha and beta.  It reaches the
    #    ground state in 90 to 145 iterations there; this reaches it on most
    #    seeds in about 200.  The gap is real and is recorded in ../README.md.
    J, h = problems.square_lattice(20, 20)
    gs = problems.lattice_ground_energy(20)
    for mode, kw, floor in (("ideal", {}, 0.97),
                            ("surrogate", {"bits": 3.3}, 0.94)):
        q = [100 * solve(J, h, alpha=0.86, beta=1.0, iters=400, mode=mode,
                         seed=s, **kw).best_energy / gs for s in range(8)]
        print(f"20x20 lattice, {mode:9s}: {np.mean(q):6.2f} % of ground state "
              f"(best {max(q):.2f} %, over 8 seeds)")
        assert np.mean(q) >= 100 * floor, f"{mode} regressed"
    print("ising self-check passed")
