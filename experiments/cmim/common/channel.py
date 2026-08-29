"""channel.py — the analogue link as a fast behavioural model.

WHY THIS EXISTS.  One iteration of the paper's largest problem is 434,176
symbols at 106 GBaud, and a benchmark is 800 to 1,000 iterations.  A circuit
simulation of one such iteration is about a million timesteps.  A thousand of
them is not going to happen, in this simulator or any other.

So the machine is modelled twice, and the two are kept honest against each
other:

    link.py     the real deck.  Newton, every device, every pole.  Used to
                CHARACTERISE: frequency response, transfer function, noise
                densities, and to validate this file.
    channel.py  the same physics in numpy, evaluated in one pass.  Used to RUN
                the iteration loop.

Nothing here is invented.  The linear response is measured from the deck, the
noise densities are measured from the deck, and `validate()` pushes the same
waveform through both and reports the difference.  The paper does the same
thing for the same reason and says so, in Supplementary S1.7.1: "we developed a
simulation model that emulates our experimental setup".

WHAT IS ASSERTED RATHER THAN MEASURED.  A small-signal sweep gives the product
of every pole in the chain, and this model needs to know which of them sit
BEFORE the modulator's nonlinearity and which after.  `PRE_POLES` is that split.
The product is measured, only the split is asserted, and `validate()` is what
would catch a wrong one.
"""
from __future__ import annotations

import numpy as np

from instruments import AWG_FS, enob_for_baud, quantise
from link import Link

# Poles ahead of the modulator's sine, in hertz.  These are `va_mod_driver`'s
# `f_3db`, the electrode capacitance against the driver's 50 ohm source in
# parallel with the modulator's 50 ohm termination, and `va_tfln_mzm`'s `f_eo`.
# `link/check.py` pins the whole five-pole cascade against the deck to 4e-5 dB,
# which is where these three came from.
PRE_POLES = (65e9, 1.0 / (2 * np.pi * 25.0 * 50e-15), 110e9)

Q_ELECTRON = 1.602176634e-19
BW_OPT = 4.77e12          # the receiver's optical noise bandwidth, Hz


def _pole_cascade(f: np.ndarray, poles) -> np.ndarray:
    h = np.ones_like(f, dtype=complex)
    for fp in poles:
        h = h / (1.0 + 1j * f / fp)
    return h


def _fft_filter(x: np.ndarray, h_of_f, fs: float) -> np.ndarray:
    n = len(x)
    f = np.fft.rfftfreq(n, d=1.0 / fs)
    return np.fft.irfft(np.fft.rfft(x) * h_of_f(f), n=n)


class Channel:
    """The CMIM feedforward path, from AWG samples to oscilloscope samples.

    Build it with `Channel.from_link()`, which measures what it needs from the
    deck.  Building it by hand is possible and is how you sweep a parameter the
    deck does not carry.
    """

    # Driver gains, referred to the AWG's full scale.  The driver is a 50 ohm
    # source into the modulator's 50 ohm termination, so the modulator sees half
    # the open-circuit swing:
    #
    #   spin channel    a_v = 1.5  ->  0.75 V  ->  pi*v/V_pi = pi/2 exactly,
    #                   the full nonlinear half-wave Methods drives channel 1
    #                   across, reached at full scale and never passed.
    #   weight channel  a_v = 0.4  ->  0.20 V  ->  0.42 rad, where the sine is
    #                   within 3 % of its tangent, which is the linear regime
    #                   Methods keeps channel 2 in.
    #
    # The output swings are set well above both, so the soft clip limits an
    # overdrive rather than costing gain at the nominal drive.
    def __init__(self, f: np.ndarray, h_total: np.ndarray, *,
                 v_pi: float = 1.5, a_v_x: float = 1.5, a_v_w: float = 0.4,
                 sw_x: float = 4.0, sw_w: float = 4.0, clip_n: float = 4.0,
                 p_in_w: float = 4.47e-3, il_dB: float = 6.0, er_dB: float = 42.0,
                 g_soa: float = 44.3, split: float = 0.5,
                 responsivity: float = 0.8, z_t: float = 300.0,
                 rx_swing: float = 2.0,
                 i_ph_n: float = 0.0, p_sig_n: float = 0.0,
                 rho_ase_n: float = 0.0, i_n_in: float = 3e-12,
                 seed: int | None = None):
        self.f, self.h_total = np.asarray(f), np.asarray(h_total)
        self.h_total = self.h_total / np.abs(self.h_total[0])
        for k, v in dict(v_pi=v_pi, a_v_x=a_v_x, a_v_w=a_v_w, sw_x=sw_x,
                         sw_w=sw_w, clip_n=clip_n, p_in_w=p_in_w, il_dB=il_dB,
                         er_dB=er_dB, g_soa=g_soa, split=split,
                         responsivity=responsivity, z_t=z_t, rx_swing=rx_swing,
                         i_ph_n=i_ph_n, p_sig_n=p_sig_n, rho_ase_n=rho_ase_n,
                         i_n_in=i_n_in).items():
            setattr(self, k, v)
        self.rng = np.random.default_rng(seed)

    # ── construction ─────────────────────────────────────────────────────
    @classmethod
    def from_link(cls, lk: Link | None = None, *, seed: int | None = None,
                  n_per_decade: int = 30, **params) -> "Channel":
        """Measure everything measurable from the deck, then build the model.

        `params` are `.param` overrides, so this is also how you set the bulk
        SOA bias: `Channel.from_link(i_noise_mA=60, p_ase_ref=-38)`.
        """
        lk = lk or Link()
        f, h = lk.ac(1e8, 3e11, n_per_decade, **params)
        npar = lk.noise_params(**params)
        op = lk.op(["v(comb_re_0)", "v(comb_im_0)",
                    "v(soao_re_0)", "v(soao_im_0)"], **params)
        lo = {k.lower(): v for k, v in op.items()}
        p_in = lo["v(comb_re_0)"] ** 2 + lo["v(comb_im_0)"] ** 2
        p_out = lo["v(soao_re_0)"] ** 2 + lo["v(soao_im_0)"] ** 2
        return cls(f, h, g_soa=p_out / p_in, seed=seed, **npar)

    # ── the physical stages ──────────────────────────────────────────────
    def _soft_clip(self, v, swing):
        u = v / swing
        return swing * u / (1.0 + np.abs(u) ** self.clip_n) ** (1.0 / self.clip_n)

    def _h_pre(self, f):
        return _pole_cascade(f, PRE_POLES)

    def _h_post(self, f):
        """Whatever is left of the measured response after the pre-filter.

        Interpolated in log magnitude and in unwrapped phase, because linear
        interpolation of a complex number across a decade of a rolling-off
        response undershoots the magnitude and wanders in phase.
        """
        mag = np.interp(f, self.f, np.log(np.abs(self.h_total)),
                        left=0.0, right=np.log(np.abs(self.h_total[-1])))
        ph = np.interp(f, self.f, np.unwrap(np.angle(self.h_total)),
                       left=0.0, right=np.unwrap(np.angle(self.h_total))[-1])
        return np.exp(mag + 1j * ph) / self._h_pre(f)

    def __call__(self, awg_x: np.ndarray, awg_w: np.ndarray,
                 baud: float, noise: bool = True,
                 quantise_rx: bool = True) -> np.ndarray:
        """Push both AWG channels through the link.  Returns scope samples."""
        n = min(len(awg_x), len(awg_w))
        x, w = np.asarray(awg_x[:n], float), np.asarray(awg_w[:n], float)

        # Drivers: gain, soft saturation, then the source-to-termination halving.
        vx = self._soft_clip(self.a_v_x * x, self.sw_x) / 2.0
        vw = self._soft_clip(self.a_v_w * w, self.sw_w) / 2.0

        # Everything ahead of the sine.
        vx = _fft_filter(vx, self._h_pre, AWG_FS)
        vw = _fft_filter(vw, self._h_pre, AWG_FS)

        # Two quadrature-biased Mach-Zehnders in series, Methods equations 6-7.
        alpha = 10.0 ** (-self.il_dB / 10.0)
        inv_er = 10.0 ** (-self.er_dB / 10.0)
        def t(v):
            return alpha * ((1.0 - inv_er)
                            * (1.0 + np.sin(np.pi * v / self.v_pi)) / 2.0 + inv_er)
        p_opt = self.p_in_w * t(vx) * t(vw) * self.split * self.g_soa

        # Detection, then the noise that is generated with the carriers.
        i_ph = self.responsivity * p_opt
        if noise:
            r2 = self.responsivity ** 2
            s = (2 * Q_ELECTRON * self.i_ph_n
                 + 4 * r2 * self.p_sig_n * self.rho_ase_n
                 + 2 * r2 * self.rho_ase_n ** 2 * BW_OPT
                 + self.i_n_in ** 2)
            # One-sided density S over [0, fs/2] is a per-sample variance of
            # S*fs/2 — the same zero-order-hold convention the simulator uses.
            i_ph = i_ph + self.rng.normal(0.0, np.sqrt(s * AWG_FS / 2.0), n)

        v = _fft_filter(self.z_t * i_ph, self._h_post, AWG_FS)
        v = self._soft_clip(v, self.rx_swing)
        if quantise_rx:
            v = quantise(v - np.mean(v), enob_for_baud(baud))
        return v


# ── the annealing noise source ──────────────────────────────────────────────
def ase_density(i_mA: float, p_ase_ref_dBm: float = -38.0,
                i_th_mA: float = 20.0, i_ref_mA: float = 100.0,
                n_exp: float = 1.5, bw_ase_nm: float = 40.0,
                lambda_nm: float = 1310.0) -> float:
    """Bulk SOA bias current to injected ASE density, W/Hz.

    The same arithmetic as `models/va_ase_source.va`, repeated here so that the
    annealing schedule can be evaluated without starting a simulator.  Keep the
    two in step: the model is the one that runs in the deck.
    """
    x = max((i_mA - i_th_mA) / (i_ref_mA - i_th_mA), 0.0)
    p = 1e-3 * 10.0 ** (p_ase_ref_dBm / 10.0) * x ** n_exp
    lam = lambda_nm * 1e-9
    return p / (299792458.0 * bw_ase_nm * 1e-9 / lam ** 2)


def anneal_current(t: int, i0_mA: float = 100.0, gamma: float = 0.0235,
                   i_th_mA: float = 20.0, n_exp: float = 1.5) -> float:
    """Bias current at iteration `t` for the paper's exponential schedule.

    Methods gives the schedule as a variance: eta^2 = eta0^2 * exp(-gamma*t),
    with eta0^2 = 18.2 mV^2 and gamma = 0.0235 per iteration, starting from a
    100 mA bias.  The variance at the oscilloscope is proportional to the
    injected ASE power, which goes as (I - I_th)^n_exp, so the current that
    realises the schedule is

        I(t) = I_th + (I0 - I_th) * exp(-gamma*t/n_exp)

    Clamped at the threshold, below which the source emits nothing.
    """
    return i_th_mA + (i0_mA - i_th_mA) * np.exp(-gamma * t / n_exp)


# ── validation ──────────────────────────────────────────────────────────────
def validate(n_sym: int = 96, baud: float = 106e9, rolloff: float = 0.2,
             seed: int = 0, oversample: tuple[int, ...] = (1, 2, 4, 8, 16)
             ) -> list[dict]:
    """Push one waveform through both models and report the difference.

    This is the check that makes the fast model usable.  Without it, the numpy
    version is a plausible story about the deck rather than a stand-in for it.
    Noise is off on both sides: a noisy comparison measures two different random
    sequences, not two different models.

    IT SWEEPS THE TIMESTEP, and that is the whole point.  A transient run at one
    step per AWG sample, 3.9 ps, does not resolve this link's own poles: the
    fastest is at 127 GHz, whose time constant is 1.25 ps, and backward Euler at
    three time constants per step damps hard.  The deck then swings 35 % less
    than it should, and comparing against it at that step would "prove" the fast
    model wrong by exactly that much.

    Refine the step and the deck walks toward the fast model.  The `scale`
    column is the amplitude ratio the two agree on, and it should approach 1.
    If it does not, the fast model is wrong.  If it does, the coarse-step run
    was.

    Read this before choosing a step for any transient in this directory:
    16 samples per AWG sample is the honest setting, and it is 16 times the
    cost.  It is also why the iteration loop uses this file and not the deck.
    """
    import dsp

    rng = np.random.default_rng(seed)
    sx = np.sign(rng.normal(size=n_sym))
    sw = rng.uniform(-1.0, 1.0, n_sym)
    ix, iw = dsp.interleave(sx, sw)
    awg_x, awg_w = dsp.tx(ix, iw, baud, rolloff, n_pilot=64)

    lk = Link()
    ch = Channel.from_link(lk, seed=seed)
    t = np.arange(len(awg_x)) / AWG_FS
    v_fast = ch(awg_x, awg_w, baud, noise=False, quantise_rx=False)

    out = []
    for over in oversample:
        t_sim, v_sim = lk.tran(t, awg_x, awg_w, tstep=1.0 / (over * AWG_FS))
        b_all = np.interp(t_sim, t, v_fast)
        a = v_sim - np.mean(v_sim)
        b = b_all - np.mean(b_all)
        m = slice(len(a) // 8, len(a) - len(a) // 8)
        best = None
        for lag in range(-6 * over, 6 * over + 1):
            bb = np.roll(b, lag)
            s = float(np.dot(a[m], bb[m]) / np.dot(bb[m], bb[m]))
            e = float(np.std(a[m] - s * bb[m]) / np.std(a[m]))
            if best is None or e < best[0]:
                best = (e, lag, s)
        out.append({"oversample": over, "n_samples": len(t_sim),
                    "rms_deck": float(np.std(a)), "rms_fast": float(np.std(b)),
                    "rel_error": best[0], "lag": best[1], "scale": best[2]})
    return out


if __name__ == "__main__":
    print("ASE density and the annealing schedule:")
    for it in (0, 25, 50, 100, 200):
        i = anneal_current(it)
        print(f"  iteration {it:4d}  bias {i:6.2f} mA  "
              f"rho {ase_density(i):.3e} W/Hz")

    print("\nfast model against the deck, no noise on either side.")
    print("The deck should walk toward the fast model as the step refines.\n")
    print(f"  {'step/AWG':>9} {'points':>8} {'deck rms':>9} {'fast rms':>9}"
          f" {'scale':>7} {'rel err':>8}")
    rows = validate()
    for r in rows:
        print(f"  {r['oversample']:9d} {r['n_samples']:8d} {r['rms_deck']:9.5f}"
              f" {r['rms_fast']:9.5f} {r['scale']:7.4f} {r['rel_error']:8.4f}")
    last = rows[-1]
    assert abs(last["scale"] - 1.0) < 0.06, "the fast model does not converge"
    assert last["rel_error"] < 0.12, "residual too large at the finest step"
    print("\nchannel self-check passed")
