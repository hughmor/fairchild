"""Netlist generator for the hypermultiplexed photonic tensor core.

Builds the architecture of Pappas et al. as analysed in the scaling manuscript
(`ptc-scaling-analysis/manuscript.md`, sec:architecture), with K = N:

    N lasers -> mux -> 1:N fan-out -> N MZMs (weights W)
             -> N x N AWGR (cyclic route)
             -> 1:S fan-out -> S*N MZMs (inputs X)
             -> demux -> N^2 S photodiodes + integrating capacitors

Receiver (j, s, k) accumulates the weight row i = (j - k) mod N against input
column (j, s), so a block of L symbols yields Y[i, s, j] = sum_l W[i,l] X[l,s,j]
-- the matrix-tensor product of eq:throughput.

The integrator is a capacitor on the photodiode anode: charge, not current, is
the observable, and a block result is the *difference* of the node voltage
across the block. No reset device, no op-amp; the leak is r_shunt * C_int.

Noise (`Hw.noise=True`) is injected as pre-sampled PWL sources, since fairchild
has no random-waveform primitive. TIA + shot noise share one independent
unit-variance node per receiver; RIN gets one node per *wavelength*, shared by
every receiver reading that laser, because RIN is a property of the source.
ASE is absent -- see README.

    python3 experiments/ptc/ptc.py --selftest
"""
from __future__ import annotations

import argparse
import math
from dataclasses import dataclass, field

import numpy as np

C0 = 299_792_458.0
Q_E = 1.602176634e-19


@dataclass
class Arch:
    """Scaling parameters. K = N throughout (manuscript sec:architecture)."""

    N: int = 4  # wavelengths == AWGR ports == first-bank modulators
    S: int = 1  # spatial fan-out after the AWGR
    L: int = 4  # symbols accumulated per readout (paper baseline L = N)
    f_B: float = 25e9  # symbol rate


@dataclass
class Hw:
    """Hardware, defaulting to scenario B of tbl:hardware_params."""

    p_laser_mW: float = 10.0  # per line, <= P_0max = 13 dBm
    lambda0_nm: float = 1550.0
    df_ghz: float = 100.0  # channel spacing; B_opt/N in the paper
    fwhm_ghz: float = 40.0  # AWG passband
    m: float = 0.8  # modulation index
    v_pi: float = 3.0
    responsivity: float = 0.8
    c_int_fF: float = 500.0
    l_mux_db: float = 0.5
    l_awgr_db: float = 2.0
    l_demux_db: float = 1.5
    l_mzm_db: float = 1.0  # per stage, excluding the 3 dB quadrature bias
    l_cpl_db: float = 1.0  # chip coupling charged per modulator stage
    xt_adj_db: float = -30.0  # AWGR port isolation, adjacent
    xt_bg_db: float = -40.0  # AWGR port isolation, background
    er_db: float = 30.0  # MZM extinction ratio
    r_shunt: float = 1e15  # integrator leak path
    # --- noise: `.options trannoise=1`, so shot and RIN come from the devices ---
    noise: bool = False
    noise_seed: int = 1
    s_i_tia: float = 18e-12  # A/sqrt(Hz), input-referred; no TIA device, see below
    rin_db_hz: float = -155.0
    # --- knobs the experiments sweep ---
    predistort: bool = True  # arcsin drive -> exactly linear MZM
    tr_frac: float = 0.25  # NRZ edge as a fraction of the symbol period
    oversample: int = 8  # transient steps per symbol


@dataclass
class Deck:
    text: str
    arch: Arch
    hw: Hw
    n_blocks: int
    t_edges: np.ndarray = field(repr=False)  # block boundaries, seconds
    step: float = 0.0
    stop: float = 0.0


def grid_nm(hw: Hw, n: int) -> np.ndarray:
    """Channel wavelengths on the device's *frequency* grid."""
    f0 = C0 / (hw.lambda0_nm * 1e-9)
    return np.array([C0 / (f0 + k * hw.df_ghz * 1e9) * 1e9 for k in range(n)])


def _drive(hw: Hw, u: np.ndarray) -> np.ndarray:
    """Symbol values in [-1, 1] -> MZM drive volts about quadrature.

    fc_mzm transmits (alpha/2)(1 + sin(pi*dV/V_pi * -1)) about V_pi/2, so a
    linear drive gives 1 + sin(m u) and an arcsin drive gives exactly 1 + m u.
    """
    arg = np.clip(hw.m * u, -1.0, 1.0)
    delta = np.arcsin(arg) if hw.predistort else arg
    return hw.v_pi / 2.0 - (hw.v_pi / math.pi) * delta


def _pwl(t_sym: float, vals: np.ndarray, tr: float) -> str:
    """NRZ piecewise-linear waveform: symbol l held over [l*T, (l+1)*T)."""
    pts = [(0.0, vals[0])]
    for l in range(1, len(vals)):
        t = l * t_sym
        pts += [(t - tr / 2, vals[l - 1]), (t + tr / 2, vals[l])]
    pts.append((len(vals) * t_sym, vals[-1]))
    return "PWL(" + " ".join(f"{t:.9e} {v:.9e}" for t, v in pts) + ")"


def _fanout(lines: list[str], src: str, m: int, tag: str, n_ch: int) -> list[str]:
    """1 -> m even split as a splitter chain (any m, not just powers of two)."""
    if m == 1:
        return [src]
    outs, rest = [], src
    for t in range(m - 1):
        out = f"{tag}{t}"
        nxt = f"{tag}r{t}" if t < m - 2 else f"{tag}{m - 1}"
        lines += [f".optical_port {out} {n_ch}", f".optical_port {nxt} {n_ch}"]
        lines.append(f"X{tag}s{t} {rest} {out} {nxt} fc_splitter r={1.0 / (m - t):.12g}")
        outs.append(out)
        rest = nxt
    outs.append(rest)
    return outs


def build(arch: Arch, hw: Hw, W: np.ndarray, X: np.ndarray, seed: int = 0) -> Deck:
    """W is (N, L_tot); X is (L_tot, S, N). Both in [-1, 1]."""
    N, S, L = arch.N, arch.S, arch.L
    l_tot = W.shape[1]
    assert W.shape == (N, l_tot) and X.shape == (l_tot, S, N)
    assert l_tot % L == 0, "symbol count must be a whole number of blocks"
    n_blocks = l_tot // L

    t_sym = 1.0 / arch.f_B
    tr = hw.tr_frac * t_sym
    dt = t_sym / hw.oversample
    # An unresolved symbol edge is indistinguishable from a slower edge: the
    # solver smears the transition over one step, so the MAC error stops
    # tracking tr and starts tracking dt. Four steps per edge keeps the ISI
    # term (below) physical rather than numerical.
    assert tr >= 4 * dt, (
        f"tr = {tr:.3g} s needs >= 4 solver steps of {dt:.3g} s; "
        f"raise oversample above {4 / hw.tr_frac:.0f}")
    t_stop = l_tot * t_sym
    lam = grid_nm(hw, N)
    mzm = (f"V_pi={hw.v_pi} alpha_dB={hw.l_mzm_db + hw.l_cpl_db:.6g}"
           f" e_r_dB={hw.er_db:.6g}")

    ln = [f"* hypermultiplexed PTC  N={N} S={S} K={N} L={L} f_B={arch.f_B:.4g}",
          f".options lambda_center_nm={hw.lambda0_nm} vntol=1e-14 reltol=1e-12",
          f".options max_step={dt:.6e} method=gear"]
    if hw.noise:
        ln.append(f".options trannoise=1 noiseseed={hw.noise_seed}"
                  " variable_step=0")

    # --- source bank -> WDM bus ---
    for k in range(N):
        ln.append(f".optical_port ch{k}")
        # RIN is injected at the source, so it propagates through both
        # modulators and arrives at every receiver of wavelength k perfectly
        # correlated -- which is what a shared laser physically does.
        rin = f" rin_db_hz={hw.rin_db_hz:.6g}" if hw.noise else ""
        ln.append(f"Xl{k} ch{k} fc_cw_laser power_mW={hw.p_laser_mW:.9g}"
                  f" wavelength_nm={lam[k]:.9f}{rin}")
    ln.append(".optical_port src %d" % N)
    ln.append("Xmux src " + " ".join(f"ch{k}" for k in range(N))
              + f" fc_mux il_db={hw.l_mux_db:.6g}")

    # --- 1:N fan-out, weight bank ---
    branches = _fanout(ln, "src", N, "fw", N)
    for i in range(N):
        ln += [f".optical_port a{i} {N}",
               f"Vw{i} vw{i} 0 " + _pwl(t_sym, _drive(hw, W[i]), tr),
               f"Xm{i} {branches[i]} a{i} vw{i} 0 fc_mzm {mzm}"]

    # --- AWGR ---
    for j in range(N):
        ln.append(f".optical_port b{j} {N}")
    ln.append("Xawgr " + " ".join(f"a{i}" for i in range(N))
              + " " + " ".join(f"b{j}" for j in range(N))
              + f" fc_awgr df_ghz={hw.df_ghz:.6g} fwhm_ghz={hw.fwhm_ghz:.6g}"
              f" lambda0_nm={hw.lambda0_nm:.9g} il_db={hw.l_awgr_db:.6g}"
              f" xt_adj_db={hw.xt_adj_db:.6g} xt_bg_db={hw.xt_bg_db:.6g}")

    # --- spatial fan-out, input bank, demux, receivers ---
    for j in range(N):
        taps = _fanout(ln, f"b{j}", S, f"fs{j}", N)
        for s in range(S):
            ln += [f".optical_port d{j}_{s} {N}",
                   f"Vx{j}_{s} vx{j}_{s} 0 " + _pwl(t_sym, _drive(hw, X[:, s, j]), tr),
                   f"Xn{j}_{s} {taps[s]} d{j}_{s} vx{j}_{s} 0 fc_mzm {mzm}"]
            for k in range(N):
                ln.append(f".optical_port p{j}_{s}_{k}")
            ln.append(f"Xdm{j}_{s} d{j}_{s} "
                      + " ".join(f"p{j}_{s}_{k}" for k in range(N))
                      + f" fc_demux il_db={hw.l_demux_db:.6g}"
                      f" fwhm_ghz={hw.fwhm_ghz:.6g} df_ghz={hw.df_ghz:.6g}"
                      f" lambda0_nm={hw.lambda0_nm:.9g}")
            for k in range(N):
                q, p = f"q{j}_{s}_{k}", f"p{j}_{s}_{k}"
                ln += [f"Xpd{j}_{s}_{k} {p} {q} 0 fc_photodetector"
                       f" responsivity={hw.responsivity:.6g}"
                       f" r_shunt={hw.r_shunt:.6g} i_dark=0",
                       f"Ci{j}_{s}_{k} {q} 0 {hw.c_int_fF:.6g}f IC=0"]
                if hw.noise:
                    ln += _tia_noise(hw, j, s, k)

    ln += [f".tran {dt:.6e} {t_stop:.6e}", ".end"]
    d = Deck("\n".join(ln) + "\n", arch, hw, n_blocks,
             np.arange(n_blocks + 1) * L * t_sym)
    d.step, d.stop = dt, t_stop
    return d


KB = 1.380649e-23


def _tia_noise(hw: Hw, j: int, s: int, k: int) -> list[str]:
    """Input-referred TIA current noise.

    Shot and RIN come from `fc_photodetector` and `fc_cw_laser` under
    `.options trannoise=1`. The TIA does not: there is no TIA device, and no
    element that spells "a current source of PSD S". A resistor of
    `R = 4kT/S_i²` would have the right PSD but loads the integrator to death
    (51 ohm across 500 fF is a 25 fs time constant).

    So: pick any convenient R, and buffer its 4kT·R volts through a B-element
    transconductance chosen to land the wanted current PSD. The resistor sets
    the spectrum, the B-element sets the amplitude, and nothing loads the
    integrating node.
    """
    r_n = 1e3
    g = hw.s_i_tia / math.sqrt(4.0 * KB * 300.15 * r_n)
    tag = f"{j}_{s}_{k}"
    return [f"Rtia{tag} ntia{tag} 0 {r_n:.6g}",
            f"Btia{tag} 0 q{tag} I={g:.9e}*V(ntia{tag})"]


def readout(deck: Deck, res) -> np.ndarray:
    """Integrator voltage differences per block -> (n_blocks, N, S, N) charge.

    Indexed [block, i, s, j] with i = (j - k) mod N: the weight row a receiver
    actually accumulated, not its wavelength slot.
    """
    N, S = deck.arch.N, deck.arch.S
    t = res.time()
    out = np.empty((deck.n_blocks, N, S, N))
    for j in range(N):
        for s in range(S):
            for k in range(N):
                v = np.interp(deck.t_edges, t, res[f"V(q{j}_{s}_{k})"])
                out[:, (j - k) % N, s, j] = np.diff(v) * deck.hw.c_int_fF * 1e-15
    return out


def decode(deck: Deck, charge: np.ndarray, cal: np.ndarray,
           W: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Charge -> MAC estimate, removing the deterministic linear MZM terms.

    p = P_det (1 + m w)(1 + m x) leaves `1`, `m w` and `m x` alongside the
    wanted `m^2 w x`. All three are known at calibration time (the weights are
    programmed and the inputs are the data), so they are subtracted rather than
    filtered -- eq:cascaded_mzm and the note under it.
    """
    N, S, L = deck.arch.N, deck.arch.S, deck.arch.L
    m = deck.hw.m
    kappa = cal / L  # charge per symbol at zero drive
    wb = W.reshape(N, -1, L).sum(axis=2)  # (N, n_blocks)
    xb = X.reshape(-1, L, S, N).sum(axis=1)  # (n_blocks, S, N)
    lin = m * (wb.T[:, :, None, None] + xb[:, None, :, :])
    return (charge / kappa[None] - L - lin) / m**2


def ideal(deck: Deck, W: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Y[b, i, s, j] = sum over the block of W[i,l] X[l,s,j]."""
    N, S, L = deck.arch.N, deck.arch.S, deck.arch.L
    w = W.reshape(N, -1, L)
    x = X.reshape(-1, L, S, N)
    return np.einsum("ibl,blsj->bisj", w, x)


def enob(err: np.ndarray, span: float = 1.0) -> float:
    """Operand-referenced ENOB (eq:enob) from the MAC residual.

    2^ENOB = FS / (sqrt(12) sigma_n) with FS = 2*span the full-scale *range* of
    the quantity measured. One MAC of operands on [-1,1] has span 1; a block of
    L accumulated MACs has span L, which is where the +0.5*log2(L) of the
    output-referred convention comes from -- do not also divide by sigma_d.
    """
    return math.log2(span / (math.sqrt(3.0) * err.std()))


def enob_theory(hw: Hw, arch: Arch, i_ph: float,
                rin_instantaneous: bool = True) -> float:
    """Block ENOB from eq:master_snr / eq:coefficients, no ASE (no SOA yet).

    The integrator is a boxcar of length 1/f_B, whose equivalent noise bandwidth
    is f_B/2 -- half the bandwidth the manuscript's per-symbol convention uses,
    hence the sqrt(2). The L-symbol accumulation then adds 0.5*log2(L).

    `rin_instantaneous` inflates the RIN term by (1 + m^2<u^2>)^2. eq:coefficients
    writes `r I_ph^2 B` against the *average* photocurrent, but RIN is
    multiplicative on the instantaneous power and two cascaded modulators raise
    <P^2>/<P>^2 by that factor. Set False to compare against the paper as
    written; the simulator does the instantaneous thing whichever you pick.
    """
    b = arch.f_B
    a0 = hw.s_i_tia**2 * b
    a1 = 2.0 * Q_E * b
    a2 = 10.0 ** (hw.rin_db_hz / 10.0) * b
    if rin_instantaneous:
        a2 *= (1.0 + hw.m**2 / 3.0) ** 2  # uniform operands, <u^2> = 1/3
    sigma = math.sqrt((a0 + a1 * i_ph + a2 * i_ph**2) / 2.0)
    per_mac = math.log2(hw.m**2 * i_ph / (math.sqrt(3.0) * sigma))
    return per_mac + 0.5 * math.log2(arch.L)


def run(deck: Deck, method: str = "tr"):
    """Solve the transient.

    `method="tr"` matters: the observable is the integral of a piecewise-linear
    photocurrent, which the trapezoidal rule reproduces exactly. BE and GEAR are
    first-order on it and leak charge at every symbol edge, which shows up as a
    fake ~1/oversample error floor in the MAC.
    """
    import fairchild as fc

    c = fc.Circuit()
    c.load_str(deck.text)
    # uic is not optional: the integrator node sees only C_int and r_shunt,
    # so a DC operating point would put I_ph * 1e15 ohm on it and never
    # converge. The physical initial state is a discharged integrator.
    # trannoise also needs a fixed step: the LTE controller would read a fresh
    # noise sample as a fast signal and shrink the step to chase it, correlating
    # step size with noise. The solver refuses the combination outright.
    return c.run("tran", step=deck.step, stop=deck.stop, method=method,
                 max_step=deck.step, uic=True, variable_step=False,
                 trannoise=deck.hw.noise, noiseseed=deck.hw.noise_seed)


def calibrate(arch: Arch, hw: Hw, W: np.ndarray, X: np.ndarray, seed: int = 0):
    """Zero-drive reference charge per block, (N, S, N). Noise sources off."""
    cal_hw = Hw(**{**hw.__dict__, "noise": False})
    deck = build(arch, cal_hw, np.zeros_like(W), np.zeros_like(X), seed)
    return readout(deck, run(deck))[0]


def simulate(arch: Arch, hw: Hw, W: np.ndarray, X: np.ndarray, seed: int = 0):
    """Build, solve, calibrate, decode. Returns (measured, ideal) MAC tensors."""
    deck = build(arch, hw, W, X, seed)
    charge = readout(deck, run(deck))
    cal = calibrate(arch, hw, W, X, seed)
    return decode(deck, charge, cal, W, X), ideal(deck, W, X)


IDEAL = dict(fwhm_ghz=0.0, l_mux_db=0, l_awgr_db=0, l_demux_db=0,
             l_mzm_db=0, l_cpl_db=0, er_db=100)


def selftest() -> None:
    arch = Arch(N=4, S=1, L=4, f_B=25e9)
    rng = np.random.default_rng(1)
    n_sym = 3 * arch.L
    W = rng.uniform(-1, 1, (arch.N, n_sym))
    X = rng.uniform(-1, 1, (n_sym, arch.S, arch.N))

    # With no loss, no passband and a predistorted drive, the only thing left
    # between the sim and the exact MAC is the symbol edge: two operands ramping
    # together across a boundary integrate to (tr/6)(dW)(dX) less than the
    # rectangular product. That is ISI, which the manuscript defers, so it must
    # scale with tr and not with the timestep.
    errs = []
    for tr_frac, os_ in [(0.08, 100), (0.04, 200), (0.02, 400)]:
        hw = Hw(**IDEAL, tr_frac=tr_frac, oversample=os_)
        got, want = simulate(arch, hw, W, X)
        errs.append(np.abs(got - want).max())
        print(f"  tr/T = {tr_frac:<5} max |Y_sim - Y_exact| = {errs[-1]:.3e}")
    assert errs[0] / errs[2] > 3.5, f"ISI error is not tracking tr: {errs}"
    assert errs[2] < 2e-2, errs[2]

    hw = Hw(**IDEAL, tr_frac=0.02, oversample=400)
    _, want = simulate(arch, hw, W, X)
    # The routing map is the claim: receiver (j,s,k) must carry weight row
    # (j-k) mod N, so permuting the weights must break the answer.
    got2, _ = simulate(arch, hw, W[::-1], X)
    assert np.abs(got2 - want).max() > 0.5, "readout index map is not being tested"

    hw_xt = Hw(**{**hw.__dict__, "fwhm_ghz": 40.0, "xt_adj_db": -15.0})
    got3, want3 = simulate(arch, hw_xt, W, X)
    e_xt = np.abs(got3 - want3).max()
    print(f"  AWGR port isolation -15 dB: max err = {e_xt:.3e}")
    assert e_xt > 10 * errs[2], (e_xt, errs[2])

    # Noise must land where eq:master_snr says it does. Use a long run so the
    # sample std is worth comparing, and read I_ph off the calibration charge
    # rather than re-deriving the link budget by hand.
    n_long = 40 * arch.L
    Wl = rng.uniform(-1, 1, (arch.N, n_long))
    Xl = rng.uniform(-1, 1, (n_long, arch.S, arch.N))
    hw_n = Hw(**{**hw.__dict__, "noise": True})
    got4, want4 = simulate(arch, hw_n, Wl, Xl)
    i_ph = calibrate(arch, hw_n, Wl, Xl)[0, 0, 0] / arch.L * arch.f_B
    got_i, want_i = simulate(arch, hw, Wl, Xl)  # same deck, noise off
    resid = (got4 - want4) - (got_i - want_i)  # strip the ISI term
    measured, predicted = enob(resid, span=arch.L), enob_theory(hw_n, arch, i_ph)
    print(f"  shot + TIA + RIN, I_ph = {i_ph * 1e6:.1f} uA:"
          f"  ENOB {measured:.2f} b vs theory {predicted:.2f} b")
    assert abs(measured - predicted) < 0.35, (measured, predicted)
    print("selftest OK")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("-N", type=int, default=4)
    ap.add_argument("-S", type=int, default=1)
    ap.add_argument("-L", type=int, default=4)
    ap.add_argument("--blocks", type=int, default=2)
    ap.add_argument("--noise", action="store_true")
    ap.add_argument("-o", "--out", help="write the netlist here instead of stdout")
    a = ap.parse_args()
    if a.selftest:
        selftest()
    else:
        arch = Arch(N=a.N, S=a.S, L=a.L)
        rng = np.random.default_rng(0)
        n = a.blocks * a.L
        deck = build(arch, Hw(noise=a.noise),
                     rng.uniform(-1, 1, (a.N, n)),
                     rng.uniform(-1, 1, (n, a.S, a.N)))
        if a.out:
            open(a.out, "w").write(deck.text)
            print(f"wrote {a.out} ({len(deck.text) / 1e3:.0f} kB)")
        else:
            print(deck.text, end="")
