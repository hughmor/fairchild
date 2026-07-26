"""rnn_drive.py — drive the hand-written giona RNN deck from Python.

Read-only on netlists/giona_rnn_perfectW.sp: `.param` values are substituted at
parse time, so they cannot be reached with set_param, and overriding the deck
TEXT is the honest way to sweep them.

  src()         build the deck source with .param overrides, an optional .tran,
                and an optional MZM-drive stimulus
  op()          solve the DC operating point, return the 8 neuron node voltages
  solve_bias()  the recurrent bias rule (see below)
  transient()   run a transient, return (t, V[8], bus[8], wall seconds)
"""
import fairchild as fc, numpy as np, re, pathlib, time
NET = pathlib.Path(__file__).resolve().parent / "netlists" / "giona_rnn_perfectW.sp"
CAL=[7.9998737e-06,8.0074122e-06,8.0149580e-06,8.0224153e-06,8.0300681e-06,8.0376334e-06,8.0452040e-06,8.0527814e-06]
N=8

def src(W, pdb, powers, tran=None, kick=None):
    s = NET.read_text()
    o = {f"radius{i+1}": CAL[i] for i in range(N)}
    o |= {f"PD_B{i+1}": float(pdb[i]) for i in range(N)}
    o |= {f"p_{i+1}": float(powers[i]) for i in range(N)}
    for i in range(N):
        for j in range(N):
            o[f"w_{i+1}{j+1}"] = float(W[i, j])
    for k, v in o.items():
        s, n = re.subn(rf"^(\.param\s+{re.escape(k)}\s*=\s*)\S+",
                       lambda m: f"{m.group(1)}{v:.9g}", s, flags=re.M)
        assert n, k
    if kick:
        ch, amp, t0, wid = kick
        tgt = f"Vd{ch} d{ch} 0 DC 0"
        assert tgt in s, tgt
        s = s.replace(tgt, f"Vd{ch} d{ch} 0 PULSE(0 {amp:.4g} {t0:.3e} 50p 50p {wid:.3e} 1)")
    if tran:
        step, stop = tran
        s = s.replace("\n.op\n", f"\n.options method=gear\n.tran {step:.3e} {stop:.3e}\n")
    return s

def op(W, pdb, powers):
    c = fc.Circuit(); c.load_str(src(W, pdb, powers)); r = c.run("op")
    return np.array([float(r[f"V(mod_cathode{i+1})"][0]) for i in range(N)])

def solve_bias(W, powers, v_star, n_active=2, iters=14, verbose=True):
    """Per-neuron PD_B so every active neuron rests at v_star, WITH weights on.

    The rule: choose the rest state for all neurons at once and the required bias
    follows in closed form, because then every T_j(I*) is a known constant — no
    iteration on the recurrence itself. What remains is only the residual algebra
    of the bias network (2k || 10k plus a clamping diode), and that is what this
    solves numerically.

    It uses a COUPLED Newton on the full n_active x n_active sensitivity matrix
    dV_i/dPD_B_j, not per-neuron secants. Independent secants work only while
    each bias mostly moves its own neuron; once the differential mode goes
    unstable (2*G > 1, the winner-take-all regime we actually want) the
    off-diagonal terms dominate and per-neuron updates chase their own tails —
    that showed up as a 45 mV residual at 60 mW/channel.

    Raises if the residual stays above 5 mV rather than returning whatever the
    last iterate happened to be: a silently-wrong operating point is worse than
    no answer (an early version "converged" to PD_B = +17 V, reverse-biasing the
    junctions so hard there was no gain at all).
    """
    pdb = np.full(N, -9.0)
    h = 0.2
    for it in range(iters):
        v = op(W, pdb, powers)
        err = v[:n_active] - v_star
        if verbose:
            print(f"    bias iter {it}: V={np.round(v[:n_active], 5)} "
                  f"max|err|={np.abs(err).max() * 1e3:7.3f} mV", flush=True)
        if np.abs(err).max() < 1e-3:
            break
        # Numerical sensitivity matrix: column j is the response of every active
        # neuron to neuron j's bias.
        J = np.zeros((n_active, n_active))
        for j in range(n_active):
            pj = pdb.copy()
            pj[j] -= h
            J[:, j] = (v[:n_active] - op(W, pj, powers)[:n_active]) / h
        try:
            step = -np.linalg.solve(J, err)
        except np.linalg.LinAlgError:
            step = -err / np.where(np.abs(np.diag(J)) > 1e-4, np.diag(J), 1 / 6)
        # Damp and step-limit: the diode's incremental stiffness changes faster
        # than a linearisation can track. Clamp to negative bias — a positive one
        # can never forward-bias the junction, so it cannot be a solution.
        step = np.clip(0.7 * step, -2.0, 2.0)
        pdb[:n_active] = np.clip(pdb[:n_active] + step, -40.0, -0.2)
    v = op(W, pdb, powers)
    resid = np.abs(v[:n_active] - v_star).max()
    if resid > 5e-3:
        raise RuntimeError(
            f"bias solve did not converge: |V - V*| = {resid * 1e3:.2f} mV, "
            f"PD_B={np.round(pdb[:n_active], 3)}, V={np.round(v[:n_active], 5)}")
    return pdb, v


def transient(W, pdb, powers, step, stop, kick):
    c = fc.Circuit(); c.load_str(src(W, pdb, powers, tran=(step, stop), kick=kick))
    t0 = time.time()
    r = c.run("tran", step=step, stop=stop, method="gear", variable_step=True)
    t = np.asarray(r.time())
    v = np.array([np.asarray(r[f"V(mod_cathode{i+1})"]) for i in range(N)])
    bus = np.array([np.asarray(r[f"V(bout8_re_{k})"])**2 + np.asarray(r[f"V(bout8_im_{k})"])**2
                    for k in range(N)]) / 1e-3
    return t, v, bus, time.time() - t0
