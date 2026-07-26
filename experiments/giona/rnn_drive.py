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

def solve_bias(W, powers, v_star, n_active=2, iters=12, verbose=True):
    """Per-neuron PD_B so every neuron rests at v_star WITH the weights applied.

    The rule: choose the rest state for all neurons at once, and the required
    bias follows — no iteration on the recurrence itself. Here the residual
    algebra (2k/10k network + clamping diode) is mopped up with a secant, the
    slope re-measured each step because the diode stiffens as it turns on.
    """
    pdb = np.full(N, -7.0)
    for it in range(iters):
        v = op(W, pdb, powers)
        err = v[:n_active] - v_star
        if verbose:
            print(f"    bias iter {it}: V={np.round(v[:n_active],5)} "
                  f"max|err|={np.abs(err).max()*1e3:7.3f} mV", flush=True)
        if np.abs(err).max() < 1e-3:
            break
        v2 = op(W, pdb - 0.25, powers)
        slope = (v - v2) / 0.25
        slope = np.where(np.abs(slope) > 1e-4, slope, 1/6)
        # Damped + step-limited: the undamped secant rings, because the diode's
        # incremental stiffness changes faster than the secant's memory of it.
        step = -0.65 * err / slope[:n_active]
        step = np.clip(step, -2.0, 2.0)
        pdb[:n_active] = pdb[:n_active] + step
    return pdb, op(W, pdb, powers)

def transient(W, pdb, powers, step, stop, kick):
    c = fc.Circuit(); c.load_str(src(W, pdb, powers, tran=(step, stop), kick=kick))
    t0 = time.time()
    r = c.run("tran", step=step, stop=stop, method="gear", variable_step=True)
    t = np.asarray(r.time())
    v = np.array([np.asarray(r[f"V(mod_cathode{i+1})"]) for i in range(N)])
    bus = np.array([np.asarray(r[f"V(bout8_re_{k})"])**2 + np.asarray(r[f"V(bout8_im_{k})"])**2
                    for k in range(N)]) / 1e-3
    return t, v, bus, time.time() - t0
