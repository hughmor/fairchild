"""Winner-take-all on the giona RNN.

W_rec = [[1,-1],[-1,1]]: row sums both zero (so no counter-bias is needed) and
eigenvalues {2, 0}. The λ=2 eigenvector is the DIFFERENTIAL mode (1,-1) and the
λ=0 one is the common mode (1,1), so once 2·G > 1 the differential mode is
unstable and the pair latches with one neuron high and the other low. Which way
it falls is decided by the input.

Protocol (the ordering matters):
  1. Solve the per-neuron bias with the input weights at ZERO, so both neurons
     rest at the same point. Solving with the input applied would have the bias
     rule cancel the very asymmetry the input is meant to supply.
  2. Apply the input weights (w13, w23) on channel 3, and present channel 3 as a
     step at t = 5 ns (its MZM starts at v_pi = off, steps to 0 V = on).
  3. Watch V1 - V2 latch. sign(w13 - w23) should pick the winner.
"""
import sys, numpy as np, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import rnn_drive as rnn
import fairchild as fc
SC = str(Path(__file__).resolve().parent / "results")
REC = np.array([[1.0,-1.0],[-1.0,1.0]])

def W_of(in_w):
    W = np.zeros((8,8)); W[:2,:2] = REC
    W[0,2], W[1,2] = in_w
    return W

def input_step(src_text):
    """Channel 3 off -> on at t=5 ns. fc_mzm: 0 V = ON, v_pi = 1 V = OFF."""
    tgt = "Vd3 d3 0 DC 0"
    assert tgt in src_text
    return src_text.replace(tgt, "Vd3 d3 0 PULSE(1 0 5n 200p 200p 1 2)")

def run(tag, p, v_star, in_w, stop=80e-9, step=200e-12):
    powers = [p, p, p] + [0.0]*5
    # 1. bias with NO input coupling
    pdb, v0 = rnn.solve_bias(W_of((0.0, 0.0)), powers, v_star, verbose=False)
    # 2/3. apply the input weights and step channel 3 on
    W = W_of(in_w)
    src = rnn.src(W, pdb, powers, tran=(step, stop))
    src = input_step(src)
    c = fc.Circuit(); c.load_str(src)
    t0 = time.time()
    r = c.run("tran", step=step, stop=stop, method="gear", variable_step=True)
    el = time.time()-t0
    t = np.asarray(r.time())
    v = np.array([np.asarray(r[f"V(mod_cathode{i+1})"]) for i in range(8)])
    bus = np.array([np.asarray(r[f"V(bout8_re_{k})"])**2 + np.asarray(r[f"V(bout8_im_{k})"])**2
                    for k in range(8)])/1e-3
    d0, d1 = float(v[0][0]-v[1][0]), float(v[0][-1]-v[1][-1])
    win = "N1" if d1 < -1e-4 else ("N2" if d1 > 1e-4 else "tie")
    print(f"{tag:20s} in_w=({in_w[0]:.2f},{in_w[1]:.2f})  V0={np.round(v0[:2],4)}  "
          f"(V1-V2): {d0*1e3:+8.3f} -> {d1*1e3:+8.3f} mV   winner={win}   "
          f"bus λ1={bus[0][-1]:.4f} λ2={bus[1][-1]:.4f} mW  {el:5.1f}s", flush=True)
    np.savez(f"{SC}/{tag}.npz", t=t, v=v, bus=bus, pdb=pdb, p=p, v_star=v_star,
             in_w=np.array(in_w))
    return d1

if __name__ == "__main__":
    for p in (30.0, 60.0):
        print(f"\n=== WTA {p:.0f} mW/channel, rest V*=-0.98 ===", flush=True)
        for in_w in ((0.8, 0.0), (0.0, 0.8), (0.5, 0.3), (0.3, 0.5), (0.4, 0.4)):
            try:
                run(f"w3_{int(p)}_{int(in_w[0]*10)}{int(in_w[1]*10)}", p, -0.98, in_w)
            except Exception as e:
                print(f"  in_w={in_w} FAILED {str(e)[:52]}", flush=True)
