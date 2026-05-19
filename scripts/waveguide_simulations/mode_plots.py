"""Plot |E|^2 mode profiles for strip, rib, and bent rib at 1550 nm."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from wg_sim import build_strip, build_rib, solve_neff

WL = 1.55

cases = [
    ("Strip 500x220",        build_strip, np.inf, 2.44),
    ("Rib straight",         build_rib,   np.inf, 2.76),
    ("Rib bent (R=8 um)",    build_rib,   8.0,    2.76),
]

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
for ax, (title, builder, radius, ng) in zip(axes, cases):
    polys, res = builder()
    modes = solve_neff(polys, res, WL, num_modes=1, radius=radius, n_guess=ng)
    m = modes[0]
    m.plot_intensity(ax=ax)
    ax.set_title(f"{title}\nn_eff = {m.n_eff.real:.4f}")
    ax.set_aspect("equal")
plt.tight_layout()
plt.savefig("/home/claude/mode_profiles.png", dpi=130)
print("Saved mode_profiles.png")
