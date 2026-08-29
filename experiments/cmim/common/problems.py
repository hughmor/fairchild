"""problems.py — the combinatorial problems, mapped to Ising couplings.

Every mapping is Methods, "Mathematical framework", and each function names the
equation it implements.  The Ising convention throughout is the paper's:

    H = sum_{i<j} J_ij s_i s_j + sum_i h_i s_i          Methods, equation 1

with `J` symmetric, zero on the diagonal, and `s` in {-1, +1}.  `J` is a
`scipy.sparse` matrix because the paper's largest problem has 41,209 spins and
205,233 non-zero couplings: dense, that is 13 GB of zeros.

Sparsity is not an optimisation here, it is the encoding.  Methods says the
flattening step "determines the positions of the zero elements of W and removes
corresponding matrix and vector elements", so the number of non-zeros IS the
number of symbols the AWG has to send, and `scipy.sparse` counts them for free.
"""
from __future__ import annotations

import pathlib

import numpy as np
import scipy.sparse as sp

# Best-known cut values, for scoring.  Both come from Table S4 and the Figure 3
# captions.
#
# The paper contradicts itself on G22: the main text says the best-known value
# is 13,352 and the Figure 3c caption says 13,359.  Only 13,359 is consistent
# with the 99.48 % it reports for its own 13,289 cuts, and 13,359 is the
# accepted literature value, so that is what is used here.
BEST_KNOWN_CUT = {"G22": 13359, "G81": 14030}


# ── benchmark graphs ────────────────────────────────────────────────────────
def square_lattice(rows: int, cols: int | None = None) -> tuple[sp.csr_matrix,
                                                                np.ndarray]:
    """Antiferromagnetic 2D square lattice.  The checkerboard is the ground state.

    THE PAPER'S SIGN DOES NOT WORK IN THE PAPER'S OWN EQUATION.  Methods says
    "assigning antiferromagnetic coupling (J_ij = -1) to each lattice edge", and
    the ground state "corresponds to alternating neighbouring spins, forming a
    checkerboard pattern".  Under Methods equation 1, H = sum J_ij s_i s_j with
    no leading minus, J_ij = -1 is minimised by s_i s_j = +1: every neighbour
    ALIGNED, a ferromagnet, the opposite of a checkerboard.

    Which convention did they mean?  Their own max-cut mapping settles it.
    Methods equation 12 sets J_ij = k_ij / 2, positive for a positive edge
    weight, and max-cut is antiferromagnetic by definition.  So positive J is
    antiferromagnetic here, equation 1 is the convention, and the "-1" in the
    lattice paragraph is a sign slip.

    This function returns J_ij = +1, because the checkerboard is what the paper
    demonstrates and what Figure 3(a) shows.  Its ground-state energy is minus
    the edge count.

    The paper runs 10x10 up to 203x203, which is 41,209 spins and 205,233
    couplings — the 205,233 in Extended Data Table 1 counts each edge twice,
    plus the local-field column, which is how a flattened W is actually sent.
    """
    cols = cols or rows
    n = rows * cols
    idx = np.arange(n).reshape(rows, cols)
    a = np.concatenate([idx[:, :-1].ravel(), idx[:-1, :].ravel()])
    b = np.concatenate([idx[:, 1:].ravel(), idx[1:, :].ravel()])
    J = sp.coo_matrix((np.ones(2 * len(a)),
                       (np.concatenate([a, b]), np.concatenate([b, a]))),
                      shape=(n, n)).tocsr()
    return J, np.zeros(n)


def lattice_ground_energy(rows: int, cols: int | None = None) -> float:
    """Minus the edge count: every edge is satisfied by the checkerboard."""
    cols = cols or rows
    return -float(rows * (cols - 1) + (rows - 1) * cols)


def read_gset(path: str | pathlib.Path) -> tuple[sp.csr_matrix, np.ndarray]:
    """Read a Gset graph and return its max-cut Ising couplings.

    The Gset format is one header line `n_vertices n_edges`, then one line per
    edge `i j w` with one-based vertex indices.  The files are the standard
    Stanford set; the paper uses G22 (2,000 vertices, 19,990 edges) and G81
    (20,000 vertices, 40,000 edges).

    The mapping is Methods equation 12: J_ij = k_ij / 2, and no local field.
    """
    path = pathlib.Path(path)
    with path.open() as fh:
        n, m = (int(v) for v in fh.readline().split()[:2])
        rows, cols, vals = [], [], []
        for line in fh:
            if not line.strip():
                continue
            i, j, w = line.split()[:3]
            rows.append(int(i) - 1)
            cols.append(int(j) - 1)
            vals.append(float(w) / 2.0)
    if len(vals) != m:
        raise ValueError(f"{path.name}: header says {m} edges, found {len(vals)}")
    J = sp.coo_matrix((vals + vals, (rows + cols, cols + rows)),
                      shape=(n, n)).tocsr()
    return J, np.zeros(n)


def cut_value(J: sp.spmatrix, sigma: np.ndarray) -> float:
    """The cut a spin configuration achieves.  Methods equation 10.

    J carries k_ij/2 and both triangles, so the sum of J over the graph is the
    total edge weight, and C = sum_{(i,j) in E} k_ij (1 - s_i s_j) / 2 becomes
    (sum(J) - s.J.s) / 2 with this bookkeeping.
    """
    s = np.asarray(sigma, dtype=float)
    return float((J.sum() - s @ (J @ s)) / 2.0)


# ── number partitioning ─────────────────────────────────────────────────────
def number_partition(values: np.ndarray, normalise: bool = True
                     ) -> tuple[sp.csr_matrix, np.ndarray]:
    """Split a set of integers into two subsets of equal sum.

    Methods equation 18: J_ij = s_i * s_j, dense and all-to-all, with no local
    field.  The constant term of equation 17 is dropped, as the paper drops it:
    a gradient descent never sees it.

    The paper draws its integers uniformly from {0, ..., 16} and goes up to
    N = 256, which is the largest all-to-all problem the AWG memory holds.

    `normalise` divides the values by the largest of them, so the couplings land
    in [0, 1].  IT IS NOT COSMETIC, and the paper does not mention it.

    Raw integers make the problem degenerate.  With values up to 16 the
    couplings reach 256, so at the paper's beta = 0.3 the update
    x_i = alpha*s_i - beta*v_i*(v.s) is dominated by the second term.  Every
    v_i is positive, so every spin takes the sign of -(v.s) at once: the whole
    vector goes uniform, the difference goes to the total sum, and the machine
    does the worst possible thing on every iteration.  The alpha*I term that
    would break the symmetry is two orders of magnitude too small to be heard.

    Normalised, the coupling term is beta*v_i*|v.s| ~ 0.3*sqrt(N)/2, which is
    the same order as alpha, and the published hyperparameters mean something.
    A converter has a full scale and the weights are written to it, so this is
    what the hardware does whether or not the paper says so.  Measured on
    N = 16: 0 runs in 10 reach the ground state raw, and it is routine
    normalised.
    """
    v = np.asarray(values, dtype=float)
    if normalise and np.max(np.abs(v)) > 0:
        v = v / np.max(np.abs(v))
    J = np.outer(v, v)
    np.fill_diagonal(J, 0.0)
    return sp.csr_matrix(J), np.zeros(len(v))


def partition_difference(values: np.ndarray, sigma: np.ndarray) -> float:
    """|sum(S+) - sum(S-)|.  Zero is a perfect partition.  Methods equation 16."""
    return abs(float(np.asarray(values, float) @ np.asarray(sigma, float)))


def random_partition_instance(n: int, hi: int = 16, rng=None) -> np.ndarray:
    """N integers drawn uniformly from {0, ..., hi}, as Methods specifies."""
    rng = rng or np.random.default_rng()
    return rng.integers(0, hi + 1, size=n).astype(float)


# ── HP lattice protein folding ──────────────────────────────────────────────
def hp_qubo(sequence: str, sites: int | None = None,
            lam_site: float = 4.0, lam_chain: float = 4.0,
            e_hh: float = -1.0) -> np.ndarray:
    """A position-encoded QUBO for 2D HP lattice folding.

    THIS IS NOT THE PAPER'S ENCODING, and the difference is not small.

    Methods says the formulation follows reference 19 (Perdomo-Ortiz et al.) and
    describes the variables as "whether a particular amino acid is present at a
    given lattice site", which is a position encoding.  That is what this
    builds.  But the paper reports 630 spins for its 30-residue sequence, and a
    position encoding of 30 residues over any useful 2D lattice needs far more
    than 630 binary variables.  630 = 30 x 21 exactly, which suggests a specific
    site count this function cannot guess, and the paper gives no construction.

    So: this returns a working, checkable HP folding QUBO that solves the same
    problem class, and it does NOT reproduce the paper's spin counts.  Treat any
    Figure 4 comparison as qualitative until the mapping is pinned down.  See
    ../SPECS.md, and reference 19 for the construction the authors used.

    Terms: one-hot per residue, one residue per site, chain adjacency, and an
    `e_hh` reward for every non-adjacent H-H pair on neighbouring sites.
    """
    seq = sequence.upper()
    if set(seq) - {"H", "P"}:
        raise ValueError("an HP sequence contains only H and P")
    n = len(seq)
    side = int(np.ceil(np.sqrt(n))) + 1
    sites = sites or side * side
    grid = np.array([(i // side, i % side) for i in range(sites)])
    adj = (np.abs(grid[:, None, :] - grid[None, :, :]).sum(axis=2) == 1)

    def v(res, site):
        return res * sites + site

    Q = np.zeros((n * sites, n * sites))
    for r in range(n):                      # exactly one site per residue
        for a in range(sites):
            Q[v(r, a), v(r, a)] -= lam_site
            for b in range(a + 1, sites):
                Q[v(r, a), v(r, b)] += 2.0 * lam_site
    for a in range(sites):                  # at most one residue per site
        for r in range(n):
            for s in range(r + 1, n):
                Q[v(r, a), v(s, a)] += 2.0 * lam_site
    for r in range(n - 1):                  # consecutive residues are adjacent
        for a in range(sites):
            for b in range(sites):
                if not adj[a, b]:
                    Q[v(r, a), v(r + 1, b)] += lam_chain
    for r in range(n):                      # the reward the model exists for
        for s in range(r + 2, n):
            if seq[r] == "H" and seq[s] == "H":
                for a in range(sites):
                    for b in range(sites):
                        if adj[a, b]:
                            Q[v(r, a), v(s, b)] += e_hh
    return Q + Q.T - np.diag(np.diag(Q))


def qubo_to_ising(Q: np.ndarray) -> tuple[sp.csr_matrix, np.ndarray]:
    """b in {0,1} to s in {-1,+1}, via s = 2b - 1.  Methods equations 14, 15.

        J_ij = Q_ij / 4
        h_i  = Q_ii / 2 + sum_j Q_ij / 4

    The transformation adds a constant energy offset.  The paper notes that a
    gradient-descent solver never sees it, so it is dropped here too.
    """
    Q = np.asarray(Q, dtype=float)
    diag = np.diag(Q).copy()
    off = Q - np.diag(diag)
    J = off / 4.0
    h = diag / 2.0 + off.sum(axis=1) / 4.0
    return sp.csr_matrix(J), h


# ── energy ──────────────────────────────────────────────────────────────────
def ising_energy(J: sp.spmatrix, h: np.ndarray, sigma: np.ndarray) -> float:
    """H = sum_{i<j} J_ij s_i s_j + sum_i h_i s_i.  Methods equation 1.

    J holds both triangles, so the pair sum is halved.
    """
    s = np.asarray(sigma, dtype=float)
    return float(0.5 * s @ (J @ s) + h @ s)


if __name__ == "__main__":
    # Each check is against a value that can be worked out on paper.
    J, h = square_lattice(10, 10)
    assert J.nnz == 2 * 180, J.nnz
    board = np.array([(-1.0) ** (i + j) for i in range(10) for j in range(10)])
    e = ising_energy(J, h, board)
    print(f"10x10 lattice: {J.nnz//2} edges, checkerboard energy {e:.0f}, "
          f"ground {lattice_ground_energy(10):.0f}")
    assert e == lattice_ground_energy(10)
    # Sabotage check: flipping one spin must cost exactly 2*degree of satisfied
    # edges.  A corner has 2 neighbours, so the energy rises by 4.
    bad = board.copy(); bad[0] = -bad[0]
    assert ising_energy(J, h, bad) - e == 4.0

    v = np.array([3.0, 1.0, 1.0, 2.0, 2.0, 1.0])
    Jp, hp = number_partition(v)
    s = np.array([1.0, -1.0, -1.0, 1.0, -1.0, -1.0])   # {3,2} vs {1,1,2,1}
    print(f"partition {v.tolist()} -> difference "
          f"{partition_difference(v, s):.0f} (a perfect split is 0)")
    assert partition_difference(v, s) == 0.0
    # The Hamiltonian must be minimal exactly when the difference is.  Summing
    # equation 17 over i<j rather than over all i != j gives
    # H = (delta^2 - sum v_i^2) / 2, so a perfect split lands at -sum(v^2)/2.
    # Checked unnormalised, where that closed form holds as written.
    Jraw, hraw = number_partition(v, normalise=False)
    assert abs(ising_energy(Jraw, hraw, s) + float(v @ v) / 2) < 1e-9
    worse = np.array([1.0, -1.0, -1.0, 1.0, 1.0, -1.0])   # difference 4
    assert partition_difference(v, worse) == 4.0
    assert ising_energy(Jraw, hraw, worse) > ising_energy(Jraw, hraw, s)
    # Normalising divides every coupling by max(v)^2 and so scales the energy,
    # which must not change which configuration wins.
    assert abs(ising_energy(Jp, hp, s) * v.max() ** 2
               - ising_energy(Jraw, hraw, s)) < 1e-9
    assert ising_energy(Jp, hp, worse) > ising_energy(Jp, hp, s)

    Q = hp_qubo("HPHP")
    Ji, hi = qubo_to_ising(Q)
    print(f"HP 'HPHP': {Q.shape[0]} binary variables -> {Ji.shape[0]} spins")
    print("problems self-check passed")
