//! `.pz` — poles and zeros of the small-signal transfer function.
//!
//! Poles are the `s` for which the linearised network has a non-trivial
//! solution with no excitation: `det(G + sC) = 0`.  Zeros are the `s` at which
//! the transfer from the card's input port to its output port vanishes.  Both
//! are generalised eigenvalue problems on matrices [`crate::ac`] already
//! assembles at the operating point.
//!
//! ## Clearing the `1/s`, without a quadratic pencil
//!
//! The AC assembly carries inductors as their admittance `1/(jωL)`, so the
//! system it describes is `G + sC + Λ/s`, not `G + sC`.  Multiplying through by
//! `s` to clear that would give `s²C + sG + Λ` — a *quadratic* eigenvalue
//! problem, and worse, one whose `n` extra roots include a cluster at the
//! origin that has to be told apart from a genuine pole at the origin (an
//! integrator has one).  Sorting real origin poles from bookkeeping ones by
//! magnitude is exactly the sort of plausible-looking rule that produces a
//! confidently wrong answer.
//!
//! So this does not clear the `1/s`.  It reintroduces each inductor current as
//! an unknown instead — the ordinary MNA branch stamp, `v₊ − v₋ = sL·i` — which
//! is linear in `s` by construction.  The pencil stays first-order, every
//! finite eigenvalue is a pole, and the non-dynamic modes come back from the QZ
//! flagged as infinite (`β = 0`) rather than having to be guessed at.  That is
//! what [`crate::ac::LBranch`] exists to carry.
//!
//! ## Zeros
//!
//! For `H(s) = lᵀ(G + sC)⁻¹b`, the Schur complement gives
//!
//! ```text
//!     det [ G + sC   b ]  =  −det(G + sC) · H(s)
//!         [   lᵀ     0 ]
//! ```
//!
//! so the zeros of `H` are the roots of that bordered determinant — the same
//! pencil, bordered by the input column and output row, and solved by the same
//! QZ.  A root that is also a pole is a pole-zero cancellation and is reported
//! by both lists, which is honest: the cancellation is a property of the
//! circuit, not an artefact to be tidied away.
//!
//! ## Why dense, and why it refuses to grow
//!
//! QZ here is dense and `O(N³)`.  That is fine for the size of circuit anyone
//! reads a pole-zero listing for and useless above it, so there is a hard
//! ceiling ([`MAX_PZ_SIZE`]) and an error naming it rather than an hour of
//! silence. A sparse shift-invert Arnoldi pass is the real answer for large
//! circuits and is its own project; a `.pz` that quietly took forty minutes
//! would be a worse outcome than one that says it will not.

use faer::Mat;
use fairchild_parser::{Netlist, PzDrive, PzWant};

use crate::ac::assemble_ac_matrices;
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::mna::CircuitTopology;
use crate::options::SimOptions;

/// Largest pencil `.pz` will factor densely.
///
/// 400 is roughly where dense QZ stops being instant on one core (~0.5 s) and
/// starts being something you wait for; twice that is eight times the work.
pub const MAX_PZ_SIZE: usize = 400;

/// One root of the characteristic equation, in rad/s.
///
/// Complex, and not folded into conjugate pairs — a listing that showed only
/// the upper half-plane would hide whether the solver actually found both.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Root {
    /// Real part, rad/s.  Negative is a stable (decaying) mode.
    pub re: f64,
    /// Imaginary part, rad/s.
    pub im: f64,
}

impl Root {
    /// The root in Hz, as `s/2π` — the axis a Bode plot is read on.
    pub fn hz(&self) -> (f64, f64) {
        let tau = std::f64::consts::TAU;
        (self.re / tau, self.im / tau)
    }
}

/// What `.pz` reports.
#[derive(Debug, Clone)]
pub struct PzResult {
    pub poles: Vec<Root>,
    pub zeros: Vec<Root>,
    /// Modes the pencil reported as infinite — the algebraic (non-dynamic)
    /// part of the system, one per unknown that carries no reactance.  Not an
    /// error, and not hidden either: if this is the whole system, the circuit
    /// has no dynamics and both lists above are legitimately empty.
    pub infinite_poles: usize,
    /// The same count for the bordered pencil the zeros came from.
    pub infinite_zeros: usize,
}

impl PzResult {
    /// Every root, labelled `pole(n)` / `zero(n)` the way ngspice names them.
    fn labelled(&self) -> Vec<(String, Root)> {
        let mut v: Vec<(String, Root)> = self
            .poles
            .iter()
            .enumerate()
            .map(|(i, r)| (format!("pole({})", i + 1), *r))
            .collect();
        v.extend(
            self.zeros
                .iter()
                .enumerate()
                .map(|(i, r)| (format!("zero({})", i + 1), *r)),
        );
        v
    }

    /// Both rad/s and Hz, because a pole is quoted either way depending on who
    /// is reading it, and converting by hand is where a factor of 2π gets lost.
    pub fn write_csv<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        writeln!(w, "root,real_rad_s,imag_rad_s,real_hz,imag_hz")?;
        for (name, r) in self.labelled() {
            let (fre, fim) = r.hz();
            writeln!(
                w,
                "{name},{:.6e},{:.6e},{:.6e},{:.6e}",
                r.re, r.im, fre, fim
            )?;
        }
        writeln!(w, "# infinite (non-dynamic) modes: {}", self.infinite_poles)?;
        // A bordered pencil always carries at least one infinite mode, so a
        // non-zero count is what says the zeros were actually asked for — as
        // opposed to a `pol` run, where an empty zero list means "not computed"
        // rather than "none".
        if self.infinite_zeros > 0 {
            writeln!(
                w,
                "# infinite modes of the bordered pencil: {}",
                self.infinite_zeros
            )?;
        }
        Ok(())
    }

    pub fn write_nutmeg<W: std::io::Write>(&self, mut w: W, title: &str) -> std::io::Result<()> {
        let rows = self.labelled();
        writeln!(w, "Title: {title}")?;
        writeln!(w, "Plotname: Pole-Zero Analysis")?;
        writeln!(w, "Flags: complex")?;
        writeln!(w, "No. Variables: {}", rows.len())?;
        writeln!(w, "No. Points: 1")?;
        writeln!(w, "Variables:")?;
        for (i, (name, _)) in rows.iter().enumerate() {
            writeln!(w, "\t{i}\t{name}\tnotype")?;
        }
        writeln!(w, "Values:")?;
        for (i, (_, r)) in rows.iter().enumerate() {
            if i == 0 {
                writeln!(w, " 0\t{:.6e},{:.6e}", r.re, r.im)?;
            } else {
                writeln!(w, "\t{:.6e},{:.6e}", r.re, r.im)?;
            }
        }
        Ok(())
    }
}

/// `.pz` at the DC operating point.
pub fn pole_zero(
    netlist: &Netlist,
    registry: &DeviceRegistry,
    opts: &SimOptions,
    in_pos: &str,
    in_neg: &str,
    out_pos: &str,
    out_neg: &str,
    drive: PzDrive,
    want: PzWant,
) -> Result<PzResult, SimError> {
    let (topo, g_mat, c_mat, _l_mat, l_branches) = assemble_ac_matrices(netlist, registry, opts)?;
    let n = topo.size;

    // Rows: the MNA unknowns, one per inductor current, and — for a voltage
    // drive — one for the driving source's own branch.  The source is part of
    // the network the poles belong to: driving the input port from a voltage
    // source shorts it for the homogeneous problem, and a circuit's poles
    // measured through a short are not the ones measured through an open.
    let n_l = l_branches.len();

    // A `vol` drive needs a voltage source across the input port.  Nearly every
    // deck carrying a `.pz` card already has one there — that is what makes the
    // port interesting — and adding a second in parallel with it would leave
    // the pencil rank-deficient, which QZ reports as garbage rather than as an
    // error. So the deck's own source is used when there is one, and its branch
    // row *is* the excitation; a new row is only added when the port is
    // undriven.
    let existing = match drive {
        PzDrive::Vol => find_vsrc_across(netlist, &topo, in_pos, in_neg),
        PzDrive::Cur => None,
    };
    let vol_row = matches!(drive, PzDrive::Vol)
        .then_some(n + n_l)
        .filter(|_| existing.is_none());
    let size = n + n_l + usize::from(vol_row.is_some());

    if size > MAX_PZ_SIZE {
        return Err(SimError::ParameterError(format!(
            "`.pz` builds a dense {size}×{size} pencil, over the {MAX_PZ_SIZE} limit. \
             Dense QZ is O(N³) and this would run for a long time without saying so. \
             Reduce the circuit, or extract the sub-block you want the poles of; a \
             sparse eigensolver for large circuits is not implemented"
        )));
    }

    let node = |name: &str| resolve_node(&topo, name);
    let (ip, ineg) = (node(in_pos)?, node(in_neg)?);
    let (op, oneg) = (node(out_pos)?, node(out_neg)?);

    let (mut g, mut c) = (Mat::<f64>::zeros(size, size), Mat::<f64>::zeros(size, size));
    for (i, row) in g_mat.iter().enumerate() {
        for (j, v) in row.iter() {
            g[(i, j)] += v;
        }
    }
    for (i, row) in c_mat.iter().enumerate() {
        for (j, v) in row.iter() {
            c[(i, j)] += v;
        }
    }

    // Each inductor as a branch: `v₊ − v₋ − sL·i = 0`, with `+i` leaving the
    // `+` node.  Identical to the `Λ/s` admittance stamp it replaces — eliminate
    // `i` from these two and `(1/sL)(v₊ − v₋)` is what lands back on the node
    // rows — but linear in `s`.
    for (b, br) in l_branches.iter().enumerate() {
        let r = n + b;
        for (nd, sign) in [(br.pos, 1.0), (br.neg, -1.0)] {
            if let Some(k) = nd {
                g[(k, r)] += sign;
                g[(r, k)] += sign;
            }
        }
        c[(r, r)] -= br.henries;
    }

    // The excitation column and the observation row.
    let mut b_col = vec![0.0; size];
    let mut l_row = vec![0.0; size];
    if let Some((row, sign)) = existing {
        b_col[row] = sign;
    } else if let Some(r) = vol_row {
        for (nd, sign) in [(ip, 1.0), (ineg, -1.0)] {
            if let Some(k) = nd {
                g[(k, r)] += sign;
                g[(r, k)] += sign;
            }
        }
        b_col[r] = 1.0;
    } else {
        // A current injected into `in_pos` and out of `in_neg`.  The overall
        // sign of `b` scales `H(s)` and cannot move a root, so the injection
        // convention is free here in a way it is not for `.tf`.
        if let Some(k) = ip {
            b_col[k] += 1.0;
        }
        if let Some(k) = ineg {
            b_col[k] -= 1.0;
        }
    }
    if let Some(k) = op {
        l_row[k] += 1.0;
    }
    if let Some(k) = oneg {
        l_row[k] -= 1.0;
    }
    if l_row.iter().all(|v| *v == 0.0) {
        return Err(SimError::ParameterError(format!(
            "`.pz` output port ({out_pos}, {out_neg}) is ground-to-ground, so the \
             transfer function is identically zero and has no zeros to report"
        )));
    }

    let (poles, infinite_poles) = if matches!(want, PzWant::Zeros) {
        (Vec::new(), 0)
    } else {
        roots(&g, &c)?
    };

    let (zeros, infinite_zeros) = if matches!(want, PzWant::Poles) {
        (Vec::new(), 0)
    } else {
        let mut gz = Mat::<f64>::zeros(size + 1, size + 1);
        let mut cz = Mat::<f64>::zeros(size + 1, size + 1);
        for i in 0..size {
            for j in 0..size {
                gz[(i, j)] = g[(i, j)];
                cz[(i, j)] = c[(i, j)];
            }
            gz[(i, size)] = b_col[i];
            gz[(size, i)] = l_row[i];
        }
        roots(&gz, &cz)?
    };

    Ok(PzResult {
        poles,
        zeros,
        infinite_poles,
        infinite_zeros,
    })
}

/// Finite roots of `det(G + sC) = 0`, and how many modes came back infinite.
///
/// `(G + sC)v = 0` is `(−G)v = s·Cv`, a generalised eigenproblem whose QZ
/// reports each eigenvalue as a ratio `α/β` rather than a number — precisely so
/// that `β = 0`, the infinite eigenvalue, is representable instead of being an
/// overflow.  `C` is singular for any real circuit (every node without a
/// capacitor to somewhere is an algebraic constraint), so that case is the
/// normal one, not an edge case.
fn roots(g: &Mat<f64>, c: &Mat<f64>) -> Result<(Vec<Root>, usize), SimError> {
    let neg_g = -g;
    let eig = faer::linalg::solvers::GeneralizedEigen::new_from_real(neg_g.as_ref(), c.as_ref())
        .map_err(|_| {
            SimError::ParameterError(
                "`.pz` eigensolve did not converge. The pencil is usually ill-conditioned \
                 because the operating point is: check that .op converges cleanly first"
                    .into(),
            )
        })?;

    let (alpha, beta) = (eig.S_a().column_vector(), eig.S_b().column_vector());
    let mut finite = Vec::new();
    let mut infinite = 0;
    for i in 0..alpha.nrows() {
        let (a, b) = (alpha[i], beta[i]);
        let (an, bn) = (a.norm(), b.norm());
        // LAPACK's test: a mode is infinite when its β vanishes *relative to
        // the α it is paired with*.  Comparing β against a fixed floor instead
        // would call a slow pole infinite on a circuit scaled in seconds and
        // finite on the same circuit scaled in nanoseconds.
        if bn <= f64::EPSILON * an || bn == 0.0 {
            infinite += 1;
            continue;
        }
        let s = a / b;
        if !s.re.is_finite() || !s.im.is_finite() {
            infinite += 1;
            continue;
        }
        finite.push(Root { re: s.re, im: s.im });
    }
    // Slowest first — the pole that sets the bandwidth is the one being looked
    // for, and it is the one nearest the origin.
    finite.sort_by(|p, q| {
        (p.re * p.re + p.im * p.im)
            .partial_cmp(&(q.re * q.re + q.im * q.im))
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(p.im.partial_cmp(&q.im).unwrap_or(std::cmp::Ordering::Equal))
    });
    Ok((finite, infinite))
}

/// A node name as a matrix row, with ground reported as `None` rather than as
/// row zero — ground has no row, and silently using one would pin the wrong
/// potential.
fn resolve_node(topo: &CircuitTopology, name: &str) -> Result<Option<usize>, SimError> {
    if name == "0" || name == "gnd" {
        return Ok(None);
    }
    topo.node_index
        .get(name)
        .copied()
        .map(Some)
        .ok_or_else(|| SimError::UnknownNode(name.to_string()))
}

/// An independent voltage source connected directly across `(pos, neg)`, as
/// `(its MNA branch row, +1 if it faces that way and −1 if reversed)`.
fn find_vsrc_across(
    netlist: &Netlist,
    topo: &CircuitTopology,
    pos: &str,
    neg: &str,
) -> Option<(usize, f64)> {
    netlist.elements.iter().find_map(|el| {
        let fairchild_parser::Element::VoltageSource {
            name,
            pos: p,
            neg: n,
            ..
        } = el
        else {
            return None;
        };
        let sign = if p == pos && n == neg {
            1.0
        } else if p == neg && n == pos {
            -1.0
        } else {
            return None;
        };
        let idx = *topo.vsrc_index.get(&name.to_lowercase())?;
        Some((topo.n_nodes() + idx, sign))
    })
}
