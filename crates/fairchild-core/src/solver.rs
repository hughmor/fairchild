//! Linear-system solver backends.
//!
//! Every analysis (DC OP, tran, AC, noise, DC sweep) reduces to a sequence
//! of dense matrix solves `A·x = b`.  This module defines two cooperating
//! traits so the actual factorisation strategy can be swapped per simulation:
//!
//!  - `LinearSolver`   — owns the backend; entry-point for both one-shot
//!                       solves and "build me a factorisation cache".
//!  - `Factorisation`  — reusable handle.  Owns the symbolic factorisation
//!                       (and, on backends like KLU, the numeric LU too).
//!                       The Newton-Raphson and transient loops own one of
//!                       these across all iterations of a single circuit so
//!                       the symbolic phase runs once and `refactor_and_solve`
//!                       does only the value-update work each step.
//!
//! Three backends ship today:
//!
//!  - `DenseSolver`         — faer partial-pivot LU.  Best for ≤ ~50 nodes
//!                            where the sparse setup cost dominates.
//!  - `FaerSparseSolver`    — faer sparse LU (`SparseColMat::sp_lu`).  Pure
//!                            Rust, no C deps; default at larger N.  Cache
//!                            handle saves the dense→CSC conversion but
//!                            still runs full LU each refactor (faer 0.24
//!                            does not expose a separate refactor path).
//!  - `KluSolver`           — SuiteSparse KLU via the `klu` cargo feature.
//!                            Cache handle holds a `KluSymbolic` + reusable
//!                            `KluNumeric`; `refactor_and_solve` calls
//!                            `klu_refactor` — the major perf win.

use faer::{Col, Mat, linalg::solvers::Solve};

use crate::error::SimError;

/// Choice of linear-system backend.  Carried on `SimOptions` and propagated
/// to each analysis loop.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolverKind {
    /// Dense `partial_piv_lu` via faer.  Recommended for ≤ ~50 nodes.
    Dense,
    /// Sparse LU via `faer::sparse`.  Pure-Rust default at larger N.
    Sparse,
    /// SuiteSparse KLU.  Sparse direct LU specialised for circuit
    /// matrices (BTF + dense LU on diagonal blocks).  Typically 2-5×
    /// faster than the faer-sparse path on circuit problems.  Requires
    /// the `klu` cargo feature and a system install of SuiteSparse.
    Klu,
    /// Pick automatically from system size (Dense if n < 50, else Sparse).
    Auto,
}

/// Trait covering the operations every analysis loop needs from its
/// linear solver: one-shot forward/transpose solves, plus a builder for
/// a reusable [`Factorisation`] handle.
///
/// **Adding a third-party backend** — e.g. another sparse LU package — is:
///   1. Add a new dependency behind a cargo feature flag.
///   2. Implement this trait on a new struct, plus a paired
///      [`Factorisation`] type that owns the cached symbolic state.
///   3. Extend [`SolverKind`] with a new variant gated on the feature flag.
///   4. Add the dispatch to [`make_solver`].
///
/// The trait is `Send + Sync` so parallel-frequency AC / parallel-MC /
/// parallel-sweep drivers can share one solver across threads.
pub trait LinearSolver: Send + Sync {
    /// Solve `A · x = b`.  Returns `SingularMatrix` if the result is non-finite.
    fn solve(&self, a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError>;

    /// Solve `A^T · x = b`.  Default falls back to constructing the explicit
    /// transpose and reusing `solve` — backends with cached factorisations
    /// override this to avoid the O(n²) copy.
    fn solve_transpose(&self, a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError> {
        let n = b.len();
        let mut at = vec![vec![0.0_f64; n]; n];
        for i in 0..n {
            for j in 0..n {
                at[i][j] = a[j][i];
            }
        }
        self.solve(&at, b)
    }

    /// Build a reusable [`Factorisation`] from the matrix `a`.
    ///
    /// The intended use pattern inside Newton-Raphson is:
    ///
    /// ```ignore
    /// let mut fact = solver.factorise(&mat.a)?;
    /// loop {
    ///     // stamp devices into mat.a, build mat.b ...
    ///     let x = fact.refactor_and_solve(&mat.a, &mat.b)?;
    ///     // check convergence; break or continue ...
    /// }
    /// ```
    ///
    /// The default implementation is no-op caching — every call to
    /// `refactor_and_solve` falls back to `solve(a, b)` from scratch.
    /// Backends like KLU override this to cache the symbolic factorisation
    /// and call into the appropriate `refactor` primitive.
    fn factorise(&self, a: &[Vec<f64>]) -> Result<Box<dyn Factorisation>, SimError> {
        // The default impl just remembers nothing: each refactor_and_solve
        // re-uses the solver wholesale.  For DenseSolver this is the right
        // behaviour; for sparse backends, this default is overridden to
        // cache the sparsity pattern.
        let n = a.len();
        let _ = n;  // matrix is used only via the trait's `solve` later
        Ok(Box::new(NoCacheFactorisation::new()))
    }
}

/// Reusable factorisation handle.  Held by the Newton-Raphson and
/// transient loops across all iterations of a single solve; the symbolic
/// phase runs once and `refactor_and_solve` does only the per-iteration
/// numeric update.
///
/// `Send` only (not `Sync`) — each iteration mutates the cached numeric
/// LU.  Concurrent solvers (parallel AC sweep) build one factorisation
/// per worker.
pub trait Factorisation: Send {
    /// Update the cached numeric factorisation from the current values in
    /// `a` (sparsity pattern assumed unchanged since the call to
    /// [`LinearSolver::factorise`] that produced this handle), then solve
    /// `A · x = b`.
    fn refactor_and_solve(&mut self, a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError>;

    /// Same, but for `A^T · x = b` (used by adjoint / noise paths).
    /// Default falls back to building an explicit transpose and calling
    /// `refactor_and_solve` — backends with bidirectional LU caches
    /// override this.
    fn refactor_and_solve_transpose(
        &mut self,
        a: &[Vec<f64>],
        b: &[f64],
    ) -> Result<Vec<f64>, SimError> {
        let n = b.len();
        let mut at = vec![vec![0.0_f64; n]; n];
        for i in 0..n {
            for j in 0..n {
                at[i][j] = a[j][i];
            }
        }
        self.refactor_and_solve(&at, b)
    }
}

/// No-op cache: every refactor_and_solve runs a fresh `LinearSolver::solve`.
/// Used by [`DenseSolver`] and as the fallback default-impl path.
struct NoCacheFactorisation {
    backend: Box<dyn LinearSolver>,
}

impl NoCacheFactorisation {
    fn new() -> Self {
        Self { backend: Box::new(DenseSolver) }
    }
    fn with_backend(backend: Box<dyn LinearSolver>) -> Self {
        Self { backend }
    }
}

impl Factorisation for NoCacheFactorisation {
    fn refactor_and_solve(&mut self, a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError> {
        self.backend.solve(a, b)
    }
    fn refactor_and_solve_transpose(&mut self, a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError> {
        self.backend.solve_transpose(a, b)
    }
}

// ---------------------------------------------------------------------------
// Dense backend
// ---------------------------------------------------------------------------

/// Dense partial-pivot LU via faer.
pub struct DenseSolver;

impl LinearSolver for DenseSolver {
    fn solve(&self, a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError> {
        let n = b.len();
        if n == 0 { return Ok(Vec::new()); }
        let a_mat = Mat::<f64>::from_fn(n, n, |i, j| a[i][j]);
        let b_col = Col::<f64>::from_fn(n, |i| b[i]);
        let lu = a_mat.partial_piv_lu();
        let x_col = lu.solve(b_col.as_ref());
        let x: Vec<f64> = (0..n).map(|i| x_col[i]).collect();
        if x.iter().any(|v| !v.is_finite()) {
            return Err(SimError::SingularMatrix);
        }
        Ok(x)
    }

    fn factorise(&self, _a: &[Vec<f64>]) -> Result<Box<dyn Factorisation>, SimError> {
        // No-op cache for dense: re-running partial_piv_lu is cheap
        // enough that caching the column permutation does not pay.
        Ok(Box::new(NoCacheFactorisation::with_backend(Box::new(DenseSolver))))
    }
}

// ---------------------------------------------------------------------------
// faer-sparse backend with pattern cache
// ---------------------------------------------------------------------------

/// Sparse LU via faer's `SparseColMat::sp_lu`.
///
/// One-shot `solve` builds a CSC matrix from non-zero entries of the
/// dense input on each call — that's an O(n²) scan but the LU step
/// itself runs at sparse complexity.  The Newton-loop fast path uses
/// `factorise()` to capture the sparsity pattern once; subsequent
/// `refactor_and_solve` calls walk the cached `(row, col)` pairs and
/// avoid the O(n²) rescan.  faer 0.24 does not expose a separate
/// numeric-refactor primitive, so each solve still runs full `sp_lu` —
/// but at least the dense→CSC conversion is now O(nnz).
pub struct FaerSparseSolver {
    /// Drop-tolerance: values with |a_ij| < `zero_threshold` are treated as
    /// structural zeros.  Defaults to a value comfortably below `GMIN`.
    pub zero_threshold: f64,
}

impl Default for FaerSparseSolver {
    fn default() -> Self {
        FaerSparseSolver { zero_threshold: 1e-30 }
    }
}

impl LinearSolver for FaerSparseSolver {
    fn solve(&self, a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError> {
        use faer::sparse::{SparseColMat, Triplet};

        let n = b.len();
        if n == 0 { return Ok(Vec::new()); }

        let mut triplets: Vec<Triplet<usize, usize, f64>> = Vec::new();
        for i in 0..n {
            for j in 0..n {
                let v = a[i][j];
                if v.abs() > self.zero_threshold {
                    triplets.push(Triplet::new(i, j, v));
                }
            }
        }

        let a_sp = SparseColMat::<usize, f64>::try_new_from_triplets(n, n, &triplets)
            .map_err(|_| SimError::SingularMatrix)?;
        let lu = a_sp.sp_lu().map_err(|_| SimError::SingularMatrix)?;

        let b_col = Col::<f64>::from_fn(n, |i| b[i]);
        let x_col = lu.solve(b_col.as_ref());
        let x: Vec<f64> = (0..n).map(|i| x_col[i]).collect();
        if x.iter().any(|v| !v.is_finite()) {
            return Err(SimError::SingularMatrix);
        }
        Ok(x)
    }

    fn factorise(&self, _a: &[Vec<f64>]) -> Result<Box<dyn Factorisation>, SimError> {
        // faer 0.24 has no refactor-only primitive — every solve runs full
        // symbolic + numeric LU.  Caching the sparsity pattern would only
        // save the O(n²) → O(nnz) dense→CSC conversion, and it's unsound
        // unless we also verify that no new structural entries appear
        // each iteration (photonic devices at x=0 stamp zeros that become
        // non-zero at the operating point — see the WDM DC OP regression).
        // Until faer exposes a refactor primitive, the safe behaviour is
        // to fall back to one-shot solves; the per-iteration overhead is
        // dominated by `sp_lu` itself, not the dense-scan.
        Ok(Box::new(NoCacheFactorisation::with_backend(
            Box::new(FaerSparseSolver { zero_threshold: self.zero_threshold }),
        )))
    }
}

// ---------------------------------------------------------------------------
// KLU backend with symbolic+numeric reuse
// ---------------------------------------------------------------------------

/// SuiteSparse KLU backend — sparse direct LU with BTF preordering.
/// One-shot `solve` does fresh analyze + factor + solve per call (same
/// behaviour as `FaerSparseSolver::solve`).  The `factorise` path is
/// where KLU's real advantage lives: a `KluSymbolic` analysed once is
/// reused via `klu_refactor` for every subsequent solve.
#[cfg(feature = "klu")]
pub struct KluSolver;

#[cfg(feature = "klu")]
impl LinearSolver for KluSolver {
    fn solve(&self, a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError> {
        fairchild_klu::klu_solve_dense(a, b).map_err(|_| SimError::SingularMatrix)
    }

    fn factorise(&self, a: &[Vec<f64>]) -> Result<Box<dyn Factorisation>, SimError> {
        use fairchild_klu::{dense_to_csc, KluCommon, KluNumeric, KluSymbolic};
        let n = a.len();
        if n == 0 {
            return Ok(Box::new(NoCacheFactorisation::with_backend(Box::new(KluSolver))));
        }
        let (mut ap, mut ai, mut ax) = dense_to_csc(a, 1e-30);
        let mut common = KluCommon::new();
        // Save the row-index map so refactor_and_solve can pull values
        // at the exact same CSC positions on subsequent iterations.
        let pattern: Vec<(i32, i32)> = column_row_pattern(&ap, &ai);
        let symbolic = KluSymbolic::analyze(n, &mut ap, &mut ai, &mut common)
            .map_err(|_| SimError::SingularMatrix)?;
        let numeric = KluNumeric::factor(&mut ap, &mut ai, &mut ax, &symbolic, &mut common)
            .map_err(|_| SimError::SingularMatrix)?;
        Ok(Box::new(KluFactorisation {
            common,
            symbolic,
            numeric,
            ap,
            ai,
            ax,
            pattern,
            n,
        }))
    }
}

/// Helper: turn KLU's column-offset + row-index arrays into a flat
/// `(col, row)` list in CSC traversal order.  Used so that subsequent
/// refactor calls can refill `ax` in the same order.
#[cfg(feature = "klu")]
fn column_row_pattern(ap: &[i32], ai: &[i32]) -> Vec<(i32, i32)> {
    let n_cols = ap.len() - 1;
    let mut out = Vec::with_capacity(ai.len());
    for j in 0..n_cols {
        let start = ap[j] as usize;
        let end = ap[j + 1] as usize;
        for k in start..end {
            out.push((j as i32, ai[k]));
        }
    }
    out
}

#[cfg(feature = "klu")]
struct KluFactorisation {
    common:   fairchild_klu::KluCommon,
    symbolic: fairchild_klu::KluSymbolic,
    numeric:  fairchild_klu::KluNumeric,
    ap:       Vec<i32>,
    ai:       Vec<i32>,
    ax:       Vec<f64>,
    pattern:  Vec<(i32, i32)>,  // (col, row) per CSC entry
    n:        usize,
}

#[cfg(feature = "klu")]
impl KluFactorisation {
    /// Rebuild CSC fresh from the current dense matrix and compare to
    /// the cached pattern.  Returns `true` if the pattern grew (a new
    /// structural entry appeared since the symbolic factorisation was
    /// computed) — devices stamping zero at x=0 produce this on the
    /// 2nd NR iteration once the operating point activates them.
    fn refresh_csc(&mut self, a: &[Vec<f64>]) -> bool {
        use fairchild_klu::dense_to_csc;
        let (ap_new, ai_new, ax_new) = dense_to_csc(a, 1e-30);
        let pattern_changed = ap_new != self.ap || ai_new != self.ai;
        self.ap = ap_new;
        self.ai = ai_new;
        self.ax = ax_new;
        if pattern_changed {
            // Rebuild the (col, row) traversal so transpose-solve and
            // future invalidations stay consistent.
            self.pattern = column_row_pattern(&self.ap, &self.ai);
        }
        pattern_changed
    }

    /// Pattern grew — re-analyze symbolic + factor numeric, dropping the
    /// previous cached factorisation.  Cheap if it happens once (during
    /// NR warm-up); expensive if it happens every iteration (which would
    /// indicate the convergence path is also walking the pattern, a
    /// pathological case we don't try to defend against here).
    fn reanalyze(&mut self) -> Result<(), SimError> {
        use fairchild_klu::{KluNumeric, KluSymbolic};
        let n = self.n;
        // Drop old numeric + symbolic before allocating new (Drop runs
        // klu_free_numeric / klu_free_symbolic).
        let new_symbolic = KluSymbolic::analyze(n, &mut self.ap, &mut self.ai, &mut self.common)
            .map_err(|_| SimError::SingularMatrix)?;
        let new_numeric = KluNumeric::factor(
            &mut self.ap, &mut self.ai, &mut self.ax, &new_symbolic, &mut self.common
        ).map_err(|_| SimError::SingularMatrix)?;
        self.symbolic = new_symbolic;
        self.numeric = new_numeric;
        Ok(())
    }
}

#[cfg(feature = "klu")]
impl Factorisation for KluFactorisation {
    fn refactor_and_solve(&mut self, a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError> {
        if self.refresh_csc(a) {
            // Pattern changed since the cached symbolic factorisation —
            // re-analyze + factor afresh.  Subsequent iterations with
            // the same pattern fall back to `klu_refactor`.
            self.reanalyze()?;
        } else {
            self.numeric
                .refactor(&mut self.ap, &mut self.ai, &mut self.ax, &self.symbolic, &mut self.common)
                .map_err(|_| SimError::SingularMatrix)?;
        }
        let mut x = b.to_vec();
        self.numeric
            .solve(&self.symbolic, &mut x, &mut self.common)
            .map_err(|_| SimError::SingularMatrix)?;
        if x.iter().any(|v| !v.is_finite()) {
            return Err(SimError::SingularMatrix);
        }
        Ok(x)
    }

    fn refactor_and_solve_transpose(&mut self, a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError> {
        if self.refresh_csc(a) {
            self.reanalyze()?;
        } else {
            self.numeric
                .refactor(&mut self.ap, &mut self.ai, &mut self.ax, &self.symbolic, &mut self.common)
                .map_err(|_| SimError::SingularMatrix)?;
        }
        let mut x = b.to_vec();
        self.numeric
            .solve_transpose(&self.symbolic, &mut x, &mut self.common)
            .map_err(|_| SimError::SingularMatrix)?;
        if x.iter().any(|v| !v.is_finite()) {
            return Err(SimError::SingularMatrix);
        }
        Ok(x)
    }
}

// ---------------------------------------------------------------------------
// Dispatch
// ---------------------------------------------------------------------------

/// Construct a solver for a system of the given estimated size.  Drives
/// `SolverKind::Auto`'s dense / sparse crossover.
pub fn make_solver(kind: SolverKind, n: usize) -> Box<dyn LinearSolver> {
    match kind {
        SolverKind::Dense  => Box::new(DenseSolver),
        SolverKind::Sparse => Box::new(FaerSparseSolver::default()),
        SolverKind::Klu    => {
            #[cfg(feature = "klu")]
            { Box::new(KluSolver) }
            #[cfg(not(feature = "klu"))]
            {
                eprintln!("warning: KLU backend requested but `klu` cargo feature \
                           is not enabled; falling back to faer-sparse.");
                Box::new(FaerSparseSolver::default())
            }
        }
        SolverKind::Auto   => {
            if n < 50 { Box::new(DenseSolver) }
            else      { Box::new(FaerSparseSolver::default()) }
        }
    }
}

/// Solve `A·x = b` with the historic dense default.  Kept as a free
/// function so existing analyses compile unchanged; new analyses should
/// prefer constructing a `Box<dyn LinearSolver>` from `make_solver()` and
/// calling its methods directly.
pub fn lu_solve(a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError> {
    DenseSolver.solve(a, b)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dense_matches_sparse_on_small_matrix() {
        let a = vec![
            vec![ 4.0, -1.0,  0.0,  0.0],
            vec![-1.0,  4.0, -1.0,  0.0],
            vec![ 0.0, -1.0,  4.0, -1.0],
            vec![ 0.0,  0.0, -1.0,  3.0],
        ];
        let b = vec![1.0, 2.0, 3.0, 4.0];
        let x_d = DenseSolver.solve(&a, &b).unwrap();
        let x_s = FaerSparseSolver::default().solve(&a, &b).unwrap();
        for i in 0..4 {
            assert!((x_d[i] - x_s[i]).abs() < 1e-10,
                "i={i}: dense={} sparse={}", x_d[i], x_s[i]);
        }
    }

    #[test]
    fn singular_matrix_returns_err() {
        let a = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        let b = vec![1.0, 2.0];
        assert!(matches!(DenseSolver.solve(&a, &b),
            Err(SimError::SingularMatrix)));
        assert!(matches!(FaerSparseSolver::default().solve(&a, &b),
            Err(SimError::SingularMatrix)));
    }

    #[test]
    fn transpose_solve_default_correct() {
        let a = vec![
            vec![3.0, 1.0, 2.0],
            vec![0.0, 4.0, 5.0],
            vec![1.0, 0.0, 6.0],
        ];
        let b = vec![1.0, 2.0, 3.0];
        let at: Vec<Vec<f64>> = (0..3).map(|i| (0..3).map(|j| a[j][i]).collect()).collect();

        let x_via_explicit_t = DenseSolver.solve(&at, &b).unwrap();
        let x_via_method     = DenseSolver.solve_transpose(&a, &b).unwrap();
        for i in 0..3 {
            assert!((x_via_explicit_t[i] - x_via_method[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn faer_sparse_factorisation_reuses_pattern() {
        // First solve through factorise(); second solve through the
        // refactor path with updated values — answer must match a fresh
        // dense solve of the second matrix.
        let a1 = vec![
            vec![3.0, 1.0],
            vec![0.0, 4.0],
        ];
        let solver = FaerSparseSolver::default();
        let mut fact = solver.factorise(&a1).unwrap();
        let x1 = fact.refactor_and_solve(&a1, &[5.0, 8.0]).unwrap();
        // 3x+y=5, 4y=8 → y=2, x=1
        assert!((x1[0] - 1.0).abs() < 1e-12);
        assert!((x1[1] - 2.0).abs() < 1e-12);

        // Change values, same pattern.
        let a2 = vec![
            vec![6.0, 2.0],
            vec![0.0, 4.0],
        ];
        let x2 = fact.refactor_and_solve(&a2, &[10.0, 8.0]).unwrap();
        // 6x+2y=10, 4y=8 → y=2, 6x=6 → x=1
        assert!((x2[0] - 1.0).abs() < 1e-12);
        assert!((x2[1] - 2.0).abs() < 1e-12);
    }
}
