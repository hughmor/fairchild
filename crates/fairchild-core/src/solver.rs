//! Linear-system solver backends.
//!
//! Every analysis (DC OP, tran, AC, noise, DC sweep) reduces to a sequence
//! of dense matrix solves `A·x = b`.  This module defines two cooperating
//! traits so the actual factorisation strategy can be swapped per simulation:
//!
//!  - `LinearSolver`   — owns the backend; entry-point for both one-shot
//!    solves and "build me a factorisation cache".
//!  - `Factorisation`  — reusable handle.  Owns the symbolic factorisation
//!    (and, on backends like KLU, the numeric LU too).
//!    The Newton-Raphson and transient loops own one of
//!    these across all iterations of a single circuit so
//!    the symbolic phase runs once and `refactor_and_solve`
//!    does only the value-update work each step.
//!
//! Three backends ship today:
//!
//!  - `DenseSolver`         — faer partial-pivot LU.  Best for ≤ ~50 nodes
//!    where the sparse setup cost dominates.
//!  - `FaerSparseSolver`    — faer sparse LU.  Pure Rust, no C deps; default
//!    at larger N.  Given a structural pattern (see
//!    `mna::Pattern`) the cache handle keeps both the
//!    CSC structure and the symbolic LU, so a refactor
//!    is a value refill plus
//!    `Lu::try_new_with_symbolic` — no dense rescan,
//!    no re-run of the column ordering.
//!  - `KluSolver`           — SuiteSparse KLU via the `klu` cargo feature.
//!    Cache handle holds a `KluSymbolic` + reusable
//!    `KluNumeric`; `refactor_and_solve` calls
//!    `klu_refactor` — the major perf win.

use faer::{linalg::solvers::Solve, Col, Mat};
use std::sync::Arc;

use crate::error::SimError;
use crate::mna::MnaMatrix;

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
        let _ = n; // matrix is used only via the trait's `solve` later
        Ok(Box::new(NoCacheFactorisation::new()))
    }

    /// [`LinearSolver::factorise`] with access to the whole matrix, including
    /// its structural sparsity pattern when the caller built one.  Sparse
    /// backends override this to cache the CSC structure and the symbolic LU;
    /// the default ignores the pattern.
    fn factorise_mat(&self, mat: &MnaMatrix) -> Result<Box<dyn Factorisation>, SimError> {
        self.factorise(&mat.a)
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

    /// [`Factorisation::refactor_and_solve`] taking the whole matrix, so a
    /// backend holding a cached CSC structure can refill values by walking the
    /// pattern instead of rescanning n² dense cells.  Default delegates.
    fn refactor_and_solve_mat(&mut self, mat: &MnaMatrix) -> Result<Vec<f64>, SimError> {
        self.refactor_and_solve(&mat.a, &mat.b)
    }

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
        Self {
            backend: Box::new(DenseSolver),
        }
    }
    fn with_backend(backend: Box<dyn LinearSolver>) -> Self {
        Self { backend }
    }
}

impl Factorisation for NoCacheFactorisation {
    fn refactor_and_solve(&mut self, a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError> {
        self.backend.solve(a, b)
    }
    fn refactor_and_solve_transpose(
        &mut self,
        a: &[Vec<f64>],
        b: &[f64],
    ) -> Result<Vec<f64>, SimError> {
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
        if n == 0 {
            return Ok(Vec::new());
        }
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
        Ok(Box::new(NoCacheFactorisation::with_backend(Box::new(
            DenseSolver,
        ))))
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
        FaerSparseSolver {
            zero_threshold: 1e-30,
        }
    }
}

impl LinearSolver for FaerSparseSolver {
    fn solve(&self, a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError> {
        use faer::sparse::{SparseColMat, Triplet};

        let n = b.len();
        if n == 0 {
            return Ok(Vec::new());
        }

        let mut triplets: Vec<Triplet<usize, usize, f64>> = Vec::new();
        for (i, row) in a.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
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
        // Without a structural pattern there is nothing sound to cache:
        // discovering the pattern from values on iteration 0 misses cells that
        // are zero at x=0 and non-zero at the solution, which photonic devices
        // do stamp (see the WDM DC OP regression).  Callers that want the fast
        // path hand over a pattern via `factorise_mat`.
        Ok(Box::new(NoCacheFactorisation::with_backend(Box::new(
            FaerSparseSolver {
                zero_threshold: self.zero_threshold,
            },
        ))))
    }

    fn factorise_mat(&self, mat: &MnaMatrix) -> Result<Box<dyn Factorisation>, SimError> {
        match mat.pattern() {
            Some(p) => Ok(Box::new(FaerSparseFactorisation::new(
                Arc::clone(p),
                self.zero_threshold,
            ))),
            None => self.factorise(&mat.a),
        }
    }
}

/// faer-sparse with a cached CSC structure and a cached symbolic LU.
///
/// Two costs disappear versus one-shot `solve`: the O(n²) dense scan that built
/// triplets every iteration, and the column ordering (colamd), which was re-run
/// every iteration even though the pattern never changes.
///
/// The structure is the *numerically* non-zero subset of the structural
/// pattern, not the pattern itself — a structural superset would hand extra
/// explicit zeros to the ordering and pay for the fill-in.  Narrowing is safe
/// because growth is detectable in O(nnz): the refill walk visits every
/// structural cell, so a cell that turns non-zero later is seen on the
/// iteration it happens and triggers one rebuild.
struct FaerSparseFactorisation {
    pattern: Arc<crate::mna::Pattern>,
    zero_threshold: f64,
    /// CSC column pointers / row indices of the active set.
    col_ptr: Vec<usize>,
    row_idx: Vec<usize>,
    /// Where in `values` each active `(row, col)` lives, indexed the same way
    /// the refill walk visits them: row-major over `pattern.cols`.
    slot: Vec<u32>,
    values: Vec<f64>,
    symbolic: Option<faer::sparse::linalg::solvers::SymbolicLu<usize>>,
}

/// `slot` entry for a structural cell that is not in the active set.
const NO_SLOT: u32 = u32::MAX;

impl FaerSparseFactorisation {
    fn new(pattern: Arc<crate::mna::Pattern>, zero_threshold: f64) -> Self {
        FaerSparseFactorisation {
            pattern,
            zero_threshold,
            col_ptr: Vec::new(),
            row_idx: Vec::new(),
            slot: Vec::new(),
            values: Vec::new(),
            symbolic: None,
        }
    }

    /// Rebuild the CSC structure and symbolic LU from the cells of `a` that are
    /// currently non-zero.  Runs on the first solve and again only if a new
    /// cell inside the structural pattern turns non-zero.
    fn rebuild(&mut self, a: &[Vec<f64>]) -> Result<(), SimError> {
        use faer::sparse::linalg::solvers::SymbolicLu;
        use faer::sparse::SymbolicSparseColMatRef;

        let n = a.len();
        let thr = self.zero_threshold;
        // Count per column first so the CSC arrays can be filled in one pass.
        let mut col_count = vec![0usize; n];
        for (i, cols) in self.pattern.cols.iter().enumerate() {
            for &j in cols {
                if a[i][j as usize].abs() > thr {
                    col_count[j as usize] += 1;
                }
            }
        }
        let mut col_ptr = vec![0usize; n + 1];
        for j in 0..n {
            col_ptr[j + 1] = col_ptr[j] + col_count[j];
        }
        let nnz = col_ptr[n];
        let mut row_idx = vec![0usize; nnz];
        let mut values = vec![0.0_f64; nnz];
        let mut fill = col_ptr.clone();
        // Row-major walk keeps `row_idx` ascending within each column, which is
        // what `new_checked` requires of a sorted CSC.
        let mut slot: Vec<u32> = Vec::with_capacity(self.pattern.nnz);
        for (i, cols) in self.pattern.cols.iter().enumerate() {
            for &j in cols {
                let v = a[i][j as usize];
                if v.abs() > thr {
                    let k = fill[j as usize];
                    row_idx[k] = i;
                    values[k] = v;
                    fill[j as usize] += 1;
                    slot.push(k as u32);
                } else {
                    slot.push(NO_SLOT);
                }
            }
        }
        let sym = SymbolicSparseColMatRef::new_checked(n, n, &col_ptr, None, &row_idx);
        self.symbolic = Some(SymbolicLu::try_new(sym).map_err(|_| SimError::SingularMatrix)?);
        self.col_ptr = col_ptr;
        self.row_idx = row_idx;
        self.slot = slot;
        self.values = values;
        Ok(())
    }

    /// Copy current values into the cached CSC.  Returns true if a structural
    /// cell outside the active set has become non-zero, meaning the caller must
    /// rebuild before factorising.
    fn refill(&mut self, a: &[Vec<f64>]) -> bool {
        let thr = self.zero_threshold;
        let mut grew = false;
        let mut k = 0usize;
        for (i, cols) in self.pattern.cols.iter().enumerate() {
            let row = &a[i];
            for &j in cols {
                let v = row[j as usize];
                let s = self.slot[k];
                if s == NO_SLOT {
                    grew |= v.abs() > thr;
                } else {
                    self.values[s as usize] = v;
                }
                k += 1;
            }
        }
        grew
    }

    fn solve_cached(&self, b: &[f64]) -> Result<Vec<f64>, SimError> {
        use faer::sparse::linalg::solvers::Lu;
        use faer::sparse::{SparseColMatRef, SymbolicSparseColMatRef};

        let n = b.len();
        let symbolic = self.symbolic.as_ref().ok_or(SimError::SingularMatrix)?;
        let sym = SymbolicSparseColMatRef::new_checked(n, n, &self.col_ptr, None, &self.row_idx);
        let mat = SparseColMatRef::<usize, f64>::new(sym, &self.values);
        let lu = Lu::try_new_with_symbolic(symbolic.clone(), mat)
            .map_err(|_| SimError::SingularMatrix)?;
        let b_col = Col::<f64>::from_fn(n, |i| b[i]);
        let x_col = lu.solve(b_col.as_ref());
        let x: Vec<f64> = (0..n).map(|i| x_col[i]).collect();
        if x.iter().any(|v| !v.is_finite()) {
            return Err(SimError::SingularMatrix);
        }
        Ok(x)
    }
}

impl Factorisation for FaerSparseFactorisation {
    fn refactor_and_solve(&mut self, a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError> {
        // `refill` indexes the cached slot map, so it is only valid once
        // `rebuild` has run; short-circuit ordering matters here.
        if self.symbolic.is_none() || self.refill(a) {
            self.rebuild(a)?;
        }
        self.solve_cached(b)
    }

    fn refactor_and_solve_transpose(
        &mut self,
        a: &[Vec<f64>],
        b: &[f64],
    ) -> Result<Vec<f64>, SimError> {
        // Adjoint paths are cold; the explicit transpose keeps this simple.
        let n = b.len();
        let mut at = vec![vec![0.0_f64; n]; n];
        for i in 0..n {
            for j in 0..n {
                at[i][j] = a[j][i];
            }
        }
        FaerSparseSolver {
            zero_threshold: self.zero_threshold,
        }
        .solve(&at, b)
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
            return Ok(Box::new(NoCacheFactorisation::with_backend(Box::new(
                KluSolver,
            ))));
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
        for &row in &ai[start..end] {
            out.push((j as i32, row));
        }
    }
    out
}

#[cfg(feature = "klu")]
struct KluFactorisation {
    common: fairchild_klu::KluCommon,
    symbolic: fairchild_klu::KluSymbolic,
    numeric: fairchild_klu::KluNumeric,
    ap: Vec<i32>,
    ai: Vec<i32>,
    ax: Vec<f64>,
    pattern: Vec<(i32, i32)>, // (col, row) per CSC entry
    n: usize,
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
            &mut self.ap,
            &mut self.ai,
            &mut self.ax,
            &new_symbolic,
            &mut self.common,
        )
        .map_err(|_| SimError::SingularMatrix)?;
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
                .refactor(
                    &mut self.ap,
                    &mut self.ai,
                    &mut self.ax,
                    &self.symbolic,
                    &mut self.common,
                )
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

    fn refactor_and_solve_transpose(
        &mut self,
        a: &[Vec<f64>],
        b: &[f64],
    ) -> Result<Vec<f64>, SimError> {
        if self.refresh_csc(a) {
            self.reanalyze()?;
        } else {
            self.numeric
                .refactor(
                    &mut self.ap,
                    &mut self.ai,
                    &mut self.ax,
                    &self.symbolic,
                    &mut self.common,
                )
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
        SolverKind::Dense => Box::new(DenseSolver),
        SolverKind::Sparse => Box::new(FaerSparseSolver::default()),
        SolverKind::Klu => {
            #[cfg(feature = "klu")]
            {
                Box::new(KluSolver)
            }
            #[cfg(not(feature = "klu"))]
            unreachable!(
                "SolverKind::Klu reached make_solver without the `klu` feature — \
                          SimOptions::set should have rejected it earlier"
            )
        }
        SolverKind::Auto => {
            if n < 50 {
                Box::new(DenseSolver)
            } else {
                Box::new(FaerSparseSolver::default())
            }
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

// ---------------------------------------------------------------------------
// Matrix equilibration (two-sided Ruiz scaling)
// ---------------------------------------------------------------------------

/// Two-sided ∞-norm (Ruiz) scaling factors `(D_r, D_c)` such that the rows and
/// columns of `diag(D_r)·A·diag(D_c)` have ∞-norm ≈ 1.  A few iterations is
/// enough in practice; empty rows/columns are left unscaled (factor 1).
fn equilibration_factors(a: &[Vec<f64>]) -> (Vec<f64>, Vec<f64>) {
    let n = a.len();
    let mut dr = vec![1.0_f64; n];
    let mut dc = vec![1.0_f64; n];
    for _ in 0..3 {
        for i in 0..n {
            let mut m = 0.0_f64;
            for j in 0..n {
                m = m.max((dr[i] * a[i][j] * dc[j]).abs());
            }
            if m > 0.0 {
                dr[i] /= m.sqrt();
            }
        }
        for j in 0..n {
            let mut m = 0.0_f64;
            for (i, row) in a.iter().enumerate() {
                m = m.max((dr[i] * row[j] * dc[j]).abs());
            }
            if m > 0.0 {
                dc[j] /= m.sqrt();
            }
        }
    }
    (dr, dc)
}

/// `A' = diag(D_r)·A·diag(D_c)`, `b' = D_r·b`.  Returns the scaled system.
fn apply_equilibration(
    a: &[Vec<f64>],
    b: &[f64],
    dr: &[f64],
    dc: &[f64],
) -> (Vec<Vec<f64>>, Vec<f64>) {
    let n = a.len();
    let mut a_s = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            a_s[i][j] = dr[i] * a[i][j] * dc[j];
        }
    }
    let b_s: Vec<f64> = (0..n).map(|i| dr[i] * b[i]).collect();
    (a_s, b_s)
}

/// LinearSolver decorator that equilibrates the system before delegating to an
/// inner backend, then unscales the solution.  Equilibration is applied to the
/// forward `solve` / `refactor_and_solve` path only; the transpose/adjoint path
/// (`.noise`) is delegated unscaled to keep the scaling unambiguous.
pub struct EquilibratedSolver {
    inner: Box<dyn LinearSolver>,
}

impl EquilibratedSolver {
    pub fn new(inner: Box<dyn LinearSolver>) -> Self {
        Self { inner }
    }
}

impl LinearSolver for EquilibratedSolver {
    fn solve(&self, a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError> {
        let (dr, dc) = equilibration_factors(a);
        let (a_s, b_s) = apply_equilibration(a, b, &dr, &dc);
        let x_s = self.inner.solve(&a_s, &b_s)?;
        Ok((0..x_s.len()).map(|j| dc[j] * x_s[j]).collect())
    }

    fn solve_transpose(&self, a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError> {
        // Adjoint path left unscaled (forward-only equilibration).
        self.inner.solve_transpose(a, b)
    }

    fn factorise(&self, a: &[Vec<f64>]) -> Result<Box<dyn Factorisation>, SimError> {
        // Scaling preserves sparsity, so the inner symbolic factorisation built
        // on the unscaled `a` stays valid across refactors.
        Ok(Box::new(EquilibratedFactorisation {
            inner: self.inner.factorise(a)?,
        }))
    }
}

struct EquilibratedFactorisation {
    inner: Box<dyn Factorisation>,
}

impl Factorisation for EquilibratedFactorisation {
    fn refactor_and_solve(&mut self, a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError> {
        let (dr, dc) = equilibration_factors(a);
        let (a_s, b_s) = apply_equilibration(a, b, &dr, &dc);
        let x_s = self.inner.refactor_and_solve(&a_s, &b_s)?;
        Ok((0..x_s.len()).map(|j| dc[j] * x_s[j]).collect())
    }

    fn refactor_and_solve_transpose(
        &mut self,
        a: &[Vec<f64>],
        b: &[f64],
    ) -> Result<Vec<f64>, SimError> {
        self.inner.refactor_and_solve_transpose(a, b)
    }
}

// ---------------------------------------------------------------------------
// Condition-number estimate (2-norm, power iteration)
// ---------------------------------------------------------------------------

/// Estimate the 2-norm condition number κ₂(A) = σ_max / σ_min.
///
/// `σ_max` comes from power iteration on `AᵀA`; `σ_min` from inverse power
/// iteration `(AᵀA)⁻¹ = A⁻¹A⁻ᵀ` using a dense LU for the inner solves. Returns
/// `None` if `A` is singular (the inverse iteration fails) or empty. This is a
/// diagnostic — it builds dense LUs, so it is opt-in (`.options cond_estimate`).
pub fn estimate_condition_2norm(a: &[Vec<f64>]) -> Option<f64> {
    let n = a.len();
    if n == 0 {
        return Some(1.0);
    }
    let mat = |v: &[f64]| -> Vec<f64> {
        (0..n)
            .map(|i| (0..n).map(|j| a[i][j] * v[j]).sum())
            .collect()
    };
    let mat_t = |v: &[f64]| -> Vec<f64> {
        (0..n)
            .map(|j| (0..n).map(|i| a[i][j] * v[i]).sum())
            .collect()
    };
    let l2 = |v: &[f64]| -> f64 { v.iter().map(|x| x * x).sum::<f64>().sqrt() };
    let normed = |v: Vec<f64>| -> Vec<f64> {
        let nrm = l2(&v);
        if nrm > 0.0 {
            v.iter().map(|x| x / nrm).collect()
        } else {
            v
        }
    };
    // Deterministic seed (no RNG): a simple non-uniform vector.
    let seed: Vec<f64> = (0..n).map(|i| 1.0 + (i % 7) as f64 * 0.13).collect();

    // σ_max² = λ_max(AᵀA): power iteration v ← AᵀA v / ‖·‖.
    let mut v = normed(seed.clone());
    let mut lam_max = 0.0;
    for _ in 0..100 {
        let av = mat(&v);
        let atav = mat_t(&av);
        lam_max = l2(&atav);
        v = normed(atav);
    }
    let sigma_max = lam_max.sqrt();

    // σ_min² = 1/λ_max((AᵀA)⁻¹): inverse power iteration w ← A⁻¹A⁻ᵀ w / ‖·‖.
    let solver = DenseSolver;
    let mut w = normed(seed);
    let mut lam_inv_max = 0.0;
    for _ in 0..100 {
        let u = solver.solve_transpose(a, &w).ok()?;
        let z = solver.solve(a, &u).ok()?;
        lam_inv_max = l2(&z);
        w = normed(z);
    }
    if lam_inv_max <= 0.0 {
        return None;
    }
    let sigma_min = 1.0 / lam_inv_max.sqrt();
    if sigma_min <= 0.0 || !sigma_max.is_finite() {
        return None;
    }
    Some(sigma_max / sigma_min)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dense_matches_sparse_on_small_matrix() {
        let a = vec![
            vec![4.0, -1.0, 0.0, 0.0],
            vec![-1.0, 4.0, -1.0, 0.0],
            vec![0.0, -1.0, 4.0, -1.0],
            vec![0.0, 0.0, -1.0, 3.0],
        ];
        let b = vec![1.0, 2.0, 3.0, 4.0];
        let x_d = DenseSolver.solve(&a, &b).unwrap();
        let x_s = FaerSparseSolver::default().solve(&a, &b).unwrap();
        for i in 0..4 {
            assert!(
                (x_d[i] - x_s[i]).abs() < 1e-10,
                "i={i}: dense={} sparse={}",
                x_d[i],
                x_s[i]
            );
        }
    }

    #[test]
    fn singular_matrix_returns_err() {
        let a = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        let b = vec![1.0, 2.0];
        assert!(matches!(
            DenseSolver.solve(&a, &b),
            Err(SimError::SingularMatrix)
        ));
        assert!(matches!(
            FaerSparseSolver::default().solve(&a, &b),
            Err(SimError::SingularMatrix)
        ));
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
        let x_via_method = DenseSolver.solve_transpose(&a, &b).unwrap();
        for i in 0..3 {
            assert!((x_via_explicit_t[i] - x_via_method[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn faer_sparse_factorisation_reuses_pattern() {
        // First solve through factorise(); second solve through the
        // refactor path with updated values — answer must match a fresh
        // dense solve of the second matrix.
        let a1 = vec![vec![3.0, 1.0], vec![0.0, 4.0]];
        let solver = FaerSparseSolver::default();
        let mut fact = solver.factorise(&a1).unwrap();
        let x1 = fact.refactor_and_solve(&a1, &[5.0, 8.0]).unwrap();
        // 3x+y=5, 4y=8 → y=2, x=1
        assert!((x1[0] - 1.0).abs() < 1e-12);
        assert!((x1[1] - 2.0).abs() < 1e-12);

        // Change values, same pattern.
        let a2 = vec![vec![6.0, 2.0], vec![0.0, 4.0]];
        let x2 = fact.refactor_and_solve(&a2, &[10.0, 8.0]).unwrap();
        // 6x+2y=10, 4y=8 → y=2, 6x=6 → x=1
        assert!((x2[0] - 1.0).abs() < 1e-12);
        assert!((x2[1] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn condition_estimate_diagonal_and_identity() {
        // Identity → κ = 1.
        let id = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let k = estimate_condition_2norm(&id).unwrap();
        assert!((k - 1.0).abs() < 1e-6, "identity κ={k}");
        // diag(1, 1000) → κ = 1000.
        let d = vec![vec![1.0, 0.0], vec![0.0, 1000.0]];
        let k = estimate_condition_2norm(&d).unwrap();
        assert!((k - 1000.0).abs() / 1000.0 < 1e-3, "diag κ={k}");
    }

    #[test]
    fn equilibration_solves_badly_scaled_system_correctly() {
        // Badly-scaled system: row 0 ~1e6, row 1 ~1e-6. Exact solution x=[1,1].
        let a = vec![vec![2.0e6, 1.0e6], vec![1.0e-6, 3.0e-6]];
        let b = vec![3.0e6, 4.0e-6];
        let eq = EquilibratedSolver::new(Box::new(DenseSolver));
        let x = eq.solve(&a, &b).unwrap();
        assert!((x[0] - 1.0).abs() < 1e-6, "x0={}", x[0]);
        assert!((x[1] - 1.0).abs() < 1e-6, "x1={}", x[1]);
        // Equilibrated answer matches the plain dense solve (scaling is exact).
        let x_plain = DenseSolver.solve(&a, &b).unwrap();
        assert!((x[0] - x_plain[0]).abs() < 1e-9);
        assert!((x[1] - x_plain[1]).abs() < 1e-9);
    }

    #[test]
    fn equilibrated_factorisation_matches_direct() {
        let a = vec![vec![3.0, 1.0], vec![0.0, 4.0]];
        let eq = EquilibratedSolver::new(Box::new(DenseSolver));
        let mut fact = eq.factorise(&a).unwrap();
        let x = fact.refactor_and_solve(&a, &[5.0, 8.0]).unwrap();
        assert!((x[0] - 1.0).abs() < 1e-12 && (x[1] - 2.0).abs() < 1e-12);
    }
}
