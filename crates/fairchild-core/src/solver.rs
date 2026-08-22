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

use crate::error::SimError;
use crate::mna::{MnaMatrix, SparseRow};

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
    fn solve(&self, a: &[SparseRow], b: &[f64]) -> Result<Vec<f64>, SimError>;

    /// Solve `A^T · x = b`.  Default falls back to constructing the explicit
    /// transpose and reusing `solve` — backends with cached factorisations
    /// override this to avoid the O(n²) copy.
    fn solve_transpose(&self, a: &[SparseRow], b: &[f64]) -> Result<Vec<f64>, SimError> {
        let at = transpose_sparse(a, b.len());
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
    fn factorise(&self, a: &[SparseRow]) -> Result<Box<dyn Factorisation>, SimError> {
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
    fn refactor_and_solve(&mut self, a: &[SparseRow], b: &[f64]) -> Result<Vec<f64>, SimError>;

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
        a: &[SparseRow],
        b: &[f64],
    ) -> Result<Vec<f64>, SimError> {
        let at = transpose_sparse(a, b.len());
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
    #[cfg(feature = "klu")]
    fn with_backend(backend: Box<dyn LinearSolver>) -> Self {
        Self { backend }
    }
}

impl Factorisation for NoCacheFactorisation {
    fn refactor_and_solve(&mut self, a: &[SparseRow], b: &[f64]) -> Result<Vec<f64>, SimError> {
        self.backend.solve(a, b)
    }
    fn refactor_and_solve_transpose(
        &mut self,
        a: &[SparseRow],
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
    fn solve(&self, a: &[SparseRow], b: &[f64]) -> Result<Vec<f64>, SimError> {
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

    fn factorise(&self, _a: &[SparseRow]) -> Result<Box<dyn Factorisation>, SimError> {
        Ok(Box::new(DenseFactorisation::new()))
    }
}

/// Dense LU that survives an unchanged matrix.
///
/// The cache here is the **numeric factors**, not the pivot order — caching the
/// pivots alone genuinely does not pay, which is what the no-op cache this
/// replaced was reasoning about. Reusing the factors is a different trade: an LU
/// is O(n³) and the triangular solves that follow are O(n²), so on a run where
/// `A` holds still (a linear circuit at fixed timestep, where the MNA matrix is
/// constant for the whole transient) this turns every step after the first into
/// solves alone.
struct DenseFactorisation {
    /// Row-major values the cached factors were built from.
    values: Vec<f64>,
    n: usize,
    lu: Option<faer::linalg::solvers::PartialPivLu<f64>>,
}

impl DenseFactorisation {
    fn new() -> Self {
        Self {
            values: Vec::new(),
            n: 0,
            lu: None,
        }
    }
}

impl Factorisation for DenseFactorisation {
    fn refactor_and_solve(&mut self, a: &[SparseRow], b: &[f64]) -> Result<Vec<f64>, SimError> {
        let n = b.len();
        if n == 0 {
            return Ok(Vec::new());
        }
        if self.n != n {
            self.n = n;
            self.values = vec![f64::NAN; n * n]; // NAN != anything, so this reads as changed
            self.lu = None;
        }

        // One walk that both refreshes the cache and notices whether it moved.
        // Same `row[j]` access the one-shot `solve` uses, so the two paths
        // cannot disagree about what the matrix is.
        let mut changed = false;
        for (row, dst) in a.iter().zip(self.values.chunks_mut(n)) {
            for (j, slot) in dst.iter_mut().enumerate() {
                let v = row[j];
                changed |= *slot != v;
                *slot = v;
            }
        }

        if changed || self.lu.is_none() {
            let a_mat = Mat::<f64>::from_fn(n, n, |i, j| self.values[i * n + j]);
            self.lu = Some(a_mat.partial_piv_lu());
        }
        let lu = self.lu.as_ref().expect("just populated");
        let b_col = Col::<f64>::from_fn(n, |i| b[i]);
        let x_col = lu.solve(b_col.as_ref());
        let x: Vec<f64> = (0..n).map(|i| x_col[i]).collect();
        if x.iter().any(|v| !v.is_finite()) {
            return Err(SimError::SingularMatrix);
        }
        Ok(x)
    }

    fn refactor_and_solve_transpose(
        &mut self,
        a: &[SparseRow],
        b: &[f64],
    ) -> Result<Vec<f64>, SimError> {
        // Deliberately uncached. The adjoint paths are cold, and routing Aᵀ
        // through the forward cache would make A and Aᵀ evict each other on
        // every alternation — slower than not caching, for no gain.
        DenseSolver.solve_transpose(a, b)
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
    fn solve(&self, a: &[SparseRow], b: &[f64]) -> Result<Vec<f64>, SimError> {
        use faer::sparse::{SparseColMat, Triplet};

        let n = b.len();
        if n == 0 {
            return Ok(Vec::new());
        }

        let mut triplets: Vec<Triplet<usize, usize, f64>> = Vec::new();
        for (i, row) in a.iter().enumerate() {
            for (j, v) in row.iter() {
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

    /// Caching is sound straight off the rows now, which is why there is no
    /// `factorise_mat` override any more.
    ///
    /// It was not sound while the matrix was dense: discovering the structure
    /// from values on iteration 0 misses cells that are zero at x=0 and
    /// non-zero at the solution, which photonic devices do stamp (see the WDM
    /// DC OP regression) — hence the separate [`crate::mna::Pattern`] that used
    /// to be threaded in. A `SparseRow` built from that pattern *allocates*
    /// those cells at zero, so the rebuild walk sees them and records them as
    /// `NO_SLOT`; the same growth detection then fires on the iteration one
    /// turns non-zero. The pattern still does its job — it just does it once,
    /// when the matrix is constructed, instead of on every solve.
    fn factorise(&self, _a: &[SparseRow]) -> Result<Box<dyn Factorisation>, SimError> {
        Ok(Box::new(FaerSparseFactorisation::new(self.zero_threshold)))
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
    zero_threshold: f64,
    /// CSC column pointers / row indices of the active set.
    col_ptr: Vec<usize>,
    row_idx: Vec<usize>,
    /// Where in `values` each active `(row, col)` lives, indexed the same way
    /// the refill walk visits them: row-major over `pattern.cols`.
    slot: Vec<u32>,
    values: Vec<f64>,
    symbolic: Option<faer::sparse::linalg::solvers::SymbolicLu<usize>>,
    /// Numeric factors of the matrix currently in `values`.  `None` means they
    /// are stale (or never computed) and the next solve must factorise.
    numeric: Option<faer::sparse::linalg::solvers::Lu<usize, f64>>,
}

/// `slot` entry for a structural cell that is not in the active set.
const NO_SLOT: u32 = u32::MAX;

/// What a [`FaerSparseFactorisation::refill`] pass found, and therefore how much
/// of the cache survives.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Refill {
    /// The structural pattern changed — both symbolic and numeric are stale.
    Rebuild,
    /// Same pattern, at least one value moved — numeric factors are stale.
    Changed,
    /// Not one value moved.  The numeric factors are still exact; skip
    /// factorising and go straight to the triangular solves.
    Unchanged,
}

impl FaerSparseFactorisation {
    fn new(zero_threshold: f64) -> Self {
        FaerSparseFactorisation {
            zero_threshold,
            col_ptr: Vec::new(),
            row_idx: Vec::new(),
            slot: Vec::new(),
            values: Vec::new(),
            symbolic: None,
            numeric: None,
        }
    }

    /// Rebuild the CSC structure and symbolic LU from the cells of `a` that are
    /// currently non-zero.  Runs on the first solve and again only if a new
    /// cell inside the structural pattern turns non-zero.
    fn rebuild(&mut self, a: &[SparseRow]) -> Result<(), SimError> {
        use faer::sparse::linalg::solvers::SymbolicLu;
        use faer::sparse::SymbolicSparseColMatRef;

        let n = a.len();
        let thr = self.zero_threshold;
        // Count per column first so the CSC arrays can be filled in one pass.
        let mut col_count = vec![0usize; n];
        for row in a {
            for (j, v) in row.iter() {
                if v.abs() > thr {
                    col_count[j] += 1;
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
        let mut slot: Vec<u32> = Vec::new();
        for (i, row) in a.iter().enumerate() {
            for (j, v) in row.iter() {
                if v.abs() > thr {
                    let k = fill[j];
                    row_idx[k] = i;
                    values[k] = v;
                    fill[j] += 1;
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
        self.numeric = None; // new pattern — any cached factors are meaningless
        Ok(())
    }

    /// Copy current values into the cached CSC, reporting what that implies for
    /// the cached factors.
    ///
    /// The `Unchanged` case is the one worth having: an LU depends only on `A`,
    /// so if not one value moved, the existing factors are still exactly right
    /// and only the triangular solves need re-running.  Detecting it costs one
    /// comparison per stored value — negligible against the factorisation it
    /// skips, and it *checks* rather than assuming, so no caller can be wrong
    /// about whether its matrix is constant.
    fn refill(&mut self, a: &[SparseRow]) -> Refill {
        let thr = self.zero_threshold;
        let mut grew = false;
        let mut changed = false;
        let mut k = 0usize;
        // Positional: `slot` was built by this same walk, so reading the row's
        // own values needs no column lookup. Indexing `row[j]` here instead
        // would binary-search every cell now that rows are stored sparse.
        for row in a {
            let (_, vals) = row.entries();
            for &v in vals {
                match self.slot.get(k).copied() {
                    Some(NO_SLOT) => grew |= v.abs() > thr,
                    Some(s) => {
                        let slot = &mut self.values[s as usize];
                        // Bit-exact: a value that moved by any amount at all
                        // invalidates the factors.  Nothing here should tolerate
                        // "close enough" — that would silently reuse factors for
                        // a matrix they were not computed from.
                        changed |= *slot != v;
                        *slot = v;
                    }
                    None => return Refill::Rebuild,
                }
                k += 1;
            }
        }
        if grew || k != self.slot.len() {
            Refill::Rebuild
        } else if changed {
            Refill::Changed
        } else {
            Refill::Unchanged
        }
    }

    /// Solve with the cached factors, computing them first if they are absent.
    fn solve_cached(&mut self, b: &[f64]) -> Result<Vec<f64>, SimError> {
        use faer::sparse::linalg::solvers::Lu;
        use faer::sparse::{SparseColMatRef, SymbolicSparseColMatRef};

        let n = b.len();
        if self.numeric.is_none() {
            let symbolic = self.symbolic.clone().ok_or(SimError::SingularMatrix)?;
            let sym =
                SymbolicSparseColMatRef::new_checked(n, n, &self.col_ptr, None, &self.row_idx);
            let mat = SparseColMatRef::<usize, f64>::new(sym, &self.values);
            self.numeric = Some(
                Lu::try_new_with_symbolic(symbolic, mat).map_err(|_| SimError::SingularMatrix)?,
            );
        }
        let lu = self.numeric.as_ref().expect("just populated");
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
    fn refactor_and_solve(&mut self, a: &[SparseRow], b: &[f64]) -> Result<Vec<f64>, SimError> {
        // `refill` indexes the cached slot map, so it is only valid once
        // `rebuild` has run; short-circuit ordering matters here.
        if self.symbolic.is_none() {
            self.rebuild(a)?;
        } else {
            match self.refill(a) {
                Refill::Rebuild => self.rebuild(a)?,
                Refill::Changed => self.numeric = None,
                // Factors still describe this exact matrix — keep them.
                Refill::Unchanged => {}
            }
        }
        self.solve_cached(b)
    }

    fn refactor_and_solve_transpose(
        &mut self,
        a: &[SparseRow],
        b: &[f64],
    ) -> Result<Vec<f64>, SimError> {
        // Adjoint paths are cold; the explicit transpose keeps this simple.
        let at = transpose_sparse(a, b.len());
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
/// behaviour as `FaerSparseSolver::solve`).  The `factorise_mat` path is
/// where KLU's advantage lives: a `KluSymbolic` analysed once is reused via
/// `klu_refactor` for every subsequent solve.
///
/// That advantage was theoretical until 2026-08-01. The cache reused the
/// symbolic factorisation correctly, then threw the win away by rebuilding the
/// CSC arrays from a **full dense O(n²) scan** on every `refactor_and_solve` —
/// and in column-major order over row-major storage, so it fell out of cache
/// as well. Measured: 41 ms per call at n=3200, against a ~120 ms total excess
/// over faer-sparse across ~3 Newton iterations. Net effect, KLU was 4.1×
/// *slower* than the pure-Rust backend it was supposed to beat.
///
/// It now uses the same structural-pattern slot map as
/// [`FaerSparseFactorisation`], which `d149a26` gave that backend and this one
/// never received. Measured after: 1.2-1.5× faster than faer-sparse, widening
/// with size (`cargo bench -p fairchild-core --features klu`).
#[cfg(feature = "klu")]
pub struct KluSolver;

#[cfg(feature = "klu")]
impl LinearSolver for KluSolver {
    fn solve(&self, a: &[SparseRow], b: &[f64]) -> Result<Vec<f64>, SimError> {
        let dense = crate::mna::CircuitTopology::to_dense(a, b.len());
        fairchild_klu::klu_solve_dense(&dense, b).map_err(|_| SimError::SingularMatrix)
    }

    /// Sparse rows carry their own structure, so this needs no `MnaMatrix` and no
    /// separate [`crate::mna::Pattern`] — which is why there is no
    /// `factorise_mat` override any more.
    fn factorise(&self, a: &[SparseRow]) -> Result<Box<dyn Factorisation>, SimError> {
        let n = a.len();
        if n == 0 {
            return Ok(Box::new(NoCacheFactorisation::with_backend(Box::new(
                KluSolver,
            ))));
        }
        let mut f = KluFactorisation::new(1e-30, n);
        f.rebuild(a)?;
        Ok(Box::new(f))
    }
}

#[cfg(feature = "klu")]
struct KluFactorisation {
    common: fairchild_klu::KluCommon,
    /// `None` until the first `rebuild`.
    symbolic: Option<fairchild_klu::KluSymbolic>,
    numeric: Option<fairchild_klu::KluNumeric>,
    zero_threshold: f64,
    ap: Vec<i32>,
    ai: Vec<i32>,
    ax: Vec<f64>,
    /// Where in `ax` each allocated cell lives, in the order the refill walk
    /// visits them — row-major over the rows' own entries, so refill is a
    /// positional copy with no lookup. `NO_SLOT` = not currently in the
    /// numerically active set.
    slot: Vec<u32>,
    n: usize,
}

#[cfg(feature = "klu")]
impl KluFactorisation {
    fn new(zero_threshold: f64, n: usize) -> Self {
        KluFactorisation {
            common: fairchild_klu::KluCommon::new(),
            symbolic: None,
            numeric: None,
            zero_threshold,
            ap: Vec::new(),
            ai: Vec::new(),
            ax: Vec::new(),
            slot: Vec::new(),
            n,
        }
    }

    /// Build the CSC arrays and the slot map from the currently non-zero cells,
    /// then run `klu_analyze` + `klu_factor`.
    ///
    /// The alternative this replaced — rebuilding CSC from a full dense scan on
    /// every solve — cost O(n²) with a column-major traversal of row-major
    /// storage, measured at 41 ms per call at n = 3200, which swamped
    /// everything `klu_refactor` was saving. Now that the matrix is stored
    /// sparse the walk is over the rows' own entries, so both the structure and
    /// the values come out in one O(nnz) pass.
    fn rebuild(&mut self, a: &[SparseRow]) -> Result<(), SimError> {
        use fairchild_klu::{KluNumeric, KluSymbolic};
        let n = self.n;
        let thr = self.zero_threshold;

        let mut col_count = vec![0i32; n];
        for row in a {
            for (j, v) in row.iter() {
                if v.abs() > thr {
                    col_count[j] += 1;
                }
            }
        }
        let mut ap = vec![0i32; n + 1];
        for j in 0..n {
            ap[j + 1] = ap[j] + col_count[j];
        }
        let nnz = ap[n] as usize;
        let mut ai = vec![0i32; nnz];
        let mut ax = vec![0.0f64; nnz];
        let mut fill: Vec<i32> = ap.clone();
        // Row-major walk keeps row indices ascending within each column, which
        // is what KLU expects of a sorted CSC.
        let mut slot: Vec<u32> = Vec::new();
        for (i, row) in a.iter().enumerate() {
            for (j, v) in row.iter() {
                if v.abs() > thr {
                    let k = fill[j] as usize;
                    ai[k] = i as i32;
                    ax[k] = v;
                    fill[j] += 1;
                    slot.push(k as u32);
                } else {
                    slot.push(NO_SLOT);
                }
            }
        }

        // Drop the old handles before allocating new ones (Drop runs
        // klu_free_numeric / klu_free_symbolic).
        self.numeric = None;
        self.symbolic = None;
        let symbolic = KluSymbolic::analyze(n, &mut ap, &mut ai, &mut self.common)
            .map_err(|_| SimError::SingularMatrix)?;
        let numeric = KluNumeric::factor(&mut ap, &mut ai, &mut ax, &symbolic, &mut self.common)
            .map_err(|_| SimError::SingularMatrix)?;
        self.ap = ap;
        self.ai = ai;
        self.ax = ax;
        self.slot = slot;
        self.symbolic = Some(symbolic);
        self.numeric = Some(numeric);
        Ok(())
    }

    /// Copy current values into the cached CSC, reporting how much of the cache
    /// survives — same contract as [`FaerSparseFactorisation::refill`].
    ///
    /// Purely positional: `slot` was built by this same walk order, so no
    /// column lookup happens here at all.
    fn refill(&mut self, a: &[SparseRow]) -> Refill {
        let thr = self.zero_threshold;
        let mut grew = false;
        let mut changed = false;
        let mut k = 0usize;
        for row in a {
            let (_, vals) = row.entries();
            for &v in vals {
                match self.slot.get(k).copied() {
                    Some(NO_SLOT) => grew |= v.abs() > thr,
                    Some(s) => {
                        let slot = &mut self.ax[s as usize];
                        changed |= *slot != v;
                        *slot = v;
                    }
                    // The row grew a column since the rebuild — only possible
                    // on a patternless matrix, and it means rebuild anyway.
                    None => return Refill::Rebuild,
                }
                k += 1;
            }
        }
        // Fewer cells than the slot map expects: same conclusion.
        if grew || k != self.slot.len() {
            Refill::Rebuild
        } else if changed {
            Refill::Changed
        } else {
            Refill::Unchanged
        }
    }

    /// Solve on the cached pattern, numerically refactoring first only when the
    /// matrix actually moved. `transpose` picks `klu_tsolve` for the adjoint
    /// path.
    fn refactor_then(
        &mut self,
        b: &[f64],
        transpose: bool,
        needs_refactor: bool,
    ) -> Result<Vec<f64>, SimError> {
        let symbolic = self.symbolic.as_ref().ok_or(SimError::SingularMatrix)?;
        let numeric = self.numeric.as_mut().ok_or(SimError::SingularMatrix)?;
        if needs_refactor {
            numeric
                .refactor(
                    &mut self.ap,
                    &mut self.ai,
                    &mut self.ax,
                    symbolic,
                    &mut self.common,
                )
                .map_err(|_| SimError::SingularMatrix)?;
        }
        let mut x = b.to_vec();
        // `klu_tsolve` solves Aᵀx = b from A's own factorisation, so the
        // adjoint path reuses the cache instead of materialising a dense
        // transpose the way it used to (and the way the faer path still does).
        let r = if transpose {
            numeric.solve_transpose(symbolic, &mut x, &mut self.common)
        } else {
            numeric.solve(symbolic, &mut x, &mut self.common)
        };
        r.map_err(|_| SimError::SingularMatrix)?;
        if x.iter().any(|v| !v.is_finite()) {
            return Err(SimError::SingularMatrix);
        }
        Ok(x)
    }

    fn solve_with(
        &mut self,
        a: &[SparseRow],
        b: &[f64],
        transpose: bool,
    ) -> Result<Vec<f64>, SimError> {
        // `rebuild` ends in a fresh `klu_factor` over the current values, so
        // both rebuild paths leave the numeric factors already valid.
        let needs_refactor = if self.symbolic.is_none() {
            self.rebuild(a)?;
            false
        } else {
            match self.refill(a) {
                Refill::Rebuild => {
                    self.rebuild(a)?;
                    false
                }
                Refill::Changed => true,
                Refill::Unchanged => false,
            }
        };
        self.refactor_then(b, transpose, needs_refactor)
    }
}

#[cfg(feature = "klu")]
impl Factorisation for KluFactorisation {
    fn refactor_and_solve(&mut self, a: &[SparseRow], b: &[f64]) -> Result<Vec<f64>, SimError> {
        self.solve_with(a, b, false)
    }

    fn refactor_and_solve_transpose(
        &mut self,
        a: &[SparseRow],
        b: &[f64],
    ) -> Result<Vec<f64>, SimError> {
        self.solve_with(a, b, true)
    }
}

// ---------------------------------------------------------------------------
// Dispatch
// ---------------------------------------------------------------------------

/// Largest system `Auto` hands to the dense backend.
///
/// Measured, both directions, best of 3 (2026-08-17). `n` is a stand-in for the
/// thing that actually decides this — how sparse the matrix is — and the two
/// circuit families bracket it: a ring oscillator's MOSFETs couple four nodes
/// each, an RC ladder is tridiagonal.
///
/// | nodes | nonlinear: cost of sparse | linear: cost of dense |
/// |---|---|---|
/// | ~7–11 | +35 % | +9…16 % |
/// | ~21–23 | +12 % | +47 % |
/// | ~43 | −3 % (sparse ahead) | +115 % |
///
/// So the crossover sits near 20, not the 50 this used to be: past ~20 the
/// penalty for guessing sparse is small and bounded while the penalty for
/// guessing dense keeps growing with sparsity. Below it, dense's fit is real
/// enough to keep. Both backends reuse numeric factors now, so this is a
/// judgement about algorithmic fit to `n` and nothing else.
const AUTO_DENSE_MAX: usize = 20;

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
            if n < AUTO_DENSE_MAX {
                Box::new(DenseSolver)
            } else {
                Box::new(FaerSparseSolver::default())
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Matrix equilibration (two-sided Ruiz scaling)
// ---------------------------------------------------------------------------

/// Two-sided ∞-norm (Ruiz) scaling factors `(D_r, D_c)` such that the rows and
/// columns of `diag(D_r)·A·diag(D_c)` have ∞-norm ≈ 1.  A few iterations is
/// enough in practice; empty rows/columns are left unscaled (factor 1).
fn equilibration_factors(a: &[SparseRow]) -> (Vec<f64>, Vec<f64>) {
    let n = a.len();
    let mut dr = vec![1.0_f64; n];
    let mut dc = vec![1.0_f64; n];
    // Row and column ∞-norms only ever see non-zero cells, so both sweeps walk
    // the sparse entries — O(nnz) per pass instead of O(n²).
    for _ in 0..3 {
        for (i, row) in a.iter().enumerate() {
            let mut m = 0.0_f64;
            for (j, v) in row.iter() {
                m = m.max((dr[i] * v * dc[j]).abs());
            }
            if m > 0.0 {
                dr[i] /= m.sqrt();
            }
        }
        let mut col_max = vec![0.0_f64; n];
        for (i, row) in a.iter().enumerate() {
            for (j, v) in row.iter() {
                col_max[j] = col_max[j].max((dr[i] * v * dc[j]).abs());
            }
        }
        for (j, m) in col_max.into_iter().enumerate() {
            if m > 0.0 {
                dc[j] /= m.sqrt();
            }
        }
    }
    (dr, dc)
}

/// `A' = diag(D_r)·A·diag(D_c)`, `b' = D_r·b`.  Returns the scaled system.
fn apply_equilibration(
    a: &[SparseRow],
    b: &[f64],
    dr: &[f64],
    dc: &[f64],
) -> (Vec<SparseRow>, Vec<f64>) {
    let n = a.len();
    // Diagonal scaling is elementwise, so the structure carries over untouched.
    let a_s: Vec<SparseRow> = a
        .iter()
        .enumerate()
        .map(|(i, row)| {
            SparseRow::from_sorted_cells(
                row.iter()
                    .map(|(j, v)| (j as u32, dr[i] * v * dc[j]))
                    .collect(),
            )
        })
        .collect();
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

/// Sparse transpose. The adjoint paths are cold, so this stays a plain rebuild
/// rather than a cached structure.
pub(crate) fn transpose_sparse(a: &[SparseRow], n: usize) -> Vec<SparseRow> {
    let mut cells: Vec<Vec<(u32, f64)>> = vec![Vec::new(); n];
    for (i, row) in a.iter().enumerate() {
        for (j, v) in row.iter() {
            cells[j].push((i as u32, v));
        }
    }
    cells
        .into_iter()
        .map(SparseRow::from_sorted_cells)
        .collect()
}

impl LinearSolver for EquilibratedSolver {
    fn solve(&self, a: &[SparseRow], b: &[f64]) -> Result<Vec<f64>, SimError> {
        let (dr, dc) = equilibration_factors(a);
        let (a_s, b_s) = apply_equilibration(a, b, &dr, &dc);
        let x_s = self.inner.solve(&a_s, &b_s)?;
        Ok((0..x_s.len()).map(|j| dc[j] * x_s[j]).collect())
    }

    fn solve_transpose(&self, a: &[SparseRow], b: &[f64]) -> Result<Vec<f64>, SimError> {
        // Adjoint path left unscaled (forward-only equilibration).
        self.inner.solve_transpose(a, b)
    }

    fn factorise(&self, a: &[SparseRow]) -> Result<Box<dyn Factorisation>, SimError> {
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
    fn refactor_and_solve(&mut self, a: &[SparseRow], b: &[f64]) -> Result<Vec<f64>, SimError> {
        let (dr, dc) = equilibration_factors(a);
        let (a_s, b_s) = apply_equilibration(a, b, &dr, &dc);
        let x_s = self.inner.refactor_and_solve(&a_s, &b_s)?;
        Ok((0..x_s.len()).map(|j| dc[j] * x_s[j]).collect())
    }

    fn refactor_and_solve_transpose(
        &mut self,
        a: &[SparseRow],
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
pub fn estimate_condition_2norm(a: &[SparseRow]) -> Option<f64> {
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
    /// These fixtures are written as dense literals because that is the
    /// readable way to write a 2×2; the solvers take sparse rows.
    fn sp(a: &[Vec<f64>]) -> Vec<SparseRow> {
        crate::mna::CircuitTopology::sparse_from_dense(a)
    }

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
        let x_d = DenseSolver.solve(&sp(&a), &b).unwrap();
        let x_s = FaerSparseSolver::default().solve(&sp(&a), &b).unwrap();
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
            DenseSolver.solve(&sp(&a), &b),
            Err(SimError::SingularMatrix)
        ));
        assert!(matches!(
            FaerSparseSolver::default().solve(&sp(&a), &b),
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

        let x_via_explicit_t = DenseSolver.solve(&sp(&at), &b).unwrap();
        let x_via_method = DenseSolver.solve_transpose(&sp(&a), &b).unwrap();
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
        let mut fact = solver.factorise(&sp(&a1)).unwrap();
        let x1 = fact.refactor_and_solve(&sp(&a1), &[5.0, 8.0]).unwrap();
        // 3x+y=5, 4y=8 → y=2, x=1
        assert!((x1[0] - 1.0).abs() < 1e-12);
        assert!((x1[1] - 2.0).abs() < 1e-12);

        // Change values, same pattern.
        let a2 = vec![vec![6.0, 2.0], vec![0.0, 4.0]];
        let x2 = fact.refactor_and_solve(&sp(&a2), &[10.0, 8.0]).unwrap();
        // 6x+2y=10, 4y=8 → y=2, 6x=6 → x=1
        assert!((x2[0] - 1.0).abs() < 1e-12);
        assert!((x2[1] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn condition_estimate_diagonal_and_identity() {
        // Identity → κ = 1.
        let id = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let k = estimate_condition_2norm(&sp(&id)).unwrap();
        assert!((k - 1.0).abs() < 1e-6, "identity κ={k}");
        // diag(1, 1000) → κ = 1000.
        let d = vec![vec![1.0, 0.0], vec![0.0, 1000.0]];
        let k = estimate_condition_2norm(&sp(&d)).unwrap();
        assert!((k - 1000.0).abs() / 1000.0 < 1e-3, "diag κ={k}");
    }

    #[test]
    fn equilibration_solves_badly_scaled_system_correctly() {
        // Badly-scaled system: row 0 ~1e6, row 1 ~1e-6. Exact solution x=[1,1].
        let a = vec![vec![2.0e6, 1.0e6], vec![1.0e-6, 3.0e-6]];
        let b = vec![3.0e6, 4.0e-6];
        let eq = EquilibratedSolver::new(Box::new(DenseSolver));
        let x = eq.solve(&sp(&a), &b).unwrap();
        assert!((x[0] - 1.0).abs() < 1e-6, "x0={}", x[0]);
        assert!((x[1] - 1.0).abs() < 1e-6, "x1={}", x[1]);
        // Equilibrated answer matches the plain dense solve (scaling is exact).
        let x_plain = DenseSolver.solve(&sp(&a), &b).unwrap();
        assert!((x[0] - x_plain[0]).abs() < 1e-9);
        assert!((x[1] - x_plain[1]).abs() < 1e-9);
    }

    #[test]
    fn equilibrated_factorisation_matches_direct() {
        let a = vec![vec![3.0, 1.0], vec![0.0, 4.0]];
        let eq = EquilibratedSolver::new(Box::new(DenseSolver));
        let mut fact = eq.factorise(&sp(&a)).unwrap();
        let x = fact.refactor_and_solve(&sp(&a), &[5.0, 8.0]).unwrap();
        assert!((x[0] - 1.0).abs() < 1e-12 && (x[1] - 2.0).abs() < 1e-12);
    }

    /// `refill` must distinguish "not one value moved" from "a value moved",
    /// because that is the only thing standing between reusing the numeric
    /// factors and reusing them when they are stale.
    #[test]
    fn refill_reports_unchanged_only_when_nothing_moved() {
        let a = vec![vec![4.0, -1.0], vec![-1.0, 3.0]];
        let mut f = FaerSparseFactorisation::new(0.0);
        f.rebuild(&sp(&a)).unwrap();

        assert_eq!(f.refill(&sp(&a)), Refill::Unchanged, "identical matrix");

        let moved = vec![vec![4.0, -1.0], vec![-1.0, 3.0 + 1e-15]];
        assert_eq!(
            f.refill(&sp(&moved)),
            Refill::Changed,
            "a 1e-15 move is still a different matrix"
        );

        // A cell outside the active set turning non-zero needs the pattern back.
        let grown = vec![vec![4.0, -1.0], vec![-1.0, 3.0]];
        let mut g = FaerSparseFactorisation::new(1.0); // threshold hides the -1s
        g.rebuild(&sp(&grown)).unwrap();
        let big = vec![vec![4.0, -2.0], vec![-2.0, 3.0]];
        assert_eq!(g.refill(&sp(&big)), Refill::Rebuild, "pattern grew");
    }

    /// Every backend that caches factors, by name, so a failure says which one.
    /// Dense is included deliberately: it caches nothing, and these two
    /// properties must hold for it too.
    fn caching_backends() -> Vec<(&'static str, Box<dyn LinearSolver>)> {
        #[allow(unused_mut)] // `mut` is only needed with the `klu` feature on
        let mut v: Vec<(&'static str, Box<dyn LinearSolver>)> = vec![
            ("dense", Box::new(DenseSolver)),
            ("faer-sparse", Box::new(FaerSparseSolver::default())),
        ];
        #[cfg(feature = "klu")]
        v.push(("klu", Box::new(KluSolver)));
        v
    }

    /// The reuse must cache **factors**, not answers: same `A`, three different
    /// right-hand sides, each answer checked against the analytic solution.
    /// Caching a solution instead would pass a single-solve test and fail here.
    #[test]
    fn reused_factors_still_solve_new_right_hand_sides() {
        // [[2,0],[0,4]] — diagonal, so x = [b0/2, b1/4] by inspection.
        let a = vec![vec![2.0, 0.0], vec![0.0, 4.0]];
        for (name, solver) in caching_backends() {
            let mut fact = solver.factorise(&sp(&a)).unwrap();
            for (b, want) in [
                ([2.0, 4.0], [1.0, 1.0]),
                ([6.0, 8.0], [3.0, 2.0]),
                ([-4.0, 2.0], [-2.0, 0.5]),
            ] {
                let x = fact.refactor_and_solve(&sp(&a), &b).unwrap();
                assert!(
                    (x[0] - want[0]).abs() < 1e-12 && (x[1] - want[1]).abs() < 1e-12,
                    "{name}: b={b:?} got {x:?} want {want:?}"
                );
            }
        }
    }

    /// And when `A` *does* change between solves, the answer must follow it.
    /// This is the failure the reuse could introduce: stale factors give a
    /// plausible number for the previous matrix.
    #[test]
    fn changed_matrix_is_refactorised() {
        let two = sp(&[vec![2.0, 0.0], vec![0.0, 2.0]]);
        let four = sp(&[vec![4.0, 0.0], vec![0.0, 4.0]]);
        for (name, solver) in caching_backends() {
            let mut fact = solver.factorise(&two).unwrap();
            let x = fact.refactor_and_solve(&two, &[2.0, 2.0]).unwrap();
            assert!((x[0] - 1.0).abs() < 1e-12, "{name}: 2x=2 => x=1, got {x:?}");

            // Same pattern, different values: x must be 0.5, not the cached 1.0.
            let x = fact.refactor_and_solve(&four, &[2.0, 2.0]).unwrap();
            assert!(
                (x[0] - 0.5).abs() < 1e-12,
                "{name}: 4x=2 => x=0.5, got {x:?}"
            );
        }
    }
}
