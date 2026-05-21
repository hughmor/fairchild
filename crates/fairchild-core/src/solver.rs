//! Linear-system solver backends.
//!
//! Every analysis (DC OP, tran, AC, noise, DC sweep) reduces to a sequence
//! of dense matrix solves `A·x = b`.  This module defines a small trait so
//! the actual factorisation strategy can be swapped per simulation:
//!
//!  - `DenseSolver`         — faer partial-pivot LU (the historic path).  Best
//!                            for ≤ ~50 nodes where the sparse setup cost
//!                            dominates.
//!  - `FaerSparseSolver`    — faer sparse LU (`SparseColMat::sp_lu`).  Pure
//!                            Rust, no C deps; the default at larger N.
//!  - (planned) KLU         — SuiteSparse FFI, behind a `klu` cargo feature.
//!
//! `SolverKind::Auto` picks dense or sparse based on system size at the
//! analysis entry point.  Adjoint (`solve_transpose`) is provided so noise /
//! sensitivity paths can re-use the same backend.

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
/// linear solver: forward solve and (for adjoint / noise) transpose solve.
///
/// **Adding a third-party backend** — e.g. KLU via SuiteSparse FFI — is:
///   1. Add a new dependency behind a cargo feature flag.
///   2. Implement this trait on a new struct (see `FaerSparseSolver` for
///      the reference shape).  Override `solve_transpose` if your
///      factorisation caches an LU so transpose-solve is free.
///   3. Extend `SolverKind` with a `Klu` variant gated on the feature flag.
///   4. Add the dispatch to `make_solver`.
///
/// The trait is `Send + Sync` so future parallel-frequency AC / parallel-MC
/// drivers can share one solver across threads.
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
}

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
}

/// Sparse LU via faer's `SparseColMat::sp_lu`.
///
/// Builds a CSC matrix from non-zero entries of the dense input on each call
/// — that's an O(n²) scan but the LU step itself runs at sparse complexity,
/// which dominates above ~50 nodes for circuit matrices.  A future
/// commit will let `stamp_netlist` emit triplets directly so this scan
/// disappears entirely.
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

        // Collect (row, col, val) triplets, dropping near-zeros so faer's
        // symbolic phase actually sees the sparsity.
        let mut triplets: Vec<Triplet<usize, usize, f64>> = Vec::new();
        for i in 0..n {
            for j in 0..n {
                let v = a[i][j];
                if v.abs() > self.zero_threshold {
                    triplets.push(Triplet::new(i, j, v));
                }
            }
        }

        // SparseColMat owns CSC storage.  Sum duplicate triplets (none here,
        // but the API is symmetric with the future triplet-stamp path).
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
}

/// SuiteSparse KLU backend — sparse direct LU with BTF preordering.
///
/// Like `FaerSparseSolver`, this currently converts the dense input to
/// CSC on every call; the symbolic / numeric factorisation split that
/// makes KLU's repeated-pattern fast path pay off is added by the
/// follow-on `LinearSolver` trait extension (factorisation cache).
#[cfg(feature = "klu")]
pub struct KluSolver;

#[cfg(feature = "klu")]
impl LinearSolver for KluSolver {
    fn solve(&self, a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError> {
        fairchild_klu::klu_solve_dense(a, b).map_err(|_| SimError::SingularMatrix)
    }
}

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
        // Small 4×4 system; dense and sparse must agree to machine eps.
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
        // Use a non-symmetric A; verify solve(A^T, b) == solve_transpose(A, b).
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
}
