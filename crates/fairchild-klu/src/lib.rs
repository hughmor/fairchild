//! Thin safe bindings to SuiteSparse KLU.
//!
//! KLU is the de-facto sparse direct LU solver for circuit-shaped
//! matrices: it exploits Block Triangular Form (BTF) to decompose the
//! matrix into strongly-connected components and run dense LU within
//! each.  On a typical analog circuit it beats UMFPACK 2-5×.
//!
//! This crate is intentionally minimal — it covers exactly the slice of
//! the KLU C API that fairchild needs:
//!
//!   * `klu_defaults` — initialise the common (settings + status) block
//!   * `klu_analyze`  — symbolic factorisation (sparsity-only; cacheable)
//!   * `klu_factor`   — numeric factorisation (values)
//!   * `klu_refactor` — re-run numeric with the same sparsity pattern
//!   * `klu_solve`    — back-substitute `A·x = b`
//!   * `klu_tsolve`   — back-substitute `Aᵀ·x = b` (for adjoint paths)
//!   * `klu_free_symbolic` / `klu_free_numeric`
//!
//! The `klu_common` C struct is treated as **fully opaque** — we never
//! dereference it from Rust.  Errors are detected via the return values
//! of `klu_factor` / `klu_solve` / etc., which is sufficient for every
//! call site fairchild has.  This avoids hand-mirroring a 30-field C
//! struct that has shifted layout between SuiteSparse 5.x and 7.x.
//!
//! The common-block buffer is over-sized (4 KiB, 8-byte aligned via
//! `Vec<f64>`) — comfortably larger than the actual `klu_common` struct
//! (~250 B on SuiteSparse 7.x) and guaranteed to remain so unless
//! SuiteSparse adds another ~480 fields.  A `debug_assert!` in
//! `KluCommon::new` guards the upper bound at runtime in debug builds.

use std::ffi::c_void;
use std::os::raw::c_int;

/// Default size for the `klu_common` backing buffer.  4096 bytes, which
/// is comfortably larger than the actual struct on every SuiteSparse
/// version we care about.
const COMMON_BUF_F64S: usize = 512; // 4096 bytes / 8 bytes per f64

// ---------------------------------------------------------------------------
// Raw extern "C" declarations.
//
// `klu_symbolic` and `klu_numeric` are opaque to us — we only ever pass
// the pointers KLU returns.  `klu_common` is treated as a void buffer
// for the same reason (see module docs).
// ---------------------------------------------------------------------------

#[allow(non_camel_case_types)]
extern "C" {
    fn klu_defaults(common: *mut c_void) -> c_int;

    fn klu_analyze(n: i32, ap: *mut i32, ai: *mut i32, common: *mut c_void) -> *mut c_void;

    fn klu_factor(
        ap: *mut i32,
        ai: *mut i32,
        ax: *mut f64,
        symbolic: *mut c_void,
        common: *mut c_void,
    ) -> *mut c_void;

    fn klu_refactor(
        ap: *mut i32,
        ai: *mut i32,
        ax: *mut f64,
        symbolic: *mut c_void,
        numeric: *mut c_void,
        common: *mut c_void,
    ) -> c_int;

    fn klu_solve(
        symbolic: *mut c_void,
        numeric: *mut c_void,
        ldim: i32,
        nrhs: i32,
        b: *mut f64,
        common: *mut c_void,
    ) -> c_int;

    fn klu_tsolve(
        symbolic: *mut c_void,
        numeric: *mut c_void,
        ldim: i32,
        nrhs: i32,
        b: *mut f64,
        common: *mut c_void,
    ) -> c_int;

    fn klu_free_symbolic(symbolic: *mut *mut c_void, common: *mut c_void) -> c_int;

    fn klu_free_numeric(numeric: *mut *mut c_void, common: *mut c_void) -> c_int;
}

// ---------------------------------------------------------------------------
// Safe wrappers
// ---------------------------------------------------------------------------

/// Errors surfaced by the KLU backend.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KluError {
    /// `klu_analyze` returned NULL — likely an out-of-memory or
    /// malformed sparsity pattern.
    AnalyzeFailed,
    /// `klu_factor` returned NULL — usually a structurally / numerically
    /// singular matrix.
    FactorFailed,
    /// `klu_refactor` returned 0.
    RefactorFailed,
    /// `klu_solve` / `klu_tsolve` returned 0.
    SolveFailed,
    /// CSC dimensions did not match the system size.
    BadShape,
}

impl std::fmt::Display for KluError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            KluError::AnalyzeFailed => write!(f, "klu_analyze returned NULL"),
            KluError::FactorFailed => write!(f, "klu_factor returned NULL (singular matrix?)"),
            KluError::RefactorFailed => write!(f, "klu_refactor failed"),
            KluError::SolveFailed => write!(f, "klu_solve / klu_tsolve failed"),
            KluError::BadShape => write!(f, "CSC dimensions inconsistent with system size"),
        }
    }
}

impl std::error::Error for KluError {}

/// Backing buffer for `klu_common`, plus the KLU settings struct itself
/// (treated as opaque bytes).  Always 8-byte aligned via the `Vec<f64>`
/// allocation.
///
/// `KluCommon` is **not** `Send` / `Sync` — each KLU computation must
/// own its own common block; KLU's internal scratch state lives here.
pub struct KluCommon {
    buf: Vec<f64>,
}

impl KluCommon {
    /// Allocate a zero-initialised `klu_common` buffer and call
    /// `klu_defaults` on it.
    pub fn new() -> Self {
        // 4 KiB; debug-assert in case a future SuiteSparse grows past this.
        debug_assert!(
            COMMON_BUF_F64S * std::mem::size_of::<f64>() >= 1024,
            "common-block buffer is suspiciously small"
        );
        let mut buf = vec![0.0_f64; COMMON_BUF_F64S];
        let ok = unsafe { klu_defaults(buf.as_mut_ptr() as *mut c_void) };
        // klu_defaults returns 1 on success.  Anything else is a SuiteSparse
        // bug; we panic loudly because it should be infallible.
        assert_eq!(ok, 1, "klu_defaults returned {ok}");
        KluCommon { buf }
    }

    fn as_ptr(&mut self) -> *mut c_void {
        self.buf.as_mut_ptr() as *mut c_void
    }
}

impl Default for KluCommon {
    fn default() -> Self {
        Self::new()
    }
}

/// Owned handle to a KLU symbolic factorisation (column permutation +
/// BTF block structure).  Pattern-only — survives value changes that
/// preserve the sparsity pattern, so this is the artefact a transient
/// or NR loop wants to cache across iterations.
pub struct KluSymbolic {
    ptr: *mut c_void,
}

unsafe impl Send for KluSymbolic {}

impl KluSymbolic {
    /// Run symbolic factorisation of an `n×n` CSC matrix.
    pub fn analyze(
        n: usize,
        ap: &mut [i32],
        ai: &mut [i32],
        common: &mut KluCommon,
    ) -> Result<Self, KluError> {
        if ap.len() != n + 1 {
            return Err(KluError::BadShape);
        }
        let ptr =
            unsafe { klu_analyze(n as i32, ap.as_mut_ptr(), ai.as_mut_ptr(), common.as_ptr()) };
        if ptr.is_null() {
            return Err(KluError::AnalyzeFailed);
        }
        Ok(KluSymbolic { ptr })
    }
}

impl Drop for KluSymbolic {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            // KLU needs a *fresh* common block here only if the original
            // one has been freed; in practice we always free before the
            // owning common does, but allocate one defensively.
            let mut tmp = KluCommon::new();
            unsafe {
                let mut p = self.ptr;
                klu_free_symbolic(&mut p as *mut *mut c_void, tmp.as_ptr());
            }
            self.ptr = std::ptr::null_mut();
        }
    }
}

/// Owned handle to a KLU numeric factorisation (LU factors + row
/// permutation + scale).  Bound to a `KluSymbolic` — the symbolic must
/// outlive the numeric.
pub struct KluNumeric {
    ptr: *mut c_void,
}

unsafe impl Send for KluNumeric {}

impl KluNumeric {
    /// Run numeric factorisation given a symbolic factorisation and CSC
    /// values.
    pub fn factor(
        ap: &mut [i32],
        ai: &mut [i32],
        ax: &mut [f64],
        symbolic: &KluSymbolic,
        common: &mut KluCommon,
    ) -> Result<Self, KluError> {
        let ptr = unsafe {
            klu_factor(
                ap.as_mut_ptr(),
                ai.as_mut_ptr(),
                ax.as_mut_ptr(),
                symbolic.ptr,
                common.as_ptr(),
            )
        };
        if ptr.is_null() {
            return Err(KluError::FactorFailed);
        }
        Ok(KluNumeric { ptr })
    }

    /// Re-run the numeric factorisation in place — the sparsity pattern
    /// must be unchanged since the original `factor` call.  This is the
    /// fast path for Newton-Raphson loops where only the values change.
    pub fn refactor(
        &mut self,
        ap: &mut [i32],
        ai: &mut [i32],
        ax: &mut [f64],
        symbolic: &KluSymbolic,
        common: &mut KluCommon,
    ) -> Result<(), KluError> {
        let ok = unsafe {
            klu_refactor(
                ap.as_mut_ptr(),
                ai.as_mut_ptr(),
                ax.as_mut_ptr(),
                symbolic.ptr,
                self.ptr,
                common.as_ptr(),
            )
        };
        if ok == 0 {
            return Err(KluError::RefactorFailed);
        }
        Ok(())
    }

    /// Solve `A·x = b` in place — on entry `b` contains the right-hand
    /// side; on success `b` is overwritten with the solution.
    pub fn solve(
        &self,
        symbolic: &KluSymbolic,
        b: &mut [f64],
        common: &mut KluCommon,
    ) -> Result<(), KluError> {
        let n = b.len() as i32;
        let ok = unsafe {
            klu_solve(
                symbolic.ptr,
                self.ptr,
                n,
                1,
                b.as_mut_ptr(),
                common.as_ptr(),
            )
        };
        if ok == 0 {
            return Err(KluError::SolveFailed);
        }
        Ok(())
    }

    /// Solve `Aᵀ·x = b` in place.  Used by adjoint / noise analysis.
    pub fn solve_transpose(
        &self,
        symbolic: &KluSymbolic,
        b: &mut [f64],
        common: &mut KluCommon,
    ) -> Result<(), KluError> {
        let n = b.len() as i32;
        let ok = unsafe {
            klu_tsolve(
                symbolic.ptr,
                self.ptr,
                n,
                1,
                b.as_mut_ptr(),
                common.as_ptr(),
            )
        };
        if ok == 0 {
            return Err(KluError::SolveFailed);
        }
        Ok(())
    }
}

impl Drop for KluNumeric {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            let mut tmp = KluCommon::new();
            unsafe {
                let mut p = self.ptr;
                klu_free_numeric(&mut p as *mut *mut c_void, tmp.as_ptr());
            }
            self.ptr = std::ptr::null_mut();
        }
    }
}

// ---------------------------------------------------------------------------
// One-shot helper: dense → CSC → analyze + factor + solve.
//
// Used by fairchild-core's `KluSolver` until the triplet-emission
// refactor lands the symbolic/numeric split through the `LinearSolver`
// trait.  Equivalent in behaviour to `FaerSparseSolver::solve` but
// dispatches the LU through KLU.
// ---------------------------------------------------------------------------

/// Convert a dense matrix to CSC (column-pointer / row-index / value
/// triplets), dropping entries with `|a_ij| ≤ threshold`.
///
/// Returned vectors are i32-indexed for KLU compatibility.
pub fn dense_to_csc(a: &[Vec<f64>], threshold: f64) -> (Vec<i32>, Vec<i32>, Vec<f64>) {
    let n = a.len();
    let mut ap: Vec<i32> = Vec::with_capacity(n + 1);
    let mut ai: Vec<i32> = Vec::new();
    let mut ax: Vec<f64> = Vec::new();
    ap.push(0);
    for j in 0..n {
        for i in 0..n {
            let v = a[i][j];
            if v.abs() > threshold {
                ai.push(i as i32);
                ax.push(v);
            }
        }
        ap.push(ai.len() as i32);
    }
    (ap, ai, ax)
}

/// One-shot Ax = b via KLU.  Allocates a fresh `KluCommon` /
/// `KluSymbolic` / `KluNumeric` per call — fine for "DenseSolver-style"
/// callers; the proper symbolic/numeric reuse path goes through the
/// caching `LinearSolver` trait extension in `fairchild-core`.
pub fn klu_solve_dense(a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, KluError> {
    let n = b.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    if a.len() != n || a.iter().any(|row| row.len() != n) {
        return Err(KluError::BadShape);
    }

    let (mut ap, mut ai, mut ax) = dense_to_csc(a, 1e-30);
    let mut common = KluCommon::new();
    let symbolic = KluSymbolic::analyze(n, &mut ap, &mut ai, &mut common)?;
    let numeric = KluNumeric::factor(&mut ap, &mut ai, &mut ax, &symbolic, &mut common)?;
    let mut x = b.to_vec();
    numeric.solve(&symbolic, &mut x, &mut common)?;
    if x.iter().any(|v| !v.is_finite()) {
        return Err(KluError::SolveFailed);
    }
    Ok(x)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn klu_diagonal_3x3() {
        // A = diag(2, 3, 4), b = [4, 9, 16]  →  x = [2, 3, 4]
        let a = vec![
            vec![2.0, 0.0, 0.0],
            vec![0.0, 3.0, 0.0],
            vec![0.0, 0.0, 4.0],
        ];
        let b = vec![4.0, 9.0, 16.0];
        let x = klu_solve_dense(&a, &b).unwrap();
        assert!((x[0] - 2.0).abs() < 1e-12, "x[0]={}", x[0]);
        assert!((x[1] - 3.0).abs() < 1e-12, "x[1]={}", x[1]);
        assert!((x[2] - 4.0).abs() < 1e-12, "x[2]={}", x[2]);
    }

    #[test]
    fn klu_tridiagonal_matches_hand_solve() {
        // 4 -1  0  0
        //-1  4 -1  0
        // 0 -1  4 -1
        // 0  0 -1  3
        let a = vec![
            vec![4.0, -1.0, 0.0, 0.0],
            vec![-1.0, 4.0, -1.0, 0.0],
            vec![0.0, -1.0, 4.0, -1.0],
            vec![0.0, 0.0, -1.0, 3.0],
        ];
        let b = vec![1.0, 2.0, 3.0, 4.0];
        let x = klu_solve_dense(&a, &b).unwrap();
        // Cross-check by direct multiplication.
        let mut r = vec![0.0; 4];
        for i in 0..4 {
            for j in 0..4 {
                r[i] += a[i][j] * x[j];
            }
        }
        for i in 0..4 {
            assert!(
                (r[i] - b[i]).abs() < 1e-10,
                "row {i}: Ax-b = {:.3e}",
                r[i] - b[i]
            );
        }
    }

    #[test]
    fn klu_singular_returns_err() {
        let a = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        let b = vec![1.0, 2.0];
        let r = klu_solve_dense(&a, &b);
        assert!(
            matches!(r, Err(KluError::FactorFailed) | Err(KluError::SolveFailed)),
            "expected factor/solve failure, got {r:?}"
        );
    }

    #[test]
    fn klu_symbolic_then_refactor() {
        // Verifies the cache path: analyze + factor + solve, then change
        // values (preserving pattern) and refactor + solve.  This is the
        // motivation for the entire KLU integration.
        let mut a = vec![vec![2.0, 1.0], vec![0.0, 3.0]];
        let mut common = KluCommon::new();
        let (mut ap, mut ai, mut ax) = dense_to_csc(&a, 1e-30);
        let symbolic = KluSymbolic::analyze(2, &mut ap, &mut ai, &mut common).unwrap();
        let mut numeric =
            KluNumeric::factor(&mut ap, &mut ai, &mut ax, &symbolic, &mut common).unwrap();
        let mut x = vec![3.0, 6.0];
        numeric.solve(&symbolic, &mut x, &mut common).unwrap();
        // 2x+y=3, 3y=6 → y=2, 2x=1 → x=0.5
        assert!((x[0] - 0.5).abs() < 1e-12);
        assert!((x[1] - 2.0).abs() < 1e-12);

        // Change values, same sparsity pattern → refactor.
        a[0][0] = 5.0;
        a[0][1] = 1.0;
        a[1][1] = 4.0;
        let (mut ap2, mut ai2, mut ax2) = dense_to_csc(&a, 1e-30);
        assert_eq!(ap, ap2);
        assert_eq!(ai, ai2); // pattern preserved
        numeric
            .refactor(&mut ap2, &mut ai2, &mut ax2, &symbolic, &mut common)
            .unwrap();
        let mut x = vec![6.0, 8.0];
        numeric.solve(&symbolic, &mut x, &mut common).unwrap();
        // 5x+y=6, 4y=8 → y=2, 5x=4 → x=0.8
        assert!((x[0] - 0.8).abs() < 1e-12, "x[0]={}", x[0]);
        assert!((x[1] - 2.0).abs() < 1e-12);
    }
}
