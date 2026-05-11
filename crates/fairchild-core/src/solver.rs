use faer::{Col, Mat, linalg::solvers::Solve};

use crate::error::SimError;

/// Solve A·x = b using partial-pivot LU (faer dense).
/// Returns Err(SingularMatrix) if the result contains NaN/Inf.
pub fn lu_solve(a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SimError> {
    let n = b.len();
    if n == 0 {
        return Ok(vec![]);
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
