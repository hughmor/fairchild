use faer::{Col, Mat, linalg::solvers::Solve};

use crate::error::SimError;
use crate::mna::MnaSystem;

/// Solve the MNA system A*x = b and return the solution vector.
pub fn solve_dc(sys: &MnaSystem) -> Result<Vec<f64>, SimError> {
    let n = sys.size;
    if n == 0 {
        return Ok(vec![]);
    }

    let a_mat = Mat::<f64>::from_fn(n, n, |i, j| sys.a[i][j]);
    let b_col = Col::<f64>::from_fn(n, |i| sys.b[i]);

    let lu = a_mat.partial_piv_lu();

    // Check for NaN/Inf which indicate a singular or ill-conditioned system.
    let x_col = lu.solve(b_col.as_ref());
    let x: Vec<f64> = (0..n).map(|i| x_col[i]).collect();
    if x.iter().any(|v| !v.is_finite()) {
        return Err(SimError::SingularMatrix);
    }

    Ok(x)
}

/// Convenience wrapper: build MNA from netlist, solve, return OpResult.
pub fn dc_op(netlist: &fairchild_parser::Netlist) -> Result<OpResult, SimError> {
    let sys = MnaSystem::build(netlist)?;
    let x = solve_dc(&sys)?;
    Ok(OpResult { sys, x })
}

/// DC operating-point result.
pub struct OpResult {
    sys: MnaSystem,
    x: Vec<f64>,
}

impl OpResult {
    pub fn node_voltage(&self, node: &str) -> Result<f64, SimError> {
        self.sys.node_voltage(node, &self.x)
    }

    pub fn vsrc_current(&self, name: &str) -> Result<f64, SimError> {
        self.sys.vsrc_current(name, &self.x)
    }

    pub fn all_voltages(&self) -> impl Iterator<Item = (&str, f64)> {
        self.sys.node_index.iter().map(|(name, &i)| (name.as_str(), self.x[i]))
    }
}
