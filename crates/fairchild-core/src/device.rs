use crate::mna::MnaMatrix;

pub const K_BOLTZMANN: f64 = 1.380649e-23;
pub const Q_ELECTRON: f64 = 1.602176634e-19;

/// Index of a terminal in the MNA solution vector; `None` → ground (excluded from matrix).
pub type NodeId = Option<usize>;

/// Simulator context passed to device model callbacks at every eval.
pub struct SimContext {
    pub temperature: f64,   // Kelvin; default 300.15 K (27 °C, SPICE TNOM)
    pub omega_0: f64,        // rad/s carrier frequency for optical ports; 0 for electrical
}

impl Default for SimContext {
    fn default() -> Self {
        SimContext { temperature: 300.15, omega_0: 0.0 }
    }
}

impl SimContext {
    /// Thermal voltage kT/q in volts.
    pub fn vt(&self) -> f64 {
        K_BOLTZMANN * self.temperature / Q_ELECTRON
    }
}

/// Bit-flags controlling which contributions `eval` should compute.
pub struct EvalFlags {
    /// Compute resistive (DC / quasi-static) contributions — I(V) and dI/dV.
    pub resistive: bool,
}

impl EvalFlags {
    pub fn dc() -> Self {
        EvalFlags { resistive: true }
    }
}

/// Unified interface for circuit device models — both OSDI-loaded and Rust-native.
///
/// Call order per Newton-Raphson iteration:
///   1. `eval`           — device evaluates physics at current `x`, caches Norton equivalent
///   2. `load_residual`  — stamps cached currents into `b`
///   3. `load_jacobian`  — stamps cached conductances into `mat.a`
pub trait Device: Send + Sync {
    fn num_terminals(&self) -> usize;

    /// Finalise model-level derived quantities (temperature-scaled params, etc.).
    /// Called once after construction before any `setup_instance`.
    fn setup_model(&mut self, ctx: &SimContext);

    /// Bind this instance to specific MNA node indices.
    /// `terminals` length must equal `num_terminals()`; `None` = connected to ground.
    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext);

    /// Evaluate device physics at operating point `x`.
    /// Results are cached internally for subsequent `load_residual` / `load_jacobian` calls.
    fn eval(&mut self, x: &[f64], flags: EvalFlags, ctx: &SimContext);

    /// Add the Norton-equivalent current source contribution to the residual vector `b`.
    fn load_residual(&self, b: &mut [f64]);

    /// Add the Norton-equivalent conductance contribution to the MNA Jacobian `mat`.
    fn load_jacobian(&self, mat: &mut MnaMatrix);
}
