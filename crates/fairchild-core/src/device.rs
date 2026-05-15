use crate::mna::MnaMatrix;

pub const K_BOLTZMANN: f64 = 1.380649e-23;
pub const Q_ELECTRON: f64 = 1.602176634e-19;

/// Index of a terminal in the MNA solution vector; `None` → ground (excluded from matrix).
pub type NodeId = Option<usize>;

/// Simulator context passed to device model callbacks at every eval.
pub struct SimContext {
    pub temperature: f64,   // Kelvin; default 300.15 K (27 °C, SPICE TNOM)
    pub omega_0: f64,        // rad/s carrier frequency for optical ports; 0 for electrical
    /// When true, device models apply junction-step limiters (pnjlim, fetlim).
    /// Mapped from `SimOptions::pnjlim` / `.options nopnjlim`.
    pub jlim_enabled: bool,
}

impl Default for SimContext {
    fn default() -> Self {
        SimContext { temperature: 300.15, omega_0: 0.0, jlim_enabled: true }
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
    /// Compute reactive (capacitive/inductive) contributions — dq/dt and dq/dV.
    pub transient: bool,
}

impl EvalFlags {
    pub fn dc() -> Self {
        EvalFlags { resistive: true, transient: false }
    }

    pub fn tran() -> Self {
        EvalFlags { resistive: true, transient: true }
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

    /// Transient residual: stamp reactive + resistive contributions for BE companion.
    /// `alpha` = h_new/h_old (1.0 for fixed-step BE).
    /// Default falls back to DC residual (correct for purely resistive devices).
    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) { self.load_residual(b); }

    /// Transient Jacobian: stamp resistive + α·reactive conductances.
    /// Default falls back to DC Jacobian.
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) { self.load_jacobian(mat); }

    /// Record the current solution as the previous-timestep reference for reactive terms.
    ///
    /// Called once per accepted timestep in the variable-step integrator.
    /// Default is a no-op (correct for resistive-only and built-in devices).
    fn commit_timestep(&mut self, _x: &[f64]) {}

    /// Set a named real-valued parameter on this device instance.
    ///
    /// Returns `true` if the parameter was found and set; `false` if not supported
    /// or if `name` is not a parameter of this device (built-in devices return false
    /// by default — their parameters are set at construction time).
    fn set_real_param(&mut self, _name: &str, _value: f64) -> bool { false }

    /// Internal noise current sources for `.noise` analysis.
    ///
    /// Called AFTER a DC-OP `eval()` so the device can read its cached
    /// bias-point quantities (e.g. Id, gm) when building the PSDs.
    ///
    /// Returns one or more `(pos, neg, S_i)` tuples where `S_i` is the
    /// one-sided current PSD in A²/Hz of an uncorrelated noise source
    /// connected between those two terminal node indices.  Default is empty
    /// (resistor thermal noise is iterated as `Element::Resistor` in
    /// `noise_analysis`, not through this hook).
    fn noise_sources(&self, _ctx: &SimContext) -> Vec<(NodeId, NodeId, f64)> {
        Vec::new()
    }

    /// How many extra MNA rows this device needs beyond its terminal nodes.
    ///
    /// Verilog-A models with potential contributions (`V(port) <+ value;`)
    /// implicitly introduce branch-flow nodes that OSDI exposes as
    /// `num_nodes − num_terminals`.  The topology allocates that many rows
    /// per device and then calls `bind_extra_nodes` with the starting index.
    /// Returns 0 by default — native Rust devices that stamp directly into
    /// terminal nodes have no extras.
    fn num_extra_nodes(&self) -> usize { 0 }

    /// Tell the device which contiguous MNA rows (starting at `first_idx`)
    /// have been allocated for its internal nodes.  Called once after
    /// construction and before the first `eval()`.  Default is a no-op.
    fn bind_extra_nodes(&mut self, _first_idx: usize) {}
}
