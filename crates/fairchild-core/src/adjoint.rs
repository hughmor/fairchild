//! Parameter sensitivity at a DC operating point, by the adjoint method.
//!
//! Given a converged operating point `x*` solving `f(x, p) = 0` and a scalar
//! output `L(x)`, the total derivative with respect to a design parameter `p` is
//!
//! ```text
//!     dL/dp = −λᵀ · ∂f/∂p        where     Jᵀ · λ = (∂L/∂x)ᵀ
//! ```
//!
//! and `J = ∂f/∂x` is the Newton Jacobian — the matrix already stamped and
//! already factorised at the last NR iteration.  So the cost of a gradient is
//! **one transposed back-substitution per output**, plus one residual
//! re-evaluation per parameter.  Neither is a nonlinear solve.  Sweeping 50
//! parameters costs one simulation, not 50.
//!
//! Two things make this worth having over "just finite-difference the whole
//! simulation":
//!
//! 1. **Cost.** Forward differences pay a full Newton solve per parameter.  The
//!    adjoint pays a back-substitution per *output*, and outputs are usually far
//!    fewer than parameters — that asymmetry is the entire point of running the
//!    system backwards.
//! 2. **Accuracy, which matters more than the cost.** Differencing a converged
//!    solution differences a quantity that is only accurate to `reltol`
//!    (1e-3 by default).  With a 1e-6 relative step that leaves roughly three
//!    significant figures of *noise* in the answer, and shrinking the step makes
//!    it worse.  The adjoint differences `f`, which is an explicit function of
//!    `(x, p)` evaluated to machine precision, so the usual `∛ε` step analysis
//!    applies and the result is good to ~1e-10 relative.
//!
//! ## What is differentiable, and what is not
//!
//! `∂f/∂p` is taken by central difference **on the residual only** — perturb the
//! parameter, re-stamp at the frozen `x*`, subtract.  That is deliberate: it
//! works for every device in the tree the day it is written, including OSDI /
//! Verilog-A models and every photonic device, without a per-device derivative
//! to write or keep correct.  For any parameter that enters the stamp linearly
//! (an `R`, `C`, `L`, a source level) the central difference of a linear
//! function is exact, so "finite difference" costs nothing there at all.
//!
//! Two paths reach a parameter, chosen by element kind:
//!
//! * `R`/`C`/`L`/`V`/`I` are stamped straight from the netlist, so the netlist
//!   copy is edited and the element stamp re-run.  No device is touched, so no
//!   device state moves.
//! * Everything else is a `Device`, and the *live, converged* instance is
//!   retuned via [`Device::set_real_param`].  Rebuilding it instead would reset
//!   the junction limiter's `v_prev` to zero, and a fresh diode re-evaluated at
//!   `x*` limits its first step — the residual would not be `f(x*)` at all.
//!
//! A parameter that reaches neither is reported as `reached = false` rather than
//! contributing a zero gradient.  A silent zero is indistinguishable from a real
//! insensitivity, and this is a numerical-optimisation API: a wrong zero stalls
//! the optimiser somewhere that looks like a stationary point.  Today that
//! covers diode and BJT model parameters (neither implements `set_real_param`;
//! see `docs/model_status.md`) and MOSFET geometry.

use fairchild_parser::{Element, Netlist};

use crate::device::Device;
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::mna::{CircuitTopology, MnaMatrix, StampPlan};
use crate::newton::{
    build_devices_with_footprints, dc_op_nr_with_devices_opts, device_element_names, residual_l2,
};
use crate::options::SimOptions;

/// A design parameter to differentiate with respect to, named as it is in the
/// netlist: the element's refdes plus a parameter name
/// [`crate::netlist_edit::set_element_param`] understands.
#[derive(Clone, Debug)]
pub struct ParamRef {
    pub element: String,
    pub param: String,
    /// Nominal value.  `None` reads it off the instance line; supply it
    /// explicitly for a parameter that only appears on a `.model` card, which
    /// the netlist getter cannot see.
    pub nominal: Option<f64>,
    /// Absolute finite-difference step.  `None` picks `∛ε·|p|`, which assumes
    /// the residual varies on the parameter's own scale.
    ///
    /// It does not always.  An optical length enters the residual through
    /// `φ = 2π·n_g·L/λ`, so `∂φ/∂L ≈ 17 rad/µm` — a step scaled to `L` moves the
    /// phase by radians, and the power gradient is a near-total cancellation
    /// between two much larger phase terms.  [`dc_sensitivity`] recovers from
    /// that by differencing at two step sizes and extrapolating, and reports
    /// what it had to remove as `Sensitivities::fd_error`; set this to pin the
    /// step instead.
    pub step: Option<f64>,
}

impl ParamRef {
    pub fn new(element: impl Into<String>, param: impl Into<String>) -> Self {
        ParamRef {
            element: element.into(),
            param: param.into(),
            nominal: None,
            step: None,
        }
    }

    /// [`ParamRef::new`] with the nominal value supplied rather than looked up.
    pub fn with_nominal(element: impl Into<String>, param: impl Into<String>, v: f64) -> Self {
        ParamRef {
            nominal: Some(v),
            ..ParamRef::new(element, param)
        }
    }

    /// Pin the finite-difference step rather than letting it be chosen.
    pub fn with_step(mut self, h: f64) -> Self {
        self.step = Some(h);
        self
    }
}

/// The scalar being differentiated.  Each one costs one transposed solve.
#[derive(Clone, Debug)]
pub enum Output {
    /// `v(node)`.
    NodeVoltage(String),
    /// `v(pos) − v(neg)` — a port voltage rather than a ground-referenced one.
    /// `.tf` and `.sens` both need it: a `.tf` output resistance is measured
    /// *across a port*, and the port is only the same thing as a node when the
    /// far side happens to be ground.
    NodeVoltageDiff { pos: String, neg: String },
    /// The branch current through a voltage source, SPICE sign convention.
    BranchCurrent(String),
    /// Optical power `re² + im²` on one channel of an optical net, in watts —
    /// the quantity a photodetector responds to and the natural objective for
    /// an optical design.  Reads the `<net>_re_<ch>` / `<net>_im_<ch>` wires
    /// that `.optical_port` creates.
    OpticalPower { net: String, channel: usize },
    /// An arbitrary `∂L/∂x`, length `topo.size`, for anything the variants
    /// above do not spell.  The reported value is `∂L/∂x · x`, which is `L`
    /// itself only when `L` is linear — supply the value separately if it
    /// is not.
    Custom(Vec<f64>),
}

impl Output {
    /// `(value at x*, ∂L/∂x at x*)`.
    pub(crate) fn seed(
        &self,
        topo: &CircuitTopology,
        x: &[f64],
    ) -> Result<(f64, Vec<f64>), SimError> {
        let mut seed = vec![0.0; topo.size];
        match self {
            Output::NodeVoltage(node) => {
                let value = topo.node_voltage(node, x)?;
                // Ground is a legitimate request with an identically zero
                // gradient; it has no row, so leave the seed at zero.
                if let Some(&i) = topo.node_index.get(node) {
                    seed[i] = 1.0;
                }
                Ok((value, seed))
            }
            Output::NodeVoltageDiff { pos, neg } => {
                let value = topo.node_voltage(pos, x)? - topo.node_voltage(neg, x)?;
                // Ground has no row, so it contributes nothing to the seed —
                // the same rule `NodeVoltage` follows, applied to both ends.
                if let Some(&i) = topo.node_index.get(pos) {
                    seed[i] += 1.0;
                }
                if let Some(&i) = topo.node_index.get(neg) {
                    seed[i] -= 1.0;
                }
                Ok((value, seed))
            }
            Output::BranchCurrent(vsrc) => {
                let value = topo.vsrc_current(vsrc, x)?;
                let i = *topo
                    .vsrc_index
                    .get(vsrc)
                    .ok_or_else(|| SimError::UnknownNode(vsrc.clone()))?;
                seed[topo.n_nodes() + i] = 1.0;
                Ok((value, seed))
            }
            Output::OpticalPower { net, channel } => {
                let re_name = format!("{net}_re_{channel}");
                let im_name = format!("{net}_im_{channel}");
                let re = topo.node_voltage(&re_name, x)?;
                let im = topo.node_voltage(&im_name, x)?;
                // ∂(re² + im²)/∂x — ground-referenced wires again contribute
                // nothing, and a dark port is exactly the zero-gradient case.
                if let Some(&i) = topo.node_index.get(&re_name) {
                    seed[i] = 2.0 * re;
                }
                if let Some(&i) = topo.node_index.get(&im_name) {
                    seed[i] = 2.0 * im;
                }
                Ok((re * re + im * im, seed))
            }
            Output::Custom(v) => {
                if v.len() != topo.size {
                    return Err(SimError::ParameterError(format!(
                        "Output::Custom has length {} but the system has {} unknowns",
                        v.len(),
                        topo.size
                    )));
                }
                let value = v.iter().zip(x.iter()).map(|(a, b)| a * b).sum();
                Ok((value, v.clone()))
            }
        }
    }
}

/// Result of [`dc_sensitivity`].
pub struct Sensitivities {
    /// The converged operating point.
    pub x: Vec<f64>,
    /// Each requested output, evaluated at `x`.
    pub values: Vec<f64>,
    /// `grad[output][param] = dL/dp`.
    pub grad: Vec<Vec<f64>>,
    /// Per parameter: did the perturbation actually move the residual?  `false`
    /// means the gradient entry is a placeholder zero, not a computed one.
    pub reached: Vec<bool>,
    /// Per parameter: the relative disagreement between the two finite-difference
    /// step sizes, before Richardson extrapolation removed it.  A conservative
    /// error bar on that column of the gradient — the extrapolated value is
    /// normally far better than this suggests, but a large value means the
    /// parameter's natural scale and the residual's disagree badly enough to be
    /// worth pinning `ParamRef::step` by hand.
    pub fd_error: Vec<f64>,
}

impl Sensitivities {
    /// The names of any parameters that could not be reached, for a caller that
    /// wants to fail loudly rather than optimise against placeholder zeros.
    pub fn unreached<'a>(&self, params: &'a [ParamRef]) -> Vec<&'a ParamRef> {
        params
            .iter()
            .zip(self.reached.iter())
            .filter(|(_, ok)| !**ok)
            .map(|(p, _)| p)
            .collect()
    }
}

/// Which mechanism carries a perturbation to the equations.
pub(crate) enum Handle {
    /// Index into `netlist.elements`; the element stamper reads its value each
    /// pass, so editing the netlist copy is enough.
    NetlistElement,
    /// Index into `devices`.
    Device(usize),
}

/// `∂(each output)/∂(each parameter)` at the DC operating point.
///
/// One nonlinear solve, then one transposed back-substitution per output and
/// two residual re-stamps per parameter.
pub fn dc_sensitivity(
    netlist: &Netlist,
    registry: &DeviceRegistry,
    opts: &SimOptions,
    outputs: &[Output],
    params: &[ParamRef],
) -> Result<Sensitivities, SimError> {
    let ctx = opts.sim_context();
    let mut topo = CircuitTopology::build_resolved(netlist, &ctx, registry);
    let (mut devices, footprints) =
        build_devices_with_footprints(netlist, &mut topo, &ctx, registry)?;
    let plan = StampPlan::new(&topo, netlist, &footprints);
    plan.resolve_device_cells(&mut devices);

    let op = dc_op_nr_with_devices_opts(netlist, &topo, &mut devices, &ctx, opts)?;
    let x = op.x;

    // Re-stamp at the converged point.  `residual_l2` leaves the matrix holding
    // that linearisation, which is exactly the Jacobian the adjoint needs — the
    // norm it returns is the convergence check, discarded here.
    let mut mat = MnaMatrix::with_pattern(topo.size, plan.pattern.clone());
    let mut work = netlist.clone();
    let _ = residual_l2(
        &mut mat,
        &topo,
        &work,
        &mut devices,
        &ctx,
        1.0,
        0.0,
        Some(&plan),
        &x,
    );

    // Repair the columns the devices freeze rather than differentiate.  This is
    // what makes the gradient a *total* derivative on an electro-optic path;
    // see `Device::frozen_jacobian_columns`.  It happens before the
    // factorisation and only touches the adjoint's copy of the matrix, so
    // Newton's iteration matrix is unaffected.
    let frozen = frozen_columns(&topo, &devices);
    let repaired = frozen
        .iter()
        .map(|&col| {
            let column =
                fd_jacobian_column(col, &mut mat, &topo, &work, &mut devices, &ctx, &plan, &x);
            (col, column)
        })
        .collect::<Vec<_>>();
    // The FD stamps above left `mat` linearised somewhere other than `x`.
    let _ = residual_l2(
        &mut mat,
        &topo,
        &work,
        &mut devices,
        &ctx,
        1.0,
        0.0,
        Some(&plan),
        &x,
    );
    for (col, column) in &repaired {
        for (row, v) in column.iter().enumerate() {
            // Only write cells that carry something, or that the stamp already
            // has an entry for — writing every row would densify the column.
            if *v != 0.0 || mat.a[row][*col] != 0.0 {
                mat.a[row][*col] = *v;
            }
        }
    }

    // The residual at the operating point itself.
    let mut f0 = vec![0.0; topo.size];
    mat.residual_into(&x, &mut f0);

    // One transposed solve per output, all sharing the factorisation.
    let solver = opts.linear_solver(topo.size);
    let mut fact = solver.factorise_mat(&mat)?;
    let mut values = Vec::with_capacity(outputs.len());
    let mut lambdas = Vec::with_capacity(outputs.len());
    for out in outputs {
        let (value, seed) = out.seed(&topo, &x)?;
        values.push(value);
        lambdas.push(fact.refactor_and_solve_transpose(&mat.a, &seed)?);
    }

    let dev_names = device_element_names(netlist);
    let mut grad = vec![vec![0.0; params.len()]; outputs.len()];
    let mut reached = vec![false; params.len()];
    let mut fd_error = vec![0.0; params.len()];
    let n = topo.size;
    let mut df = vec![0.0; n];
    let mut coarse = vec![0.0; n];
    let mut scratch = (vec![0.0; n], vec![0.0; n]);

    for (pi, p) in params.iter().enumerate() {
        let Some((handle, nominal)) = resolve(p, &work, &dev_names) else {
            continue;
        };
        let h = p.step.unwrap_or_else(|| default_step(nominal));

        // Two step sizes, then Richardson.  A central difference has an O(h²)
        // error term, so `(4·D(h/2) − D(h))/3` cancels it and leaves O(h⁴).
        // That is what rescues a parameter whose natural scale is not the
        // residual's: an optical length moves the propagation phase by
        // ~17 rad/µm, so the power gradient is a near-total cancellation
        // between two much larger phase terms and a step sized to the length
        // is far too coarse.  Their disagreement is also the honest error bar,
        // reported as `fd_error`.
        let mut sample = |value: f64, out: &mut [f64]| {
            apply(&handle, p, value, &mut work, &mut devices);
            let _ = residual_l2(
                &mut mat,
                &topo,
                &work,
                &mut devices,
                &ctx,
                1.0,
                0.0,
                Some(&plan),
                &x,
            );
            mat.residual_into(&x, out);
        };
        for (step, out) in [(h, &mut coarse), (0.5 * h, &mut df)] {
            sample(nominal + step, &mut scratch.0);
            sample(nominal - step, &mut scratch.1);
            let inv = 1.0 / (2.0 * step);
            for ((o, a), b) in out.iter_mut().zip(scratch.0.iter()).zip(scratch.1.iter()) {
                *o = (a - b) * inv;
            }
        }
        apply(&handle, p, nominal, &mut work, &mut devices);

        let mut moved = false;
        let (mut spread, mut mag) = (0.0_f64, 0.0_f64);
        for (d, c) in df.iter_mut().zip(coarse.iter()) {
            let fine = *d;
            spread = spread.max((fine - c).abs());
            *d = (4.0 * fine - c) / 3.0;
            mag = mag.max(d.abs());
            moved |= *d != 0.0;
        }
        fd_error[pi] = if mag > 0.0 { spread / mag } else { 0.0 };

        if !moved {
            // The perturbation changed nothing in the equations.  Either the
            // device ignored `set_real_param` or the name is wrong; both are
            // reportable, neither is a gradient.
            continue;
        }
        reached[pi] = true;

        for (oi, lambda) in lambdas.iter().enumerate() {
            grad[oi][pi] = -lambda
                .iter()
                .zip(df.iter())
                .map(|(l, d)| l * d)
                .sum::<f64>();
        }
    }

    // `f0` is the residual at the operating point; a converged solve leaves it
    // near zero and it plays no part in a central difference.  Kept as the
    // caller-visible sanity handle rather than dropped silently.
    debug_assert!(
        f0.iter().all(|v| v.is_finite()),
        "non-finite residual at the operating point"
    );

    Ok(Sensitivities {
        x,
        values,
        grad,
        reached,
        fd_error,
    })
}

/// Every column whose `∂f/∂x` has to be re-derived numerically for the adjoint.
///
/// Declared by the device that froze it, and only that.  λ columns used to be
/// added here from the netlist as well: every optical device froze λ, because
/// `∂φ/∂λ = φ/λ` is of order 1e9 per metre and differentiating against it does
/// not converge, and asking fourteen models to say so would have been fourteen
/// copies of one fact.  That rule also only ever matched a λ net named
/// `_wl_<k>` by an `.optical_port`, so it silently missed every PCell that
/// hand-wires its bundle.  Both problems are gone with the rows: λ is resolved
/// before the solve, no device reads a λ column, and there is no λ column.
pub(crate) fn frozen_columns(topo: &CircuitTopology, devices: &[Box<dyn Device>]) -> Vec<usize> {
    let mut cols: Vec<usize> = Vec::new();
    for dev in devices {
        cols.extend(dev.frozen_jacobian_columns());
    }
    cols.retain(|&c| c < topo.size);
    cols.sort_unstable();
    cols.dedup();
    cols
}

/// One column of `∂f/∂x` by central difference, leaving `mat` stamped wherever
/// the last probe put it — the caller re-stamps at `x` afterwards.
#[allow(clippy::too_many_arguments)]
fn fd_jacobian_column(
    col: usize,
    mat: &mut MnaMatrix,
    topo: &CircuitTopology,
    netlist: &Netlist,
    devices: &mut Vec<Box<dyn Device>>,
    ctx: &crate::device::SimContext,
    plan: &StampPlan,
    x: &[f64],
) -> Vec<f64> {
    let n = topo.size;
    // A λ wire is ~1.55e-6 and a node voltage is O(1); scale to whichever the
    // column actually holds, with a floor so a column sitting at zero still
    // gets probed.
    let h = 6.0554545e-6_f64 * x[col].abs().max(1e-3);
    let mut probe = x.to_vec();
    let mut f_plus = vec![0.0; n];
    let mut f_minus = vec![0.0; n];

    probe[col] = x[col] + h;
    let _ = residual_l2(
        mat,
        topo,
        netlist,
        devices,
        ctx,
        1.0,
        0.0,
        Some(plan),
        &probe,
    );
    mat.residual_into(&probe, &mut f_plus);

    probe[col] = x[col] - h;
    let _ = residual_l2(
        mat,
        topo,
        netlist,
        devices,
        ctx,
        1.0,
        0.0,
        Some(plan),
        &probe,
    );
    mat.residual_into(&probe, &mut f_minus);

    let scale = 1.0 / (2.0 * h);
    f_plus
        .iter()
        .zip(f_minus.iter())
        .map(|(p, m)| (p - m) * scale)
        .collect()
}

// ---------------------------------------------------------------------------
// Jacobian completeness
// ---------------------------------------------------------------------------

/// One entry of the stamped Jacobian that disagrees with `∂f/∂x`.
#[derive(Clone, Debug)]
pub struct JacobianMismatch {
    pub row: usize,
    pub col: usize,
    pub stamped: f64,
    pub numeric: f64,
    /// Whether this column is one the adjoint already knows to re-derive — a λ
    /// wire, or a column some device declared via
    /// [`Device::frozen_jacobian_columns`].  Those are expected and repaired.
    /// **An undeclared mismatch is a bug**: either a wrong stamp, or a device
    /// that freezes a coefficient without saying so, which makes every gradient
    /// through it silently zero.
    pub frozen: bool,
}

impl JacobianMismatch {
    /// Size of the disagreement, relative to the larger of the two entries.
    pub fn rel_err(&self) -> f64 {
        (self.stamped - self.numeric).abs()
            / self
                .stamped
                .abs()
                .max(self.numeric.abs())
                .max(f64::MIN_POSITIVE)
    }
}

/// Compare the Jacobian the devices stamp against `∂f/∂x` taken numerically.
///
/// **This is a precondition of the adjoint method, not a nicety.**  `dL/dp =
/// −λᵀ·∂f/∂p` is only the total derivative if `Jᵀλ = ∂L/∂x` was solved with the
/// *true* `J = ∂f/∂x`.  A device is free to converge Newton with an approximate
/// Jacobian — successive substitution on a frozen coefficient reaches the same
/// fixed point, just linearly instead of quadratically — and the forward answer
/// is right either way.  The gradient is not: every path through a missing block
/// silently contributes zero.
///
/// Costs `2n` re-stamps, so this is a test and diagnostic tool, not something to
/// run inside an optimiser.  `atol` filters entries too small to compare
/// meaningfully; `rtol` is the acceptance threshold on the rest.
pub fn jacobian_check(
    netlist: &Netlist,
    registry: &DeviceRegistry,
    opts: &SimOptions,
    x: &[f64],
    rtol: f64,
    atol: f64,
) -> Result<Vec<JacobianMismatch>, SimError> {
    let ctx = opts.sim_context();
    let mut topo = CircuitTopology::build_resolved(netlist, &ctx, registry);
    let (mut devices, footprints) =
        build_devices_with_footprints(netlist, &mut topo, &ctx, registry)?;
    let plan = StampPlan::new(&topo, netlist, &footprints);
    plan.resolve_device_cells(&mut devices);
    let n = topo.size;
    if x.len() != n {
        return Err(SimError::ParameterError(format!(
            "jacobian_check: x has length {} but the system has {n} unknowns",
            x.len()
        )));
    }

    let mut mat = MnaMatrix::with_pattern(n, plan.pattern.clone());
    let frozen = frozen_columns(&topo, &devices);
    let stamp = |mat: &mut MnaMatrix, devices: &mut Vec<Box<dyn Device>>, at: &[f64]| {
        let _ = residual_l2(
            mat,
            &topo,
            netlist,
            devices,
            &ctx,
            1.0,
            0.0,
            Some(&plan),
            at,
        );
    };

    // Settle the devices at `x` before snapshotting anything.  A freshly built
    // device carries no history, so its first eval limits: a MOSFET handed
    // V_gs = 2 V with `v_prev = 0` stamps a `g_m` for the limited operating
    // point, not the real one, and the snapshot would disagree with `∂f/∂x` for
    // reasons that have nothing to do with a frozen coefficient.  Limiting is a
    // contraction, so re-evaluating at a fixed `x` converges in a few passes.
    let mut settle = vec![0.0; n];
    let mut previous = vec![f64::INFINITY; n];
    for _ in 0..32 {
        stamp(&mut mat, &mut devices, x);
        mat.residual_into(x, &mut settle);
        if settle
            .iter()
            .zip(previous.iter())
            .all(|(a, b)| (a - b).abs() <= 1e-15 * a.abs().max(1.0))
        {
            break;
        }
        previous.copy_from_slice(&settle);
    }
    let stamped = crate::mna::CircuitTopology::to_dense(&mat.a, n);

    let mut f_plus = vec![0.0; n];
    let mut f_minus = vec![0.0; n];
    let mut probe = x.to_vec();
    let mut out = Vec::new();

    for col in 0..n {
        let h = 6.0554545e-6_f64 * x[col].abs().max(1e-3);
        probe[col] = x[col] + h;
        stamp(&mut mat, &mut devices, &probe);
        mat.residual_into(&probe, &mut f_plus);
        probe[col] = x[col] - h;
        stamp(&mut mat, &mut devices, &probe);
        mat.residual_into(&probe, &mut f_minus);
        probe[col] = x[col];

        for row in 0..n {
            let numeric = (f_plus[row] - f_minus[row]) / (2.0 * h);
            let s = stamped[row][col];
            if s.abs().max(numeric.abs()) < atol {
                continue;
            }
            let m = JacobianMismatch {
                row,
                col,
                stamped: s,
                numeric,
                frozen: frozen.contains(&col),
            };
            if m.rel_err() > rtol {
                out.push(m);
            }
        }
    }
    Ok(out)
}

/// `∛ε·|p|` — the central-difference optimum when the residual varies on the
/// parameter's own scale, balancing `O(h²)` truncation against `O(ε/h)` roundoff.
///
/// ponytail: the zero-nominal fallback is a guess.  A parameter that is
/// genuinely zero has no scale of its own to borrow; pass `step` explicitly if
/// 1e-9 is wrong for yours.
pub(crate) fn default_step(nominal: f64) -> f64 {
    if nominal == 0.0 {
        1e-9
    } else {
        6.0554545e-6 * nominal.abs()
    }
}

/// Find how to reach `p`, and its nominal value.
pub(crate) fn resolve(
    p: &ParamRef,
    netlist: &Netlist,
    dev_names: &[String],
) -> Option<(Handle, f64)> {
    let el_lc = p.element.to_lowercase();
    let nominal = match p.nominal {
        Some(v) => v,
        None => crate::netlist_edit::get_element_param(netlist, &p.element, &p.param)?,
    };

    let is_netlist_stamped = netlist.elements.iter().any(|el| {
        matches!(
            el,
            Element::Resistor { .. }
                | Element::Capacitor { .. }
                | Element::Inductor { .. }
                | Element::VoltageSource { .. }
                | Element::CurrentSource { .. }
        ) && element_name(el).is_some_and(|n| n.to_lowercase() == el_lc)
    });
    if is_netlist_stamped {
        return Some((Handle::NetlistElement, nominal));
    }

    let idx = dev_names
        .iter()
        .position(|n| n.to_lowercase() == el_lc)
        .map(Handle::Device)?;
    Some((idx, nominal))
}

pub(crate) fn apply(
    handle: &Handle,
    p: &ParamRef,
    value: f64,
    netlist: &mut Netlist,
    devices: &mut [Box<dyn Device>],
) {
    match handle {
        Handle::NetlistElement => {
            crate::netlist_edit::set_element_param(netlist, &p.element, &p.param, value);
        }
        Handle::Device(i) => {
            devices[*i].set_real_param(&p.param, value);
        }
    }
}

/// The name of an element the netlist stamps directly — the set
/// [`crate::netlist_edit::set_element_param`] can retune without touching a
/// device.  `None` for everything else, which is what makes it the right
/// collision check for a synthetic probe element.
pub(crate) fn element_name(el: &Element) -> Option<&str> {
    match el {
        Element::Resistor { name, .. }
        | Element::Capacitor { name, .. }
        | Element::Inductor { name, .. }
        | Element::VoltageSource { name, .. }
        | Element::CurrentSource { name, .. } => Some(name),
        _ => None,
    }
}
