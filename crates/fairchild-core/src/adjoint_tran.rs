//! Parameter sensitivity of a **transient** run, by the discrete adjoint.
//!
//! [`crate::adjoint`] differentiates one operating point.  This differentiates a
//! whole waveform: given a scalar built from the trajectory — a value at one
//! timepoint, an integral, a mismatch against a target — it returns `dL/dp` for
//! every design parameter, at a cost that does not grow with the number of
//! *outputs* being traded off, and that never differentiates a converged solve.
//!
//! # The recursion
//!
//! Each timestep is a constraint.  Writing the resistive part as `F` and the
//! branch current the integrator produces as `i_k`, step `k` solves
//!
//! ```text
//!     G_k = F(x_k, p) + i_k = 0,      i_k = α_k·(q(x_k) − q(x_{k−1})) + β_k·i_{k−1}
//! ```
//!
//! with `(α, β)` reading straight off [`crate::reactive::charge_current`]:
//! Backward Euler is `α = 1/h, β = 0`, Trapezoidal is `α = 2/h, β = −1`.  That
//! `β` is the whole difficulty — under TR the current at step `k` depends on
//! *every* earlier step, so `∂G_j/∂x_k` is non-zero for all `j ≥ k` and a naive
//! backward pass is `O(n²)` in the number of timesteps.
//!
//! Carrying one extra co-state collapses it.  With `ū_k` the objective's
//! sensitivity to the stored current `i_k`,
//!
//! ```text
//!     ū_k = λ_k + β_{k+1}·ū_{k+1}
//!     A_kᵀ·λ_k = ∂L/∂x_k + (α_{k+1} − β_{k+1}·α_k)·J_q(x_k)ᵀ·ū_{k+1}
//! ```
//!
//! and the left-hand side comes out as exactly `A_k = ∂F/∂x_k + α_k·J_q(x_k)`,
//! which is **the matrix the forward solve already stamped**.  So the backward
//! pass is one transposed solve per timestep against a matrix that already
//! exists, and the whole method needs just one thing the forward pass does not
//! already produce: the charge Jacobian `J_q = ∂q/∂x`.
//!
//! The `t = 0` operating point is a constraint like any other, with `α_0 = 0`
//! and the DC Jacobian in place of `A_0` — so a transient gradient contains a DC
//! adjoint at its tail.  Under `.tran ... UIC` the initial state is imposed
//! rather than solved, and that term drops out.
//!
//! # Getting `J_q` without a per-device derivative
//!
//! `α` is the *only* way `h` enters a stamped transient Jacobian, and every
//! charge term enters through `α`.  So stamping the same point twice at two step
//! sizes and differencing isolates the charge block exactly:
//!
//! ```text
//!     J_q = (A(h) − A(2h)) / (α(h) − α(2h))
//! ```
//!
//! No new `Device` method, no per-model derivative to keep correct, and it
//! covers OSDI / Verilog-A `ddt` charges and bias-dependent junction caps the
//! day it is written — the same reasoning that makes `∂f/∂p` a residual
//! difference in [`crate::adjoint`].
//!
//! # `∂G_k/∂p` is a replay, not a re-solve
//!
//! A parameter does not only enter step `k`'s stamp: it enters the *history*
//! that step `k` integrates from.  Perturbing `C` changes `q(x_{k−1})` even
//! though `x_{k−1}` is held fixed, and under TR it changes the stored `i_{k−1}`
//! as well.  Rather than enumerate those paths — which is where a hand-derived
//! transient adjoint goes wrong — this re-walks the whole run with the parameter
//! nudged and **every `x_k` frozen at its solved value**, accumulating
//! `Σ_k λ_kᵀ·G_k` as it goes.  Nothing is solved; each step is one stamp.  Every
//! history path is then included by construction, because the history is
//! propagated by the same code that propagated it forward.
//!
//! The honest cost, therefore, is `O(n_params)` **restamps** — not the textbook
//! "independent of parameter count".  What the adjoint still buys is the far
//! more valuable half: the accuracy.  Differencing a converged transient
//! differences a quantity known only to `reltol`, and the error compounds
//! step over step; differencing the residual at a frozen trajectory does not.
//!
//! # Not covered yet
//!
//! * **Inductors**, netlist or device-declared.  A Norton-companion inductor
//!   scales as `1/α` rather than `α` and carries its flux as hidden state, so it
//!   fits neither the `J_q` extraction nor the co-state recursion above.  It is
//!   rejected loudly rather than silently mis-differentiated.
//! * **Variable step.**  This drives the fixed-step integrator
//!   ([`TranStepper`]), so `gear2_h_prev` is always absent and BDF-2 demotes to
//!   BE exactly as it does for a plain `.tran`.

use fairchild_parser::Netlist;

use crate::adjoint::{apply, default_step, frozen_columns, resolve, Handle, Output, ParamRef};
use crate::device::ReactiveKind;
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::mna::{CircuitTopology, SparseRow};
use crate::options::SimOptions;
use crate::reactive::conductance;
use crate::tran::IntegratorMode;
use crate::tran_step::TranStepper;

/// A transient run, plus everything the backward pass needs to differentiate it.
///
/// Built by [`TranAdjoint::run`]; ask it for gradients with
/// [`TranAdjoint::gradient`], as many times and against as many objectives as
/// you like — the expensive part is the forward pass, and it is already done.
pub struct TranAdjoint {
    netlist: Netlist,
    opts: SimOptions,
    step: f64,
    topo: CircuitTopology,
    /// Accepted timepoints, `t[0] = 0`.
    t: Vec<f64>,
    /// The solution at each timepoint.
    x: Vec<Vec<f64>>,
    /// `A_k = ∂G_k/∂x_k`.  `a[0]` is the DC Jacobian, or empty under UIC.
    a: Vec<Vec<SparseRow>>,
    /// `J_q(x_k) = ∂q/∂x` at each timepoint.
    jq: Vec<Vec<SparseRow>>,
    /// `α_k` from the integrator; `alpha[0] = 0` — the initial condition
    /// integrates nothing.
    alpha: Vec<f64>,
    /// `β_k`: `−1` where step `k` integrated with Trapezoidal, `0` otherwise.
    beta: Vec<f64>,
    /// The method each step actually used.  Trapezoidal takes its first step
    /// with Backward Euler, and the adjoint has to agree with what happened.
    mode: Vec<IntegratorMode>,
    /// Whether `t = 0` was solved (a constraint that carries `p`) or imposed
    /// by `.ic` under UIC (one that does not).
    dc_start: bool,
    /// Columns whose stamped `∂f/∂x` is frozen and has to be re-derived — see
    /// [`crate::device::Device::frozen_jacobian_columns`].
    frozen: Vec<usize>,
}

/// Result of [`TranAdjoint::gradient`].
pub struct TranSensitivities {
    /// `dL/dp`, one entry per requested parameter.
    pub grad: Vec<f64>,
    /// Per parameter: did the perturbation actually move any residual?  `false`
    /// means the entry is a placeholder zero, not a computed gradient — the same
    /// contract as [`crate::adjoint::Sensitivities::reached`], and for the same
    /// reason: a wrong zero stalls an optimiser at a point that looks stationary.
    pub reached: Vec<bool>,
    /// Per parameter: relative disagreement between the two finite-difference
    /// step sizes, before Richardson extrapolation removed it.  A conservative
    /// error bar on that gradient entry.
    pub fd_error: Vec<f64>,
}

impl TranSensitivities {
    /// Parameters that could not be reached, for a caller that would rather fail
    /// than optimise against placeholder zeros.
    pub fn unreached<'a>(&self, params: &'a [ParamRef]) -> Vec<&'a ParamRef> {
        params
            .iter()
            .zip(self.reached.iter())
            .filter(|(_, ok)| !**ok)
            .map(|(p, _)| p)
            .collect()
    }
}

impl TranAdjoint {
    /// Run the transient, capturing the per-timestep state the adjoint needs.
    ///
    /// Same integrator, same step, same answers as
    /// [`crate::tran::tran_nr_with_registry_opts`] — this drives the identical
    /// [`TranStepper`] and only reads more out of it.
    pub fn run(
        netlist: &Netlist,
        registry: &DeviceRegistry,
        opts: &SimOptions,
        step: f64,
        stop: f64,
    ) -> Result<Self, SimError> {
        let mut st = TranStepper::new(netlist.clone(), registry, opts, step)?;
        reject_inductance(netlist, &st)?;

        let frozen = frozen_columns(netlist, st.topology(), st.devices());
        let mut out = TranAdjoint {
            netlist: netlist.clone(),
            opts: opts.clone(),
            step: st.step_size(),
            topo: st.topology().clone(),
            t: Vec::new(),
            x: Vec::new(),
            a: Vec::new(),
            jq: Vec::new(),
            alpha: Vec::new(),
            beta: Vec::new(),
            mode: Vec::new(),
            dc_start: !opts.uic,
            frozen,
        };

        // ── t = 0 ──────────────────────────────────────────────────────────
        // The initial condition is a constraint whose Jacobian is the DC one,
        // and whose residual carries `p` — unless UIC imposed it, in which case
        // it carries nothing and needs no co-state.
        let x0 = st.solution().to_vec();
        let mode0 = st.step_mode();
        let a0 = if out.dc_start {
            out.stamp_and_repair(&mut st, &x0, true)
        } else {
            Vec::new()
        };
        // `J_q(x_0)` still matters: step 1 integrates `q(x_0)`, so `x_0` reaches
        // the objective through the charge even when it reaches nothing else.
        let (_, jq0) = out.capture(&mut st, 0.0, mode0, &x0);
        out.push(0.0, x0, a0, jq0, 0.0, 0.0, mode0);

        // ── the run ────────────────────────────────────────────────────────
        let eps = out.step * 1e-9;
        let mut t = 0.0;
        while t < stop - eps {
            let t_next = t + out.step;
            let mode = st.step_mode();
            st.solve_at(t_next)?;
            st.commit(t_next);
            let xk = st.solution().to_vec();

            // Between `commit` and `advance_history` is the only window where
            // the companions still describe *this* step.
            let (ak, jqk) = out.capture(&mut st, t_next, mode, &xk);
            // The probing above left the devices linearised at a probe point and
            // the companions built for the second step size; the history advance
            // reads both, so put them back first.
            st.prepare(t_next, mode, out.step);
            st.stamp_at(&xk);
            st.advance_history();

            let alpha = alpha_of(mode, out.step);
            let beta = beta_of(mode);
            out.push(t_next, xk, ak, jqk, alpha, beta, mode);
            t = t_next;
        }

        Ok(out)
    }

    /// Accepted timepoints, in seconds.  `time()[0]` is always 0.
    pub fn time(&self) -> &[f64] {
        &self.t
    }

    /// The MNA solution at each timepoint.
    pub fn trajectory(&self) -> &[Vec<f64>] {
        &self.x
    }

    /// One output's value at every timepoint.
    ///
    /// The forward half of a gradient loop: read the waveform, build a loss from
    /// it, hand the loss's derivative back through [`TranAdjoint::weighted`].
    pub fn signal(&self, out: &Output) -> Result<Vec<f64>, SimError> {
        self.x
            .iter()
            .map(|x| Ok(out.seed(&self.topo, x)?.0))
            .collect()
    }

    /// `(L, ∂L/∂x_k for every k)` for `L = Σ_k weights[k]·out(x_k)`.
    ///
    /// `weights` has one entry per timepoint, so the common objectives are just
    /// choices of weight: a one in the last slot is the final value, `dt`
    /// everywhere is a time integral, and `2·(v_k − target_k)` — which needs the
    /// trajectory, hence [`TranAdjoint::trajectory`] — is a waveform match.
    pub fn weighted(
        &self,
        out: &Output,
        weights: &[f64],
    ) -> Result<(f64, Vec<Vec<f64>>), SimError> {
        if weights.len() != self.t.len() {
            return Err(SimError::ParameterError(format!(
                "weights has length {} but the run has {} timepoints",
                weights.len(),
                self.t.len()
            )));
        }
        let mut value = 0.0;
        let mut seeds = Vec::with_capacity(self.t.len());
        for (x, w) in self.x.iter().zip(weights.iter()) {
            let (v, mut s) = out.seed(&self.topo, x)?;
            value += w * v;
            for e in &mut s {
                *e *= w;
            }
            seeds.push(s);
        }
        Ok((value, seeds))
    }

    /// `dL/dp` for every parameter, given `∂L/∂x_k` at every timepoint.
    ///
    /// One transposed solve per timestep, then four residual replays per
    /// parameter — a central difference at two step sizes, Richardson
    /// extrapolated, as [`crate::adjoint::dc_sensitivity`] does it and for the
    /// same reason: a step scaled to the parameter is not always scaled to the
    /// quantity being differenced.  A parameter whose default step turns out to
    /// be noise-limited costs a few more rounds while the step is walked up; one
    /// with [`ParamRef::step`] pinned costs exactly four.
    pub fn gradient(
        &self,
        registry: &DeviceRegistry,
        seeds: &[Vec<f64>],
        params: &[ParamRef],
    ) -> Result<TranSensitivities, SimError> {
        if seeds.len() != self.t.len() {
            return Err(SimError::ParameterError(format!(
                "seeds has {} timepoints but the run has {}",
                seeds.len(),
                self.t.len()
            )));
        }
        if let Some(bad) = seeds.iter().find(|s| s.len() != self.topo.size) {
            return Err(SimError::ParameterError(format!(
                "a seed has length {} but the system has {} unknowns",
                bad.len(),
                self.topo.size
            )));
        }

        let lambdas = self.backward(seeds)?;
        // The replay at the nominal parameter values does not depend on which
        // parameter is being differentiated, so one is enough — and any
        // perturbation that reaches nothing reproduces it bit for bit.
        let (_, nominal_residual) = self.replay(
            registry,
            &Handle::NetlistElement,
            &ParamRef::new("", ""),
            0.0,
            &lambdas,
        )?;

        let mut grad = vec![0.0; params.len()];
        let mut reached = vec![false; params.len()];
        let mut fd_error = vec![0.0; params.len()];
        let names = crate::newton::device_element_names(&self.netlist);

        for (pi, p) in params.iter().enumerate() {
            let Some((handle, nominal)) = resolve(p, &self.netlist, &names) else {
                continue;
            };
            // Choose the step by measurement, not by rule.  `∛ε·|p|` is the
            // central-difference optimum only when the differenced quantity
            // varies on the parameter's own scale, and `Σ_k λ_kᵀ·G_k` often does
            // not: a 200 fF capacitor whose whole contribution to an optical
            // objective arrives through co-state rows carrying 1e9-scale
            // wavelength terms has a noise floor four orders above what the rule
            // assumes, and at the rule's step the difference is pure roundoff —
            // 6 % wrong, quietly.
            //
            // So walk the step up while the two-size disagreement keeps
            // improving.  That disagreement is already computed, already the
            // reported error bar, and is exactly the thing being minimised;
            // there is nothing else to tune it against.
            let pinned = p.step;
            let mut h = pinned.unwrap_or_else(|| default_step(nominal));
            // Past a percent of the parameter, genuine curvature dominates and a
            // larger step stops helping — so that is where the walk stops.
            let ceiling = 1e-2 * nominal.abs();
            let mut best: Option<(f64, f64)> = None;

            loop {
                let mut d = [0.0_f64; 2];
                for (i, step) in [h, 0.5 * h].into_iter().enumerate() {
                    let (plus, mp) = self.replay(registry, &handle, p, nominal + step, &lambdas)?;
                    let (minus, mm) =
                        self.replay(registry, &handle, p, nominal - step, &lambdas)?;
                    d[i] = (plus - minus) / (2.0 * step);
                    reached[pi] |= mp != nominal_residual || mm != nominal_residual;
                }
                let (coarse, fine) = (d[0], d[1]);
                // A central difference has an O(h²) error term, so
                // `(4·D(h/2) − D(h))/3` cancels it and leaves O(h⁴).
                let value = (4.0 * fine - coarse) / 3.0;
                let spread = if value != 0.0 {
                    (fine - coarse).abs() / value.abs()
                } else {
                    0.0
                };
                if best.is_none_or(|(s, _)| spread < s) {
                    best = Some((spread, value));
                }
                if pinned.is_some() || spread <= STEP_TARGET || 4.0 * h > ceiling {
                    break;
                }
                h *= 4.0;
            }

            let (spread, value) = best.expect("the loop runs at least once");
            fd_error[pi] = spread;
            // `dL/dp = −Σ_k λ_kᵀ·∂G_k/∂p`, and the replay differenced exactly
            // that sum.
            grad[pi] = -value;
        }

        Ok(TranSensitivities {
            grad,
            reached,
            fd_error,
        })
    }

    // ── forward-pass internals ────────────────────────────────────────────

    #[allow(clippy::too_many_arguments)]
    fn push(
        &mut self,
        t: f64,
        x: Vec<f64>,
        a: Vec<SparseRow>,
        jq: Vec<SparseRow>,
        alpha: f64,
        beta: f64,
        mode: IntegratorMode,
    ) {
        self.t.push(t);
        self.x.push(x);
        self.a.push(a);
        self.jq.push(jq);
        self.alpha.push(alpha);
        self.beta.push(beta);
        self.mode.push(mode);
    }

    /// `(A_k, J_q(x_k))` at a frozen `x`, by stamping the same point at two step
    /// sizes.  Leaves the stepper linearised at the second one — every caller
    /// re-stamps afterwards.
    ///
    /// The charge Jacobian is differenced from the **raw** stamps, before the
    /// frozen columns are re-derived, and that ordering is load-bearing.  A
    /// re-derived column is a finite difference in its own right; differencing
    /// two of them cancels four significant figures and then divides by
    /// `Δα ≈ 1/2h`, which for an optical wavelength wire — where `∂f/∂λ` runs to
    /// 1e9 — turns 1e-16 of roundoff into a spurious 1e-15 farad.  Against a
    /// real 100 fF that is a 1 % error in the gradient, and it was one.
    ///
    /// Nothing is lost by taking the raw pair: a frozen coefficient is frozen at
    /// the previous iterate, so it does not depend on `h` and differences to
    /// exactly zero, while the reactive companions that share those columns are
    /// stamped normally and difference correctly.
    ///
    /// ponytail: what that does forgo is a *charge* coupling hidden inside a
    /// frozen coefficient.  No built-in freezes one — every declaration in the
    /// tree is a resistive electro-optic term — and there is no cheap way to
    /// recover it without a `∂q/∂x` hook on `Device`.  Add the hook if a model
    /// ever needs it.
    fn capture(
        &self,
        st: &mut TranStepper,
        t: f64,
        mode: IntegratorMode,
        x: &[f64],
    ) -> (Vec<SparseRow>, Vec<SparseRow>) {
        st.prepare(t, mode, self.step);
        st.stamp_at(x);
        let raw = st.matrix().a.clone();
        st.prepare(t, mode, 2.0 * self.step);
        st.stamp_at(x);
        let d_alpha = alpha_of(mode, self.step) - alpha_of(mode, 2.0 * self.step);
        let jq = scaled_difference(&raw, &st.matrix().a, 1.0 / d_alpha);

        st.prepare(t, mode, self.step);
        (self.stamp_and_repair(st, x, false), jq)
    }

    /// Stamp at `x` and re-derive the frozen columns, as
    /// [`crate::adjoint::dc_sensitivity`] does — a frozen block contributes a
    /// silent zero to every gradient path through it, and an electro-optic
    /// modulator is nothing but such a path.
    fn stamp_and_repair(&self, st: &mut TranStepper, x: &[f64], dc: bool) -> Vec<SparseRow> {
        let repaired: Vec<(usize, Vec<f64>)> = self
            .frozen
            .iter()
            .map(|&c| (c, fd_column(st, c, x, dc)))
            .collect();
        // The probes above moved the linearisation; put it back at `x`.
        if dc {
            st.stamp_dc_at(x);
        } else {
            st.stamp_at(x);
        }
        let mut a = st.matrix().a.clone();
        for (col, column) in &repaired {
            for (row, v) in column.iter().enumerate() {
                // Only touch cells that carry something or already exist —
                // writing every row would densify the column.
                if *v != 0.0 || a[row][*col] != 0.0 {
                    a[row][*col] = *v;
                }
            }
        }
        a
    }

    // ── backward pass ─────────────────────────────────────────────────────

    /// `λ_k` for every timepoint, walked backwards.
    fn backward(&self, seeds: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, SimError> {
        let n = self.topo.size;
        let last = self.t.len() - 1;
        let solver = self.opts.linear_solver(n);
        let mut lambdas = vec![Vec::new(); self.t.len()];
        // ū_{k+1}: the objective's sensitivity to the current stored at step
        // k+1.  Zero past the end of the run.
        let mut u_next = vec![0.0; n];
        let mut rhs = vec![0.0; n];
        let mut scratch = vec![0.0; n];

        for k in (0..=last).rev() {
            rhs.copy_from_slice(&seeds[k]);
            if k < last {
                // (α_{k+1} − β_{k+1}·α_k)·J_q(x_k)ᵀ·ū_{k+1}
                let c = self.alpha[k + 1] - self.beta[k + 1] * self.alpha[k];
                if c != 0.0 {
                    transpose_mul(&self.jq[k], &u_next, &mut scratch);
                    for (r, s) in rhs.iter_mut().zip(scratch.iter()) {
                        *r += c * s;
                    }
                }
            }

            if k == 0 && !self.dc_start {
                // UIC: `x_0` is imposed, so there is no constraint to attach a
                // co-state to and nothing downstream reads `λ_0`.
                lambdas[0] = vec![0.0; n];
                break;
            }

            // ponytail: one symbolic factorisation per timestep.  The frozen-
            // column repair can insert a cell the stamped pattern lacks, so the
            // per-step patterns are not provably identical and a single cached
            // handle cannot be reused blindly.  Prove that and this becomes one
            // `factorise` plus N refactors.
            let mut fact = solver.factorise(&self.a[k])?;
            let lambda = fact.refactor_and_solve_transpose(&self.a[k], &rhs)?;

            // ū_k = λ_k + β_{k+1}·ū_{k+1}
            let b = if k < last { self.beta[k + 1] } else { 0.0 };
            for (u, l) in u_next.iter_mut().zip(lambda.iter()) {
                *u = l + b * *u;
            }
            lambdas[k] = lambda;
        }
        Ok(lambdas)
    }

    /// `Σ_k λ_kᵀ·G_k` with one parameter held at `value` and every `x_k` frozen.
    ///
    /// Re-walks the run without solving anything: each step is one stamp, and
    /// the reactive history is propagated by the same code that propagated it
    /// forward, so every path a parameter takes into the history is included
    /// without being enumerated.
    fn replay(
        &self,
        registry: &DeviceRegistry,
        handle: &Handle,
        p: &ParamRef,
        value: f64,
        lambdas: &[Vec<f64>],
    ) -> Result<(f64, f64), SimError> {
        let mut work = self.netlist.clone();
        if matches!(handle, Handle::NetlistElement) {
            // Before the stepper is built: the `t = 0` reactive history is
            // seeded from these values.
            apply(handle, p, value, &mut work, &mut []);
        }
        let mut st = TranStepper::seeded(work, registry, &self.opts, self.step, &self.x[0])?;
        if let Handle::Device(i) = handle {
            st.set_device_param(*i, &p.param, value);
        }

        let mut acc = 0.0;
        // Sum of squares of every residual seen, so a parameter that moved
        // nothing produces a bit-identical replay and can be reported as
        // unreached rather than as a computed zero.
        let mut moved = 0.0;
        let mut r = vec![0.0; self.topo.size];
        if self.dc_start {
            st.stamp_dc_at(&self.x[0]);
            st.matrix().residual_into(&self.x[0], &mut r);
            acc += dot(&lambdas[0], &r);
            moved += dot(&r, &r);
        }
        for (k, x) in self.x.iter().enumerate().skip(1) {
            st.prepare(self.t[k], self.mode[k], self.step);
            st.stamp_at(x);
            st.matrix().residual_into(x, &mut r);
            acc += dot(&lambdas[k], &r);
            moved += dot(&r, &r);
            st.force_solution(x);
            st.commit(self.t[k]);
            st.advance_history();
        }
        Ok((acc, moved))
    }

    /// The row layout of the captured trajectory.
    pub fn topology(&self) -> &CircuitTopology {
        &self.topo
    }
}

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

/// Two-step-size disagreement below which the finite difference is taken to be
/// step-limited rather than noise-limited, and the walk stops.
const STEP_TARGET: f64 = 1e-6;

/// `α` for a step: the factor multiplying `∂q/∂x` in the stamped Jacobian.
/// Taken from [`crate::reactive::conductance`] at unit capacitance so it cannot
/// drift from what the integrator actually stamped.
fn alpha_of(mode: IntegratorMode, h: f64) -> f64 {
    conductance(ReactiveKind::Capacitor, 1.0, mode, h, None)
}

/// `β`: the coefficient on the *previous* step's stored current.  Only
/// Trapezoidal carries one, which is why only Trapezoidal needs a second
/// co-state.
fn beta_of(mode: IntegratorMode) -> f64 {
    match mode {
        IntegratorMode::Trapezoidal => -1.0,
        _ => 0.0,
    }
}

/// One column of `∂f/∂x` by central difference at a frozen trajectory point.
fn fd_column(st: &mut TranStepper, col: usize, x: &[f64], dc: bool) -> Vec<f64> {
    let n = x.len();
    // A λ wire is ~1.55e-6 and a node voltage is O(1); scale to whichever the
    // column holds, with a floor so a column sitting at zero still gets probed.
    let h = 6.0554545e-6_f64 * x[col].abs().max(1e-3);
    let mut probe = x.to_vec();
    let mut f_plus = vec![0.0; n];
    let mut f_minus = vec![0.0; n];

    for (delta, out) in [(h, &mut f_plus), (-h, &mut f_minus)] {
        probe[col] = x[col] + delta;
        if dc {
            st.stamp_dc_at(&probe);
        } else {
            st.stamp_at(&probe);
        }
        st.matrix().residual_into(&probe, out);
    }

    let scale = 1.0 / (2.0 * h);
    f_plus
        .iter()
        .zip(f_minus.iter())
        .map(|(p, m)| (p - m) * scale)
        .collect()
}

/// `(a − b)·inv`, keeping `a`'s sparsity.  `b` is stamped from the same plan at
/// a different step size, so it cannot carry a cell `a` lacks.
fn scaled_difference(a: &[SparseRow], b: &[SparseRow], inv: f64) -> Vec<SparseRow> {
    a.iter()
        .zip(b.iter())
        .map(|(ra, rb)| {
            let (cols, vals) = ra.entries();
            let cells: Vec<(u32, f64)> = cols
                .iter()
                .zip(vals.iter())
                .filter_map(|(&j, &v)| {
                    let d = (v - rb[j as usize]) * inv;
                    (d != 0.0).then_some((j, d))
                })
                .collect();
            SparseRow::from_sorted_cells(cells)
        })
        .collect()
}

/// `out = Aᵀ·u`, over CSR rows.
fn transpose_mul(a: &[SparseRow], u: &[f64], out: &mut [f64]) {
    out.fill(0.0);
    for (i, row) in a.iter().enumerate() {
        let ui = u[i];
        if ui == 0.0 {
            continue;
        }
        for (j, v) in row.iter() {
            out[j] += v * ui;
        }
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Inductance is not differentiable here yet — say so rather than return a
/// plausible wrong number.
fn reject_inductance(netlist: &Netlist, st: &TranStepper) -> Result<(), SimError> {
    use fairchild_parser::Element;
    let named = netlist.elements.iter().find_map(|el| match el {
        Element::Inductor { name, .. } | Element::CoupledInductors { name, .. } => Some(name),
        _ => None,
    });
    let from_device = st
        .devices()
        .iter()
        .flat_map(|d| d.reactive_branches())
        .any(|b| b.kind == ReactiveKind::Inductor);
    if named.is_some() || from_device {
        return Err(SimError::ParameterError(format!(
            "transient adjoint: inductance is not differentiable yet ({}).  A \
             Norton-companion inductor scales as 1/α rather than α and carries \
             its flux as hidden state, so it fits neither the charge-Jacobian \
             extraction nor the co-state recursion",
            named.map_or("a device-declared branch".to_string(), |n| format!("'{n}'"))
        )));
    }
    Ok(())
}
