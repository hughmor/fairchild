//! Companion-model history for reactive branches — one implementation, shared
//! by every integrator and every source of reactance.
//!
//! # Why this module exists
//!
//! The transient code has three independent concerns:
//!
//! 1. **Where reactance comes from** — netlist `C`, netlist `L`, coupled `K`
//!    groups, and `Device::reactive_branches()`.
//! 2. **How history is represented** and turned into a companion model.
//! 3. **Step control** — fixed `h`, or LTE-driven variable `h` with a
//!    predictor and accept/reject.
//!
//! It used to be organised by (3): each integrator carried its own answers to
//! (1) and (2). That multiplies out to a matrix of `source × integrator ×
//! method` cells, every one hand-written, and a bug is any cell nobody filled
//! in. Four were empty — device branches dropped entirely under variable step,
//! then integrated at the wrong order under both, then coupled inductors
//! diverging under variable step. Nothing caught them because nothing *forced*
//! a new source or a new method to be handled everywhere.
//!
//! So (1) and (2) live here instead, and integrators keep only (3).
//!
//! # Raw state, not companion pairs
//!
//! History is stored as **physical** quantities:
//!
//! | | `state` | `aux` |
//! |---|---|---|
//! | capacitor | `v_C` | `i_C` |
//! | inductor | `i_L` | `v_L` |
//!
//! That choice is what lets one representation serve both integrators. Storing
//! `(G_eq, I_hist)` instead — as the fixed-step path once did — bakes the step
//! size into the state: the Trapezoidal recursion
//! `I_hist' = 2·G·v − I_hist` is only valid while `h` holds still, which is why
//! the variable-step integrator could never offer TR and had to keep its own
//! parallel history. From `(state, aux)` every method's companion is derivable
//! at any `h`, so that split disappears and TR works under variable step too.

use indexmap::IndexMap;

use fairchild_parser::{Element, Netlist};

use crate::device::{Device, ReactiveBranchSpec, ReactiveKind};
use crate::mna::{cap_companion, cap_companion_gear2, ind_companion, ind_companion_gear2};
use crate::mna::{CircuitTopology, MnaMatrix};
use crate::tran::{coupled_inductor_currents, IntegratorMode};

/// Physical history of one reactive branch, as of the last accepted timepoint.
#[derive(Clone, Copy, Debug)]
pub struct BranchHistory {
    /// C (farads) or L (henries). Stored rather than re-read from the device
    /// because it is bias-dependent, and after a rejected step the device's
    /// cache reflects the rejected iterate.
    pub value: f64,
    /// `v_C` for a capacitor, `i_L` for an inductor.
    pub state: f64,
    /// `i_C` for a capacitor, `v_L` for an inductor. Only Trapezoidal reads it;
    /// zero from a DC operating point, where a cap carries no current.
    pub aux: f64,
    /// `state` one timepoint further back. `None` until two steps have been
    /// accepted, which is what gates BDF-2 for this branch.
    pub state_prev2: Option<f64>,
}

impl BranchHistory {
    /// Seed from an operating point: `state` from the solution, no current
    /// through a capacitor (or, by the established convention, none through an
    /// inductor either), no second-order history yet.
    pub fn seeded(kind: ReactiveKind, value: f64, across: f64) -> Self {
        BranchHistory {
            value,
            state: match kind {
                ReactiveKind::Capacitor => across,
                // Inductors start from zero current — the convention the
                // built-in `ind_i` seeding has always used.
                ReactiveKind::Inductor => 0.0,
            },
            aux: 0.0,
            state_prev2: None,
        }
    }
}

/// Norton companion `(G_eq, I_hist)` for one branch at step size `h`.
///
/// This is the only place an integration method is interpreted. A new method is
/// one arm here, and no caller can forget it because `mode` is a parameter of
/// the single function everyone calls.
///
/// `gear2_h_prev` is `Some(h_prev)` only when BDF-2 is permitted for this step
/// — mode is GEAR, no recent rejection, and the step ratio is sane. GEAR
/// without it demotes to BE, which is ordinary order control.
pub fn companion(
    kind: ReactiveKind,
    hist: &BranchHistory,
    mode: IntegratorMode,
    h: f64,
    gear2_h_prev: Option<f64>,
) -> (f64, f64) {
    let value = hist.value;
    if let (IntegratorMode::Gear, Some(h_prev), Some(prev2)) =
        (mode, gear2_h_prev, hist.state_prev2)
    {
        return match kind {
            ReactiveKind::Capacitor => cap_companion_gear2(value, h, h_prev, hist.state, prev2),
            ReactiveKind::Inductor => ind_companion_gear2(value, h, h_prev, hist.state, prev2),
        };
    }
    match (kind, mode) {
        // TR from physical state.  i_C = G·v − I_hist with G = 2C/h and
        // I_hist = G·v_n + i_C(t_n) — no dependence on the previous h, which is
        // the whole point of storing `aux`.
        (ReactiveKind::Capacitor, IntegratorMode::Trapezoidal) => {
            let g = 2.0 * value / h;
            (g, g * hist.state + hist.aux)
        }
        // i_L = G·v + I_hist with G = h/2L and I_hist = i_L(t_n) + G·v_L(t_n).
        (ReactiveKind::Inductor, IntegratorMode::Trapezoidal) => {
            let g = h / (2.0 * value.max(1e-30));
            (g, hist.state + g * hist.aux)
        }
        (ReactiveKind::Capacitor, _) => cap_companion(value, h, hist.state),
        (ReactiveKind::Inductor, _) => ind_companion(value, h, hist.state),
    }
}

/// The same method interpretation as [`companion`], expressed on charge rather
/// than on a branch value.
///
/// [`companion`] is parametrised by a C or an L, which presumes the reactance is
/// a two-terminal branch with one number describing it. A Verilog-A device has
/// no such number: its `ddt` contributions form a charge *vector* `q(x)` whose
/// Jacobian `∂q_i/∂x_j` need not be symmetric — transcapacitance is exactly what
/// a BSIM-class model is made of. This form takes `(q, i)` instead of
/// `(value, v)`, so it applies **per row** and serves both.
///
/// Returns `(i_n, dq_scale)`: the branch current at this iterate, and the factor
/// such that `dq_scale · ∂q/∂x` is its Jacobian. `dq_scale` is deliberately the
/// same quantity [`conductance`] returns for `value = 1`, so the two forms
/// cannot drift apart — [`charge_current_agrees_with_scalar_companion`] pins
/// that.
///
/// `q_prev2` is `Some` only when a second charge history exists; BDF-2 also
/// needs `gear2_h_prev`, and demotes to BE without either, which is ordinary
/// order control.
///
/// (This `(q, i)` formulation is also the prerequisite noted on [`companion`]
/// for offering real Trapezoidal on the variable-step path — it carries no
/// `h`-dependent history term.)
pub fn charge_current(
    mode: IntegratorMode,
    h: f64,
    gear2_h_prev: Option<f64>,
    q_new: f64,
    q_prev: f64,
    q_prev2: Option<f64>,
    i_prev: f64,
) -> (f64, f64) {
    if let (IntegratorMode::Gear, Some(h_prev), Some(q2)) = (mode, gear2_h_prev, q_prev2) {
        // Variable-step BDF-2, matching `cap_companion_gear2` term for term
        // under q = C·v.
        let rho = h / h_prev;
        let denom = h * (1.0 + rho);
        let scale = (1.0 + 2.0 * rho) / denom;
        let i = scale * q_new - (1.0 + rho) / h * q_prev + (rho * rho) / denom * q2;
        return (i, scale);
    }
    match mode {
        // TR: (q_n − q_{n−1})·2/h = i_n + i_{n−1}.
        IntegratorMode::Trapezoidal => {
            let scale = 2.0 / h;
            (scale * (q_new - q_prev) - i_prev, scale)
        }
        // BE, and GEAR before it has the history for second order.
        _ => {
            let scale = 1.0 / h;
            (scale * (q_new - q_prev), scale)
        }
    }
}

/// One charge branch's history, for a device that stamps its own reactance.
///
/// [`BranchHistory`] serves devices that *declare* a reactive branch and let
/// the shared stamper integrate it. A device that stamps its own — the diode's
/// `Cj + TT·gd`, the MOSFET's Meyer and depletion caps, a Verilog-A `ddt` —
/// only ever receives `alpha`, which cannot express anything but Backward
/// Euler. This is the state those devices need to build any other method's
/// companion, plus the two operations on it, so they cannot each grow their own
/// answer to "what does `method` mean" — the failure mode this module exists to
/// prevent.
#[derive(Clone, Copy, Debug, Default)]
pub struct ChargeHistory {
    q_prev: f64,
    /// `None` until two timepoints have been accepted; gates BDF-2, exactly as
    /// [`BranchHistory::state_prev2`] does.
    q_prev2: Option<f64>,
    /// Branch current at the last accepted timepoint. Only Trapezoidal reads
    /// it; zero from a DC operating point, where a capacitor carries none.
    i_prev: f64,
}

impl ChargeHistory {
    /// Companion for this branch at the current Newton iterate.
    ///
    /// `q_new` is the branch charge there; `cv` is `C·v` at the same point —
    /// the term the Jacobian stamp itself contributes to the residual, which
    /// the history current has to cancel. (For a linear cap the two are equal;
    /// for a depletion cap they are not, which is the whole reason both are
    /// arguments.)
    ///
    /// Returns `(i_hist, scale)`: the current to stamp from the branch's
    /// positive node to its negative one, and the factor multiplying `C` in the
    /// Jacobian. Taking both from one call is what keeps residual and Jacobian
    /// linearised about the same point.
    ///
    /// `disc` is `None` outside the transient loop — a caller with no method to
    /// offer, where the `alpha` it passes is Backward Euler by definition.
    pub fn companion(
        &self,
        disc: Option<crate::device::Discretisation>,
        alpha: f64,
        q_new: f64,
        cv: f64,
    ) -> (f64, f64) {
        let (i_n, scale) = match disc {
            Some(d) => charge_current(
                d.mode,
                d.h,
                d.gear2_h_prev,
                q_new,
                self.q_prev,
                self.q_prev2,
                self.i_prev,
            ),
            None => (alpha * (q_new - self.q_prev), alpha),
        };
        (scale * cv - i_n, scale)
    }

    /// The Jacobian factor alone, for the stamp that has no charge to hand.
    ///
    /// [`conductance`] is linear in the branch value, so evaluating it at 1.0
    /// gives exactly the scalar — and it is the same quantity
    /// [`charge_current`] returns as `scale`, so the two stamps cannot drift.
    pub fn scale(disc: Option<crate::device::Discretisation>, alpha: f64) -> f64 {
        match disc {
            Some(d) => conductance(ReactiveKind::Capacitor, 1.0, d.mode, d.h, d.gear2_h_prev),
            None => alpha,
        }
    }

    /// Roll past an accepted step, `q_now` being the branch charge at the
    /// converged solution.
    ///
    /// The current *entering* this step becomes Trapezoidal's history for the
    /// next one, and must be computed before `q_prev` is overwritten.
    pub fn advance(&mut self, disc: Option<crate::device::Discretisation>, q_now: f64) {
        if let Some(d) = disc {
            let (i_n, _) = charge_current(
                d.mode,
                d.h,
                d.gear2_h_prev,
                q_now,
                self.q_prev,
                self.q_prev2,
                self.i_prev,
            );
            self.i_prev = i_n;
        }
        self.q_prev2 = Some(self.q_prev);
        self.q_prev = q_now;
    }
}

/// Conductance alone, for a branch whose `value` is re-queried mid-Newton.
///
/// Device-declared branches have bias-dependent C, so `G_eq` must be recomputed
/// from the current iterate while `I_hist` stays pinned to the accepted
/// timepoint. That pairing is the charge-conserving form:
/// `i_C = C_new·v_new/h − C_old·v_old/h`.
///
/// Both GEAR-2 companions are linear in the branch value, so the BDF-2
/// coefficient factors out exactly and this stays consistent with [`companion`].
pub fn conductance(
    kind: ReactiveKind,
    value: f64,
    mode: IntegratorMode,
    h: f64,
    gear2_h_prev: Option<f64>,
) -> f64 {
    let v = value.max(1e-30);
    if let (IntegratorMode::Gear, Some(h_prev)) = (mode, gear2_h_prev) {
        let rho = h / h_prev;
        return match kind {
            ReactiveKind::Capacitor => v * (1.0 + 2.0 * rho) / (h * (1.0 + rho)),
            ReactiveKind::Inductor => h * (1.0 + rho) / (v * (1.0 + 2.0 * rho)),
        };
    }
    match (kind, mode) {
        (ReactiveKind::Capacitor, IntegratorMode::Trapezoidal) => 2.0 * v / h,
        (ReactiveKind::Capacitor, _) => v / h,
        (ReactiveKind::Inductor, IntegratorMode::Trapezoidal) => h / (2.0 * v),
        (ReactiveKind::Inductor, _) => h / v,
    }
}

/// Roll one branch's history past an accepted step.
///
/// Takes no `mode` and no `h`: everything method-specific is already inside the
/// `stamped` companion that was actually used, so the update is
/// method-agnostic. That is why adding a method cannot break the advance.
///
/// `across` is the branch voltage from the converged solution; `new_value` the
/// branch's C or L at the new operating point.
pub fn advance(
    hist: &mut BranchHistory,
    kind: ReactiveKind,
    stamped: (f64, f64),
    across: f64,
    new_value: f64,
) {
    let (g_eq, i_hist) = stamped;
    hist.state_prev2 = Some(hist.state);
    match kind {
        ReactiveKind::Capacitor => {
            hist.state = across;
            hist.aux = g_eq * across - i_hist; // i_C = G·v − I_hist
        }
        ReactiveKind::Inductor => {
            hist.state = g_eq * across + i_hist; // i_L = G·v + I_hist
            hist.aux = across;
        }
    }
    hist.value = new_value;
}

/// Override an inductor's advanced current — for a coupled pair, whose current
/// comes from the mutual form rather than the standalone one.
pub fn set_inductor_current(hist: &mut BranchHistory, current: f64, across: f64) {
    hist.state = current;
    hist.aux = across;
}

/// All reactive history in a circuit, plus the companion buffers the stampers
/// read. Both integrators own one of these and nothing else about reactance.
///
/// The buffers are members rather than return values so the per-step rebuild
/// allocates nothing: `build` clears and refills in place.
/// One netlist `C` or `L`, with everything the per-step work needs resolved up
/// front. No name lookups and no netlist scans inside the loop — those cost
/// O(elements²) per step, which is exactly the shape of regression a shared
/// abstraction is prone to hiding.
struct ElemBranch {
    kind: ReactiveKind,
    /// Resolved MNA row for each terminal; `None` is ground.
    pos: Option<usize>,
    neg: Option<usize>,
    /// Slot in `cap_state` / `ind_state`, so companions are written by index.
    slot: usize,
    hist: BranchHistory,
}

/// A `K` element, with both windings pre-resolved to `elems` indices.
struct CoupledPair {
    a: usize,
    b: usize,
    coupling: f64,
}

/// All reactive history in a circuit, plus the companion buffers the stampers
/// read. Both integrators own one of these and nothing else about reactance.
///
/// The buffers are members rather than return values, and are keyed once at
/// construction: the per-step rebuild overwrites values by index, so it neither
/// allocates nor hashes.
pub struct ReactiveState {
    /// Netlist `C` / `L`, in netlist order.
    elems: Vec<ElemBranch>,
    /// `Device::reactive_branches()` history, `[device][branch]`.
    devs: Vec<Vec<BranchHistory>>,
    /// Coupled inductor pairs, pre-resolved.
    coupled: Vec<CoupledPair>,
    /// Companion for each netlist capacitor at the current step.
    pub cap_state: IndexMap<String, (f64, f64)>,
    /// Companion for each netlist inductor at the current step.
    pub ind_state: IndexMap<String, (f64, f64)>,
    /// Companion for each device branch at the current step, `[device][branch]`.
    pub dev_state: Vec<Vec<(f64, f64)>>,
}

impl ReactiveState {
    /// Seed every branch — netlist and device alike — from an operating point.
    /// Seed every reactive branch's history from the operating point `x`.
    ///
    /// `devices` is `&mut` for one reason, and it is load-bearing: a
    /// bias-dependent capacitance reports its value through
    /// `reactive_branches()` from whatever its last `eval` cached, and the
    /// transient builds **fresh** device instances that the DC solve never
    /// touched. Seeding without evaluating first therefore recorded the
    /// constructor's default — `c_j_cached: 0.0` for a depletion junction —
    /// so `q_prev` came out zero while the first step's companion used the real
    /// `C_j(V_op)`. The branch then had to acquire its entire charge in one
    /// timestep: measured on a PN phase shifter reverse-biased at −2.93 V
    /// through 25 Ω, the anode jumped to −0.26 V on the first 1 ps step and
    /// took ~10 ps to recover, identically under BE, TR and GEAR — which is
    /// what says it is the initial charge and not the integrator.
    ///
    /// Evaluating here rather than at the call sites is deliberate: there are
    /// two of them, they are the only places that can get this wrong, and
    /// `&mut` makes it impossible to call this without the eval happening.
    pub fn new(
        netlist: &Netlist,
        topo: &CircuitTopology,
        devices: &mut [Box<dyn Device>],
        ctx: &crate::device::SimContext,
        x: &[f64],
    ) -> Self {
        for dev in devices.iter_mut() {
            dev.eval(x, crate::device::EvalFlags::dc(), ctx);
        }
        let devices: &[Box<dyn Device>] = devices;
        let mut elems = Vec::new();
        let mut cap_state = IndexMap::new();
        let mut ind_state = IndexMap::new();
        let mut by_name: IndexMap<String, usize> = IndexMap::new();

        for el in &netlist.elements {
            let (kind, name, pos, neg, value) = match el {
                Element::Capacitor {
                    name,
                    pos,
                    neg,
                    capacitance,
                } => (ReactiveKind::Capacitor, name, pos, neg, *capacitance),
                Element::Inductor {
                    name,
                    pos,
                    neg,
                    inductance,
                } => (ReactiveKind::Inductor, name, pos, neg, *inductance),
                _ => continue,
            };
            let (p, n) = (
                topo.node_index.get(pos).copied(),
                topo.node_index.get(neg).copied(),
            );
            let across = node_diff(p, n, x);
            let slot = match kind {
                ReactiveKind::Capacitor => {
                    cap_state.insert(name.clone(), (0.0, 0.0));
                    cap_state.len() - 1
                }
                ReactiveKind::Inductor => {
                    ind_state.insert(name.clone(), (0.0, 0.0));
                    ind_state.len() - 1
                }
            };
            by_name.insert(name.clone(), elems.len());
            elems.push(ElemBranch {
                kind,
                pos: p,
                neg: n,
                slot,
                hist: BranchHistory::seeded(kind, value, across),
            });
        }

        let coupled = netlist
            .elements
            .iter()
            .filter_map(|el| {
                let Element::CoupledInductors {
                    l1, l2, coupling, ..
                } = el
                else {
                    return None;
                };
                Some(CoupledPair {
                    a: *by_name.get(l1)?,
                    b: *by_name.get(l2)?,
                    coupling: *coupling,
                })
            })
            .collect();

        let devs: Vec<Vec<BranchHistory>> = devices
            .iter()
            .map(|dev| {
                dev.reactive_branches()
                    .iter()
                    .map(|br| BranchHistory::seeded(br.kind, br.value, branch_voltage(br, x)))
                    .collect()
            })
            .collect();
        let dev_state = devs.iter().map(|d| vec![(0.0, 0.0); d.len()]).collect();

        ReactiveState {
            elems,
            devs,
            coupled,
            cap_state,
            ind_state,
            dev_state,
        }
    }

    /// Rebuild every companion for step size `h` under `mode`.
    ///
    /// Call once per attempted step, before the Newton loop. Safe to call again
    /// with a different `h` — that is what a rejected step needs, and the whole
    /// reason history is physical rather than companion-shaped.
    pub fn build(
        &mut self,
        devices: &[Box<dyn Device>],
        mode: IntegratorMode,
        h: f64,
        gear2_h_prev: Option<f64>,
    ) {
        for e in &self.elems {
            let c = companion(e.kind, &e.hist, mode, h, gear2_h_prev);
            let slot = match e.kind {
                ReactiveKind::Capacitor => self.cap_state.get_index_mut(e.slot),
                ReactiveKind::Inductor => self.ind_state.get_index_mut(e.slot),
            };
            if let Some((_, v)) = slot {
                *v = c;
            }
        }
        for (d, dev) in devices.iter().enumerate() {
            for (b, br) in dev.reactive_branches().iter().enumerate() {
                self.dev_state[d][b] = companion(br.kind, &self.devs[d][b], mode, h, gear2_h_prev);
            }
        }
    }

    /// Roll all history forward past an accepted solution.
    ///
    /// Needs neither `mode` nor `h`: the companions in `cap_state` / `ind_state`
    /// / `dev_state` are the ones that were actually stamped, so everything
    /// method-specific is already folded into them.
    pub fn accept(&mut self, devices: &[Box<dyn Device>], x: &[f64]) {
        // Coupled pairs need what was stamped, before the advance overwrites it.
        let stamped: Vec<(f64, f64)> = if self.coupled.is_empty() {
            Vec::new()
        } else {
            self.elems.iter().map(|e| self.stamped_of(e)).collect()
        };

        for e in &mut self.elems {
            let stamped_e = match e.kind {
                ReactiveKind::Capacitor => self.cap_state[e.slot],
                ReactiveKind::Inductor => self.ind_state[e.slot],
            };
            let across = node_diff(e.pos, e.neg, x);
            let value = e.hist.value; // netlist values are constant
            advance(&mut e.hist, e.kind, stamped_e, across, value);
        }

        // Coupled inductors: the standalone advance above missed the mutual
        // term, so redo those two currents from the mutual companion.  Left as a
        // post-pass rather than a branch group — a group would need matrix-valued
        // companions throughout, which is the next simplification, not this one.
        for k in &self.coupled {
            let (ea, eb) = (&self.elems[k.a], &self.elems[k.b]);
            let (val_a, val_b) = (ea.hist.value, eb.hist.value);
            let (vla, vlb) = (node_diff(ea.pos, ea.neg, x), node_diff(eb.pos, eb.neg, x));
            let ((g_eq_a, i_hist_a), (_, i_hist_b)) = (stamped[k.a], stamped[k.b]);
            if let Some((ia, ib)) = coupled_inductor_currents(
                val_a,
                val_b,
                k.coupling,
                g_eq_a * val_a, // conductance scale: h under BE, 1/α under BDF-2
                vla,
                vlb,
                i_hist_a,
                i_hist_b,
            ) {
                set_inductor_current(&mut self.elems[k.a].hist, ia, vla);
                set_inductor_current(&mut self.elems[k.b].hist, ib, vlb);
            }
        }

        for (d, dev) in devices.iter().enumerate() {
            for (b, br) in dev.reactive_branches().iter().enumerate() {
                let across = branch_voltage(br, x);
                advance(
                    &mut self.devs[d][b],
                    br.kind,
                    self.dev_state[d][b],
                    across,
                    br.value,
                );
            }
        }
    }

    fn stamped_of(&self, e: &ElemBranch) -> (f64, f64) {
        match e.kind {
            ReactiveKind::Capacitor => self.cap_state[e.slot],
            ReactiveKind::Inductor => self.ind_state[e.slot],
        }
    }

    /// Whether BDF-2 has enough history everywhere. BDF-2 on some branches and
    /// BE on others inside one step is an order mismatch, so it is all or none.
    pub fn gear2_ready(&self) -> bool {
        self.elems.iter().all(|e| e.hist.state_prev2.is_some())
            && self
                .devs
                .iter()
                .all(|d| d.iter().all(|h| h.state_prev2.is_some()))
    }
}

/// Voltage between two resolved MNA rows; `None` is ground.
fn node_diff(pos: Option<usize>, neg: Option<usize>, x: &[f64]) -> f64 {
    let vp = pos.map_or(0.0, |i| x[i]);
    let vn = neg.map_or(0.0, |i| x[i]);
    vp - vn
}

/// Stamp the companion of every device-declared reactive branch into `mat`.
///
/// Called inside the Newton loop, after each device's `load_jacobian_tran`, so
/// the device's cached `value` reflects the current iterate.  `G_eq` is
/// recomputed from that fresh value while `I_hist` stays pinned to the accepted
/// timepoint — the charge-conserving pairing for a bias-dependent C.
pub fn stamp_device_branches(
    devices: &[Box<dyn Device>],
    dev_state: &[Vec<(f64, f64)>],
    mat: &mut MnaMatrix,
    x: &[f64],
    h: f64,
    mode: IntegratorMode,
    gear2_h_prev: Option<f64>,
) {
    for (d, dev) in devices.iter().enumerate() {
        for (b, br) in dev.reactive_branches().iter().enumerate() {
            let (_g_old, i_hist) = dev_state[d][b];
            // `g_val` is the residual's coefficient: the branch carries
            // `q = value·v`, so its current is `g_val·v − i_hist`.
            let g_val = conductance(br.kind, br.value, mode, h, gear2_h_prev);
            // `g_jac` is `∂(current)/∂v`, which for a bias-dependent value picks
            // up a second term.  `conductance` is linear in `value` for every
            // mode, so the extra derivative scales the same way and can go
            // through the same helper.
            let g_jac = if br.dvalue_dstate == 0.0 {
                g_val
            } else {
                let v = branch_voltage(br, x);
                g_val + v * conductance(br.kind, br.dvalue_dstate, mode, h, gear2_h_prev)
            };
            // Norton: stamp the true derivative into `A`, and put the difference
            // back on the RHS so the *residual* `A·x − b` is unchanged. Getting
            // this wrong would move the answer, not just the iteration count.
            let v0 = branch_voltage(br, x);
            let comp = (g_jac - g_val) * v0;
            // For an inductor the current is i = G_eq·v + I_hist, so the history
            // adds rather than subtracts.
            let i_sign = match br.kind {
                ReactiveKind::Capacitor => 1.0,
                ReactiveKind::Inductor => -1.0,
            };
            if let Some(p) = br.pos {
                mat.a[p][p] += g_jac;
                if let Some(n) = br.neg {
                    mat.a[p][n] -= g_jac;
                }
                mat.b[p] += i_sign * i_hist + comp;
            }
            if let Some(n) = br.neg {
                mat.a[n][n] += g_jac;
                if let Some(p) = br.pos {
                    mat.a[n][p] -= g_jac;
                }
                mat.b[n] -= i_sign * i_hist + comp;
            }
        }
    }
}

/// Voltage across a device-declared reactive branch, from a solution vector.
pub fn branch_voltage(br: &ReactiveBranchSpec, x: &[f64]) -> f64 {
    match (br.pos, br.neg) {
        (Some(p), Some(n)) => x[p] - x[n],
        (Some(p), None) => x[p],
        (None, Some(n)) => -x[n],
        (None, None) => 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `charge_current` must be the *same* method interpretation as
    /// `companion`, not a second one that happens to look similar. For a
    /// constant `C` the two describe one physical branch, so with `q = C·v`
    /// they have to agree on both the current and the Jacobian — under every
    /// mode, including a non-uniform BDF-2 step.
    ///
    /// This is the test that stops the Verilog-A path and the native path
    /// drifting into two integrators again.
    #[test]
    fn charge_current_agrees_with_scalar_companion() {
        let c = 2.5e-12;
        let cases: [(IntegratorMode, f64, Option<f64>); 5] = [
            (IntegratorMode::BackwardEuler, 1e-9, None),
            (IntegratorMode::Trapezoidal, 1e-9, None),
            (IntegratorMode::Gear, 1e-9, None), // no prev2 yet -> demotes to BE
            (IntegratorMode::Gear, 1e-9, Some(1e-9)), // uniform BDF-2
            (IntegratorMode::Gear, 4e-10, Some(1.3e-9)), // non-uniform BDF-2
        ];

        for (mode, h, gear2_h_prev) in cases {
            // A history with distinct values in every slot, so a dropped or
            // swapped term cannot pass by coincidence.
            let (v_prev, v_prev2, i_prev) = (0.37, 0.11, 4.2e-5);
            let hist = BranchHistory {
                value: c,
                state: v_prev,
                aux: i_prev,
                state_prev2: Some(v_prev2),
            };
            let v_new = 0.83;

            // Voltage form: i_C = G_eq·v − I_hist.
            let (g_eq, i_hist) = companion(ReactiveKind::Capacitor, &hist, mode, h, gear2_h_prev);
            let i_voltage_form = g_eq * v_new - i_hist;

            // Charge form, same branch expressed as q = C·v.
            let (i_charge_form, scale) = charge_current(
                mode,
                h,
                gear2_h_prev,
                c * v_new,
                c * v_prev,
                Some(c * v_prev2),
                i_prev,
            );

            let tol = 1e-9 * i_voltage_form.abs().max(1e-12);
            assert!(
                (i_charge_form - i_voltage_form).abs() <= tol,
                "{mode:?} h={h:e} gear2={gear2_h_prev:?}: current \
                 {i_charge_form:e} != {i_voltage_form:e}"
            );
            // And the Jacobians: scale·∂q/∂v = scale·C must equal G_eq.
            assert!(
                (scale * c - g_eq).abs() <= 1e-9 * g_eq.abs(),
                "{mode:?} h={h:e}: scale·C {:e} != G_eq {g_eq:e}",
                scale * c
            );
            // `conductance` is documented as `scale` at value = 1; hold it to that.
            let cond = conductance(ReactiveKind::Capacitor, 1.0, mode, h, gear2_h_prev);
            assert!(
                (cond - scale).abs() <= 1e-9 * scale.abs(),
                "{mode:?} h={h:e}: conductance(1.0) {cond:e} != scale {scale:e}"
            );
        }
    }

    /// The property the whole module rests on: a companion built from physical
    /// state is the same one the incremental recursion would have produced, and
    /// stays correct when `h` changes — which the incremental form cannot.
    #[test]
    fn tr_from_raw_state_survives_a_step_size_change() {
        let (c, h1, h2) = (1e-9, 1e-6, 3.7e-7);
        let mut hist = BranchHistory::seeded(ReactiveKind::Capacitor, c, 0.0);

        // One TR step at h1 towards 1 V, then keep integrating at a different h.
        let mut v = 0.0;
        for (i, h) in [h1, h1, h2, h2, h1].into_iter().enumerate() {
            let (g, i_hist) = companion(
                ReactiveKind::Capacitor,
                &hist,
                IntegratorMode::Trapezoidal,
                h,
                None,
            );
            // Solve the one-node RC: (1/R + g)·v = 1/R·1 + i_hist, R = 1k.
            let g_r = 1.0 / 1e3;
            v = (g_r * 1.0 + i_hist) / (g_r + g);
            advance(&mut hist, ReactiveKind::Capacitor, (g, i_hist), v, c);
            assert!(
                v.is_finite() && (0.0..=1.0).contains(&v),
                "step {i} produced {v}"
            );
            // i_C must stay consistent with the charge balance: i_C = (1−v)/R.
            let expected_i = (1.0 - v) / 1e3;
            assert!(
                (hist.aux - expected_i).abs() < 1e-12,
                "step {i}: tracked i_C {} vs KCL {}",
                hist.aux,
                expected_i
            );
        }
        assert!(v > 0.0, "the cap should have charged");
    }

    /// BE ignores `aux`, so a BE companion must depend only on `value`, `state`
    /// and `h` — the guarantee that switching representation changed nothing for
    /// the methods that were already correct.
    #[test]
    fn be_companion_matches_the_direct_formula() {
        let mut hist = BranchHistory::seeded(ReactiveKind::Capacitor, 2e-12, 0.75);
        hist.aux = 1234.0; // must be ignored
        let got = companion(
            ReactiveKind::Capacitor,
            &hist,
            IntegratorMode::BackwardEuler,
            5e-9,
            None,
        );
        assert_eq!(got, cap_companion(2e-12, 5e-9, 0.75));
    }

    #[test]
    fn conductance_agrees_with_the_full_companion() {
        for kind in [ReactiveKind::Capacitor, ReactiveKind::Inductor] {
            let hist = BranchHistory {
                value: 3e-9,
                state: 0.4,
                aux: 0.1,
                state_prev2: Some(0.3),
            };
            for mode in [
                IntegratorMode::BackwardEuler,
                IntegratorMode::Trapezoidal,
                IntegratorMode::Gear,
            ] {
                for gear2 in [None, Some(7e-7)] {
                    let (g, _) = companion(kind, &hist, mode, 1e-6, gear2);
                    let g2 = conductance(kind, hist.value, mode, 1e-6, gear2);
                    assert!(
                        (g - g2).abs() <= g.abs() * 1e-12,
                        "{kind:?} {mode:?} gear2={gear2:?}: {g} vs {g2}"
                    );
                }
            }
        }
    }

    /// GEAR must demote to BE until two timepoints exist, or BDF-2 reads a
    /// `state_prev2` that was never set.
    #[test]
    fn gear_demotes_to_be_without_history() {
        let hist = BranchHistory::seeded(ReactiveKind::Capacitor, 1e-12, 0.5);
        assert!(hist.state_prev2.is_none());
        let got = companion(
            ReactiveKind::Capacitor,
            &hist,
            IntegratorMode::Gear,
            1e-9,
            Some(1e-9),
        );
        assert_eq!(got, cap_companion(1e-12, 1e-9, 0.5));
    }
}
