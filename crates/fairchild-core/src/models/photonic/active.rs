//! Active optical devices: an [`OpticalSegment`] (the optical core) driven by a
//! [`PhotonicActiveModel`] (the physics that perturbs the mode and contributes
//! to the electrical equations).
//!
//! This is the collapse target for the phase-shifter / modulator family: one
//! generic [`ActiveOpticalDevice`] parameterised by a swappable drive model,
//! instead of one bespoke `Device` impl per (drive × electrical-feature)
//! combination. New device physics — free-carrier dispersion, thermo-optic,
//! Pockels, photoconductive back-action, an externally-tabulated model — is a
//! new `PhotonicActiveModel`, not a new copy of the optical stamp loop.
//!
//! `docs/photonic-models.md` documents every device built this way, including
//! the tier tables for the phase-shifter families.

use super::segment::{OpticalSegment, PerChannel};
use super::{dB_per_cm_to_neper_per_m, stamp_resistor};

/// A depletion junction's capacitance **and its stored charge** at bias `v`.
///
/// ```text
///   C(v) = C_j0·(1 − v/V_bi)^−m          q(v) = ∫₀^v C dv'
///        = C_j0·V_bi/(1−m)·[1 − (1 − v/V_bi)^(1−m)]
/// ```
///
/// The charge is the point of this function. `C(v)·v` is not it, and using
/// `C(v)·v` as the integrator's state makes a `m_j = 0.5` junction come out 2.3x
/// too fast on a 2 V step — see `ReactiveBranchSpec::charge`. The capacitance
/// alone is still what the Jacobian and the small-signal analyses want, because
/// there it is `∂q/∂v`.
///
/// Past `V_bi/2` both switch to the tangent line, so `C` stays finite into
/// forward bias and `q` stays its exact integral: `q(v) = q(knee) +
/// C_knee·Δ + ½·(dC/dv)·Δ²`.
fn junction_cap_and_charge(v: f64, c_j0: f64, v_bi: f64, m_j: f64) -> (f64, f64) {
    let v_knee = 0.5 * v_bi;
    // q(v) for the power law, exact for m ≠ 1 (m_j is clamped below 1).
    let q_law = |v: f64| c_j0 * v_bi / (1.0 - m_j) * (1.0 - (1.0 - v / v_bi).powf(1.0 - m_j));
    if v < v_knee {
        let c = c_j0 / (1.0 - v / v_bi).powf(m_j);
        (c, q_law(v))
    } else {
        let c_knee = c_j0 / (1.0 - v_knee / v_bi).powf(m_j);
        let dc_dv = c_knee * m_j / (v_bi - v_knee);
        let d = v - v_knee;
        // Linear past the knee, so the slope is the tangent's own and the
        // charge is that line's integral.
        (
            c_knee + dc_dv * d,
            q_law(v_knee) + c_knee * d + 0.5 * dc_dv * d * d,
        )
    }
}
use crate::device::{Device, EvalFlags, NodeId, ReactiveBranchSpec, ReactiveKind, SimContext};
use crate::mna::MnaMatrix;
use fairchild_parser::{EvalContext, Expr};

/// Per-segment optical perturbation produced by a [`PhotonicActiveModel`] at the
/// current operating point:
///   - `dn_eff` — effective-index change (a λ-dependent phase `2π·Δn_eff·L/λ`,
///     e.g. free-carrier dispersion or Pockels via index);
///   - `dphi` — a λ-independent direct phase (e.g. a calibrated thermo-optic
///     `φ_th = π·P/P_π`, applied identically to every WDM channel);
///   - `dalpha_neper_m` — excess propagation loss.
///
/// All are mechanism-agnostic and additive — the segment never asks how they
/// were produced — so any drive physics, and any composition of drives, reduces
/// to this triple.
// No longer `Copy`: a per-channel value owns a Vec. Every consumer clones
// explicitly, which keeps the allocation visible at the call site.
#[derive(Clone, Debug, Default)]
pub struct OpticalPerturbation {
    pub dn_eff: PerChannel,
    pub dphi: f64,
    pub dalpha_neper_m: PerChannel,
}

/// The active physics of an optical device: how electrical/thermal state and
/// the optical field in the segment couple.
pub trait PhotonicActiveModel: Send + Sync {
    /// Electrical/thermal terminals beyond the optical bundle (anode/cathode,
    /// heater±, electrode±, …).
    fn num_electrical_terminals(&self) -> usize;

    /// Internal MNA rows the model needs (series-R node, RC thermal node, …).
    fn num_internal_nodes(&self) -> usize {
        0
    }

    /// Finalise band-centre / temperature-derived defaults (mirrors
    /// `Device::setup_model`). Default no-op.
    fn setup_model(&mut self, _ctx: &SimContext) {}

    /// Bind the model to its electrical terminal nodes (length =
    /// `num_electrical_terminals`).
    fn set_terminals(&mut self, electrical: &[NodeId]);

    /// Bind the model's internal MNA rows to a contiguous block starting at
    /// `first_idx`. Default no-op (no internal nodes).
    fn bind_internal(&mut self, _first_idx: usize) {}

    /// Evaluate at the current iterate: read electrical node voltages (from `x`
    /// via the bound terminals) and the per-channel optical `intensity_w`
    /// (enables photoconductive / detection back-action), cache electrical
    /// contributions, and return the optical perturbation the segment applies.
    fn eval(&mut self, x: &[f64], intensity_w: &[f64], ctx: &SimContext) -> OpticalPerturbation;

    /// Stamp the electrical/thermal contributions (junction conductance,
    /// resistor, electrode cap, photocurrent, heat source) into the Jacobian.
    fn stamp(&self, mat: &mut MnaMatrix);

    /// Transient Jacobian stamp (`alpha` = BE/GEAR coefficient). Default falls
    /// back to the DC stamp — correct for models whose electrical contributions
    /// are bias-only (PN, injection) and whose reactances flow through
    /// `reactive_branches`. Override for device-owned discretised state (e.g. a
    /// thermal-RC node whose BE state equation differs from its DC form).
    fn stamp_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.stamp(mat);
    }

    /// Stamp electrical RHS contributions (history sources, fixed currents).
    /// Default no-op.
    fn stamp_residual(&self, _b: &mut [f64]) {}

    /// Transient RHS stamp (`alpha` = BE/GEAR coefficient). Default falls back
    /// to the DC residual. Override alongside `stamp_tran` for device-owned
    /// state.
    fn stamp_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.stamp_residual(b);
    }

    /// Linear/bias-dependent reactive branches (junction C_j, electrode cap)
    /// for the transient integrator and frequency-domain analyses.
    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> {
        Vec::new()
    }

    /// Record post-converged state for the next timestep. Default no-op.
    fn commit(&mut self, _x: &[f64]) {}

    /// Set a model parameter (already lower-cased). Returns whether recognised.
    fn set_param(&mut self, name: &str, value: f64) -> bool;
}

/// A generic optical device: a passive/active [`OpticalSegment`] driven by a
/// [`PhotonicActiveModel`]. Terminal layout: the segment's `2·wpc·N` optical
/// bundle wires first, then the model's electrical terminals.
pub struct ActiveOpticalDevice {
    seg: OpticalSegment,
    model: Box<dyn PhotonicActiveModel>,
    /// The drive's electrical terminals, kept so `eval` can finite-difference
    /// the perturbation against each one (see the Newton cross-terms in
    /// `OpticalSegment::stamp`). Numerical rather than analytic so every drive —
    /// including `ExprDrive`, whose constitutive law is user text — gets exact
    /// Jacobian coupling for free; the cost is a few extra scalar model evals
    /// per iteration, negligible against one linear solve.
    elec_nodes: Vec<NodeId>,
    /// Scratch for the finite difference, reused to avoid per-iteration allocs.
    /// Per control node: `(d Δn/dV, d Δφ/dV, d Δα/dV)`. The two per-channel
    /// members stay `Uniform` for every drive whose effect is shared, so the
    /// common case allocates nothing.
    dpert_dv: Vec<(PerChannel, f64, PerChannel)>,
    /// Smallest valid terminal count, recorded when `setup_instance` is handed a
    /// count it cannot use. Without it `num_terminals()` would report the
    /// unconfigured segment's 0 optical wires, and the caller's error message
    /// would quote a nonsense expectation.
    min_terminals: Option<usize>,
}

impl ActiveOpticalDevice {
    pub fn new(seg: OpticalSegment, model: Box<dyn PhotonicActiveModel>) -> Self {
        ActiveOpticalDevice {
            seg,
            model,
            elec_nodes: Vec::new(),
            dpert_dv: Vec::new(),
            min_terminals: None,
        }
    }
}

impl Device for ActiveOpticalDevice {
    fn num_terminals(&self) -> usize {
        // A refused setup leaves the segment with zero optical wires, so report
        // the smallest count that would have worked instead of a misleading 0+N.
        if let Some(min) = self.min_terminals {
            return min;
        }
        self.seg.num_optical_wires() + self.model.num_electrical_terminals()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.seg.setup_model(ctx);
        self.model.setup_model(ctx);
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        let n_elec = self.model.num_electrical_terminals();
        let stride = 2 * wpc;
        // A mis-wired instance used to panic here, which crosses the pyo3
        // boundary as a PanicException and tells the user nothing about which
        // element is wrong. Leave the device unconfigured instead; the terminal
        // count is checked against num_terminals() in build_devices_with_
        // footprints, which can name the element and fail cleanly.
        if terminals.len() < stride + n_elec || !(terminals.len() - n_elec).is_multiple_of(stride) {
            self.min_terminals = Some(stride + n_elec);
            return;
        }
        let optical_len = terminals.len() - n_elec;
        self.seg.setup_instance(&terminals[..optical_len], ctx);
        self.model.set_terminals(&terminals[optical_len..]);
        self.elec_nodes = terminals[optical_len..].to_vec();
        self.seg.set_control_nodes(&self.elec_nodes);
        self.dpert_dv = vec![(PerChannel::zero(), 0.0, PerChannel::zero()); self.elec_nodes.len()];
    }

    fn num_extra_nodes(&self) -> usize {
        self.seg.num_aux_branches() + self.model.num_internal_nodes()
    }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        self.seg.bind_branches(first_idx);
        self.model
            .bind_internal(first_idx + self.seg.num_aux_branches());
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        let lc = name.to_lowercase();
        // Route to both: geometry params land on the segment, drive params on
        // the model. A few (wl_ref) are consumed by both, so don't early-return.
        let s = self.seg.set_param(&lc, value);
        let m = self.model.set_param(&lc, value);
        s || m
    }

    fn lambda_routing(&self) -> Vec<(usize, usize)> {
        self.seg.lambda_routing()
    }

    fn set_resolved_lambda(&mut self, per_terminal: &[f64]) {
        self.seg.set_resolved_lambda(per_terminal);
    }

    fn eval(&mut self, x: &[f64], flags: EvalFlags, ctx: &SimContext) {
        let intensity = self.seg.channel_intensities(x);
        let pert = self.model.eval(x, &intensity, ctx);
        // Finite-difference the perturbation against each control voltage so the
        // segment can stamp ∂(optical output)/∂V. Step is relative to the local
        // voltage with an absolute floor, the usual compromise between
        // truncation and cancellation error.
        if !self.elec_nodes.is_empty() {
            let mut xp = x.to_vec();
            for (i, node) in self.elec_nodes.iter().enumerate() {
                let Some(j) = *node else {
                    self.dpert_dv[i] = (PerChannel::zero(), 0.0, PerChannel::zero());
                    continue;
                };
                let v = x[j];
                let h = 1e-6_f64.max(v.abs() * 1e-6);
                xp[j] = v + h;
                let pp = self.model.eval(&xp, &intensity, ctx);
                xp[j] = v;
                self.dpert_dv[i] = (
                    pp.dn_eff.diff_by(&pert.dn_eff, h),
                    (pp.dphi - pert.dphi) / h,
                    pp.dalpha_neper_m.diff_by(&pert.dalpha_neper_m, h),
                );
            }
            // The probing evals left the model's cached electrical state at the
            // perturbed point; restore it at x before anything stamps.
            let _ = self.model.eval(x, &intensity, ctx);
        }
        // Group-delay line: engaged in transient when the `waveguide_delay`
        // option is on and the segment has a finite group delay τ_g = L·n_g/c.
        // The lumped model delays the input envelope by τ_g and applies the
        // current transmission (which carries Δn_eff/Δφ/Δα from the drive at
        // time t). DC/AC and the default instantaneous path are unaffected; a
        // zero-length segment (e.g. fc_thermal_ps) has τ_g = 0 and stays
        // instantaneous regardless.
        let delay_active = flags.transient && ctx.waveguide_delay && self.seg.tau_g_s() > 0.0;
        self.seg.refresh_with_sens(
            x,
            pert.dn_eff.clone(),
            pert.dphi,
            pert.dalpha_neper_m.clone(),
            delay_active,
            ctx,
            &self.dpert_dv,
        );
    }

    fn load_residual(&self, b: &mut [f64]) {
        self.seg.stamp_residual(b);
        self.model.stamp_residual(b);
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        self.seg.stamp(mat);
        self.model.stamp(mat);
    }

    fn load_residual_tran(&self, b: &mut [f64], alpha: f64) {
        self.seg.stamp_residual(b);
        self.model.stamp_residual_tran(b, alpha);
    }

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, alpha: f64) {
        self.seg.stamp(mat);
        self.model.stamp_tran(mat, alpha);
    }

    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> {
        self.model.reactive_branches()
    }

    fn requested_max_timestep(&self) -> Option<f64> {
        self.seg.requested_max_timestep()
    }

    fn ac_stamps(&self, omega: f64) -> Vec<crate::device::AcStamp> {
        self.seg.ac_stamps(omega)
    }

    fn commit_timestep(&mut self, x: &[f64]) {
        self.seg.commit(x);
        self.model.commit(x);
    }

    /// The segment stamps true `∂(out)/∂V` cross-terms whenever its drive model
    /// supplies `∂perturbation/∂V`, and falls back to a frozen coefficient when
    /// it does not.  A `PhotonicActiveModel` that reports no derivatives is
    /// therefore invisible to the adjoint, so say so in exactly that case.
    fn frozen_jacobian_columns(&self) -> Vec<usize> {
        self.seg.unstamped_control_columns()
    }
}

// ────────────────────────────────────────────────────────────────────────
// Device constructors (defaults live here, where the dB→Neper helper lives)
// ────────────────────────────────────────────────────────────────────────

/// Build the `fc_pn_ps` device — a linear free-carrier PN-junction phase
/// shifter. SOI rib waveguide PN modulator section (R = 8 µm bent) defaults;
/// `dn_dv` follows from `V_pi_L = 0.015 V·m` (≈ 15 V·π at the 1-mm default
/// length). `pin_at_ref = false`: physical absolute propagation phase, so ring
/// resonances depend on L.
pub fn pn_phase_shifter() -> ActiveOpticalDevice {
    let seg = OpticalSegment::new(1e-3, 2.7654, 4.02, dB_per_cm_to_neper_per_m(20.0));
    ActiveOpticalDevice::new(seg, Box::new(PnDrive::new()))
}

/// Build the `fc_pn_ps_cap` device — the PN phase shifter plus a bias-dependent
/// depletion junction capacitance (and optional `da/dV` reverse-bias loss).
/// Same SOI-rib optics as `fc_pn_ps`.
pub fn pn_phase_shifter_cap() -> ActiveOpticalDevice {
    let seg = OpticalSegment::new(1e-3, 2.7654, 4.02, dB_per_cm_to_neper_per_m(20.0));
    ActiveOpticalDevice::new(seg, Box::new(PnDrive::with_depletion_cap()))
}

/// Build the `fc_pn_ps_inj` device — a forward-bias carrier-injection PN phase
/// shifter (Shockley diode + diffusion cap + exponential injection Δn/Δα).
pub fn pn_phase_shifter_inj() -> ActiveOpticalDevice {
    let seg = OpticalSegment::new(1e-3, 2.7654, 4.02, dB_per_cm_to_neper_per_m(1.0));
    ActiveOpticalDevice::new(seg, Box::new(Injection::new()))
}

/// Build the `fc_pn_th_ps_inj` device — the injection PN phase shifter with a
/// metal heater bolted on.
pub fn pn_thermal_phase_shifter_inj() -> ActiveOpticalDevice {
    let seg = OpticalSegment::new(1e-3, 2.7654, 4.02, dB_per_cm_to_neper_per_m(1.0));
    ActiveOpticalDevice::new(seg, Box::new(WithHeater::new(Box::new(Injection::new()))))
}

/// Build the `fc_thermal_ps` device — a metal-heater thermo-optic phase shifter
/// (no PN junction). A pure phase rotation `φ = π·P/P_π`, modelled on a
/// zero-length segment (no propagation phase, no loss).
pub fn thermal_phase_shifter() -> ActiveOpticalDevice {
    let seg = OpticalSegment::new(0.0, 2.445, 4.19, 0.0);
    ActiveOpticalDevice::new(seg, Box::new(Heater::new()))
}

/// Build the `fc_pn_ps_full` device — the L3 PN modulator (depletion +
/// injection + TPA + static self-heating + series resistance).
pub fn pn_phase_shifter_full() -> ActiveOpticalDevice {
    let seg = OpticalSegment::new(1e-3, 2.7654, 4.02, dB_per_cm_to_neper_per_m(1.0));
    ActiveOpticalDevice::new(seg, Box::new(FullPnDrive::new()))
}

/// Build the `fc_pn_th_ps_full` device — the L3 PN modulator with a metal
/// heater bolted on.
pub fn pn_thermal_phase_shifter_full() -> ActiveOpticalDevice {
    let seg = OpticalSegment::new(1e-3, 2.7654, 4.02, dB_per_cm_to_neper_per_m(1.0));
    ActiveOpticalDevice::new(seg, Box::new(WithHeater::new(Box::new(FullPnDrive::new()))))
}

/// Build the `fc_thermal_ps_rc` device — a metal-heater phase shifter with a
/// first-order thermal RC (filtered heater power). Zero-length segment.
pub fn thermal_rc_phase_shifter() -> ActiveOpticalDevice {
    let seg = OpticalSegment::new(0.0, 2.445, 4.19, 0.0);
    ActiveOpticalDevice::new(seg, Box::new(HeaterRc::new()))
}

/// Build the `fc_pn_th_ps` device — a PN phase shifter with a metal heater
/// bolted on (Δn_eff from V_pn and φ_th from Joule power sum). Terminal order:
/// optical bundle, then anode, cathode, heat_p, heat_n.
pub fn pn_thermal_phase_shifter() -> ActiveOpticalDevice {
    let seg = OpticalSegment::new(1e-3, 2.7654, 4.02, dB_per_cm_to_neper_per_m(20.0));
    ActiveOpticalDevice::new(seg, Box::new(WithHeater::new(Box::new(PnDrive::new()))))
}

/// Build the `fc_pn_th_ps_cap` device — depletion-mode PN (C_j(V) + reverse-bias
/// FCA) plus a metal heater. Calibrated defaults from a lateral-PN modulator
/// extraction (5e17/5e17, 300 K): distinct dn_dv, loss, and C_j from the plain
/// `fc_pn_ps_cap`.
pub fn pn_thermal_phase_shifter_cap() -> ActiveOpticalDevice {
    let seg = OpticalSegment::new(1e-3, 2.7654, 4.02, 29.78); // 2.59 dB/cm in Np/m
    let mut pn = PnDrive::with_depletion_cap();
    pn.set_param("dn_dv", 5.024e-5);
    pn.set_param("c_j0", 1.375e-13);
    pn.set_param("v_bi", 0.917);
    pn.set_param("da_dv", 7.83);
    // g_pn (1e-3) and m_j (0.5) already match the PnDrive defaults.
    ActiveOpticalDevice::new(seg, Box::new(WithHeater::new(Box::new(pn))))
}

// ────────────────────────────────────────────────────────────────────────
// Drive models
// ────────────────────────────────────────────────────────────────────────

/// Free-carrier PN-junction drive (the `fc_pn_ps` / `fc_pn_ps_cap` physics):
/// a single shared junction conductance `g_pn` between anode and cathode, a
/// linear electro-optic index change `Δn_eff = (dn/dV)·V_pn`, an optional
/// bias-dependent depletion capacitance `C_j(V_pn)` (active only when
/// `c_j0 > 0`), and an optional linear reverse-bias free-carrier loss
/// `Δα = (dα/dV)·max(0, −V_pn)`.
///
/// `c_j0 = 0` / `da_dv = 0` reproduce `fc_pn_ps` exactly (no cap, no loss
/// change); a non-zero `c_j0` makes it `fc_pn_ps_cap`. This is the LEVEL
/// collapse in miniature — one drive model, parameters select the electrical
/// sophistication. Forward-injection and full (TPA/self-heat) regimes are
/// separate models, not flags here (they are alternative / superset physics).
pub struct PnDrive {
    dn_dv: f64,
    g_pn: f64,
    /// Design wavelength for the `V_pi_L → dn_dv` conversion. Tracks the
    /// segment's `wl_ref_m` (both default to the band centre and both update on
    /// a `wl_ref` param), kept here so the conversion needs no segment handle.
    wl_ref_m: f64,
    // Optional depletion junction cap (inactive when c_j0 == 0).
    c_j0: f64,
    v_bi: f64,
    m_j: f64,
    // Optional linear reverse-bias FCA loss (Np/m per volt of reverse bias).
    da_dv: f64,
    // Per-eval cache: C_j(V_pn) at the current iterate.
    c_j_cached: f64,
    /// The junction's stored charge at the cached bias, `∫C dv`.
    ///
    /// This is the state the integrator advances, and it is *not* `C·V`. See
    /// `ReactiveBranchSpec::charge`. Carrying it also removes the need for
    /// `dC/dV`: the Jacobian of a charge branch is `α·dq/dv = α·C`.
    q_j_cached: f64,
    anode: NodeId,
    cathode: NodeId,
}

impl PnDrive {
    /// The `fc_pn_ps` drive: linear EO only (no junction cap, no loss change).
    pub fn new() -> Self {
        PnDrive {
            dn_dv: 1.55e-6 / (2.0 * 0.015),
            g_pn: 1e-3,
            wl_ref_m: 1.55e-6,
            c_j0: 0.0,
            v_bi: 0.7,
            m_j: 0.5,
            da_dv: 0.0,
            c_j_cached: 0.0,
            q_j_cached: 0.0,
            anode: None,
            cathode: None,
        }
    }

    /// The `fc_pn_ps_cap` drive: adds the depletion junction capacitance.
    pub fn with_depletion_cap() -> Self {
        PnDrive {
            c_j0: 20e-15,
            c_j_cached: 20e-15,
            q_j_cached: 0.0,
            ..Self::new()
        }
    }
}

impl Default for PnDrive {
    fn default() -> Self {
        Self::new()
    }
}

impl PhotonicActiveModel for PnDrive {
    fn num_electrical_terminals(&self) -> usize {
        2 // anode, cathode
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wl_ref_m = ctx.lambda_center_m;
    }

    fn set_terminals(&mut self, electrical: &[NodeId]) {
        self.anode = electrical[0];
        self.cathode = electrical[1];
    }

    fn eval(&mut self, x: &[f64], _intensity_w: &[f64], _ctx: &SimContext) -> OpticalPerturbation {
        let v_a = self.anode.map_or(0.0, |i| x[i]);
        let v_c = self.cathode.map_or(0.0, |i| x[i]);
        let v_pn = v_a - v_c;

        // Depletion C_j(V_pn) and its stored charge, with a linear tangent past
        // V_bi/2 to stay finite and keep the NR Jacobian smooth (only
        // meaningful when c_j0 > 0).
        if self.c_j0 > 0.0 {
            let (c, q) = junction_cap_and_charge(v_pn, self.c_j0, self.v_bi, self.m_j);
            self.c_j_cached = c;
            self.q_j_cached = q;
        }

        // Reverse-bias FCA loss: Δα = (dα/dV)·max(0, −V_pn).
        let v_rev = (-v_pn).max(0.0);
        OpticalPerturbation {
            dn_eff: (self.dn_dv * v_pn).into(),
            dphi: 0.0,
            dalpha_neper_m: (self.da_dv * v_rev).into(),
        }
    }

    fn stamp(&self, mat: &mut MnaMatrix) {
        stamp_resistor(mat, self.anode, self.cathode, self.g_pn);
    }

    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> {
        if self.c_j0 <= 0.0 {
            return Vec::new();
        }
        vec![ReactiveBranchSpec {
            kind: ReactiveKind::Capacitor,
            pos: self.anode,
            neg: self.cathode,
            value: self.c_j_cached,
            dvalue_dstate: 0.0,
            charge: Some(self.q_j_cached),
        }]
    }

    fn set_param(&mut self, name: &str, value: f64) -> bool {
        match name {
            "dn_dv" => {
                self.dn_dv = value;
                true
            }
            "g_pn" => {
                self.g_pn = value;
                true
            }
            "v_pi_l" => {
                // dn_dv such that 2π·L·dn_dv·Vπ/λ = π at V = Vπ ⇒
                // dn_dv = λ_ref / (2·V_pi_L).
                if value > 0.0 {
                    self.dn_dv = self.wl_ref_m / (2.0 * value);
                }
                true
            }
            "c_j0" => {
                self.c_j0 = value.max(0.0);
                true
            }
            "v_bi" => {
                self.v_bi = value.max(1e-3);
                true
            }
            "m_j" => {
                self.m_j = value.clamp(0.0, 0.99);
                true
            }
            "da_dv" => {
                self.da_dv = value;
                true
            }
            // wl_ref is also consumed by the segment; cache it for v_pi_l.
            "wl_ref_m" | "lambda_ref_m" => {
                self.wl_ref_m = value;
                false // let the segment record it as the authoritative copy
            }
            "wl_ref_nm" | "lambda_ref_nm" => {
                self.wl_ref_m = value * 1e-9;
                false
            }
            _ => false,
        }
    }
}

/// Forward-bias carrier-injection PN drive (the `fc_pn_ps_inj` physics).
/// Forward bias only (V_pn ∈ [0, ~0.8 V]); does not model depletion. Differs
/// from [`PnDrive`]:
///   - nonlinear Shockley diode `I = I_s·(exp(V/(n·V_T)) − 1)`, Norton-stamped
///     at the operating point (`g_d` Jacobian + `i_eq` residual);
///   - diffusion capacitance `C_d = τ_carrier·g_d` (replaces depletion C_j);
///   - exponential injection `Δn_eff = −K_inj·(e−1)` and `Δα = K_α·(e−1)`,
///     where `e = exp(V/(n·V_T))` (Soref–Bennett carrier-density form).
pub struct Injection {
    i_sat: f64,
    n_diode: f64,
    tau_carrier: f64,
    dn_dv_inj: f64,
    da_dv_inj: f64,
    anode: NodeId,
    cathode: NodeId,
    g_d_cached: f64,
    i_eq_cached: f64,
    c_d_cached: f64,
    /// `dC_d/dV` at the cached bias; see `ReactiveBranchSpec::dvalue_dstate`.
    dc_d_dv_cached: f64,
}

impl Injection {
    pub fn new() -> Self {
        Injection {
            i_sat: 1e-12,
            n_diode: 1.05,
            tau_carrier: 10e-9,
            dn_dv_inj: 1.311e-4,
            da_dv_inj: 150.0,
            anode: None,
            cathode: None,
            g_d_cached: 1e-9,
            i_eq_cached: 0.0,
            c_d_cached: 0.0,
            dc_d_dv_cached: 0.0,
        }
    }
}

impl Default for Injection {
    fn default() -> Self {
        Self::new()
    }
}

impl PhotonicActiveModel for Injection {
    fn num_electrical_terminals(&self) -> usize {
        2 // anode, cathode
    }

    fn set_terminals(&mut self, electrical: &[NodeId]) {
        self.anode = electrical[0];
        self.cathode = electrical[1];
    }

    fn eval(&mut self, x: &[f64], _intensity_w: &[f64], ctx: &SimContext) -> OpticalPerturbation {
        let v_a = self.anode.map_or(0.0, |i| x[i]);
        let v_c = self.cathode.map_or(0.0, |i| x[i]);
        let v_pn = v_a - v_c;
        let vt = ctx.vt() * self.n_diode;
        let arg = (v_pn / vt).clamp(-40.0, 40.0);
        let e = arg.exp();
        let i_diode = self.i_sat * (e - 1.0);
        self.g_d_cached = self.i_sat * e / vt;
        self.i_eq_cached = i_diode - self.g_d_cached * v_pn;
        self.c_d_cached = self.tau_carrier * self.g_d_cached;
        // g_d is Shockley, so dg_d/dv = g_d/vt (vt already carries n_diode) and
        // the diffusion cap inherits that slope exactly.
        self.dc_d_dv_cached = self.c_d_cached / vt;
        // Normalised injected carrier density (≥ 0; reverse bias contributes 0).
        let inj = (e - 1.0).max(0.0);
        OpticalPerturbation {
            dn_eff: (-self.dn_dv_inj * inj).into(), // more carriers → lower index
            dphi: 0.0,
            dalpha_neper_m: (self.da_dv_inj * inj).into(),
        }
    }

    fn stamp(&self, mat: &mut MnaMatrix) {
        stamp_resistor(mat, self.anode, self.cathode, self.g_d_cached);
    }

    fn stamp_residual(&self, b: &mut [f64]) {
        if let Some(a) = self.anode {
            b[a] -= self.i_eq_cached;
        }
        if let Some(c) = self.cathode {
            b[c] += self.i_eq_cached;
        }
    }

    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> {
        vec![ReactiveBranchSpec {
            kind: ReactiveKind::Capacitor,
            pos: self.anode,
            neg: self.cathode,
            // Still the `q = C·v` state, which is wrong for the depletion
            // half of this capacitance in the same way `PnDrive`'s was — the
            // device comes out faster than it is. Not converted with it because
            // the charge here is a function of the *internal* junction voltage,
            // which the series resistance separates from the branch voltage, so
            // the branch cannot report `q(v_branch)` without solving for it.
            value: self.c_d_cached,
            dvalue_dstate: self.dc_d_dv_cached,
            charge: None,
        }]
    }

    fn set_param(&mut self, name: &str, value: f64) -> bool {
        match name {
            "i_sat" | "is" => {
                self.i_sat = value;
                true
            }
            "n_diode" | "n" => {
                self.n_diode = value;
                true
            }
            "tau_carrier" | "tau" => {
                self.tau_carrier = value;
                true
            }
            "dn_dv_inj" | "dn_dv" => {
                self.dn_dv_inj = value;
                true
            }
            "da_dv_inj" | "da_dv" => {
                self.da_dv_inj = value;
                true
            }
            _ => false,
        }
    }
}

/// "Full" L3 PN drive (the `fc_pn_ps_full` physics): both bias regimes
/// (reverse depletion + forward injection), depletion + diffusion capacitance,
/// reverse/forward free-carrier absorption, two-photon absorption (TPA), static
/// thermal self-heating from absorbed optical power, and an optional ohmic
/// series resistance (the junction voltage is solved implicitly each iterate).
///
/// This is the model that exercises the optical→thermal/carrier **back-action**:
/// `eval` reads the segment's optical `intensity_w` to compute TPA loss and
/// self-heating Δn — the in-tree proof that the abstraction admits back-action
/// (no new stub device needed). It caches `length_m`/`alpha_neper_m` (the
/// segment geometry its self-heating needs) the way [`PnDrive`] caches `wl_ref`.
pub struct FullPnDrive {
    // Reverse / depletion.
    dn_dv_rev: f64,
    da_dv_rev: f64,
    c_j0: f64,
    v_bi: f64,
    m_j: f64,
    // Forward / injection.
    i_sat: f64,
    n_diode: f64,
    tau_carrier: f64,
    dn_dv_inj: f64,
    da_dv_inj: f64,
    // Current-parametrized injection: Δn = −dn_di·I_fwd, Δα = da_di·I_fwd
    // (per ampere of forward diode current). Physically equivalent to the
    // (e−1)-factor form above — I_fwd = i_sat·(e−1) — but decoupled from
    // (i_sat, n_diode), so a fit can pin the diode from the measured IV and
    // then fit these two as well-conditioned linear coefficients. Default 0
    // (off); use EITHER this or dn_dv_inj/da_dv_inj, not both.
    dn_di: f64,
    da_di: f64,
    // Series resistance.
    r_series: f64,
    // TPA + thermal.
    beta_tpa_m_per_w: f64,
    a_eff_m2: f64,
    r_th_k_per_w: f64,
    dn_dt: f64,
    // Cached segment geometry needed for the self-heating power budget.
    length_m: f64,
    alpha_neper_m: f64,
    anode: NodeId,
    cathode: NodeId,
    g_pn_cached: f64,
    i_eq_cached: f64,
    c_eff_cached: f64,
    /// `dC_eff/dV` at the cached bias; see `ReactiveBranchSpec::dvalue_dstate`.
    dc_eff_dv_cached: f64,
}

impl FullPnDrive {
    pub fn new() -> Self {
        FullPnDrive {
            dn_dv_rev: 5.024e-5,
            da_dv_rev: 7.83,
            c_j0: 1.375e-13,
            v_bi: 0.917,
            m_j: 0.5,
            i_sat: 1e-12,
            n_diode: 1.05,
            tau_carrier: 10e-9,
            dn_dv_inj: 1.311e-4,
            da_dv_inj: 150.0,
            dn_di: 0.0,
            da_di: 0.0,
            r_series: 0.0,
            beta_tpa_m_per_w: 7.9e-12,
            a_eff_m2: 1.257e-13,
            r_th_k_per_w: 0.0,
            dn_dt: 1.86e-4,
            length_m: 1e-3,
            alpha_neper_m: dB_per_cm_to_neper_per_m(1.0),
            anode: None,
            cathode: None,
            g_pn_cached: 1e-9,
            i_eq_cached: 0.0,
            c_eff_cached: 1.375e-13,
            dc_eff_dv_cached: 0.0,
        }
    }
}

impl Default for FullPnDrive {
    fn default() -> Self {
        Self::new()
    }
}

impl PhotonicActiveModel for FullPnDrive {
    fn num_electrical_terminals(&self) -> usize {
        2 // anode, cathode
    }

    fn set_terminals(&mut self, electrical: &[NodeId]) {
        self.anode = electrical[0];
        self.cathode = electrical[1];
    }

    fn eval(&mut self, x: &[f64], intensity_w: &[f64], ctx: &SimContext) -> OpticalPerturbation {
        let v_a = self.anode.map_or(0.0, |i| x[i]);
        let v_c = self.cathode.map_or(0.0, |i| x[i]);
        let v_pn = v_a - v_c;
        let vt = ctx.vt() * self.n_diode;

        // Junction voltage: implicit when R_s > 0 (solve V_j + R_s·I_d = V_pn).
        let v_junc = if self.r_series <= 0.0 {
            v_pn
        } else {
            let mut vj = v_pn;
            for _ in 0..50 {
                let arg = (vj / vt).clamp(-40.0, 40.0);
                let e = arg.exp();
                let id = self.i_sat * (e - 1.0);
                let f = vj + self.r_series * id - v_pn;
                let gd = self.i_sat * e / vt;
                let delta = f / (1.0 + self.r_series * gd);
                vj -= delta;
                if delta.abs() < 1e-12 {
                    break;
                }
            }
            vj
        };

        let arg = (v_junc / vt).clamp(-40.0, 40.0);
        let e = arg.exp();
        let i_diode = self.i_sat * (e - 1.0);
        let g_d = self.i_sat * e / vt;
        // Norton stamp accounting for series resistance.
        let g_eff = g_d / (1.0 + g_d * self.r_series);
        self.g_pn_cached = g_eff.max(1e-15);
        self.i_eq_cached = i_diode - g_eff * v_pn;

        // Depletion C_j(V_j) (linear past the V_bi/2 knee) + diffusion C_d.
        let (c_j_v, dc_j_dvj) = {
            let v_knee = 0.5 * self.v_bi;
            if v_junc < v_knee {
                let c = self.c_j0 / (1.0 - v_junc / self.v_bi).powf(self.m_j);
                (c, c * self.m_j / (self.v_bi - v_junc))
            } else {
                let c_knee = self.c_j0 / (1.0 - v_knee / self.v_bi).powf(self.m_j);
                let dc_dv = c_knee * self.m_j / (self.v_bi - v_knee);
                (c_knee + dc_dv * (v_junc - v_knee), dc_dv)
            }
        };
        self.c_eff_cached = c_j_v + self.tau_carrier * g_d;
        // Both terms are functions of the *junction* voltage, and the branch is
        // stamped across the terminals.  With `r_series = 0` those are the same
        // node pair and this is exact; with series resistance they differ by
        // `dv_junc/dv_pn = 1/(1 + g_d·r_series)`, which is applied here so the
        // derivative is against the branch's own state either way.
        let dvj_dv = 1.0 / (1.0 + g_d * self.r_series);
        self.dc_eff_dv_cached = (dc_j_dvj + self.tau_carrier * g_d / vt) * dvj_dv;

        // Optical back-action, over the whole bus and both directions
        // (`channel_intensities` counts both). A shared effect must be driven by
        // the shared total: reading channel 0 alone under-counted by 1/N with
        // every channel lit, and gave *exactly zero* whenever channel 0 was the
        // dark one.
        let total: f64 = intensity_w.iter().copied().map(|i| i.max(0.0)).sum();
        let inj = (e - 1.0).max(0.0);
        let i_fwd = i_diode.max(0.0);
        let v_rev = (-v_junc).max(0.0);
        // Loss common to every channel: reverse + forward free-carrier
        // absorption. Forward FCA has two equivalent parametrizations (see the
        // field docs).
        let dalpha_shared = self.da_dv_rev * v_rev + self.da_dv_inj * inj + self.da_di * i_fwd;

        // TPA is PER CHANNEL, because cross-TPA between distinct frequencies is
        // twice self-TPA:
        //     α_TPA,j = β/A_eff · (I_j + 2·Σ_{k≠j} I_k) = β/A_eff · (2·Σ_k I_k − I_j)
        // A single channel gives 2·I − I = I, so this reduces exactly to
        // self-TPA and the one-channel answer is unchanged. Before there was
        // anywhere to put a per-channel Δα this had to be the
        // no-cross-enhancement bound `β/A_eff·Σ_k I_k`, under-estimating the
        // loss by up to a factor of two on a fully loaded bus.
        let tpa_of =
            |i_j: f64| self.beta_tpa_m_per_w * (2.0 * total - i_j.max(0.0)) / self.a_eff_m2;

        // Absorbed power is what heats the device, and it is additive over
        // channels with each channel paying its own loss. One temperature, so
        // the resulting Δn is shared.
        let p_abs: f64 = intensity_w
            .iter()
            .map(|&i_j| {
                let i_j = i_j.max(0.0);
                (self.alpha_neper_m + dalpha_shared + tpa_of(i_j)) * self.length_m * i_j
            })
            .sum();
        let dn_self = self.dn_dt * self.r_th_k_per_w * p_abs;
        // All three index changes are λ-dependent ⇒ fold into dn_eff.
        let dn_eff = self.dn_dv_rev * v_junc - self.dn_dv_inj * inj - self.dn_di * i_fwd + dn_self;

        // One channel (or none) needs no vector: keep the scalar path allocation
        // free, since this runs per Newton iteration and again per FD probe.
        let dalpha = if intensity_w.len() <= 1 {
            PerChannel::Uniform(dalpha_shared + tpa_of(total))
        } else {
            PerChannel::Each(
                intensity_w
                    .iter()
                    .map(|&i_j| dalpha_shared + tpa_of(i_j))
                    .collect(),
            )
        };

        OpticalPerturbation {
            dn_eff: dn_eff.into(),
            dphi: 0.0,
            dalpha_neper_m: dalpha,
        }
    }

    fn stamp(&self, mat: &mut MnaMatrix) {
        stamp_resistor(mat, self.anode, self.cathode, self.g_pn_cached);
    }

    fn stamp_residual(&self, b: &mut [f64]) {
        if let Some(a) = self.anode {
            b[a] -= self.i_eq_cached;
        }
        if let Some(c) = self.cathode {
            b[c] += self.i_eq_cached;
        }
    }

    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> {
        vec![ReactiveBranchSpec {
            kind: ReactiveKind::Capacitor,
            pos: self.anode,
            neg: self.cathode,
            // Still the `q = C·v` state, which is wrong for the depletion
            // half of this capacitance in the same way `PnDrive`'s was — the
            // device comes out faster than it is. Not converted with it because
            // the charge here is a function of the *internal* junction voltage,
            // which the series resistance separates from the branch voltage, so
            // the branch cannot report `q(v_branch)` without solving for it.
            value: self.c_eff_cached,
            dvalue_dstate: self.dc_eff_dv_cached,
            charge: None,
        }]
    }

    fn set_param(&mut self, name: &str, value: f64) -> bool {
        match name {
            "dn_dv_rev" | "dn_dv" => {
                self.dn_dv_rev = value;
                true
            }
            "da_dv_rev" | "da_dv" => {
                self.da_dv_rev = value;
                true
            }
            "c_j0" => {
                self.c_j0 = value.max(0.0);
                true
            }
            "v_bi" => {
                self.v_bi = value.max(1e-3);
                true
            }
            "m_j" => {
                self.m_j = value.clamp(0.0, 0.99);
                true
            }
            "i_sat" | "is" => {
                self.i_sat = value;
                true
            }
            "n_diode" | "n" => {
                self.n_diode = value;
                true
            }
            "tau_carrier" | "tau" => {
                self.tau_carrier = value;
                true
            }
            "dn_dv_inj" => {
                self.dn_dv_inj = value;
                true
            }
            "da_dv_inj" => {
                self.da_dv_inj = value;
                true
            }
            "dn_di" | "dn_di_per_a" => {
                self.dn_di = value;
                true
            }
            "da_di" | "da_di_per_a" => {
                self.da_di = value;
                true
            }
            "beta_tpa" | "beta_tpa_m_per_w" => {
                self.beta_tpa_m_per_w = value;
                true
            }
            "a_eff_m2" => {
                self.a_eff_m2 = value.max(1e-20);
                true
            }
            "a_eff_um2" => {
                self.a_eff_m2 = value.max(1e-8) * 1e-12;
                true
            }
            "r_th" | "r_th_k_per_w" => {
                self.r_th_k_per_w = value.max(0.0);
                true
            }
            "dn_dt" => {
                self.dn_dt = value;
                true
            }
            "r_series" => {
                self.r_series = value.max(0.0);
                true
            }
            // Geometry the self-heating needs; segment stays authoritative.
            "l_um" => {
                self.length_m = value * 1e-6;
                false
            }
            "l_m" | "length" => {
                self.length_m = value;
                false
            }
            "alpha_db_cm" => {
                self.alpha_neper_m = dB_per_cm_to_neper_per_m(value);
                false
            }
            _ => false,
        }
    }
}

/// Resistive metal-heater drive (thermo-optic). A heater resistor `1/R_heater`
/// between `heat_p`/`heat_n`; Joule power `P = V²/R` produces a calibrated,
/// wavelength-independent phase `φ_th = π·P/P_π` applied to every channel — the
/// `fc_thermal_ps` physics. No optical loss, no index dispersion (it is a pure
/// phase rotation; modelled on a zero-length segment when standalone).
pub struct Heater {
    r_heater: f64,
    p_pi: f64,
    heat_p: NodeId,
    heat_n: NodeId,
}

impl Heater {
    pub fn new() -> Self {
        Heater {
            r_heater: 1000.0,
            p_pi: 10e-3,
            heat_p: None,
            heat_n: None,
        }
    }
}

impl Default for Heater {
    fn default() -> Self {
        Self::new()
    }
}

impl PhotonicActiveModel for Heater {
    fn num_electrical_terminals(&self) -> usize {
        2 // heat_p, heat_n
    }

    fn set_terminals(&mut self, electrical: &[NodeId]) {
        self.heat_p = electrical[0];
        self.heat_n = electrical[1];
    }

    fn eval(&mut self, x: &[f64], _intensity_w: &[f64], _ctx: &SimContext) -> OpticalPerturbation {
        let v = self.heat_p.map_or(0.0, |i| x[i]) - self.heat_n.map_or(0.0, |i| x[i]);
        let p = v * v / self.r_heater;
        OpticalPerturbation {
            dn_eff: PerChannel::zero(),
            dphi: std::f64::consts::PI * p / self.p_pi,
            dalpha_neper_m: PerChannel::zero(),
        }
    }

    fn stamp(&self, mat: &mut MnaMatrix) {
        stamp_resistor(mat, self.heat_p, self.heat_n, 1.0 / self.r_heater);
    }

    fn set_param(&mut self, name: &str, value: f64) -> bool {
        match name {
            "r_heater" | "r" => {
                self.r_heater = value;
                true
            }
            // `p_pi_th` is the alias the combined PN+thermal devices used.
            "p_pi" | "p_pi_w" | "p_pi_th" => {
                self.p_pi = value;
                true
            }
            _ => false,
        }
    }
}

/// Metal-heater drive with a first-order thermal RC (the `fc_thermal_ps_rc`
/// physics). The optical phase tracks the *filtered* heater power `T(t)` rather
/// than the instantaneous Joule dissipation: `dT/dt = (P − T)/τ_th`, with `T` a
/// device-owned state on one internal MNA row. At steady state `T = P`, so the
/// phase reduces to the L1 `φ = π·P/P_π`. This is the "path B" pattern (device
/// stamps its own discretised state equation), so it overrides `stamp_tran` /
/// `stamp_residual_tran` for the BE form.
pub struct HeaterRc {
    r_heater: f64,
    p_pi: f64,
    tau_th: f64,
    heat_p: NodeId,
    heat_n: NodeId,
    t_state_idx: Option<usize>,
    t_old: f64,
    v_h_op: f64,
}

impl HeaterRc {
    pub fn new() -> Self {
        HeaterRc {
            r_heater: 1000.0,
            p_pi: 10e-3,
            tau_th: 10e-6,
            heat_p: None,
            heat_n: None,
            t_state_idx: None,
            t_old: 0.0,
            v_h_op: 0.0,
        }
    }

    fn stamp_heater_resistor(&self, mat: &mut MnaMatrix) {
        stamp_resistor(mat, self.heat_p, self.heat_n, 1.0 / self.r_heater);
    }
}

impl Default for HeaterRc {
    fn default() -> Self {
        Self::new()
    }
}

impl PhotonicActiveModel for HeaterRc {
    fn num_electrical_terminals(&self) -> usize {
        2 // heat_p, heat_n
    }

    fn num_internal_nodes(&self) -> usize {
        1 // the T(t) state row
    }

    fn set_terminals(&mut self, electrical: &[NodeId]) {
        self.heat_p = electrical[0];
        self.heat_n = electrical[1];
    }

    fn bind_internal(&mut self, first_idx: usize) {
        self.t_state_idx = Some(first_idx);
    }

    fn eval(&mut self, x: &[f64], _intensity_w: &[f64], _ctx: &SimContext) -> OpticalPerturbation {
        self.v_h_op = self.heat_p.map_or(0.0, |i| x[i]) - self.heat_n.map_or(0.0, |i| x[i]);
        // Phase driven by the filtered heater power state T (power-equivalent W).
        let t = self.t_state_idx.map_or(0.0, |i| x[i]);
        OpticalPerturbation {
            dn_eff: PerChannel::zero(),
            dphi: std::f64::consts::PI * t / self.p_pi,
            dalpha_neper_m: PerChannel::zero(),
        }
    }

    fn stamp(&self, mat: &mut MnaMatrix) {
        self.stamp_heater_resistor(mat);
        // DC state row: T − P_lin(V_h) = 0  (T = P at steady state).
        if let Some(t) = self.t_state_idx {
            mat.a[t][t] += 1.0;
            let two_vop_over_r = 2.0 * self.v_h_op / self.r_heater;
            if let Some(hp) = self.heat_p {
                mat.a[t][hp] -= two_vop_over_r;
            }
            if let Some(hn) = self.heat_n {
                mat.a[t][hn] += two_vop_over_r;
            }
        }
    }

    fn stamp_tran(&self, mat: &mut MnaMatrix, alpha: f64) {
        self.stamp_heater_resistor(mat);
        // BE state row: T·(α + 1/τ) − 2·V_h_op·V_h/(R·τ) = …
        if let Some(t) = self.t_state_idx {
            let inv_tau = 1.0 / self.tau_th;
            mat.a[t][t] += alpha + inv_tau;
            let two_vop_over_r = 2.0 * self.v_h_op / self.r_heater;
            if let Some(hp) = self.heat_p {
                mat.a[t][hp] -= two_vop_over_r * inv_tau;
            }
            if let Some(hn) = self.heat_n {
                mat.a[t][hn] += two_vop_over_r * inv_tau;
            }
        }
    }

    fn stamp_residual(&self, b: &mut [f64]) {
        if let Some(t) = self.t_state_idx {
            b[t] -= self.v_h_op * self.v_h_op / self.r_heater;
        }
    }

    fn stamp_residual_tran(&self, b: &mut [f64], alpha: f64) {
        if let Some(t) = self.t_state_idx {
            let inv_tau = 1.0 / self.tau_th;
            let p_op = self.v_h_op * self.v_h_op / self.r_heater;
            b[t] += self.t_old * alpha - p_op * inv_tau;
        }
    }

    fn commit(&mut self, x: &[f64]) {
        if let Some(t) = self.t_state_idx {
            self.t_old = x[t];
        }
    }

    fn set_param(&mut self, name: &str, value: f64) -> bool {
        match name {
            "r_heater" | "r" => {
                self.r_heater = value;
                true
            }
            "p_pi" | "p_pi_w" => {
                self.p_pi = value;
                true
            }
            "tau_th" | "tau" => {
                self.tau_th = value.max(1e-30);
                true
            }
            _ => false,
        }
    }
}

/// Compose any drive model with a metal heater bolted on. The heater's
/// `heat_p`/`heat_n` terminals follow the inner model's electrical terminals;
/// Δn/Δφ/Δα sum, electrical stamps and reactive branches concatenate. This is
/// the orthogonal "+ optional heater" axis — `WithHeater::new(PnDrive::new())`
/// is `fc_pn_th_ps`, `WithHeater::new(PnDrive::with_depletion_cap())` is
/// `fc_pn_th_ps_cap`, etc. — applicable to any future drive (Pockels, …).
pub struct WithHeater {
    inner: Box<dyn PhotonicActiveModel>,
    heater: Heater,
}

impl WithHeater {
    pub fn new(inner: Box<dyn PhotonicActiveModel>) -> Self {
        WithHeater {
            inner,
            heater: Heater::new(),
        }
    }
}

impl PhotonicActiveModel for WithHeater {
    fn num_electrical_terminals(&self) -> usize {
        self.inner.num_electrical_terminals() + self.heater.num_electrical_terminals()
    }

    fn num_internal_nodes(&self) -> usize {
        self.inner.num_internal_nodes() + self.heater.num_internal_nodes()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.inner.setup_model(ctx);
        self.heater.setup_model(ctx);
    }

    fn set_terminals(&mut self, electrical: &[NodeId]) {
        let ni = self.inner.num_electrical_terminals();
        self.inner.set_terminals(&electrical[..ni]);
        self.heater.set_terminals(&electrical[ni..]);
    }

    fn bind_internal(&mut self, first_idx: usize) {
        self.inner.bind_internal(first_idx);
        self.heater
            .bind_internal(first_idx + self.inner.num_internal_nodes());
    }

    fn eval(&mut self, x: &[f64], intensity_w: &[f64], ctx: &SimContext) -> OpticalPerturbation {
        let a = self.inner.eval(x, intensity_w, ctx);
        let b = self.heater.eval(x, intensity_w, ctx);
        OpticalPerturbation {
            dn_eff: a.dn_eff.add(&b.dn_eff),
            dphi: a.dphi + b.dphi,
            dalpha_neper_m: a.dalpha_neper_m.add(&b.dalpha_neper_m),
        }
    }

    fn stamp(&self, mat: &mut MnaMatrix) {
        self.inner.stamp(mat);
        self.heater.stamp(mat);
    }

    fn stamp_tran(&self, mat: &mut MnaMatrix, alpha: f64) {
        self.inner.stamp_tran(mat, alpha);
        self.heater.stamp_tran(mat, alpha);
    }

    fn stamp_residual(&self, b: &mut [f64]) {
        self.inner.stamp_residual(b);
        self.heater.stamp_residual(b);
    }

    fn stamp_residual_tran(&self, b: &mut [f64], alpha: f64) {
        self.inner.stamp_residual_tran(b, alpha);
        self.heater.stamp_residual_tran(b, alpha);
    }

    fn reactive_branches(&self) -> Vec<ReactiveBranchSpec> {
        let mut v = self.inner.reactive_branches();
        v.extend(self.heater.reactive_branches());
        v
    }

    fn commit(&mut self, x: &[f64]) {
        self.inner.commit(x);
        self.heater.commit(x);
    }

    fn set_param(&mut self, name: &str, value: f64) -> bool {
        let i = self.inner.set_param(name, value);
        let h = self.heater.set_param(name, value);
        i || h
    }
}

// ────────────────────────────────────────────────────────────────────────
// Tier-1: declarative expression-driven model (no recompile)
// ────────────────────────────────────────────────────────────────────────

/// Variable environment for a device constitutive map: `V` = junction/terminal
/// bias, `T` = temperature (K), `lambda` = centre wavelength (m). Node/branch
/// references are not meaningful here and read 0.
struct VarCtx {
    v: f64,
    t: f64,
    lambda: f64,
}

impl EvalContext for VarCtx {
    fn node_voltage(&self, _: &str) -> f64 {
        0.0
    }
    fn branch_current(&self, _: &str) -> f64 {
        0.0
    }
    fn time(&self) -> f64 {
        0.0
    }
    fn variable(&self, name: &str) -> f64 {
        match name {
            "v" | "vpn" => self.v,
            "t" | "temp" | "temperature" => self.t,
            "lambda" | "wl" | "lam" => self.lambda,
            _ => 0.0,
        }
    }
}

/// A drive whose constitutive map is **declarative** — parsed expressions over
/// `(V, T, lambda)` rather than hard-coded Rust. This is the Tier-1
/// runtime-loadable model: a designer writes
/// `.model myps fc_phase_shifter_expr dneff="-3.1e-5*V - 1.2e-5*V*V"
/// dalpha="8.0" g_pn=1m`, re-runs, and gets new physics with no recompile. The
/// expressions are parsed once at setup and evaluated per NR-iterate.
///
/// Covers the closed-form-map case — any drive expressible as `Δn(V)` and
/// `Δα(V)`, which is most PN, thermal and Pockels maps. Stateful physics
/// (carrier rate equations, lookup tables) still needs a compiled
/// `PhotonicActiveModel`, because an expression has nowhere to keep state.
pub struct ExprDrive {
    dneff: Option<Expr>,
    dalpha: Option<Expr>,
    g_pn: f64,
    anode: NodeId,
    cathode: NodeId,
    lambda_center_m: f64,
}

impl ExprDrive {
    /// Build from already-parsed constitutive expressions.
    pub fn new(dneff: Option<Expr>, dalpha: Option<Expr>, g_pn: f64) -> Self {
        ExprDrive {
            dneff,
            dalpha,
            g_pn,
            anode: None,
            cathode: None,
            lambda_center_m: 1.55e-6,
        }
    }
}

impl PhotonicActiveModel for ExprDrive {
    fn num_electrical_terminals(&self) -> usize {
        2 // anode, cathode
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.lambda_center_m = ctx.lambda_center_m;
    }

    fn set_terminals(&mut self, electrical: &[NodeId]) {
        self.anode = electrical[0];
        self.cathode = electrical[1];
    }

    fn eval(&mut self, x: &[f64], _intensity_w: &[f64], ctx: &SimContext) -> OpticalPerturbation {
        let v = self.anode.map_or(0.0, |i| x[i]) - self.cathode.map_or(0.0, |i| x[i]);
        let env = VarCtx {
            v,
            t: ctx.temperature,
            lambda: self.lambda_center_m,
        };
        OpticalPerturbation {
            dn_eff: self.dneff.as_ref().map_or(0.0, |e| e.eval(&env)).into(),
            dphi: 0.0,
            dalpha_neper_m: self.dalpha.as_ref().map_or(0.0, |e| e.eval(&env)).into(),
        }
    }

    fn stamp(&self, mat: &mut MnaMatrix) {
        stamp_resistor(mat, self.anode, self.cathode, self.g_pn);
    }

    fn set_param(&mut self, name: &str, value: f64) -> bool {
        match name {
            "g_pn" => {
                self.g_pn = value;
                true
            }
            _ => false,
        }
    }
}

/// Build an `fc_phase_shifter_expr` device from parsed constitutive expressions
/// (`dneff`, `dalpha`) and a junction conductance. Optics default to the SOI-rib
/// PN baseline and are overridden by the card's numeric params.
pub fn expr_phase_shifter(
    dneff: Option<Expr>,
    dalpha: Option<Expr>,
    g_pn: f64,
) -> ActiveOpticalDevice {
    let seg = OpticalSegment::new(1e-3, 2.7654, 4.02, dB_per_cm_to_neper_per_m(1.0));
    ActiveOpticalDevice::new(seg, Box::new(ExprDrive::new(dneff, dalpha, g_pn)))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Group delay on a phase shifter: with `waveguide_delay` on, the output
    /// envelope lags the input by τ_g = L·n_g/c (the segment's `DelayLine`);
    /// with it off, the path is instantaneous (homogeneous residual).
    #[test]
    fn pn_ps_group_delay_engages_under_option() {
        // pn_phase_shifter defaults: L = 1 mm, n_g = 4.02.
        let tau = 1e-3 * 4.02 / 299_792_458.0;
        let ctx = |delay: bool, t: f64| SimContext {
            waveguide_delay: delay,
            time_s: t,
            ..SimContext::default()
        };
        // 8 terminals: in re/im/λ, out re/im/λ, anode, cathode; branches at 8.
        let terminals: Vec<NodeId> = (0..8).map(Some).collect();
        let build = |delay: bool| {
            let mut d = pn_phase_shifter();
            d.setup_model(&ctx(delay, 0.0));
            d.setup_instance(&terminals, &ctx(delay, 0.0));
            d.bind_extra_nodes(8);
            d
        };

        // ── Delay ON: drive in_re=1 for t=0 and t=τ, then query at t=2τ where
        // the delayed input (t−τ = τ) is the recorded 1.0. ──
        let mut on = build(true);
        let mut x_hi = vec![0.0; 11];
        x_hi[0] = 1.0; // in_re
        on.eval(&x_hi, EvalFlags::tran(), &ctx(true, 0.0));
        on.commit_timestep(&x_hi);
        on.eval(&x_hi, EvalFlags::tran(), &ctx(true, tau));
        on.commit_timestep(&x_hi);
        let x_lo = vec![0.0; 11]; // current input zero; output reflects delayed past
        on.eval(&x_lo, EvalFlags::tran(), &ctx(true, 2.0 * tau));
        let mut b = vec![0.0; 11];
        on.load_residual(&mut b);
        let out_mag = (b[8] * b[8] + b[9] * b[9]).sqrt();
        // The delayed input was 1.0; the output is it attenuated by the
        // propagation loss over L: t_amp = exp(-α·L/2), α from 20 dB/cm.
        let t_amp = (-dB_per_cm_to_neper_per_m(20.0) * 1e-3 / 2.0).exp();
        assert!(
            (out_mag - t_amp).abs() < 1e-6,
            "delay-on: output should be the delayed input × t_amp ({t_amp:.4}); got {out_mag:.4}"
        );

        // ── Delay OFF: instantaneous path — output lives in the Jacobian, the
        // residual is homogeneous. ──
        let mut off = build(false);
        off.eval(&x_hi, EvalFlags::tran(), &ctx(false, 0.0));
        off.commit_timestep(&x_hi);
        off.eval(&x_lo, EvalFlags::tran(), &ctx(false, 2.0 * tau));
        let mut b_off = vec![0.0; 11];
        off.load_residual(&mut b_off);
        assert!(
            b_off.iter().all(|v| v.abs() < 1e-12),
            "delay-off: residual must be homogeneous, got {b_off:?}"
        );
    }

    #[test]
    fn pn_drive_dn_eff_is_linear_in_vpn() {
        let mut m = PnDrive::new();
        m.dn_dv = 1e-4;
        m.set_terminals(&[Some(0), Some(1)]); // anode=node0, cathode=node1
        let x = [0.7, 0.2]; // V_pn = 0.5
        let p = m.eval(&x, &[], &SimContext::default());
        assert!(
            (p.dn_eff.at(0) - 1e-4 * 0.5).abs() < 1e-18,
            "dn_eff={:?}",
            p.dn_eff
        );
        assert_eq!(p.dalpha_neper_m, PerChannel::zero());
    }

    #[test]
    fn v_pi_l_sets_dn_dv_from_wl_ref() {
        let mut m = PnDrive::new();
        m.wl_ref_m = 1.55e-6;
        assert!(m.set_param("v_pi_l", 0.015));
        // dn_dv = λ_ref / (2·V_pi_L).
        assert!((m.dn_dv - 1.55e-6 / (2.0 * 0.015)).abs() < 1e-20);
    }
}
