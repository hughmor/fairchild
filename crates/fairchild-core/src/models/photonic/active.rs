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
//! See `_notes/optical_abstraction_design.md` for the full design and the
//! future device classes this admits.

use super::segment::OpticalSegment;
use super::{dB_per_cm_to_neper_per_m, stamp_resistor};
use crate::device::{Device, EvalFlags, NodeId, ReactiveBranchSpec, ReactiveKind, SimContext};
use crate::mna::MnaMatrix;

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
#[derive(Clone, Copy, Debug, Default)]
pub struct OpticalPerturbation {
    pub dn_eff: f64,
    pub dphi: f64,
    pub dalpha_neper_m: f64,
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
}

impl ActiveOpticalDevice {
    pub fn new(seg: OpticalSegment, model: Box<dyn PhotonicActiveModel>) -> Self {
        ActiveOpticalDevice { seg, model }
    }
}

impl Device for ActiveOpticalDevice {
    fn num_terminals(&self) -> usize {
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
        assert!(
            terminals.len() >= stride + n_elec && (terminals.len() - n_elec).is_multiple_of(stride),
            "ActiveOpticalDevice: terminal count must be {stride}·N + {n_elec} (wpc={wpc}); got {}",
            terminals.len()
        );
        let optical_len = terminals.len() - n_elec;
        self.seg.setup_instance(&terminals[..optical_len], ctx);
        self.model.set_terminals(&terminals[optical_len..]);
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

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, ctx: &SimContext) {
        let intensity = self.seg.channel_intensities(x);
        let pert = self.model.eval(x, &intensity, ctx);
        // Delay-on-phase-shifters is a future step; today's classes are
        // instantaneous, so the segment runs with the delay line disengaged.
        self.seg
            .refresh(x, pert.dn_eff, pert.dphi, pert.dalpha_neper_m, false, ctx);
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

    fn commit_timestep(&mut self, x: &[f64]) {
        self.seg.commit(x);
        self.model.commit(x);
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
            anode: None,
            cathode: None,
        }
    }

    /// The `fc_pn_ps_cap` drive: adds the depletion junction capacitance.
    pub fn with_depletion_cap() -> Self {
        PnDrive {
            c_j0: 20e-15,
            c_j_cached: 20e-15,
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

        // Depletion C_j(V_pn) with a linear tangent past V_bi/2 to stay finite
        // and keep the NR Jacobian smooth (only meaningful when c_j0 > 0).
        if self.c_j0 > 0.0 {
            let v_knee = 0.5 * self.v_bi;
            self.c_j_cached = if v_pn < v_knee {
                self.c_j0 / (1.0 - v_pn / self.v_bi).powf(self.m_j)
            } else {
                let c_knee = self.c_j0 / (1.0 - v_knee / self.v_bi).powf(self.m_j);
                let dc_dv = c_knee * self.m_j / (self.v_bi - v_knee);
                c_knee + dc_dv * (v_pn - v_knee)
            };
        }

        // Reverse-bias FCA loss: Δα = (dα/dV)·max(0, −V_pn).
        let v_rev = (-v_pn).max(0.0);
        OpticalPerturbation {
            dn_eff: self.dn_dv * v_pn,
            dphi: 0.0,
            dalpha_neper_m: self.da_dv * v_rev,
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
        // Normalised injected carrier density (≥ 0; reverse bias contributes 0).
        let inj = (e - 1.0).max(0.0);
        OpticalPerturbation {
            dn_eff: -self.dn_dv_inj * inj, // more carriers → lower index
            dphi: 0.0,
            dalpha_neper_m: self.da_dv_inj * inj,
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
            value: self.c_d_cached,
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
        let c_j_v = {
            let v_knee = 0.5 * self.v_bi;
            if v_junc < v_knee {
                self.c_j0 / (1.0 - v_junc / self.v_bi).powf(self.m_j)
            } else {
                let c_knee = self.c_j0 / (1.0 - v_knee / self.v_bi).powf(self.m_j);
                let dc_dv = c_knee * self.m_j / (self.v_bi - v_knee);
                c_knee + dc_dv * (v_junc - v_knee)
            }
        };
        self.c_eff_cached = c_j_v + self.tau_carrier * g_d;

        // Optical intensity (channel 0) drives TPA + self-heating back-action.
        let intensity = intensity_w.first().copied().unwrap_or(0.0).max(0.0);
        let inj = (e - 1.0).max(0.0);
        let v_rev = (-v_junc).max(0.0);
        let alpha_tpa = self.beta_tpa_m_per_w * intensity / self.a_eff_m2;
        // Extra loss beyond the segment's base α: reverse + forward FCA + TPA.
        let dalpha = self.da_dv_rev * v_rev + self.da_dv_inj * inj + alpha_tpa;
        // Absorbed power uses the TOTAL loss (segment base α + dalpha).
        let alpha_total = self.alpha_neper_m + dalpha;
        let p_abs = alpha_total * self.length_m * intensity;
        let dn_self = self.dn_dt * self.r_th_k_per_w * p_abs;
        // All three index changes are λ-dependent ⇒ fold into dn_eff.
        let dn_eff = self.dn_dv_rev * v_junc - self.dn_dv_inj * inj + dn_self;

        OpticalPerturbation {
            dn_eff,
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
            value: self.c_eff_cached,
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
            dn_eff: 0.0,
            dphi: std::f64::consts::PI * p / self.p_pi,
            dalpha_neper_m: 0.0,
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
            dn_eff: 0.0,
            dphi: std::f64::consts::PI * t / self.p_pi,
            dalpha_neper_m: 0.0,
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
            dn_eff: a.dn_eff + b.dn_eff,
            dphi: a.dphi + b.dphi,
            dalpha_neper_m: a.dalpha_neper_m + b.dalpha_neper_m,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pn_drive_dn_eff_is_linear_in_vpn() {
        let mut m = PnDrive::new();
        m.dn_dv = 1e-4;
        m.set_terminals(&[Some(0), Some(1)]); // anode=node0, cathode=node1
        let x = [0.7, 0.2]; // V_pn = 0.5
        let p = m.eval(&x, &[], &SimContext::default());
        assert!((p.dn_eff - 1e-4 * 0.5).abs() < 1e-18, "dn_eff={}", p.dn_eff);
        assert_eq!(p.dalpha_neper_m, 0.0);
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
