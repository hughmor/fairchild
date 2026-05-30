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
use crate::device::{Device, EvalFlags, NodeId, ReactiveBranchSpec, SimContext};
use crate::mna::MnaMatrix;

/// Per-segment optical perturbation produced by a [`PhotonicActiveModel`] at the
/// current operating point: an effective-index change `Δn_eff` (added to the
/// segment index) and an excess loss `Δα` (Neper/m, added to the propagation
/// loss). Both are mechanism-agnostic — the segment never asks how they were
/// produced — so any drive physics reduces to this pair.
#[derive(Clone, Copy, Debug, Default)]
pub struct OpticalPerturbation {
    pub dn_eff: f64,
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

    /// Stamp electrical RHS contributions (history sources, fixed currents).
    /// Default no-op.
    fn stamp_residual(&self, _b: &mut [f64]) {}

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
            .refresh(x, pert.dn_eff, pert.dalpha_neper_m, false, ctx);
    }

    fn load_residual(&self, b: &mut [f64]) {
        self.seg.stamp_residual(b);
        self.model.stamp_residual(b);
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        self.seg.stamp(mat);
        self.model.stamp(mat);
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
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

// ────────────────────────────────────────────────────────────────────────
// Drive models
// ────────────────────────────────────────────────────────────────────────

/// Linear free-carrier PN-junction drive (the `fc_pn_ps` physics): a single
/// shared junction conductance `g_pn` between anode and cathode, and a linear
/// electro-optic index change `Δn_eff = (dn/dV)·V_pn`. No loss change, no
/// junction capacitance (those are added by richer drive models / flags).
pub struct PnDrive {
    dn_dv: f64,
    g_pn: f64,
    /// Design wavelength for the `V_pi_L → dn_dv` conversion. Tracks the
    /// segment's `wl_ref_m` (both default to the band centre and both update on
    /// a `wl_ref` param), kept here so the conversion needs no segment handle.
    wl_ref_m: f64,
    anode: NodeId,
    cathode: NodeId,
}

impl PnDrive {
    pub fn new() -> Self {
        PnDrive {
            dn_dv: 1.55e-6 / (2.0 * 0.015),
            g_pn: 1e-3,
            wl_ref_m: 1.55e-6,
            anode: None,
            cathode: None,
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
        OpticalPerturbation {
            dn_eff: self.dn_dv * v_pn,
            dalpha_neper_m: 0.0,
        }
    }

    fn stamp(&self, mat: &mut MnaMatrix) {
        stamp_resistor(mat, self.anode, self.cathode, self.g_pn);
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
                // dn_dv = λ_ref / (2·Vπ·L) — but Vπ·L is the product, so
                // dn_dv = λ_ref / (2·V_pi_L).
                if value > 0.0 {
                    self.dn_dv = self.wl_ref_m / (2.0 * value);
                }
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
