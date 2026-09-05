//! `fc_tw_ps` — travelling-wave phase shifter.
//!
//! A lumped phase shifter says the whole electrode is at one voltage. That is
//! true while the device is short against the RF wavelength and false exactly
//! where a travelling-wave modulator earns its name. This device is the
//! interleaved ladder instead: `N` optical slices, `N` electrode sections, and
//! slice `i` driven by node `i` of the electrode.
//!
//! ```text
//!   rf_in ──[T]──┬──[T]──┬── … ──┬──[T]── rf_out
//!                │       │       │
//!   opt_in ─[seg]┴─[seg]─┴─ … ───┴─[seg]─ opt_out
//! ```
//!
//! Velocity mismatch, termination ripple, and the bandwidth collapse when the
//! RF is launched against the light are **not modelled here**. They emerge from
//! the topology, because the RF and the optical envelope accumulate delay at
//! different rates down the same ladder. `tests/native/travelling_wave_ladder.rs`
//! pins each of them against its closed form on a hand-written ladder, which is
//! the same construction this device builds for you.
//!
//! # Why `N` is not a parameter
//!
//! The requirement is that the RF voltage be roughly constant across one slice,
//! so `Δz ≪ λ_RF = c/(f_max·n_m)`. Asking a user for a slice count is asking
//! them to do that arithmetic, and to redo it whenever they change the length or
//! the electrode. So the card takes `f_max` — the top of the band they care
//! about — and the device solves for
//!
//! ```text
//!   N = ceil(slices_per_wave · L · n_m · f_max / c)
//! ```
//!
//! with `slices_per_wave` defaulting to 10. Convergence is `O(Δz²)`, so the
//! honest check is to raise `slices_per_wave` until the answer stops moving;
//! `n_slices` overrides the count outright for exactly that sweep.
//!
//! # What it does not model
//!
//! * **RF loss.** The electrode sections are lossless (`T`), so there is no
//!   conductor or dielectric loss and no skin effect. A real electrode rolls
//!   off from `√f` attenuation as well as from walk-off, and this device will
//!   be optimistic above the frequency where that dominates.
//! * **A junction capacitance per slice.** The electrode is unloaded, so its
//!   `n_m` and `Z0` are whatever the card says rather than the loaded values a
//!   real segmented electrode has. Give the *loaded* numbers.
//! * **Bidirectional light.** `wpc = 5` is refused: the backward wave would
//!   need its own ladder, and leaving it undriven would be worse.

use super::segment::{OpticalSegment, PerChannel};
use crate::device::{AcStamp, Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;
use crate::models::tline::NativeTLine;

/// Speed of light, as everywhere else in this module.
use super::C0;

pub struct NativeTwPhaseShifter {
    // ── geometry and electrode ───────────────────────────────────────────────
    length_m: f64,
    n_m: f64,
    z0: f64,
    f_max: f64,
    slices_per_wave: f64,
    n_slices_override: Option<usize>,
    /// Optical index change per volt, from `v_pi_l` (V·m) at the reference
    /// wavelength: a differential `V_pi` over length `L` is a π phase shift.
    dn_dv: f64,
    v_pi_l: f64,
    /// Geometry handed to every slice. Held as a template because the slices do
    /// not exist until `validate`, which is the first point at which every
    /// parameter has been applied.
    template: OpticalSegment,

    // ── structure, built in `validate` ───────────────────────────────────────
    n_channels: usize,
    wpc: usize,
    segs: Vec<OpticalSegment>,
    lines: Vec<NativeTLine>,
    /// External optical bundle wires: `in` block then `out` block.
    optical: Vec<NodeId>,
    rf_in: NodeId,
    rf_out: NodeId,
    /// Electrode nodes `e0 … eN`; the ends are the terminals, the interior is
    /// allocated.
    elec: Vec<NodeId>,
    /// The per-slice drive voltage at the current iterate.
    v_slice: Vec<f64>,
    /// Whether the group delay is engaged (the run-level option).
    delay_option: bool,
    min_terminals: Option<usize>,
}

impl Default for NativeTwPhaseShifter {
    fn default() -> Self {
        Self::new()
    }
}

impl NativeTwPhaseShifter {
    pub fn new() -> Self {
        let mut template = OpticalSegment::new(3e-3, 2.445, 4.19, 0.0);
        template.pin_at_ref = true;
        NativeTwPhaseShifter {
            length_m: 3e-3,
            n_m: 4.2,
            z0: 35.0,
            f_max: 50e9,
            slices_per_wave: 10.0,
            n_slices_override: None,
            // 1.2 V·cm, a typical depletion shifter.
            dn_dv: 0.0,
            v_pi_l: 0.012,
            template,
            n_channels: 0,
            wpc: 3,
            segs: Vec::new(),
            lines: Vec::new(),
            optical: Vec::new(),
            rf_in: None,
            rf_out: None,
            elec: Vec::new(),
            v_slice: Vec::new(),
            delay_option: false,
            min_terminals: None,
        }
    }

    /// Slices needed to keep the RF voltage roughly constant across one.
    ///
    /// At least two: one slice is a lumped device with extra steps, and the
    /// whole point of this element is that it is not one.
    fn slice_count(&self) -> usize {
        if let Some(n) = self.n_slices_override {
            return n.max(1);
        }
        let waves = self.slices_per_wave * self.length_m * self.n_m * self.f_max / C0;
        (waves.ceil() as usize).clamp(2, 512)
    }

    /// `Δn_eff` per volt at the reference wavelength, from `v_pi_l`.
    ///
    /// `V_pi·L = v_pi_l` and a π phase over `L` is `Δn_eff = λ/(2L)`, so
    /// `dn/dV = λ/(2·v_pi_l)` — independent of the length, which is what makes
    /// it a per-slice constant.
    fn refresh_dn_dv(&mut self) {
        self.dn_dv = self.template.wl_ref_m / (2.0 * self.v_pi_l);
    }

    /// Optical wires per bundle end.
    fn wires_per_end(&self) -> usize {
        self.wpc * self.n_channels
    }

    /// Interior optical rows: `re` and `im` for every channel at every internal
    /// boundary. The λ wire gets none — a wavelength is resolved before the
    /// solve, so an interior λ row would read zero while resolution says
    /// 1.55 µm.
    fn interior_optical_rows(&self) -> usize {
        self.segs.len().saturating_sub(1) * 2 * self.n_channels
    }
}

impl Device for NativeTwPhaseShifter {
    fn num_terminals(&self) -> usize {
        if let Some(min) = self.min_terminals {
            return min;
        }
        2 * self.wires_per_end() + 2
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.template.setup_model(ctx);
        self.wpc = ctx.wires_per_channel();
        self.delay_option = ctx.waveguide_delay;
        self.refresh_dn_dv();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        // `2·wpc·N + 2`: the optical bundle at both ends, plus the electrode's.
        if terminals.len() < 2 * wpc + 2 || !(terminals.len() - 2).is_multiple_of(2 * wpc) {
            self.min_terminals = Some(2 * wpc + 2);
            return;
        }
        self.n_channels = (terminals.len() - 2) / (2 * wpc);
        self.optical = terminals[..terminals.len() - 2].to_vec();
        self.rf_in = terminals[terminals.len() - 2];
        self.rf_out = terminals[terminals.len() - 1];
    }

    /// Build the ladder.
    ///
    /// Not in `setup_instance`, because the slice count depends on `f_max`,
    /// `n_m` and the length, and instance parameters are applied *after*
    /// `setup_instance`. `validate` is the first point where every parameter is
    /// in and the row count has not yet been asked for.
    fn validate(&mut self) -> Result<(), String> {
        if self.min_terminals.is_some() {
            return Err(format!(
                "expects {} terminals: the optical bundle in, the bundle out, then \
                 rf_in and rf_out",
                2 * self.wpc + 2
            ));
        }
        if self.wpc == 5 {
            return Err(
                "bidirectional mode is not supported: the backward wave needs \
                        its own ladder, and leaving it undriven would be worse"
                    .into(),
            );
        }
        if !(self.length_m > 0.0 && self.n_m > 0.0 && self.z0 > 0.0 && self.f_max > 0.0) {
            return Err(format!(
                "needs positive l_um, n_m, z0 and f_max (got {:.3e} m, {}, {} Ω, {:.3e} Hz)",
                self.length_m, self.n_m, self.z0, self.f_max
            ));
        }
        self.refresh_dn_dv();
        let n = self.slice_count();
        let dz = self.length_m / n as f64;
        let td = self.n_m * dz / C0;
        self.segs = (0..n)
            .map(|_| {
                let mut s = self.template.clone();
                s.length_m = dz;
                s.refresh_tau();
                s
            })
            .collect();
        self.lines = (0..n).map(|_| NativeTLine::new(self.z0, td)).collect();
        self.v_slice = vec![0.0; n];
        Ok(())
    }

    /// Interior electrode nodes, interior optical wires, and every child's own
    /// rows — all out of the one pool `push_device` allocates.
    fn num_extra_nodes(&self) -> usize {
        let n = self.segs.len();
        if n == 0 {
            return 0;
        }
        let interior_elec = n - 1;
        let seg_branches = n * (self.wpc - 1) * self.n_channels;
        let line_branches = n * 2;
        interior_elec + self.interior_optical_rows() + seg_branches + line_branches
    }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        let n = self.segs.len();
        if n == 0 {
            return;
        }
        let mut next = first_idx;
        let mut take = |count: usize| {
            let start = next;
            next += count;
            start
        };

        // Electrode: the ends are terminals, the interior is ours.
        self.elec = Vec::with_capacity(n + 1);
        self.elec.push(self.rf_in);
        let elec_first = take(n - 1);
        for i in 0..n - 1 {
            self.elec.push(Some(elec_first + i));
        }
        self.elec.push(self.rf_out);

        // Optical: build each slice's wire list, using the external bundle at
        // the two ends and fresh rows in between. λ wires are `None` inside —
        // see `interior_optical_rows`.
        let per_end = self.wires_per_end();
        let interior_first = take(self.interior_optical_rows());
        let interior_wire = |boundary: usize, k: usize, quad: usize| {
            Some(interior_first + (boundary * self.n_channels + k) * 2 + quad)
        };
        for i in 0..n {
            let mut wires = Vec::with_capacity(2 * per_end);
            for k in 0..self.n_channels {
                for q in 0..self.wpc {
                    wires.push(if i == 0 {
                        self.optical[self.wpc * k + q]
                    } else if q + 1 == self.wpc {
                        None // interior λ wire: resolved, never solved
                    } else {
                        interior_wire(i - 1, k, q)
                    });
                }
            }
            for k in 0..self.n_channels {
                for q in 0..self.wpc {
                    wires.push(if i + 1 == n {
                        self.optical[per_end + self.wpc * k + q]
                    } else if q + 1 == self.wpc {
                        None
                    } else {
                        interior_wire(i, k, q)
                    });
                }
            }
            let ctx = SimContext::default();
            self.segs[i].setup_instance(&wires, &ctx);
            self.segs[i].set_control_nodes(&[self.elec[i]]);
        }

        // Children's own rows.
        for seg in &mut self.segs {
            let first = take((self.wpc - 1) * self.n_channels);
            seg.bind_branches(first);
        }
        for (i, line) in self.lines.iter_mut().enumerate() {
            let ctx = SimContext::default();
            line.setup_instance(&[self.elec[i], None, self.elec[i + 1], None], &ctx);
            let first = take(2);
            line.bind_extra_nodes(first);
        }
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "l_um" => {
                self.length_m = value * 1e-6;
                true
            }
            "l_m" | "length" => {
                self.length_m = value;
                true
            }
            "n_m" | "n_rf" => {
                self.n_m = value;
                true
            }
            "z0" | "z_0" => {
                self.z0 = value;
                true
            }
            "f_max" | "fmax" => {
                self.f_max = value;
                true
            }
            "slices_per_wave" => {
                self.slices_per_wave = value.max(1.0);
                true
            }
            "n_slices" => {
                self.n_slices_override = Some(value.round().max(1.0) as usize);
                true
            }
            "v_pi_l" => {
                self.v_pi_l = value;
                self.refresh_dn_dv();
                true
            }
            other => {
                // Everything else is segment geometry, applied to the template
                // the slices are cut from. The per-slice length is ours.
                if other == "l_um" || other == "l_m" {
                    return false;
                }
                let took = self.template.set_param(other, value);
                if took {
                    self.refresh_dn_dv();
                }
                took
            }
        }
    }

    fn lambda_routing(&self) -> Vec<(usize, usize)> {
        let per_end = self.wires_per_end();
        let lam = self.wpc - 1;
        (0..self.n_channels)
            .map(|k| (self.wpc * k + lam, per_end + self.wpc * k + lam))
            .collect()
    }

    fn set_resolved_lambda(&mut self, per_terminal: &[f64]) {
        // Every slice sees the same wavelength: they are one waveguide cut into
        // pieces, and the interior λ wires deliberately do not exist.
        for seg in &mut self.segs {
            seg.set_resolved_lambda(per_terminal);
        }
    }

    fn eval(&mut self, x: &[f64], flags: EvalFlags, ctx: &SimContext) {
        let delay_active = flags.transient && self.delay_option;
        for i in 0..self.segs.len() {
            let v = self.elec[i].map_or(0.0, |j| x[j]);
            self.v_slice[i] = v;
            let dn = PerChannel::Uniform(self.dn_dv * v);
            let engaged = delay_active && self.segs[i].tau_g_s() > 0.0;
            self.segs[i].refresh_with_sens(
                x,
                dn,
                0.0,
                PerChannel::zero(),
                engaged,
                ctx,
                // One control node per slice, with the analytic derivative.
                &[(PerChannel::Uniform(self.dn_dv), 0.0, PerChannel::zero())],
            );
        }
        for line in &mut self.lines {
            line.eval(x, EvalFlags { ..flags }, ctx);
        }
    }

    fn load_residual(&self, b: &mut [f64]) {
        for seg in &self.segs {
            seg.stamp_residual(b);
        }
        for line in &self.lines {
            line.load_residual(b);
        }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        for seg in &self.segs {
            seg.stamp(mat);
        }
        for line in &self.lines {
            line.load_jacobian(mat);
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], alpha: f64) {
        for seg in &self.segs {
            seg.stamp_residual(b);
        }
        for line in &self.lines {
            line.load_residual_tran(b, alpha);
        }
    }

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, alpha: f64) {
        for seg in &self.segs {
            seg.stamp(mat);
        }
        for line in &self.lines {
            line.load_jacobian_tran(mat, alpha);
        }
    }

    fn ac_stamps(&self, omega: f64) -> Vec<AcStamp> {
        let mut out = Vec::new();
        for seg in &self.segs {
            out.extend(seg.ac_stamps(omega));
        }
        for line in &self.lines {
            out.extend(line.ac_stamps(omega));
        }
        out
    }

    /// The tightest step any part of the ladder needs.
    ///
    /// Both delays are per *slice*, so this shrinks as the device is cut
    /// finer — which is the price of resolving a travelling wave, and the
    /// reason `f_max` is a parameter rather than a guess.
    fn requested_max_timestep(&self) -> Option<f64> {
        let optical = self
            .segs
            .first()
            .filter(|_| self.delay_option)
            .map(|s| s.tau_g_s() / 2.0)
            .filter(|t| *t > 0.0);
        let electrical = self.lines.first().and_then(|l| l.requested_max_timestep());
        match (optical, electrical) {
            (Some(a), Some(b)) => Some(a.min(b)),
            (a, b) => a.or(b),
        }
    }

    fn commit_timestep(&mut self, x: &[f64]) {
        for seg in &mut self.segs {
            seg.commit(x);
        }
        for line in &mut self.lines {
            line.commit_timestep(x);
        }
    }

    fn frozen_jacobian_columns(&self) -> Vec<usize> {
        self.segs
            .iter()
            .flat_map(|s| s.unstamped_control_columns())
            .collect()
    }
}
