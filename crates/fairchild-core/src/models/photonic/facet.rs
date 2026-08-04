use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

use super::stamp_potential_eq;

// ────────────────────────────────────────────────────────────────────────
// Single-port facet: terminator ↔ partial reflector ↔ mirror
// ────────────────────────────────────────────────────────────────────────

/// One optical port whose forward field is split three ways: reflected back
/// into the same port, transmitted out of the model, absorbed.
///
/// ```text
///   ── A_fw ──►│                     R + T + L = 1
///   ◄─ A_bw ───│  ─── T ──►  (leaves the simulation)
/// ```
///
/// It is a chip facet, which is why it is named for one: an AR-coated facet is
/// a terminator (`R = 0`, the default), a cleaved facet is `R ≈ 0.3`, an HR
/// coating is a mirror (`R = 1`), and a fibre-coupled edge is mostly `T`.  One
/// device covers all of them because from inside the circuit they differ only
/// in how much light comes back.
///
/// **Only `reflectance` changes the answer.**  `transmittance` and `loss` are
/// bookkeeping: light that leaves through either one is gone either way, and
/// there is no second port for it to arrive at.  They exist so the budget is
/// written down and checked — `reflectance=0.9 transmittance=0.5` is a typo the
/// device should catch, not average away.
///
/// Set any one, any two, or all three; the unset ones absorb the remainder
/// (`loss` first).  Setting all three requires them to sum to 1.
///
/// | Parameter | Default | Meaning |
/// |---|---|---|
/// | `reflectance` / `r` | 0 | Power fraction returned into the port. |
/// | `transmittance` / `t` | 0 | Power fraction leaving the model. |
/// | `loss` | remainder | Power fraction absorbed. |
/// | `phase_deg` | 0 | Phase added on reflection (180° for a metal mirror). |
///
/// **Needs `.options enable_bidirectional=1`** for any non-zero reflectance —
/// a unidirectional bundle has no backward wire to drive.  Under unidirectional
/// propagation the device is a pure terminator and a non-zero `reflectance`
/// is a hard error rather than a silent no-op.
///
/// Bundle-aware: `wpc·N` terminals for `N` WDM channels, with one budget shared
/// across channels.  Wavelength-dependent reflectance (a DBR, a coated facet
/// near its design wavelength) would want `spectrum.rs`; this device is flat.
///
/// This is an **end cap**: light arrives on the port's forward wires and leaves
/// on its backward ones.  A Fabry-Pérot cavity needs a partial mirror that
/// couples an outside port to an inside one in both directions — two ports, not
/// one — so it is a different device, closer in shape to `fc_dcoupler`.
pub struct NativeFacet {
    reflectance: Option<f64>,
    transmittance: Option<f64>,
    loss: Option<f64>,
    phase_deg: f64,
    /// `√R`, resolved on the first `eval` — the budget cannot be checked at
    /// `setup_instance` because the registry sets parameters after it.
    rho: f64,
    checked: bool,
    n_channels: usize,
    wpc: usize,
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,
}

impl Default for NativeFacet {
    fn default() -> Self {
        Self::new()
    }
}

impl NativeFacet {
    pub fn new() -> Self {
        Self {
            reflectance: None,
            transmittance: None,
            loss: None,
            phase_deg: 0.0,
            rho: 0.0,
            checked: false,
            n_channels: 0,
            wpc: 3,
            nodes: Vec::new(),
            branches: Vec::new(),
        }
    }

    /// Resolve `R`, `T`, `L` from whichever were given and check they are a
    /// power budget.  Panics on an over-unity or negative one — a facet that
    /// quietly normalised its own numbers would hide the typo it exists to
    /// catch.
    fn resolve(&mut self) {
        let given: Vec<(&str, f64)> = [
            ("reflectance", self.reflectance),
            ("transmittance", self.transmittance),
            ("loss", self.loss),
        ]
        .into_iter()
        .filter_map(|(n, v)| v.map(|v| (n, v)))
        .collect();
        for (name, v) in &given {
            assert!(
                (0.0..=1.0).contains(v),
                "fc_facet: {name}={v} is not a power fraction in [0, 1]"
            );
        }
        let sum: f64 = given.iter().map(|(_, v)| v).sum();
        match given.len() {
            3 => assert!(
                (sum - 1.0).abs() < 1e-9,
                "fc_facet: reflectance + transmittance + loss = {sum}, must be 1 \
                 (leave one out and it takes the remainder)"
            ),
            _ => assert!(
                sum <= 1.0 + 1e-9,
                "fc_facet: {} sum to {sum} > 1 — no power left for the rest",
                given
                    .iter()
                    .map(|(n, _)| *n)
                    .collect::<Vec<_>>()
                    .join(" + ")
            ),
        }
        let r = match (self.reflectance, self.transmittance, self.loss) {
            (Some(r), _, _) => r,
            // R is the remainder only when the other two pin it; otherwise a
            // facet that was told nothing about reflection does not reflect.
            (None, Some(t), Some(l)) => 1.0 - t - l,
            _ => 0.0,
        };
        self.rho = r.max(0.0).sqrt();
        assert!(
            self.rho == 0.0 || self.wpc == 5,
            "fc_facet: reflectance={r} needs a backward wire; \
             set `.options enable_bidirectional=1` (or use reflectance=0 as a terminator)"
        );
        self.checked = true;
    }
}

impl Device for NativeFacet {
    fn num_terminals(&self) -> usize {
        self.nodes.len()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        self.wpc = wpc;
        assert!(
            !terminals.is_empty() && terminals.len().is_multiple_of(wpc),
            "fc_facet: terminal count must be {wpc}·N for N ≥ 1 channels; got {}",
            terminals.len()
        );
        self.n_channels = terminals.len() / wpc;
        self.nodes = terminals.to_vec();
        // Two driven potentials per channel (re_bw, im_bw), and none at all
        // under unidirectional propagation — where the device is a terminator
        // and the port's only wires are already driven from upstream.
        self.branches = vec![None; if wpc == 5 { 2 * self.n_channels } else { 0 }];
    }

    fn num_extra_nodes(&self) -> usize {
        self.branches.len()
    }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for (i, b) in self.branches.iter_mut().enumerate() {
            *b = Some(first_idx + i);
        }
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "reflectance" | "r" | "r_power" => self.reflectance = Some(value),
            "transmittance" | "t" | "t_power" => self.transmittance = Some(value),
            "loss" | "l" | "loss_power" => self.loss = Some(value),
            "phase_deg" | "phi_deg" => self.phase_deg = value,
            _ => return false,
        }
        self.checked = false;
        true
    }

    fn eval(&mut self, _x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        if !self.checked {
            self.resolve();
        }
    }

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        if self.branches.is_empty() {
            return;
        }
        // A_bw = ρ·e^(−jφ)·A_fw, the same sign convention `OpticalSegment` uses
        // for propagation phase, so a facet and a length of waveguide compose
        // the way the algebra says they should.
        let phi = self.phase_deg.to_radians();
        let (c, s) = (self.rho * phi.cos(), self.rho * phi.sin());
        for k in 0..self.n_channels {
            let base = self.wpc * k;
            let re_fw = self.nodes[base];
            let im_fw = self.nodes[base + 1];
            let re_bw = self.nodes[base + 2];
            let im_bw = self.nodes[base + 3];
            stamp_potential_eq(
                mat,
                &self.branches,
                2 * k,
                re_bw,
                &[(re_fw, -c), (im_fw, -s)],
            );
            stamp_potential_eq(
                mat,
                &self.branches,
                2 * k + 1,
                im_bw,
                &[(re_fw, s), (im_fw, -c)],
            );
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }
}
