//! Parameter sensitivity of an **`.ac` sweep**, by the adjoint method.
//!
//! [`crate::adjoint`] differentiates one operating point and
//! [`crate::adjoint_tran`] a whole waveform.  This differentiates a frequency
//! response: given a scalar built from the sweep — a mismatch against a target
//! passband, the power at one wavelength, an integral over a band — it returns
//! `dL/dp` for every design parameter at a cost that does not grow with the
//! number of frequencies being traded off.
//!
//! It is the shape a photonic design problem actually takes.  "Put the
//! resonance here", "flatten this passband", "hit this free spectral range" are
//! all least-squares fits of `|H(f)|²` against a target across a sweep, and the
//! parameters are ring radii, coupling coefficients and phase-shifter biases.
//!
//! # The system
//!
//! Each frequency is an independent linear constraint.  `.ac` solves the real
//! block form of `(G + jB)·V = b` with `B = ωC − L/ω`:
//!
//! ```text
//!     A_k·V_k = b,     A_k = [ G  −B_k ]      V_k = [ Re V ]
//!                            [ B_k   G ]            [ Im V ]
//! ```
//!
//! so for `L = Σ_k w_k·g(V_k)`,
//!
//! ```text
//!     A_kᵀ·λ_k = w_k·∇g(V_k)          dL/dp = −Σ_k λ_kᵀ·∂(A_k·V_k − b)/∂p
//! ```
//!
//! One transposed solve per frequency, and the parameter derivative is taken
//! the way both other adjoints take it: re-assemble at a perturbed parameter
//! with `V_k` **frozen** and difference the residual.  Every path a parameter
//! takes into `G`, `C`, `L` or the excitation is then included without being
//! enumerated — which matters here more than in DC, because an optical length
//! reaches the answer through `G` *and* through the operating point the
//! small-signal matrices were linearised about.
//!
//! # Why the transpose and not the conjugate transpose
//!
//! The real block form of `M = G + jB` has `[G −B; B G]`, whose plain transpose
//! is `[Gᵀ Bᵀ; −Bᵀ Gᵀ]` — the real block form of `Mᴴ`, not `Mᵀ`.  For a
//! magnitude objective that distinction does not survive `|·|²`, and
//! `crate::noise` relies on exactly that.  Here it does matter, because a
//! gradient carries a sign: the seed below is written in the real block basis
//! and solved against the real block transpose, so the pairing is
//! `⟨λ, r⟩` in `R^2n` throughout and no complex conjugation is ever implied.

use fairchild_parser::Netlist;

use crate::ac::{assemble_ac, AcSystem};
use crate::adjoint::{default_step, resolve, ParamRef};
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::mna::CircuitTopology;
use crate::newton::device_element_names;
use crate::options::SimOptions;

/// A scalar read off one frequency point.
#[derive(Clone, Debug)]
pub enum AcOutput {
    /// `|V(node)|²` — power, and the objective a filter fit wants.  Smooth
    /// everywhere including at a null, which `|V|` is not.
    MagSquared { node: String },
    /// `Re V(node)`.
    Real { node: String },
    /// `Im V(node)`.
    Imag { node: String },
}

impl AcOutput {
    fn node(&self) -> &str {
        match self {
            AcOutput::MagSquared { node } | AcOutput::Real { node } | AcOutput::Imag { node } => {
                node
            }
        }
    }

    /// The scalar and its gradient `∂g/∂V_k` in the real block basis.
    fn seed(&self, topo: &CircuitTopology, v: &[f64]) -> Result<(f64, Vec<f64>), SimError> {
        let size = topo.size;
        let idx = *topo
            .node_index
            .get(self.node())
            .ok_or_else(|| SimError::ParameterError(format!("unknown node '{}'", self.node())))?;
        let (re, im) = (v[idx], v[size + idx]);
        let mut g = vec![0.0; 2 * size];
        let value = match self {
            AcOutput::MagSquared { .. } => {
                g[idx] = 2.0 * re;
                g[size + idx] = 2.0 * im;
                re * re + im * im
            }
            AcOutput::Real { .. } => {
                g[idx] = 1.0;
                re
            }
            AcOutput::Imag { .. } => {
                g[size + idx] = 1.0;
                im
            }
        };
        Ok((value, g))
    }
}

/// A captured `.ac` sweep, ready to differentiate.
pub struct AcAdjoint {
    netlist: Netlist,
    opts: SimOptions,
    ac_source: Option<String>,
    sys: AcSystem,
    freqs: Vec<f64>,
    /// The solved `[Re V; Im V]` at each frequency.
    v: Vec<Vec<f64>>,
}

/// `dL/dp` for one objective over a sweep.
#[derive(Clone, Debug)]
pub struct AcSensitivities {
    pub grad: Vec<f64>,
    /// Whether the perturbation changed the equations at all.  A `false` is a
    /// misspelled parameter or a device that ignored `set_real_param` — both
    /// reportable, neither a gradient.
    pub reached: Vec<bool>,
    /// Two-step-size disagreement, as the honest error bar on each entry.
    pub fd_error: Vec<f64>,
}

impl AcAdjoint {
    /// Solve the sweep and keep everything the backward pass needs.
    pub fn run(
        netlist: &Netlist,
        registry: &DeviceRegistry,
        opts: &SimOptions,
        freqs: &[f64],
        ac_source: Option<&str>,
    ) -> Result<Self, SimError> {
        let sys = assemble_ac(netlist, ac_source, registry, opts)?;
        let solver = opts.linear_solver(2 * sys.topo.size);
        let mut v = Vec::with_capacity(freqs.len());
        for &f in freqs {
            let (a2, rhs) = sys.at(f);
            v.push(solver.solve(&CircuitTopology::sparse_from_dense(&a2), &rhs)?);
        }
        Ok(Self {
            netlist: netlist.clone(),
            opts: opts.clone(),
            ac_source: ac_source.map(str::to_string),
            sys,
            freqs: freqs.to_vec(),
            v,
        })
    }

    pub fn freqs(&self) -> &[f64] {
        &self.freqs
    }

    /// The objective's value per frequency, without weighting.
    pub fn response(&self, out: &AcOutput) -> Result<Vec<f64>, SimError> {
        self.v
            .iter()
            .map(|v| out.seed(&self.sys.topo, v).map(|(value, _)| value))
            .collect()
    }

    /// `L = Σ_k w_k·g(V_k)` and the per-frequency seeds `∂L/∂V_k`.
    ///
    /// Weights are the caller's: `[0,…,1,…,0]` picks one frequency, `df`
    /// integrates a band, and `2·(g_k − target_k)·df` is a least-squares fit
    /// against a target response — which is the one that makes this useful.
    pub fn weighted(
        &self,
        out: &AcOutput,
        weights: &[f64],
    ) -> Result<(f64, Vec<Vec<f64>>), SimError> {
        if weights.len() != self.freqs.len() {
            return Err(SimError::ParameterError(format!(
                "weights has {} entries but the sweep has {}",
                weights.len(),
                self.freqs.len()
            )));
        }
        let mut total = 0.0;
        let mut seeds = Vec::with_capacity(self.freqs.len());
        for (v, &w) in self.v.iter().zip(weights) {
            let (value, mut g) = out.seed(&self.sys.topo, v)?;
            total += w * value;
            for e in g.iter_mut() {
                *e *= w;
            }
            seeds.push(g);
        }
        Ok((total, seeds))
    }

    /// `dL/dp` for each parameter, from one transposed solve per frequency.
    pub fn gradient(
        &self,
        registry: &DeviceRegistry,
        seeds: &[Vec<f64>],
        params: &[ParamRef],
    ) -> Result<AcSensitivities, SimError> {
        if seeds.len() != self.freqs.len() {
            return Err(SimError::ParameterError(format!(
                "seeds has {} timepoints but the sweep has {}",
                seeds.len(),
                self.freqs.len()
            )));
        }
        let size = self.sys.topo.size;
        let solver = self.opts.linear_solver(2 * size);

        // Backward pass: A_kᵀ λ_k = ∂L/∂V_k, one per frequency.
        let mut lambdas = Vec::with_capacity(self.freqs.len());
        for (k, &f) in self.freqs.iter().enumerate() {
            let (a2, _) = self.sys.at(f);
            let mut at = vec![vec![0.0f64; 2 * size]; 2 * size];
            for (i, row) in a2.iter().enumerate() {
                for (j, val) in row.iter().enumerate() {
                    at[j][i] = *val;
                }
            }
            lambdas.push(solver.solve(&CircuitTopology::sparse_from_dense(&at), &seeds[k])?);
        }

        let names = device_element_names(&self.netlist);
        let mut grad = vec![0.0; params.len()];
        let mut reached = vec![false; params.len()];
        let mut fd_error = vec![0.0; params.len()];

        for (pi, p) in params.iter().enumerate() {
            // `resolve` is used only for the nominal value. Its `Handle` says
            // whether the DC path would mutate a live device in place, which is
            // irrelevant here: this re-assembles the whole system from the
            // netlist per sample, so every parameter — device instance params
            // included — is set on the netlist and picked up by that rebuild.
            let Some((_, nominal)) = resolve(p, &self.netlist, &names) else {
                continue;
            };
            let h = p.step.unwrap_or_else(|| default_step(nominal));

            // Two step sizes, then Richardson — the same construction the DC
            // path uses, and for the same reason: a parameter whose natural
            // scale is not the residual's (an optical length moves the
            // propagation phase by ~17 rad/µm) needs the O(h²) term cancelled
            // or the difference is dominated by it.
            let mut d = [0.0_f64; 2];
            let mut moved = false;
            for (i, step) in [h, 0.5 * h].into_iter().enumerate() {
                let plus = self.replay(registry, p, nominal + step, &lambdas)?;
                let minus = self.replay(registry, p, nominal - step, &lambdas)?;
                d[i] = (plus.0 - minus.0) / (2.0 * step);
                moved |= plus.1 || minus.1;
            }
            let value = (4.0 * d[1] - d[0]) / 3.0;
            fd_error[pi] = if value != 0.0 {
                (d[1] - d[0]).abs() / value.abs()
            } else {
                0.0
            };
            if !moved {
                continue;
            }
            reached[pi] = true;
            grad[pi] = -value;
        }

        Ok(AcSensitivities {
            grad,
            reached,
            fd_error,
        })
    }

    /// `Σ_k λ_kᵀ·(A_k·V_k − b)` with one parameter held at `value` and every
    /// `V_k` frozen, plus whether anything moved.
    ///
    /// Re-assembling the whole system per sample is the expensive part — it
    /// redoes the operating point — but it is also what makes the gradient a
    /// *total* derivative: a parameter that moves the DC bias moves the
    /// small-signal `G` the sweep was linearised about, and enumerating that
    /// path by hand is exactly the kind of thing that silently comes out zero.
    fn replay(
        &self,
        registry: &DeviceRegistry,
        p: &ParamRef,
        value: f64,
        lambdas: &[Vec<f64>],
    ) -> Result<(f64, bool), SimError> {
        let mut work = self.netlist.clone();
        if !crate::netlist_edit::set_element_param(&mut work, &p.element, &p.param, value) {
            // Nothing to differentiate: report it as unreached rather than
            // returning a confident zero.
            return Ok((0.0, false));
        }
        let sys = assemble_ac(&work, self.ac_source.as_deref(), registry, &self.opts)?;
        let size = sys.topo.size;
        if size != self.sys.topo.size {
            return Err(SimError::ParameterError(
                "the perturbed netlist has a different row count; a parameter that \
                 changes the topology cannot be differentiated this way"
                    .into(),
            ));
        }

        let mut acc = 0.0;
        let mut moved = false;
        for (k, &f) in self.freqs.iter().enumerate() {
            let (a2, rhs) = sys.at(f);
            let v = &self.v[k];
            for (i, row) in a2.iter().enumerate() {
                let r: f64 = row.iter().zip(v.iter()).map(|(a, x)| a * x).sum::<f64>() - rhs[i];
                if r != 0.0 {
                    moved = true;
                }
                acc += lambdas[k][i] * r;
            }
        }
        Ok((acc, moved))
    }
}
