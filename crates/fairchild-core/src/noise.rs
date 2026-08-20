//! Small-signal noise analysis.
//!
//! For each swept frequency:
//!   1. Build the same complex linearized system A(jω) = G + jB used by AC.
//!   2. Solve the adjoint problem  A^T · λ = e_out  to obtain the transfer
//!      impedance from any internal current injection to the observation
//!      node V(out_pos) − V(out_neg).
//!   3. For every uncorrelated noise source k at nodes (p,n) with one-sided
//!      current PSD S_ik(f) [A²/Hz]:
//!      S_V_out_k(f) = |λ[p] − λ[n]|² · S_ik(f).
//!   4. Sum the contributions and divide by the squared signal-path gain to
//!      get the input-referred PSD.
//!
//! This commit ships the full plumbing (parser → CLI/Python) plus resistor
//! thermal noise (4kT/R).  Diode shot noise and MOSFET channel noise land in
//! the next commit once the Device trait grows a `noise_sources()` hook.

use crate::warn_user;
use indexmap::IndexMap;

use fairchild_parser::{Element, Netlist};

use crate::device::{Device, EvalFlags, NodeId, ReactiveKind, SimContext};
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::mna::{
    stamp_2port_by_id, stamp_netlist_scaled, stamp_passive_2port, CircuitTopology, RowFloor,
};
use crate::newton::build_devices;
use crate::options::SimOptions;
use crate::solver::LinearSolver;

/// Boltzmann's constant in J/K.
const KB: f64 = 1.380649e-23;

// ───────────────────────────────────────────────────────────────────────────
// The one enumeration of "what is a noise source in this circuit"
// ───────────────────────────────────────────────────────────────────────────

/// Every noise generator in a circuit, resolved to MNA rows.
///
/// `.noise` and transient noise are two consumers of the same physics, and a
/// source that exists in only one of them is a silent inconsistency — the
/// frequency-domain answer and the time-domain answer would simply disagree,
/// with nothing to notice. So the list lives here once, exactly as
/// `crate::reactive` owns "what is reactive". `transient_noise_agrees_with_the
/// _noise_analysis` is the regression that holds the two together.
///
/// Resistor taps are resolved once at construction because `4kT/R` and the node
/// indices are static. Device sources are asked for on every visit: shot noise
/// follows the bias, and in a transient the bias moves.
pub struct NoiseSources {
    /// `(pos, neg, 4kT/R)` for every linear resistor with `R > 0`.
    resistors: Vec<(NodeId, NodeId, f64)>,
}

impl NoiseSources {
    pub fn build(netlist: &Netlist, topo: &CircuitTopology, temp_k: f64) -> Self {
        let four_kt = 4.0 * KB * temp_k;
        let resistors = netlist
            .elements
            .iter()
            .filter_map(|el| match el {
                Element::Resistor {
                    pos,
                    neg,
                    resistance,
                    ..
                } if *resistance > 0.0 => Some((
                    topo.node_index.get(pos).copied(),
                    topo.node_index.get(neg).copied(),
                    four_kt / resistance,
                )),
                _ => None,
            })
            .collect();
        Self { resistors }
    }

    /// Visit `(taps, psd)` for every generator at the devices' current bias and
    /// at frequency `freq` (Hz).
    ///
    /// `psd` is one-sided, in A²/Hz for an injection into a node row and in the
    /// square of the enforced potential's unit for a branch row. Taps belonging
    /// to one generator are driven by one random process and must be combined
    /// coherently by the caller.
    ///
    /// Every native generator here is flat, so `freq` changes nothing for them.
    /// It exists because an OSDI model may call `flicker_noise()`, and a hook
    /// that cannot express frequency dependence would silently flatten it.
    pub fn for_each(
        &self,
        devices: &[Box<dyn Device>],
        ctx: &SimContext,
        freq: f64,
        mut f: impl FnMut(&[(NodeId, NodeId, f64)], f64),
    ) {
        for &(p, n, s_i) in &self.resistors {
            f(&[(p, n, 1.0)], s_i);
        }
        for dev in devices {
            for (p, n, s_i) in dev.noise_sources(ctx, freq) {
                f(&[(p, n, 1.0)], s_i);
            }
            for src in dev.correlated_noise_sources(ctx) {
                f(&src.taps, src.psd);
            }
        }
    }
}

// ───────────────────────────────────────────────────────────────────────────
// Transient noise
// ───────────────────────────────────────────────────────────────────────────

/// splitmix64 — the whole random-number dependency.
///
/// A noise injector needs a reproducible stream of gaussians, not a
/// cryptographic one, and splitmix64 clears BigCrush in five lines. Pulling in
/// `rand` for this would add a dependency tree to the simulator so a test can
/// draw a bell curve.
struct Rng {
    state: u64,
    /// Box-Muller produces two normals per pair of uniforms; keep the spare.
    spare: Option<f64>,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Self {
            state: seed,
            spare: None,
        }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform on [0, 1) with 53 bits of mantissa.
    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Standard normal, Box-Muller.
    fn normal(&mut self) -> f64 {
        if let Some(v) = self.spare.take() {
            return v;
        }
        // `ln(0)` is −inf and would poison the whole step; the uniform is on
        // [0, 1) so zero is reachable, once in 2^53 draws.
        let u1 = self.uniform().max(f64::MIN_POSITIVE);
        let u2 = self.uniform();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = std::f64::consts::TAU * u2;
        self.spare = Some(r * theta.sin());
        r * theta.cos()
    }
}

/// Time-domain noise: the same generators `.noise` reports as PSDs, realised as
/// random currents injected at every timestep.
///
/// # The amplitude
///
/// A zero-order-held random sequence of variance `σ²` at interval `h` has a
/// one-sided PSD of `2σ²h` below its Nyquist frequency, so a generator with PSD
/// `S` is realised by drawing
///
/// ```text
/// i_n = √(S / 2h) · N(0, 1)
/// ```
///
/// and holding it for the step. The consistency check: integrating `S` over the
/// resolved band `[0, 1/2h]` gives `S/2h`, which is exactly `σ²`.
///
/// # Bandwidth is set by the timestep — and that is usually fine
///
/// The injected noise is white up to `1/2h` and absent above it. Real thermal
/// and shot noise are white far past any timestep worth using, so this is a
/// truncation. It does not bias observable quantities, because any transient
/// that resolves the circuit at all has its Nyquist frequency well above the
/// circuit's own bandwidth, and the circuit filters the difference away: an RC
/// low-pass settles at `kT/C` for any `h ≪ RC`, independent of `h`. What *is*
/// step-dependent is a measurement with no bandwidth limit of its own — the
/// voltage across a bare resistor divider has variance `S_V/2h` and always
/// will, in any simulator, because unbounded-bandwidth white noise has infinite
/// power.
///
/// ponytail: no `noisefmax`. Decoupling the noise bandwidth from the timestep
/// means holding each sample for several steps, which is ~8 lines — worth it
/// only if someone needs an unfiltered node to be step-independent.
///
/// # What is missing
///
/// Flicker (1/f) and RTS noise, which `.noise` does not model either — this
/// injects exactly the generators that analysis knows about, no more. The
/// sequence is white and gaussian per generator, and independent generators are
/// independent here as they are there.
pub struct TransientNoise {
    sources: NoiseSources,
    rng: Rng,
    scale: f64,
    /// RHS contribution for the current step. Frozen across Newton iterations —
    /// see [`TransientNoise::draw`].
    rhs: Vec<f64>,
    /// The flatness probe in `draw` runs once, not once per timestep.
    flatness_checked: bool,
}

impl TransientNoise {
    /// `None` unless `.options trannoise=1`, so the cost is zero when off and
    /// the caller has one branch instead of a flag test per step.
    pub fn new(netlist: &Netlist, topo: &CircuitTopology, opts: &SimOptions) -> Option<Self> {
        if !opts.trannoise {
            return None;
        }
        Some(Self {
            sources: NoiseSources::build(netlist, topo, opts.temp_k),
            rng: Rng::new(opts.noiseseed),
            scale: opts.noisescale,
            rhs: vec![0.0; topo.size],
            flatness_checked: false,
        })
    }

    /// Draw one sample per generator for a step of length `h`, at the bias the
    /// devices are currently evaluated at.
    ///
    /// **Once per timestep, never per Newton iteration.** Redrawing inside the
    /// loop changes the equation being solved between iterations, so Newton is
    /// chasing a target that moves as fast as it does and simply never
    /// converges. Freezing the sample also makes it what it physically is: one
    /// noise realisation held across the step.
    pub fn draw(&mut self, devices: &[Box<dyn Device>], ctx: &SimContext, h: f64) {
        // A held i.i.d. sample sequence is white by construction, so there is
        // exactly one density it can realise. Probe mid-band of what this step
        // resolves — [0, 1/2h], so 1/4h — and check once that the generators
        // really are flat across it. Every native source is; an OSDI model
        // calling `flicker_noise()` is not, and silently flattening it is the
        // shape of bug this file exists to avoid.
        let f_mid = 1.0 / (4.0 * h);
        if !self.flatness_checked {
            self.flatness_checked = true;
            self.warn_if_not_flat(devices, ctx, h);
        }
        self.rhs.fill(0.0);
        let rng = &mut self.rng;
        let scale = self.scale;
        let rhs = &mut self.rhs;
        self.sources.for_each(devices, ctx, f_mid, |taps, psd| {
            if psd <= 0.0 {
                return;
            }
            let amp = scale * (psd / (2.0 * h)).sqrt() * rng.normal();
            for &(p, n, w) in taps {
                // Same sign convention as `stamp_current_source_at`, mirrored:
                // a positive sample drives current INTO `pos`, which is the
                // direction `.noise`'s `λ[p] − λ[n]` transfer assumes.
                if let Some(i) = p {
                    rhs[i] += amp * w;
                }
                if let Some(i) = n {
                    rhs[i] -= amp * w;
                }
            }
        });
    }

    /// Warn once if any generator's density varies across the resolved band.
    ///
    /// A held sample sequence is white, so a sloped density cannot be realised
    /// by this injector at all — the honest options are to say so or to refuse,
    /// and a warning naming the ratio lets the user decide whether the band
    /// they care about is flat enough. `.noise` handles such a source exactly;
    /// it is only the time-domain injection that cannot.
    fn warn_if_not_flat(&self, devices: &[Box<dyn Device>], ctx: &SimContext, h: f64) {
        let (f_lo, f_hi) = (1.0 / (20.0 * h), 1.0 / (2.2 * h));
        let mut lo = Vec::new();
        self.sources
            .for_each(devices, ctx, f_lo, |_, psd| lo.push(psd));
        let mut worst = 1.0_f64;
        let mut i = 0;
        self.sources.for_each(devices, ctx, f_hi, |_, psd| {
            if let Some(&a) = lo.get(i) {
                if a > 0.0 && psd > 0.0 {
                    worst = worst.max((psd / a).max(a / psd));
                }
            }
            i += 1;
        });
        if worst > 1.05 {
            warn_user!(
                "a noise generator varies by {worst:.2}x across the band this \
                 timestep resolves ({:.3e} to {:.3e} Hz), but transient noise injects a \
                 white sample sequence and can only realise one density — it uses the \
                 mid-band value. `.noise` handles the frequency dependence exactly; for \
                 the time domain, either shorten the step so the band sits where the \
                 density is flat, or treat the result as the mid-band approximation it is.",
                f_lo,
                f_hi
            );
        }
    }

    /// This step's RHS contribution, to be added after every netlist stamp.
    pub fn rhs(&self) -> &[f64] {
        &self.rhs
    }
}

/// Result of a `.noise` sweep.  All PSDs are one-sided (V²/Hz).
pub struct NoiseResult {
    pub freq: Vec<f64>,
    /// Total output-referred voltage noise PSD at each frequency (V²/Hz).
    pub onoise_psd: Vec<f64>,
    /// Input-referred PSD = onoise / |H(f)|² where H is the small-signal gain
    /// from `input_src` to (out_pos, out_neg).  NaN at frequencies where the
    /// transfer function is too small to invert reliably.
    pub inoise_psd: Vec<f64>,
}

/// One-sided PSD (V²/Hz) → amplitude density (V/√Hz), the form every rawfile
/// reader expects under `voltage-density`.
///
/// A tiny negative PSD is numerical dust and clamps to zero; a NaN stays NaN.
/// The distinction matters: `inoise_psd` is deliberately NaN where the transfer
/// function is too small to invert, and reporting that as `0` would claim a
/// noiseless input — a plausible-looking number in place of "not computable".
fn amplitude_density(psd: f64) -> f64 {
    if psd.is_nan() {
        f64::NAN
    } else {
        psd.max(0.0).sqrt()
    }
}

impl NoiseResult {
    /// Write a CSV with columns `freq_hz, onoise_v2hz, onoise_vrthz, inoise_v2hz, inoise_vrthz`.
    pub fn write_csv<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        writeln!(
            w,
            "freq_hz,onoise_v2hz,onoise_vrthz,inoise_v2hz,inoise_vrthz"
        )?;
        for i in 0..self.freq.len() {
            let on = self.onoise_psd[i];
            let in_ = self.inoise_psd[i];
            writeln!(
                w,
                "{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}",
                self.freq[i],
                on,
                amplitude_density(on),
                in_,
                amplitude_density(in_)
            )?;
        }
        Ok(())
    }

    /// Write the noise sweep as an ngspice-compatible Nutmeg ASCII rawfile.
    ///
    /// **Values are amplitude densities in V/√Hz, not the PSDs in V²/Hz that
    /// this struct stores.**  That is not a choice: the Nutmeg variable type
    /// `voltage-density` means V/√Hz, and `onoise_spectrum` / `inoise_spectrum`
    /// are the names every rawfile reader expects those units under.  Emitting
    /// the stored PSD under these names would be wrong by a square, and wrong
    /// in the way that looks plausible — a reader has no way to tell V/√Hz from
    /// V²/Hz by inspection, so nothing downstream would report a fault.  The
    /// CSV writer keeps both, since its column names say which is which.
    pub fn write_nutmeg<W: std::io::Write>(&self, mut w: W, title: &str) -> std::io::Result<()> {
        let n_pts = self.freq.len();
        writeln!(w, "Title: {title}")?;
        writeln!(w, "Plotname: Noise Spectral Density Curves")?;
        writeln!(w, "Flags: real")?;
        writeln!(w, "No. Variables: 3")?;
        writeln!(w, "No. Points: {n_pts}")?;
        writeln!(w, "Variables:")?;
        writeln!(w, "\t0\tfrequency\tfrequency")?;
        writeln!(w, "\t1\tonoise_spectrum\tvoltage-density")?;
        writeln!(w, "\t2\tinoise_spectrum\tvoltage-density")?;
        writeln!(w, "Values:")?;
        for i in 0..n_pts {
            // Point index on the first variable's line only, as in the other
            // analyses' writers.
            writeln!(w, " {i}\t{:.6e}", self.freq[i])?;
            writeln!(w, "\t{:.6e}", amplitude_density(self.onoise_psd[i]))?;
            writeln!(w, "\t{:.6e}", amplitude_density(self.inoise_psd[i]))?;
        }
        Ok(())
    }
}

/// Run a `.noise` analysis.
///
/// `freqs` — already-expanded sweep points (use the same `freq_decade` /
///   `freq_linear` / `freq_oct` helpers as `ac::ac_analysis`).
/// `out_pos`, `out_neg` — observation node pair (`"0"` for ground).
/// `input_src` — name of the voltage source that defines the input port; used
///   only for input-referred PSD via the signal-path gain.
pub fn noise_analysis(
    netlist: &Netlist,
    freqs: &[f64],
    out_pos: &str,
    out_neg: &str,
    input_src: &str,
    registry: &DeviceRegistry,
    opts: &SimOptions,
) -> Result<NoiseResult, SimError> {
    crate::connectivity::check_connectivity(netlist)?;
    let ctx = opts.sim_context();
    let mut topo = CircuitTopology::build(netlist);
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();

    // DC operating point.
    let mut devices = build_devices(netlist, &mut topo, &ctx, registry)?;
    let n_nodes = topo.n_nodes();
    let size = topo.size;
    let dc_solver = opts.linear_solver(size);
    let noise_solver = opts.linear_solver(2 * size);
    let x0 = run_dc_op(&topo, netlist, &mut devices, &ctx, opts, &*dc_solver)?;

    // Real (G) and imaginary-coefficient (C, L⁻¹) parts of Y(jω) = G + j(ωC − L⁻¹/ω).
    let mut g_mat = stamp_netlist_scaled(
        &topo,
        netlist,
        1.0,
        &empty,
        &empty,
        crate::mna::InductorDc::Reactive,
    )
    .a;
    for dev in devices.iter_mut() {
        // tran() flags so device small-signal cap caches populate (see ac.rs).
        dev.eval(&x0, EvalFlags::tran(), &ctx);
        let mut tmp = crate::mna::MnaMatrix::zeros(size);
        dev.load_jacobian(&mut tmp);
        // `tmp.a` is sparse now, so this walks only the cells the device
        // actually stamped instead of the whole row.
        for (g_row, t_row) in g_mat.iter_mut().zip(tmp.a.iter()) {
            for (j, t) in t_row.iter() {
                g_row[j] += t;
            }
        }
    }
    topo.stamp_gmin(&mut g_mat, opts.gmin, RowFloor::GminOnly);
    // ponytail: dense G/C/L and a dense 2n×2n adjoint system, same trade-off
    // and same upgrade path as `ac.rs`. Tracked as task #12.
    let mut c_mat = vec![vec![0.0f64; size]; size];
    let mut l_mat = vec![vec![0.0f64; size]; size];
    for el in &netlist.elements {
        match el {
            Element::Capacitor {
                pos,
                neg,
                capacitance,
                ..
            } => stamp_passive_2port(&mut c_mat, &topo.node_index, pos, neg, *capacitance),
            Element::Inductor {
                pos,
                neg,
                inductance,
                ..
            } => stamp_passive_2port(&mut l_mat, &topo.node_index, pos, neg, 1.0 / *inductance),
            _ => {}
        }
    }
    // Device-internal small-signal reactances (diode Cj, MOSFET caps, photonic
    // parasitics) — previously absent from noise; see ac.rs for the rationale.
    for dev in devices.iter() {
        for r in dev.small_signal_reactances() {
            match r.kind {
                ReactiveKind::Capacitor => stamp_2port_by_id(&mut c_mat, r.pos, r.neg, r.value),
                ReactiveKind::Inductor if r.value != 0.0 => {
                    stamp_2port_by_id(&mut l_mat, r.pos, r.neg, 1.0 / r.value)
                }
                ReactiveKind::Inductor => {}
            }
        }
        // Devices whose reactance is a general ∂q/∂x matrix rather than a set
        // of two-terminal branches (OSDI/Verilog-A) stamp it themselves.
        dev.load_reactive_jacobian(&mut c_mat);
    }

    // Locate the named input source so we can compute its signal-path gain.
    let input_vsrc_idx = topo.vsrc_index.get(input_src).copied().ok_or_else(|| {
        SimError::ParameterError(format!("noise: input source '{input_src}' not found"))
    })?;

    let out_pos_idx = if out_pos == "0" {
        None
    } else {
        topo.node_index.get(out_pos).copied()
    };
    let out_neg_idx = if out_neg == "0" {
        None
    } else {
        topo.node_index.get(out_neg).copied()
    };
    if out_pos == "0" && out_neg == "0" {
        return Err(SimError::ParameterError(
            "noise: output node cannot be ground".into(),
        ));
    }

    let sources = NoiseSources::build(netlist, &topo, opts.temp_k);

    let mut result = NoiseResult {
        freq: freqs.to_vec(),
        onoise_psd: Vec::with_capacity(freqs.len()),
        inoise_psd: Vec::with_capacity(freqs.len()),
    };

    for &f in freqs {
        let omega = 2.0 * std::f64::consts::PI * f;

        // B = ωC − L⁻¹/ω (purely real coefficient matrix; imaginary part of Y).
        let mut b_mat = vec![vec![0.0f64; size]; size];
        for i in 0..size {
            for j in 0..size {
                b_mat[i][j] = omega * c_mat[i][j] - l_mat[i][j] / omega;
            }
        }

        // Forward 2N×2N block system for the signal-path gain H(f):
        //   [ G  -B ] [V_re]   [b_re]
        //   [ B   G ] [V_im] = [b_im]
        // where the RHS injects 1 V at the input source's branch row.
        let n2 = 2 * size;
        let mut a_fwd = vec![vec![0.0f64; n2]; n2];
        for i in 0..size {
            for j in 0..size {
                a_fwd[i][j] = g_mat[i][j];
                a_fwd[i][size + j] = -b_mat[i][j];
                a_fwd[size + i][j] = b_mat[i][j];
                a_fwd[size + i][size + j] = g_mat[i][j];
            }
        }
        let mut rhs_fwd = vec![0.0f64; n2];
        rhs_fwd[n_nodes + input_vsrc_idx] = 1.0; // unit AC amplitude on V source
        let v_fwd = noise_solver.solve(&CircuitTopology::sparse_from_dense(&a_fwd), &rhs_fwd)?;
        let v_re_fwd = &v_fwd[..size];
        let v_im_fwd = &v_fwd[size..];
        let h_re = pick(v_re_fwd, out_pos_idx) - pick(v_re_fwd, out_neg_idx);
        let h_im = pick(v_im_fwd, out_pos_idx) - pick(v_im_fwd, out_neg_idx);
        let h_mag_sq = h_re * h_re + h_im * h_im;

        // Adjoint solve.  Real-block transpose of [G -B; B G] is
        //   [ G^T   B^T ]
        //   [-B^T   G^T ]
        // which, by the advisor's note, represents M^H rather than M^T.  For
        // |λ[p]−λ[n]|² in the noise sum this is fine — conjugation preserves
        // magnitude — but we MUST NOT compare λ values directly to an A^T λ
        // solution from a different convention.
        let mut a_adj = vec![vec![0.0f64; n2]; n2];
        for i in 0..size {
            for j in 0..size {
                a_adj[i][j] = g_mat[j][i];
                a_adj[i][size + j] = b_mat[j][i];
                a_adj[size + i][j] = -b_mat[j][i];
                a_adj[size + i][size + j] = g_mat[j][i];
            }
        }
        let mut rhs_adj = vec![0.0f64; n2];
        // e_out: +1 at out_pos, -1 at out_neg, both in the real block (imag
        // RHS is zero — we observe a real-valued node voltage).
        if let Some(i) = out_pos_idx {
            rhs_adj[i] += 1.0;
        }
        if let Some(i) = out_neg_idx {
            rhs_adj[i] -= 1.0;
        }
        let lam = noise_solver.solve(&CircuitTopology::sparse_from_dense(&a_adj), &rhs_adj)?;
        let lam_re = &lam[..size];
        let lam_im = &lam[size..];

        // Every generator in the circuit, from the one enumeration both
        // analyses share.  A source's taps sum BEFORE the magnitude is taken,
        // so a multi-tap generator interferes with itself instead of adding in
        // quadrature — laser RIN is the case; see `CorrelatedNoise`.
        let mut s_v_out = 0.0_f64;
        sources.for_each(&devices, &ctx, f, |taps, psd| {
            let mut z_re = 0.0;
            let mut z_im = 0.0;
            for &(p_idx, n_idx, w) in taps {
                z_re += w * (pick(lam_re, p_idx) - pick(lam_re, n_idx));
                z_im += w * (pick(lam_im, p_idx) - pick(lam_im, n_idx));
            }
            s_v_out += psd * (z_re * z_re + z_im * z_im);
        });

        let s_v_in = if h_mag_sq > 1e-30 {
            s_v_out / h_mag_sq
        } else {
            f64::NAN
        };
        result.onoise_psd.push(s_v_out);
        result.inoise_psd.push(s_v_in);
    }

    Ok(result)
}

/// Pick a node value from a vector indexed by `NodeId`-style Option<usize>.
fn pick(v: &[f64], idx: Option<usize>) -> f64 {
    match idx {
        Some(i) => v[i],
        None => 0.0,
    }
}

/// DC OP shared with ac.rs (kept private here to avoid coupling the modules).
fn run_dc_op(
    topo: &CircuitTopology,
    netlist: &Netlist,
    devices: &mut [Box<dyn Device>],
    ctx: &crate::device::SimContext,
    opts: &SimOptions,
    solver: &dyn LinearSolver,
) -> Result<Vec<f64>, SimError> {
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();
    let n_nodes = topo.n_nodes();
    // Not every unknown is a volt — see `crate::tolerance`.
    let tol = crate::tolerance::Tolerances::build(netlist, topo, opts);
    let mut x = vec![0.0f64; topo.size];

    for _ in 0..opts.itl1 {
        let mut mat = stamp_netlist_scaled(
            topo,
            netlist,
            1.0,
            &empty,
            &empty,
            crate::mna::InductorDc::Reactive,
        );
        for dev in devices.iter_mut() {
            dev.eval(&x, EvalFlags::dc(), ctx);
            dev.load_residual(&mut mat.b);
            dev.load_jacobian(&mut mat);
        }
        topo.stamp_gmin(&mut mat.a, opts.gmin, RowFloor::PinEmptyRows);
        let x_new = solver.solve(&mat.a, &mat.b)?;
        let max_dv = x_new
            .iter()
            .zip(x.iter())
            .take(n_nodes)
            .map(|(n, o)| (n - o).abs())
            .fold(0.0f64, f64::max);
        let x_next: Vec<f64> = if max_dv > opts.vmax {
            let scale = opts.vmax / max_dv;
            x.iter()
                .zip(x_new.iter())
                .map(|(o, n)| o + scale * (n - o))
                .collect()
        } else {
            x_new
        };
        let converged = tol.converged(&x_next, &x);
        x = x_next;
        if converged {
            return Ok(x);
        }
    }
    Err(SimError::NoConvergence { iters: opts.itl1 })
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

    /// Resistor-divider thermal noise.  At any frequency:
    ///   S_V_out = 4kT·(R1 ‖ R2)
    /// for a 1-V AC input driving R1→out and R2→0.  Independent of frequency
    /// in this purely resistive circuit.
    #[test]
    fn resistor_divider_thermal_noise() {
        // The current parser doesn't have an `AC` source tag; AC analysis
        // treats every voltage source with a nonzero DC value as the
        // excitation.  noise() routes its own RHS injection by name, so the
        // DC value here is irrelevant for the test.
        let net = parse_spice(
            "* thermal\nV1 in 0 DC 1\nR1 in out 1k\nR2 out 0 1k\n\
             .noise V(out) V1 DEC 1 1k 1k\n.end\n",
        )
        .unwrap();
        let mut registry = crate::device_registry::DeviceRegistry::new();
        registry.register_builtin_models(&net.models);
        let opts = SimOptions::default();
        let r = noise_analysis(&net, &[1e3], "out", "0", "v1", &registry, &opts).unwrap();

        let kt = KB * opts.temp_k;
        let expected = 4.0 * kt * 500.0; // 4kT·(R1||R2)
        let s = r.onoise_psd[0];
        let rel = (s - expected).abs() / expected;
        assert!(
            rel < 0.01,
            "S_V_out={s:.3e} expected={expected:.3e} rel={rel:.3e}"
        );
    }

    /// Diode shot noise.  A current source biases a diode to ~1 mA so its
    /// internal Id is set by the source, not by an exponential.  At low f
    /// (cap negligible) the diode's small-signal resistance r_d = V_T / Id
    /// terminates the noise current into the output node.  The expected
    /// output PSD is:
    ///     S_V_out = (R_load || r_d)² · (2qId + 4kT/R_load · (r_d/(R_load+r_d))²)
    /// We just verify shot dominates: S_V_out > 4kT/R_load · z² and the
    /// magnitude is roughly 2qId · (R_load || r_d)² within 30%.
    #[test]
    fn diode_shot_noise_dominates_at_high_bias() {
        let src = "* diode shot\n\
                   Vbias bias 0 DC 1\n\
                   Ib   0 b  1m\n\
                   D1   b 0  myd\n\
                   .model myd D (Is=1e-14 N=1)\n\
                   .noise V(b) Vbias DEC 1 1k 1k\n.end\n";
        let net = parse_spice(src).unwrap();
        let mut registry = crate::device_registry::DeviceRegistry::new();
        registry.register_builtin_models(&net.models);
        let opts = SimOptions::default();
        let r = noise_analysis(&net, &[1e3], "b", "0", "vbias", &registry, &opts).unwrap();

        // The output is the diode anode; its small-signal resistance is the
        // only impedance to ground at that node aside from gmin.  λ at the
        // diode terminals therefore ≈ r_d (= V_T / Id ≈ 25.85 Ω at 27 °C,
        // Id=1mA), so output PSD ≈ 2qId · r_d².
        const Q: f64 = 1.602176634e-19;
        let id = 1e-3;
        let vt = KB * opts.temp_k / Q;
        let r_d = vt / id;
        let expected = 2.0 * Q * id * r_d * r_d;
        let s = r.onoise_psd[0];
        assert!(
            (s - expected).abs() / expected < 0.1,
            "diode shot S_V_out={s:.3e} expected≈{expected:.3e}"
        );
    }

    /// The Nutmeg writer must emit **amplitude density (V/√Hz)**, not the
    /// V²/Hz PSD it stores, because that is what `voltage-density` /
    /// `onoise_spectrum` mean to a rawfile reader.
    ///
    /// Anchored on the analytic value — √(4kT·500Ω) ≈ 2.879 nV/√Hz — rather
    /// than on `onoise_psd`, so the assertion still holds if both the solver
    /// and the writer are wrong together.  The second half is the real guard:
    /// it fails if anyone emits the PSD under these names, which no unit tag
    /// in the file would reveal.
    #[test]
    fn write_nutmeg_emits_amplitude_density_not_psd() {
        let net = parse_spice(
            "* thermal\nV1 in 0 DC 1\nR1 in out 1k\nR2 out 0 1k\n\
             .noise V(out) V1 DEC 1 1k 1k\n.end\n",
        )
        .unwrap();
        let mut registry = crate::device_registry::DeviceRegistry::new();
        registry.register_builtin_models(&net.models);
        let opts = SimOptions::default();
        let r = noise_analysis(&net, &[1e3], "out", "0", "v1", &registry, &opts).unwrap();

        let mut buf = Vec::new();
        r.write_nutmeg(&mut buf, "thermal test").unwrap();
        let s = String::from_utf8(buf).unwrap();

        assert!(
            s.contains("Plotname: Noise Spectral Density Curves"),
            "plotname: {s}"
        );
        assert!(s.contains("Flags: real"), "flags: {s}");
        assert!(s.contains("frequency\tfrequency"), "freq var: {s}");
        assert!(
            s.contains("onoise_spectrum\tvoltage-density"),
            "onoise var: {s}"
        );
        assert!(
            s.contains("inoise_spectrum\tvoltage-density"),
            "inoise var: {s}"
        );

        // Values block: point index + frequency, then onoise, then inoise.
        let values = s.split("Values:").nth(1).unwrap();
        let rows: Vec<&str> = values.lines().filter(|l| !l.trim().is_empty()).collect();
        assert_eq!(rows.len(), 3, "one point => three lines: {values}");
        let emitted: f64 = rows[1].trim().parse().unwrap();

        let psd = 4.0 * KB * opts.temp_k * 500.0; // analytic 4kT·(R1||R2)
        let expected_vrthz = psd.sqrt(); // ≈ 2.879e-9
        assert!(
            (emitted - expected_vrthz).abs() / expected_vrthz < 0.01,
            "emitted={emitted:.4e} expected≈{expected_vrthz:.4e} V/√Hz"
        );
        // Emitting the PSD instead would be ~9 orders of magnitude off.
        assert!(
            (emitted - psd).abs() / expected_vrthz > 0.5,
            "emitted the V²/Hz PSD ({psd:.4e}) where V/√Hz was required"
        );
    }

    /// A NaN input-referred PSD must survive to the output as NaN, not become
    /// zero.  `inoise_psd` is NaN by design where the transfer function is too
    /// small to invert, and `0` there would read as a noiseless input.
    #[test]
    fn amplitude_density_preserves_nan_and_clamps_negatives() {
        assert!(amplitude_density(f64::NAN).is_nan());
        assert_eq!(amplitude_density(-1e-30), 0.0, "numerical dust clamps");
        assert_eq!(amplitude_density(4.0), 2.0);
    }
}
