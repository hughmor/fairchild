use super::spectrum::{AwgSpectrum, ChannelGrid, SpectrumTable};
use super::{stamp_potential_eq, C0};
use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

// ────────────────────────────────────────────────────────────────────────
// Arrayed-waveguide grating router (`fc_awgr`)
// ────────────────────────────────────────────────────────────────────────

/// N×N cyclic arrayed-waveguide grating router.
///
/// Input port `i` carrying channel `k` is routed to output port
/// `j = (i + k) mod N`, still in channel slot `k` — the cyclic wavelength
/// shift that makes an AWGR an all-to-all interconnect: every output port
/// receives exactly one wavelength from every input port.
///
/// # The model
///
/// Adopting the convention that **channel slot index ≡ wavelength index**, the
/// whole device is one complex matrix per slot:
///
/// ```text
///   out_j[k] = Σ_i  t_ij(λ_{i,k}) · in_i[k]
/// ```
///
/// with `t_ij` the field transmission from input `i` to output `j`. The ideal
/// router is the permutation `t_ij = 1 ⟺ (j − i) mod N == k`.
///
/// That single form captures both crosstalk mechanisms, and it is worth being
/// explicit about why, because it is the reason this device is representable
/// and a crosstalking 1×N demux is not:
///
/// - **Wrong-port crosstalk** (light from an unintended input reaching this
///   output) arrives at the *same* wavelength, so it lands in the same slot and
///   sums coherently with the intended signal. Correct, and it is the dominant
///   term in an AWG penalty budget.
/// - **Wrong-wavelength crosstalk** lands in *its own* slot and stays separate.
///
/// Neither path ever adds two different carriers into one complex envelope,
/// which would inject a spurious DC term where the physical 100 GHz beat note
/// would have been filtered out by any real photodiode.
///
/// The transmission is a *static* coefficient evaluated at each channel's
/// carrier — the exact narrowband limit. See `docs/photonic-models.md` for the
/// error bound and for what this deliberately cannot model (sideband shaping,
/// PM→AM off a detuned slope, channel skew).
///
/// # Terminals — `2·wpc·N²`
///
/// All `N` input ports, then all `N` output ports; each port contributes its
/// `N` channels in order, each channel `wpc` wires:
///
/// ```text
///   in_0[0..N] … in_{N−1}[0..N]   out_0[0..N] … out_{N−1}[0..N]
/// ```
///
/// `N = √(len / (2·wpc))` and must come out exact. Declare it with vector
/// ports and the width follows the netlist:
///
/// ```text
///   .optical_port in0 8
///   … in1 … in7, out0 … out7 …
///   Xr in0 in1 … in7 out0 … out7 fc_awgr df_ghz=100 fwhm_ghz=40 il_db=3
/// ```
///
/// Sixteen port tokens is a lot to type; wrap it in a `.subckt` PCell.
///
/// # Modes
///
/// Selected by which parameters are present, not by a mode string:
///
/// - **ideal** (nothing set) — the exact cyclic permutation, lossless. Free:
///   one coupling per output slot.
/// - **gauss** (`fwhm_ghz` set) — super-Gaussian passbands on a periodic
///   frequency grid, with crosstalk floors. See [`AwgSpectrum`].
/// - **table** (`.model … fc_awgr sfile="…"`) — measured `N×N` spectra
///   interpolated per channel. See [`SpectrumTable`].
///
/// # Parameters (gauss mode)
///
/// | name | default | meaning |
/// |---|---|---|
/// | `lambda0_nm` | `.options lambda_center_nm` | passband centre of the `(j−i) mod N == 0` pairs |
/// | `df_ghz` | 100 | channel spacing (a **frequency** grid — AWGs are periodic in f) |
/// | `fsr_ghz` | `N·df_ghz` | free spectral range; the default is the cyclic condition |
/// | `fwhm_ghz` | — | power FWHM; **a positive value selects gauss mode**, 0 stays ideal |
/// | `shape_p` | 1 | super-Gaussian order: 1 = Gaussian, 2–4 = flat-top |
/// | `il_db` | 3 (gauss) / 0 (ideal) | peak insertion loss |
/// | `il_tilt_db` | 0 | extra loss at the outermost channel vs the centre one |
/// | `xt_adj_db` | −30 | adjacent-channel crosstalk floor |
/// | `xt_bg_db` | −40 | non-adjacent crosstalk floor |
/// | `dlambda_dt_pm_per_k` | 0 | grid thermal drift (silica ≈ 11, SOI ≈ 80) |
/// | `t_nom_k` | 300.15 | reference temperature for the drift |
/// | `lambda_src` | 0 | which input port's λ tags the outputs mirror; −1 = device grid |
///
/// # Using it as a mux or demux
///
/// A demux **is** this device with `N−1` input ports left dark; a mux is it
/// with `N−1` output ports left dark. That is not an analogy — it is the same
/// physical star-coupler-plus-array used three ways, and it is why
/// `fc_mux`/`fc_demux` do not model crosstalk themselves: a demux with
/// single-channel output ports cannot represent it (two wavelengths would have
/// to share one envelope). Dark ports contribute nothing regardless of their λ
/// tags, so nothing special is needed.
///
/// # Transmission phase
///
/// Analytic modes produce **purely real** transmission — every crosstalk term
/// adds in phase, which is the pessimistic bound on coherent crosstalk. A
/// synthetic random-phase mode for sampling the crosstalk-penalty distribution
/// is not implemented; see `docs/photonic-models.md` §4, "Not implemented, deliberately". Table
/// mode does honour measured phase from `t_<i>_<j>_deg` columns.
///
/// # Not supported
///
/// Bidirectional propagation (`wpc = 5`) — the backward fields would need the
/// transposed routing, and `setup_instance` rejects it rather than leaving the
/// backward wires undriven. Latency (`tau_s`) is likewise not modelled; an AWG's
/// few-ps array transit is far below any timestep this simulator runs at.
pub struct NativeAwgr {
    n: usize,
    wpc: usize,
    /// Smallest terminal count that would have worked, recorded when
    /// `setup_instance` refuses; see the other photonic devices.
    min_terminals: Option<usize>,
    nodes: Vec<NodeId>,
    branches: Vec<Option<usize>>,

    // ── Spectral parameters (raw; resolved lazily, once N is known) ──
    lambda0_m: f64,
    df_hz: f64,
    fsr_hz: Option<f64>,
    fwhm_hz: Option<f64>,
    shape_p: f64,
    il_db: Option<f64>,
    il_tilt_db: f64,
    xt_adj_db: f64,
    xt_bg_db: f64,
    dlambda_dt_m_per_k: f64,
    t_nom_k: f64,
    lambda_src: i64,
    /// Measured-table mode, when built from a `.model … sfile="…"` card.
    table: Option<SpectrumTable>,

    // ── Cached evaluation state ──
    /// λ read from every input slot last rebuild, `[i·N + k]`.
    lam_cache: Vec<f64>,
    /// Compiled stamp list: (branch slot, output node, weighted inputs).
    rows: Vec<StampRow>,
    /// λ branch rows driven to a constant (slot, wavelength) — the fallback
    /// when there is no input tag to mirror.
    lambda_rhs: Vec<(usize, f64)>,
    built: bool,
    warned_gain: bool,
}

/// Branch rows per output channel slot: out_re, out_im, out_λ.
const BRANCHES_PER_SLOT: usize = 3;

/// One compiled output equation: which branch row, which output node, and the
/// weighted input nodes feeding it.
type StampRow = (usize, NodeId, Vec<(NodeId, f64)>);

impl Default for NativeAwgr {
    fn default() -> Self {
        Self::new()
    }
}

impl NativeAwgr {
    pub fn new() -> Self {
        Self {
            n: 0,
            wpc: 3,
            min_terminals: None,
            nodes: Vec::new(),
            branches: Vec::new(),
            lambda0_m: 1.55e-6,
            df_hz: 100e9,
            fsr_hz: None,
            fwhm_hz: None,
            shape_p: 1.0,
            il_db: None,
            il_tilt_db: 0.0,
            xt_adj_db: -30.0,
            xt_bg_db: -40.0,
            dlambda_dt_m_per_k: 0.0,
            t_nom_k: 300.15,
            lambda_src: 0,
            table: None,
            lam_cache: Vec::new(),
            rows: Vec::new(),
            lambda_rhs: Vec::new(),
            built: false,
            warned_gain: false,
        }
    }

    /// Port count, known after `setup_instance`. The registrar needs it to
    /// size a measured table against the netlist rather than the file.
    pub fn n_ports(&self) -> usize {
        self.n
    }

    /// Switch to measured-table mode. Used by the `.model … fc_awgr sfile="…"`
    /// registrar, which is the only path that can carry a string parameter.
    pub fn set_table(&mut self, table: SpectrumTable) {
        self.table = Some(table);
        self.built = false;
    }

    /// Wire index of input port `i`, slot `k`, wire `w`.
    #[inline]
    fn in_wire(&self, i: usize, k: usize, w: usize) -> NodeId {
        self.nodes[self.wpc * self.n * i + self.wpc * k + w]
    }

    /// Wire index of output port `j`, slot `k`, wire `w`.
    #[inline]
    fn out_wire(&self, j: usize, k: usize, w: usize) -> NodeId {
        self.nodes[self.wpc * self.n * (self.n + j) + self.wpc * k + w]
    }

    /// Read input port `i` slot `k`'s λ wire, bootstrapped to the reference
    /// wavelength when undriven (matching `OpticalSegment::lambda_of`).
    fn lambda_in(&self, x: &[f64], i: usize, k: usize) -> f64 {
        match self.in_wire(i, k, self.wpc - 1) {
            Some(idx) if x[idx].abs() > 1e-9 => x[idx],
            _ => self.lambda0_m,
        }
    }

    /// The channel grid implied by the current parameters.
    fn grid(&self) -> ChannelGrid {
        let df = self.df_hz;
        ChannelGrid {
            f0_hz: if self.lambda0_m > 0.0 {
                C0 / self.lambda0_m
            } else {
                0.0
            },
            df_hz: df,
            fsr_hz: self.fsr_hz.unwrap_or(self.n as f64 * df),
            fwhm_hz: self.fwhm_hz.unwrap_or(0.4 * df),
            shape_p: self.shape_p,
            dlambda_dt_m_per_k: self.dlambda_dt_m_per_k,
            t_nom_k: self.t_nom_k,
        }
    }

    /// Peak insertion loss in dB. Ideal mode is a *perfect* router and
    /// defaults to lossless; a modelled passband defaults to a realistic 3 dB,
    /// which also keeps the crosstalk floors from summing past unity and
    /// tripping the energy warning on every run.
    fn il_db_or_default(&self) -> f64 {
        self.il_db
            .unwrap_or(if self.fwhm_hz.is_some() { 3.0 } else { 0.0 })
    }

    fn spectrum(&self) -> AwgSpectrum {
        AwgSpectrum {
            grid: self.grid(),
            xt_adj: 10f64.powf(self.xt_adj_db / 10.0),
            xt_bg: 10f64.powf(self.xt_bg_db / 10.0),
            il_peak: 10f64.powf(-self.il_db_or_default() / 10.0),
            il_tilt: 10f64.powf(-self.il_tilt_db / 10.0),
        }
    }

    /// Field transmission from input `i` to output `j` for channel slot `k`,
    /// whose carrier sits at `lambda`, in the currently selected mode.
    fn t_ij(
        &self,
        spec: &AwgSpectrum,
        lambda: f64,
        i: usize,
        j: usize,
        k: usize,
        t_k: f64,
    ) -> (f64, f64) {
        if let Some(table) = &self.table {
            return table.amp(lambda, i, j);
        }
        // Channel index this port pair passes: the cyclic shift (j − i) mod N.
        let m = (j + self.n - i) % self.n;
        match self.fwhm_hz {
            // gauss: a super-Gaussian passband floored by the crosstalk spec,
            // evaluated at the carrier actually present on that input.
            Some(_) => (spec.power(lambda, m, self.n, t_k).sqrt(), 0.0),
            // ideal: a pure permutation on slot indices. Deliberately *not* a
            // grid lookup — an ideal router should route whatever comb it is
            // handed, not silently go dark because the lasers sit off the
            // device's nominal grid.
            None => {
                let amp = if m == k {
                    10f64.powf(-self.il_db_or_default() / 20.0)
                } else {
                    0.0
                };
                (amp, 0.0)
            }
        }
    }

    /// Recompute every coefficient and recompile the stamp list.
    fn rebuild(&mut self, ctx: &SimContext) {
        let (n, wpc) = (self.n, self.wpc);
        let lam_w = wpc - 1;
        let t_k = ctx.temperature;
        let spec = self.spectrum();
        self.rows.clear();
        self.lambda_rhs.clear();
        // Per-input energy accounting for the passivity warning.
        let mut out_power = vec![0.0f64; n * n];

        for j in 0..n {
            for k in 0..n {
                let slot = BRANCHES_PER_SLOT * (j * n + k);
                let mut re_ins: Vec<(NodeId, f64)> = Vec::with_capacity(2 * n);
                let mut im_ins: Vec<(NodeId, f64)> = Vec::with_capacity(2 * n);
                for i in 0..n {
                    let lambda = self.lam_cache[i * n + k];
                    let (a, b) = self.t_ij(&spec, lambda, i, j, k, t_k);
                    out_power[i * n + k] += a * a + b * b;
                    let (in_re, in_im) = (self.in_wire(i, k, 0), self.in_wire(i, k, 1));
                    // stamp_potential_eq writes V(out) + Σ c·V(in) = 0, so the
                    // couplings carry the negated matrix entries:
                    //   out_re = Σ (a·in_re − b·in_im)
                    //   out_im = Σ (b·in_re + a·in_im)
                    re_ins.push((in_re, -a));
                    re_ins.push((in_im, b));
                    im_ins.push((in_re, -b));
                    im_ins.push((in_im, -a));
                }
                self.rows.push((slot, self.out_wire(j, k, 0), re_ins));
                self.rows.push((slot + 1, self.out_wire(j, k, 1), im_ins));

                // λ tag: mirror input port `lambda_src`'s tag for this slot,
                // falling back to the device's own grid wavelength when that
                // wire does not exist in the netlist.
                //
                // The test is deliberately *structural* (does the node exist)
                // rather than "is it driven": the first rebuild happens on the
                // zero initial guess, where nothing is driven yet, so a
                // value-based test would freeze every output onto the grid and
                // throw away the input comb's detuning.
                let src_wire = usize::try_from(self.lambda_src)
                    .ok()
                    .filter(|&s| s < n)
                    .and_then(|s| self.in_wire(s, k, lam_w));
                let out_lam = self.out_wire(j, k, lam_w);
                match src_wire {
                    Some(_) => self.rows.push((slot + 2, out_lam, vec![(src_wire, -1.0)])),
                    None => {
                        self.rows.push((slot + 2, out_lam, Vec::new()));
                        self.lambda_rhs
                            .push((slot + 2, self.grid().lambda_center(k, t_k)));
                    }
                }
            }
        }
        self.built = true;
        self.warn_if_energy_created(&out_power);
    }

    /// A superposed passband and crosstalk floor can total more than unity —
    /// physically the floor is stolen *from* the passband, but the parameters
    /// are independent, so nothing enforces it. Warn once rather than assert:
    /// the excess is tiny in any realistic setup (`il_db = 3` with −30 dB
    /// floors totals 0.504), and refusing to run would be obnoxious.
    fn warn_if_energy_created(&mut self, out_power: &[f64]) {
        if self.warned_gain {
            return;
        }
        let Some((idx, &p)) = out_power
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
        else {
            return;
        };
        if p > 1.0 + 1e-9 {
            self.warned_gain = true;
            eprintln!(
                "warning: fc_awgr input {} channel {} sends {:.4}× its power to the outputs \
                 (>1). Raise il_db or lower xt_adj_db/xt_bg_db; in a feedback path this grows \
                 without bound.",
                idx / self.n,
                idx % self.n,
                p
            );
        }
    }
}

impl Device for NativeAwgr {
    fn num_terminals(&self) -> usize {
        if let Some(min) = self.min_terminals {
            return min;
        }
        self.nodes.len()
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        self.wpc = ctx.wires_per_channel();
        self.lambda0_m = ctx.lambda_center_m;
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        let wpc = ctx.wires_per_channel();
        if wpc != 3 {
            // The caller reports a terminal-count mismatch, which is true but
            // not the cause: with 5 wires per channel the netlist simply cannot
            // supply the 2·3·N² this device wants. Name the real fix here.
            eprintln!(
                "warning: fc_awgr does not support bidirectional propagation (wpc={wpc}); \
                 the backward-travelling fields would need the transposed routing. Drop \
                 `.options enable_bidirectional=1`. The terminal-count error that follows \
                 is a consequence of it."
            );
            self.min_terminals = Some(2 * 3);
            return;
        }
        self.wpc = wpc;
        self.lambda0_m = ctx.lambda_center_m;
        let len = terminals.len();
        // 2 sides × N ports × N channels × wpc wires.
        let per_side = len as f64 / (2.0 * wpc as f64);
        let n = per_side.sqrt().round() as usize;
        if n == 0 || 2 * wpc * n * n != len {
            // Not 2·wpc·N² for any N — most often because the ports do not all
            // declare the same channel count as the port count.
            self.min_terminals = Some(2 * wpc);
            return;
        }
        self.min_terminals = None;
        self.n = n;
        self.nodes = terminals.to_vec();
        self.branches = vec![None; BRANCHES_PER_SLOT * n * n];
        self.lam_cache = vec![f64::NAN; n * n];
        self.built = false;
    }

    fn num_extra_nodes(&self) -> usize {
        self.branches.len()
    }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..self.branches.len() {
            self.branches[i] = Some(first_idx + i);
        }
    }

    fn stamp_pairs(&self) -> Option<Vec<(usize, usize)>> {
        // Declared explicitly because the default footprint is a clique over
        // every row the device owns, and this device owns 9N² of them: at N=8
        // that is a 332k-entry clique standing in for a true footprint of ~6k,
        // and it grows as N⁴ rather than N³.
        let (n, wpc) = (self.n, self.wpc);
        let mut pairs = Vec::with_capacity(12 * n * n * n);
        for j in 0..n {
            for k in 0..n {
                let base = BRANCHES_PER_SLOT * (j * n + k);
                // Output wire w is driven by branch row `base + w` (re, im, λ).
                for w in 0..BRANCHES_PER_SLOT {
                    let (Some(row), Some(out)) = (self.branches[base + w], self.out_wire(j, k, w))
                    else {
                        continue;
                    };
                    // The branch row's own potential equation, both directions.
                    pairs.push((row, out));
                    pairs.push((out, row));
                    pairs.push((row, row));
                }
                // Every input this output slot can draw from. Both field rows
                // touch both field wires — a complex coefficient mixes re and
                // im, so declaring only the diagonal would silently drop the
                // off-diagonal cells the moment a transmission has phase.
                for i in 0..n {
                    for rowoff in 0..2 {
                        let Some(row) = self.branches[base + rowoff] else {
                            continue;
                        };
                        for w in 0..2 {
                            if let Some(inp) = self.in_wire(i, k, w) {
                                pairs.push((row, inp));
                            }
                        }
                    }
                    if let (Some(row), Some(inp)) =
                        (self.branches[base + 2], self.in_wire(i, k, wpc - 1))
                    {
                        pairs.push((row, inp));
                    }
                }
            }
        }
        pairs.sort_unstable();
        pairs.dedup();
        Some(pairs)
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        let hit = match name.to_lowercase().as_str() {
            "lambda0_nm" | "wavelength_nm" => {
                self.lambda0_m = value * 1e-9;
                true
            }
            "lambda0_m" => {
                self.lambda0_m = value;
                true
            }
            "df_ghz" | "spacing_ghz" => {
                self.df_hz = value * 1e9;
                true
            }
            "fsr_ghz" => {
                self.fsr_hz = Some(value * 1e9);
                true
            }
            "fwhm_ghz" | "bw_ghz" => {
                // A non-positive width means "no passband modelled" — i.e.
                // ideal mode — rather than a zero-width filter that passes
                // nothing. Lets a generated netlist pass fwhm_ghz=0 to mean
                // "leave it ideal" instead of going silently dark.
                self.fwhm_hz = (value > 0.0).then_some(value * 1e9);
                true
            }
            "shape_p" => {
                self.shape_p = value.max(0.1);
                true
            }
            "il_db" => {
                self.il_db = Some(value);
                true
            }
            "il_tilt_db" => {
                self.il_tilt_db = value;
                true
            }
            "xt_adj_db" => {
                self.xt_adj_db = value;
                true
            }
            "xt_bg_db" => {
                self.xt_bg_db = value;
                true
            }
            "dlambda_dt_pm_per_k" => {
                self.dlambda_dt_m_per_k = value * 1e-12;
                true
            }
            "t_nom_k" => {
                self.t_nom_k = value;
                true
            }
            "lambda_src" => {
                self.lambda_src = value.round() as i64;
                true
            }
            _ => false,
        };
        if hit {
            self.built = false; // parameters changed → coefficients are stale
        }
        hit
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, ctx: &SimContext) {
        // Coefficients depend only on the input λ tags (and on parameters,
        // which invalidate directly). Lasers are CW, so after the first NR
        // iteration this is a comparison and nothing else.
        let n = self.n;
        let mut changed = !self.built;
        for i in 0..n {
            for k in 0..n {
                let lam = self.lambda_in(x, i, k);
                let slot = i * n + k;
                // The cache starts as NaN, so the first pass always rebuilds.
                let prev = self.lam_cache[slot];
                if !prev.is_finite() || (lam - prev).abs() > 1e-18 {
                    self.lam_cache[slot] = lam;
                    changed = true;
                }
            }
        }
        if changed {
            self.rebuild(ctx);
        }
    }

    fn load_residual(&self, b: &mut [f64]) {
        for &(slot, lambda) in &self.lambda_rhs {
            if let Some(row) = self.branches[slot] {
                b[row] += lambda;
            }
        }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        for (slot, out, ins) in &self.rows {
            stamp_potential_eq(mat, &self.branches, *slot, *out, ins);
        }
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The routing rule, checked directly on the coefficient builder: input `i`
    /// channel `k` must reach output `(i + k) mod N` and nowhere else.
    #[test]
    fn ideal_mode_is_the_cyclic_permutation() {
        let n = 4;
        let mut d = NativeAwgr::new();
        d.n = n;
        d.lambda0_m = 1.55e-6;
        let spec = d.spectrum();
        let grid = d.grid();
        for k in 0..n {
            let lambda = grid.lambda_center(k, 300.15);
            for i in 0..n {
                let mut hits = Vec::new();
                for j in 0..n {
                    let (a, _) = d.t_ij(&spec, lambda, i, j, k, 300.15);
                    if a.abs() > 1e-12 {
                        hits.push(j);
                    }
                }
                assert_eq!(
                    hits,
                    vec![(i + k) % n],
                    "input {i} channel {k} should reach only output {}",
                    (i + k) % n
                );
            }
        }
    }

    /// Every output port must receive exactly one wavelength from every input
    /// port — the defining property of a cyclic router.
    #[test]
    fn every_output_sees_every_input_exactly_once() {
        let n = 8;
        let mut d = NativeAwgr::new();
        d.n = n;
        d.lambda0_m = 1.55e-6;
        let spec = d.spectrum();
        let grid = d.grid();
        for j in 0..n {
            let mut sources: Vec<usize> = Vec::new();
            for k in 0..n {
                let lambda = grid.lambda_center(k, 300.15);
                for i in 0..n {
                    if d.t_ij(&spec, lambda, i, j, k, 300.15).0.abs() > 1e-12 {
                        sources.push(i);
                    }
                }
            }
            sources.sort_unstable();
            assert_eq!(sources, (0..n).collect::<Vec<_>>(), "output port {j}");
        }
    }

    /// Gauss mode at the grid centres must reproduce the same permutation, with
    /// the crosstalk terms sitting on their specified floors.
    #[test]
    fn gauss_mode_hits_the_specified_insertion_loss_and_floors() {
        let n = 4;
        let mut d = NativeAwgr::new();
        d.n = n;
        d.lambda0_m = 1.55e-6;
        d.fwhm_hz = Some(40e9);
        d.il_db = Some(3.0);
        d.xt_adj_db = -30.0;
        d.xt_bg_db = -40.0;
        let spec = d.spectrum();
        let grid = d.grid();
        let c = 1usize; // probe at the centre of channel 1
        let lambda = grid.lambda_center(c, 300.15);
        for i in 0..n {
            for j in 0..n {
                let (a, b) = d.t_ij(&spec, lambda, i, j, c, 300.15);
                let p = a * a + b * b;
                // How far this pair's own passband sits from the probe, in
                // channels, folded into the FSR.
                let m = (j + n - i) % n;
                let off = super::super::spectrum::wrap_half(c as f64 - m as f64, n as f64).abs();
                // 3 dB is 0.50119, not 0.5 — the whole point of dB params.
                let want = 10f64.powf(-3.0 / 10.0)
                    * match off.round() as i64 {
                        0 => 1.0,  // routed: peak, i.e. the full −3 dB
                        1 => 1e-3, // adjacent-channel floor
                        _ => 1e-4, // background floor
                    };
                assert!(
                    (p / want - 1.0).abs() < 1e-6,
                    "in {i} → out {j} (m={m}, {off} channels off): got {p:.3e}, want {want:.3e}"
                );
            }
        }
    }

    /// Two crosstalk contributions landing in the same output slot must add as
    /// fields, not as powers — the coherent convention this device is built on.
    #[test]
    fn same_slot_contributions_are_coherent() {
        let n = 4;
        let mut d = NativeAwgr::new();
        d.n = n;
        d.lambda0_m = 1.55e-6;
        d.fwhm_hz = Some(40e9);
        d.xt_bg_db = -40.0;
        d.xt_adj_db = -40.0;
        let spec = d.spectrum();
        let lambda = d.grid().lambda_center(0, 300.15);
        // Inputs 1 and 2 both leak into output 0 at channel 0 (only input 0 is
        // routed there). Equal real amplitudes → the pair carries 4× one's power.
        let a1 = d.t_ij(&spec, lambda, 1, 0, 0, 300.15).0;
        let a2 = d.t_ij(&spec, lambda, 2, 0, 0, 300.15).0;
        assert!((a1 - a2).abs() < 1e-15, "{a1} vs {a2}");
        let coherent = (a1 + a2).powi(2);
        assert!(
            (coherent / (4.0 * a1 * a1) - 1.0).abs() < 1e-12,
            "fields must add, not powers"
        );
    }

    /// Ideal mode must not consult the grid: an off-grid comb still routes.
    /// The alternative (matching λ against the nominal grid) fails silently and
    /// completely the moment a user's lasers sit somewhere else.
    #[test]
    fn ideal_mode_routes_a_comb_that_is_nowhere_near_the_nominal_grid() {
        let n = 4;
        let mut d = NativeAwgr::new();
        d.n = n;
        d.lambda0_m = 1.55e-6;
        let spec = d.spectrum();
        for k in 0..n {
            for i in 0..n {
                let (a, _) = d.t_ij(&spec, 1.31e-6, i, (i + k) % n, k, 300.15);
                assert!((a - 1.0).abs() < 1e-12, "O-band comb must still route");
            }
        }
    }
}
