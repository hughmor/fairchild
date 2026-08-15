//! Shared spectral-response model for wavelength-selective photonic devices.
//!
//! Every filter in the tree — AWG router, WDM mux/demux — is described by the
//! same three things: a periodic channel grid, a super-Gaussian passband shape,
//! and a crosstalk floor. This module owns that description so the devices only
//! own their port wiring.
//!
//! # Why the grid is in frequency
//!
//! An arrayed-waveguide grating is periodic in *optical frequency*, not
//! wavelength: its free spectral range is set by the path-length increment
//! `ΔL` between adjacent arms, `FSR = c / (n_g·ΔL)`, a constant in Hz. Channel
//! spacings on the ITU grid are likewise defined in GHz. Working in Hz makes
//! the FSR wrap one `rem_euclid` instead of a transcendental, and makes a
//! cyclic router exactly cyclic rather than cyclic to first order.
//!
//! # Narrowband validity
//!
//! These functions return a *static* transmission evaluated at a channel's
//! carrier. That is the exact narrowband limit of the true baseband
//! convolution `B_k = (h_k ⊛ A_k)`, with relative field error
//! `≈ 2·ln2·(B/FWHM)²` for a centred carrier of modulation bandwidth `B`
//! (1.4 % at `B = FWHM/10`), plus a first-order amplitude tilt
//! `≈ 4·ln2·Δ·B/FWHM²` when the carrier is detuned by `Δ`. It therefore
//! reproduces insertion loss, detuning penalty and crosstalk exactly, and
//! cannot reproduce sideband shaping (ISI), PM→AM conversion off a detuned
//! passband slope, or differential channel skew. See `docs/photonic-models.md` §4.

use super::C0;

/// Fold `x` into the half-open interval `(−period/2, +period/2]`.
///
/// Used to reduce an optical-frequency offset into one free spectral range,
/// which is what makes a passband repeat every FSR.
#[inline]
pub fn wrap_half(x: f64, period: f64) -> f64 {
    if period > 0.0 && x.is_finite() {
        (x + 0.5 * period).rem_euclid(period) - 0.5 * period
    } else {
        x
    }
}

/// Super-Gaussian power transmission, normalised to 1 at `delta = 0` and
/// 0.5 at `delta = ±fwhm/2` for every order.
///
/// `p = 1` is the ordinary Gaussian `exp(−4·ln2·δ²/FWHM²)` that a
/// waveguide-input AWG produces (the passband is the overlap integral of two
/// near-Gaussian mode fields). `p ≥ 2` flattens the top and steepens the
/// skirts, which is the shape of an MMI- or parabolic-horn-input AWG.
#[inline]
pub fn super_gaussian(delta: f64, fwhm: f64, p: f64) -> f64 {
    if fwhm <= 0.0 || fwhm.is_nan() {
        // Zero width: a delta function. Only exact centre transmits.
        return if delta == 0.0 { 1.0 } else { 0.0 };
    }
    let u = (2.0 * delta / fwhm).abs();
    let e = -std::f64::consts::LN_2 * u.powf(2.0 * p.max(0.1));
    if e < -700.0 {
        0.0 // exp() would underflow to 0 anyway; skip it and stay finite
    } else {
        e.exp()
    }
}

/// A periodic WDM channel grid with a super-Gaussian passband per channel.
#[derive(Clone, Debug)]
pub struct ChannelGrid {
    /// Optical frequency of channel index 0 at `t_nom_k` (Hz).
    pub f0_hz: f64,
    /// Channel spacing (Hz).
    pub df_hz: f64,
    /// Free spectral range (Hz). For a cyclic N×N router this is `N·df_hz`.
    pub fsr_hz: f64,
    /// Power-transmission FWHM of one passband (Hz).
    pub fwhm_hz: f64,
    /// Super-Gaussian order — 1 = Gaussian, 2–4 = flat-top.
    pub shape_p: f64,
    /// Thermo-optic grid drift (m/K). Silica AWG ≈ 11 pm/K, SOI ≈ 80 pm/K.
    pub dlambda_dt_m_per_k: f64,
    /// Reference temperature for the drift (K).
    pub t_nom_k: f64,
}

impl ChannelGrid {
    /// Centre frequency of channel `m` at temperature `t_k`.
    ///
    /// Thermal drift is specified the way datasheets quote it — as a wavelength
    /// shift of the anchor channel — and applied by moving the anchor, leaving
    /// the channel spacing alone. Both halves of that are physical: the AWG
    /// resonance `m·λ = n_eff·ΔL` moves with `dn_eff/dT`, while the spacing is
    /// set by the diffraction order and does not meaningfully drift. Moving the
    /// anchor in wavelength (rather than converting to a rigid `Δf` via
    /// `c/λ₀²`) also makes the quoted pm/K exact instead of correct only to
    /// first order in `Δλ/λ₀`.
    pub fn f_center(&self, m: usize, t_k: f64) -> f64 {
        let mut f0 = self.f0_hz;
        if self.dlambda_dt_m_per_k != 0.0 && f0 > 0.0 {
            let lambda0 = C0 / f0 + self.dlambda_dt_m_per_k * (t_k - self.t_nom_k);
            f0 = if lambda0 > 0.0 { C0 / lambda0 } else { 0.0 };
        }
        f0 + (m as f64) * self.df_hz
    }

    /// Wavelength (m) of channel `m` at temperature `t_k` — the label an
    /// output port carries when it has no input λ tag to mirror.
    pub fn lambda_center(&self, m: usize, t_k: f64) -> f64 {
        let f = self.f_center(m, t_k);
        if f > 0.0 {
            C0 / f
        } else {
            0.0
        }
    }
}

/// Full AWG-style response: a channel grid plus the crosstalk floors that a
/// real device's array phase errors impose on the Gaussian tail.
///
/// Without the floors this model is wildly optimistic — a Gaussian tail three
/// channels out is below −1000 dB, whereas a fabricated AWG floors at −25 to
/// −35 dB adjacent and a few dB below that non-adjacent. Those two numbers are
/// what every AWG datasheet quotes, so they are what this takes.
#[derive(Clone, Debug)]
pub struct AwgSpectrum {
    pub grid: ChannelGrid,
    /// Adjacent-channel crosstalk floor, linear power (e.g. 1e-3 for −30 dB).
    pub xt_adj: f64,
    /// Non-adjacent (background) crosstalk floor, linear power.
    pub xt_bg: f64,
    /// Peak insertion loss, linear power (e.g. 0.5 for 3 dB).
    pub il_peak: f64,
    /// Extra loss at the outermost channel relative to the centre one, linear
    /// power. Models the star coupler's far-field envelope roll-off, which is
    /// why edge channels of a real AWG are the lossy ones.
    pub il_tilt: f64,
}

impl AwgSpectrum {
    /// Power transmission (0..1) for a port pair whose ideal channel index is
    /// `m`, evaluated at wavelength `lambda_m`, on an `n_ch`-channel grid.
    ///
    /// Returns 0 for a non-physical wavelength, which is what an undriven
    /// (dark) port's λ wire reads before any laser has driven it.
    pub fn power(&self, lambda_m: f64, m: usize, n_ch: usize, t_k: f64) -> f64 {
        if lambda_m <= 0.0 || !lambda_m.is_finite() {
            return 0.0;
        }
        let f = C0 / lambda_m;
        if !self.in_band(f, n_ch, t_k) {
            return 0.0;
        }
        let delta = wrap_half(f - self.grid.f_center(m, t_k), self.grid.fsr_hz);
        let shape = super_gaussian(delta, self.grid.fwhm_hz, self.grid.shape_p);
        // How many channel slots away from this pair's passband we are. The
        // floor is what the Gaussian tail saturates onto, so it must not lift
        // the in-band response — hence no floor at offset 0.
        let n_off = if self.grid.df_hz > 0.0 {
            (delta / self.grid.df_hz).round().abs() as i64
        } else {
            0
        };
        let floor = match n_off {
            0 => 0.0,
            1 => self.xt_adj,
            _ => self.xt_bg,
        };
        self.il(m, n_ch) * shape.max(floor)
    }

    /// Whether `f` is inside the device's usable band.
    ///
    /// The FSR periodicity is real but **not unbounded**: the star coupler's
    /// far-field envelope rolls the response off after a few free spectral
    /// ranges, and a physical AWG passes nothing an octave away. Modelling the
    /// periodicity as infinite is not merely optimistic, it wrecks the solver —
    /// Newton's early iterates put λ at ~1e-8 m, and an unbounded wrap folds
    /// that straight back onto a passband, so the coefficients thrash from one
    /// iteration to the next and the line search collapses to minimum steps.
    ///
    /// The window is one FSR beyond each end of the channel grid, which leaves
    /// the adjacent replicas usable (a cyclic router is often operated one FSR
    /// up) while making nonsense wavelengths honestly dark. An aperiodic filter
    /// bank (`fsr_hz = 0`) gets the grid span as its margin instead.
    fn in_band(&self, f: f64, n_ch: usize, t_k: f64) -> bool {
        let margin = if self.grid.fsr_hz > 0.0 {
            self.grid.fsr_hz
        } else {
            (n_ch.max(1) as f64) * self.grid.df_hz
        };
        let lo = self.grid.f_center(0, t_k) - margin;
        let hi = self.grid.f_center(n_ch.saturating_sub(1), t_k) + margin;
        f >= lo && f <= hi
    }

    /// Peak transmission of channel `m`, including the port non-uniformity.
    fn il(&self, m: usize, n_ch: usize) -> f64 {
        if self.il_tilt >= 1.0 || n_ch < 2 {
            return self.il_peak;
        }
        let mid = (n_ch - 1) as f64 / 2.0;
        let u = ((m as f64) - mid).abs() / mid; // 0 at grid centre, 1 at the edges
        self.il_peak * self.il_tilt.powf(u)
    }
}

// ─── Measured-table mode ────────────────────────────────────────────────────

/// An `N×N` grid of measured transmission spectra, interpolated at run time.
///
/// CSV layout — column 0 is `wavelength_nm`, then one column per port pair:
///
/// ```text
/// wavelength_nm,t_0_0_db,t_0_1_db,…,t_0_0_deg,…
/// 1545.00,-3.10,-31.4,…,0.0,…
/// ```
///
/// `t_<in>_<out>_db` is required for every pair; `t_<in>_<out>_deg` is
/// optional and defaults to zero phase. Rows may arrive in any order and are
/// sorted on load. Missing pairs read as −∞ dB (dark), so a partially measured
/// device is usable without editing the file.
///
/// Interpolation is **linear in dB and in unwrapped degrees**, clamped to the
/// endpoint values outside the measured span. Linear beats cubic here: a
/// spline rings on the steep passband skirts and can overshoot into negative
/// power, which then shows up as a NaN amplitude three devices downstream.
#[derive(Clone, Debug)]
pub struct SpectrumTable {
    n: usize,
    /// Ascending wavelengths (m).
    lambda_m: Vec<f64>,
    /// `[pair][point]` power dB, pair index = `i·n + j`.
    db: Vec<Vec<f64>>,
    /// `[pair][point]` unwrapped phase in radians.
    phase: Vec<Vec<f64>>,
}

impl SpectrumTable {
    /// Parse the CSV described on [`SpectrumTable`]. `n` is the port count,
    /// taken from the device's terminal arithmetic rather than the file, so a
    /// file that disagrees with the netlist is an error rather than a resize.
    pub fn from_csv(text: &str, n: usize) -> Result<Self, String> {
        let mut lines = text
            .lines()
            .map(str::trim)
            .filter(|l| !l.is_empty() && !l.starts_with('#') && !l.starts_with('*'));
        let header = lines.next().ok_or("empty S-parameter file")?;
        // Map each column index to (pair, is_phase). Column 0 is wavelength.
        let mut cols: Vec<Option<(usize, bool)>> = Vec::new();
        for (c, raw) in header.split(',').enumerate() {
            let name = raw.trim().to_lowercase();
            if c == 0 {
                if !name.starts_with("wavelength") && !name.starts_with("lambda") {
                    return Err(format!(
                        "first column must be wavelength_nm; got '{}'",
                        raw.trim()
                    ));
                }
                cols.push(None);
                continue;
            }
            cols.push(parse_pair_column(&name, n)?);
        }
        let npair = n * n;
        let mut lambda_m: Vec<f64> = Vec::new();
        let mut db: Vec<Vec<f64>> = vec![Vec::new(); npair];
        let mut phase: Vec<Vec<f64>> = vec![Vec::new(); npair];
        for (row, line) in lines.enumerate() {
            let fields: Vec<&str> = line.split(',').collect();
            let parse = |s: &str| -> Result<f64, String> {
                s.trim()
                    .parse::<f64>()
                    .map_err(|_| format!("row {}: '{}' is not a number", row + 2, s.trim()))
            };
            lambda_m.push(parse(fields.first().copied().unwrap_or(""))? * 1e-9);
            // Absent pairs stay dark; absent phase stays zero.
            for v in db.iter_mut() {
                v.push(f64::NEG_INFINITY);
            }
            for v in phase.iter_mut() {
                v.push(0.0);
            }
            let last = lambda_m.len() - 1;
            for (c, field) in fields.iter().enumerate().skip(1) {
                let Some(Some((pair, is_phase))) = cols.get(c) else {
                    continue;
                };
                let v = parse(field)?;
                if *is_phase {
                    phase[*pair][last] = v.to_radians();
                } else {
                    db[*pair][last] = v;
                }
            }
        }
        if lambda_m.len() < 2 {
            return Err("need at least two wavelength points to interpolate".into());
        }
        // Sort by wavelength, carrying every column along.
        let mut order: Vec<usize> = (0..lambda_m.len()).collect();
        order.sort_by(|&a, &b| lambda_m[a].total_cmp(&lambda_m[b]));
        let permute = |v: &Vec<f64>| -> Vec<f64> { order.iter().map(|&i| v[i]).collect() };
        let lambda_sorted: Vec<f64> = order.iter().map(|&i| lambda_m[i]).collect();
        let db: Vec<Vec<f64>> = db.iter().map(permute).collect();
        let mut phase: Vec<Vec<f64>> = phase.iter().map(permute).collect();
        for p in phase.iter_mut() {
            unwrap_phase(p);
        }
        Ok(SpectrumTable {
            n,
            lambda_m: lambda_sorted,
            db,
            phase,
        })
    }

    /// Field amplitude `(re, im)` from input port `i` to output port `j` at
    /// `lambda_m`, interpolated from the table.
    pub fn amp(&self, lambda_m: f64, i: usize, j: usize) -> (f64, f64) {
        if i >= self.n || j >= self.n || lambda_m <= 0.0 || lambda_m.is_nan() {
            return (0.0, 0.0);
        }
        // Outside the measured span the value is held (a laser a hair past the
        // last sweep point should not fall off a cliff), but only for a 10 %
        // margin — beyond that it is dark. Holding indefinitely would report a
        // full passband at the nonsense wavelengths Newton visits on its way
        // out of the zero initial guess, which stalls the line search.
        let span = self.lambda_m[self.lambda_m.len() - 1] - self.lambda_m[0];
        if lambda_m < self.lambda_m[0] - 0.1 * span
            || lambda_m > self.lambda_m[self.lambda_m.len() - 1] + 0.1 * span
        {
            return (0.0, 0.0);
        }
        let pair = i * self.n + j;
        let (lo, hi, w) = self.bracket(lambda_m);
        let db = lerp(self.db[pair][lo], self.db[pair][hi], w);
        if !db.is_finite() {
            return (0.0, 0.0); // NEG_INFINITY dB — an unmeasured pair
        }
        let ph = lerp(self.phase[pair][lo], self.phase[pair][hi], w);
        let mag = 10f64.powf(db / 20.0); // power dB → field amplitude
        (mag * ph.cos(), mag * ph.sin())
    }

    /// Indices bracketing `lambda_m` and the interpolation weight, clamped at
    /// both ends (hold, never extrapolate — extrapolating dB linearly off the
    /// end of a passband skirt reaches +30 dB of gain within a nanometre).
    fn bracket(&self, lambda_m: f64) -> (usize, usize, f64) {
        let xs = &self.lambda_m;
        match xs.partition_point(|&x| x < lambda_m) {
            0 => (0, 0, 0.0),
            k if k >= xs.len() => (xs.len() - 1, xs.len() - 1, 0.0),
            k => {
                let (a, b) = (xs[k - 1], xs[k]);
                let w = if b > a { (lambda_m - a) / (b - a) } else { 0.0 };
                (k - 1, k, w)
            }
        }
    }
}

/// `"t_2_5_db"` → `Some((2·n + 5, false))`. Unrecognised columns are ignored
/// (returns `Ok(None)`) so extra bookkeeping columns don't break the file.
fn parse_pair_column(name: &str, n: usize) -> Result<Option<(usize, bool)>, String> {
    let Some(rest) = name.strip_prefix("t_") else {
        return Ok(None);
    };
    let parts: Vec<&str> = rest.split('_').collect();
    if parts.len() != 3 {
        return Ok(None);
    }
    let is_phase = match parts[2] {
        "db" => false,
        "deg" => true,
        _ => return Ok(None),
    };
    let (Ok(i), Ok(j)) = (parts[0].parse::<usize>(), parts[1].parse::<usize>()) else {
        return Ok(None);
    };
    if i >= n || j >= n {
        return Err(format!(
            "column 't_{i}_{j}_…' is outside the {n}×{n} port range implied by the netlist"
        ));
    }
    Ok(Some((i * n + j, is_phase)))
}

/// Remove 2π jumps in place, so linear interpolation across a wrap doesn't
/// sweep the phasor the long way round.
fn unwrap_phase(p: &mut [f64]) {
    let two_pi = 2.0 * std::f64::consts::PI;
    for i in 1..p.len() {
        let d = p[i] - p[i - 1];
        p[i] -= two_pi * (d / two_pi).round();
    }
}

#[inline]
fn lerp(a: f64, b: f64, w: f64) -> f64 {
    a + (b - a) * w
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn super_gaussian_is_half_at_half_fwhm_for_every_order() {
        for p in [1.0, 2.0, 3.0, 5.0] {
            let t = super_gaussian(0.5, 1.0, p);
            assert!((t - 0.5).abs() < 1e-12, "p={p}: {t}");
            assert!((super_gaussian(0.0, 1.0, p) - 1.0).abs() < 1e-12);
        }
        // p=1 must be the plain Gaussian.
        let d = 0.31;
        let want = (-4.0 * std::f64::consts::LN_2 * d * d).exp();
        assert!((super_gaussian(d, 1.0, 1.0) - want).abs() < 1e-12);
        // Higher order = flatter top, steeper skirt.
        assert!(super_gaussian(0.25, 1.0, 4.0) > super_gaussian(0.25, 1.0, 1.0));
        assert!(super_gaussian(0.9, 1.0, 4.0) < super_gaussian(0.9, 1.0, 1.0));
    }

    #[test]
    fn wrap_half_folds_into_a_symmetric_window() {
        assert!((wrap_half(0.0, 10.0) - 0.0).abs() < 1e-12);
        assert!((wrap_half(11.0, 10.0) - 1.0).abs() < 1e-12);
        assert!((wrap_half(-11.0, 10.0) + 1.0).abs() < 1e-12);
        assert!((wrap_half(-4.0, 10.0) + 4.0).abs() < 1e-12);
        assert!(wrap_half(4.0, 10.0).abs() <= 5.0);
    }

    /// The crosstalk floor must catch the tail without lifting the passband.
    #[test]
    fn crosstalk_floor_applies_outside_the_passband_only() {
        let s = AwgSpectrum {
            grid: ChannelGrid {
                f0_hz: 193.4e12,
                df_hz: 100e9,
                fsr_hz: 800e9,
                fwhm_hz: 40e9,
                shape_p: 1.0,
                dlambda_dt_m_per_k: 0.0,
                t_nom_k: 300.15,
            },
            xt_adj: 1e-3,
            xt_bg: 1e-4,
            il_peak: 1.0,
            il_tilt: 1.0,
        };
        let lam = |f: f64| C0 / f;
        // Dead centre: unity, floor did not intrude.
        assert!((s.power(lam(193.4e12), 0, 8, 300.15) - 1.0).abs() < 1e-12);
        // One channel off: the Gaussian tail there is ~1e-30, so we see the floor.
        assert!((s.power(lam(193.5e12), 0, 8, 300.15) - 1e-3).abs() < 1e-9);
        // Three channels off: background floor.
        assert!((s.power(lam(193.7e12), 0, 8, 300.15) - 1e-4).abs() < 1e-9);
        // One FSR away is the same passband again (this is what makes it cyclic).
        assert!((s.power(lam(193.4e12 + 800e9), 0, 8, 300.15) - 1.0).abs() < 1e-12);
        // Non-physical λ is dark, not NaN — dark ports read λ = 0.
        assert_eq!(s.power(0.0, 0, 8, 300.15), 0.0);
    }

    #[test]
    fn fwhm_lands_at_the_three_db_points() {
        let s = AwgSpectrum {
            grid: ChannelGrid {
                f0_hz: 193.4e12,
                df_hz: 100e9,
                fsr_hz: 800e9,
                fwhm_hz: 40e9,
                shape_p: 1.0,
                dlambda_dt_m_per_k: 0.0,
                t_nom_k: 300.15,
            },
            xt_adj: 0.0,
            xt_bg: 0.0,
            il_peak: 0.5,
            il_tilt: 1.0,
        };
        let t = s.power(C0 / (193.4e12 + 20e9), 0, 8, 300.15);
        assert!((t - 0.25).abs() < 1e-12, "half of the 3 dB peak: {t}");
    }

    #[test]
    fn thermal_drift_moves_the_grid_by_the_quoted_pm_per_k() {
        let mut g = ChannelGrid {
            f0_hz: C0 / 1550e-9,
            df_hz: 100e9,
            fsr_hz: 800e9,
            fwhm_hz: 40e9,
            shape_p: 1.0,
            dlambda_dt_m_per_k: 11e-12,
            t_nom_k: 300.0,
        };
        let shifted = g.lambda_center(0, 310.0);
        assert!(
            (shifted - (1550e-9 + 110e-12)).abs() < 1e-18,
            "10 K × 11 pm/K = 110 pm exactly; got {} nm",
            shifted * 1e9
        );
        g.dlambda_dt_m_per_k = 0.0;
        assert!((g.lambda_center(0, 400.0) - 1550e-9).abs() < 1e-18);
    }

    #[test]
    fn edge_channels_carry_the_non_uniformity_tilt() {
        let s = AwgSpectrum {
            grid: ChannelGrid {
                f0_hz: 193.4e12,
                df_hz: 100e9,
                fsr_hz: 800e9,
                fwhm_hz: 40e9,
                shape_p: 1.0,
                dlambda_dt_m_per_k: 0.0,
                t_nom_k: 300.15,
            },
            xt_adj: 0.0,
            xt_bg: 0.0,
            il_peak: 1.0,
            il_tilt: 0.5, // 3 dB extra at the outermost channel
        };
        // n_ch = 5 → centre index 2 is untilted, indices 0 and 4 take the full hit.
        let centre = s.power(C0 / s.grid.f_center(2, 300.15), 2, 5, 300.15);
        let edge = s.power(C0 / s.grid.f_center(0, 300.15), 0, 5, 300.15);
        assert!((centre - 1.0).abs() < 1e-12, "{centre}");
        assert!((edge - 0.5).abs() < 1e-12, "{edge}");
    }

    #[test]
    fn csv_table_round_trips_magnitude_and_phase() {
        let csv = "\
wavelength_nm,t_0_0_db,t_0_1_db,t_1_0_db,t_1_1_db,t_0_0_deg
1550.0,-3.0,-30.0,-30.0,-3.0,90.0
1551.0,-6.0,-30.0,-30.0,-6.0,90.0
";
        let t = SpectrumTable::from_csv(csv, 2).unwrap();
        // On a grid point: −3 dB power → 0.7079 field, at +90° → purely imaginary.
        let (re, im) = t.amp(1550e-9, 0, 0);
        assert!(re.abs() < 1e-12, "re={re}");
        assert!((im - 10f64.powf(-3.0 / 20.0)).abs() < 1e-12, "im={im}");
        // Halfway: linear in dB → −4.5 dB.
        let (_, im) = t.amp(1550.5e-9, 0, 0);
        assert!((im - 10f64.powf(-4.5 / 20.0)).abs() < 1e-12);
        // Just off the end of the measured span: held, not extrapolated (dB
        // extrapolated linearly off a passband skirt reaches gain in a nm).
        let (_, im_lo) = t.amp(1549.95e-9, 0, 0);
        assert!((im_lo - 10f64.powf(-3.0 / 20.0)).abs() < 1e-12);
        let (_, im_hi) = t.amp(1551.05e-9, 0, 0);
        assert!((im_hi - 10f64.powf(-6.0 / 20.0)).abs() < 1e-12);
        // Far outside it, dark — holding a passband value out to arbitrary
        // wavelengths is what stalls Newton on its way up from x = 0.
        assert_eq!(t.amp(1500e-9, 0, 0), (0.0, 0.0));
        assert_eq!(t.amp(1600e-9, 0, 0), (0.0, 0.0));
        // A pair with no column at all is dark rather than an error.
        let t1 = SpectrumTable::from_csv("wavelength_nm,t_0_0_db\n1550.0,-3.0\n1551.0,-3.0\n", 2)
            .unwrap();
        assert_eq!(t1.amp(1550e-9, 1, 1), (0.0, 0.0));
    }

    #[test]
    fn csv_rows_may_arrive_out_of_order() {
        let a = SpectrumTable::from_csv("wavelength_nm,t_0_0_db\n1551.0,-6.0\n1550.0,-3.0\n", 1)
            .unwrap();
        let (re, _) = a.amp(1550.5e-9, 0, 0);
        assert!((re - 10f64.powf(-4.5 / 20.0)).abs() < 1e-12, "re={re}");
    }

    #[test]
    fn csv_rejects_a_pair_index_outside_the_netlist_port_count() {
        let err = SpectrumTable::from_csv("wavelength_nm,t_0_4_db\n1550.0,-3.0\n1551.0,-3.0\n", 2)
            .unwrap_err();
        assert!(err.contains("2×2"), "{err}");
    }

    #[test]
    fn phase_unwrap_keeps_interpolation_on_the_short_arc() {
        // 179° → −179° is a 2° step, not a 358° sweep.
        let csv = "\
wavelength_nm,t_0_0_db,t_0_0_deg
1550.0,0.0,179.0
1551.0,0.0,-179.0
";
        let t = SpectrumTable::from_csv(csv, 1).unwrap();
        let (re, im) = t.amp(1550.5e-9, 0, 0);
        assert!(re < -0.999, "should stay near 180°, got re={re} im={im}");
        assert!(im.abs() < 0.02, "im={im}");
    }
}
