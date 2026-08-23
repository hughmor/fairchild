//! Solver tuning knobs.
//!
//! `SimOptions` is the single struct that owns every numerical parameter the
//! solver consumes: tolerances, iteration limits, integration method, source-
//! stepping parameters, etc.  CLI flags, the `.options` directive in a netlist,
//! and Python kwargs all merge into one of these and pass it into the analysis
//! entry points.

use crate::warn_user;
use fairchild_parser::Netlist;

use crate::device::SimContext;
use crate::solver::{make_solver, LinearSolver, SolverKind};
use crate::tran::IntegratorMode;

/// Numerical options consumed by every analysis entry point.
///
/// Defaults match the historic hardcoded SPICE-standard constants.  Override
/// fields individually before passing to `dc_op_nr_with_options` /
/// `tran_nr_with_options` / `ac_analysis_with_options`.
#[derive(Debug, Clone)]
pub struct SimOptions {
    // ── convergence tolerances ─────────────────────────────────────────────
    /// Relative tolerance on Newton update (typical 1e-3).
    pub reltol: f64,
    /// Absolute current tolerance (A); the convergence floor for
    /// voltage-source branch-current rows.
    pub abstol: f64,
    /// Absolute node-voltage tolerance (V); convergence floor for voltage NR.
    pub vntol: f64,
    /// Absolute temperature tolerance (K) for thermal rows.
    ///
    /// A thermal node's potential is kelvin, so `vntol` — a microvolt — is not a
    /// convergence bound on it at all, it is a demand for eight digits nobody
    /// asked for. 1 mK is tight against any self-heating worth modelling and
    /// loose enough that Newton stops.
    pub temptol: f64,
    /// Maximum allowed |Δv| per NR iteration before damping (V).
    pub vmax: f64,
    /// Minimum conductance added to every diagonal entry (S).
    pub gmin: f64,

    // ── iteration limits ────────────────────────────────────────────────────
    /// Max NR iterations per DC-OP solve (ngspice ITL1).
    pub itl1: usize,
    /// Max NR iterations per transient timestep (ngspice ITL4).
    pub itl4: usize,
    /// Max rejected transient steps before bailing (no ngspice equivalent).
    pub max_rejections: usize,

    // ── transient integration ──────────────────────────────────────────────
    /// Integration method.  `BackwardEuler` is unconditionally stable;
    /// `Trapezoidal` is second-order but can ring on discontinuities.
    pub method: IntegratorMode,
    /// Maximum allowed step size (s).  `f64::INFINITY` means "use whatever the
    /// solver decides up to the .tran step argument."
    pub max_step: f64,
    /// Discard transient output before this time — `.tran`'s third argument.
    /// The run still integrates from 0.
    pub tstart: f64,

    // ── convergence aids ───────────────────────────────────────────────────
    /// Initial extra GMIN added during GMIN-stepping (S).
    pub gmin_max: f64,
    /// Number of source-stepping increments to try before giving up.
    pub srcsteps: usize,

    // ── environment ────────────────────────────────────────────────────────
    /// Circuit temperature (K).  Currently informational; future device models
    /// will read it for thermal voltage etc.
    pub temp_k: f64,

    // ── transient initial conditions ───────────────────────────────────────
    /// If true, skip the DC operating point at t=0 and instead seed every
    /// node voltage from `.ic` directives (zero where unspecified).  Equivalent
    /// to the `UIC` keyword on a `.tran` line.
    pub uic: bool,

    // ── Newton-Raphson aids ────────────────────────────────────────────────
    /// Apply junction-step limiters (pnjlim for diodes/BJTs, fetlim for
    /// MOSFETs).  On by default — matches ngspice/hspice.  Disable via
    /// `.options nopnjlim` (or CLI `--no-pnjlim`).
    pub pnjlim: bool,

    // ── linear-system backend ──────────────────────────────────────────────
    /// LU factorisation backend.  `Auto` picks dense for small systems, then
    /// KLU when the `klu` feature is compiled in, `faer-sparse` otherwise
    /// (see `solver::make_solver`).  Override with
    /// `.options solver=sparse|dense|klu|auto` or CLI `--solver`.
    pub solver: SolverKind,

    /// Centre wavelength of the photonic band of interest (m).  Used by
    /// photonic devices as the bootstrap λ for the initial NR iterate, before
    /// a laser has driven the actual λ wire.  Override with `.options
    /// lambda_center_nm=1310` (O-band) or `lambda_center_m=…`.  Default
    /// 1.55 µm (C-band).
    pub lambda_center_m: f64,

    /// Run the netlist sanity-check preflight pass before analysis begins.
    /// On by default — emits warnings to stderr for obvious-but-fatal
    /// netlist errors (R=0, duplicate refdes, fc_* zero-param, etc.).
    /// Silence with `.options nosanitycheck=1` or `sanity_check=0`.
    pub sanity_check: bool,

    /// Emit diagnostic notes about solver progress to stderr.  Off by default.
    /// When on, the analysis entry points (`dc_op_*`, `tran_*`) print:
    ///   - MNA matrix size, NNZ, sparsity, diagonal magnitude spread (once,
    ///     before NR begins);
    ///   - which convergence phase ran (direct NR, source-stepping, gmin-
    ///     stepping) and which one ultimately succeeded;
    ///   - on NR non-convergence: the top-5 rows of the residual vector with
    ///     node / source-branch names, and the device contributing the
    ///     largest residual to each.
    ///     Set via `.options verbose=1`, CLI `--verbose`, or pyo3 `verbose=True`.
    pub verbose: bool,

    /// Enable bidirectional optical propagation.  When off (default),
    /// optical bundles carry (re, im, λ) per channel — light only flows
    /// forward in the direction the device's port topology implies.  When
    /// on, bundles carry (re_fw, im_fw, re_bw, im_bw, λ) per channel, and
    /// every photonic device stamps independent forward + backward paths;
    /// circulators, terminators, and reflective devices become meaningful.
    /// Set via `.options enable_bidirectional=1` or `--opt
    /// enable_bidirectional=1`.
    pub bidirectional_propagation: bool,

    /// Use the LTE-controlled variable-step transient solver instead of the
    /// default fixed-step solver.  The variable-step solver adapts the
    /// internal timestep to keep local truncation error below the
    /// `reltol`/`vntol` budget; `step` becomes the initial (and maximum)
    /// timestep rather than the fixed stride.  Set via `.options
    /// variable_step=1`, CLI `--variable-step`, or Python `variable_step=True`.
    pub variable_step: bool,

    /// Inject device and resistor noise as random currents during `.tran`,
    /// turning the PSDs `.noise` reports into a time-domain waveform. Off by
    /// default: a transient is expected to be reproducible and deterministic,
    /// and every golden in the tree depends on it being so.
    ///
    /// Fixed step only — see `crate::noise::TransientNoise`. Set via `.options
    /// trannoise=1`, or Python `trannoise=True`.
    pub trannoise: bool,

    /// Seed for the transient-noise generator. The same seed gives the same
    /// waveform, so a noisy run is still a reproducible one; sweep it to get
    /// independent trials for a BER or Monte-Carlo estimate.
    pub noiseseed: u64,

    /// Multiplier on every injected noise AMPLITUDE (not power). `2.0` gives
    /// 4× the noise power everywhere, which is the usual trick for pulling a
    /// deep-BER eye closure into a simulation short enough to run.
    pub noisescale: f64,

    /// Model the group delay of optical waveguides (and any device exposing a
    /// group delay τ_g) as a true delay line: the output optical envelope is
    /// the input envelope delayed by τ_g = L·n_g/c, reconstructed from a
    /// per-channel history buffer.  When off (default) the waveguide applies an
    /// instantaneous transmission `exp(-αL/2)·exp(-jβL)` — correct for DC and
    /// steady-state spectra, and cheaper, but it ignores the finite transit
    /// time that matters at modulation bandwidths comparable to 1/τ_g.  For an
    /// electrical transmission line the delay is intrinsic and always modelled.
    /// This flag governs every `OpticalSegment`-based device — the waveguide AND
    /// the active phase shifters / modulators (which carry the same τ_g = L·n_g/c
    /// over their length); a zero-length segment (e.g. `fc_thermal_ps`) stays
    /// instantaneous regardless.  Set via `.options waveguide_delay=1`
    /// (aliases `optical_delay`, `wg_delay`) or `--opt waveguide_delay=1`.
    pub waveguide_delay: bool,

    /// Estimate the 2-norm condition number κ(A) of the MNA matrix at the start
    /// of the DC operating point and print it as a diagnostic.  Useful for
    /// quantifying ill-conditioning (κ ≫ 1e8 is a red flag) before deciding
    /// whether to enable equilibration.  Costs one extra factorisation + a few
    /// power-iteration solves, so it is opt-in.  Set via `.options
    /// cond_estimate=1` or `--opt cond_estimate=1`.
    pub cond_estimate: bool,

    /// Apply two-sided (row/column) equilibration to the MNA matrix before LU
    /// factorisation and unscale the solution afterwards — a readability-neutral
    /// way to improve numerical conditioning of badly-scaled systems (it is
    /// transparent to device code).  Applies to the forward solve only (DC /
    /// transient); the adjoint/transpose path used by `.noise` is left
    /// unscaled.  Set via `.options equilibrate=1` or `--opt equilibrate=1`.
    pub equilibrate: bool,
}

impl Default for SimOptions {
    fn default() -> Self {
        SimOptions {
            reltol: 1e-3,
            abstol: 1e-12,
            vntol: 1e-6,
            temptol: 1e-3,
            vmax: 0.5,
            gmin: 1e-12,
            itl1: 150,
            itl4: 150,
            max_rejections: 30,
            method: IntegratorMode::Trapezoidal,
            max_step: f64::INFINITY,
            tstart: 0.0,
            gmin_max: 1.0,
            srcsteps: 10,
            temp_k: 300.15,
            uic: false,
            pnjlim: true,
            solver: SolverKind::Auto,
            lambda_center_m: 1.55e-6,
            bidirectional_propagation: false,
            verbose: false,
            sanity_check: true,
            variable_step: false,
            trannoise: false,
            noiseseed: 1,
            noisescale: 1.0,
            waveguide_delay: false,
            cond_estimate: false,
            equilibrate: false,
        }
    }
}

impl SimOptions {
    /// Construct options by starting from defaults and folding in every
    /// `.options KEY=VAL` token parsed from the netlist.
    ///
    /// Unrecognised keys are silently dropped (the caller may collect them by
    /// calling `set` directly).  Order of keys in the netlist is preserved, so
    /// later `.options` lines override earlier ones — same semantics as ngspice.
    pub fn from_netlist(netlist: &Netlist) -> Self {
        let mut opts = SimOptions::default();
        // `.temp` precedes `.options temp=…` only as a default; `.options
        // temp=…` written after `.temp` still wins, matching ngspice's
        // last-token-wins rule.
        if let Some(&t_k) = netlist.temps.first() {
            opts.temp_k = t_k;
        }
        for (k, v) in &netlist.options {
            // `set` already reports whether it recognised the key; discarding
            // that made `.options trtol=7` — a real ngspice option — a silent
            // no-op, indistinguishable from one that took effect.
            if !opts.set(k, v) {
                warn_user!(
                    ".options '{k}' is not recognised and has no effect \
                     (see docs/spice_support.md)"
                );
            }
        }
        opts
    }

    /// Apply one `.tran` card's option-bearing fields — `tstart`, `tmax` and
    /// `UIC` — as a unit.
    ///
    /// Call this only from the code about to run *that* card. It used to live in
    /// [`SimOptions::from_netlist`], which every frontend calls whether or not
    /// it is running the deck's analyses: a caller supplying its own `step` and
    /// `stop` still inherited the card's `tmax`, so half a card applied and
    /// nobody chose that. A deck with two `.tran` cards was worse — both runs
    /// got the tightest `tmax` of the two. Either the card is taken whole or it
    /// is not taken.
    pub fn apply_tran_card(&mut self, tstart: f64, tmax: Option<f64>, uic: bool) {
        if tstart > 0.0 {
            self.tstart = tstart;
        }
        if let Some(tmax) = tmax {
            // Tightest constraint wins, so this can never loosen a step ceiling
            // the user set through `.options maxstep` or a CLI flag.
            self.max_step = self.max_step.min(tmax);
        }
        if uic {
            self.uic = true;
        }
    }

    /// Build a `SimContext` consistent with these options.  Threads temperature
    /// and the junction-limiter flag into every device-eval callback.
    pub fn sim_context(&self) -> SimContext {
        SimContext {
            temperature: self.temp_k,
            omega_0: 0.0,
            jlim_enabled: self.pnjlim,
            lambda_center_m: self.lambda_center_m,
            bidirectional_propagation: self.bidirectional_propagation,
            waveguide_delay: self.waveguide_delay,
            time_s: 0.0,
            // Set per step by the transient loops; meaningless in DC/AC.
            discretisation: None,
        }
    }

    /// Build the linear solver matching this options' `solver` choice, sized
    /// for an `n`-row system.  Used by every analysis entry point.  When
    /// `equilibrate` is set, the chosen backend is wrapped in an
    /// [`EquilibratedSolver`] that two-sided-scales the system before LU.
    pub fn linear_solver(&self, n: usize) -> Box<dyn LinearSolver> {
        let base = make_solver(self.solver, n);
        if self.equilibrate {
            Box::new(crate::solver::EquilibratedSolver::new(base))
        } else {
            base
        }
    }

    /// Apply a single `.options KEY=VALUE` token, returning `true` if recognised.
    ///
    /// Used by the parser's `.options` directive and by the CLI's `--opt KEY=VAL`
    /// flag.  Unrecognised keys return `false` so the caller can warn.
    pub fn set(&mut self, key: &str, value: &str) -> bool {
        let key_lc = key.to_lowercase();
        match key_lc.as_str() {
            "reltol" => self.reltol = parse_num(value).unwrap_or(self.reltol),
            "abstol" => self.abstol = parse_num(value).unwrap_or(self.abstol),
            "vntol" => self.vntol = parse_num(value).unwrap_or(self.vntol),
            "temptol" => self.temptol = parse_num(value).unwrap_or(self.temptol),
            "lambdatol" => {
                // Retired rather than dropped: `.options lambdatol=…` in an old
                // deck must say why it no longer applies instead of falling
                // through to "not recognised", which reads like a typo. λ is
                // resolved before the solve, so there is no λ row left to
                // converge and no tolerance for one.
                warn_user!(
                    ".options 'lambdatol' no longer applies: a wavelength is \
                     resolved before the solve rather than solved for, so no \
                     matrix row carries one (see docs/photonic-models.md)"
                );
            }
            "vmax" => self.vmax = parse_num(value).unwrap_or(self.vmax),
            "gmin" => self.gmin = parse_num(value).unwrap_or(self.gmin),
            "itl1" => self.itl1 = parse_int(value).unwrap_or(self.itl1),
            "itl4" => self.itl4 = parse_int(value).unwrap_or(self.itl4),
            "maxstep" | "max_step" => self.max_step = parse_num(value).unwrap_or(self.max_step),
            "tstart" => self.tstart = parse_num(value).unwrap_or(self.tstart),
            "gmin_max" | "gminmax" => self.gmin_max = parse_num(value).unwrap_or(self.gmin_max),
            "srcsteps" | "srcmax" => self.srcsteps = parse_int(value).unwrap_or(self.srcsteps),
            "temp" => self.temp_k = parse_num(value).unwrap_or(self.temp_k) + 273.15,
            "tnom" => self.temp_k = parse_num(value).unwrap_or(self.temp_k) + 273.15,
            "lambda_center_nm" => {
                self.lambda_center_m = parse_num(value)
                    .map(|nm| nm * 1e-9)
                    .unwrap_or(self.lambda_center_m);
            }
            "lambda_center_m" => {
                self.lambda_center_m = parse_num(value).unwrap_or(self.lambda_center_m);
            }
            "verbose" => {
                self.verbose = matches!(
                    value.to_lowercase().as_str(),
                    "" | "1" | "true" | "yes" | "on"
                );
            }
            "sanity_check" | "sanitycheck" => {
                self.sanity_check = matches!(
                    value.to_lowercase().as_str(),
                    "" | "1" | "true" | "yes" | "on"
                );
            }
            "nosanitycheck" | "no_sanity_check" => {
                // ngspice-style bare keyword: `.options nosanitycheck` disables.
                let off = matches!(
                    value.to_lowercase().as_str(),
                    "" | "1" | "true" | "yes" | "on"
                );
                self.sanity_check = !off;
            }
            "enable_bidirectional" | "bidirectional" | "bidirectional_propagation" => {
                self.bidirectional_propagation = matches!(
                    value.to_lowercase().as_str(),
                    "" | "1" | "true" | "yes" | "on"
                );
            }
            "variable_step" | "variablestep" => {
                self.variable_step = matches!(
                    value.to_lowercase().as_str(),
                    "" | "1" | "true" | "yes" | "on"
                );
            }
            "trannoise" | "tran_noise" | "transient_noise" => {
                self.trannoise = matches!(
                    value.to_lowercase().as_str(),
                    "" | "1" | "true" | "yes" | "on"
                );
            }
            "noiseseed" | "noise_seed" => {
                self.noiseseed = parse_num(value).map_or(self.noiseseed, |v| v.max(0.0) as u64);
            }
            "noisescale" | "noise_scale" => {
                self.noisescale = parse_num(value).map_or(self.noisescale, |v| v.max(0.0));
            }
            "waveguide_delay" | "wg_delay" | "optical_delay" => {
                self.waveguide_delay = matches!(
                    value.to_lowercase().as_str(),
                    "" | "1" | "true" | "yes" | "on"
                );
            }
            "cond_estimate" | "estimate_condition_number" | "condest" => {
                self.cond_estimate = matches!(
                    value.to_lowercase().as_str(),
                    "" | "1" | "true" | "yes" | "on"
                );
            }
            "equilibrate" | "equilibration" | "scale" => {
                self.equilibrate = matches!(
                    value.to_lowercase().as_str(),
                    "" | "1" | "true" | "yes" | "on"
                );
            }
            "max_rejections" => {
                self.max_rejections = parse_int(value).unwrap_or(self.max_rejections)
            }
            "method" => match value.to_lowercase().as_str() {
                "be" | "backwardeuler" | "gear1" => self.method = IntegratorMode::BackwardEuler,
                "tr" | "trap" | "trapezoidal" => self.method = IntegratorMode::Trapezoidal,
                "gear" | "gear2" | "bdf2" => self.method = IntegratorMode::Gear,
                _ => return false,
            },
            "uic" => {
                self.uic = matches!(value.to_lowercase().as_str(), "1" | "true" | "yes" | "on");
            }
            "pnjlim" => {
                self.pnjlim = matches!(
                    value.to_lowercase().as_str(),
                    "" | "1" | "true" | "yes" | "on"
                );
            }
            "nopnjlim" => {
                // ngspice-style bare keyword: `.options nopnjlim` disables limiting.
                // Bare token comes in with empty value, treated as "on".
                let off = matches!(
                    value.to_lowercase().as_str(),
                    "" | "1" | "true" | "yes" | "on"
                );
                self.pnjlim = !off;
            }
            "solver" => {
                self.solver = match value.to_lowercase().as_str() {
                    "dense" => SolverKind::Dense,
                    "sparse" | "faer-sparse" | "faer_sparse" => SolverKind::Sparse,
                    "klu" | "suitesparse" | "suitesparse-klu" => {
                        // Reject at options-parse time when KLU feature is absent so
                        // callers (Python, CLI) can surface a clear error rather than
                        // silently falling back to faer-sparse.
                        if cfg!(not(feature = "klu")) {
                            return false;
                        }
                        SolverKind::Klu
                    }
                    "auto" => SolverKind::Auto,
                    _ => return false,
                };
            }
            _ => return false,
        }
        true
    }
}

/// Parse a SPICE-style number with optional suffix (k, meg, m, u, n, p, f).
fn parse_num(s: &str) -> Option<f64> {
    let s_lc = s.to_lowercase();
    let (num, mult) = if let Some(n) = s_lc.strip_suffix("meg") {
        (n, 1e6)
    } else if let Some(n) = s_lc.strip_suffix('k') {
        (n, 1e3)
    } else if let Some(n) = s_lc.strip_suffix('m') {
        (n, 1e-3)
    } else if let Some(n) = s_lc.strip_suffix('u') {
        (n, 1e-6)
    } else if let Some(n) = s_lc.strip_suffix('n') {
        (n, 1e-9)
    } else if let Some(n) = s_lc.strip_suffix('p') {
        (n, 1e-12)
    } else if let Some(n) = s_lc.strip_suffix('f') {
        (n, 1e-15)
    } else {
        (s_lc.as_str(), 1.0)
    };
    num.parse::<f64>().ok().map(|v| v * mult)
}

fn parse_int(s: &str) -> Option<usize> {
    s.parse::<usize>()
        .ok()
        .or_else(|| parse_num(s).map(|v| v as usize))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `solver=klu` without the `klu` feature is **refused**, not quietly served
    /// by faer-sparse — a caller who asked for a specific factorisation and got
    /// another one has no way to tell.
    ///
    /// This behaviour exists only in the default build, never the
    /// `--features klu` one. (Since #76 the two builds also differ in what
    /// `Auto` dispatches to — KLU when compiled in — so the default-features
    /// CI run covers the faer-sparse paths as well as this refusal.)
    #[test]
    #[cfg(not(feature = "klu"))]
    fn solver_klu_is_refused_when_it_is_not_compiled_in() {
        let mut o = SimOptions::default();
        assert!(
            !o.set("solver", "klu"),
            "`solver=klu` must be refused without the feature, not fall back"
        );
        assert!(!o.set("solver", "suitesparse"));
        assert_eq!(o.solver, SolverKind::Auto, "the option must not have moved");
        // The other names still work, so the refusal is about KLU and not about
        // `solver` having stopped parsing.
        assert!(o.set("solver", "sparse"));
        assert_eq!(o.solver, SolverKind::Sparse);
    }

    /// The other half, so the pair cannot both be satisfied by refusing
    /// everything: with the feature, the same string is accepted.
    #[test]
    #[cfg(feature = "klu")]
    fn solver_klu_is_accepted_when_it_is_compiled_in() {
        let mut o = SimOptions::default();
        assert!(o.set("solver", "klu"));
        assert_eq!(o.solver, SolverKind::Klu);
    }

    #[test]
    fn defaults_match_legacy_constants() {
        let o = SimOptions::default();
        assert_eq!(o.reltol, 1e-3);
        assert_eq!(o.vntol, 1e-6);
        assert_eq!(o.vmax, 0.5);
        assert_eq!(o.gmin, 1e-12);
        assert_eq!(o.itl1, 150);
        assert!(matches!(o.method, IntegratorMode::Trapezoidal));
    }

    #[test]
    fn set_recognised_key() {
        let mut o = SimOptions::default();
        assert!(o.set("reltol", "1e-5"));
        assert_eq!(o.reltol, 1e-5);
        assert!(o.set("gmin", "1p"));
        assert!((o.gmin - 1e-12).abs() < 1e-18);
        assert!(o.set("itl1", "300"));
        assert_eq!(o.itl1, 300);
    }

    #[test]
    fn set_method() {
        let mut o = SimOptions::default();
        assert!(o.set("method", "be"));
        assert!(matches!(o.method, IntegratorMode::BackwardEuler));
        assert!(o.set("method", "tr"));
        assert!(matches!(o.method, IntegratorMode::Trapezoidal));
    }

    #[test]
    fn unknown_key_returns_false() {
        let mut o = SimOptions::default();
        assert!(!o.set("not_a_key", "1"));
    }

    #[test]
    fn temp_converts_celsius_to_kelvin() {
        let mut o = SimOptions::default();
        assert!(o.set("temp", "27"));
        assert!((o.temp_k - 300.15).abs() < 1e-6);
    }

    #[test]
    fn uic_flag() {
        let mut o = SimOptions::default();
        assert!(o.set("uic", "1"));
        assert!(o.uic);
        assert!(o.set("uic", "0"));
        assert!(!o.uic);
    }

    #[test]
    fn temp_directive_seeds_temp_k() {
        // .temp 75 should set temp_k to 75 + 273.15 K, threaded through the
        // SimContext so device evals see the right thermal voltage.
        let net =
            fairchild_parser::parse_spice("* tempd\nV1 in 0 DC 0.7\nR1 in out 1k\n.temp 75\n.op\n")
                .unwrap();
        assert_eq!(net.temps.len(), 1);
        assert!((net.temps[0] - 348.15).abs() < 1e-9);
        let opts = SimOptions::from_netlist(&net);
        assert!((opts.temp_k - 348.15).abs() < 1e-9);
        assert!((opts.sim_context().temperature - 348.15).abs() < 1e-9);
    }

    #[test]
    fn temp_list_parsed_for_sweep() {
        let net =
            fairchild_parser::parse_spice("* tempd\nV1 in 0 DC 1\n.temp -40 27 85 125\n.op\n")
                .unwrap();
        assert_eq!(net.temps.len(), 4);
        assert!((net.temps[0] - 233.15).abs() < 1e-9);
        assert!((net.temps[3] - 398.15).abs() < 1e-9);
    }

    #[test]
    fn pnjlim_flag_default_on() {
        let o = SimOptions::default();
        assert!(
            o.pnjlim,
            "pnjlim must default to on (matches ngspice/hspice)"
        );
        let ctx = o.sim_context();
        assert!(ctx.jlim_enabled);
    }

    #[test]
    fn nopnjlim_bare_keyword_disables() {
        let mut o = SimOptions::default();
        assert!(o.set("nopnjlim", "")); // bare token
        assert!(!o.pnjlim);
        // also via pnjlim=0
        let mut o2 = SimOptions::default();
        assert!(o2.set("pnjlim", "0"));
        assert!(!o2.pnjlim);
    }

    #[test]
    fn from_netlist_picks_up_options_directive() {
        let net = fairchild_parser::parse_spice(
            "* opts test\nV1 in 0 DC 1\nR1 in out 1k\n\
             .options reltol=1e-5 gmin=1p method=be itl1=300\n.op\n",
        )
        .unwrap();
        let o = SimOptions::from_netlist(&net);
        assert!((o.reltol - 1e-5).abs() < 1e-18);
        assert!((o.gmin - 1e-12).abs() < 1e-18);
        assert_eq!(o.itl1, 300);
        assert!(matches!(o.method, IntegratorMode::BackwardEuler));
    }

    #[test]
    fn from_netlist_later_overrides_earlier() {
        let net = fairchild_parser::parse_spice(
            "* opts test\nV1 in 0 DC 1\nR1 in out 1k\n\
             .options reltol=1e-3\n.options reltol=1e-7\n.op\n",
        )
        .unwrap();
        let o = SimOptions::from_netlist(&net);
        assert!(
            (o.reltol - 1e-7).abs() < 1e-18,
            "expected last value, got {}",
            o.reltol
        );
    }

    // ── `.tran` cards are applied per run, not folded into every options set ──

    #[test]
    fn from_netlist_leaves_tran_card_alone() {
        // A caller that supplies its own step and stop must get none of the
        // card. `from_netlist` used to fold tstart and tmax in, so the deck's
        // 0.1 ps ceiling silently clamped a run timed entirely from Python.
        let net = fairchild_parser::parse_spice(
            "* tran card
V1 in 0 PULSE(0 1 0 1n 1n 1u 2u)
R1 in out 1k
             C1 out 0 1n
.tran 1p 5n 2n 0.1p
",
        )
        .unwrap();
        let o = SimOptions::from_netlist(&net);
        assert!(
            o.max_step.is_infinite(),
            "card tmax reached options that asked for no card: {}",
            o.max_step
        );
        assert_eq!(o.tstart, 0.0, "card tstart reached options unasked");
    }

    #[test]
    fn apply_tran_card_takes_the_whole_card() {
        let net = fairchild_parser::parse_spice(
            "* tran card
V1 in 0 DC 1
R1 in out 1k
C1 out 0 1n
             .tran 1p 5n 2n 0.1p UIC
",
        )
        .unwrap();
        let fairchild_parser::Analysis::Tran {
            tstart, tmax, uic, ..
        } = net.analyses[0]
        else {
            panic!("expected a .tran card, got {:?}", net.analyses);
        };
        let mut o = SimOptions::from_netlist(&net);
        o.apply_tran_card(tstart, tmax, uic);
        assert_eq!(o.tstart, 2e-9);
        assert!((o.max_step - 1e-13).abs() < 1e-20);
        assert!(o.uic, "UIC on the card must reach the run");
    }

    #[test]
    fn each_tran_card_applies_only_to_its_own_run() {
        // Two cards in one deck: the coarse run must not inherit the fine run's
        // step ceiling. Folding both into one options set gave every run
        // min(tmax) over the whole deck.
        let net = fairchild_parser::parse_spice(
            "* two trans
V1 in 0 DC 1
R1 in out 1k
C1 out 0 1n
             .tran 1n 100n 0 2n
.tran 1p 5n 0 0.1p
",
        )
        .unwrap();
        let base = SimOptions::from_netlist(&net);
        let mut per_card = Vec::new();
        for a in &net.analyses {
            if let fairchild_parser::Analysis::Tran {
                tstart, tmax, uic, ..
            } = *a
            {
                let mut o = base.clone();
                o.apply_tran_card(tstart, tmax, uic);
                per_card.push(o.max_step);
            }
        }
        assert_eq!(per_card.len(), 2);
        assert!((per_card[0] - 2e-9).abs() < 1e-18, "got {}", per_card[0]);
        assert!((per_card[1] - 1e-13).abs() < 1e-20, "got {}", per_card[1]);
    }

    #[test]
    fn options_maxstep_survives_a_looser_card() {
        // Tightest constraint wins: a card may lower the ceiling, never raise it.
        let net = fairchild_parser::parse_spice(
            "* tighter option
V1 in 0 DC 1
R1 in out 1k
             .options maxstep=1e-13
.tran 1p 5n 0 1p
",
        )
        .unwrap();
        let mut o = SimOptions::from_netlist(&net);
        o.apply_tran_card(0.0, Some(1e-12), false);
        assert!((o.max_step - 1e-13).abs() < 1e-20, "got {}", o.max_step);
    }
}
