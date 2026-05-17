//! Solver tuning knobs.
//!
//! `SimOptions` is the single struct that owns every numerical parameter the
//! solver consumes: tolerances, iteration limits, integration method, source-
//! stepping parameters, etc.  CLI flags, the `.options` directive in a netlist,
//! and Python kwargs all merge into one of these and pass it into the analysis
//! entry points.

use fairchild_parser::Netlist;

use crate::device::SimContext;
use crate::solver::{LinearSolver, SolverKind, make_solver};
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
    /// Absolute current tolerance (A); used by current-domain residuals.
    pub abstol: f64,
    /// Absolute node-voltage tolerance (V); convergence floor for voltage NR.
    pub vntol: f64,
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
    /// LU factorisation backend.  `Auto` picks dense for ≤50 nodes and
    /// `faer-sparse` above.  Override with `.options solver=sparse|dense|auto`
    /// or CLI `--solver`.
    pub solver: SolverKind,

    /// Centre wavelength of the photonic band of interest (m).  Used by
    /// photonic devices as the bootstrap λ for the initial NR iterate, before
    /// a laser has driven the actual λ wire.  Override with `.options
    /// lambda_center_nm=1310` (O-band) or `lambda_center_m=…`.  Default
    /// 1.55 µm (C-band).
    pub lambda_center_m: f64,

    /// Emit diagnostic notes about solver progress to stderr.  Off by default.
    /// When on, the analysis entry points (`dc_op_*`, `tran_*`) print:
    ///   - MNA matrix size, NNZ, sparsity, diagonal magnitude spread (once,
    ///     before NR begins);
    ///   - which convergence phase ran (direct NR, source-stepping, gmin-
    ///     stepping) and which one ultimately succeeded;
    ///   - on NR non-convergence: the top-5 rows of the residual vector with
    ///     node / source-branch names, and the device contributing the
    ///     largest residual to each.
    /// Set via `.options verbose=1`, CLI `--verbose`, or pyo3 `verbose=True`.
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
}

impl Default for SimOptions {
    fn default() -> Self {
        SimOptions {
            reltol:         1e-3,
            abstol:         1e-12,
            vntol:          1e-6,
            vmax:           0.5,
            gmin:           1e-12,
            itl1:           150,
            itl4:           150,
            max_rejections: 30,
            method:         IntegratorMode::Trapezoidal,
            max_step:       f64::INFINITY,
            gmin_max:       1.0,
            srcsteps:       10,
            temp_k:         300.15,
            uic:            false,
            pnjlim:         true,
            solver:         SolverKind::Auto,
            lambda_center_m: 1.55e-6,
            bidirectional_propagation: false,
            verbose:        false,
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
            opts.set(k, v);
        }
        opts
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
        }
    }

    /// Build the linear solver matching this options' `solver` choice, sized
    /// for an `n`-row system.  Used by every analysis entry point.
    pub fn linear_solver(&self, n: usize) -> Box<dyn LinearSolver> {
        make_solver(self.solver, n)
    }

    /// Apply a single `.options KEY=VALUE` token, returning `true` if recognised.
    ///
    /// Used by the parser's `.options` directive and by the CLI's `--opt KEY=VAL`
    /// flag.  Unrecognised keys return `false` so the caller can warn.
    pub fn set(&mut self, key: &str, value: &str) -> bool {
        let key_lc = key.to_lowercase();
        match key_lc.as_str() {
            "reltol"  => self.reltol  = parse_num(value).unwrap_or(self.reltol),
            "abstol"  => self.abstol  = parse_num(value).unwrap_or(self.abstol),
            "vntol"   => self.vntol   = parse_num(value).unwrap_or(self.vntol),
            "vmax"    => self.vmax    = parse_num(value).unwrap_or(self.vmax),
            "gmin"    => self.gmin    = parse_num(value).unwrap_or(self.gmin),
            "itl1"    => self.itl1    = parse_int(value).unwrap_or(self.itl1),
            "itl4"    => self.itl4    = parse_int(value).unwrap_or(self.itl4),
            "maxstep" | "max_step" => self.max_step = parse_num(value).unwrap_or(self.max_step),
            "gmin_max" | "gminmax" => self.gmin_max = parse_num(value).unwrap_or(self.gmin_max),
            "srcsteps" | "srcmax"  => self.srcsteps = parse_int(value).unwrap_or(self.srcsteps),
            "temp"    => self.temp_k = parse_num(value).unwrap_or(self.temp_k) + 273.15,
            "tnom"    => self.temp_k = parse_num(value).unwrap_or(self.temp_k) + 273.15,
            "lambda_center_nm" => {
                self.lambda_center_m = parse_num(value).map(|nm| nm * 1e-9)
                    .unwrap_or(self.lambda_center_m);
            }
            "lambda_center_m"  => {
                self.lambda_center_m = parse_num(value).unwrap_or(self.lambda_center_m);
            }
            "verbose" => {
                self.verbose = matches!(value.to_lowercase().as_str(),
                    "" | "1" | "true" | "yes" | "on");
            }
            "enable_bidirectional" | "bidirectional" | "bidirectional_propagation" => {
                self.bidirectional_propagation = matches!(value.to_lowercase().as_str(),
                    "" | "1" | "true" | "yes" | "on");
            }
            "max_rejections" => self.max_rejections = parse_int(value).unwrap_or(self.max_rejections),
            "method" => {
                match value.to_lowercase().as_str() {
                    "be" | "backwardeuler" | "gear1" => self.method = IntegratorMode::BackwardEuler,
                    "tr" | "trap" | "trapezoidal"    => self.method = IntegratorMode::Trapezoidal,
                    "gear" | "gear2" | "bdf2"        => self.method = IntegratorMode::Gear,
                    _ => return false,
                }
            }
            "uic" => {
                self.uic = matches!(value.to_lowercase().as_str(),
                    "1" | "true" | "yes" | "on");
            }
            "pnjlim" => {
                self.pnjlim = matches!(value.to_lowercase().as_str(),
                    "" | "1" | "true" | "yes" | "on");
            }
            "nopnjlim" => {
                // ngspice-style bare keyword: `.options nopnjlim` disables limiting.
                // Bare token comes in with empty value, treated as "on".
                let off = matches!(value.to_lowercase().as_str(),
                    "" | "1" | "true" | "yes" | "on");
                self.pnjlim = !off;
            }
            "solver" => {
                self.solver = match value.to_lowercase().as_str() {
                    "dense"  => SolverKind::Dense,
                    "sparse" | "faer-sparse" | "faer_sparse" => SolverKind::Sparse,
                    "auto"   => SolverKind::Auto,
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
    let (num, mult) = if let Some(n) = s_lc.strip_suffix("meg") { (n, 1e6) }
        else if let Some(n) = s_lc.strip_suffix('k') { (n, 1e3) }
        else if let Some(n) = s_lc.strip_suffix('m') { (n, 1e-3) }
        else if let Some(n) = s_lc.strip_suffix('u') { (n, 1e-6) }
        else if let Some(n) = s_lc.strip_suffix('n') { (n, 1e-9) }
        else if let Some(n) = s_lc.strip_suffix('p') { (n, 1e-12) }
        else if let Some(n) = s_lc.strip_suffix('f') { (n, 1e-15) }
        else { (s_lc.as_str(), 1.0) };
    num.parse::<f64>().ok().map(|v| v * mult)
}

fn parse_int(s: &str) -> Option<usize> {
    s.parse::<usize>().ok().or_else(|| parse_num(s).map(|v| v as usize))
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let net = fairchild_parser::parse_spice(
            "* tempd\nV1 in 0 DC 0.7\nR1 in out 1k\n.temp 75\n.op\n.end\n"
        ).unwrap();
        assert_eq!(net.temps.len(), 1);
        assert!((net.temps[0] - 348.15).abs() < 1e-9);
        let opts = SimOptions::from_netlist(&net);
        assert!((opts.temp_k - 348.15).abs() < 1e-9);
        assert!((opts.sim_context().temperature - 348.15).abs() < 1e-9);
    }

    #[test]
    fn temp_list_parsed_for_sweep() {
        let net = fairchild_parser::parse_spice(
            "* tempd\nV1 in 0 DC 1\n.temp -40 27 85 125\n.op\n.end\n"
        ).unwrap();
        assert_eq!(net.temps.len(), 4);
        assert!((net.temps[0] - 233.15).abs() < 1e-9);
        assert!((net.temps[3] - 398.15).abs() < 1e-9);
    }

    #[test]
    fn pnjlim_flag_default_on() {
        let o = SimOptions::default();
        assert!(o.pnjlim, "pnjlim must default to on (matches ngspice/hspice)");
        let ctx = o.sim_context();
        assert!(ctx.jlim_enabled);
    }

    #[test]
    fn nopnjlim_bare_keyword_disables() {
        let mut o = SimOptions::default();
        assert!(o.set("nopnjlim", ""));   // bare token
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
             .options reltol=1e-5 gmin=1p method=be itl1=300\n.op\n.end\n"
        ).unwrap();
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
             .options reltol=1e-3\n.options reltol=1e-7\n.op\n.end\n"
        ).unwrap();
        let o = SimOptions::from_netlist(&net);
        assert!((o.reltol - 1e-7).abs() < 1e-18, "expected last value, got {}", o.reltol);
    }
}
