use std::fs;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use clap::{Parser, ValueEnum};
use rayon::prelude::*;

use fairchild_core::{
    ac_analysis_opts, dc_op_nr_with_registry_opts, dc_sweep_with_registry_opts,
    evaluate_measurements,
    freq_decade, freq_linear, freq_oct,
    tran_nr_with_registry_opts, tran_nr_with_registry_var_opts, DeviceRegistry, SimOptions,
};
use fairchild_osdi::OsdiLibrary;
use fairchild_parser::{check_disciplines, parse_spice_file, AcVariation, Analysis, Element, Netlist};

#[derive(Parser)]
#[command(
    name = "fairchild",
    version,
    about = "Open-source time-domain electro-optic circuit simulator",
    long_about = "Fairchild simulates analog circuits containing both electronic and photonic \
                  components in the same Newton-Raphson loop.  Supports DC, transient, and \
                  small-signal AC analyses; loads Verilog-A models compiled with OpenVAF \
                  via the OSDI v0.4 interface."
)]
struct Cli {
    /// Input SPICE netlist file
    #[arg(short, long)]
    file: PathBuf,

    /// Output format
    #[arg(long, default_value = "csv")]
    format: Format,

    /// Output file (default: stdout)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Comma-separated list of signals to include in output.
    /// Example: --probe "V(out),V(in),I(V1)"
    /// Applies to CSV output only; nutmeg always outputs all signals.
    #[arg(long, value_name = "SIGNAL,...")]
    probe: Option<String>,

    /// Override a circuit parameter.  Format: ELEMENT.PARAM=VALUE
    /// Example: --param "Xcoupler.kappa_0=0.05" --param "Rload.resistance=2e3"
    /// Can be specified multiple times.
    #[arg(long = "param", value_name = "ELEMENT.PARAM=VALUE")]
    params: Vec<String>,

    /// Parse and discipline-check the netlist, then exit without simulating.
    /// Exit code 0 if valid, 1 if any errors found.
    #[arg(long)]
    check: bool,

    /// List node names parsed from the netlist, then exit.
    #[arg(long)]
    list_nodes: bool,

    /// List model cards (.model statements) parsed from the netlist, then exit.
    #[arg(long)]
    list_models: bool,

    /// Print simulation progress and iteration counts.
    #[arg(long, short)]
    verbose: bool,

    /// Suppress all warning messages.
    #[arg(long, short)]
    quiet: bool,

    // ── solver tuning knobs (overlay onto netlist `.options`) ────────────
    /// Override an arbitrary solver option.  Format: KEY=VALUE.  Layered on
    /// top of any `.options` directives in the netlist.  Can be repeated.
    ///
    /// Recognised keys: reltol, abstol, vntol, gmin, vmax, itl1, itl4,
    /// maxstep, gminmax, srcsteps, method (be|tr|gear), uic, temp,
    /// variable_step.
    ///
    /// Example: --opt reltol=1e-5 --opt method=gear --opt variable_step=1
    #[arg(long = "opt", value_name = "KEY=VALUE")]
    options: Vec<String>,

    /// Convenience flag: relative Newton tolerance.
    #[arg(long, value_name = "VALUE")]
    reltol: Option<String>,

    /// Convenience flag: minimum diagonal conductance (S).
    #[arg(long, value_name = "VALUE")]
    gmin: Option<String>,

    /// Convenience flag: transient integration method.
    ///   be    — Backward Euler (BDF-1, robust)
    ///   tr    — Trapezoidal Rule (2nd-order, can ring)
    ///   gear  — GEAR / BDF-2 (2nd-order, L-stable, no ringing)
    #[arg(long, value_name = "be|tr|gear")]
    method: Option<String>,

    /// Convenience flag: maximum transient step size (s).
    #[arg(long = "maxstep", value_name = "VALUE")]
    max_step: Option<String>,

    /// Disable junction-step limiters (pnjlim for diodes/BJTs, fetlim for
    /// MOSFETs).  Equivalent to `.options nopnjlim`.  Limiters are ON by
    /// default — turning them off occasionally helps diagnose convergence
    /// failures by surfacing the raw NR step.
    #[arg(long = "no-pnjlim")]
    no_pnjlim: bool,

    /// Linear-system backend.
    ///   dense  — faer partial-pivot LU (recommended for ≤ ~50 nodes)
    ///   sparse — faer sparse LU (pure Rust, default at larger N)
    ///   klu    — SuiteSparse KLU (requires `klu` cargo feature + system
    ///            install of suite-sparse; 2–5× faster on circuit matrices)
    ///   auto   — pick from system size (default)
    #[arg(long, value_name = "dense|sparse|klu|auto")]
    solver: Option<String>,

    /// Use the LTE-controlled variable-step transient solver.  `step` in the
    /// netlist becomes the initial / maximum timestep rather than a fixed
    /// stride.  Equivalent to `.options variable_step=1`.
    #[arg(long)]
    variable_step: bool,

    /// Bundle all `.alter` × `.temp` corner outputs into the single
    /// `--output` file (with `# alter=…` / `# temp_c=…` header lines),
    /// preserving the historic concatenated layout.
    ///
    /// Default behaviour when `--output` is given: write one file per
    /// corner (e.g. `out.alter_pvtfast.temp_-40c.csv`) and run the
    /// corners in parallel.  With `--single-output`, runs are serial
    /// so headers and rows stay in deterministic order.
    #[arg(long)]
    single_output: bool,
}

#[derive(Clone, ValueEnum)]
enum Format {
    Csv,
    Nutmeg,
}

// ── Parameter override helpers ─────────────────────────────────────────────

/// Apply CLI `--param` overrides to a netlist in-place.
///
/// Format: "ELEMENT.PARAM=VALUE" (case-insensitive element and param names).
/// Supports XOsdi elements, Resistor, Capacitor, Inductor.
fn apply_params(netlist: &mut Netlist, overrides: &[String], quiet: bool) {
    for raw in overrides {
        let (lhs, rhs) = match raw.split_once('=') {
            Some(pair) => pair,
            None => {
                if !quiet { eprintln!("warning: --param '{raw}': expected ELEMENT.PARAM=VALUE, skipping"); }
                continue;
            }
        };
        let value: f64 = match rhs.parse() {
            Ok(v) => v,
            Err(_) => {
                if !quiet { eprintln!("warning: --param '{raw}': cannot parse value '{rhs}', skipping"); }
                continue;
            }
        };
        let (elem_name, param_name) = match lhs.split_once('.') {
            Some(pair) => pair,
            None => {
                if !quiet { eprintln!("warning: --param '{raw}': expected ELEMENT.PARAM, skipping"); }
                continue;
            }
        };
        let elem_name_lc  = elem_name.to_lowercase();
        let param_name_lc = param_name.to_lowercase();

        let mut applied = false;
        for el in &mut netlist.elements {
            match el {
                Element::XOsdi { name, params, .. } if name.to_lowercase() == elem_name_lc => {
                    if let Some(slot) = params.iter_mut().find(|(k, _)| k.to_lowercase() == param_name_lc) {
                        slot.1 = value;
                    } else {
                        params.push((param_name_lc.clone(), value));
                    }
                    applied = true;
                    break;
                }
                Element::Resistor { name, resistance, .. }
                    if name.to_lowercase() == elem_name_lc
                       && (param_name_lc == "resistance" || param_name_lc == "value" || param_name_lc == "r") =>
                {
                    *resistance = value;
                    applied = true;
                    break;
                }
                Element::Capacitor { name, capacitance, .. }
                    if name.to_lowercase() == elem_name_lc
                       && (param_name_lc == "capacitance" || param_name_lc == "value" || param_name_lc == "c") =>
                {
                    *capacitance = value;
                    applied = true;
                    break;
                }
                Element::Inductor { name, inductance, .. }
                    if name.to_lowercase() == elem_name_lc
                       && (param_name_lc == "inductance" || param_name_lc == "value" || param_name_lc == "l") =>
                {
                    *inductance = value;
                    applied = true;
                    break;
                }
                Element::VoltageSource { name, waveform, .. }
                    if name.to_lowercase() == elem_name_lc
                       && (param_name_lc == "dc" || param_name_lc == "value" || param_name_lc == "v") =>
                {
                    *waveform = fairchild_parser::Waveform::Dc(value);
                    applied = true;
                    break;
                }
                Element::CurrentSource { name, waveform, .. }
                    if name.to_lowercase() == elem_name_lc
                       && (param_name_lc == "dc" || param_name_lc == "value" || param_name_lc == "i") =>
                {
                    *waveform = fairchild_parser::Waveform::Dc(value);
                    applied = true;
                    break;
                }
                _ => {}
            }
        }
        if !applied && !quiet {
            eprintln!("warning: --param '{raw}': element '{elem_name}' not found or param not applicable");
        }
    }
}

// ── Probe / column filtering ───────────────────────────────────────────────

/// Parse a comma-separated probe string into normalised signal names.
fn parse_probe(s: &str) -> Vec<String> {
    s.split(',').map(|t| t.trim().to_lowercase()).filter(|t| !t.is_empty()).collect()
}

/// Filter a CSV string to keep only the header columns that match `probes`.
///
/// Column names are lowercased before matching.  The first column (analysis /
/// time) is always kept.  Returns the full CSV if `probes` is empty.
fn filter_csv(csv: &str, probes: &[String]) -> String {
    if probes.is_empty() {
        return csv.to_string();
    }
    let mut out = String::new();
    let mut keep_cols: Option<Vec<usize>> = None;

    for line in csv.lines() {
        if keep_cols.is_none() {
            // Header row — determine which columns to keep
            let headers: Vec<&str> = line.split(',').collect();
            let cols: Vec<usize> = headers
                .iter()
                .enumerate()
                .filter(|(i, h)| *i == 0 || probes.iter().any(|p| p == &h.to_lowercase()))
                .map(|(i, _)| i)
                .collect();
            // Emit filtered header
            let header_out: Vec<&str> = cols.iter().map(|&i| headers[i]).collect();
            out.push_str(&header_out.join(","));
            out.push('\n');
            keep_cols = Some(cols);
        } else {
            let cols = keep_cols.as_ref().unwrap();
            let fields: Vec<&str> = line.split(',').collect();
            let row: Vec<&str> = cols.iter().filter_map(|&i| fields.get(i).copied()).collect();
            out.push_str(&row.join(","));
            out.push('\n');
        }
    }
    out
}

// ── SimOptions builder: netlist .options + CLI flags ───────────────────────

/// Build a `SimOptions` by starting from netlist `.options` and applying the
/// CLI-specified overlays.  Unknown keys emit a warning (unless `--quiet`).
fn build_options(netlist: &Netlist, cli: &Cli) -> SimOptions {
    let mut opts = SimOptions::from_netlist(netlist);

    let mut apply = |key: &str, value: &str| {
        if !opts.set(key, value) && !cli.quiet {
            eprintln!("warning: unknown solver option '{key}={value}'");
        }
    };

    if let Some(v) = &cli.reltol      { apply("reltol",        v); }
    if let Some(v) = &cli.gmin       { apply("gmin",          v); }
    if let Some(v) = &cli.method     { apply("method",        v); }
    if let Some(v) = &cli.max_step   { apply("maxstep",       v); }
    if cli.no_pnjlim                 { apply("pnjlim",        "0"); }
    if let Some(v) = &cli.solver     { apply("solver",        v); }
    if cli.variable_step             { apply("variable_step", "1"); }

    for raw in &cli.options {
        if let Some((k, v)) = raw.split_once('=') {
            apply(k.trim(), v.trim().trim_matches('"').trim_matches('\''));
        } else if !cli.quiet {
            eprintln!("warning: --opt '{raw}': expected KEY=VALUE, skipping");
        }
    }

    if cli.verbose { opts.verbose = true; }
    opts
}

// ── OSDI registry builder ─────────────────────────────────────────────────

/// Load built-in models + any OSDI shared libraries listed in the netlist.
/// Relative `.osdi` paths are resolved against `netlist_dir`.
fn build_registry(netlist: &Netlist, netlist_dir: Option<&PathBuf>, quiet: bool) -> DeviceRegistry {
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_diodes(&netlist.models);
    registry.register_builtin_mosfets(&netlist.models);
    registry.register_builtin_bjts(&netlist.models);

    // Photonic-model authoring guidance: as of the B-phase refactor, native
    // Rust photonic devices (`fc_waveguide`, `fc_dcoupler`, `fc_splitter`,
    // `fc_photodetector`, `fc_thermal_ps`, `fc_pn_ps`) are the recommended
    // path.  Surface a one-shot info note so users with `.osdi` photonic
    // models know there's a faster, cleaner alternative.
    if !quiet && !netlist.osdi_paths.is_empty() {
        let photonic_count = netlist.osdi_paths.iter()
            .filter(|p| p.contains("photonic")
                || p.contains("waveguide")
                || p.contains("mrr")
                || p.contains("mzi")
                || p.contains("laser"))
            .count();
        if photonic_count > 0 {
            eprintln!(
                "info: {} OSDI photonic library/libraries loaded — note that native Rust devices \
                 (fc_waveguide etc.) are now the recommended path; see fairchild-osdi crate \
                 docs for the deprecation rationale.",
                photonic_count
            );
        }
    }

    for osdi_path in &netlist.osdi_paths {
        let path = if std::path::Path::new(osdi_path).is_absolute() {
            PathBuf::from(osdi_path)
        } else if let Some(dir) = netlist_dir {
            dir.join(osdi_path)
        } else {
            PathBuf::from(osdi_path)
        };

        let lib = unsafe { OsdiLibrary::open(&path) }.unwrap_or_else(|e| {
            eprintln!("error: cannot load OSDI library '{}': {e}", path.display());
            std::process::exit(1);
        });
        let lib = Arc::new(lib);
        lib.register_into(&mut registry);
    }

    registry
}

// ── Main ──────────────────────────────────────────────────────────────────

fn main() {
    let cli = Cli::parse();

    let mut netlist = parse_spice_file(&cli.file).unwrap_or_else(|e| {
        eprintln!("error: {e}");
        std::process::exit(1);
    });

    // Apply --param overrides before any structural checks
    if !cli.params.is_empty() {
        apply_params(&mut netlist, &cli.params, cli.quiet);
    }

    // Discipline check (always performed; emits errors + exits on failure)
    if let Err(e) = check_disciplines(&netlist) {
        eprintln!("error: discipline mismatch: {e}");
        std::process::exit(1);
    }

    // --check: validate only
    if cli.check {
        if !cli.quiet {
            let n_el = netlist.elements.len();
            let n_an = netlist.analyses.len();
            eprintln!("ok: {} element(s), {} analysis/analyses, disciplines clean", n_el, n_an);
        }
        std::process::exit(0);
    }

    // --list-nodes: enumerate nodes and exit
    if cli.list_nodes {
        // Build a registry for the netlist_dir (needed to load OSDI libs for topology)
        let netlist_dir_tmp = cli.file.parent().map(|p| p.to_path_buf());
        let reg_tmp = build_registry(&netlist, netlist_dir_tmp.as_ref(), cli.quiet);
        let opts_tmp = build_options(&netlist, &cli);
        let result = dc_op_nr_with_registry_opts(&netlist, &reg_tmp, &opts_tmp).unwrap_or_else(|e| {
            eprintln!("error: cannot build topology: {e}");
            std::process::exit(1);
        });
        let mut nodes: Vec<&str> = result.topo.node_index.keys().map(|s| s.as_str()).collect();
        nodes.sort_unstable();
        for n in nodes { println!("V({n})"); }
        for n in result.topo.vsrc_index.keys() { println!("I({n})"); }
        std::process::exit(0);
    }

    // --list-models: print model cards
    if cli.list_models {
        if netlist.models.is_empty() {
            println!("(no .model cards found)");
        }
        for m in &netlist.models {
            let params: Vec<String> = m.params.iter().map(|(k, v)| format!("{k}={v}")).collect();
            println!(".model {} {} {}", m.name, m.kind, params.join(" "));
        }
        for path in &netlist.osdi_paths {
            println!(".osdi {path}");
        }
        std::process::exit(0);
    }

    let title = netlist.title.clone();
    let probe_list: Vec<String> = cli.probe.as_deref().map(parse_probe).unwrap_or_default();

    // Build device registry: built-in models + OSDI shared libraries.
    // Constructed inside the .alter loop below so model overrides take effect.
    let netlist_dir = cli.file.parent().map(|p| p.to_path_buf());

    // Merge netlist `.options` + CLI flag overrides into a single SimOptions.
    let opts = build_options(&netlist, &cli);
    if cli.verbose {
        eprintln!("info: solver options: reltol={:e} gmin={:e} method={:?} itl1={} itl4={}",
            opts.reltol, opts.gmin, opts.method, opts.itl1, opts.itl4);
    }

    // .temp <T1> [<T2> ...] sweep: re-run every analysis once per temperature.
    // Empty `temps` ⇒ single pass at whatever `opts.temp_k` already is.
    let temp_sweep: Vec<f64> = if netlist.temps.len() > 1 {
        netlist.temps.clone()
    } else {
        vec![opts.temp_k]
    };
    // .alter sweep: base run + one re-run per .alter block.
    let mut alter_runs: Vec<(String, Netlist)> = vec![("base".into(), netlist.clone())];
    for block in &netlist.alters {
        let mut patched = netlist.clone();
        patched.apply_alter(block);
        alter_runs.push((block.label.clone(), patched));
    }

    // Flatten the (alter × temp) grid into a list of corners.  Each
    // corner is an independent simulation; we either run them serially
    // into one shared writer (`--single-output`, no `--output`, only
    // one corner, or `--verbose`) or in parallel into per-corner files.
    let mut corners: Vec<Corner> = Vec::with_capacity(alter_runs.len() * temp_sweep.len());
    for (ai, (alter_label, alter_netlist)) in alter_runs.iter().enumerate() {
        for (ti, &temp_k) in temp_sweep.iter().enumerate() {
            let mut corner_opts = opts.clone();
            corner_opts.temp_k = temp_k;
            corners.push(Corner {
                alter_idx: ai,
                temp_idx: ti,
                alter_label: alter_label.clone(),
                temp_k,
                netlist: alter_netlist.clone(),
                opts: corner_opts,
            });
        }
    }

    let n_alters = alter_runs.len();
    let n_temps = temp_sweep.len();
    let n_corners = corners.len();

    // Dispatch: file-per-corner only when (a) an output path is given,
    // (b) we have more than one corner, (c) the user hasn't asked for
    // the single-output bundle, and (d) verbose is off (otherwise the
    // interleaved per-corner NR diagnostics would be unreadable).
    let parallel_eligible = cli.output.is_some()
        && n_corners > 1
        && !cli.single_output
        && !cli.verbose;

    let ran_something = if parallel_eligible {
        run_corners_parallel(
            &corners, n_alters, n_temps,
            netlist_dir.as_ref(), &probe_list, &title, &cli,
        )
    } else {
        run_corners_serial(
            &corners, n_alters, n_temps,
            netlist_dir.as_ref(), &probe_list, &title, &cli,
        )
    };

    if !ran_something && !cli.quiet {
        eprintln!("warning: no analyses found in netlist (add .op, .tran, or .ac)");
    }
}

// ---------------------------------------------------------------------------
// Corner sweep — `.alter` × `.temp` grid
// ---------------------------------------------------------------------------

/// One leaf of the `.alter` × `.temp` grid: a fully-resolved netlist
/// plus a `SimOptions` carrying the per-corner temperature.
struct Corner {
    alter_idx:  usize,
    temp_idx:   usize,
    alter_label: String,
    temp_k:     f64,
    netlist:    Netlist,
    opts:       SimOptions,
}

/// Derive a per-corner output path by suffixing the base `--output`
/// path with `.alter_<label>` and `.temp_<C>c` where the corresponding
/// sweep is non-trivial.  Single-corner runs pass through unchanged.
fn corner_path(base: &Path, alter_label: &str, n_alters: usize, temp_k: f64, n_temps: usize) -> PathBuf {
    if n_alters <= 1 && n_temps <= 1 {
        return base.to_path_buf();
    }
    let stem = base.file_stem().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default();
    let ext  = base.extension().map(|e| e.to_string_lossy().into_owned()).unwrap_or_default();
    let parent = base.parent().unwrap_or_else(|| Path::new("."));
    let mut name = stem;
    if n_alters > 1 {
        name.push_str(".alter_");
        name.push_str(&sanitize_label(alter_label));
    }
    if n_temps > 1 {
        let temp_c = temp_k - 273.15;
        // `{:+.0}` keeps the sign so `-40c` and `27c` are immediately
        // distinguishable, while `{:.0}` collapses 26.9°C to "27c".
        name.push_str(&format!(".temp_{:.0}c", temp_c));
    }
    if !ext.is_empty() {
        name.push('.');
        name.push_str(&ext);
    }
    parent.join(name)
}

/// Replace whitespace and path separators in `.alter` labels so they
/// produce valid filename suffixes.
fn sanitize_label(s: &str) -> String {
    s.chars().map(|c| match c {
        ' ' | '\t' | '/' | '\\' | ':' => '_',
        _ => c,
    }).collect()
}

fn open_writer(path: &Path) -> Box<dyn Write + Send> {
    let f = fs::File::create(path).unwrap_or_else(|e| {
        eprintln!("error: cannot create {path:?}: {e}");
        std::process::exit(1);
    });
    Box::new(BufWriter::new(f))
}

/// Run every corner sequentially into a single shared writer.  Used
/// when `--single-output` is given, when there's only one corner,
/// when no `--output` was specified (writes go to stdout), or under
/// `--verbose` (per-corner diagnostics would otherwise interleave).
fn run_corners_serial(
    corners: &[Corner], n_alters: usize, n_temps: usize,
    netlist_dir: Option<&PathBuf>, probe_list: &[String], title: &str, cli: &Cli,
) -> bool {
    let mut w: Box<dyn Write> = match &cli.output {
        Some(path) => open_writer(path),
        None       => Box::new(BufWriter::new(io::stdout())),
    };
    let mut ran_something = false;
    let mut last_alter: Option<usize> = None;
    let mut registry: Option<DeviceRegistry> = None;
    for corner in corners {
        // Rebuild registry on every new alter block (model overrides may
        // differ); reuse across temperature sweep points within a block.
        if last_alter != Some(corner.alter_idx) {
            registry = Some(build_registry(&corner.netlist, netlist_dir, cli.quiet));
            last_alter = Some(corner.alter_idx);
        }
        if n_alters > 1 {
            writeln!(w, "# alter={} (block {}/{})",
                corner.alter_label, corner.alter_idx + 1, n_alters)
                .unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
            if cli.verbose {
                eprintln!("info: .alter block {}/{}: '{}'",
                    corner.alter_idx + 1, n_alters, corner.alter_label);
            }
        }
        if n_temps > 1 {
            writeln!(w, "# temp_c={:.3} (point {}/{})",
                corner.temp_k - 273.15, corner.temp_idx + 1, n_temps)
                .unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
            if cli.verbose {
                eprintln!("info: temperature sweep point {}/{}: {:.2} °C",
                    corner.temp_idx + 1, n_temps, corner.temp_k - 273.15);
            }
        }
        if run_corner_analyses(
            corner, registry.as_ref().unwrap(),
            probe_list, title, cli, &mut w,
        ) {
            ran_something = true;
        }
    }
    ran_something
}

/// Run every corner on its own thread, writing to per-corner files
/// derived from the base `--output` path.  Each thread builds its own
/// device registry so OSDI library handles aren't shared (the OSDI
/// loader uses `dlopen` which is process-global, but each registry
/// owns its own `OsdiLibrary` handles).
fn run_corners_parallel(
    corners: &[Corner], n_alters: usize, n_temps: usize,
    netlist_dir: Option<&PathBuf>, probe_list: &[String], title: &str, cli: &Cli,
) -> bool {
    let base = cli.output.as_ref().expect("--output required for parallel mode");
    let netlist_dir = netlist_dir.cloned();
    let probe_list = probe_list.to_vec();
    let title = title.to_string();
    let cli_quiet = cli.quiet;
    let cli_verbose = cli.verbose;
    let cli_format = cli.format.clone();

    let ran: Vec<bool> = corners.par_iter()
        .map(|corner| {
            let out_path = corner_path(base, &corner.alter_label, n_alters, corner.temp_k, n_temps);
            let mut w: Box<dyn Write> = open_writer(&out_path);
            let registry = build_registry(&corner.netlist, netlist_dir.as_ref(), cli_quiet);
            // Build a synthetic Cli reference for run_corner_analyses; we
            // only need a few fields, so pass them explicitly via a
            // miniature struct rather than threading the whole Cli through.
            let ctx = CornerCtx {
                verbose: cli_verbose,
                format: &cli_format,
            };
            run_corner_analyses_ctx(corner, &registry, &probe_list, &title, &ctx, &mut w)
        })
        .collect();
    ran.into_iter().any(|x| x)
}

/// Minimal subset of CLI flags that `run_corner_analyses_ctx` actually
/// needs.  Used so the parallel path can pass small `Copy`/borrowed
/// fields per-worker instead of an entire `Cli` (which holds owned
/// `String` / `Vec<String>` that would force cloning per corner).
struct CornerCtx<'a> {
    verbose: bool,
    format:  &'a Format,
}

/// Convenience wrapper for the serial path that has a `&Cli` available.
fn run_corner_analyses(
    corner: &Corner, registry: &DeviceRegistry,
    probe_list: &[String], title: &str, cli: &Cli,
    w: &mut dyn Write,
) -> bool {
    let ctx = CornerCtx { verbose: cli.verbose, format: &cli.format };
    run_corner_analyses_ctx(corner, registry, probe_list, title, &ctx, w)
}

/// Run every `Analysis` declared on this corner's netlist, writing
/// results to `w`.  Returns `true` if anything was emitted.
fn run_corner_analyses_ctx(
    corner: &Corner, registry: &DeviceRegistry,
    probe_list: &[String], title: &str, ctx: &CornerCtx,
    w: &mut dyn Write,
) -> bool {
    let netlist = &corner.netlist;
    let opts = &corner.opts;
    let mut ran_something = false;
    for analysis in &netlist.analyses {
        match analysis {
            Analysis::Op => {
                if ctx.verbose { eprintln!("info: running DC operating-point analysis..."); }
                let t0 = Instant::now();
                let result = dc_op_nr_with_registry_opts(netlist, registry, opts).unwrap_or_else(|e| {
                    eprintln!("error: DC op failed: {e}");
                    std::process::exit(1);
                });
                if ctx.verbose {
                    eprintln!("info: DC op converged in {} iteration(s) [{:.1} ms]",
                        result.iters, t0.elapsed().as_secs_f64() * 1000.0);
                }
                match ctx.format {
                    Format::Csv => {
                        let mut buf = Vec::new();
                        result.write_csv(&mut buf).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                        let csv = String::from_utf8_lossy(&buf);
                        let filtered = filter_csv(&csv, probe_list);
                        w.write_all(filtered.as_bytes()).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                    }
                    Format::Nutmeg => {
                        result.write_nutmeg(&mut *w, title).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                    }
                }
                ran_something = true;
            }

            Analysis::Tran { step, stop } => {
                if ctx.verbose {
                    let mode = if opts.variable_step { "variable-step" } else { "fixed-step" };
                    eprintln!("info: running transient analysis (step={step:.2e} stop={stop:.2e} method={:?} {mode})...", opts.method);
                }
                let t0 = Instant::now();
                let result = if opts.variable_step {
                    tran_nr_with_registry_var_opts(netlist, *step, *stop, registry, opts)
                } else {
                    tran_nr_with_registry_opts(netlist, *step, *stop, registry, opts)
                }.unwrap_or_else(|e| {
                    eprintln!("error: tran failed: {e}");
                    std::process::exit(1);
                });
                if ctx.verbose {
                    eprintln!("info: transient complete: {} time-points [{:.1} ms]",
                        result.time.len(), t0.elapsed().as_secs_f64() * 1000.0);
                }
                match ctx.format {
                    Format::Csv => {
                        let mut buf = Vec::new();
                        result.write_csv(&mut buf).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                        let csv = String::from_utf8_lossy(&buf);
                        let filtered = filter_csv(&csv, probe_list);
                        w.write_all(filtered.as_bytes()).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                    }
                    Format::Nutmeg => {
                        result.write_nutmeg(&mut *w, title).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                    }
                }
                if !netlist.measurements.is_empty() {
                    let ms = evaluate_measurements(&netlist.measurements, &result);
                    for m in ms {
                        eprintln!("{:<24} = {:.6e}", m.name, m.value);
                    }
                }
                ran_something = true;
            }

            Analysis::Dc { src, start, stop, step, nested } => {
                if ctx.verbose {
                    let extra = match nested {
                        Some(n) => format!(" × {} {}..{}", n.src, n.start, n.stop),
                        None    => String::new(),
                    };
                    eprintln!("info: running DC sweep on {src} ({start}..{stop} step={step}){extra}...");
                }
                let t0 = Instant::now();
                let nested_arg = nested.as_ref().map(|n| (n.src.as_str(), n.start, n.stop, n.step));
                let result = dc_sweep_with_registry_opts(
                    netlist, src, *start, *stop, *step, nested_arg, registry, opts
                ).unwrap_or_else(|e| {
                    eprintln!("error: DC sweep failed: {e}");
                    std::process::exit(1);
                });
                if ctx.verbose {
                    eprintln!("info: DC sweep complete: {} point(s) [{:.1} ms]",
                        result.n_points(), t0.elapsed().as_secs_f64() * 1000.0);
                }
                match ctx.format {
                    Format::Csv => {
                        let mut buf = Vec::new();
                        result.write_csv(&mut buf).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                        let csv = String::from_utf8_lossy(&buf);
                        let filtered = filter_csv(&csv, probe_list);
                        w.write_all(filtered.as_bytes()).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                    }
                    Format::Nutmeg => {
                        result.write_nutmeg(&mut *w, title).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                    }
                }
                ran_something = true;
            }

            Analysis::Ac { variation, points, fstart, fstop } => {
                if ctx.verbose { eprintln!("info: running AC analysis ({fstart:.2e}–{fstop:.2e} Hz, {points} pts)..."); }
                let t0 = Instant::now();
                let freqs = match variation {
                    AcVariation::Dec => freq_decade(*fstart, *fstop, *points),
                    AcVariation::Oct => freq_oct(*fstart, *fstop, *points),
                    AcVariation::Lin => freq_linear(*fstart, *fstop, *points),
                };
                let result = ac_analysis_opts(netlist, &freqs, None, registry, opts).unwrap_or_else(|e| {
                    eprintln!("error: AC analysis failed: {e}");
                    std::process::exit(1);
                });
                if ctx.verbose {
                    eprintln!("info: AC analysis complete: {} frequency points [{:.1} ms]",
                        freqs.len(), t0.elapsed().as_secs_f64() * 1000.0);
                }
                match ctx.format {
                    Format::Csv => {
                        let mut buf = Vec::new();
                        result.write_csv(&mut buf).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                        let csv = String::from_utf8_lossy(&buf);
                        let filtered = filter_csv(&csv, probe_list);
                        w.write_all(filtered.as_bytes()).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                    }
                    Format::Nutmeg => {
                        result.write_nutmeg(&mut *w, title).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                    }
                }
                ran_something = true;
            }

            Analysis::Noise { out_pos, out_neg, input_src, variation, points, fstart, fstop } => {
                if ctx.verbose {
                    eprintln!("info: running noise analysis V({out_pos},{out_neg}) on {input_src} \
                              ({fstart:.2e}–{fstop:.2e} Hz, {points} pts)...");
                }
                let t0 = Instant::now();
                let freqs = match variation {
                    AcVariation::Dec => freq_decade(*fstart, *fstop, *points),
                    AcVariation::Oct => freq_oct(*fstart, *fstop, *points),
                    AcVariation::Lin => freq_linear(*fstart, *fstop, *points),
                };
                let result = fairchild_core::noise_analysis(
                    netlist, &freqs, out_pos, out_neg, input_src, registry, opts,
                ).unwrap_or_else(|e| {
                    eprintln!("error: noise analysis failed: {e}");
                    std::process::exit(1);
                });
                if ctx.verbose {
                    eprintln!("info: noise analysis complete: {} pts [{:.1} ms]",
                        freqs.len(), t0.elapsed().as_secs_f64() * 1000.0);
                }
                let mut buf = Vec::new();
                result.write_csv(&mut buf).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                w.write_all(&buf).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                ran_something = true;
            }
        }
    }
    ran_something
}
