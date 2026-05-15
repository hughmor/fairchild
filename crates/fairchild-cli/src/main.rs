use std::fs;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use clap::{Parser, ValueEnum};

use fairchild_core::{
    ac_analysis_opts, dc_op_nr_with_registry_opts, dc_sweep_with_registry_opts,
    evaluate_measurements,
    freq_decade, freq_linear, freq_oct,
    tran_nr_with_registry_opts, DeviceRegistry, SimOptions,
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
    /// maxstep, gminmax, srcsteps, method (be|tr|gear), uic, temp.
    ///
    /// Example: --opt reltol=1e-5 --opt method=gear
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

    if let Some(v) = &cli.reltol   { apply("reltol",  v); }
    if let Some(v) = &cli.gmin     { apply("gmin",    v); }
    if let Some(v) = &cli.method   { apply("method",  v); }
    if let Some(v) = &cli.max_step { apply("maxstep", v); }
    if cli.no_pnjlim               { apply("pnjlim",  "0"); }

    for raw in &cli.options {
        if let Some((k, v)) = raw.split_once('=') {
            apply(k.trim(), v.trim().trim_matches('"').trim_matches('\''));
        } else if !cli.quiet {
            eprintln!("warning: --opt '{raw}': expected KEY=VALUE, skipping");
        }
    }

    opts
}

// ── OSDI registry builder ─────────────────────────────────────────────────

/// Load built-in models + any OSDI shared libraries listed in the netlist.
/// Relative `.osdi` paths are resolved against `netlist_dir`.
fn build_registry(netlist: &Netlist, netlist_dir: Option<&PathBuf>, _quiet: bool) -> DeviceRegistry {
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_diodes(&netlist.models);
    registry.register_builtin_mosfets(&netlist.models);

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

    // Build device registry: built-in models + OSDI shared libraries
    let netlist_dir = cli.file.parent().map(|p| p.to_path_buf());
    let registry = build_registry(&netlist, netlist_dir.as_ref(), cli.quiet);

    // Merge netlist `.options` + CLI flag overrides into a single SimOptions.
    let opts = build_options(&netlist, &cli);
    if cli.verbose {
        eprintln!("info: solver options: reltol={:e} gmin={:e} method={:?} itl1={} itl4={}",
            opts.reltol, opts.gmin, opts.method, opts.itl1, opts.itl4);
    }

    let writer: Box<dyn Write> = match &cli.output {
        Some(path) => {
            let f = fs::File::create(path).unwrap_or_else(|e| {
                eprintln!("error: cannot create {:?}: {e}", path);
                std::process::exit(1);
            });
            Box::new(BufWriter::new(f))
        }
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let mut w = writer;

    let mut ran_something = false;

    for analysis in &netlist.analyses {
        match analysis {
            Analysis::Op => {
                if cli.verbose { eprintln!("info: running DC operating-point analysis..."); }
                let t0 = Instant::now();
                let result = dc_op_nr_with_registry_opts(&netlist, &registry, &opts).unwrap_or_else(|e| {
                    eprintln!("error: DC op failed: {e}");
                    std::process::exit(1);
                });
                if cli.verbose {
                    eprintln!("info: DC op converged in {} iteration(s) [{:.1} ms]",
                        result.iters, t0.elapsed().as_secs_f64() * 1000.0);
                }
                match cli.format {
                    Format::Csv => {
                        let mut buf = Vec::new();
                        result.write_csv(&mut buf).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                        let csv = String::from_utf8_lossy(&buf);
                        let filtered = filter_csv(&csv, &probe_list);
                        w.write_all(filtered.as_bytes()).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                    }
                    Format::Nutmeg => {
                        result.write_nutmeg(&mut w, &title).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                    }
                }
                ran_something = true;
            }

            Analysis::Tran { step, stop } => {
                if cli.verbose { eprintln!("info: running transient analysis (step={step:.2e} stop={stop:.2e} method={:?})...", opts.method); }
                let t0 = Instant::now();
                let result = tran_nr_with_registry_opts(&netlist, *step, *stop, &registry, &opts).unwrap_or_else(|e| {
                    eprintln!("error: tran failed: {e}");
                    std::process::exit(1);
                });
                if cli.verbose {
                    eprintln!("info: transient complete: {} time-points [{:.1} ms]",
                        result.time.len(), t0.elapsed().as_secs_f64() * 1000.0);
                }
                match cli.format {
                    Format::Csv => {
                        let mut buf = Vec::new();
                        result.write_csv(&mut buf).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                        let csv = String::from_utf8_lossy(&buf);
                        let filtered = filter_csv(&csv, &probe_list);
                        w.write_all(filtered.as_bytes()).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                    }
                    Format::Nutmeg => {
                        result.write_nutmeg(&mut w, &title).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                    }
                }
                // Evaluate `.measure` directives on the just-finished tran run.
                if !netlist.measurements.is_empty() {
                    let ms = evaluate_measurements(&netlist.measurements, &result);
                    for m in ms {
                        eprintln!("{:<24} = {:.6e}", m.name, m.value);
                    }
                }
                ran_something = true;
            }

            Analysis::Dc { src, start, stop, step, nested } => {
                if cli.verbose {
                    let extra = match nested {
                        Some(n) => format!(" × {} {}..{}", n.src, n.start, n.stop),
                        None    => String::new(),
                    };
                    eprintln!("info: running DC sweep on {src} ({start}..{stop} step={step}){extra}...");
                }
                let t0 = Instant::now();
                let nested_arg = nested.as_ref().map(|n| (n.src.as_str(), n.start, n.stop, n.step));
                let result = dc_sweep_with_registry_opts(
                    &netlist, src, *start, *stop, *step, nested_arg, &registry, &opts
                ).unwrap_or_else(|e| {
                    eprintln!("error: DC sweep failed: {e}");
                    std::process::exit(1);
                });
                if cli.verbose {
                    eprintln!("info: DC sweep complete: {} point(s) [{:.1} ms]",
                        result.n_points(), t0.elapsed().as_secs_f64() * 1000.0);
                }
                match cli.format {
                    Format::Csv => {
                        let mut buf = Vec::new();
                        result.write_csv(&mut buf).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                        let csv = String::from_utf8_lossy(&buf);
                        let filtered = filter_csv(&csv, &probe_list);
                        w.write_all(filtered.as_bytes()).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                    }
                    Format::Nutmeg => {
                        result.write_nutmeg(&mut w, &title).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                    }
                }
                ran_something = true;
            }

            Analysis::Ac { variation, points, fstart, fstop } => {
                if cli.verbose { eprintln!("info: running AC analysis ({fstart:.2e}–{fstop:.2e} Hz, {points} pts)..."); }
                let t0 = Instant::now();
                let freqs = match variation {
                    AcVariation::Dec => freq_decade(*fstart, *fstop, *points),
                    AcVariation::Oct => freq_oct(*fstart, *fstop, *points),
                    AcVariation::Lin => freq_linear(*fstart, *fstop, *points),
                };
                let result = ac_analysis_opts(&netlist, &freqs, None, &registry, &opts).unwrap_or_else(|e| {
                    eprintln!("error: AC analysis failed: {e}");
                    std::process::exit(1);
                });
                if cli.verbose {
                    eprintln!("info: AC analysis complete: {} frequency points [{:.1} ms]",
                        freqs.len(), t0.elapsed().as_secs_f64() * 1000.0);
                }
                match cli.format {
                    Format::Csv => {
                        let mut buf = Vec::new();
                        result.write_csv(&mut buf).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                        let csv = String::from_utf8_lossy(&buf);
                        let filtered = filter_csv(&csv, &probe_list);
                        w.write_all(filtered.as_bytes()).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                    }
                    Format::Nutmeg => {
                        result.write_nutmeg(&mut w, &title).unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                    }
                }
                ran_something = true;
            }
        }
    }

    if !ran_something && !cli.quiet {
        eprintln!("warning: no analyses found in netlist (add .op, .tran, or .ac)");
    }
}
