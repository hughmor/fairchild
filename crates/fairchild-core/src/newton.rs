use indexmap::IndexMap;

use fairchild_parser::{Element, Netlist};

use std::sync::Arc;

use crate::behavioral::BehavioralDevice;
use crate::connectivity::check_connectivity;
use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::mna::{stamp_netlist_scaled, CircuitTopology};
use crate::options::SimOptions;
use crate::solver::LinearSolver;

/// Result of a nonlinear DC operating-point solve.
pub struct NrResult {
    pub topo: CircuitTopology,
    pub x: Vec<f64>,
    pub iters: usize,
}

impl NrResult {
    pub fn node_voltage(&self, node: &str) -> Result<f64, SimError> {
        self.topo.node_voltage(node, &self.x)
    }

    pub fn vsrc_current(&self, name: &str) -> Result<f64, SimError> {
        self.topo.vsrc_current(name, &self.x)
    }

    pub fn all_voltages(&self) -> impl Iterator<Item = (&str, f64)> {
        self.topo.node_index.iter().map(|(name, &i)| (name.as_str(), self.x[i]))
    }

    /// Write the DC operating point as an ngspice-compatible Nutmeg ASCII rawfile.
    pub fn write_nutmeg<W: std::io::Write>(&self, mut w: W, title: &str) -> std::io::Result<()> {
        let n_nodes = self.topo.n_nodes();
        let n_vars = n_nodes + self.topo.vsrc_index.len();
        writeln!(w, "Title: {title}")?;
        writeln!(w, "Plotname: Operating Point")?;
        writeln!(w, "Flags: real")?;
        writeln!(w, "No. Variables: {n_vars}")?;
        writeln!(w, "No. Points: 1")?;
        writeln!(w, "Variables:")?;
        for (idx, name) in self.topo.node_index.keys().enumerate() {
            writeln!(w, "\t{idx}\tv({name})\tvoltage")?;
        }
        for (i, name) in self.topo.vsrc_index.keys().enumerate() {
            writeln!(w, "\t{}\ti({name})\tcurrent", n_nodes + i)?;
        }
        writeln!(w, "Values:")?;
        let mut first = true;
        for &idx in self.topo.node_index.values() {
            if first {
                writeln!(w, " 0\t{:.6e}", self.x[idx])?;
                first = false;
            } else {
                writeln!(w, "\t{:.6e}", self.x[idx])?;
            }
        }
        for &idx in self.topo.vsrc_index.values() {
            writeln!(w, "\t{:.6e}", self.x[n_nodes + idx])?;
        }
        Ok(())
    }

    /// Write DC operating point as a two-row CSV (header + one data row).
    pub fn write_csv<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        write!(w, "analysis")?;
        for name in self.topo.node_index.keys() {
            write!(w, ",V({name})")?;
        }
        let n_nodes = self.topo.n_nodes();
        for name in self.topo.vsrc_index.keys() {
            write!(w, ",I({name})")?;
        }
        writeln!(w)?;
        write!(w, "dc_op")?;
        for &idx in self.topo.node_index.values() {
            write!(w, ",{:.6e}", self.x[idx])?;
        }
        for &idx in self.topo.vsrc_index.values() {
            write!(w, ",{:.6e}", self.x[n_nodes + idx])?;
        }
        writeln!(w)
    }
}

/// Return refdes/instance names for each device in the order that
/// `build_devices` produces them.  Used by the verbose failure-reporter
/// to attribute residual rows to a named device.
pub fn build_device_names(netlist: &Netlist) -> Vec<String> {
    let mut names: Vec<String> = Vec::new();
    for el in &netlist.elements {
        match el {
            Element::Diode { name, model_name, .. } => {
                names.push(format!("{name} ({model_name})"));
            }
            Element::Mosfet { name, model_name, .. } => {
                names.push(format!("{name} ({model_name})"));
            }
            Element::Behavioral { name, kind, .. } => {
                names.push(format!("{name} ({kind:?})"));
            }
            Element::XOsdi { name, model_name, .. } => {
                names.push(format!("{name} ({model_name})"));
            }
            _ => {}
        }
    }
    names
}

/// Build device instances from the elements in a netlist via the registry.
pub fn build_devices(
    netlist: &Netlist,
    topo: &mut CircuitTopology,
    ctx: &SimContext,
    registry: &DeviceRegistry,
) -> Result<Vec<Box<dyn Device>>, SimError> {
    let mut devices: Vec<Box<dyn Device>> = Vec::new();
    let topo_arc = Arc::new(topo.clone());
    for el in &netlist.elements {
        match el {
            Element::Diode { anode, cathode, model_name, .. } => {
                let factory = registry.get(model_name)
                    .ok_or_else(|| SimError::UnknownModel(model_name.clone()))?;
                let pos: NodeId = topo.node_index.get(anode).copied();
                let neg: NodeId = topo.node_index.get(cathode).copied();
                devices.push(factory(&[pos, neg], ctx));
            }
            Element::Mosfet { drain, gate, source, bulk, model_name, params, .. } => {
                let d: NodeId = topo.node_index.get(drain).copied();
                let g: NodeId = topo.node_index.get(gate).copied();
                let s: NodeId = topo.node_index.get(source).copied();
                let b: NodeId = topo.node_index.get(bulk).copied();
                if let Some(dev) = registry.build_mosfet(model_name, params, &[d, g, s, b], ctx) {
                    devices.push(dev);
                } else {
                    let factory = registry.get(model_name)
                        .ok_or_else(|| SimError::UnknownModel(model_name.clone()))?;
                    devices.push(factory(&[d, g, s, b], ctx));
                }
            }
            Element::Behavioral { name, pos, neg, kind, expr } => {
                let dev = BehavioralDevice::build(
                    Arc::clone(&topo_arc),
                    name, pos, neg, *kind, expr.clone(),
                );
                devices.push(Box::new(dev));
            }
            Element::XOsdi { nets, model_name, params, .. } => {
                let factory = registry.get(model_name)
                    .ok_or_else(|| SimError::UnknownModel(model_name.clone()))?;
                let terminals: Vec<NodeId> = nets.iter()
                    .map(|net| topo.node_index.get(net).copied())
                    .collect();
                let mut dev = factory(&terminals, ctx);
                let expected = dev.num_terminals();
                if terminals.len() != expected {
                    eprintln!(
                        "warning: XOsdi '{model_name}': netlist provides {} net(s) but model \
                         expects {expected} terminal(s); extra terminals default to ground",
                        terminals.len()
                    );
                }
                for (name, value) in params {
                    dev.set_real_param(name, *value);
                }
                // OSDI models that use direct potential contributions
                // (`V(port) <+ ...`) declare internal flow-branch nodes:
                // OpenVAF surfaces these as `num_nodes − num_terminals`.
                // Allocate MNA rows for them so the OSDI runtime has real
                // slots to stamp into instead of running past mna_nodes.
                let extras = dev.num_extra_nodes();
                if extras > 0 {
                    let first = topo.allocate_extra_rows(extras);
                    dev.bind_extra_nodes(first);
                }
                devices.push(dev);
            }
            _ => {}
        }
    }
    Ok(devices)
}

/// Report MNA matrix size, NNZ, sparsity, and diagonal magnitude spread.
/// Called once at the start of `dc_op_nr_with_registry_opts` when verbose
/// is enabled — gives a sense of problem size and rough conditioning
/// before NR even starts.
/// Returns true if any Jacobian or residual cell at x=0 is NaN/Inf — a
/// signal that a device is stamping a non-finite value, which will make
/// LU factorisation fail immediately.  The caller uses this to gate the
/// (more expensive) per-device finiteness validator.
fn report_matrix_stats(
    opts: &SimOptions,
    topo: &CircuitTopology,
    netlist: &Netlist,
    devices: &mut [Box<dyn Device>],
    ctx: &SimContext,
) -> bool {
    if !opts.verbose { return false; }
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();
    let x0 = vec![0.0f64; topo.size];
    let mut mat = stamp_netlist_scaled(topo, netlist, 1.0, &empty, &empty);
    for dev in devices.iter_mut() {
        dev.eval(&x0, EvalFlags::dc(), ctx);
        dev.load_residual(&mut mat.b);
        dev.load_jacobian(&mut mat);
    }
    let n = topo.size;
    let mut nnz = 0usize;
    let mut diag_max = 0.0_f64;
    let mut diag_min = f64::INFINITY;
    let mut nonfinite = false;
    for i in 0..n {
        for j in 0..n {
            let v = mat.a[i][j];
            if v != 0.0 { nnz += 1; }
            if !v.is_finite() { nonfinite = true; }
        }
        if !mat.b[i].is_finite() { nonfinite = true; }
        let d = mat.a[i][i].abs();
        if d.is_finite() && d > 0.0 {
            if d > diag_max { diag_max = d; }
            if d < diag_min { diag_min = d; }
        }
    }
    let total = (n as f64) * (n as f64);
    let sparsity_pct = 100.0 * (1.0 - (nnz as f64) / total.max(1.0));
    eprintln!("info: MNA size: {n} rows ({} nodes + {} vsrc + {} extras); \
               nnz={nnz} ({sparsity_pct:.1}% sparse)",
        topo.n_nodes(), topo.vsrc_index.len(),
        n - topo.n_nodes() - topo.vsrc_index.len());
    if diag_min < f64::INFINITY {
        eprintln!("info: MNA diagonal magnitude spread (mixed units; large spread \
                   often indicates poor scaling): min={:.2e} max={:.2e} ratio={:.2e}",
            diag_min, diag_max, diag_max / diag_min.max(1e-300));
    }
    if nonfinite {
        eprintln!("info: MNA contains non-finite values at x=0 — running \
                   per-device finiteness validator to identify the offending \
                   device(s).");
    }
    nonfinite
}

/// Locate non-finite rows in the cumulative MNA stamp at iterate `x` and
/// cross-reference each with the netlist elements that touch the matching
/// net.  Cheaper than a per-device stamping pass — important for very
/// large circuits where allocating a fresh MnaMatrix per device would
/// cost N² × N_devices.
///
/// Walks the cumulative stamp once (O(N²)), then for each bad row scans
/// the netlist's elements (O(N_elements × max_nets_per_element)) to find
/// every X-element / R / L / C / V / I / D / M that mentions the offending
/// net.  The user gets a short list of suspect refdes per bad net.
fn validate_devices_finite(
    topo: &CircuitTopology,
    netlist: &Netlist,
    devices: &mut [Box<dyn Device>],
    ctx: &SimContext,
    x: &[f64],
) {
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();
    let mut mat = stamp_netlist_scaled(topo, netlist, 1.0, &empty, &empty);
    for dev in devices.iter_mut() {
        dev.eval(x, EvalFlags::dc(), ctx);
        dev.load_residual(&mut mat.b);
        dev.load_jacobian(&mut mat);
    }

    // Identify bad rows once (single pass; checks any non-finite cell in
    // row r OR a non-finite residual at b[r]).
    let n = topo.size;
    let mut bad_rows: Vec<usize> = Vec::new();
    for r in 0..n {
        let mut row_bad = !mat.b[r].is_finite();
        if !row_bad {
            for c in 0..n {
                if !mat.a[r][c].is_finite() {
                    row_bad = true;
                    break;
                }
            }
        }
        if row_bad { bad_rows.push(r); }
    }
    if bad_rows.is_empty() {
        eprintln!("info: validator found NO non-finite rows in the cumulative \
                   device stamp (inf must come from the linear stamp or LU).");
        return;
    }

    // Build a reverse lookup: net_index → net_name.
    let n_nodes = topo.n_nodes();
    let mut node_name_by_idx: Vec<&str> = vec![""; n_nodes];
    for (name, &i) in &topo.node_index {
        if i < n_nodes { node_name_by_idx[i] = name.as_str(); }
    }
    let mut vsrc_name_by_idx: Vec<&str> = vec![""; topo.vsrc_index.len()];
    for (name, &i) in &topo.vsrc_index {
        if i < vsrc_name_by_idx.len() { vsrc_name_by_idx[i] = name.as_str(); }
    }

    eprintln!("info: {} MNA row(s) contain non-finite stamps at x=0:",
        bad_rows.len());
    const MAX_BAD: usize = 10;
    for &r in bad_rows.iter().take(MAX_BAD) {
        let label = if r < n_nodes {
            format!("v({}) — node {r}", node_name_by_idx[r])
        } else if r < n_nodes + topo.vsrc_index.len() {
            format!("i({}) — vsrc branch", vsrc_name_by_idx[r - n_nodes])
        } else {
            format!("x[{r}] — device-internal (likely an OSDI flow-branch row)")
        };
        eprintln!("info:   row {r}: {label}");

        // Cross-reference netlist elements that touch the named net.
        if r < n_nodes {
            let net_name_lc = node_name_by_idx[r].to_lowercase();
            let mut culprits: Vec<String> = Vec::new();
            let mut suspects: Vec<String> = Vec::new();
            for el in &netlist.elements {
                let (refdes, touches, bad) = element_touches(el, &net_name_lc);
                if !touches { continue; }
                match bad {
                    Some(hint) => culprits.push(format!("{refdes} [{hint}]")),
                    None => suspects.push(refdes),
                }
            }
            if !culprits.is_empty() {
                eprintln!("info:     LIKELY CULPRIT(s): {}", culprits.join(", "));
            }
            if !suspects.is_empty() {
                let trunc = if suspects.len() > 8 {
                    format!("{} (and {} more)",
                        suspects[..8].join(", "), suspects.len() - 8)
                } else { suspects.join(", ") };
                eprintln!("info:     other elements on this net: {trunc}");
            }
        }
    }
    if bad_rows.len() > MAX_BAD {
        eprintln!("info:   ... and {} more non-finite row(s) (truncated)",
            bad_rows.len() - MAX_BAD);
    }
}

/// Return (label, touches, bad_stamp_hint) — whether `el` references the
/// given lowercase net on any of its terminals, plus a hint string if the
/// element's value alone would stamp a non-finite contribution (e.g. an
/// R=0 resistor stamps 1/0 = ∞).  Used by the finiteness validator to
/// surface the exact bad element rather than just the list of suspects.
fn element_touches(el: &Element, net_lc: &str) -> (String, bool, Option<String>) {
    let net_match = |s: &str| s.to_lowercase() == net_lc;
    let bad_value = |v: f64, what: &str| -> Option<String> {
        if v == 0.0 {
            Some(format!("{what}=0 → stamps 1/0 = ∞"))
        } else if !v.is_finite() {
            Some(format!("{what}={v} (non-finite)"))
        } else {
            None
        }
    };
    match el {
        Element::Resistor { name, pos, neg, resistance } => {
            (name.clone(), net_match(pos) || net_match(neg),
             bad_value(*resistance, "R"))
        }
        Element::Capacitor { name, pos, neg, capacitance } => {
            (name.clone(), net_match(pos) || net_match(neg),
             if !capacitance.is_finite() { Some(format!("C={capacitance} (non-finite)")) } else { None })
        }
        Element::Inductor { name, pos, neg, inductance } => {
            (name.clone(), net_match(pos) || net_match(neg),
             if !inductance.is_finite() { Some(format!("L={inductance} (non-finite)")) } else { None })
        }
        Element::VoltageSource { name, pos, neg, .. }
        | Element::CurrentSource { name, pos, neg, .. } => {
            (name.clone(), net_match(pos) || net_match(neg), None)
        }
        Element::Diode { name, anode, cathode, .. } => {
            (name.clone(), net_match(anode) || net_match(cathode), None)
        }
        Element::Mosfet { name, drain, gate, source, bulk, .. } => {
            (name.clone(),
             net_match(drain) || net_match(gate) || net_match(source) || net_match(bulk),
             None)
        }
        Element::Behavioral { name, pos, neg, .. } => {
            (name.clone(), net_match(pos) || net_match(neg), None)
        }
        Element::XOsdi { name, nets, model_name, .. } => {
            let hits = nets.iter().any(|n| net_match(n));
            (format!("{name} ({model_name})"), hits, None)
        }
    }
}

/// Reverse-lookup an MNA row index → human-readable name.
fn row_name(topo: &CircuitTopology, r: usize) -> String {
    let n_nodes = topo.n_nodes();
    if r < n_nodes {
        for (name, &i) in &topo.node_index {
            if i == r { return format!("v({name})"); }
        }
    } else if r < n_nodes + topo.vsrc_index.len() {
        for (name, &i) in &topo.vsrc_index {
            if i + n_nodes == r { return format!("i({name})"); }
        }
    } else {
        return format!("x[{r}] (device-internal)");
    }
    format!("x[{r}]")
}

/// Print the top-K rows of the residual vector with names and the
/// dominant-contributing device for each, plus the offending iterate's
/// most-changed nodes.
fn report_failure(
    phase: &str,
    opts: &SimOptions,
    topo: &CircuitTopology,
    netlist: &Netlist,
    devices: &mut [Box<dyn Device>],
    dev_names: &[String],
    ctx: &SimContext,
    x: &[f64],
    source_scale: f64,
    gmin_extra: f64,
) {
    if !opts.verbose { return; }
    // Recompute the residual at the failed iterate `x` so we can attribute
    // rows to devices.
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();
    let mut mat = stamp_netlist_scaled(topo, netlist, source_scale, &empty, &empty);
    for dev in devices.iter_mut() {
        dev.eval(x, EvalFlags::dc(), ctx);
        dev.load_residual(&mut mat.b);
    }
    let b = mat.b.clone();

    // Score each row's dominant device by probing each device into a scratch
    // residual.  O(N_devices × N_rows) but only runs on failure.
    let mut row_owner: Vec<Option<usize>> = vec![None; b.len()];
    let mut row_owner_mag: Vec<f64> = vec![0.0; b.len()];
    let mut scratch = vec![0.0f64; b.len()];
    for (di, dev) in devices.iter_mut().enumerate() {
        scratch.iter_mut().for_each(|v| *v = 0.0);
        dev.eval(x, EvalFlags::dc(), ctx);
        dev.load_residual(&mut scratch);
        for (r, &v) in scratch.iter().enumerate() {
            if v.abs() > row_owner_mag[r] {
                row_owner_mag[r] = v.abs();
                row_owner[r] = Some(di);
            }
        }
    }

    let mut idx: Vec<usize> = (0..b.len()).collect();
    idx.sort_by(|&a, &c| b[c].abs().partial_cmp(&b[a].abs()).unwrap());

    let l2: f64 = b.iter().map(|v| v * v).sum::<f64>().sqrt();
    eprintln!("info: NR did NOT converge in {phase} (residual L2 = {l2:.3e}, \
               source_scale={source_scale:.3}, gmin_extra={gmin_extra:.2e})");
    eprintln!("info: top 5 residual rows:");
    for &r in idx.iter().take(5) {
        let owner = match row_owner[r] {
            Some(d) if d < dev_names.len() => dev_names[d].as_str(),
            _ => "(linear stamp)",
        };
        eprintln!("info:   {:>4}  {:<35}  b={:>11.3e}  x={:>11.3e}  dom: {}",
            r, row_name(topo, r), b[r], x[r], owner);
    }
}

/// Compute the L2 norm of the nonlinear residual f(x) = J(x)·x − b(x)
/// at iterate `x`.  Requires a full eval + stamp at `x` — only call when
/// running line search on a clamped step.
///
/// The Norton-equivalent MNA stamp returns (J, b) such that the linearised
/// solve J·x_{k+1} = b yields the Newton step x_{k+1} − x_k = −J⁻¹·f(x_k);
/// equivalently f(x) = J(x)·x − b(x), zero at the true operating point.
fn residual_l2(
    topo: &CircuitTopology,
    netlist: &Netlist,
    devices: &mut [Box<dyn Device>],
    ctx: &SimContext,
    opts: &SimOptions,
    source_scale: f64,
    gmin_extra: f64,
    x: &[f64],
) -> f64 {
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();
    let mut mat = stamp_netlist_scaled(topo, netlist, source_scale, &empty, &empty);
    for dev in devices.iter_mut() {
        dev.eval(x, EvalFlags::dc(), ctx);
        dev.load_residual(&mut mat.b);
        dev.load_jacobian(&mut mat);
    }
    let n_nodes = topo.n_nodes();
    for i in 0..n_nodes {
        mat.a[i][i] += opts.gmin + gmin_extra;
    }
    let vsrc_end = n_nodes + topo.vsrc_index.len();
    for i in vsrc_end..topo.size {
        mat.a[i][i] += opts.gmin + gmin_extra;
    }
    let n = topo.size;
    let mut sumsq = 0.0_f64;
    for i in 0..n {
        let mut row = 0.0_f64;
        let a_row = &mat.a[i];
        for j in 0..n {
            row += a_row[j] * x[j];
        }
        let fi = row - mat.b[i];
        sumsq += fi * fi;
    }
    sumsq.sqrt()
}

/// Core Newton-Raphson loop at a fixed source scale and gmin.
///
/// `dev_names` and `phase` are only used to emit verbose diagnostics on
/// non-convergence.  Pass an empty slice + "" to suppress.
fn nr_inner(
    topo: &CircuitTopology,
    netlist: &Netlist,
    devices: &mut [Box<dyn Device>],
    ctx: &SimContext,
    opts: &SimOptions,
    solver: &dyn LinearSolver,
    mut x: Vec<f64>,
    source_scale: f64,
    gmin_extra: f64,
    dev_names: &[String],
    phase: &str,
) -> Result<Vec<f64>, SimError> {
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();
    let n_nodes = topo.n_nodes();

    for _ in 0..opts.itl1 {
        let mut mat = stamp_netlist_scaled(topo, netlist, source_scale, &empty, &empty);

        for dev in devices.iter_mut() {
            dev.eval(&x, EvalFlags::dc(), ctx);
            dev.load_residual(&mut mat.b);
            dev.load_jacobian(&mut mat);
        }

        for i in 0..n_nodes {
            mat.a[i][i] += opts.gmin + gmin_extra;
        }
        // Stamp gmin on OSDI internal-node rows too.  Their default state
        // when the device emits a degenerate contribution (e.g.
        // `OWL(out_lambda) <+ OWL(in_lambda)` where both ports map to the
        // same circuit net) is a zero row — singular without gmin.  Skip
        // voltage-source aux rows, which need a clean diagonal for their
        // own KCL equation.
        let vsrc_end = n_nodes + topo.vsrc_index.len();
        for i in vsrc_end..topo.size {
            mat.a[i][i] += opts.gmin + gmin_extra;
        }

        let x_new = match solver.solve(&mat.a, &mat.b) {
            Ok(v) => v,
            Err(e) => {
                if opts.verbose && !phase.is_empty() {
                    eprintln!("info: linear solve failed in {phase}: {e}");
                    report_failure(phase, opts, topo, netlist, devices, dev_names,
                                   ctx, &x, source_scale, gmin_extra);
                }
                return Err(e);
            }
        };

        let max_dv = x_new.iter().zip(x.iter()).take(n_nodes)
            .map(|(n, o)| (n - o).abs())
            .fold(0.0f64, f64::max);

        // Damping path. When the proposed Newton step is small enough that
        // every node's update is below the `vmax` clamp threshold, take it
        // verbatim — happy-path circuits incur zero extra cost. Otherwise
        // clamp the step to `vmax` (the existing trust-region bound) and
        // run Armijo backtracking on the L2 residual norm: pick the
        // largest α ∈ {1, 1/2, 1/4, 1/8, 1/16} for which
        //     ‖f(x + α·Δ)‖ ≤ (1 − c·α)·‖f(x)‖
        // with c = 1e-4. If no α in the budget satisfies it, fall through
        // with α_min — the next iteration's stamp will see the partial
        // step and try again.
        //
        // Residual is f(x) = J(x)·x − b(x); evaluating it requires a full
        // restamp at the trial point, so the line search costs up to 5
        // extra eval+stamp passes on each clamped iteration.
        let x_next: Vec<f64> = if max_dv > opts.vmax {
            let scale = opts.vmax / max_dv;
            // Clamped Newton step: at α=1 this is the existing vmax-clamped
            // update; Armijo lets us back off when the residual would grow.
            let delta: Vec<f64> = x.iter().zip(x_new.iter())
                .map(|(o, n)| scale * (n - o)).collect();

            let f_prev = residual_l2(topo, netlist, devices, ctx, opts,
                                     source_scale, gmin_extra, &x);
            const C_ARMIJO: f64 = 1e-4;
            const ALPHA_MIN: f64 = 1.0 / 16.0;
            let mut alpha = 1.0_f64;
            let mut x_trial: Vec<f64>;
            loop {
                x_trial = x.iter().zip(delta.iter())
                    .map(|(o, d)| o + alpha * d).collect();
                let f_trial = residual_l2(topo, netlist, devices, ctx, opts,
                                          source_scale, gmin_extra, &x_trial);
                if f_trial <= (1.0 - C_ARMIJO * alpha) * f_prev || alpha <= ALPHA_MIN {
                    break;
                }
                alpha *= 0.5;
            }
            x_trial
        } else {
            x_new
        };

        let converged = x_next.iter().zip(x.iter())
            .all(|(n, o)| (n - o).abs() < opts.vntol + opts.reltol * n.abs());

        x = x_next;
        if converged {
            return Ok(x);
        }
    }
    if opts.verbose && !phase.is_empty() {
        report_failure(phase, opts, topo, netlist, devices, dev_names, ctx,
                       &x, source_scale, gmin_extra);
    }
    Err(SimError::NoConvergence { iters: opts.itl1 })
}

/// Source-stepping homotopy: ramp sources from 0 → full value in adaptive increments.
fn source_stepping(
    topo: &CircuitTopology,
    netlist: &Netlist,
    devices: &mut [Box<dyn Device>],
    ctx: &SimContext,
    opts: &SimOptions,
    solver: &dyn LinearSolver,
    x0: Vec<f64>,
    dev_names: &[String],
) -> Result<Vec<f64>, SimError> {
    let mut x = x0;
    let mut scale = 0.0_f64;
    let mut ds = 1.0_f64 / (opts.srcsteps as f64).max(1.0);
    let min_ds = 1e-6_f64;

    while scale < 1.0 {
        let next = (scale + ds).min(1.0);
        // Suppress per-step failure reports during the inner stepping
        // loop — most failures are recoverable by halving ds.  Only the
        // outer give-up below is worth reporting.
        match nr_inner(topo, netlist, devices, ctx, opts, solver, x.clone(),
                       next, 0.0, &[], "") {
            Ok(x_new) => {
                x = x_new;
                scale = next;
                ds = (ds * 2.0).min(2.0 * (1.0 / opts.srcsteps.max(1) as f64));
            }
            Err(_) => {
                ds *= 0.5;
                if ds < min_ds {
                    if opts.verbose {
                        eprintln!("info: source-stepping gave up at scale={scale:.3} \
                                   (ds shrank to {ds:.2e} < {min_ds:.0e}); \
                                   replaying last failing step for diagnosis");
                        // Replay the last failing step with reporting on.
                        let _ = nr_inner(topo, netlist, devices, ctx, opts, solver,
                                         x.clone(), (scale + ds*2.0).min(1.0), 0.0,
                                         dev_names,
                                         &format!("source-stepping @ scale={:.3}",
                                                  (scale + ds*2.0).min(1.0)));
                    }
                    return Err(SimError::NoConvergence { iters: opts.itl1 });
                }
            }
        }
    }
    Ok(x)
}

/// GMIN stepping: add a large artificial conductance to all nodes, then ramp it down.
fn gmin_stepping(
    topo: &CircuitTopology,
    netlist: &Netlist,
    devices: &mut [Box<dyn Device>],
    ctx: &SimContext,
    opts: &SimOptions,
    solver: &dyn LinearSolver,
    dev_names: &[String],
) -> Result<Vec<f64>, SimError> {
    let mut gmin_extra = opts.gmin_max;
    let target = opts.gmin;
    let mut x = vec![0.0f64; topo.size];

    loop {
        // Quiet during the inner ramp; only the outer failure is reported.
        match nr_inner(topo, netlist, devices, ctx, opts, solver, x.clone(),
                       1.0, gmin_extra, &[], "") {
            Ok(x_new) => {
                x = x_new;
                if gmin_extra <= target { break; }
                gmin_extra = (gmin_extra * 0.1).max(target);
            }
            Err(_) => {
                if opts.verbose {
                    eprintln!("info: gmin-stepping FAILED at gmin_extra={gmin_extra:.2e} \
                               (target gmin={target:.2e}); replaying for diagnosis");
                    let _ = nr_inner(topo, netlist, devices, ctx, opts, solver,
                                     x.clone(), 1.0, gmin_extra, dev_names,
                                     &format!("gmin-stepping @ gmin_extra={gmin_extra:.2e}"));
                }
                return Err(SimError::NoConvergence { iters: opts.itl1 });
            }
        }
    }
    Ok(x)
}

/// Build an initial `x` vector for DC NR, seeding from any `.nodeset` entries
/// in the netlist.  Unknown nodes are silently ignored.
fn build_x0_from_nodeset(netlist: &Netlist, topo: &CircuitTopology) -> Vec<f64> {
    let mut x0 = vec![0.0f64; topo.size];
    for (name, value) in &netlist.nodeset {
        if let Some(&i) = topo.node_index.get(name) {
            x0[i] = *value;
        }
    }
    x0
}

/// DC operating-point with explicit `SimOptions`.
///
/// Convergence strategy (in order):
///   1. Direct Newton-Raphson from the `.nodeset` seed (or x=0).
///   2. Source stepping: ramp sources from 0 → full value.
///   3. GMIN stepping: add large diagonal conductance, ramp to standard GMIN.
pub fn dc_op_nr_with_registry_opts(
    netlist: &Netlist,
    registry: &DeviceRegistry,
    opts: &SimOptions,
) -> Result<NrResult, SimError> {
    if opts.sanity_check {
        crate::sanity::check_netlist_sanity(netlist);
    }
    check_connectivity(netlist)?;
    let ctx = opts.sim_context();
    let mut topo = CircuitTopology::build(netlist);

    let mut devices = build_devices(netlist, &mut topo, &ctx, registry)?;
    let dev_names = build_device_names(netlist);
    let x0 = build_x0_from_nodeset(netlist, &topo);
    let solver = opts.linear_solver(topo.size);

    if opts.verbose {
        let nonfinite = report_matrix_stats(opts, &topo, netlist, &mut devices, &ctx);
        if nonfinite {
            let x_probe = vec![0.0f64; topo.size];
            validate_devices_finite(&topo, netlist, &mut devices, &ctx, &x_probe);
        }
        eprintln!("info: DC OP: trying direct Newton-Raphson from \
                   {} seed...",
            if !netlist.nodeset.is_empty() {"nodeset"} else {"x=0"});
    }
    if let Ok(x) = nr_inner(&topo, netlist, &mut devices, &ctx, opts, &*solver,
                            x0.clone(), 1.0, 0.0, &dev_names, "direct NR (DC OP)") {
        if opts.verbose { eprintln!("info: DC OP: direct NR succeeded"); }
        return Ok(NrResult { topo, x, iters: 1 });
    }

    if opts.verbose { eprintln!("info: DC OP: direct NR failed; trying source-stepping..."); }
    if let Ok(x) = source_stepping(&topo, netlist, &mut devices, &ctx, opts, &*solver,
                                   x0, &dev_names) {
        if opts.verbose { eprintln!("info: DC OP: source-stepping succeeded"); }
        return Ok(NrResult { topo, x, iters: 2 });
    }

    if opts.verbose { eprintln!("info: DC OP: source-stepping failed; trying gmin-stepping..."); }
    match gmin_stepping(&topo, netlist, &mut devices, &ctx, opts, &*solver, &dev_names) {
        Ok(x) => {
            if opts.verbose { eprintln!("info: DC OP: gmin-stepping succeeded"); }
            Ok(NrResult { topo, x, iters: 3 })
        }
        Err(e) => Err(e),
    }
}

/// DC operating-point with options taken from any `.options` directives in
/// the netlist (defaults where unspecified).
///
/// This is the recommended entry point for CLI/Python callers: it honours
/// `.options reltol=… gmin=… method=…` automatically.
pub fn dc_op_nr_with_registry(
    netlist: &Netlist,
    registry: &DeviceRegistry,
) -> Result<NrResult, SimError> {
    dc_op_nr_with_registry_opts(netlist, registry, &SimOptions::from_netlist(netlist))
}

/// DC operating-point with pre-built devices (for sweeps / parametric analysis).
pub fn dc_op_nr_with_devices_opts(
    netlist: &Netlist,
    topo: &CircuitTopology,
    devices: &mut Vec<Box<dyn Device>>,
    ctx: &SimContext,
    opts: &SimOptions,
) -> Result<NrResult, SimError> {
    let x0 = vec![0.0f64; topo.size];
    let solver = opts.linear_solver(topo.size);
    let dev_names = build_device_names(netlist);

    if let Ok(x) = nr_inner(topo, netlist, devices, ctx, opts, &*solver,
                            x0.clone(), 1.0, 0.0, &dev_names, "direct NR (DC OP)") {
        return Ok(NrResult { topo: topo.clone(), x, iters: 1 });
    }

    if let Ok(x) = source_stepping(topo, netlist, devices, ctx, opts, &*solver,
                                   x0, &dev_names) {
        return Ok(NrResult { topo: topo.clone(), x, iters: 2 });
    }

    match gmin_stepping(topo, netlist, devices, ctx, opts, &*solver, &dev_names) {
        Ok(x) => Ok(NrResult { topo: topo.clone(), x, iters: 3 }),
        Err(e) => Err(e),
    }
}

/// DC operating-point with pre-built devices, default options.
pub fn dc_op_nr_with_devices(
    netlist: &Netlist,
    topo: &CircuitTopology,
    devices: &mut Vec<Box<dyn Device>>,
    ctx: &SimContext,
) -> Result<NrResult, SimError> {
    dc_op_nr_with_devices_opts(netlist, topo, devices, ctx, &SimOptions::default())
}

/// DC operating-point using only built-in models, default options.
pub fn dc_op_nr(netlist: &Netlist) -> Result<NrResult, SimError> {
    dc_op_nr_opts(netlist, &SimOptions::default())
}

/// DC operating-point using only built-in models, with explicit options.
pub fn dc_op_nr_opts(netlist: &Netlist, opts: &SimOptions) -> Result<NrResult, SimError> {
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_diodes(&netlist.models);
    registry.register_builtin_mosfets(&netlist.models);
    dc_op_nr_with_registry_opts(netlist, &registry, opts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

    #[test]
    fn linear_circuit_converges() {
        let net = parse_spice(
            "* divider\nV1 in 0 DC 1.0\nR1 in out 1k\nR2 out 0 1k\n.op\n.end\n",
        ).unwrap();
        let r = dc_op_nr(&net).unwrap();
        assert!(r.iters <= 5, "linear circuit should converge fast, took {}", r.iters);
        let v = r.node_voltage("out").unwrap();
        assert!((v - 0.5).abs() < 1e-6, "v(out)={v}");
    }

    #[test]
    fn pnjlim_helps_stiff_diode() {
        // 20 V across a tiny series R into a diode: an easy case for an
        // ideal solver, a stiff one without junction limiting.  With pnjlim
        // ON (default) this must converge in well under ITL1=150 iters.
        let net = parse_spice(
            "* stiff diode\nVdd a 0 DC 20\nR1 a b 100\nD1 b 0 myd\n\
             .model myd D (Is=1e-15 N=1)\n.op\n.end\n"
        ).unwrap();
        let r = dc_op_nr(&net).unwrap();
        // The pn-junction equation gives V(b) ≈ 18·V_T·ln(I/Is) ≈ 0.7 V,
        // and most of the supply drops across R1.
        let vb = r.node_voltage("b").unwrap();
        assert!(vb > 0.5 && vb < 1.0, "V(b)={vb:.4} should be ~0.7V");
        assert!(r.iters < 50, "convergence took {} iters with pnjlim on", r.iters);
    }

    #[test]
    fn current_source_biased_diode() {
        let net = parse_spice(
            "* Diode bias\nIb 0 b 1m\nD1 b 0 myd\n.model myd D (Is=1e-14 N=1)\n.op\n.end\n",
        ).unwrap();
        let r = dc_op_nr(&net).unwrap();
        let vb = r.node_voltage("b").unwrap();
        let vt = 1.380649e-23 * 300.15 / 1.602176634e-19;
        let expected = vt * (1e-3_f64 / 1e-14_f64 + 1.0).ln();
        let tol = 1e-4 * expected;
        assert!(
            (vb - expected).abs() < tol,
            "V(b)={vb:.6e}  expected={expected:.6e}  diff={:.2e}",
            (vb - expected).abs()
        );
    }

    #[test]
    fn series_rd_circuit() {
        let net = parse_spice(
            "* R-D series\nVdd a 0 DC 5\nR1 a b 10k\nD1 b 0 myd\n\
             .model myd D (Is=1e-14 N=1)\n.op\n.end\n",
        ).unwrap();
        let r = dc_op_nr(&net).unwrap();
        let vb = r.node_voltage("b").unwrap();
        assert!(vb > 0.5 && vb < 0.8, "V(b) out of expected range: {vb}");
        let vt = 1.380649e-23 * 300.15 / 1.602176634e-19;
        let ir = (5.0 - vb) / 10e3;
        let id = 1e-14 * ((vb / vt).exp() - 1.0);
        assert!((ir - id).abs() < 1e-8, "KCL error: ir={ir:.4e} id={id:.4e}");
    }

    #[test]
    fn write_csv_dc_op() {
        let net = parse_spice(
            "* divider\nV1 in 0 DC 2.0\nR1 in out 1k\nR2 out 0 1k\n.op\n.end\n",
        ).unwrap();
        let r = dc_op_nr(&net).unwrap();
        let mut buf = Vec::new();
        r.write_csv(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("analysis,"), "header: {s}");
        assert!(s.contains("V(out)"), "should contain V(out): {s}");
        assert!(s.contains("dc_op"), "should have dc_op row: {s}");
        assert!(s.contains("1.000000e0") || s.contains("1.000000e+0"), "V(out)≈1V missing: {s}");
    }

    #[test]
    fn write_nutmeg_dc_op() {
        let net = parse_spice(
            "* divider\nV1 in 0 DC 2.0\nR1 in out 1k\nR2 out 0 1k\n.op\n.end\n",
        ).unwrap();
        let r = dc_op_nr(&net).unwrap();
        let mut buf = Vec::new();
        r.write_nutmeg(&mut buf, "divider test").unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("Plotname: Operating Point"), "plotname: {s}");
        assert!(s.contains("Flags: real"), "flags: {s}");
        assert!(s.contains("v(out)\tvoltage"), "v(out): {s}");
        assert!(s.contains("No. Points: 1"), "single point: {s}");
    }

    #[test]
    fn tighter_reltol_takes_more_iterations() {
        // A circuit with mild nonlinearity. With opts.reltol=1e-10 NR should still
        // converge but use more iterations than at the default 1e-3.
        let net = parse_spice(
            "* R-D\nVdd a 0 DC 5\nR1 a b 10k\nD1 b 0 myd\n\
             .model myd D (Is=1e-14 N=1)\n.op\n.end\n",
        ).unwrap();
        let mut opts = SimOptions::default();
        opts.reltol = 1e-10;
        opts.vntol  = 1e-12;
        let r = dc_op_nr_opts(&net, &opts).unwrap();
        let vb = r.node_voltage("b").unwrap();
        assert!(vb > 0.5 && vb < 0.8);
    }
}
