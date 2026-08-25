use crate::warn_user;
use indexmap::IndexMap;

use fairchild_parser::{Element, Netlist};

use std::sync::Arc;

use crate::behavioral::BehavioralDevice;
use crate::connectivity::check_connectivity;
use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::device_registry::{DeviceRegistry, ParamSet};
use crate::error::SimError;
use crate::mna::{stamp_netlist_scaled, CircuitTopology, Footprint, RowFloor};
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

    /// Every voltage-like signal: solved node voltages, then the λ labels
    /// (resolved before the solve, see `CircuitTopology::lambda_signals`) —
    /// so enumeration and by-name probing name the same set (#71).
    pub fn all_voltages(&self) -> impl Iterator<Item = (&str, f64)> {
        self.topo
            .node_index
            .iter()
            .map(|(name, &i)| (name.as_str(), self.x[i]))
            .chain(self.topo.lambda_signals())
    }

    /// Write the DC operating point as an ngspice-compatible Nutmeg ASCII rawfile.
    ///
    /// The λ labels ride along as voltage variables after the solved nodes:
    /// ngspice would have carried these nets, and this file format is the one
    /// with consumers outside this repo (#71). The variable count and the
    /// `Variables:` block are two statements of one number, so both are
    /// derived from the same three sets.
    pub fn write_nutmeg<W: std::io::Write>(&self, mut w: W, title: &str) -> std::io::Result<()> {
        let n_nodes = self.topo.n_nodes();
        let lambda = self.topo.lambda_signals();
        let n_vars = n_nodes + lambda.len() + self.topo.vsrc_index.len();
        writeln!(w, "Title: {title}")?;
        writeln!(w, "Plotname: Operating Point")?;
        writeln!(w, "Flags: real")?;
        writeln!(w, "No. Variables: {n_vars}")?;
        writeln!(w, "No. Points: 1")?;
        writeln!(w, "Variables:")?;
        let mut idx = 0usize;
        for name in self.topo.node_index.keys() {
            writeln!(w, "\t{idx}\tv({name})\tvoltage")?;
            idx += 1;
        }
        for (name, _) in &lambda {
            writeln!(w, "\t{idx}\tv({name})\tvoltage")?;
            idx += 1;
        }
        for name in self.topo.vsrc_index.keys() {
            writeln!(w, "\t{idx}\ti({name})\tcurrent")?;
            idx += 1;
        }
        writeln!(w, "Values:")?;
        let values = self
            .topo
            .node_index
            .values()
            .map(|&idx| self.x[idx])
            .chain(lambda.iter().map(|&(_, wl)| wl))
            .chain(
                self.topo
                    .vsrc_index
                    .values()
                    .map(|&idx| self.x[n_nodes + idx]),
            );
        for (k, v) in values.enumerate() {
            if k == 0 {
                writeln!(w, " 0\t{v:.6e}")?;
            } else {
                writeln!(w, "\t{v:.6e}")?;
            }
        }
        Ok(())
    }

    /// Write DC operating point as a two-row CSV (header + one data row).
    /// Column order: solved node voltages, λ labels, branch currents.
    pub fn write_csv<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        let lambda = self.topo.lambda_signals();
        write!(w, "analysis")?;
        for name in self.topo.node_index.keys() {
            write!(w, ",V({name})")?;
        }
        for (name, _) in &lambda {
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
        for &(_, wl) in &lambda {
            write!(w, ",{wl:.6e}")?;
        }
        for &idx in self.topo.vsrc_index.values() {
            write!(w, ",{:.6e}", self.x[n_nodes + idx])?;
        }
        writeln!(w)
    }
}

/// Return refdes/instance names for each device in the order that
/// `build_devices` produces them.  Used by the verbose failure-reporter to
/// attribute residual rows to a named device, and by `crate::adjoint` to map a
/// `ParamRef`'s element name onto a live device.
///
/// **This must stay arm-for-arm parallel with `build_devices_with_footprints`.**
/// It was not: switches and transmission lines became devices and never got a
/// name here, so every device after one in the list was attributed to the wrong
/// element.  `device_element_names` below is the same walk over element names
/// alone, and `names_stay_parallel_to_devices` pins both against a netlist that
/// exercises every arm.
pub fn build_device_names(netlist: &Netlist) -> Vec<String> {
    let mut names: Vec<String> = Vec::new();
    for el in &netlist.elements {
        match el {
            Element::Diode {
                name, model_name, ..
            }
            | Element::Mosfet {
                name, model_name, ..
            }
            | Element::Bjt {
                name, model_name, ..
            }
            | Element::XOsdi {
                name, model_name, ..
            }
            | Element::VoltageSwitch {
                name, model_name, ..
            }
            | Element::CurrentSwitch {
                name, model_name, ..
            } => {
                names.push(format!("{name} ({model_name})"));
            }
            Element::Behavioral { name, kind, .. } => {
                names.push(format!("{name} ({kind:?})"));
            }
            Element::TransmissionLine { name, .. } => {
                names.push(format!("{name} (T)"));
            }
            _ => {}
        }
    }
    names
}

/// The bare element name for each device, same order as `build_devices`.
pub fn device_element_names(netlist: &Netlist) -> Vec<String> {
    let mut names: Vec<String> = Vec::new();
    for el in &netlist.elements {
        match el {
            Element::Diode { name, .. }
            | Element::Mosfet { name, .. }
            | Element::Bjt { name, .. }
            | Element::XOsdi { name, .. }
            | Element::VoltageSwitch { name, .. }
            | Element::CurrentSwitch { name, .. }
            | Element::Behavioral { name, .. }
            | Element::TransmissionLine { name, .. } => names.push(name.clone()),
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
    build_devices_with_footprints(netlist, topo, ctx, registry).map(|(d, _)| d)
}

/// Name the element a construction failure belongs to.
///
/// A factory's `Err` describes what is wrong with the *device* — "reflectance +
/// transmittance = 1.4, must be ≤ 1" — because that is all a device knows. Which
/// element and which model is the caller's half, and this is the only place that
/// has both.
fn attribute(
    built: Result<Box<dyn Device>, String>,
    name: &str,
    model_name: &str,
) -> Result<Box<dyn Device>, SimError> {
    built.map_err(|e| SimError::ParameterError(format!("{name} ('{model_name}'): {e}")))
}

/// Register a built device: allocate its extra MNA rows, bind them, and record
/// the structural footprint `mna::Pattern` needs — every row/column it can
/// stamp into.
fn push_device(
    devices: &mut Vec<Box<dyn Device>>,
    foot: &mut Vec<Footprint>,
    topo: &mut CircuitTopology,
    terminals: &[NodeId],
    mut dev: Box<dyn Device>,
) {
    let mut rows: Vec<usize> = terminals.iter().flatten().copied().collect();
    let extras = dev.num_extra_nodes();
    let mut extra_first = None;
    if extras > 0 {
        let first = topo.allocate_extra_rows(extras);
        dev.bind_extra_nodes(first);
        rows.extend(first..first + extras);
        extra_first = Some(first);
    }
    // Which of this device's rows carry kelvin. Asked after the extras are
    // allocated, since a self-heating node is usually one of them, and resolved
    // here because this is the only point where a device's own node numbering
    // and its MNA rows are both in hand. A terminal tied to ground is `None` and
    // is skipped: row 0 does not exist, and a temperature clamped to ambient has
    // no unknown to bound.
    for k in dev.thermal_nodes() {
        let row = match k.checked_sub(terminals.len()) {
            None => terminals[k],
            Some(off) => extra_first.filter(|_| off < extras).map(|f| f + off),
        };
        if let Some(r) = row {
            topo.thermal_rows.push(r);
        }
    }
    rows.extend(dev.extra_stamp_rows());
    rows.sort_unstable();
    rows.dedup();
    // Asked only after the extra rows are bound, so a device can name them.
    foot.push(match dev.stamp_pairs() {
        Some(pairs) => Footprint::Pairs(pairs),
        None => Footprint::Clique(rows),
    });
    devices.push(dev);
}

/// Look up a switch's `.model` card and build the device from it.
///
/// Shared by the `S` and `W` arms of `build_devices_with_footprints`, which
/// differ only in what they bind as the control.
fn build_switch(
    registry: &DeviceRegistry,
    model_name: &str,
    inst_name: &str,
    initial_on: bool,
) -> Result<crate::models::Switch, SimError> {
    let (is_current, params) = registry
        .switch_cards
        .get(model_name)
        .ok_or_else(|| SimError::UnknownModel(model_name.to_string()))?;
    crate::models::Switch::from_model_params(*is_current, params, initial_on)
        .map(|(dev, _)| dev)
        .map_err(|e| SimError::ParameterError(format!("{inst_name}: {e}")))
}

/// A device list paired with each device's structural MNA footprint.
pub type DevicesWithFootprints = (Vec<Box<dyn Device>>, Vec<Footprint>);

/// [`build_devices`] plus each device's structural footprint, in device order.
/// The hot NR / transient loops need the footprints to build a sparsity
/// pattern; everything else can use the simpler wrapper.
pub fn build_devices_with_footprints(
    netlist: &Netlist,
    topo: &mut CircuitTopology,
    ctx: &SimContext,
    registry: &DeviceRegistry,
) -> Result<DevicesWithFootprints, SimError> {
    let mut devices: Vec<Box<dyn Device>> = Vec::new();
    let mut foot: Vec<Footprint> = Vec::new();
    // Every auxiliary row allocated below lands at or after this index, in
    // device order — `check_exclusive_potential_drivers` walks the same order
    // to attribute a row back to the device that owns it.
    let extras_base = topo.size;
    let topo_arc = Arc::new(topo.clone());
    // Which element each device came from, so the λ pass below can hand a
    // device the wavelengths on its own nets. Filled by watching the vector
    // grow rather than threaded through `push_device`'s dozen call sites — one
    // list, kept true by construction rather than by a promise.
    let mut dev_elem: Vec<usize> = Vec::new();
    for (elem_idx, el) in netlist.elements.iter().enumerate() {
        match el {
            Element::Diode {
                name,
                anode,
                cathode,
                model_name,
                params,
            } => {
                let factory = registry
                    .get(model_name)
                    .ok_or_else(|| SimError::UnknownModel(model_name.clone()))?;
                let pos: NodeId = topo.node_index.get(anode).copied();
                let neg: NodeId = topo.node_index.get(cathode).copied();
                let ps = ParamSet::new(params);
                let dev = attribute(factory(&[pos, neg], &ps, ctx), name, model_name)?;
                // A diode instance parameter used to reach the netlist and stop:
                // `D1 a k dm area=2` parsed, changed nothing, and said nothing.
                // AREA is honoured now; anything else gets named here.
                for p in ps.unconsumed() {
                    match &p.card {
                        Some(card) => warn_user!(
                            "{name}: .model '{card}' parameter '{}' is not honoured by \
                             '{model_name}' and was dropped",
                            p.key
                        ),
                        None => warn_user!(
                            "{name} ('{model_name}'): instance parameter '{}' is not \
                             honoured by this model and was dropped",
                            p.key
                        ),
                    }
                }
                push_device(&mut devices, &mut foot, topo, &[pos, neg], dev);
            }
            Element::Mosfet {
                name,
                drain,
                gate,
                source,
                bulk,
                model_name,
                params,
            } => {
                let d: NodeId = topo.node_index.get(drain).copied();
                let g: NodeId = topo.node_index.get(gate).copied();
                let s: NodeId = topo.node_index.get(source).copied();
                let b: NodeId = topo.node_index.get(bulk).copied();
                if let Some(dev) =
                    registry.build_mosfet(model_name, name, params, &[d, g, s, b], ctx)
                {
                    push_device(&mut devices, &mut foot, topo, &[d, g, s, b], dev);
                } else {
                    let factory = registry
                        .get(model_name)
                        .ok_or_else(|| SimError::UnknownModel(model_name.clone()))?;
                    let dev = attribute(
                        factory(&[d, g, s, b], &ParamSet::new(params), ctx),
                        name,
                        model_name,
                    )?;
                    push_device(&mut devices, &mut foot, topo, &[d, g, s, b], dev);
                }
            }
            Element::Bjt {
                name,
                collector,
                base,
                emitter,
                substrate,
                model_name,
                params,
            } => {
                let c: NodeId = topo.node_index.get(collector).copied();
                let b: NodeId = topo.node_index.get(base).copied();
                let e: NodeId = topo.node_index.get(emitter).copied();
                let s: NodeId = topo.node_index.get(substrate).copied(); // typically ground
                let dev = registry
                    .build_bjt(model_name, name, params, &[c, b, e, s], ctx)
                    .ok_or_else(|| SimError::UnknownModel(model_name.clone()))?;
                // RB/RC/RE series resistances declare internal nodes (one per
                // non-zero resistance); push_device allocates and binds them.
                push_device(&mut devices, &mut foot, topo, &[c, b, e, s], dev);
            }
            Element::Behavioral {
                name,
                pos,
                neg,
                kind,
                expr,
            } => {
                let dev = BehavioralDevice::build(
                    Arc::clone(&topo_arc),
                    name,
                    pos,
                    neg,
                    *kind,
                    expr.clone(),
                );
                // Terminals empty on purpose: a B-element reports its whole
                // footprint (terminals, aux row, referenced nodes) itself.
                push_device(&mut devices, &mut foot, topo, &[], Box::new(dev));
            }
            Element::XOsdi {
                name,
                nets,
                model_name,
                params,
            } => {
                let factory = registry
                    .get(model_name)
                    .ok_or_else(|| SimError::UnknownModel(model_name.clone()))?;
                let terminals: Vec<NodeId> = nets
                    .iter()
                    .map(|net| topo.node_index.get(net).copied())
                    .collect();
                // build() applies the instance params (and the model-card
                // defaults baked into the factory). It also tracks which params
                // the device consumed, so we can warn about typos.
                let ps = ParamSet::new(params);
                let dev = attribute(factory(&terminals, &ps, ctx), name, model_name)?;
                let expected = dev.num_terminals();
                if terminals.len() != expected {
                    // Used to be a warning that grounded the missing terminals
                    // and carried on, which is how a mis-wired photonic device
                    // reached its own assert! and panicked out through pyo3 —
                    // or worse, silently simulated a circuit nobody drew. A
                    // wrong port count is a netlist error, so say so.
                    return Err(SimError::ParameterError(format!(
                        "{name}: '{model_name}' expects {expected} terminal(s) but the \
                         netlist gives {}. Check the element's port order against the \
                         model's card — for photonic devices the optical ports come \
                         first, then the electrical ones.",
                        terminals.len()
                    )));
                }
                for p in ps.unconsumed() {
                    match &p.card {
                        Some(card) => warn_user!(
                            "{name}: .model '{card}' parameter '{}' is unknown to \
                             '{model_name}' and was ignored",
                            p.key
                        ),
                        None => warn_user!(
                            "'{model_name}' instance: unknown parameter '{}' ignored",
                            p.key
                        ),
                    }
                }
                // OSDI models that use direct potential contributions
                // (`V(port) <+ ...`) declare internal flow-branch nodes:
                // OpenVAF surfaces these as `num_nodes − num_terminals`.
                // push_device allocates MNA rows for them so the OSDI runtime
                // has real slots to stamp into instead of running past
                // mna_nodes.
                push_device(&mut devices, &mut foot, topo, &terminals, dev);
            }
            Element::VoltageSwitch {
                name,
                pos,
                neg,
                ctrl_pos,
                ctrl_neg,
                model_name,
                initial_on,
            } => {
                let node = |n: &fairchild_parser::NodeName| topo.node_index.get(n).copied();
                let terms = [node(pos), node(neg), node(ctrl_pos), node(ctrl_neg)];
                let dev = build_switch(registry, model_name, name, *initial_on)?;
                let mut dev: Box<dyn Device> = Box::new(dev);
                dev.setup_model(ctx);
                dev.setup_instance(&terms, ctx);
                push_device(&mut devices, &mut foot, topo, &terms, dev);
            }
            Element::CurrentSwitch {
                name,
                pos,
                neg,
                ctrl_vsrc,
                model_name,
                initial_on,
            } => {
                let node = |n: &fairchild_parser::NodeName| topo.node_index.get(n).copied();
                let terms = [node(pos), node(neg)];
                let mut dev = build_switch(registry, model_name, name, *initial_on)?;
                // Resolving the controlling branch row is the whole reason a
                // switch is not a plain registry factory: nothing below
                // `build_devices` knows `vsrc_index`.
                let vname = ctrl_vsrc.to_lowercase();
                let idx = topo.vsrc_index.get(&vname).copied().ok_or_else(|| {
                    SimError::ParameterError(format!(
                        "{name}: controlling source '{vname}' is not a voltage source in this netlist"
                    ))
                })?;
                dev.set_control(crate::models::SwitchControl::Current {
                    row: Some(topo.n_nodes() + idx),
                });
                let mut dev: Box<dyn Device> = Box::new(dev);
                dev.setup_model(ctx);
                dev.setup_instance(&terms, ctx);
                push_device(&mut devices, &mut foot, topo, &terms, dev);
            }
            Element::TransmissionLine {
                a_pos,
                a_neg,
                b_pos,
                b_neg,
                z0,
                td,
                ..
            } => {
                let term = |n: &fairchild_parser::NodeName| topo.node_index.get(n).copied();
                let terms = [term(a_pos), term(a_neg), term(b_pos), term(b_neg)];
                let mut dev: Box<dyn Device> =
                    Box::new(crate::models::tline::NativeTLine::new(*z0, *td));
                dev.setup_model(ctx);
                dev.setup_instance(&terms, ctx);
                // Two branch-current rows (i1, i2), allocated by push_device.
                push_device(&mut devices, &mut foot, topo, &terms, dev);
            }
            _ => {}
        }
        while dev_elem.len() < devices.len() {
            dev_elem.push(elem_idx);
        }
    }
    apply_resolved_lambda(netlist, &mut devices, &dev_elem, topo, ctx, registry);
    check_exclusive_potential_drivers(&devices, netlist, topo, extras_base)?;
    Ok((devices, foot))
}

/// Hand every device the wavelength resolved for each of its terminals.
///
/// λ is a label, not a state: it is routed from sources, never computed, so it
/// is resolved once here rather than solved for (see [`crate::lambda`]). Doing
/// it inside the builder is what makes it unforgettable — a caller that
/// assembles devices by hand cannot end up with a photonic device evaluating
/// its phase at the wrong colour.
///
/// A terminal no λ net reached takes the band centre, which is exactly what an
/// undriven λ wire used to bootstrap to.
fn apply_resolved_lambda(
    netlist: &Netlist,
    devices: &mut [Box<dyn Device>],
    dev_elem: &[usize],
    topo: &CircuitTopology,
    ctx: &SimContext,
    registry: &DeviceRegistry,
) {
    // Normally the topology already resolved λ — it had to, to know which nets
    // were not rows. Resolve again only for a caller that built its topology
    // with `CircuitTopology::build`, so a hand-assembled device list still
    // evaluates at the right wavelength instead of the band centre.
    let owned;
    let map = if topo.lambda.is_empty() {
        owned = crate::lambda::resolve(netlist, ctx, registry);
        &owned
    } else {
        &topo.lambda
    };
    let mut per_terminal: Vec<f64> = Vec::new();
    for (dev, &ei) in devices.iter_mut().zip(dev_elem) {
        let Some(Element::XOsdi { nets, .. }) = netlist.elements.get(ei) else {
            continue;
        };
        per_terminal.clear();
        per_terminal.extend(
            nets.iter()
                .map(|n| map.get(n).unwrap_or(ctx.lambda_center_m)),
        );
        dev.set_resolved_lambda(&per_terminal);
    }
}

/// Refuse a netlist in which two devices pin the same node's potential.
///
/// A device that drives a node through an auxiliary branch row is asserting
/// exclusive ownership of it. When a second device asserts the same node, the
/// two columns differ only by the `gmin` on their own diagonals, so the block
/// is rank-deficient in a way LU never reports: the factorisation succeeds and
/// returns a `gmin`-weighted average of the two assertions, with no error and
/// no warning. That has shipped twice — a laser driving its backward wire to
/// zero while the waveguide drove the returning field onto it (4x low), and
/// `fc_mux`/`fc_demux` stamping the backward path in the forward direction.
///
/// The ownership is not declared anywhere, so this reads it off the stamp
/// rather than off a second list that could disagree with it: every driven
/// potential goes through `models::photonic::stamp_potential_eq` (and the OSDI
/// runtime's equivalent), which writes exactly `a[row][node] = a[node][row] = 1`
/// for its own auxiliary `row`. Only the owning device writes its own rows and
/// columns, so that cell pair identifies both the node and the owner.
///
/// One device pinning a node twice is a bug in that device, not in the netlist,
/// and is left to the device's own tests.
fn check_exclusive_potential_drivers(
    devices: &[Box<dyn Device>],
    netlist: &Netlist,
    topo: &CircuitTopology,
    extras_base: usize,
) -> Result<(), SimError> {
    if topo.size == extras_base {
        return Ok(());
    }
    // Stamped before any `eval`, so the coefficients are whatever the devices
    // were constructed with. That does not matter: the `1` marking the driven
    // node is written unconditionally, independent of every coefficient.
    let mut scratch = crate::mna::MnaMatrix::zeros(topo.size);
    for dev in devices {
        dev.load_jacobian(&mut scratch);
    }
    let mut owner: Vec<usize> = vec![usize::MAX; topo.size - extras_base];
    let mut next = extras_base;
    for (i, dev) in devices.iter().enumerate() {
        let extras = dev.num_extra_nodes();
        owner[next - extras_base..next + extras - extras_base].fill(i);
        next += extras;
    }
    let n_nodes = topo.n_nodes();
    let mut pinned: IndexMap<usize, usize> = IndexMap::new();
    let mut names: Option<Vec<String>> = None;
    for row in extras_base..topo.size {
        for (node, v) in scratch.a[row].iter() {
            if node >= n_nodes || v != 1.0 || scratch.a[node][row] != 1.0 {
                continue;
            }
            let dev = owner[row - extras_base];
            let Some(&first) = pinned.get(&node) else {
                pinned.insert(node, dev);
                continue;
            };
            if first == dev {
                continue;
            }
            let names = names.get_or_insert_with(|| build_device_names(netlist));
            let name = |i: usize| names.get(i).cloned().unwrap_or_else(|| "?".into());
            return Err(SimError::OverdrivenNode {
                node: topo
                    .node_index
                    .get_index(node)
                    .map(|(n, _)| n.clone())
                    .unwrap_or_else(|| format!("row {node}")),
                first: name(first),
                second: name(dev),
            });
        }
    }
    Ok(())
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
    if !opts.verbose && !opts.cond_estimate {
        return false;
    }
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();
    let x0 = vec![0.0f64; topo.size];
    let mut mat = stamp_netlist_scaled(
        topo,
        netlist,
        1.0,
        &empty,
        &empty,
        crate::mna::InductorDc::Short,
    );
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
            if v != 0.0 {
                nnz += 1;
            }
            if !v.is_finite() {
                nonfinite = true;
            }
        }
        if !mat.b[i].is_finite() {
            nonfinite = true;
        }
        let d = mat.a[i][i].abs();
        if d.is_finite() && d > 0.0 {
            if d > diag_max {
                diag_max = d;
            }
            if d < diag_min {
                diag_min = d;
            }
        }
    }
    // Structurally empty rows and columns, named. A singular matrix is
    // otherwise almost impossible to localise from the outside: the failure
    // surfaces as an unsatisfied *linear* row somewhere else entirely, because
    // the solve returns garbage for the whole system. gmin is stamped first so
    // node rows that only float are not reported — what is left is a row no
    // device ever wrote to, or an unknown nothing depends on.
    {
        let mut a = mat.a.clone();
        topo.stamp_gmin(&mut a, opts.gmin.max(1e-12), RowFloor::PinEmptyRows);
        let mut empty_rows: Vec<usize> = Vec::new();
        let mut col_nz = vec![0usize; n];
        // Rows are sparse: iterating yields the stored (column, value) pairs, so
        // count structural entries directly rather than scanning n columns.
        for (i, row) in a.iter().enumerate() {
            let mut row_nz = 0usize;
            for (j, v) in row.iter() {
                if v != 0.0 {
                    row_nz += 1;
                    col_nz[j] += 1;
                }
            }
            if row_nz == 0 {
                empty_rows.push(i);
            }
        }
        let empty_cols: Vec<usize> = (0..n).filter(|&j| col_nz[j] == 0).collect();
        for (what, rows) in [("row", &empty_rows), ("column", &empty_cols)] {
            if rows.is_empty() {
                continue;
            }
            eprintln!(
                "info: {} structurally empty MNA {what}(s) — the matrix is \
                       singular and every solve from it is meaningless:",
                rows.len()
            );
            for &r in rows.iter().take(12) {
                eprintln!("info:   {r:>6}  {}", row_name(topo, r));
            }
            if rows.len() > 12 {
                eprintln!("info:   … and {} more", rows.len() - 12);
            }
        }
    }
    let total = (n as f64) * (n as f64);
    let sparsity_pct = 100.0 * (1.0 - (nnz as f64) / total.max(1.0));
    eprintln!(
        "info: MNA size: {n} rows ({} nodes + {} vsrc + {} extras); \
               nnz={nnz} ({sparsity_pct:.1}% sparse)",
        topo.n_nodes(),
        topo.vsrc_index.len(),
        n - topo.n_nodes() - topo.vsrc_index.len()
    );
    if diag_min < f64::INFINITY {
        eprintln!(
            "info: MNA diagonal magnitude spread (mixed units; large spread \
                   often indicates poor scaling): min={:.2e} max={:.2e} ratio={:.2e}",
            diag_min,
            diag_max,
            diag_max / diag_min.max(1e-300)
        );
    }
    if opts.cond_estimate {
        // Estimate κ on the matrix the solver actually factorises (with the
        // gmin floor on node / internal rows), so floating nodes don't show as
        // spuriously singular.  Work on a copy to leave `mat` untouched.
        let mut a_est = mat.a.clone();
        topo.stamp_gmin(&mut a_est, opts.gmin, RowFloor::PinEmptyRows);
        match crate::solver::estimate_condition_2norm(&a_est) {
            Some(k) => eprintln!(
                "info: estimated 2-norm condition number κ(A) ≈ {k:.3e} \
                       (κ ≫ 1e8 indicates ill-conditioning — try .options equilibrate=1)"
            ),
            None => {
                eprintln!("info: condition-number estimate unavailable (matrix singular at x=0)")
            }
        }
    }
    if nonfinite {
        eprintln!(
            "info: MNA contains non-finite values at x=0 — running \
                   per-device finiteness validator to identify the offending \
                   device(s)."
        );
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
    let mut mat = stamp_netlist_scaled(
        topo,
        netlist,
        1.0,
        &empty,
        &empty,
        crate::mna::InductorDc::Short,
    );
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
        if row_bad {
            bad_rows.push(r);
        }
    }
    if bad_rows.is_empty() {
        eprintln!(
            "info: validator found NO non-finite rows in the cumulative \
                   device stamp (inf must come from the linear stamp or LU)."
        );
        return;
    }

    // Build a reverse lookup: net_index → net_name.
    let n_nodes = topo.n_nodes();
    let mut node_name_by_idx: Vec<&str> = vec![""; n_nodes];
    for (name, &i) in &topo.node_index {
        if i < n_nodes {
            node_name_by_idx[i] = name.as_str();
        }
    }
    let mut vsrc_name_by_idx: Vec<&str> = vec![""; topo.vsrc_index.len()];
    for (name, &i) in &topo.vsrc_index {
        if i < vsrc_name_by_idx.len() {
            vsrc_name_by_idx[i] = name.as_str();
        }
    }

    eprintln!(
        "info: {} MNA row(s) contain non-finite stamps at x=0:",
        bad_rows.len()
    );
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
                if !touches {
                    continue;
                }
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
                    format!(
                        "{} (and {} more)",
                        suspects[..8].join(", "),
                        suspects.len() - 8
                    )
                } else {
                    suspects.join(", ")
                };
                eprintln!("info:     other elements on this net: {trunc}");
            }
        }
    }
    if bad_rows.len() > MAX_BAD {
        eprintln!(
            "info:   ... and {} more non-finite row(s) (truncated)",
            bad_rows.len() - MAX_BAD
        );
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
        Element::Resistor {
            name,
            pos,
            neg,
            resistance,
        } => (
            name.clone(),
            net_match(pos) || net_match(neg),
            bad_value(*resistance, "R"),
        ),
        Element::Capacitor {
            name,
            pos,
            neg,
            capacitance,
        } => (
            name.clone(),
            net_match(pos) || net_match(neg),
            if !capacitance.is_finite() {
                Some(format!("C={capacitance} (non-finite)"))
            } else {
                None
            },
        ),
        Element::Inductor {
            name,
            pos,
            neg,
            inductance,
        } => (
            name.clone(),
            net_match(pos) || net_match(neg),
            if !inductance.is_finite() {
                Some(format!("L={inductance} (non-finite)"))
            } else {
                None
            },
        ),
        Element::VoltageSource { name, pos, neg, .. }
        | Element::CurrentSource { name, pos, neg, .. } => {
            (name.clone(), net_match(pos) || net_match(neg), None)
        }
        Element::Diode {
            name,
            anode,
            cathode,
            ..
        } => (name.clone(), net_match(anode) || net_match(cathode), None),
        Element::Mosfet {
            name,
            drain,
            gate,
            source,
            bulk,
            ..
        } => (
            name.clone(),
            net_match(drain) || net_match(gate) || net_match(source) || net_match(bulk),
            None,
        ),
        Element::Bjt {
            name,
            collector,
            base,
            emitter,
            substrate,
            ..
        } => (
            name.clone(),
            net_match(collector) || net_match(base) || net_match(emitter) || net_match(substrate),
            None,
        ),
        Element::Behavioral { name, pos, neg, .. } => {
            (name.clone(), net_match(pos) || net_match(neg), None)
        }
        Element::XOsdi {
            name,
            nets,
            model_name,
            ..
        } => {
            let hits = nets.iter().any(|n| net_match(n));
            (format!("{name} ({model_name})"), hits, None)
        }
        Element::CoupledInductors { name, .. } => {
            // K elements reference inductor names, not net names directly.
            (name.clone(), false, None)
        }
        Element::VoltageSwitch {
            name,
            pos,
            neg,
            ctrl_pos,
            ctrl_neg,
            model_name,
            ..
        } => (
            format!("{name} ({model_name})"),
            net_match(pos) || net_match(neg) || net_match(ctrl_pos) || net_match(ctrl_neg),
            None,
        ),
        Element::CurrentSwitch {
            name,
            pos,
            neg,
            model_name,
            ..
        } => (
            format!("{name} ({model_name})"),
            net_match(pos) || net_match(neg),
            None,
        ),
        Element::TransmissionLine {
            name,
            a_pos,
            a_neg,
            b_pos,
            b_neg,
            ..
        } => (
            name.clone(),
            net_match(a_pos) || net_match(a_neg) || net_match(b_pos) || net_match(b_neg),
            None,
        ),
    }
}

/// Reverse-lookup an MNA row index → human-readable name.
fn row_name(topo: &CircuitTopology, r: usize) -> String {
    let n_nodes = topo.n_nodes();
    if r < n_nodes {
        for (name, &i) in &topo.node_index {
            if i == r {
                return format!("v({name})");
            }
        }
    } else if r < n_nodes + topo.vsrc_index.len() {
        for (name, &i) in &topo.vsrc_index {
            if i + n_nodes == r {
                return format!("i({name})");
            }
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
    if !opts.verbose {
        return;
    }
    // Recompute the residual at the failed iterate `x` so we can attribute
    // rows to devices.
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();
    let mut mat = stamp_netlist_scaled(
        topo,
        netlist,
        source_scale,
        &empty,
        &empty,
        crate::mna::InductorDc::Short,
    );
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
    eprintln!(
        "info: NR did NOT converge in {phase} (residual L2 = {l2:.3e}, \
               source_scale={source_scale:.3}, gmin_extra={gmin_extra:.2e})"
    );
    // 5 rows is too few on a large circuit: the voltage-source KVL rows carry
    // the biggest |b| and crowd out every device row, which is exactly the
    // wrong picture when the sources are fine and a device is the problem.
    let show = if opts.verbose { 20 } else { 5 };
    eprintln!("info: top {show} residual rows:");
    for &r in idx.iter().take(show) {
        let owner = match row_owner[r] {
            Some(d) if d < dev_names.len() => dev_names[d].as_str(),
            _ => "(linear stamp)",
        };
        eprintln!(
            "info:   {:>4}  {:<35}  b={:>11.3e}  x={:>11.3e}  dom: {}",
            r,
            row_name(topo, r),
            b[r],
            x[r],
            owner
        );
    }
}

/// Compute the L2 norm of the nonlinear residual f(x) = J(x)·x − b(x)
/// at iterate `x`.  Requires a full eval + stamp at `x` — only call when
/// running line search on a clamped step.
///
/// The Norton-equivalent MNA stamp returns (J, b) such that the linearised
/// solve J·x_{k+1} = b yields the Newton step x_{k+1} − x_k = −J⁻¹·f(x_k);
/// equivalently f(x) = J(x)·x − b(x), zero at the true operating point.
/// ‖f(x)‖₂ at an arbitrary trial point, for the Armijo line search.
///
/// Stamps into a caller-owned `scratch` matrix so the line search neither
/// allocates an n×n matrix per trial nor walks the full dense row for the
/// matvec — both were O(n²) per call, and this is the hottest path in a large
/// photonic solve.  ‖f(x)‖ at the *current* iterate needs none of this: the
/// NR loop's own matrix is already stamped there, so call
/// `MnaMatrix::residual_norm` on it directly.
/// Leaves `scratch` holding the linearisation at `x` — `crate::adjoint`
/// relies on that to read the Jacobian back out without a second stamp.
#[allow(clippy::too_many_arguments)]
pub(crate) fn residual_l2(
    scratch: &mut crate::mna::MnaMatrix,
    topo: &CircuitTopology,
    netlist: &Netlist,
    devices: &mut [Box<dyn Device>],
    ctx: &SimContext,
    opts: &SimOptions,
    source_scale: f64,
    gmin_extra: f64,
    plan: Option<&crate::mna::StampPlan>,
    x: &[f64],
) -> f64 {
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();
    crate::mna::stamp_netlist_scaled_in_place(
        scratch,
        topo,
        netlist,
        source_scale,
        &empty,
        &empty,
        plan,
        crate::mna::InductorDc::Short,
    );
    for dev in devices.iter_mut() {
        dev.set_source_scale(source_scale);
        dev.eval(x, EvalFlags::dc(), ctx);
        dev.load_residual(&mut scratch.b);
        dev.load_jacobian(scratch);
    }
    topo.stamp_gmin(
        &mut scratch.a,
        opts.gmin + gmin_extra,
        RowFloor::PinEmptyRows,
    );
    scratch.residual_norm(x)
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
    plan: Option<&crate::mna::StampPlan>,
) -> Result<Vec<f64>, SimError> {
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();
    let n_nodes = topo.n_nodes();
    // Not every unknown is a volt — see `crate::tolerance`.  Built here rather
    // than per iteration; `topo.size` is settled by the time any solver runs.
    let tol = crate::tolerance::Tolerances::build(topo, opts);

    // Sparsity pattern is fixed across this NR loop — devices stamp the
    // same matrix positions each iteration, only the values change.  The
    // factorisation cache captures the symbolic LU once on iteration 0
    // and reuses it for every subsequent solve: KLU's `klu_refactor` fast
    // path, or faer's `Lu::try_new_with_symbolic`.
    let mut fact: Option<Box<dyn crate::solver::Factorisation>> = None;

    // Reuse one MnaMatrix across iterations — saves N+1 heap allocations
    // per NR step versus rebuilding from scratch each time.  With a plan
    // attached it also carries the structural pattern, which turns the
    // per-iteration clear and the solver's dense→CSC scan into O(nnz).
    let mut mat = match plan {
        Some(p) => crate::mna::MnaMatrix::with_pattern(topo.size, p.pattern.clone()),
        None => crate::mna::MnaMatrix::zeros(topo.size),
    };
    // Second matrix for the Armijo trial points; allocated only if the line
    // search actually runs (it needs a clamped step first).
    let mut trial: Option<crate::mna::MnaMatrix> = None;
    let mut first_stamp = true;
    let mut warned_clamp = false;

    for _iter in 0..opts.itl1 {
        crate::mna::stamp_netlist_scaled_in_place(
            &mut mat,
            topo,
            netlist,
            source_scale,
            &empty,
            &empty,
            plan,
            crate::mna::InductorDc::Short,
        );

        for dev in devices.iter_mut() {
            // Independent-source devices (lasers) ramp with the homotopy too;
            // no-op for everything else.
            dev.set_source_scale(source_scale);
            dev.eval(&x, EvalFlags::dc(), ctx);
            dev.load_residual(&mut mat.b);
            dev.load_jacobian(&mut mat);
        }

        // gmin floor on node + device-internal rows (skips vsource aux rows);
        // `gmin_extra` is the homotopy step. See CircuitTopology::stamp_gmin.
        topo.stamp_gmin(&mut mat.a, opts.gmin + gmin_extra, RowFloor::PinEmptyRows);

        if first_stamp {
            // The pattern is a structural superset, so a stamped cell outside
            // it means some device coupled nodes it was never handed.  Cheap
            // insurance against a silently wrong sparse solve.
            #[cfg(debug_assertions)]
            mat.debug_assert_covers();
            first_stamp = false;
        }

        let solve_result = if let Some(f) = fact.as_mut() {
            f.refactor_and_solve_mat(&mat)
        } else {
            // First iteration of this NR loop: build the cache.
            match solver.factorise_mat(&mat) {
                Ok(mut f) => {
                    let x = f.refactor_and_solve_mat(&mat);
                    fact = Some(f);
                    x
                }
                Err(e) => Err(e),
            }
        };
        let x_new = match solve_result {
            Ok(v) => v,
            Err(e) => {
                if opts.verbose && !phase.is_empty() {
                    eprintln!("info: linear solve failed in {phase}: {e}");
                    report_failure(
                        phase,
                        opts,
                        topo,
                        netlist,
                        devices,
                        dev_names,
                        ctx,
                        &x,
                        source_scale,
                        gmin_extra,
                    );
                }
                return Err(e);
            }
        };

        // `vmax` is a limit in VOLTS, so only rows that carry volts may set it.
        //
        // λ wires used to be exempted here: they carry metres (~1.55e-6), and a
        // volt-scaled clamp once shrank every λ in a circuit to 1e-19 m because
        // one under-constrained heater node wanted 1e12 V. λ is no longer an
        // unknown at all (see `crate::lambda`) — but thermal rows are, and they
        // fail the same way in the other direction: a device settling 40 K above
        // ambient asks for a 40-unit step, which sets `max_dv` and scales
        // *every* electrical unknown down by 80×. The circuit still converges,
        // and takes a great many more iterations to do it for no reason.
        //
        // So the exemption is by unit, not by name: a row `Tolerances` bounds in
        // something other than volts does not get a vote on a volt-scaled trust
        // region. Its own equation still constrains it — it is excluded from
        // setting the clamp, not from being clamped.
        // A slice scan, not a set: `thermal_rows` holds one entry per thermal
        // node in the circuit — empty in every non-thermal deck — and this runs
        // inside the Newton loop, where building a hash set per iteration would
        // cost more than the search it replaces.
        let thermal = topo.thermal_rows.as_slice();
        let mut max_dv = 0.0f64;
        let mut max_dv_row = 0usize;
        for i in 0..n_nodes {
            if thermal.contains(&i) {
                continue;
            }
            let d = (x_new[i] - x[i]).abs();
            if d > max_dv {
                max_dv = d;
                max_dv_row = i;
            }
        }

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
        // Residual is f(x) = J(x)·x − b(x); evaluating it needs a restamp at
        // each trial point, so the line search costs up to 5 extra eval+stamp
        // passes on each clamped iteration.  ‖f(x)‖ itself is free — `mat` is
        // already stamped there.
        //
        // Declared per iteration: a value carried over from an earlier clamped
        // step would block convergence for the rest of the solve.
        let mut armijo_fell_back = false;
        let x_next: Vec<f64> = if max_dv > opts.vmax {
            let scale = opts.vmax / max_dv;
            // Clamped Newton step: at α=1 this is the existing vmax-clamped
            // update; Armijo lets us back off when the residual would grow.
            let delta: Vec<f64> = x
                .iter()
                .zip(x_new.iter())
                .map(|(o, n)| scale * (n - o))
                .collect();
            // A scale this small means the step carries no information: the
            // offending row is almost certainly under-constrained rather than
            // merely fast. Say which row, because "did not converge" sends the
            // reader looking at the physics instead.
            if scale < 1e-6 && !warned_clamp {
                warned_clamp = true;
                // Describe the symptom and let the reader draw the conclusion:
                // the first version of this asserted "no resistive path to
                // ground", which was wrong on the very circuit that prompted it
                // (the path existed; the inductors on it were being stamped open).
                warn_user!(
                    "node '{}' wants to move {max_dv:.3e} V in one Newton \
                     step, so the vmax={:.3e} trust region shrinks every unknown by \
                     {scale:.3e} — including unknowns in other units. A step this \
                     small carries almost no information, so the solve will crawl \
                     or stall. I/gmin ({:.3e} V at gmin={:.1e}) is what a node \
                     carrying no conductance at all would show.",
                    row_name(topo, max_dv_row),
                    opts.vmax,
                    1.0 / opts.gmin,
                    opts.gmin
                );
            }

            // `mat` still holds J(x) and b(x) from the stamp at the top of this
            // iteration, so ‖f(x)‖ is free — no restamp.
            let scratch = trial.get_or_insert_with(|| match plan {
                Some(p) => crate::mna::MnaMatrix::with_pattern(topo.size, p.pattern.clone()),
                None => crate::mna::MnaMatrix::zeros(topo.size),
            });
            // Measured the same way as the trial residuals below, deliberately.
            // Reading it off `mat` instead would be free — it is already
            // stamped at `x` — but device limiters (pnjlim) carry state across
            // evals, so a residual from the loop's stamp and one from a fresh
            // restamp are not the same quantity, and Armijo would be comparing
            // two different functions.  The saving here comes from making each
            // evaluation O(nnz) rather than from skipping one.
            let f_prev = residual_l2(
                scratch,
                topo,
                netlist,
                devices,
                ctx,
                opts,
                source_scale,
                gmin_extra,
                plan,
                &x,
            );
            const C_ARMIJO: f64 = 1e-4;
            const ALPHA_MIN: f64 = 1.0 / 16.0;
            let mut alpha = 1.0_f64;
            let mut x_trial: Vec<f64>;
            loop {
                armijo_fell_back = alpha <= ALPHA_MIN;
                x_trial = x
                    .iter()
                    .zip(delta.iter())
                    .map(|(o, d)| o + alpha * d)
                    .collect();
                let f_trial = residual_l2(
                    scratch,
                    topo,
                    netlist,
                    devices,
                    ctx,
                    opts,
                    source_scale,
                    gmin_extra,
                    plan,
                    &x_trial,
                );
                if f_trial <= (1.0 - C_ARMIJO * alpha) * f_prev || alpha <= ALPHA_MIN {
                    break;
                }
                alpha *= 0.5;
            }
            x_trial
        } else {
            x_new
        };

        // Never call it converged on a step the line search could not justify.
        //
        // When no alpha satisfies Armijo the loop falls through at ALPHA_MIN, so
        // the step becomes a CONSTANT vmax/16 — independent of the Newton
        // direction's magnitude. The iterate then marches at fixed velocity, and
        // the relative test `|dx| < abstol + reltol*|x|` is satisfied as soon as
        // |x| passes (vmax/16)/reltol ~ 31 V. That is not convergence, it is the
        // tolerance catching up with a stalled walk, and it returned a +56% wrong
        // operating point as a success: a photodetector shunt fed by a current
        // source read 31.25 V where the answer is 19.985 V, at every iteration
        // limit, because the stopping point depends on reltol rather than on the
        // circuit. Refusing to converge here turns a confidently wrong answer
        // into an honest failure that the homotopy can then try to fix.
        // NOTE: this test can still pass on a step the `vmax` trust region cut
        // short, and then it means nothing — `|dx|` is the clamp, not the
        // distance to the solution, and `abstol + reltol·|x|` grows with `|x|`
        // until the two meet. A node heading for hundreds of volts stops
        // wherever the walk happens to be, and is reported as converged. An
        // ideal VCCS into 1 kΩ — linear, one solve from the answer — reads
        // 0.0502 V on a node the deck pins at 0.1 V. `.options vmax=1e5` gives
        // the right answer, which is the diagnosis.
        //
        // Adding `&& max_dv <= opts.vmax` here does fix it, and costs a 15×
        // slowdown on this repository's own test suite: circuits that were
        // reaching a right answer on a clamped step now exhaust `itl1` and fall
        // into homotopy at every timestep. The honest fix is a residual-based
        // convergence test rather than a step-based one, which is a solver
        // change too broad to make from here. Tracked in #90.
        let converged = tol.converged(&x_next, &x) && !armijo_fell_back;

        x = x_next;
        if converged {
            return Ok(x);
        }
    }
    if opts.verbose && !phase.is_empty() {
        report_failure(
            phase,
            opts,
            topo,
            netlist,
            devices,
            dev_names,
            ctx,
            &x,
            source_scale,
            gmin_extra,
        );
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
    plan: Option<&crate::mna::StampPlan>,
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
        match nr_inner(
            topo,
            netlist,
            devices,
            ctx,
            opts,
            solver,
            x.clone(),
            next,
            0.0,
            &[],
            "",
            plan,
        ) {
            Ok(x_new) => {
                x = x_new;
                scale = next;
                ds = (ds * 2.0).min(2.0 * (1.0 / opts.srcsteps.max(1) as f64));
            }
            Err(_) => {
                ds *= 0.5;
                if ds < min_ds {
                    if opts.verbose {
                        eprintln!(
                            "info: source-stepping gave up at scale={scale:.3} \
                                   (ds shrank to {ds:.2e} < {min_ds:.0e}); \
                                   replaying last failing step for diagnosis"
                        );
                        // Replay the last failing step with reporting on.
                        let _ = nr_inner(
                            topo,
                            netlist,
                            devices,
                            ctx,
                            opts,
                            solver,
                            x.clone(),
                            (scale + ds * 2.0).min(1.0),
                            0.0,
                            dev_names,
                            &format!("source-stepping @ scale={:.3}", (scale + ds * 2.0).min(1.0)),
                            plan,
                        );
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
    plan: Option<&crate::mna::StampPlan>,
) -> Result<Vec<f64>, SimError> {
    let mut gmin_extra = opts.gmin_max;
    let target = opts.gmin;
    let mut x = vec![0.0f64; topo.size];

    loop {
        // Quiet during the inner ramp; only the outer failure is reported.
        match nr_inner(
            topo,
            netlist,
            devices,
            ctx,
            opts,
            solver,
            x.clone(),
            1.0,
            gmin_extra,
            &[],
            "",
            plan,
        ) {
            Ok(x_new) => {
                x = x_new;
                if gmin_extra <= target {
                    break;
                }
                gmin_extra = (gmin_extra * 0.1).max(target);
            }
            Err(e) => {
                // A singular matrix is not a convergence problem, and reporting
                // it as one sends the user hunting for a bias point that cannot
                // exist. Two voltage sources in parallel with different values,
                // or a voltage-source loop, reach here after every homotopy
                // stage has failed the same way — continuation cannot rescue a
                // topology that has no solution, so pass the real diagnosis up.
                let singular = matches!(e, SimError::SingularMatrix);
                if opts.verbose {
                    eprintln!(
                        "info: gmin-stepping FAILED at gmin_extra={gmin_extra:.2e} \
                               (target gmin={target:.2e}); replaying for diagnosis"
                    );
                    let _ = nr_inner(
                        topo,
                        netlist,
                        devices,
                        ctx,
                        opts,
                        solver,
                        x.clone(),
                        1.0,
                        gmin_extra,
                        dev_names,
                        &format!("gmin-stepping @ gmin_extra={gmin_extra:.2e}"),
                        plan,
                    );
                }
                if singular {
                    return Err(SimError::SingularMatrix);
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
    let mut topo = CircuitTopology::build_resolved(netlist, &ctx, registry);

    let (mut devices, footprints) =
        build_devices_with_footprints(netlist, &mut topo, &ctx, registry)?;
    let dev_names = build_device_names(netlist);
    let x0 = build_x0_from_nodeset(netlist, &topo);
    let solver = opts.linear_solver(topo.size);
    // Built after build_devices: extra device rows have to be allocated first,
    // since they widen topo.size.
    let plan = crate::mna::StampPlan::new(&topo, netlist, &footprints);
    plan.resolve_device_cells(&mut devices);

    if opts.verbose {
        let nonfinite = report_matrix_stats(opts, &topo, netlist, &mut devices, &ctx);
        if nonfinite {
            let x_probe = vec![0.0f64; topo.size];
            validate_devices_finite(&topo, netlist, &mut devices, &ctx, &x_probe);
        }
        eprintln!(
            "info: DC OP: trying direct Newton-Raphson from \
                   {} seed...",
            if !netlist.nodeset.is_empty() {
                "nodeset"
            } else {
                "x=0"
            }
        );
    }
    if let Ok(x) = nr_inner(
        &topo,
        netlist,
        &mut devices,
        &ctx,
        opts,
        &*solver,
        x0.clone(),
        1.0,
        0.0,
        &dev_names,
        "direct NR (DC OP)",
        Some(&plan),
    ) {
        if opts.verbose {
            eprintln!("info: DC OP: direct NR succeeded");
        }
        return Ok(NrResult { topo, x, iters: 1 });
    }

    if opts.verbose {
        eprintln!("info: DC OP: direct NR failed; trying source-stepping...");
    }
    if let Ok(x) = source_stepping(
        &topo,
        netlist,
        &mut devices,
        &ctx,
        opts,
        &*solver,
        x0,
        &dev_names,
        Some(&plan),
    ) {
        if opts.verbose {
            eprintln!("info: DC OP: source-stepping succeeded");
        }
        return Ok(NrResult { topo, x, iters: 2 });
    }

    if opts.verbose {
        eprintln!("info: DC OP: source-stepping failed; trying gmin-stepping...");
    }
    match gmin_stepping(
        &topo,
        netlist,
        &mut devices,
        &ctx,
        opts,
        &*solver,
        &dev_names,
        Some(&plan),
    ) {
        Ok(x) => {
            if opts.verbose {
                eprintln!("info: DC OP: gmin-stepping succeeded");
            }
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
    // No sparsity pattern on this path: the caller handed over pre-built
    // devices, so the per-device terminal footprints the pattern needs are
    // gone.  Falls back to the dense clear + scan.  Sweeps that care can call
    // `dc_op_nr_with_registry_opts` per point instead — they already run in
    // parallel across points.
    let plan: Option<&crate::mna::StampPlan> = None;

    if let Ok(x) = nr_inner(
        topo,
        netlist,
        devices,
        ctx,
        opts,
        &*solver,
        x0.clone(),
        1.0,
        0.0,
        &dev_names,
        "direct NR (DC OP)",
        plan,
    ) {
        return Ok(NrResult {
            topo: topo.clone(),
            x,
            iters: 1,
        });
    }

    if let Ok(x) = source_stepping(
        topo, netlist, devices, ctx, opts, &*solver, x0, &dev_names, plan,
    ) {
        return Ok(NrResult {
            topo: topo.clone(),
            x,
            iters: 2,
        });
    }

    match gmin_stepping(
        topo, netlist, devices, ctx, opts, &*solver, &dev_names, plan,
    ) {
        Ok(x) => Ok(NrResult {
            topo: topo.clone(),
            x,
            iters: 3,
        }),
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
    registry.register_builtin_models(&netlist.models);
    dc_op_nr_with_registry_opts(netlist, &registry, opts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

    #[test]
    fn linear_circuit_converges() {
        let net =
            parse_spice("* divider\nV1 in 0 DC 1.0\nR1 in out 1k\nR2 out 0 1k\n.op\n").unwrap();
        let r = dc_op_nr(&net).unwrap();
        assert!(
            r.iters <= 5,
            "linear circuit should converge fast, took {}",
            r.iters
        );
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
             .model myd D (Is=1e-15 N=1)\n.op\n",
        )
        .unwrap();
        let r = dc_op_nr(&net).unwrap();
        // The pn-junction equation gives V(b) ≈ 18·V_T·ln(I/Is) ≈ 0.7 V,
        // and most of the supply drops across R1.
        let vb = r.node_voltage("b").unwrap();
        assert!(vb > 0.5 && vb < 1.0, "V(b)={vb:.4} should be ~0.7V");
        assert!(
            r.iters < 50,
            "convergence took {} iters with pnjlim on",
            r.iters
        );
    }

    #[test]
    fn current_source_biased_diode() {
        let net =
            parse_spice("* Diode bias\nIb 0 b 1m\nD1 b 0 myd\n.model myd D (Is=1e-14 N=1)\n.op\n")
                .unwrap();
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
             .model myd D (Is=1e-14 N=1)\n.op\n",
        )
        .unwrap();
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
        let net =
            parse_spice("* divider\nV1 in 0 DC 2.0\nR1 in out 1k\nR2 out 0 1k\n.op\n").unwrap();
        let r = dc_op_nr(&net).unwrap();
        let mut buf = Vec::new();
        r.write_csv(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("analysis,"), "header: {s}");
        assert!(s.contains("V(out)"), "should contain V(out): {s}");
        assert!(s.contains("dc_op"), "should have dc_op row: {s}");
        assert!(
            s.contains("1.000000e0") || s.contains("1.000000e+0"),
            "V(out)≈1V missing: {s}"
        );
    }

    #[test]
    fn write_nutmeg_dc_op() {
        let net =
            parse_spice("* divider\nV1 in 0 DC 2.0\nR1 in out 1k\nR2 out 0 1k\n.op\n").unwrap();
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
             .model myd D (Is=1e-14 N=1)\n.op\n",
        )
        .unwrap();
        let opts = SimOptions {
            reltol: 1e-10,
            vntol: 1e-12,
            ..Default::default()
        };
        let r = dc_op_nr_opts(&net, &opts).unwrap();
        let vb = r.node_voltage("b").unwrap();
        assert!(vb > 0.5 && vb < 0.8);
    }
}
