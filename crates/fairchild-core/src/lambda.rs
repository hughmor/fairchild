//! Resolving λ before the solve, so it need not be solved for.
//!
//! A wavelength is a label. Measured across every photonic deck in the tree,
//! each λ row reads exactly a source's wavelength or the band-centre value an
//! undriven wire bootstraps to, and never anything computed
//! (`tests/lambda_is_a_label.rs` pins that). The solver has been doing label
//! propagation through a linear subsystem embedded in the matrix.
//!
//! Leaving it there costs more than the rows it occupies — 848 of 2840 unknowns
//! on the giona 8-neuron RNN, 30 % of the matrix. `vntol = 1e-6` is meaningless
//! against a 1.55e-6 quantity, so λ needed its own tolerance class. A
//! volt-scaled trust region once scaled λ from 1.55 µm to 1e-19 m and destroyed
//! every optical phase in the circuit, so λ had to be excluded from that too.
//! `LambdaSelect` latches which input to mirror from values that are still
//! settling. And no device may differentiate against λ — propagation phase is
//! thousands of radians, so `∂φ/∂λ = φ/λ` is of order 1e9 per metre — which
//! `adjoint.rs` already encodes by treating every λ column as frozen.
//!
//! # Why routing is declared, not inferred
//!
//! An earlier attempt read the routing off the assembled Jacobian, whose λ rows
//! do encode it. That cannot survive the change it was meant to enable: the
//! matrix is what is going away. Worse, it was not even where the constraint
//! lives — `V(out_λ) − V(in_λ) = 0` is stamped into a device *branch* row, and
//! the λ node rows carry only KCL over those branches, so replacing a λ node row
//! severs the KCL and leaves nodes held by nothing but `gmin`.
//!
//! So each device declares how a label moves through it
//! ([`Device::lambda_routing`]) and where one originates
//! ([`Device::lambda_emitted`]). Same shape as bundle arity: knowledge that used
//! to be inferred from a structure a device happens to build, moved next to the
//! device that knows it.
//!
//! # What resolution does not decide
//!
//! A source whose wavelength depends on circuit state — laser chirp, thermal
//! drift of an emitter — is out of scope, and is out of scope today too: the
//! objection recorded against chirp is that a drive-dependent λ wire is one
//! "every downstream device would then see move", which is exactly the Jacobian
//! that will not converge. When it is wanted, the mechanism is an outer
//! iteration around this pass.

use std::collections::HashMap;

use crate::device::Device;
use fairchild_parser::{Element, Netlist};

/// Every λ net's wavelength, in metres.
#[derive(Debug, Default, Clone)]
pub struct LambdaMap {
    by_net: HashMap<String, f64>,
    /// Nets that carry a label but no source ever reached, in netlist order.
    unreached: Vec<String>,
}

impl LambdaMap {
    /// Channel λ for a net, if resolution reached it.
    pub fn get(&self, net: &str) -> Option<f64> {
        self.by_net.get(net).copied()
    }

    pub fn len(&self) -> usize {
        self.by_net.len()
    }

    pub fn is_empty(&self) -> bool {
        self.by_net.is_empty()
    }

    /// λ nets no source reached — a dark branch of the circuit. Not an error:
    /// an undriven optical input is legitimate and bootstraps to the band
    /// centre. Exposed so a caller can say which, rather than guessing.
    pub fn unreached(&self) -> &[String] {
        &self.unreached
    }
}

/// The net name each of an element's terminals connects to.
fn terminal_nets(el: &Element) -> Option<&[String]> {
    match el {
        Element::XOsdi { nets, .. } => Some(nets),
        _ => None,
    }
}

/// Resolve every λ net by propagating each source's wavelength along the
/// routing its devices declare.
///
/// Label propagation to a fixpoint, so a cycle terminates: a ring's add/drop
/// loop carries the same λ all the way round, and revisiting a net with the
/// value it already has changes nothing. Where two different wavelengths reach
/// one net the first wins and the disagreement is *not* silently averaged — the
/// net is recorded so a caller can report it, since two sources at different
/// wavelengths on one wire is a deck bug rather than a physical mixture.
///
/// `fallback` is what an unreached net takes — `lambda_center_m`, matching the
/// bootstrap an undriven λ wire gets today.
pub fn resolve(netlist: &Netlist, devices: &[Box<dyn Device>], fallback: f64) -> LambdaMap {
    // Element order and device order agree: `build_devices` walks the elements.
    // Zip rather than index so a mismatch truncates instead of panicking.
    let mut edges: Vec<(String, String)> = Vec::new();
    let mut seeds: Vec<(String, f64)> = Vec::new();
    // Every λ net in the deck, not only the ones some routing mentions. A dark
    // port is still a λ net and must still get an answer — an AWGR's unused
    // input, a ring's dark add. Asking the netlist rather than the declarations
    // is what makes the map total.
    let mut lambda_nets: Vec<String> = netlist
        .optical_nets
        .iter()
        .filter(|n| fairchild_parser::is_lambda_wire(n))
        .cloned()
        .collect();

    for (el, dev) in netlist
        .elements
        .iter()
        .filter(|e| matches!(e, Element::XOsdi { .. }))
        .zip(devices.iter())
    {
        let Some(nets) = terminal_nets(el) else {
            continue;
        };
        for (from, to) in dev.lambda_routing() {
            if let (Some(a), Some(b)) = (nets.get(from), nets.get(to)) {
                edges.push((a.clone(), b.clone()));
                lambda_nets.push(a.clone());
                lambda_nets.push(b.clone());
            }
        }
        for (t, wl) in dev.lambda_emitted() {
            if let Some(n) = nets.get(t) {
                seeds.push((n.clone(), wl));
                lambda_nets.push(n.clone());
            }
        }
    }

    let mut by_net: HashMap<String, f64> = HashMap::new();
    for (net, wl) in seeds {
        by_net.entry(net).or_insert(wl);
    }

    // Fixpoint. Bounded by the number of edges: each round settles at least one
    // more net or stops, so a chain of E devices needs at most E rounds.
    for _ in 0..=edges.len() {
        let mut changed = false;
        for (from, to) in &edges {
            if let Some(&wl) = by_net.get(from) {
                if !by_net.contains_key(to) {
                    by_net.insert(to.clone(), wl);
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }

    lambda_nets.sort_unstable();
    lambda_nets.dedup();
    let unreached: Vec<String> = lambda_nets
        .iter()
        .filter(|n| !by_net.contains_key(*n))
        .cloned()
        .collect();
    for n in &unreached {
        by_net.insert(n.clone(), fallback);
    }

    LambdaMap { by_net, unreached }
}
