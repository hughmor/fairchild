//! Resolving λ before the solve, so it need not be solved for.
//!
//! A wavelength is a label. Measured across every photonic deck in the tree,
//! each λ row reads exactly a source's wavelength or the band-centre value an
//! undriven wire bootstraps to, and never anything computed
//! (`tests/lambda_is_a_label.rs` pins that). The solver has been doing label
//! propagation through a linear subsystem embedded in the matrix.
//!
//! Leaving it there cost more than the rows it occupied — 864 of 2840 unknowns
//! on the giona 8-neuron RNN, 30 % of the matrix. `vntol = 1e-6` is meaningless
//! against a 1.55e-6 quantity, so λ needed its own tolerance class. A
//! volt-scaled trust region once scaled λ from 1.55 µm to 1e-19 m and destroyed
//! every optical phase in the circuit, so λ had to be excluded from that too.
//! `LambdaSelect` latched which input to mirror from values that were still
//! settling. And no device could differentiate against λ — propagation phase is
//! thousands of radians, so `∂φ/∂λ = φ/λ` is of order 1e9 per metre — which
//! `adjoint.rs` encoded by treating every λ column as frozen. All four are
//! gone with the rows.
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
//! ([`crate::device::Device::lambda_routing`]) and where one originates
//! ([`crate::device::Device::lambda_emitted`]). Same shape as bundle arity: knowledge that used
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

use crate::device::SimContext;
use crate::device_registry::{DeviceRegistry, ParamSet};
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

    /// Every net this map speaks for. These are the nets that must *not* become
    /// MNA rows: their value is decided here instead of solved for.
    pub fn nets(&self) -> impl Iterator<Item = &str> {
        self.by_net.keys().map(String::as_str)
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
/// An unreached net takes `ctx.lambda_center_m`, matching the bootstrap an
/// undriven λ wire gets today.
///
/// # Why this builds its own devices
///
/// Only a device knows where its labels go, and the answer depends on its
/// bundle width — so the declarations have to be *asked of instances*, not
/// looked up in a table keyed by model name. But the row layout depends on the
/// answer, so this cannot wait for the real instances: they are bound to the
/// matrix this pass decides the shape of. So each `X` element gets a throwaway
/// instance with every terminal unbound. `lambda_routing`/`lambda_emitted` are
/// pure functions of a device's shape and parameters, neither of which a node
/// index changes, so the throwaway answers exactly as the real one would.
///
/// An unknown model is skipped rather than reported: `build_devices` raises that
/// error, with the element name, moments later.
pub fn resolve(netlist: &Netlist, ctx: &SimContext, registry: &DeviceRegistry) -> LambdaMap {
    let mut edges: Vec<(String, String)> = Vec::new();
    let mut seeds: Vec<(String, f64)> = Vec::new();
    // Every λ net in the deck, not only the ones some routing mentions: a dark
    // port is still a λ net and must still get an answer — an AWGR's unused
    // input, a ring's dark add. `lambda_terminals` is what makes it total.
    let mut lambda_nets: Vec<String> = Vec::new();

    // Models that declared nothing at all, paired with the nets they touch, so
    // the pass can say afterwards whether any of them sits on a lit path.
    let mut silent: Vec<(&str, &str, &[String])> = Vec::new();
    for el in &netlist.elements {
        let Element::XOsdi {
            name,
            nets,
            model_name,
            params,
        } = el
        else {
            continue;
        };
        let Some(factory) = registry.get(model_name) else {
            continue;
        };
        let unbound = vec![None; nets.len()];
        let dev = factory(&unbound, &ParamSet::new(params), ctx);
        let declared = dev.lambda_terminals();
        if declared.is_empty() {
            silent.push((name, model_name, nets));
        }
        for t in declared {
            if let Some(n) = nets.get(t) {
                lambda_nets.push(n.clone());
            }
        }
        for (from, to) in dev.lambda_routing() {
            if let (Some(a), Some(b)) = (nets.get(from), nets.get(to)) {
                edges.push((a.clone(), b.clone()));
            }
        }
        for (t, wl) in dev.lambda_emitted() {
            if let Some(n) = nets.get(t) {
                seeds.push((n.clone(), wl));
            }
        }
    }

    // A voltage source on a λ net is a deck author naming a wavelength by hand
    // — the idiom every hand-wired bundle in the tree and the test suite uses,
    // and the only way to label a bundle whose light comes from outside the
    // deck. It has to seed resolution: silently substituting the band centre
    // for a wire someone explicitly drove to 1551 nm is exactly the class of
    // wrong answer this codebase refuses to ship. Only nets some device already
    // agreed are λ nets qualify, so an ordinary supply is never mistaken for
    // one.
    let declared: std::collections::HashSet<&str> =
        lambda_nets.iter().map(String::as_str).collect();
    for el in &netlist.elements {
        let Element::VoltageSource {
            pos,
            waveform: fairchild_parser::Waveform::Dc(v),
            ..
        } = el
        else {
            continue;
        };
        if *v > 0.0 && declared.contains(pos.as_str()) {
            seeds.push((pos.clone(), *v));
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
        by_net.insert(n.clone(), ctx.lambda_center_m);
    }

    warn_about_silent_models(&silent, &by_net, &unreached);
    LambdaMap { by_net, unreached }
}

/// Name any model that sits on a path whose wavelength is known and says
/// nothing about what it does with the label.
///
/// A fixed-port Verilog-A optical model is the case: it has `in_wl` / `out_wl`
/// ports and used to carry the tag through the matrix with
/// `OWL(out_wl) <+ OWL(in_wl)`, which worked while λ was an unknown. It is not
/// one any more, so that contribution goes to a node the matrix does not have,
/// and everything downstream resolves to the band centre — a wrong wavelength
/// with no diagnostic. Nothing here can infer which of its ports is an input,
/// so the only honest answer is to say so.
///
/// Once per model name, and only when resolution actually reached one of its
/// nets: a model on a dark branch has nothing to get wrong. A terminator (a
/// photodetector written in Verilog-A) is a false positive, and is accepted as
/// the price of not being silent about the case that is genuinely broken.
fn warn_about_silent_models(
    silent: &[(&str, &str, &[String])],
    by_net: &HashMap<String, f64>,
    unreached: &[String],
) {
    if silent.is_empty() {
        return;
    }
    // Only a model on a LIT path is worth mentioning. One sitting entirely on
    // dark ports routes nothing anyone can observe, and warning about it would
    // train the reader to ignore the message that matters.
    let dark: std::collections::HashSet<&str> = unreached.iter().map(String::as_str).collect();
    let mut said: Vec<&str> = Vec::new();
    for (inst, model, nets) in silent {
        if said.contains(model) {
            continue;
        }
        let lit: Vec<&str> = nets
            .iter()
            .filter(|n| by_net.contains_key(*n) && !dark.contains(n.as_str()))
            .map(String::as_str)
            .collect();
        if lit.is_empty() {
            continue;
        }
        said.push(model);
        crate::warn_user!(
            "X{inst} ('{model}') is wired to {} — an optical net whose wavelength is \
             known — but declares no wavelength routing, so nothing downstream of it \
             inherits a label and will be evaluated at the band centre instead of the \
             colour actually present. A Verilog-A model written against the bundle-port \
             dialect (`optical_bundle`) declares this for you; a fixed-port one cannot, \
             because nothing can tell which of its ports is an input. Either port it to \
             the dialect, or set the wavelength on the instances downstream.",
            lit.join(", ")
        );
    }
}
