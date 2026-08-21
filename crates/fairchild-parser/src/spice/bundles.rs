use crate::{Element, ParseError};

pub(super) fn scan_bidirectional(main_lines: &[(usize, String)]) -> bool {
    let mut bidir = false;
    for (_, line) in main_lines {
        let trimmed = line.trim();
        let lc = trimmed.to_lowercase();
        if !lc.starts_with(".options") && !lc.starts_with(".option") {
            continue;
        }
        for tok in trimmed.split_whitespace().skip(1) {
            let (k, v) = if let Some((k, v)) = tok.split_once('=') {
                (k.to_lowercase(), v.to_string())
            } else {
                (tok.to_lowercase(), String::new())
            };
            if k == "enable_bidirectional"
                || k == "bidirectional"
                || k == "bidirectional_propagation"
            {
                bidir = matches!(v.to_lowercase().as_str(), "" | "1" | "true" | "yes" | "on");
            }
        }
    }
    bidir
}

/// How a device handles WDM bundles, from the *parser's* perspective.  Drives
/// `expand_bundle_ports`: whether to replicate the X-element per channel or
/// flatten the bundle into one instance.
///
/// This is the single source of truth for WDM dispatch.  When you add a new
/// native photonic device that should be bundle-aware (shared electrical
/// across channels), add its model name to `bundle_arity_for`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum BundleArity {
    /// Single-channel device: has no WDM semantics of its own.  Connecting it
    /// to a 1-channel bundle expands to that channel's wires; connecting it to
    /// a multi-channel bundle is an error.
    ///
    /// The parser used to silently replicate such an instance into N parallel
    /// copies, one per channel.  That is wrong whenever any of the device's
    /// ports is electrical, which for a subcircuit is most of the time: eight
    /// copies of an MRM PCell meant eight PN junctions and eight heaters
    /// wired to the same two electrical nodes, drawing 8x the current, with no
    /// diagnostic.  Nothing in the tree needs replication — a device that
    /// should carry a WDM bundle is `Aware` and handles all N channels in one
    /// instance — so the behaviour is gone and the case is now reported.
    Scalar,
    /// Bundle-aware — parser flattens *every* referenced bundle into its
    /// underlying wires and emits a SINGLE X-element with the combined
    /// terminal vector.  The device's `setup_instance` derives the channel
    /// count from `terminals.len()`.
    ///
    /// All bundles referenced by one such instance must agree on their channel
    /// count — optical AND electrical alike, since a device that takes one
    /// control wire per WDM channel is only well-defined when the two widths
    /// match.  This is a *per-instance* rule, not a global one: a netlist may
    /// freely mix a 4-channel optical bus with an 8-wire electrical bus as long
    /// as no single `Aware` instance straddles both.  Enforced in
    /// [`expand_bundle_ports`], which is the only layer that knows each port's
    /// declared width — by the time the device sees a flat terminal vector the
    /// grouping is no longer recoverable.
    Aware,
    /// Bundle bridge (e.g. `fc_mux`, `fc_demux`) — like `Aware`, but skips
    /// the matching-channel-count check (one side intentionally has N
    /// channels while the other has N single-channel bundles).
    Bridge,
}

/// What an instance would look like under each dispatch, so an oracle can
/// decide by shape rather than by name.
///
/// `flattened` and `single` are terminal counts: what the X-element would carry
/// if every referenced bundle were expanded to all its wires, and what it would
/// carry if each contributed one channel.  A model with a fixed terminal count
/// — anything loaded from an OSDI descriptor — is placed by comparing the two
/// against `num_terminals`, which is the rule a `.subckt` instance already
/// follows.
#[derive(Clone, Copy, Debug)]
pub struct ArityQuery<'a> {
    pub model_name: &'a str,
    pub flattened: usize,
    pub single: usize,
}

/// Who decides a model's WDM dispatch.
///
/// The parser cannot ask the device registry directly — `fairchild-core`
/// depends on `fairchild-parser`, so the dependency runs the wrong way — which
/// is why [`bundle_arity_for`] existed as a hand-maintained copy of what the
/// registry knows.  Two lists of one fact drifted, repeatedly: five tier names
/// were missing from it, and a `.model`-card-named instance is looked up under
/// the *card's* name, so no card-based device could ever be found there at all.
///
/// Accepting a policy instead of owning one inverts that without inverting the
/// crate dependency.  `None` means "no opinion" and falls back to the static
/// table, so a parse with no registry behaves exactly as before.
pub trait ArityOracle {
    fn arity(&self, q: &ArityQuery) -> Option<BundleArity>;
}

/// The historic hard-coded table, as an oracle.  What `parse_spice` uses when
/// no registry is available: parser unit tests, and any tool that only wants a
/// syntax tree.
pub struct StaticArity;

impl ArityOracle for StaticArity {
    fn arity(&self, q: &ArityQuery) -> Option<BundleArity> {
        Some(bundle_arity_for(q.model_name))
    }
}

/// Maximally permissive, for the first of two passes.
///
/// Pass one exists only to harvest `.model` cards and model-file paths so the
/// registry can be built; its expansion is thrown away.  It must therefore not
/// fail on the very case the second pass is there to get right — a card-named
/// device on a wide bundle.  `Bridge` flattens and skips the channel-count
/// agreement check, so almost any deck survives to pass two, where the real
/// oracle produces the real answer and any real error.
pub struct PermissiveArity;

impl ArityOracle for PermissiveArity {
    fn arity(&self, q: &ArityQuery) -> Option<BundleArity> {
        match bundle_arity_for(q.model_name) {
            BundleArity::Scalar => Some(BundleArity::Bridge),
            known => Some(known),
        }
    }
}

/// Return the WDM dispatch policy for a model name.
///
/// No longer the authority — see [`ArityOracle`].  This is the fallback for a
/// parse with no registry behind it, and it can only ever recognise a name
/// written literally on an X-line.
pub fn bundle_arity_for(model_name: &str) -> BundleArity {
    match model_name.to_lowercase().as_str() {
        // Bundle bridges — N single-channel ports ↔ one N-channel bus.
        "fc_mux" | "fc_demux" => BundleArity::Bridge,
        // Bundle-aware — every photonic device is bundle-aware.  Pure-
        // optical devices (waveguide, splitter, dcoupler, grating coupler)
        // run independent per-channel propagation; devices with electrical
        // state (pn_ps, thermal_ps, photodetector) share one physical
        // electrical interface across all N channels.  WDM is the rule,
        // not the exception.
        "fc_waveguide" | "fc_splitter" | "fc_dcoupler" | "fc_grating_coupler" | "fc_pn_ps"
        | "fc_pn_ps_cap" | "fc_pn_th_ps" | "fc_thermal_ps" | "fc_thermal_ps_rc" | "fc_mzm"
        | "fc_photodetector" | "fc_circulator" | "fc_optical_2x2" | "fc_awgr" | "fc_facet"
        // The tier names. Each is the same device as the family name above with
        // a LEVEL selected, and `docs/photonic-models.md` documents both
        // spellings as first-class — `fc_pn_ps_cap` and `fc_thermal_ps_rc` were
        // already here, so the list always meant to carry them and these five
        // were simply missed. The omission was a plain disagreement between two
        // lists: `fc_pn_ps_full` was refused on a WDM bus that `fc_pn_ps LEVEL=4`
        // — the identical device — accepted. Guarded by
        // `every_registered_photonic_model_declares_its_arity` in fairchild-core.
        //
        // NOTE this whole match only ever sees a name the user wrote on an
        // X-line, so it cannot help a `.model`-card-named device: the lookup key
        // is the card's name, not its kind. That is #52, and it is why a
        // card-named `fc_awgr` (the only route to table mode) is refused on the
        // very bundles it exists to route.
        | "fc_pn_ps_inj" | "fc_pn_ps_full" | "fc_pn_th_ps_cap" | "fc_pn_th_ps_inj"
        | "fc_pn_th_ps_full" => BundleArity::Aware,
        // `fc_cw_laser` / `fc_driven_laser` deliberately stay Scalar — a single laser source
        // produces one wavelength.  Combine multiple lasers via `fc_mux` for
        // WDM operation.  All non-photonic devices (R, C, L, D, MOSFETs)
        // also Scalar.
        _ => BundleArity::Scalar,
    }
}

/// Expand any tokens in an XOsdi element's net list that match a declared
/// `.optical_port` / `.electrical_port`.  Always returns exactly one element:
/// `Aware` / `Bridge` devices take the whole flattened bus, `Scalar` devices
/// take a single channel (and are rejected on a wider one).  Non-XOsdi
/// elements are passed through unchanged.
///
/// This is also where the per-instance channel-count agreement is enforced for
/// `Aware` devices — see [`BundleArity::Aware`].
/// `subckt_ports` is the declared port count of the `.subckt` this element
/// names, when it names one — it selects flatten-vs-replicate semantics for
/// subcircuit instances (see the arithmetic below).
pub(super) fn expand_bundle_ports(
    el: Element,
    ports: &[crate::BundlePort],
    subckt_ports: Option<usize>,
    lineno: usize,
    oracle: &dyn ArityOracle,
) -> Result<Vec<Element>, ParseError> {
    let Element::XOsdi {
        name,
        nets,
        model_name,
        params,
    } = el
    else {
        return Ok(vec![el]);
    };
    // Build a per-token (matching_port_index | None) map; collect each
    // referenced port's channel count to detect mismatches.
    let mut port_refs: Vec<Option<usize>> = Vec::with_capacity(nets.len());
    let mut channel_counts: Vec<usize> = Vec::new();
    for net in &nets {
        match ports.iter().position(|p| p.name == *net) {
            Some(i) => {
                port_refs.push(Some(i));
                channel_counts.push(ports[i].channels);
            }
            None => port_refs.push(None),
        }
    }
    if port_refs.iter().all(|r| r.is_none()) {
        // No port references — return the element unchanged.
        return Ok(vec![Element::XOsdi {
            name,
            nets,
            model_name,
            params,
        }]);
    }
    // Dispatch by WDM policy.  See `BundleArity` for semantics.  The oracle
    // decides; the static table is only the fallback when nothing better is
    // available (see `ArityOracle`).  Both candidate terminal counts go with
    // the query so an oracle holding a fixed-arity descriptor can place the
    // instance by shape instead of by name.
    let widths = |per_port: &dyn Fn(&crate::BundlePort) -> usize| -> usize {
        nets.iter()
            .zip(port_refs.iter())
            .map(|(_, r)| match r {
                Some(i) => per_port(&ports[*i]),
                None => 1,
            })
            .sum()
    };
    let query = ArityQuery {
        model_name: &model_name,
        flattened: widths(&|p| p.channels * p.wires_per_channel()),
        single: widths(&|p| p.wires_per_channel()),
    };
    let arity = oracle
        .arity(&query)
        .unwrap_or_else(|| bundle_arity_for(&model_name));

    // Helper: flatten every referenced bundle into its underlying wires, in
    // declaration order, into one combined net vector.
    let flatten = || -> Vec<String> {
        let mut flat: Vec<String> = Vec::with_capacity(nets.len() * 3);
        for (net, port_ref) in nets.iter().zip(port_refs.iter()) {
            if let Some(port_idx) = port_ref {
                let port = &ports[*port_idx];
                for ch in 0..port.channels {
                    flat.extend(port.wires_for_channel(ch));
                }
            } else {
                flat.push(net.clone());
            }
        }
        flat
    };

    if arity == BundleArity::Bridge {
        return Ok(vec![Element::XOsdi {
            name,
            nets: flatten(),
            model_name,
            params,
        }]);
    }
    // Validate consistent channel count when more than one port is involved.
    // A bundle-aware device serves all its channels from one instance, so a
    // 4-channel optical bus alongside a 2-wire control bus has no meaning —
    // and by the time the device sees a flat terminal vector, the grouping is
    // unrecoverable. Reject it here, naming each port so the fix is obvious.
    let max_n = channel_counts.iter().copied().max().unwrap_or(1);
    if channel_counts.iter().any(|&n| n != max_n) {
        let detail = port_refs
            .iter()
            .flatten()
            .map(|&i| {
                let p = &ports[i];
                let kind = if p.is_optical() {
                    "optical"
                } else {
                    "electrical"
                };
                format!("{}({kind}, {} ch)", p.name, p.channels)
            })
            .collect::<Vec<_>>()
            .join(", ");
        return Err(ParseError::Syntax {
            line: lineno,
            msg: format!(
                "X{name}: every bundle port on one instance must declare the same \
                 channel count; got {detail}"
            ),
        });
    }
    if arity == BundleArity::Aware {
        return Ok(vec![Element::XOsdi {
            name,
            nets: flatten(),
            model_name,
            params,
        }]);
    }
    // BundleArity::Scalar.  A `.subckt` instance arrives here too — a subckt
    // name is never in the bundle-aware table — but a subckt can legitimately
    // be written to carry a whole bus, and it says so by declaring that many
    // ports.  So: accept it as one instance over the flattened bus when the
    // port count matches, and otherwise require a single-channel connection.
    let flat_width = flatten().len();
    let per_ch_width: usize = port_refs
        .iter()
        .map(|r| match r {
            Some(i) => ports[*i].wires_per_channel(),
            None => 1,
        })
        .sum();
    if let Some(n_ports) = subckt_ports {
        if n_ports == flat_width {
            return Ok(vec![Element::XOsdi {
                name,
                nets: flatten(),
                model_name,
                params,
            }]);
        }
        if n_ports != per_ch_width {
            return Err(ParseError::Syntax {
                line: lineno,
                msg: format!(
                    "X{name}: subckt '{model_name}' declares {n_ports} ports, but this \
                     instance's bundle references expand to {flat_width} wires (the whole \
                     {max_n}-channel bus). Give the subckt {flat_width} ports to carry the \
                     bus, or connect it to a single channel."
                ),
            });
        }
    }
    if max_n > 1 {
        let detail = port_refs
            .iter()
            .flatten()
            .map(|&i| format!("{}({} ch)", ports[i].name, ports[i].channels))
            .collect::<Vec<_>>()
            .join(", ");
        return Err(ParseError::Syntax {
            line: lineno,
            msg: format!(
                "X{name}: '{model_name}' has no WDM semantics but is connected to a \
                 multi-channel bundle: {detail}. One instance cannot serve {max_n} \
                 channels, and replicating it would duplicate any electrical port \
                 {max_n}× onto the same nodes. Connect a single channel, use \
                 fc_demux/fc_mux to split and recombine, or — if this model really is \
                 bundle-aware — give it {} terminals so the whole bus fits, since a \
                 model's terminal count is what decides this.",
                query.flattened
            ),
        });
    }
    // Single channel: expand the bundle references to that channel's wires.
    let mut expanded_nets: Vec<String> = Vec::new();
    for (net, port_ref) in nets.iter().zip(port_refs.iter()) {
        if let Some(port_idx) = port_ref {
            expanded_nets.extend(ports[*port_idx].wires_for_channel(0));
        } else {
            expanded_nets.push(net.clone());
        }
    }
    Ok(vec![Element::XOsdi {
        name,
        nets: expanded_nets,
        model_name,
        params,
    }])
}
