use super::common::{canon_node, expand_bus_vectors, parse_value};
use super::element::parse_element;
use crate::{Element, ParseError, Waveform};

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
/// `expand_optical_ports`: whether to replicate the X-element per channel or
/// flatten the bundle into one instance.
///
/// This is the single source of truth for WDM dispatch.  When you add a new
/// native photonic device that should be bundle-aware (shared electrical
/// across channels), add its model name to `bundle_arity_for`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum BundleArity {
    /// Pure-electrical or per-channel-only device — parser replicates the
    /// X-element into N parallel instances when bundles are connected (one
    /// per channel).  Each replica gets the channel's underlying 3 wires.
    Scalar,
    /// Bundle-aware — parser flattens *every* referenced bundle into its
    /// underlying wires and emits a SINGLE X-element with the combined
    /// terminal vector.  The device's `setup_instance` derives the channel
    /// count from `terminals.len()`.  All referenced bundles must agree on
    /// their channel count.
    Aware,
    /// Bundle bridge (e.g. `fc_mux`, `fc_demux`) — like `Aware`, but skips
    /// the matching-channel-count check (one side intentionally has N
    /// channels while the other has N single-channel bundles).
    Bridge,
}

/// Return the WDM dispatch policy for a model name.  Centralises the
/// hard-coded list of bundle-aware native photonics so the parser doesn't
/// scatter `if model == "fc_pn_ps" || ...` chains.
///
/// Future extension: when external (PDK) devices need bundle-awareness, this
/// can be backed by a registration hook.  Today every bundle-aware model is
/// shipped in `fairchild-core`, so a static `match` suffices.
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
        | "fc_photodetector" | "fc_circulator" => BundleArity::Aware,
        // `fc_cw_laser` deliberately stays Scalar — a single laser source
        // produces one wavelength.  Combine multiple lasers via `fc_mux` for
        // WDM operation.  All non-photonic devices (R, C, L, D, MOSFETs)
        // also Scalar.
        _ => BundleArity::Scalar,
    }
}

/// Expand any tokens in an XOsdi element's net list that match a declared
/// `.optical_port`.  Returns one XOsdi per channel (most ports have channels
/// = 1, so most lines return a single element).  Non-XOsdi elements are
/// passed through unchanged.
pub(super) fn expand_optical_ports(
    el: Element,
    ports: &[crate::OpticalPort],
    lineno: usize,
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
    // Dispatch by WDM policy.  See `BundleArity` for semantics.
    let arity = bundle_arity_for(&model_name);

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
    let max_n = channel_counts.iter().copied().max().unwrap_or(1);
    if channel_counts.iter().any(|&n| n != max_n) {
        return Err(ParseError::Syntax {
            line: lineno,
            msg: format!(
                "X{name}: bundle ports must have matching channel counts, got {:?}",
                channel_counts
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
    // BundleArity::Scalar — replicate per channel.
    let mut out = Vec::with_capacity(max_n);
    for ch in 0..max_n {
        let mut expanded_nets: Vec<String> = Vec::new();
        for (net, port_ref) in nets.iter().zip(port_refs.iter()) {
            if let Some(port_idx) = port_ref {
                expanded_nets.extend(ports[*port_idx].wires_for_channel(ch));
            } else {
                expanded_nets.push(net.clone());
            }
        }
        let new_name = if max_n > 1 {
            format!("{name}_ch{ch}")
        } else {
            name.clone()
        };
        out.push(Element::XOsdi {
            name: new_name,
            nets: expanded_nets,
            model_name: model_name.clone(),
            params: params.clone(),
        });
    }
    Ok(out)
}
