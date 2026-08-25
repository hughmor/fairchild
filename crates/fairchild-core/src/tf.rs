//! `.tf` — small-signal transfer function about the DC operating point.
//!
//! Three numbers: the gain from an input source to an output, the resistance
//! the source sees looking into the circuit, and the resistance the output port
//! presents looking back.  All three are derivatives of the operating point
//! with respect to a source level, so all three are [`crate::adjoint`] calls —
//! this module contributes no linear algebra of its own.
//!
//! ## Why this is not its own solver
//!
//! The obvious implementation stamps the DC Jacobian, factorises it, and hits
//! it with three unit excitations.  That is what the machinery in `adjoint`
//! already does, and a second copy would be a second opinion about which
//! Jacobian is the right one — the operating-point Jacobian has a repair pass
//! for frozen columns (see [`crate::device::Device::frozen_jacobian_columns`])
//! that a fresh implementation would not know to apply, so its answers would
//! differ from `.sens`'s on exactly the circuits where the difference matters.
//!
//! ## Output resistance, and the probe source
//!
//! Output resistance is `∂V_port/∂I_injected`, and `I_injected` is not a
//! parameter of the user's deck — so this adds one: a zero-valued current
//! source across the output port, on a *copy* of the netlist.  A zero-valued
//! current source stamps nothing, adds no node and adds no MNA row, so the
//! operating point it is added to is the operating point without it; the
//! adjoint's central difference on its level then reads back exactly the
//! injection derivative.  This is a real element going through the real
//! stamping path, which is the point: it cannot disagree with how the solver
//! treats a current source, because it *is* how the solver treats one.
//!
//! ## Signs
//!
//! Two negations below look arbitrary and are not.  With `f = A·x − b`:
//!
//! * A current source `I p n` stamps `b[p] −= I`, `b[n] += I` — SPICE's
//!   convention that positive current leaves `p` *inside* the source.  Driving
//!   a resistor `R` that way gives `V(p,n) = −I·R`, so
//!   `R_out = −∂V_port/∂I_probe`.
//! * A voltage source's branch row holds the current leaving its `+` node
//!   through the source, which is SPICE's `I(V)` and is *negative* for a source
//!   delivering power.  Driving `R` gives `I(V) = −E/R`, so
//!   `R_in = −1/(∂I(V_in)/∂E)`.
//!
//! Both are pinned against ngspice in `tests/ngspice/ngspice_tf_pz_golden.rs`
//! rather than argued from here, because a sign convention agreed with only
//! itself is exactly the kind of self-consistent wrong answer this codebase has
//! shipped before.

use fairchild_parser::{Element, Netlist, OutVar, Waveform};

use crate::adjoint::{dc_sensitivity, element_name, Output, ParamRef};
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::options::SimOptions;

/// What `.tf` reports.
#[derive(Debug, Clone)]
pub struct TfResult {
    /// The output's value at the operating point, so the gain has something to
    /// be a gain *about*.
    pub out_value: f64,
    /// `∂out/∂in`.  Dimensionless for a voltage-in/voltage-out pair; Ω, S or
    /// dimensionless otherwise, following the two ports' units.
    pub gain: f64,
    /// Resistance seen looking into the circuit from the input source's
    /// terminals, in Ω.  Infinite when the input port draws no current — which
    /// is reported as `f64::INFINITY`, not as a large number.
    pub r_in: f64,
    /// Resistance seen looking back into the circuit from the output port, in Ω.
    pub r_out: f64,
}

impl TfResult {
    /// The report as `(label, value)` pairs, in the order ngspice prints them.
    ///
    /// One list, so the CSV and the rawfile cannot come to disagree about which
    /// three numbers a `.tf` is.
    pub fn rows(&self, out_label: &str, in_label: &str) -> Vec<(String, f64)> {
        vec![
            ("transfer_function".to_string(), self.gain),
            (format!("output_impedance_at_{out_label}"), self.r_out),
            (format!("{in_label}#input_impedance"), self.r_in),
            (out_label.to_string(), self.out_value),
        ]
    }

    pub fn write_csv<W: std::io::Write>(
        &self,
        mut w: W,
        out_label: &str,
        in_label: &str,
    ) -> std::io::Result<()> {
        writeln!(w, "quantity,value")?;
        for (name, v) in self.rows(out_label, in_label) {
            writeln!(w, "{name},{v:.6e}")?;
        }
        Ok(())
    }

    pub fn write_nutmeg<W: std::io::Write>(
        &self,
        mut w: W,
        title: &str,
        out_label: &str,
        in_label: &str,
    ) -> std::io::Result<()> {
        let rows = self.rows(out_label, in_label);
        writeln!(w, "Title: {title}")?;
        writeln!(w, "Plotname: Transfer Function")?;
        writeln!(w, "Flags: real")?;
        writeln!(w, "No. Variables: {}", rows.len())?;
        writeln!(w, "No. Points: 1")?;
        writeln!(w, "Variables:")?;
        for (i, (name, _)) in rows.iter().enumerate() {
            writeln!(w, "\t{i}\t{name}\tnotype")?;
        }
        writeln!(w, "Values:")?;
        for (i, (_, v)) in rows.iter().enumerate() {
            if i == 0 {
                writeln!(w, " 0\t{v:.6e}")?;
            } else {
                writeln!(w, "\t{v:.6e}")?;
            }
        }
        Ok(())
    }
}

/// Name of the input source and which kind it is, resolved from the deck.
enum InputPort {
    /// A voltage source: the input current is its branch current.
    Voltage { name: String },
    /// A current source: the input voltage is across its terminals.
    Current {
        name: String,
        pos: String,
        neg: String,
    },
}

/// `.tf <out> <input_src>` at the DC operating point.
pub fn transfer_function(
    netlist: &Netlist,
    registry: &DeviceRegistry,
    opts: &SimOptions,
    out: &OutVar,
    input_src: &str,
) -> Result<TfResult, SimError> {
    let port = find_input(netlist, input_src)?;

    // The output port, and the probe that measures its resistance.  A current
    // output is already a port — the voltage source it flows through — so it
    // needs no probe: sweeping that source's own level gives the same
    // derivative.
    let mut work = netlist.clone();
    let probe = match out {
        OutVar::NodeVoltage { pos, neg } => Some(add_probe(&mut work, pos, neg)),
        OutVar::BranchCurrent(_) => None,
    };

    // One adjoint call: two outputs (the user's, and whatever reads the input
    // port's response) against two parameters (the input source's level, and
    // the probe's).  That is two transposed solves and four re-stamps sharing
    // one factorisation of one operating point.
    let out_adj = adjoint_output(out);
    let in_adj = match &port {
        InputPort::Voltage { name } => Output::BranchCurrent(name.clone()),
        InputPort::Current { pos, neg, .. } => Output::NodeVoltageDiff {
            pos: pos.clone(),
            neg: neg.clone(),
        },
    };
    let outputs = vec![out_adj, in_adj];

    let in_name = match &port {
        InputPort::Voltage { name } | InputPort::Current { name, .. } => name.clone(),
    };
    let mut params = vec![ParamRef {
        element: in_name.clone(),
        param: "dc".into(),
        nominal: None,
        step: None,
    }];
    // The second parameter is whichever knob the output resistance falls out
    // of: the probe's level for a voltage output, the output source's own level
    // for a current one.  Both land in `params[1]`, except when the output
    // source *is* the input source (`.tf i(vin) vin`), where the derivative
    // wanted is already `params[0]` and asking twice would run the same
    // perturbation a second time.
    let r_out_param = match (&probe, out) {
        (Some(p), _) => {
            params.push(ParamRef {
                element: p.clone(),
                param: "dc".into(),
                // The probe is zero-valued, and `∛ε·|0|` is zero — a step that
                // perturbs nothing and reports the parameter unreached.  1 mA
                // is small enough to stay inside the linearisation and large
                // enough to move the residual well clear of `abstol`.
                nominal: Some(0.0),
                step: Some(1e-3),
            });
            1
        }
        (None, OutVar::BranchCurrent(name)) if name.to_lowercase() == in_name => 0,
        (None, OutVar::BranchCurrent(name)) => {
            params.push(ParamRef {
                element: name.clone(),
                param: "dc".into(),
                nominal: None,
                step: None,
            });
            1
        }
        (None, _) => unreachable!("a node-voltage output always gets a probe"),
    };

    let sens = dc_sensitivity(&work, registry, opts, &outputs, &params)?;

    if !sens.reached[0] {
        return Err(SimError::ParameterError(format!(
            "'{in_name}' does not respond to a level change, so there is no transfer \
             function through it. A `.tf` input must be an independent V or I source \
             carrying a DC level (e.g. `{in_name} in 0 DC 0`)"
        )));
    }

    let gain = sens.grad[0][0];
    // ∂(input port response)/∂(source level), from which the input resistance
    // follows — reciprocally for a voltage drive, directly for a current one.
    let d_in = sens.grad[1][0];
    let r_in = match &port {
        InputPort::Voltage { .. } => {
            if d_in == 0.0 {
                f64::INFINITY
            } else {
                -1.0 / d_in
            }
        }
        InputPort::Current { .. } => -d_in,
    };

    let d_out = sens.grad[0][r_out_param];
    let r_out = if probe.is_some() {
        // `∂V_port/∂I_probe`, negated for the injection convention above.
        -d_out
    } else if d_out == 0.0 {
        // A current output's port is its own voltage source: `∂I/∂E` there is
        // the port conductance, so its reciprocal is the resistance looking
        // back in — the same reciprocal, and the same sign, as `r_in`.
        f64::INFINITY
    } else {
        -1.0 / d_out
    };

    Ok(TfResult {
        out_value: sens.values[0],
        gain,
        r_in,
        r_out,
    })
}

/// The deck's `OutVar` as an adjoint `Output`.
fn adjoint_output(out: &OutVar) -> Output {
    match out {
        OutVar::NodeVoltage { pos, neg } => Output::NodeVoltageDiff {
            pos: pos.clone(),
            neg: neg.clone(),
        },
        OutVar::BranchCurrent(name) => Output::BranchCurrent(name.clone()),
    }
}

/// Locate the named input source, or say what was found instead.
fn find_input(netlist: &Netlist, name: &str) -> Result<InputPort, SimError> {
    let want = name.to_lowercase();
    for el in &netlist.elements {
        match el {
            Element::VoltageSource { name: n, .. } if n.to_lowercase() == want => {
                return Ok(InputPort::Voltage { name: want })
            }
            Element::CurrentSource {
                name: n, pos, neg, ..
            } if n.to_lowercase() == want => {
                return Ok(InputPort::Current {
                    name: want,
                    pos: pos.clone(),
                    neg: neg.clone(),
                })
            }
            _ => {}
        }
    }
    Err(SimError::ParameterError(format!(
        "'{name}' is not an independent source in this deck, so `.tf` has nothing to \
         take a transfer function from. Name a V or I source; a transfer from a \
         dependent source or a device terminal is not what `.tf` means"
    )))
}

/// Add a zero-valued current source across `(pos, neg)` and return its name.
///
/// Named so it cannot collide with anything a deck can spell, then uniquified
/// anyway — a probe that silently landed on top of a user's element would
/// re-target the whole analysis without saying so.
fn add_probe(work: &mut Netlist, pos: &str, neg: &str) -> String {
    let mut name = "i@tf_probe".to_string();
    while work
        .elements
        .iter()
        .any(|el| element_name(el).is_some_and(|n| n.to_lowercase() == name))
    {
        name.push('_');
    }
    work.elements.push(Element::CurrentSource {
        name: name.clone(),
        pos: pos.to_string(),
        neg: neg.to_string(),
        waveform: Waveform::Dc(0.0),
        ac: None,
    });
    name
}
