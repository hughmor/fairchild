//! `.sens` — DC parameter sensitivity from a deck card.
//!
//! The computation is [`crate::adjoint::dc_sensitivity`] and has been for some
//! time; what this module adds is the two things a *card* needs that a Python
//! caller supplied for itself: which parameters "all of them" means, and a
//! report that says which of them the adjoint could not reach.
//!
//! ## This is not ngspice's `.sens`
//!
//! ngspice's `.sens` perturbs each parameter and re-solves — one full
//! nonlinear solve per parameter, differencing a result that is only converged
//! to `reltol`.  The adjoint gets every parameter from one transposed solve,
//! and differences the *residual*, which is an explicit function evaluated to
//! machine precision.  The numbers here are the same quantity computed better,
//! not an approximation of ngspice's; see [`crate::adjoint`] for why the error
//! is ~1e-10 relative rather than ~1e-3.
//!
//! The one visible consequence: a parameter the adjoint cannot reach is
//! reported as unreached rather than as a zero.  A wrong zero and a real
//! insensitivity look identical in a table of numbers, and this codebase has
//! shipped that mistake before.

use fairchild_parser::{Element, Netlist, OutVar, ParamName};

use crate::adjoint::{dc_sensitivity, Output, ParamRef, Sensitivities};
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::options::SimOptions;

/// One row of a `.sens` report.
#[derive(Debug, Clone)]
pub struct SensRow {
    /// `element.param`, as the report names it.
    pub name: String,
    /// The parameter's value at the operating point.
    pub nominal: f64,
    /// `∂out/∂p`, in output-units per parameter-unit.
    pub sensitivity: f64,
    /// `∂out/∂p · p` — the change in the output for a 100 % change in the
    /// parameter, which is what makes two parameters of different units
    /// comparable.  ngspice calls this the normalised sensitivity.
    pub normalised: f64,
    /// False when perturbing the parameter did not move the residual at all.
    /// `sensitivity` is then a placeholder, not a measurement.
    pub reached: bool,
    /// Relative disagreement between the two finite-difference step sizes — a
    /// conservative error bar on this row.  See [`Sensitivities::fd_error`].
    pub fd_error: f64,
}

/// What `.sens` reports.
#[derive(Debug, Clone)]
pub struct SensResult {
    /// The output's value at the operating point.
    pub out_value: f64,
    /// One row per requested parameter, in the order the card named them (or
    /// netlist order, for the bare `.sens v(out)` form).
    pub rows: Vec<SensRow>,
}

impl SensResult {
    /// The rows whose parameter the adjoint could not reach, for a caller that
    /// would rather fail than read placeholder zeros as insensitivity.
    pub fn unreached(&self) -> Vec<&SensRow> {
        self.rows.iter().filter(|r| !r.reached).collect()
    }

    /// `reached` is a column of its own rather than a footnote: a reader
    /// scanning for the biggest number has to be able to see that a zero was
    /// never computed.
    pub fn write_csv<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        writeln!(w, "param,nominal,sensitivity,normalised,reached,fd_error")?;
        for r in &self.rows {
            writeln!(
                w,
                "{},{:.6e},{:.6e},{:.6e},{},{:.2e}",
                r.name, r.nominal, r.sensitivity, r.normalised, r.reached, r.fd_error
            )?;
        }
        Ok(())
    }

    /// Nutmeg carries the sensitivities alone — it has one value per variable
    /// and no room for the nominal, the error bar or the reached flag.  The CSV
    /// is the complete report; this is the interchange form.
    pub fn write_nutmeg<W: std::io::Write>(&self, mut w: W, title: &str) -> std::io::Result<()> {
        writeln!(w, "Title: {title}")?;
        writeln!(w, "Plotname: Sensitivity Analysis")?;
        writeln!(w, "Flags: real")?;
        writeln!(w, "No. Variables: {}", self.rows.len())?;
        writeln!(w, "No. Points: 1")?;
        writeln!(w, "Variables:")?;
        for (i, r) in self.rows.iter().enumerate() {
            writeln!(w, "\t{i}\t{}\tnotype", r.name)?;
        }
        writeln!(w, "Values:")?;
        for (i, r) in self.rows.iter().enumerate() {
            if i == 0 {
                writeln!(w, " 0\t{:.6e}", r.sensitivity)?;
            } else {
                writeln!(w, "\t{:.6e}", r.sensitivity)?;
            }
        }
        Ok(())
    }
}

/// `.sens <out> [<element>[.<param>] …]` at the DC operating point.
///
/// An empty `params` means every element value in the deck — see
/// [`default_params`] for what that set is and why it stops where it does.
pub fn sensitivity(
    netlist: &Netlist,
    registry: &DeviceRegistry,
    opts: &SimOptions,
    out: &OutVar,
    params: &[ParamName],
) -> Result<SensResult, SimError> {
    let requested: Vec<ParamRef> = if params.is_empty() {
        default_params(netlist)
    } else {
        params
            .iter()
            .map(|p| ParamRef {
                element: p.element.clone(),
                // A bare `r1` means "r1's value"; `set_element_param` accepts
                // `value` for every element it can retune, so there is no
                // per-element table to keep in step here.
                param: p.param.clone().unwrap_or_else(|| "value".into()),
                nominal: None,
                step: None,
            })
            .collect()
    };

    if requested.is_empty() {
        return Err(SimError::ParameterError(
            "`.sens` found no parameters to differentiate: this deck has no R, C, L, V \
             or I element. Name the parameters explicitly (`.sens v(out) m1.w`) if the \
             wanted ones live on a device"
                .into(),
        ));
    }

    let output = match out {
        OutVar::NodeVoltage { pos, neg } => Output::NodeVoltageDiff {
            pos: pos.clone(),
            neg: neg.clone(),
        },
        OutVar::BranchCurrent(name) => Output::BranchCurrent(name.clone()),
    };

    let sens: Sensitivities = dc_sensitivity(netlist, registry, opts, &[output], &requested)?;

    let rows = requested
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let nominal = p
                .nominal
                .or_else(|| crate::netlist_edit::get_element_param(netlist, &p.element, &p.param))
                .unwrap_or(f64::NAN);
            let sensitivity = sens.grad[0][i];
            SensRow {
                name: format!("{}.{}", p.element, p.param),
                nominal,
                sensitivity,
                normalised: sensitivity * nominal,
                reached: sens.reached[i],
                fd_error: sens.fd_error[i],
            }
        })
        .collect();

    Ok(SensResult {
        out_value: sens.values[0],
        rows,
    })
}

/// Every parameter a bare `.sens v(out)` differentiates: the value of each R,
/// C, L, V and I in the deck, in netlist order.
///
/// It stops there on purpose.  Those five are the elements the netlist stamps
/// directly, so their perturbation is exact and always reaches the residual.
/// Device model parameters reach it only where the model implements
/// [`crate::device::Device::set_real_param`], which most do not (see
/// `docs/model_status.md`) — sweeping them all by default would fill the report
/// with unreached rows and bury the ones that mean something.  Name a device
/// parameter explicitly and it is differentiated, and honestly reported if it
/// cannot be.
pub fn default_params(netlist: &Netlist) -> Vec<ParamRef> {
    netlist
        .elements
        .iter()
        .filter_map(|el| {
            let name = match el {
                Element::Resistor { name, .. }
                | Element::Capacitor { name, .. }
                | Element::Inductor { name, .. }
                | Element::VoltageSource { name, .. }
                | Element::CurrentSource { name, .. } => name,
                _ => return None,
            };
            Some(ParamRef {
                element: name.clone(),
                param: "value".into(),
                nominal: None,
                step: None,
            })
        })
        .collect()
}
