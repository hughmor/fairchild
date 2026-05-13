//! PyO3 bindings for the fairchild electro-optic circuit simulator.
//!
//! Exposes `Circuit` and `SimResult` to Python.

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use numpy::PyArray1;

use fairchild_parser::{parse_spice, Element, Netlist};
use fairchild_core::{
    dc_op_nr_with_registry, tran_nr_with_registry_tr, DeviceRegistry, NrResult, TranResult, SimError,
};
use fairchild_osdi::OsdiLibrary;

// ---------------------------------------------------------------------------
// Error conversion
// ---------------------------------------------------------------------------

fn sim_err(e: SimError) -> PyErr {
    PyRuntimeError::new_err(e.to_string())
}

fn parse_err(e: fairchild_parser::ParseError) -> PyErr {
    PyRuntimeError::new_err(e.to_string())
}

// ---------------------------------------------------------------------------
// apply_overrides: patch a cloned Netlist with parameter overrides
// ---------------------------------------------------------------------------

/// Apply a map of `"element_name.param_name" -> value` overrides to `netlist`.
///
/// Matching is case-insensitive for both the element name and the param name.
fn apply_overrides(netlist: &mut Netlist, overrides: &HashMap<String, f64>) {
    for (key, &value) in overrides {
        // Split on the first '.' only.
        let dot = match key.find('.') {
            Some(i) => i,
            None => continue,
        };
        let el_name = key[..dot].to_lowercase();
        let param_name = key[dot + 1..].to_lowercase();

        for el in &mut netlist.elements {
            match el {
                Element::Resistor { name, resistance, .. } => {
                    if name.to_lowercase() == el_name {
                        if param_name == "resistance" || param_name == "value" {
                            *resistance = value;
                        }
                    }
                }
                Element::Capacitor { name, capacitance, .. } => {
                    if name.to_lowercase() == el_name {
                        if param_name == "capacitance" || param_name == "value" {
                            *capacitance = value;
                        }
                    }
                }
                Element::Inductor { name, inductance, .. } => {
                    if name.to_lowercase() == el_name {
                        if param_name == "inductance" || param_name == "value" {
                            *inductance = value;
                        }
                    }
                }
                Element::XOsdi { name, params, .. } => {
                    if name.to_lowercase() == el_name {
                        // Update existing param or push a new one.
                        if let Some(entry) = params.iter_mut().find(|(k, _)| k.to_lowercase() == param_name) {
                            entry.1 = value;
                        } else {
                            // Preserve original capitalisation style by just using the param name
                            // as supplied (after the dot).
                            params.push((key[dot + 1..].to_string(), value));
                        }
                    }
                }
                Element::Mosfet { name, params, .. } => {
                    if name.to_lowercase() == el_name {
                        if let Some(entry) = params.iter_mut().find(|(k, _)| k.to_lowercase() == param_name) {
                            entry.1 = value;
                        } else {
                            params.push((key[dot + 1..].to_string(), value));
                        }
                    }
                }
                // VoltageSource, CurrentSource, Diode — no scalar param override supported yet.
                _ => {}
            }
        }
    }
}

// ---------------------------------------------------------------------------
// SimResult
// ---------------------------------------------------------------------------

/// Result of a single simulation run.
///
/// For DC (`.op`), `time()` returns an empty array and indexed signals return
/// 1-element arrays.  For transient, all arrays have the same length as `time()`.
#[pyclass]
pub struct SimResult {
    inner: SimResultInner,
}

enum SimResultInner {
    Dc(NrResult),
    Tran(TranResult),
}

#[pymethods]
impl SimResult {
    /// 1-D numpy array of time points (empty for DC).
    fn time<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        match &self.inner {
            SimResultInner::Dc(_) => PyArray1::from_vec_bound(py, vec![]),
            SimResultInner::Tran(r) => PyArray1::from_vec_bound(py, r.time.clone()),
        }
    }

    /// Return the waveform for `key`.
    ///
    /// Accepts `"V(node)"` or `"I(vsrc)"` (case-insensitive for the name inside
    /// the parentheses).  For DC results this is a 1-element array.
    fn __getitem__<'py>(&self, py: Python<'py>, key: &str) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let key_lc = key.to_lowercase();

        // Parse "V(name)" or "I(name)".
        if let Some(rest) = key_lc.strip_prefix("v(").and_then(|s| s.strip_suffix(')')) {
            let node = rest;
            match &self.inner {
                SimResultInner::Dc(r) => {
                    let v = r.node_voltage(node).map_err(sim_err)?;
                    Ok(PyArray1::from_vec_bound(py, vec![v]))
                }
                SimResultInner::Tran(r) => {
                    // Ground node is always 0.
                    if node == "0" || node == "gnd" {
                        let zeros = vec![0.0f64; r.time.len()];
                        return Ok(PyArray1::from_vec_bound(py, zeros));
                    }
                    let series = r.node_voltages.get(node)
                        .ok_or_else(|| PyRuntimeError::new_err(format!("unknown node '{node}'")))?;
                    Ok(PyArray1::from_vec_bound(py, series.clone()))
                }
            }
        } else if let Some(rest) = key_lc.strip_prefix("i(").and_then(|s| s.strip_suffix(')')) {
            let vsrc = rest;
            match &self.inner {
                SimResultInner::Dc(r) => {
                    let i = r.vsrc_current(vsrc).map_err(sim_err)?;
                    Ok(PyArray1::from_vec_bound(py, vec![i]))
                }
                SimResultInner::Tran(r) => {
                    let series = r.vsrc_currents.get(vsrc)
                        .ok_or_else(|| PyRuntimeError::new_err(format!("unknown vsrc '{vsrc}'")))?;
                    Ok(PyArray1::from_vec_bound(py, series.clone()))
                }
            }
        } else {
            Err(PyRuntimeError::new_err(format!(
                "unrecognised signal key '{}'; use 'V(node)' or 'I(vsrc)'",
                key
            )))
        }
    }

    /// List of available signal names (e.g. `["V(out)", "I(v1)"]`).
    fn signals(&self) -> Vec<String> {
        match &self.inner {
            SimResultInner::Dc(r) => {
                let mut sigs: Vec<String> = r.topo.node_index.keys()
                    .map(|n| format!("V({n})"))
                    .collect();
                sigs.extend(r.topo.vsrc_index.keys().map(|n| format!("I({n})")));
                sigs
            }
            SimResultInner::Tran(r) => {
                let mut sigs: Vec<String> = r.node_voltages.keys()
                    .map(|n| format!("V({n})"))
                    .collect();
                sigs.extend(r.vsrc_currents.keys().map(|n| format!("I({n})")));
                sigs
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Circuit
// ---------------------------------------------------------------------------

/// A parsed circuit that can run simulations.
///
/// ```python
/// import fairchild as fc
///
/// ckt = fc.Circuit()
/// ckt.load("path/to/netlist.sp")
/// ckt.set_param("Xlaser", "power_mW", 2.0)
/// result = ckt.run("op")
/// ```
#[pyclass]
pub struct Circuit {
    /// The parsed netlist (None until `load()` is called).
    netlist: Option<Netlist>,
    /// Directory containing the netlist file, for resolving relative .osdi paths.
    netlist_dir: Option<PathBuf>,
    /// Accumulated parameter overrides: "element_lower.param_lower" -> value.
    overrides: HashMap<String, f64>,
}

#[pymethods]
impl Circuit {
    #[new]
    pub fn new() -> Self {
        Circuit {
            netlist: None,
            netlist_dir: None,
            overrides: HashMap::new(),
        }
    }

    /// Load a SPICE netlist from `path`.
    pub fn load(&mut self, path: &str) -> PyResult<()> {
        let p = PathBuf::from(path);
        let src = std::fs::read_to_string(&p)
            .map_err(|e| PyRuntimeError::new_err(format!("cannot read '{path}': {e}")))?;
        let netlist = parse_spice(&src).map_err(parse_err)?;
        // Store the directory so we can resolve relative .osdi paths later.
        self.netlist_dir = p.parent().map(|d| d.to_path_buf());
        self.netlist = Some(netlist);
        Ok(())
    }

    /// Load a netlist from a SPICE string (no file path, so relative .osdi paths won't resolve).
    pub fn load_str(&mut self, src: &str) -> PyResult<()> {
        let netlist = parse_spice(src).map_err(parse_err)?;
        self.netlist_dir = None;
        self.netlist = Some(netlist);
        Ok(())
    }

    /// Override a parameter on an element before the next `run()`.
    ///
    /// Parameters:
    ///   element: Element name (e.g. `"Xlaser"`).
    ///   param:   Parameter name (e.g. `"power_mW"`).
    ///   value:   New value as a float.
    pub fn set_param(&mut self, element: &str, param: &str, value: f64) {
        let key = format!("{}.{}", element.to_lowercase(), param.to_lowercase());
        self.overrides.insert(key, value);
    }

    /// Run a simulation analysis.
    ///
    /// Parameters:
    ///   analysis: `"op"` for DC operating point or `"tran"` for transient.
    ///   kwargs:   For transient: `stop` (float, seconds) and `step` (float, seconds).
    ///
    /// Returns a `SimResult`.
    #[pyo3(signature = (analysis, **kwargs))]
    pub fn run(&self, analysis: &str, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<SimResult> {
        let netlist = self.netlist.as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("no netlist loaded; call load() first"))?;

        // Clone the netlist so overrides don't accumulate in self.
        let mut nl = netlist.clone();
        apply_overrides(&mut nl, &self.overrides);

        let registry = build_registry(&nl, self.netlist_dir.as_ref())?;

        match analysis.to_lowercase().as_str() {
            "op" | "dc" => {
                let result = dc_op_nr_with_registry(&nl, &registry).map_err(sim_err)?;
                Ok(SimResult { inner: SimResultInner::Dc(result) })
            }
            "tran" | "transient" => {
                let (stop, step) = parse_tran_kwargs(kwargs)?;
                let result = tran_nr_with_registry_tr(&nl, step, stop, &registry).map_err(sim_err)?;
                Ok(SimResult { inner: SimResultInner::Tran(result) })
            }
            other => Err(PyRuntimeError::new_err(format!(
                "unknown analysis '{}'; use 'op' or 'tran'",
                other
            ))),
        }
    }

    /// Run a parametric sweep.
    ///
    /// Runs `run(analysis, **kwargs)` once per value in `values`, each time
    /// setting `param` (e.g. `"Xlaser.power_mW"`) to that value.
    ///
    /// Returns a list of `SimResult` objects.
    #[pyo3(signature = (param, values, analysis, **kwargs))]
    pub fn sweep(
        &self,
        param: &str,
        values: Vec<f64>,
        analysis: &str,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Vec<SimResult>> {
        let netlist = self.netlist.as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("no netlist loaded; call load() first"))?;

        let mut results = Vec::with_capacity(values.len());

        for val in values {
            let mut nl = netlist.clone();
            // Apply existing overrides first, then the sweep override.
            apply_overrides(&mut nl, &self.overrides);
            let sweep_override: HashMap<String, f64> =
                [(param.to_lowercase(), val)].into_iter().collect();
            apply_overrides(&mut nl, &sweep_override);

            let registry = build_registry(&nl, self.netlist_dir.as_ref())?;

            let result = match analysis.to_lowercase().as_str() {
                "op" | "dc" => {
                    let r = dc_op_nr_with_registry(&nl, &registry).map_err(sim_err)?;
                    SimResult { inner: SimResultInner::Dc(r) }
                }
                "tran" | "transient" => {
                    let (stop, step) = parse_tran_kwargs(kwargs)?;
                    let r = tran_nr_with_registry_tr(&nl, step, stop, &registry).map_err(sim_err)?;
                    SimResult { inner: SimResultInner::Tran(r) }
                }
                other => {
                    return Err(PyRuntimeError::new_err(format!(
                        "unknown analysis '{}'; use 'op' or 'tran'",
                        other
                    )));
                }
            };
            results.push(result);
        }

        Ok(results)
    }
}

// ---------------------------------------------------------------------------
// Helper: build a DeviceRegistry with built-ins + OSDI models
// ---------------------------------------------------------------------------

/// Build a `DeviceRegistry` for `netlist`.
///
/// Registers built-in diodes and MOSFETs, then loads any `.osdi` shared
/// libraries referenced by `.osdi` directives in the netlist.  Relative paths
/// are resolved against `netlist_dir` (if provided).
///
/// Returns the registry and keeps the loaded `OsdiLibrary` objects alive
/// inside `Arc`s — they remain alive because `register_into` clones the Arc
/// into the registry's factory closures.
fn build_registry(netlist: &Netlist, netlist_dir: Option<&PathBuf>) -> PyResult<DeviceRegistry> {
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

        let lib = unsafe { OsdiLibrary::open(&path) }
            .map_err(|e| PyRuntimeError::new_err(format!(
                "failed to load OSDI library '{}': {e}",
                path.display()
            )))?;
        let lib = Arc::new(lib);
        lib.register_into(&mut registry);
    }

    Ok(registry)
}

// ---------------------------------------------------------------------------
// Helper: parse stop/step from keyword arguments
// ---------------------------------------------------------------------------

fn parse_tran_kwargs(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<(f64, f64)> {
    let mut stop: Option<f64> = None;
    let mut step: Option<f64> = None;

    if let Some(kw) = kwargs {
        if let Some(v) = kw.get_item("stop")? {
            stop = Some(v.extract::<f64>()?);
        }
        if let Some(v) = kw.get_item("step")? {
            step = Some(v.extract::<f64>()?);
        }
    }

    let stop = stop.ok_or_else(|| PyRuntimeError::new_err("tran requires 'stop' keyword argument"))?;
    let step = step.ok_or_else(|| PyRuntimeError::new_err("tran requires 'step' keyword argument"))?;

    Ok((stop, step))
}

// ---------------------------------------------------------------------------
// Module entry point
// ---------------------------------------------------------------------------

#[pymodule]
fn fairchild(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Circuit>()?;
    m.add_class::<SimResult>()?;
    Ok(())
}
