//! PyO3 bindings for the fairchild electro-optic circuit simulator.
//!
//! Exposes `Circuit`, `SimResult`, and `WaveformSource` to Python.

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use numpy::{PyArray1, PyReadonlyArray1};

use fairchild_parser::{parse_spice, parse_spice_file, AcVariation, Element, Netlist, Waveform};
use fairchild_core::{
    ac_analysis_opts, dc_op_nr_with_registry_opts, dc_sweep_with_registry_opts,
    evaluate_measurements,
    freq_decade, freq_linear, freq_oct,
    tran_nr_with_registry_opts, tran_nr_with_registry_var_opts,
    AcResult, DcSweepResult, DeviceRegistry, NrResult, SimError, SimOptions, TranResult,
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

fn apply_overrides(netlist: &mut Netlist, overrides: &HashMap<String, f64>) {
    for (key, &value) in overrides {
        let dot = match key.find('.') {
            Some(i) => i,
            None => continue,
        };
        let el_name = key[..dot].to_lowercase();
        let param_name = key[dot + 1..].to_lowercase();

        for el in &mut netlist.elements {
            match el {
                Element::Resistor { name, resistance, .. } => {
                    if name.to_lowercase() == el_name
                        && (param_name == "resistance" || param_name == "value")
                    {
                        *resistance = value;
                    }
                }
                Element::Capacitor { name, capacitance, .. } => {
                    if name.to_lowercase() == el_name
                        && (param_name == "capacitance" || param_name == "value")
                    {
                        *capacitance = value;
                    }
                }
                Element::Inductor { name, inductance, .. } => {
                    if name.to_lowercase() == el_name
                        && (param_name == "inductance" || param_name == "value")
                    {
                        *inductance = value;
                    }
                }
                Element::VoltageSource { name, waveform, .. } => {
                    if name.to_lowercase() == el_name
                        && (param_name == "dc" || param_name == "value" || param_name == "v")
                    {
                        *waveform = Waveform::Dc(value);
                    }
                }
                Element::CurrentSource { name, waveform, .. } => {
                    if name.to_lowercase() == el_name
                        && (param_name == "dc" || param_name == "value" || param_name == "i")
                    {
                        *waveform = Waveform::Dc(value);
                    }
                }
                Element::XOsdi { name, params, .. } => {
                    if name.to_lowercase() == el_name {
                        if let Some(entry) = params.iter_mut().find(|(k, _)| k.to_lowercase() == param_name) {
                            entry.1 = value;
                        } else {
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
                _ => {}
            }
        }
    }
}

// ---------------------------------------------------------------------------
// WaveformSource — numpy t/v arrays → Waveform::Pwl
// ---------------------------------------------------------------------------

/// A piecewise-linear waveform built from numpy time and value arrays.
///
/// Use with `Circuit.set_source()` to inject arbitrary numpy waveforms
/// (e.g. pre-computed eye pattern data or measured signal) into a circuit.
///
/// ```python
/// import numpy as np, fairchild as fc
/// t = np.linspace(0, 10e-9, 1000)
/// v = np.sin(2 * np.pi * 1e9 * t)
/// ckt.set_source("Vin", fc.WaveformSource(t, v))
/// ```
#[pyclass]
pub struct WaveformSource {
    points: Vec<(f64, f64)>,
}

#[pymethods]
impl WaveformSource {
    /// Create a PWL waveform from time and value arrays.
    ///
    /// Both arrays must be 1-D and the same length. Time values must be
    /// non-decreasing.
    #[new]
    pub fn new(
        t: PyReadonlyArray1<f64>,
        v: PyReadonlyArray1<f64>,
    ) -> PyResult<Self> {
        let t = t.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        let v = v.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        if t.len() != v.len() {
            return Err(PyRuntimeError::new_err(
                "WaveformSource: t and v must have the same length",
            ));
        }
        let points: Vec<(f64, f64)> = t.iter().copied().zip(v.iter().copied()).collect();
        Ok(Self { points })
    }

    fn __repr__(&self) -> String {
        format!("WaveformSource({} points, t=[{:.3e}..{:.3e}])",
            self.points.len(),
            self.points.first().map(|(t, _)| *t).unwrap_or(0.0),
            self.points.last().map(|(t, _)| *t).unwrap_or(0.0),
        )
    }
}

// ---------------------------------------------------------------------------
// SimResult
// ---------------------------------------------------------------------------

enum SimResultInner {
    Dc(NrResult),
    Tran(TranResult),
    Ac(AcResult),
    DcSweep(DcSweepResult),
}

/// Result of a simulation run.
///
/// For DC (`.op`), `time()` returns an empty array and indexed signals return
/// 1-element arrays.  For transient, all arrays have the same length as `time()`.
/// For AC, `freq()` returns the frequency array; signals return complex magnitude
/// arrays.
#[pyclass]
pub struct SimResult {
    inner: SimResultInner,
    /// `.measure` scalars produced from this run.  Empty for analyses that
    /// don't support measurements (DC OP, AC, DC sweep).
    measurements: Vec<(String, f64)>,
}

#[pymethods]
impl SimResult {
    /// 1-D numpy array of time points (empty for DC and AC).
    fn time<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        match &self.inner {
            SimResultInner::Tran(r) => PyArray1::from_vec_bound(py, r.time.clone()),
            _ => PyArray1::from_vec_bound(py, vec![]),
        }
    }

    /// 1-D numpy array of frequencies in Hz (empty for DC and transient).
    fn freq<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        match &self.inner {
            SimResultInner::Ac(r) => PyArray1::from_vec_bound(py, r.freq.clone()),
            _ => PyArray1::from_vec_bound(py, vec![]),
        }
    }

    /// 1-D numpy array of sweep values (only meaningful for DC sweeps).
    fn sweep<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        match &self.inner {
            SimResultInner::DcSweep(r) => PyArray1::from_vec_bound(py, r.outer.values.clone()),
            _ => PyArray1::from_vec_bound(py, vec![]),
        }
    }

    /// Return the waveform for `key`.
    ///
    /// DC/transient: accepts `"V(node)"` or `"I(vsrc)"`.
    ///
    /// AC: accepts `"V(node)"` (magnitude), `"V(node).mag"` (same),
    /// `"V(node).phase"` (degrees), `"V(node).re"`, `"V(node).im"`.
    fn __getitem__<'py>(&self, py: Python<'py>, key: &str) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let key_lc = key.to_lowercase();

        match &self.inner {
            SimResultInner::Dc(r) => {
                self.get_dc_or_tran_signal(py, &key_lc, r, None)
            }
            SimResultInner::Tran(r_tran) => {
                self.get_tran_signal(py, &key_lc, r_tran)
            }
            SimResultInner::Ac(r) => {
                self.get_ac_signal(py, &key_lc, r)
            }
            SimResultInner::DcSweep(r) => {
                self.get_dc_sweep_signal(py, &key_lc, r)
            }
        }
    }

    /// List of available signal names.
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
            SimResultInner::Ac(r) => {
                r.voltages.keys().map(|n| format!("V({n})")).collect()
            }
            SimResultInner::DcSweep(r) => {
                let mut sigs: Vec<String> = r.node_voltages.keys()
                    .map(|n| format!("V({n})"))
                    .collect();
                sigs.extend(r.vsrc_currents.keys().map(|n| format!("I({n})")));
                sigs
            }
        }
    }

    /// True if this is a DC operating-point result.
    #[getter]
    fn is_dc(&self) -> bool { matches!(&self.inner, SimResultInner::Dc(_)) }

    /// True if this is a transient result.
    #[getter]
    fn is_tran(&self) -> bool { matches!(&self.inner, SimResultInner::Tran(_)) }

    /// True if this is an AC sweep result.
    #[getter]
    fn is_ac(&self) -> bool { matches!(&self.inner, SimResultInner::Ac(_)) }

    /// True if this is a DC sweep result.
    #[getter]
    fn is_dc_sweep(&self) -> bool { matches!(&self.inner, SimResultInner::DcSweep(_)) }

    /// Return all `.measure` scalar values produced from this run.
    ///
    /// Returns a Python dict mapping measurement name → value.  Empty for
    /// analyses that don't support measurements (DC OP, AC, DC sweep).
    fn measurements<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new_bound(py);
        for (k, v) in &self.measurements {
            d.set_item(k, *v)?;
        }
        Ok(d)
    }
}

impl SimResult {
    fn get_dc_or_tran_signal<'py>(
        &self,
        py: Python<'py>,
        key: &str,
        r: &NrResult,
        _dummy: Option<()>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        if let Some(node) = key.strip_prefix("v(").and_then(|s| s.strip_suffix(')')) {
            let v = r.node_voltage(node).map_err(sim_err)?;
            Ok(PyArray1::from_vec_bound(py, vec![v]))
        } else if let Some(vsrc) = key.strip_prefix("i(").and_then(|s| s.strip_suffix(')')) {
            let i = r.vsrc_current(vsrc).map_err(sim_err)?;
            Ok(PyArray1::from_vec_bound(py, vec![i]))
        } else {
            Err(PyRuntimeError::new_err(format!(
                "unrecognised signal key '{key}'; use 'V(node)' or 'I(vsrc)'"
            )))
        }
    }

    fn get_tran_signal<'py>(
        &self,
        py: Python<'py>,
        key: &str,
        r: &TranResult,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        if let Some(node) = key.strip_prefix("v(").and_then(|s| s.strip_suffix(')')) {
            if node == "0" || node == "gnd" {
                return Ok(PyArray1::from_vec_bound(py, vec![0.0f64; r.time.len()]));
            }
            let series = r.node_voltages.get(node)
                .ok_or_else(|| PyRuntimeError::new_err(format!("unknown node '{node}'")))?;
            Ok(PyArray1::from_vec_bound(py, series.clone()))
        } else if let Some(vsrc) = key.strip_prefix("i(").and_then(|s| s.strip_suffix(')')) {
            let series = r.vsrc_currents.get(vsrc)
                .ok_or_else(|| PyRuntimeError::new_err(format!("unknown vsrc '{vsrc}'")))?;
            Ok(PyArray1::from_vec_bound(py, series.clone()))
        } else {
            Err(PyRuntimeError::new_err(format!(
                "unrecognised signal key '{key}'; use 'V(node)' or 'I(vsrc)'"
            )))
        }
    }

    fn get_dc_sweep_signal<'py>(
        &self,
        py: Python<'py>,
        key: &str,
        r: &DcSweepResult,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        if let Some(node) = key.strip_prefix("v(").and_then(|s| s.strip_suffix(')')) {
            if node == "0" || node == "gnd" {
                let n = r.n_points();
                return Ok(PyArray1::from_vec_bound(py, vec![0.0f64; n]));
            }
            let series = r.node_voltages.get(node)
                .ok_or_else(|| PyRuntimeError::new_err(format!("unknown node '{node}'")))?;
            Ok(PyArray1::from_vec_bound(py, series.clone()))
        } else if let Some(vsrc) = key.strip_prefix("i(").and_then(|s| s.strip_suffix(')')) {
            let series = r.vsrc_currents.get(vsrc)
                .ok_or_else(|| PyRuntimeError::new_err(format!("unknown vsrc '{vsrc}'")))?;
            Ok(PyArray1::from_vec_bound(py, series.clone()))
        } else {
            Err(PyRuntimeError::new_err(format!(
                "unrecognised signal key '{key}'; use 'V(node)' or 'I(vsrc)'"
            )))
        }
    }

    fn get_ac_signal<'py>(
        &self,
        py: Python<'py>,
        key: &str,
        r: &AcResult,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        // Parse "V(node)", "V(node).mag", "V(node).phase", "V(node).re", "V(node).im"
        let (node, suffix) = if let Some(rest) = key.strip_prefix("v(") {
            if let Some(dot_idx) = rest.find(").") {
                (&rest[..dot_idx], Some(&rest[dot_idx + 2..]))
            } else if let Some(bare) = rest.strip_suffix(')') {
                (bare, None)
            } else {
                return Err(PyRuntimeError::new_err(format!("malformed signal key '{key}'")));
            }
        } else {
            return Err(PyRuntimeError::new_err(format!(
                "AC result: unrecognised key '{key}'; use 'V(node)', 'V(node).mag', .phase, .re, .im"
            )));
        };

        let voltages = r.voltages.get(node)
            .ok_or_else(|| PyRuntimeError::new_err(format!("unknown AC node '{node}'")))?;

        let data: Vec<f64> = match suffix.unwrap_or("mag") {
            "mag" | "magnitude" => {
                voltages.iter().map(|(re, im)| (re * re + im * im).sqrt()).collect()
            }
            "phase" => {
                voltages.iter().map(|(re, im)| im.atan2(*re).to_degrees()).collect()
            }
            "re" | "real" => {
                voltages.iter().map(|(re, _)| *re).collect()
            }
            "im" | "imag" | "imaginary" => {
                voltages.iter().map(|(_, im)| *im).collect()
            }
            "db" => {
                voltages.iter().map(|(re, im)| {
                    20.0 * (re * re + im * im).sqrt().max(1e-300).log10()
                }).collect()
            }
            other => {
                return Err(PyRuntimeError::new_err(format!(
                    "AC suffix '{other}' not recognised; use mag, phase, re, im, db"
                )));
            }
        };

        Ok(PyArray1::from_vec_bound(py, data))
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
/// print(result["V(ph_a)"])
/// ```
#[pyclass]
pub struct Circuit {
    netlist: Option<Netlist>,
    netlist_dir: Option<PathBuf>,
    overrides: HashMap<String, f64>,
    source_overrides: HashMap<String, Vec<(f64, f64)>>,
}

#[pymethods]
impl Circuit {
    #[new]
    pub fn new() -> Self {
        Circuit {
            netlist: None,
            netlist_dir: None,
            overrides: HashMap::new(),
            source_overrides: HashMap::new(),
        }
    }

    /// Load a SPICE netlist from `path`, resolving `.include` directives
    /// relative to the file's parent directory.
    pub fn load(&mut self, path: &str) -> PyResult<()> {
        let p = PathBuf::from(path);
        let netlist = parse_spice_file(&p).map_err(parse_err)?;
        self.netlist_dir = p.parent().map(|d| d.to_path_buf());
        self.netlist = Some(netlist);
        Ok(())
    }

    /// Load a netlist from a SPICE string.
    pub fn load_str(&mut self, src: &str) -> PyResult<()> {
        let netlist = parse_spice(src).map_err(parse_err)?;
        self.netlist_dir = None;
        self.netlist = Some(netlist);
        Ok(())
    }

    /// Override a parameter on an element before the next `run()`.
    pub fn set_param(&mut self, element: &str, param: &str, value: f64) {
        let key = format!("{}.{}", element.to_lowercase(), param.to_lowercase());
        self.overrides.insert(key, value);
    }

    /// Inject a numpy waveform as the source for a voltage or current source.
    ///
    /// The source named `name` (e.g. `"Vin"`, `"V1"`) will have its waveform
    /// replaced with a PWL interpolation of the provided `WaveformSource`.
    /// Call this before `run()`.
    pub fn set_source(&mut self, name: &str, source: &WaveformSource) {
        self.source_overrides.insert(name.to_lowercase(), source.points.clone());
    }

    /// Run a simulation analysis.
    ///
    /// Parameters:
    ///   analysis: `"op"` for DC, `"tran"` for transient, `"ac"` for AC sweep.
    ///
    ///   For `"tran"`: `stop` (s) and `step` (s) are required.
    ///   Optional: `variable_step=True` enables LTE-controlled variable-step.
    ///
    ///   For `"ac"`: `fstart` (Hz), `fstop` (Hz), `points` (int, default 20),
    ///   `variation` (`"dec"`, `"oct"`, `"lin"`, default `"dec"`),
    ///   `src` (excitation source name, default `None` = first V source).
    ///
    ///   Solver options (apply to all analyses): `reltol`, `abstol`, `vntol`,
    ///   `vmax`, `gmin`, `itl1`, `itl4`, `maxstep`,
    ///   `method` (`"be"` | `"tr"` | `"gear"`), `uic`, `temp` (°C).
    ///   These overlay any `.options` directives from the netlist.
    #[pyo3(signature = (analysis, **kwargs))]
    pub fn run(&self, analysis: &str, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<SimResult> {
        let netlist = self.netlist.as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("no netlist loaded; call load() first"))?;

        let mut nl = netlist.clone();
        apply_overrides(&mut nl, &self.overrides);
        apply_source_overrides(&mut nl, &self.source_overrides);

        let registry = build_registry(&nl, self.netlist_dir.as_ref())?;
        let opts = build_sim_options(&nl, kwargs)?;

        let analysis_lc = analysis.to_lowercase();
        // "dc" with a src kwarg is a sweep; without one it's an op-point alias.
        let is_dc_sweep = analysis_lc == "dc_sweep"
            || (analysis_lc == "dc"
                && kwargs.and_then(|kw| kw.get_item("src").ok().flatten()).is_some());

        if is_dc_sweep {
            let p = parse_dc_kwargs(kwargs)?;
            let nested_arg = p.nested.as_ref()
                .map(|(s, a, b, st)| (s.as_str(), *a, *b, *st));
            let result = dc_sweep_with_registry_opts(
                &nl, &p.src, p.start, p.stop, p.step, nested_arg, &registry, &opts
            ).map_err(sim_err)?;
            return Ok(SimResult {
                inner: SimResultInner::DcSweep(result),
                measurements: Vec::new(),
            });
        }

        match analysis_lc.as_str() {
            "op" | "dc" => {
                let result = dc_op_nr_with_registry_opts(&nl, &registry, &opts).map_err(sim_err)?;
                Ok(SimResult {
                    inner: SimResultInner::Dc(result),
                    measurements: Vec::new(),
                })
            }
            "tran" | "transient" => {
                let (stop, step, variable_step) = parse_tran_kwargs(kwargs)?;
                let result = if variable_step {
                    tran_nr_with_registry_var_opts(&nl, step, stop, &registry, &opts)
                } else {
                    tran_nr_with_registry_opts(&nl, step, stop, &registry, &opts)
                }.map_err(sim_err)?;
                let measurements = evaluate_measurements(&nl.measurements, &result)
                    .into_iter().map(|m| (m.name, m.value)).collect();
                Ok(SimResult { inner: SimResultInner::Tran(result), measurements })
            }
            "ac" => {
                let (freqs, src) = parse_ac_kwargs(kwargs)?;
                let result = ac_analysis_opts(&nl, &freqs, src.as_deref(), &registry, &opts).map_err(sim_err)?;
                Ok(SimResult { inner: SimResultInner::Ac(result), measurements: Vec::new() })
            }
            other => Err(PyRuntimeError::new_err(format!(
                "unknown analysis '{}'; use 'op', 'tran', 'ac', or 'dc_sweep'",
                other
            ))),
        }
    }

    /// Run a parametric sweep over scalar element parameters.
    ///
    /// Calls `run(analysis, **kwargs)` once per value in `values`, each time
    /// setting `param` (e.g. `"Xlaser.power_mW"`) to that value.
    ///
    /// Returns a list of `SimResult` objects, one per sweep point.
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
            apply_overrides(&mut nl, &self.overrides);
            apply_source_overrides(&mut nl, &self.source_overrides);
            let sweep_override: HashMap<String, f64> =
                [(param.to_lowercase(), val)].into_iter().collect();
            apply_overrides(&mut nl, &sweep_override);

            let registry = build_registry(&nl, self.netlist_dir.as_ref())?;
            let opts = build_sim_options(&nl, kwargs)?;

            let result = match analysis.to_lowercase().as_str() {
                "op" | "dc" => {
                    let r = dc_op_nr_with_registry_opts(&nl, &registry, &opts).map_err(sim_err)?;
                    SimResult { inner: SimResultInner::Dc(r), measurements: Vec::new() }
                }
                "tran" | "transient" => {
                    let (stop, step, variable_step) = parse_tran_kwargs(kwargs)?;
                    let r = if variable_step {
                        tran_nr_with_registry_var_opts(&nl, step, stop, &registry, &opts)
                    } else {
                        tran_nr_with_registry_opts(&nl, step, stop, &registry, &opts)
                    }.map_err(sim_err)?;
                    let measurements = evaluate_measurements(&nl.measurements, &r)
                        .into_iter().map(|m| (m.name, m.value)).collect();
                    SimResult { inner: SimResultInner::Tran(r), measurements }
                }
                "ac" => {
                    let (freqs, src) = parse_ac_kwargs(kwargs)?;
                    let r = ac_analysis_opts(&nl, &freqs, src.as_deref(), &registry, &opts).map_err(sim_err)?;
                    SimResult { inner: SimResultInner::Ac(r), measurements: Vec::new() }
                }
                other => {
                    return Err(PyRuntimeError::new_err(format!(
                        "unknown analysis '{}'; use 'op', 'tran', or 'ac'", other
                    )));
                }
            };
            results.push(result);
        }

        Ok(results)
    }
}

// ---------------------------------------------------------------------------
// Helper: apply source (WaveformSource) overrides
// ---------------------------------------------------------------------------

fn apply_source_overrides(netlist: &mut Netlist, overrides: &HashMap<String, Vec<(f64, f64)>>) {
    for (name_lc, points) in overrides {
        for el in &mut netlist.elements {
            match el {
                Element::VoltageSource { name, waveform, .. } if name.to_lowercase() == *name_lc => {
                    *waveform = Waveform::Pwl { points: points.clone() };
                }
                Element::CurrentSource { name, waveform, .. } if name.to_lowercase() == *name_lc => {
                    *waveform = Waveform::Pwl { points: points.clone() };
                }
                _ => {}
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Helper: build a DeviceRegistry
// ---------------------------------------------------------------------------

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
                "failed to load OSDI library '{}': {e}", path.display()
            )))?;
        let lib = Arc::new(lib);
        lib.register_into(&mut registry);
    }

    Ok(registry)
}

// ---------------------------------------------------------------------------
// Helper: parse tran kwargs
// ---------------------------------------------------------------------------

fn parse_tran_kwargs(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<(f64, f64, bool)> {
    let mut stop: Option<f64> = None;
    let mut step: Option<f64> = None;
    let mut variable_step = false;

    if let Some(kw) = kwargs {
        if let Some(v) = kw.get_item("stop")? { stop = Some(v.extract::<f64>()?); }
        if let Some(v) = kw.get_item("step")? { step = Some(v.extract::<f64>()?); }
        if let Some(v) = kw.get_item("variable_step")? {
            variable_step = v.extract::<bool>()?;
        }
    }

    let stop = stop.ok_or_else(|| PyRuntimeError::new_err("tran requires 'stop' kwarg (seconds)"))?;
    let step = step.ok_or_else(|| PyRuntimeError::new_err("tran requires 'step' kwarg (seconds)"))?;
    Ok((stop, step, variable_step))
}

// ---------------------------------------------------------------------------
// Helper: build SimOptions from netlist + Python kwargs
// ---------------------------------------------------------------------------

/// Build a `SimOptions` by starting from the netlist's `.options` directives
/// and overlaying any Python kwargs the user passed.  Kwargs that are not
/// recognised as solver options are silently ignored (they may be analysis
/// kwargs like `stop`/`step`/`fstart`).
fn build_sim_options(netlist: &Netlist, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<SimOptions> {
    let mut opts = SimOptions::from_netlist(netlist);
    if let Some(kw) = kwargs {
        // Map of Python kwarg name -> SimOptions key (most are identical).
        const OPTION_KEYS: &[&str] = &[
            "reltol", "abstol", "vntol", "vmax", "gmin",
            "itl1", "itl4", "maxstep", "max_step", "gmin_max", "srcsteps",
            "method", "uic", "temp", "pnjlim",
        ];
        for key in OPTION_KEYS {
            if let Some(v) = kw.get_item(key)? {
                // Accept either a raw number or a string token.
                let value_str: String = if let Ok(s) = v.extract::<String>() {
                    s
                } else if let Ok(f) = v.extract::<f64>() {
                    format!("{f:e}")
                } else if let Ok(i) = v.extract::<i64>() {
                    i.to_string()
                } else if let Ok(b) = v.extract::<bool>() {
                    if b { "1".into() } else { "0".into() }
                } else {
                    return Err(PyRuntimeError::new_err(format!(
                        "kwarg '{key}': expected number, string, or bool"
                    )));
                };
                if !opts.set(key, &value_str) {
                    return Err(PyRuntimeError::new_err(format!(
                        "unknown solver option '{key}={value_str}'"
                    )));
                }
            }
        }
    }
    Ok(opts)
}

// ---------------------------------------------------------------------------
// Helper: parse dc-sweep kwargs
// ---------------------------------------------------------------------------

#[derive(Default)]
struct DcKwargs {
    src:    String,
    start:  f64,
    stop:   f64,
    step:   f64,
    nested: Option<(String, f64, f64, f64)>,
}

fn parse_dc_kwargs(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<DcKwargs> {
    let kw = kwargs.ok_or_else(|| PyRuntimeError::new_err(
        "dc requires src, start, stop, step kwargs"))?;

    let get = |k: &str| -> PyResult<_> {
        kw.get_item(k)?.ok_or_else(|| PyRuntimeError::new_err(
            format!("dc requires '{k}' kwarg")))
    };

    let src   = get("src")?.extract::<String>()?;
    let start = get("start")?.extract::<f64>()?;
    let stop  = get("stop")?.extract::<f64>()?;
    let step  = get("step")?.extract::<f64>()?;

    // Optional nested second sweep: src2, start2, stop2, step2 (all four required if any).
    let nested = match kw.get_item("src2")? {
        Some(v) => {
            let src2   = v.extract::<String>()?;
            let start2 = get("start2")?.extract::<f64>()?;
            let stop2  = get("stop2")?.extract::<f64>()?;
            let step2  = get("step2")?.extract::<f64>()?;
            Some((src2, start2, stop2, step2))
        }
        None => None,
    };

    Ok(DcKwargs { src, start, stop, step, nested })
}

// ---------------------------------------------------------------------------
// Helper: parse ac kwargs
// ---------------------------------------------------------------------------

fn parse_ac_kwargs(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<(Vec<f64>, Option<String>)> {
    let mut fstart: Option<f64> = None;
    let mut fstop: Option<f64> = None;
    let mut points: usize = 20;
    let mut variation = AcVariation::Dec;
    let mut src: Option<String> = None;

    if let Some(kw) = kwargs {
        if let Some(v) = kw.get_item("fstart")? { fstart = Some(v.extract::<f64>()?); }
        if let Some(v) = kw.get_item("fstop")?  { fstop  = Some(v.extract::<f64>()?); }
        if let Some(v) = kw.get_item("points")? { points = v.extract::<usize>()?; }
        if let Some(v) = kw.get_item("src")?    {
            let s: Option<String> = v.extract()?;
            src = s;
        }
        if let Some(v) = kw.get_item("variation")? {
            let var: String = v.extract()?;
            variation = match var.to_lowercase().as_str() {
                "dec" => AcVariation::Dec,
                "oct" => AcVariation::Oct,
                "lin" => AcVariation::Lin,
                other => return Err(PyRuntimeError::new_err(format!(
                    "unknown AC variation '{other}'; use 'dec', 'oct', or 'lin'"
                ))),
            };
        }
    }

    let fstart = fstart.ok_or_else(|| PyRuntimeError::new_err("ac requires 'fstart' kwarg (Hz)"))?;
    let fstop  = fstop .ok_or_else(|| PyRuntimeError::new_err("ac requires 'fstop' kwarg (Hz)"))?;

    let freqs = match variation {
        AcVariation::Dec => freq_decade(fstart, fstop, points),
        AcVariation::Oct => freq_oct(fstart, fstop, points),
        AcVariation::Lin => freq_linear(fstart, fstop, points),
    };

    Ok((freqs, src))
}

// ---------------------------------------------------------------------------
// Module entry point
// ---------------------------------------------------------------------------

#[pymodule]
fn fairchild(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Circuit>()?;
    m.add_class::<SimResult>()?;
    m.add_class::<WaveformSource>()?;
    Ok(())
}
