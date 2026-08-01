//! PyO3 bindings for the fairchild electro-optic circuit simulator.
//!
//! Exposes `Circuit`, `SimResult`, and `WaveformSource` to Python.

// pyo3's #[pymethods] macro expansion triggers this lint on PyResult<T> return
// types as a false positive; suppress it for the whole crate.
#![allow(clippy::useless_conversion)]

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use fairchild_core::{
    ac_analysis_opts, dc_op_nr_with_registry_opts, dc_sweep_with_registry_opts,
    evaluate_measurements, freq_decade, freq_linear, freq_oct, tran_nr_with_registry_opts,
    tran_nr_with_registry_var_opts, AcResult, DcSweepResult, DeviceRegistry, NrResult, Output,
    ParamRef, SimError, SimOptions, TranAdjoint, TranResult,
};
#[cfg(feature = "osdi")]
use fairchild_osdi::OsdiLibrary;
use fairchild_parser::{parse_spice, parse_spice_file, AcVariation, Netlist};

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
        // Keys arrive as "element.param" from `Circuit.set_param`.
        let Some(dot) = key.find('.') else { continue };
        fairchild_core::set_element_param(netlist, &key[..dot], &key[dot + 1..], value);
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
    pub fn new(t: PyReadonlyArray1<f64>, v: PyReadonlyArray1<f64>) -> PyResult<Self> {
        let t = t
            .as_slice()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        let v = v
            .as_slice()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        if t.len() != v.len() {
            return Err(PyRuntimeError::new_err(
                "WaveformSource: t and v must have the same length",
            ));
        }
        let points: Vec<(f64, f64)> = t.iter().copied().zip(v.iter().copied()).collect();
        Ok(Self { points })
    }

    fn __repr__(&self) -> String {
        format!(
            "WaveformSource({} points, t=[{:.3e}..{:.3e}])",
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
    Noise(fairchild_core::NoiseResult),
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
#[allow(clippy::useless_conversion)]
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
            SimResultInner::Noise(r) => PyArray1::from_vec_bound(py, r.freq.clone()),
            _ => PyArray1::from_vec_bound(py, vec![]),
        }
    }

    /// Output-referred voltage noise PSD in V²/Hz (only meaningful for noise).
    fn onoise<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        match &self.inner {
            SimResultInner::Noise(r) => PyArray1::from_vec_bound(py, r.onoise_psd.clone()),
            _ => PyArray1::from_vec_bound(py, vec![]),
        }
    }

    /// Input-referred voltage noise PSD in V²/Hz (only meaningful for noise).
    fn inoise<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        match &self.inner {
            SimResultInner::Noise(r) => PyArray1::from_vec_bound(py, r.inoise_psd.clone()),
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
            SimResultInner::Dc(r) => self.get_dc_or_tran_signal(py, &key_lc, r, None),
            SimResultInner::Tran(r_tran) => self.get_tran_signal(py, &key_lc, r_tran),
            SimResultInner::Ac(r) => self.get_ac_signal(py, &key_lc, r),
            SimResultInner::DcSweep(r) => self.get_dc_sweep_signal(py, &key_lc, r),
            SimResultInner::Noise(r) => {
                // Recognise the same keys as the CSV output.
                let arr = match key_lc.as_str() {
                    "onoise" | "v(onoise)" | "onoise_psd" => r.onoise_psd.clone(),
                    "inoise" | "v(inoise)" | "inoise_psd" => r.inoise_psd.clone(),
                    "onoise_vrthz" => r.onoise_psd.iter().map(|x| x.max(0.0).sqrt()).collect(),
                    "inoise_vrthz" => r.inoise_psd.iter().map(|x| x.max(0.0).sqrt()).collect(),
                    other => {
                        return Err(PyRuntimeError::new_err(format!(
                            "noise result: unknown key '{other}'; use 'onoise', 'inoise', \
                         'onoise_vrthz', or 'inoise_vrthz'"
                        )))
                    }
                };
                Ok(PyArray1::from_vec_bound(py, arr))
            }
        }
    }

    /// List of available signal names.
    fn signals(&self) -> Vec<String> {
        match &self.inner {
            SimResultInner::Dc(r) => {
                let mut sigs: Vec<String> = r
                    .topo
                    .node_index
                    .keys()
                    .map(|n| format!("V({n})"))
                    .collect();
                sigs.extend(r.topo.vsrc_index.keys().map(|n| format!("I({n})")));
                sigs
            }
            SimResultInner::Tran(r) => {
                let mut sigs: Vec<String> =
                    r.node_voltages.keys().map(|n| format!("V({n})")).collect();
                sigs.extend(r.vsrc_currents.keys().map(|n| format!("I({n})")));
                sigs
            }
            SimResultInner::Ac(r) => r.voltages.keys().map(|n| format!("V({n})")).collect(),
            SimResultInner::DcSweep(r) => {
                let mut sigs: Vec<String> =
                    r.node_voltages.keys().map(|n| format!("V({n})")).collect();
                sigs.extend(r.vsrc_currents.keys().map(|n| format!("I({n})")));
                sigs
            }
            SimResultInner::Noise(_) => {
                vec![
                    "onoise".into(),
                    "inoise".into(),
                    "onoise_vrthz".into(),
                    "inoise_vrthz".into(),
                ]
            }
        }
    }

    /// True if this is a DC operating-point result.
    #[getter]
    fn is_dc(&self) -> bool {
        matches!(&self.inner, SimResultInner::Dc(_))
    }

    /// True if this is a transient result.
    #[getter]
    fn is_tran(&self) -> bool {
        matches!(&self.inner, SimResultInner::Tran(_))
    }

    /// True if this is an AC sweep result.
    #[getter]
    fn is_ac(&self) -> bool {
        matches!(&self.inner, SimResultInner::Ac(_))
    }

    /// True if this is a DC sweep result.
    #[getter]
    fn is_dc_sweep(&self) -> bool {
        matches!(&self.inner, SimResultInner::DcSweep(_))
    }

    /// Return all `.measure` scalar values produced from this run.
    ///
    /// Returns a Python dict mapping measurement name → value.  Empty for
    /// analyses that don't support measurements (DC OP, AC, DC sweep).
    #[allow(clippy::useless_conversion)]
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
            let series = r
                .node_voltages
                .get(node)
                .ok_or_else(|| PyRuntimeError::new_err(format!("unknown node '{node}'")))?;
            Ok(PyArray1::from_vec_bound(py, series.clone()))
        } else if let Some(vsrc) = key.strip_prefix("i(").and_then(|s| s.strip_suffix(')')) {
            let series = r
                .vsrc_currents
                .get(vsrc)
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
            let series = r
                .node_voltages
                .get(node)
                .ok_or_else(|| PyRuntimeError::new_err(format!("unknown node '{node}'")))?;
            Ok(PyArray1::from_vec_bound(py, series.clone()))
        } else if let Some(vsrc) = key.strip_prefix("i(").and_then(|s| s.strip_suffix(')')) {
            let series = r
                .vsrc_currents
                .get(vsrc)
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
                return Err(PyRuntimeError::new_err(format!(
                    "malformed signal key '{key}'"
                )));
            }
        } else {
            return Err(PyRuntimeError::new_err(format!(
                "AC result: unrecognised key '{key}'; use 'V(node)', 'V(node).mag', .phase, .re, .im"
            )));
        };

        let voltages = r
            .voltages
            .get(node)
            .ok_or_else(|| PyRuntimeError::new_err(format!("unknown AC node '{node}'")))?;

        let data: Vec<f64> = match suffix.unwrap_or("mag") {
            "mag" | "magnitude" => voltages
                .iter()
                .map(|(re, im)| (re * re + im * im).sqrt())
                .collect(),
            "phase" => voltages
                .iter()
                .map(|(re, im)| im.atan2(*re).to_degrees())
                .collect(),
            "re" | "real" => voltages.iter().map(|(re, _)| *re).collect(),
            "im" | "imag" | "imaginary" => voltages.iter().map(|(_, im)| *im).collect(),
            "db" => voltages
                .iter()
                .map(|(re, im)| 20.0 * (re * re + im * im).sqrt().max(1e-300).log10())
                .collect(),
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

impl Default for Circuit {
    fn default() -> Self {
        Self::new()
    }
}

#[pymethods]
#[allow(clippy::useless_conversion)]
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
    #[allow(clippy::useless_conversion)]
    pub fn load(&mut self, path: &str) -> PyResult<()> {
        let p = PathBuf::from(path);
        let netlist = parse_spice_file(&p).map_err(parse_err)?;
        self.netlist_dir = p.parent().map(|d| d.to_path_buf());
        self.netlist = Some(netlist);
        Ok(())
    }

    /// Load a netlist from a SPICE string.
    #[allow(clippy::useless_conversion)]
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

    /// List the `.alter` block labels declared in the loaded netlist.
    ///
    /// Pass one of these as the `alter` kwarg to `run()` to apply that
    /// block on top of the base netlist for that call.
    pub fn alter_labels(&self) -> Vec<String> {
        self.netlist
            .as_ref()
            .map(|n| n.alters.iter().map(|b| b.label.clone()).collect())
            .unwrap_or_default()
    }

    /// Inject a numpy waveform as the source for a voltage or current source.
    ///
    /// The source named `name` (e.g. `"Vin"`, `"V1"`) will have its waveform
    /// replaced with a PWL interpolation of the provided `WaveformSource`.
    /// Call this before `run()`.
    pub fn set_source(&mut self, name: &str, source: &WaveformSource) {
        self.source_overrides
            .insert(name.to_lowercase(), source.points.clone());
    }

    /// Run a simulation analysis.
    ///
    /// Parameters:
    ///   analysis: `"op"` for DC, `"tran"` for transient, `"ac"` for AC sweep.
    ///
    ///   For `"tran"`: `stop` (s) and `step` (s) are required.
    ///
    ///   For `"ac"`: `fstart` (Hz), `fstop` (Hz), `points` (int, default 20),
    ///   `variation` (`"dec"`, `"oct"`, `"lin"`, default `"dec"`),
    ///   `src` (excitation source name, default `None` = first V source).
    ///
    ///   Solver options (apply to all analyses): `reltol`, `abstol`, `vntol`,
    ///   `lambdatol`, `vmax`, `gmin`, `itl1`, `itl4`, `maxstep`,
    ///   `method` (`"be"` | `"tr"` | `"gear"`), `uic`, `temp` (°C),
    ///   `solver` (`"dense"` | `"sparse"` | `"auto"` | `"klu"`),
    ///   `variable_step` (bool, enables LTE-controlled variable-step transient),
    ///   `waveguide_delay` (bool), `cond_estimate` (bool, prints κ(A)),
    ///   `equilibrate` (bool, two-sided matrix scaling before LU),
    ///   `lambda_center_nm`, `enable_bidirectional`, `sanity_check`, and more.
    ///   These overlay any `.options` directives from the netlist.
    ///   Unrecognised kwargs raise `RuntimeError` immediately.
    #[pyo3(signature = (analysis, **kwargs))]
    #[allow(clippy::useless_conversion)]
    pub fn run(
        &self,
        py: Python<'_>,
        analysis: &str,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<SimResult> {
        let netlist = self
            .netlist
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("no netlist loaded; call load() first"))?;

        let mut nl = netlist.clone();
        apply_overrides(&mut nl, &self.overrides);
        apply_source_overrides(&mut nl, &self.source_overrides);

        // Apply a `.alter` block if the user requested one via `alter=...`.
        if let Some(kw) = kwargs {
            if let Some(v) = kw.get_item("alter")? {
                let label: String = v.extract()?;
                let block = nl
                    .alters
                    .iter()
                    .find(|b| b.label == label)
                    .cloned()
                    .ok_or_else(|| {
                        PyRuntimeError::new_err(format!(
                            "no .alter block named '{label}'; available: {:?}",
                            nl.alters.iter().map(|b| &b.label).collect::<Vec<_>>()
                        ))
                    })?;
                nl.apply_alter(&block);
            }
        }

        let registry = build_registry(&nl, self.netlist_dir.as_ref())?;
        let opts = build_sim_options(&nl, kwargs)?;

        let analysis_lc = analysis.to_lowercase();
        // "dc" with a src kwarg is a sweep; without one it's an op-point alias.
        let is_dc_sweep = analysis_lc == "dc_sweep"
            || (analysis_lc == "dc"
                && kwargs
                    .and_then(|kw| kw.get_item("src").ok().flatten())
                    .is_some());

        if is_dc_sweep {
            let p = parse_dc_kwargs(kwargs)?;
            let nested_arg = p
                .nested
                .as_ref()
                .map(|(s, a, b, st)| (s.as_str(), *a, *b, *st));
            // GIL released around every solve: they touch no Python objects,
            // and holding it serialises Python threads that want to run
            // independent simulations concurrently (a finite-difference
            // Jacobian over bias points, a corner sweep).
            let result = py
                .allow_threads(|| {
                    dc_sweep_with_registry_opts(
                        &nl, &p.src, p.start, p.stop, p.step, nested_arg, &registry, &opts,
                    )
                })
                .map_err(sim_err)?;
            return Ok(SimResult {
                inner: SimResultInner::DcSweep(result),
                measurements: Vec::new(),
            });
        }

        match analysis_lc.as_str() {
            "op" | "dc" => {
                let result = py
                    .allow_threads(|| dc_op_nr_with_registry_opts(&nl, &registry, &opts))
                    .map_err(sim_err)?;
                Ok(SimResult {
                    inner: SimResultInner::Dc(result),
                    measurements: Vec::new(),
                })
            }
            "tran" | "transient" => {
                let (stop, step) = parse_tran_kwargs(kwargs)?;
                let result = py
                    .allow_threads(|| {
                        if opts.variable_step {
                            tran_nr_with_registry_var_opts(&nl, step, stop, &registry, &opts)
                        } else {
                            tran_nr_with_registry_opts(&nl, step, stop, &registry, &opts)
                        }
                    })
                    .map_err(sim_err)?;
                let measurements = evaluate_measurements(&nl.measurements, &result)
                    .into_iter()
                    .map(|m| (m.name, m.value))
                    .collect();
                Ok(SimResult {
                    inner: SimResultInner::Tran(result),
                    measurements,
                })
            }
            "ac" => {
                let (freqs, src) = parse_ac_kwargs(kwargs)?;
                let result = py
                    .allow_threads(|| {
                        ac_analysis_opts(&nl, &freqs, src.as_deref(), &registry, &opts)
                    })
                    .map_err(sim_err)?;
                Ok(SimResult {
                    inner: SimResultInner::Ac(result),
                    measurements: Vec::new(),
                })
            }
            "noise" => {
                let (freqs, out_pos, out_neg, input_src) = parse_noise_kwargs(kwargs)?;
                let result = py
                    .allow_threads(|| {
                        fairchild_core::noise_analysis(
                            &nl, &freqs, &out_pos, &out_neg, &input_src, &registry, &opts,
                        )
                    })
                    .map_err(sim_err)?;
                Ok(SimResult {
                    inner: SimResultInner::Noise(result),
                    measurements: Vec::new(),
                })
            }
            other => Err(PyRuntimeError::new_err(format!(
                "unknown analysis '{}'; use 'op', 'tran', 'ac', 'noise', or 'dc_sweep'",
                other
            ))),
        }
    }

    /// Run a transient that can be differentiated, and return the run.
    ///
    /// `probes` maps a name you choose to what it reads:
    ///
    /// ```python
    /// run = ckt.tran_adjoint(step=10e-12, stop=2e-9,
    ///                        probes={"v": "pout",                  # node voltage
    ///                                "p": ("power", "out0", 0)})   # optical power, W
    /// y = run.probes["v"]                                # (K,) numpy
    /// g = run.backward({"v": 2 * (y - target)},          # dL/dy per timepoint
    ///                  ["Xmzm.V_pi", "Rd.r"])            # -> (2,) numpy
    /// ```
    ///
    /// `params` optionally sets `"element.param"` values for this run only,
    /// which is what an optimiser's inner loop wants — it leaves the `set_param`
    /// overrides and the netlist on disk alone.
    ///
    /// Everything else — solver options, integration method — is the same
    /// kwargs as `run()`.  The step is fixed: the adjoint drives the fixed-step
    /// integrator, so `variable_step` does not apply.
    #[pyo3(signature = (probes, **kwargs))]
    pub fn tran_adjoint(
        &self,
        probes: &Bound<'_, PyDict>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<TranAdjointRun> {
        let netlist = self
            .netlist
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("no netlist loaded; call load() first"))?;

        let (stop, step) = parse_tran_kwargs(kwargs)?;

        let mut nl = netlist.clone();
        apply_overrides(&mut nl, &self.overrides);
        apply_source_overrides(&mut nl, &self.source_overrides);
        if let Some(kw) = kwargs {
            if let Some(p) = kw.get_item("params")? {
                let per_run: HashMap<String, f64> = p.extract()?;
                let lowered = per_run
                    .into_iter()
                    .map(|(k, v)| (k.to_lowercase(), v))
                    .collect();
                apply_overrides(&mut nl, &lowered);
            }
        }

        let registry = build_registry(&nl, self.netlist_dir.as_ref())?;
        let opts = build_sim_options(&nl, kwargs)?;

        let declared: Vec<(String, Output)> = probes
            .iter()
            .map(|(k, v)| Ok((k.extract::<String>()?, parse_probe(&v)?)))
            .collect::<PyResult<_>>()?;
        if declared.is_empty() {
            return Err(PyRuntimeError::new_err(
                "tran_adjoint needs at least one probe; nothing else can be differentiated",
            ));
        }

        let inner = TranAdjoint::run(&nl, &registry, &opts, step, stop).map_err(sim_err)?;
        Ok(TranAdjointRun {
            inner,
            registry,
            probes: declared,
        })
    }

    /// Run a parametric sweep over scalar element parameters.
    ///
    /// Calls `run(analysis, **kwargs)` once per value in `values`, each time
    /// setting `param` (e.g. `"Xlaser.power_mW"`) to that value.
    ///
    /// Returns a list of `SimResult` objects, one per sweep point.
    #[pyo3(signature = (param, values, analysis, **kwargs))]
    #[allow(clippy::useless_conversion)]
    pub fn sweep(
        &self,
        param: &str,
        values: Vec<f64>,
        analysis: &str,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Vec<SimResult>> {
        let netlist = self
            .netlist
            .as_ref()
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
                    SimResult {
                        inner: SimResultInner::Dc(r),
                        measurements: Vec::new(),
                    }
                }
                "tran" | "transient" => {
                    let (stop, step) = parse_tran_kwargs(kwargs)?;
                    let r = if opts.variable_step {
                        tran_nr_with_registry_var_opts(&nl, step, stop, &registry, &opts)
                    } else {
                        tran_nr_with_registry_opts(&nl, step, stop, &registry, &opts)
                    }
                    .map_err(sim_err)?;
                    let measurements = evaluate_measurements(&nl.measurements, &r)
                        .into_iter()
                        .map(|m| (m.name, m.value))
                        .collect();
                    SimResult {
                        inner: SimResultInner::Tran(r),
                        measurements,
                    }
                }
                "ac" => {
                    let (freqs, src) = parse_ac_kwargs(kwargs)?;
                    let r = ac_analysis_opts(&nl, &freqs, src.as_deref(), &registry, &opts)
                        .map_err(sim_err)?;
                    SimResult {
                        inner: SimResultInner::Ac(r),
                        measurements: Vec::new(),
                    }
                }
                other => {
                    return Err(PyRuntimeError::new_err(format!(
                        "unknown analysis '{}'; use 'op', 'tran', or 'ac'",
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
// Helper: apply source (WaveformSource) overrides
// ---------------------------------------------------------------------------

fn apply_source_overrides(netlist: &mut Netlist, overrides: &HashMap<String, Vec<(f64, f64)>>) {
    for (name_lc, points) in overrides {
        fairchild_core::set_source_pwl(netlist, name_lc, points.clone());
    }
}

// ---------------------------------------------------------------------------
// Helper: build a DeviceRegistry
// ---------------------------------------------------------------------------

fn build_registry(netlist: &Netlist, netlist_dir: Option<&PathBuf>) -> PyResult<DeviceRegistry> {
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&netlist.models);

    #[cfg(feature = "osdi")]
    for osdi_path in &netlist.osdi_paths {
        let path = if std::path::Path::new(osdi_path).is_absolute() {
            PathBuf::from(osdi_path)
        } else if let Some(dir) = netlist_dir {
            dir.join(osdi_path)
        } else {
            PathBuf::from(osdi_path)
        };

        let lib = unsafe { OsdiLibrary::open(&path) }.map_err(|e| {
            PyRuntimeError::new_err(format!(
                "failed to load OSDI library '{}': {e}",
                path.display()
            ))
        })?;
        let lib = Arc::new(lib);
        lib.register_into(&mut registry);
    }
    // `.model <card> <module> (...)` cards naming a descriptor we just loaded.
    // After the libraries, so the descriptors exist to alias.
    #[cfg(feature = "osdi")]
    registry.register_loaded_model_cards(&netlist.models);

    #[cfg(not(feature = "osdi"))]
    if !netlist.osdi_paths.is_empty() {
        return Err(PyRuntimeError::new_err(format!(
            "netlist references .osdi files but this build was compiled without OSDI support; \
             rebuild with --features osdi or remove the .osdi references"
        )));
    }

    Ok(registry)
}

// ---------------------------------------------------------------------------
// Helper: parse tran kwargs
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

    let stop =
        stop.ok_or_else(|| PyRuntimeError::new_err("tran requires 'stop' kwarg (seconds)"))?;
    let step =
        step.ok_or_else(|| PyRuntimeError::new_err("tran requires 'step' kwarg (seconds)"))?;
    Ok((stop, step))
}

// ---------------------------------------------------------------------------
// Helper: build SimOptions from netlist + Python kwargs
// ---------------------------------------------------------------------------

/// Build a `SimOptions` by starting from the netlist's `.options` directives
/// and overlaying any Python kwargs the user passed.
///
/// Analysis-specific kwargs (`stop`, `step`, `fstart`, `src`, etc.) are
/// silently skipped — they are consumed by the `parse_*_kwargs` helpers.
/// Every other kwarg is forwarded to `SimOptions::set`; an unrecognised or
/// unavailable key raises `PyRuntimeError` so misspellings surface immediately.
fn build_sim_options(
    netlist: &Netlist,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<SimOptions> {
    let mut opts = SimOptions::from_netlist(netlist);
    if let Some(kw) = kwargs {
        // Keywords consumed by parse_tran/ac/dc/noise_kwargs or run() itself.
        // Any kwarg NOT in this list is forwarded to SimOptions::set().
        const SKIP: &[&str] = &[
            "alter", // .alter block selector
            "stop",
            "step", // tran / dc sweep
            "src",
            "start",
            "src2", // dc sweep
            "start2",
            "stop2",
            "step2", // dc sweep nested
            "fstart",
            "fstop",
            "points",
            "variation", // ac / noise
            "out_pos",
            "out_neg",
            "out",    // noise
            "params", // tran_adjoint per-run parameter overrides
        ];
        for (k, v) in kw.iter() {
            let key: String = k.extract()?;
            if SKIP.contains(&key.as_str()) {
                continue;
            }
            // Accept bool before numeric to avoid True → 1.0 path.
            let value_str: String = if let Ok(b) = v.extract::<bool>() {
                if b {
                    "1".into()
                } else {
                    "0".into()
                }
            } else if let Ok(s) = v.extract::<String>() {
                s
            } else if let Ok(i) = v.extract::<i64>() {
                i.to_string()
            } else if let Ok(f) = v.extract::<f64>() {
                format!("{f:e}")
            } else {
                return Err(PyRuntimeError::new_err(format!(
                    "kwarg '{key}': expected number, string, or bool"
                )));
            };
            if !opts.set(&key, &value_str) {
                if key == "solver" {
                    return Err(PyRuntimeError::new_err(format!(
                        "unsupported solver '{value_str}' — valid values: dense, sparse, auto; \
                         use 'klu' only when fairchild is built with `--features klu` \
                         and SuiteSparse is installed (`brew install suite-sparse`)"
                    )));
                }
                return Err(PyRuntimeError::new_err(format!(
                    "unrecognised option '{key}={value_str}'"
                )));
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
    src: String,
    start: f64,
    stop: f64,
    step: f64,
    nested: Option<(String, f64, f64, f64)>,
}

fn parse_dc_kwargs(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<DcKwargs> {
    let kw = kwargs
        .ok_or_else(|| PyRuntimeError::new_err("dc requires src, start, stop, step kwargs"))?;

    let get = |k: &str| -> PyResult<_> {
        kw.get_item(k)?
            .ok_or_else(|| PyRuntimeError::new_err(format!("dc requires '{k}' kwarg")))
    };

    let src = get("src")?.extract::<String>()?;
    let start = get("start")?.extract::<f64>()?;
    let stop = get("stop")?.extract::<f64>()?;
    let step = get("step")?.extract::<f64>()?;

    // Optional nested second sweep: src2, start2, stop2, step2 (all four required if any).
    let nested = match kw.get_item("src2")? {
        Some(v) => {
            let src2 = v.extract::<String>()?;
            let start2 = get("start2")?.extract::<f64>()?;
            let stop2 = get("stop2")?.extract::<f64>()?;
            let step2 = get("step2")?.extract::<f64>()?;
            Some((src2, start2, stop2, step2))
        }
        None => None,
    };

    Ok(DcKwargs {
        src,
        start,
        stop,
        step,
        nested,
    })
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
        if let Some(v) = kw.get_item("fstart")? {
            fstart = Some(v.extract::<f64>()?);
        }
        if let Some(v) = kw.get_item("fstop")? {
            fstop = Some(v.extract::<f64>()?);
        }
        if let Some(v) = kw.get_item("points")? {
            points = v.extract::<usize>()?;
        }
        if let Some(v) = kw.get_item("src")? {
            let s: Option<String> = v.extract()?;
            src = s;
        }
        if let Some(v) = kw.get_item("variation")? {
            let var: String = v.extract()?;
            variation = match var.to_lowercase().as_str() {
                "dec" => AcVariation::Dec,
                "oct" => AcVariation::Oct,
                "lin" => AcVariation::Lin,
                other => {
                    return Err(PyRuntimeError::new_err(format!(
                        "unknown AC variation '{other}'; use 'dec', 'oct', or 'lin'"
                    )))
                }
            };
        }
    }

    let fstart =
        fstart.ok_or_else(|| PyRuntimeError::new_err("ac requires 'fstart' kwarg (Hz)"))?;
    let fstop = fstop.ok_or_else(|| PyRuntimeError::new_err("ac requires 'fstop' kwarg (Hz)"))?;

    let freqs = match variation {
        AcVariation::Dec => freq_decade(fstart, fstop, points),
        AcVariation::Oct => freq_oct(fstart, fstop, points),
        AcVariation::Lin => freq_linear(fstart, fstop, points),
    };

    Ok((freqs, src))
}

fn parse_noise_kwargs(
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<(Vec<f64>, String, String, String)> {
    let kw = kwargs.ok_or_else(|| {
        PyRuntimeError::new_err(
            "noise requires kwargs: out (or out_pos+out_neg), src, fstart, fstop",
        )
    })?;

    let mut out_pos: Option<String> = None;
    let mut out_neg: String = "0".to_string();
    if let Some(v) = kw.get_item("out_pos")? {
        out_pos = Some(v.extract()?);
    }
    if let Some(v) = kw.get_item("out_neg")? {
        out_neg = v.extract()?;
    }
    if out_pos.is_none() {
        if let Some(v) = kw.get_item("out")? {
            out_pos = Some(v.extract()?);
        }
    }
    let out_pos = out_pos.ok_or_else(|| {
        PyRuntimeError::new_err("noise: missing 'out' (or 'out_pos') kwarg — the observation node")
    })?;

    let input_src: String = kw
        .get_item("src")?
        .ok_or_else(|| PyRuntimeError::new_err("noise: missing 'src' kwarg"))?
        .extract()?;

    let mut fstart: Option<f64> = None;
    let mut fstop: Option<f64> = None;
    let mut points: usize = 20;
    let mut variation = AcVariation::Dec;
    if let Some(v) = kw.get_item("fstart")? {
        fstart = Some(v.extract()?);
    }
    if let Some(v) = kw.get_item("fstop")? {
        fstop = Some(v.extract()?);
    }
    if let Some(v) = kw.get_item("points")? {
        points = v.extract()?;
    }
    if let Some(v) = kw.get_item("variation")? {
        let var: String = v.extract()?;
        variation = match var.to_lowercase().as_str() {
            "dec" => AcVariation::Dec,
            "oct" => AcVariation::Oct,
            "lin" => AcVariation::Lin,
            other => {
                return Err(PyRuntimeError::new_err(format!(
                    "unknown variation '{other}'; use 'dec', 'oct', or 'lin'"
                )))
            }
        };
    }
    let fstart = fstart.ok_or_else(|| PyRuntimeError::new_err("noise: missing 'fstart' kwarg"))?;
    let fstop = fstop.ok_or_else(|| PyRuntimeError::new_err("noise: missing 'fstop' kwarg"))?;

    let freqs = match variation {
        AcVariation::Dec => freq_decade(fstart, fstop, points),
        AcVariation::Oct => freq_oct(fstart, fstop, points),
        AcVariation::Lin => freq_linear(fstart, fstop, points),
    };

    Ok((
        freqs,
        out_pos.to_lowercase(),
        out_neg.to_lowercase(),
        input_src.to_lowercase(),
    ))
}

// ---------------------------------------------------------------------------
// Transient adjoint
// ---------------------------------------------------------------------------

/// Turn a Python probe spec into an [`Output`].
///
/// A bare string is the common case — a node voltage — and the tuple forms
/// spell the rest: `("node", n)`, `("current", vsrc)`, `("power", net, ch)`.
fn parse_probe(spec: &Bound<'_, PyAny>) -> PyResult<Output> {
    const FORMS: &str = "a probe is a node name, or a tuple: ('node', name) | \
                         ('current', vsrc) | ('power', net, channel)";
    if let Ok(node) = spec.extract::<String>() {
        return Ok(Output::NodeVoltage(node));
    }
    let items: Vec<Bound<'_, PyAny>> =
        spec.extract().map_err(|_| PyRuntimeError::new_err(FORMS))?;
    let kind: String = items
        .first()
        .ok_or_else(|| PyRuntimeError::new_err(FORMS))?
        .extract()
        .map_err(|_| PyRuntimeError::new_err(FORMS))?;
    match (kind.as_str(), items.len()) {
        ("node", 2) => Ok(Output::NodeVoltage(items[1].extract()?)),
        ("current", 2) => Ok(Output::BranchCurrent(items[1].extract()?)),
        ("power", 3) => Ok(Output::OpticalPower {
            net: items[1].extract()?,
            channel: items[2].extract()?,
        }),
        _ => Err(PyRuntimeError::new_err(format!(
            "probe spec ('{kind}', ...) with {} entries is not one of the known forms — {FORMS}",
            items.len()
        ))),
    }
}

/// Split `"Xmzm.V_pi"` into the element and parameter a [`ParamRef`] needs.
fn parse_param(name: &str) -> PyResult<ParamRef> {
    let (element, param) = name.split_once('.').ok_or_else(|| {
        PyRuntimeError::new_err(format!(
            "parameter '{name}' must be written 'element.param', e.g. 'Xmzm.V_pi'"
        ))
    })?;
    Ok(ParamRef::new(element, param))
}

/// A transient run that can be differentiated.
///
/// Holds the trajectory plus the per-timestep state the backward pass needs, so
/// any number of objectives can be differentiated against one forward run.
///
/// Not thread-safe: it owns a device registry, which may hold `dlopen`ed OSDI
/// libraries bound to the thread that loaded them.
#[pyclass(unsendable)]
pub struct TranAdjointRun {
    inner: TranAdjoint,
    registry: DeviceRegistry,
    /// Probe name → what it reads, in the order they were declared.
    probes: Vec<(String, Output)>,
}

#[pymethods]
impl TranAdjointRun {
    /// Accepted timepoints, in seconds.  `time[0]` is always 0.
    #[getter]
    fn time<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice_bound(py, self.inner.time())
    }

    /// Each declared probe's waveform, as `{name: array of len(time)}`.
    #[getter]
    fn probes<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let out = PyDict::new_bound(py);
        for (name, probe) in &self.probes {
            let signal = self.inner.signal(probe).map_err(sim_err)?;
            out.set_item(name, PyArray1::from_vec_bound(py, signal))?;
        }
        Ok(out)
    }

    /// `dL/dp` for each named parameter, given `dL/d(probe)` at every timepoint.
    ///
    /// `cotangents` maps a probe name to an array of `len(time)` — the
    /// derivative of your loss with respect to that probe's value at each
    /// timepoint.  Build the loss from `probes` in numpy however you like; this
    /// only needs its derivative.  Probes you leave out contribute nothing.
    ///
    /// `params` are `"element.param"` strings, the same spelling `set_param`
    /// takes.  Returns one gradient per parameter, in the order given.
    ///
    /// Raises if a parameter reaches nothing in the equations — a silent zero
    /// there is indistinguishable from a real insensitivity, and would stall an
    /// optimiser at a point that looks stationary.  Warns if the finite
    /// difference behind `∂G/∂p` could not be made accurate, rather than
    /// returning a number that looks as good as the rest.
    fn backward<'py>(
        &self,
        py: Python<'py>,
        cotangents: &Bound<'py, PyDict>,
        params: Vec<String>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let n_t = self.inner.time().len();
        let n_x = self.inner.topology().size;
        let mut seeds = vec![vec![0.0; n_x]; n_t];

        for (key, value) in cotangents.iter() {
            let name: String = key.extract()?;
            let probe = self
                .probes
                .iter()
                .find(|(n, _)| *n == name)
                .map(|(_, p)| p)
                .ok_or_else(|| {
                    PyRuntimeError::new_err(format!(
                        "no probe named '{name}'; declared probes are {:?}",
                        self.probes.iter().map(|(n, _)| n).collect::<Vec<_>>()
                    ))
                })?;
            let weights: Vec<f64> = value
                .extract::<PyReadonlyArray1<f64>>()?
                .as_array()
                .to_vec();
            if weights.len() != n_t {
                return Err(PyRuntimeError::new_err(format!(
                    "cotangent for '{name}' has length {} but the run has {n_t} timepoints",
                    weights.len()
                )));
            }
            // `weighted` evaluates ∂probe/∂x at each x_k and scales it; several
            // probes just add, which is the chain rule and nothing more.
            let (_, contribution) = self.inner.weighted(probe, &weights).map_err(sim_err)?;
            for (seed, add) in seeds.iter_mut().zip(contribution.iter()) {
                for (s, a) in seed.iter_mut().zip(add.iter()) {
                    *s += a;
                }
            }
        }

        let refs: Vec<ParamRef> = params
            .iter()
            .map(|p| parse_param(p))
            .collect::<PyResult<_>>()?;
        let s = self
            .inner
            .gradient(&self.registry, &seeds, &refs)
            .map_err(sim_err)?;

        let unreached: Vec<&String> = params
            .iter()
            .zip(s.reached.iter())
            .filter(|(_, ok)| !**ok)
            .map(|(p, _)| p)
            .collect();
        if !unreached.is_empty() {
            return Err(PyRuntimeError::new_err(format!(
                "these parameters reach nothing in the equations, so their gradient is a \
                 placeholder zero rather than a computed one: {unreached:?}.  Either the name \
                 is wrong, or the model does not accept that parameter at runtime (see \
                 docs/model_status.md)"
            )));
        }

        let shaky: Vec<(&String, f64)> = params
            .iter()
            .zip(s.fd_error.iter())
            .filter(|(_, e)| **e > 1e-3)
            .map(|(p, e)| (p, *e))
            .collect();
        if !shaky.is_empty() {
            let warnings = py.import_bound("warnings")?;
            warnings.call_method1(
                "warn",
                (format!(
                    "fairchild: the gradient for {shaky:?} could not be resolved to better than \
                     the relative error shown.  The parameter's scale and the objective's \
                     disagree; pass a step explicitly if you know a better one."
                ),),
            )?;
        }

        Ok(PyArray1::from_vec_bound(py, s.grad))
    }
}

// ---------------------------------------------------------------------------
// Module entry point
// ---------------------------------------------------------------------------

#[pymodule]
fn fairchild(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Circuit>()?;
    m.add_class::<SimResult>()?;
    m.add_class::<WaveformSource>()?;
    m.add_class::<TranAdjointRun>()?;
    Ok(())
}
