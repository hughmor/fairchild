//! PyO3 bindings for the fairchild electro-optic circuit simulator.
//!
//! Exposes `Circuit`, `SimResult`, and `WaveformSource` to Python.

// pyo3's #[pymethods] macro expansion triggers this lint on PyResult<T> return
// types as a false positive; suppress it for the whole crate.
#![allow(clippy::useless_conversion)]

use std::collections::HashMap;
use std::path::PathBuf;

use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::{PyComplex, PyDict};
use rayon::prelude::*;

use fairchild_core::adjoint::dc_sensitivity;
use fairchild_core::adjoint_ac::{AcAdjoint, AcOutput};
use fairchild_core::{
    ac_analysis_opts, dc_op_nr_with_registry_opts, dc_sweep_with_registry_opts,
    evaluate_measurements, freq_decade, freq_linear, freq_oct, tran_nr_with_registry_opts,
    tran_nr_with_registry_var_opts, AcResult, DcSweepResult, DeviceRegistry, NrResult, Output,
    ParamRef, SimError, SimOptions, TranAdjoint, TranResult,
};
use fairchild_parser::{
    parse_spice_file_with_arity, parse_spice_with_arity, AcVariation, Analysis, ArityOracle,
    Netlist, OutVar, ParamName, PermissiveArity, PzDrive, PzWant,
};

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

/// Apply `"element.param" -> value` overrides, refusing any that matched nothing.
///
/// Both halves of that refusal used to be silent: a key with no dot was
/// `continue`d, and `set_element_param`'s `bool` was discarded, so a typo in
/// either half returned the un-overridden answer and reported success. That is
/// the shape of bug this codebase exists to not ship — an optimiser driving a
/// no-op sees a flat objective, converges instantly, and hands back its own
/// starting values as a fit (fairchild issue #56).
///
/// The key is deliberately an *element* parameter and not a `.param` name. A
/// `.param` is substituted textually during parse, before model cards and
/// element values exist, so there is nothing left to reach by the time a
/// netlist is in hand; overriding one means re-parsing, which is a different
/// operation with a different cost. The error says so rather than guessing.
fn apply_overrides(netlist: &mut Netlist, overrides: &HashMap<String, f64>) -> PyResult<()> {
    for (key, &value) in overrides {
        let hit = match key.find('.') {
            Some(dot) => {
                fairchild_core::set_element_param(netlist, &key[..dot], &key[dot + 1..], value)
            }
            None => false,
        };
        if !hit {
            return Err(PyRuntimeError::new_err(format!(
                "no element parameter '{key}' to override. Expected 'element.param', \
                 e.g. 'R1.value' or 'Xmzm.v_pi', naming an element that exists in the \
                 loaded netlist. If '{key}' is a `.param` in the deck, it cannot be \
                 overridden here — a `.param` is substituted at parse time, so change \
                 the deck text and re-load instead."
            )));
        }
    }
    Ok(())
}

/// The `params=` kwarg: `"element.param"` overrides for this call only.
///
/// Shared by every entry point that takes it, so one of them cannot quietly
/// gain or lose the behaviour — which is exactly how `run()` came to accept the
/// kwarg and drop it.
fn apply_param_kwarg(netlist: &mut Netlist, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
    let Some(kw) = kwargs else { return Ok(()) };
    let Some(p) = kw.get_item("params")? else {
        return Ok(());
    };
    // Not lowercased: `set_element_param` folds case on both halves itself, and
    // the key as written is what belongs in the error when it matches nothing.
    let per_run: HashMap<String, f64> = p.extract()?;
    apply_overrides(netlist, &per_run)
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
                // λ labels are probeable by name, so they are listed (#71).
                sigs.extend(
                    r.topo
                        .lambda_signals()
                        .iter()
                        .map(|(n, _)| format!("V({n})")),
                );
                sigs.extend(r.topo.vsrc_index.keys().map(|n| format!("I({n})")));
                sigs
            }
            SimResultInner::Tran(r) => {
                let mut sigs: Vec<String> =
                    r.node_voltages.keys().map(|n| format!("V({n})")).collect();
                sigs.extend(r.lambda.keys().map(|n| format!("V({n})")));
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
            if let Some(series) = r.node_voltages.get(node) {
                return Ok(PyArray1::from_vec_bound(py, series.clone()));
            }
            // A λ label is constant over the run; materialise it as a series so
            // the caller's array arithmetic works the same as for any node.
            if let Some(&wl) = r.lambda.get(node) {
                return Ok(PyArray1::from_vec_bound(py, vec![wl; r.time.len()]));
            }
            Err(PyRuntimeError::new_err(format!("unknown node '{node}'")))
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

/// Rust-side helpers — deliberately outside `#[pymethods]`, which would try to
/// expose them to Python and cannot, since they traffic in `Netlist` and
/// `SimOptions`.
impl Circuit {
    /// The netlist, registry and options every analysis entry point starts
    /// from: the loaded deck with `set_param` overrides, injected sources, a
    /// `params=` kwarg and a chosen `.alter` block applied, in that order.
    ///
    /// One copy, shared by `run` and the three small-signal reports.  Four
    /// copies of "apply the overrides, then the params, then the alter" is four
    /// chances for one entry point to honour something the others do not, which
    /// is how `params=` came to be silently dropped once already (#56).
    fn prepare(
        &self,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<(Netlist, DeviceRegistry, SimOptions)> {
        let netlist = self
            .netlist
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("no netlist loaded; call load() first"))?;

        let mut nl = netlist.clone();
        apply_overrides(&mut nl, &self.overrides)?;
        apply_source_overrides(&mut nl, &self.source_overrides);
        // `params=` means here what it means on `tran_adjoint`: per-call
        // `element.param` overrides that leave `set_param` and the deck alone.
        // It used to be accepted and dropped, which made a fitting loop return
        // its initial guess and call it converged (fairchild issue #56).
        apply_param_kwarg(&mut nl, kwargs)?;

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
        Ok((nl, registry, opts))
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
        let dir = p.parent().map(|d| d.to_path_buf());
        // See `two_pass` — WDM dispatch is the registry's answer, and the
        // registry comes from a parse.
        let netlist = two_pass(dir.as_ref(), |oracle| {
            parse_spice_file_with_arity(&p, oracle).map_err(parse_err)
        })?;
        self.netlist_dir = dir;
        self.netlist = Some(netlist);
        Ok(())
    }

    /// Load a netlist from a SPICE string.
    #[allow(clippy::useless_conversion)]
    pub fn load_str(&mut self, src: &str) -> PyResult<()> {
        let netlist = two_pass(None, |oracle| {
            parse_spice_with_arity(src, oracle).map_err(parse_err)
        })?;
        self.netlist_dir = None;
        self.netlist = Some(netlist);
        Ok(())
    }

    /// Override a parameter on an element before the next `run()`.
    ///
    /// The name is not checked here — it is checked at the next `run()`, which
    /// raises if it matches no element in the loaded netlist. A typo used to be
    /// discarded, so the run returned the un-overridden answer and reported
    /// success (fairchild issue #56).
    ///
    /// For a one-call override that does not persist, pass
    /// `run(..., params={"element.param": value})` instead.
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

    /// The analyses the deck declares, in deck order, as a list of dicts.
    ///
    /// A deck *declares* what could be run; it never runs anything on its own
    /// here. `run(kind)` with no analysis kwargs adopts the matching card whole;
    /// this is how you see what that card says, and what to pass if you want to
    /// override it.
    ///
    /// ```python
    /// ckt.analyses
    /// # [{'kind': 'op'},
    /// #  {'kind': 'tran', 'step': 1e-12, 'stop': 5e-09, 'tstart': 0.0,
    /// #   'tmax': 1e-13, 'uic': False}]
    /// ```
    #[getter]
    pub fn analyses<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
        let Some(netlist) = self.netlist.as_ref() else {
            return Ok(Vec::new());
        };
        netlist
            .analyses
            .iter()
            .map(|a| {
                let d = PyDict::new_bound(py);
                match a {
                    Analysis::Op => d.set_item("kind", "op")?,
                    Analysis::Tran {
                        step,
                        stop,
                        tstart,
                        tmax,
                        uic,
                    } => {
                        d.set_item("kind", "tran")?;
                        d.set_item("step", step)?;
                        d.set_item("stop", stop)?;
                        d.set_item("tstart", tstart)?;
                        d.set_item("tmax", tmax)?;
                        d.set_item("uic", uic)?;
                    }
                    Analysis::Ac {
                        variation,
                        points,
                        fstart,
                        fstop,
                    } => {
                        d.set_item("kind", "ac")?;
                        d.set_item("variation", variation_name(*variation))?;
                        d.set_item("points", points)?;
                        d.set_item("fstart", fstart)?;
                        d.set_item("fstop", fstop)?;
                    }
                    Analysis::Dc {
                        src,
                        start,
                        stop,
                        step,
                        nested,
                    } => {
                        d.set_item("kind", "dc")?;
                        d.set_item("src", src)?;
                        d.set_item("start", start)?;
                        d.set_item("stop", stop)?;
                        d.set_item("step", step)?;
                        if let Some(n) = nested {
                            d.set_item("src2", &n.src)?;
                            d.set_item("start2", n.start)?;
                            d.set_item("stop2", n.stop)?;
                            d.set_item("step2", n.step)?;
                        }
                    }
                    Analysis::Noise {
                        out_pos,
                        out_neg,
                        input_src,
                        variation,
                        points,
                        fstart,
                        fstop,
                    } => {
                        d.set_item("kind", "noise")?;
                        d.set_item("out_pos", out_pos)?;
                        d.set_item("out_neg", out_neg)?;
                        d.set_item("src", input_src)?;
                        d.set_item("variation", variation_name(*variation))?;
                        d.set_item("points", points)?;
                        d.set_item("fstart", fstart)?;
                        d.set_item("fstop", fstop)?;
                    }
                    Analysis::Tf { out, input_src } => {
                        d.set_item("kind", "tf")?;
                        d.set_item("out", outvar_str(out))?;
                        d.set_item("src", input_src)?;
                    }
                    Analysis::Sens { out, params } => {
                        d.set_item("kind", "sens")?;
                        d.set_item("out", outvar_str(out))?;
                        let names: Vec<String> = params.iter().map(param_str).collect();
                        d.set_item("params", names)?;
                    }
                    Analysis::Pz {
                        in_pos,
                        in_neg,
                        out_pos,
                        out_neg,
                        drive,
                        want,
                    } => {
                        d.set_item("kind", "pz")?;
                        d.set_item("in_pos", in_pos)?;
                        d.set_item("in_neg", in_neg)?;
                        d.set_item("out_pos", out_pos)?;
                        d.set_item("out_neg", out_neg)?;
                        d.set_item(
                            "drive",
                            match drive {
                                PzDrive::Vol => "vol",
                                PzDrive::Cur => "cur",
                            },
                        )?;
                        d.set_item(
                            "want",
                            match want {
                                PzWant::Poles => "pol",
                                PzWant::Zeros => "zer",
                                PzWant::Both => "pz",
                            },
                        )?;
                    }
                }
                Ok(d)
            })
            .collect()
    }

    /// Every corner the deck declares, as `{'alter': …, 'temp_c': …}`.
    ///
    /// The `.alter` × `.temp` grid, base run first. Paired with `run()`'s
    /// `alter=` and `temp=` kwargs this is the loop-it-yourself route, for when
    /// you want to do something between corners rather than collect them all:
    ///
    /// ```python
    /// for c in ckt.corners:
    ///     r = ckt.run("tran", alter=c["alter"], temp=c["temp_c"])
    ///     print(c["alter"], c["temp_c"], r["V(out)"].max())
    /// ```
    ///
    /// Always at least one entry (`alter = "base"`), so a deck declaring no
    /// corners needs no special case. `run_all()` covers the same grid when you
    /// just want every result.
    #[getter]
    pub fn corners<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
        let Some(netlist) = self.netlist.as_ref() else {
            return Ok(Vec::new());
        };
        let opts = SimOptions::from_netlist(netlist);
        fairchild_core::expand_corners(netlist, &opts)
            .corners
            .iter()
            .map(|c| {
                let d = PyDict::new_bound(py);
                d.set_item("alter", &c.alter_label)?;
                d.set_item("temp_c", c.temp_c())?;
                Ok(d)
            })
            .collect()
    }

    /// Run every analysis the deck declares, at every corner it declares.
    ///
    /// This is what the CLI does with the same file. `run()` runs one named
    /// analysis at one corner; this runs the whole deck — the `.alter` × `.temp`
    /// grid crossed with the deck's analysis list, in deck order, corners in
    /// parallel.
    ///
    /// ```python
    /// for row in ckt.run_all():
    ///     print(row["alter"], row["temp_c"], row["kind"])
    ///     row["result"]        # SimResult, or a report dict for tf/sens/pz
    /// ```
    ///
    /// Each row carries the corner it came from because the results are
    /// otherwise indistinguishable — two `.tran` rows from a two-corner deck
    /// are the same analysis at different temperatures, and nothing in a
    /// `SimResult` says which.
    ///
    /// A deck declaring no analyses returns an empty list rather than guessing
    /// at one: `run_all` runs what the deck asked for, and a deck that asked
    /// for nothing gets nothing. Solver-option kwargs (`reltol`, `method`, …)
    /// apply to every corner; `alter=` and `temp=` do not belong here — pick a
    /// single corner with `run()` instead.
    #[pyo3(signature = (**kwargs))]
    pub fn run_all<'py>(
        &self,
        py: Python<'py>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Vec<Bound<'py, PyDict>>> {
        for reject in ["alter", "temp"] {
            if kwargs.is_some_and(|kw| kw.get_item(reject).ok().flatten().is_some()) {
                return Err(PyRuntimeError::new_err(format!(
                    "run_all() runs every corner, so {reject}= would contradict it. \
                     Use run(analysis, {reject}=…) for a single corner, or ckt.corners \
                     to loop over them yourself"
                )));
            }
        }

        let (nl, _registry, opts) = self.prepare(kwargs)?;
        let grid = fairchild_core::expand_corners(&nl, &opts);
        let netlist_dir = self.netlist_dir.clone();

        // Corners are independent simulations, so they fan out. The registry is
        // built inside the worker rather than shared: it holds factory closures
        // that are not `Sync`, and building one is cheap next to a solve.
        //
        // Everything here is plain Rust — no Python object is touched until the
        // fan-out has finished, which is what makes releasing the GIL sound.
        let per_corner: Vec<Result<Vec<(String, AnyResult)>, String>> = py.allow_threads(|| {
            grid.corners
                .par_iter()
                .map(|corner| run_one_corner(corner, netlist_dir.as_ref()))
                .collect()
        });

        let mut rows = Vec::new();
        for (corner, outcome) in grid.corners.iter().zip(per_corner) {
            let results = outcome.map_err(PyRuntimeError::new_err)?;
            for (kind, result) in results {
                let d = PyDict::new_bound(py);
                d.set_item("alter", &corner.alter_label)?;
                d.set_item("temp_c", corner.temp_c())?;
                d.set_item("kind", &kind)?;
                d.set_item("result", result.into_py_object(py)?)?;
                rows.push(d);
            }
        }
        Ok(rows)
    }

    /// Small-signal transfer function about the DC operating point.
    ///
    /// ```python
    /// ckt.tf()                          # the deck's .tf card, whole
    /// ckt.tf(out="v(out)", src="Vin")   # explicit; works with no card
    /// # {'gain': 0.75, 'r_in': 4000.0, 'r_out': 750.0, 'out_value': 0.75}
    /// ```
    ///
    /// `out` is `v(node)`, `v(node,ref)` or `i(vsrc)`; `src` names an
    /// independent V or I source.  Pass neither and the deck's `.tf` card is
    /// adopted whole; pass either and the card is not used at all — the same
    /// rule `run()` follows, so the numbers in a result always come from one
    /// place.
    #[pyo3(signature = (**kwargs))]
    pub fn tf<'py>(
        &self,
        py: Python<'py>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let (nl, registry, opts) = self.prepare(kwargs)?;
        let (out, src) = tf_args(&nl, kwargs)?;
        let r = py
            .allow_threads(|| fairchild_core::transfer_function(&nl, &registry, &opts, &out, &src))
            .map_err(sim_err)?;
        tf_dict(py, &r)
    }

    /// DC parameter sensitivity of one output, by the adjoint method.
    ///
    /// ```python
    /// ckt.sens()                                    # the deck's .sens card
    /// ckt.sens(out="v(out)")                        # every element value
    /// ckt.sens(out="v(out)", params=["r1", "m1.w"]) # named
    /// # [{'param': 'r1.value', 'nominal': 1000.0, 'sensitivity': -1.875e-4,
    /// #   'normalised': -0.1875, 'reached': True, 'fd_error': 3.6e-11}, …]
    /// ```
    ///
    /// **Check `reached`.** A parameter the adjoint could not perturb reports
    /// `sensitivity = 0.0` with `reached = False`, and a real insensitivity
    /// reports the same zero with `reached = True`.  Optimising against the
    /// first without looking stalls somewhere that looks like a stationary
    /// point.
    #[pyo3(signature = (**kwargs))]
    pub fn sens<'py>(
        &self,
        py: Python<'py>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Vec<Bound<'py, PyDict>>> {
        let (nl, registry, opts) = self.prepare(kwargs)?;
        let (out, params) = sens_args(&nl, kwargs)?;
        let r = py
            .allow_threads(|| fairchild_core::sensitivity(&nl, &registry, &opts, &out, &params))
            .map_err(sim_err)?;
        sens_rows(py, &r)
    }

    /// Poles and zeros of the small-signal transfer function, in rad/s.
    ///
    /// ```python
    /// ckt.pz()                                             # the deck's .pz card
    /// ckt.pz(in_pos="in", out_pos="out", drive="vol")      # explicit
    /// # {'poles': [(-50000.0+998749.2j), (-50000.0-998749.2j)],
    /// #  'zeros': [], 'infinite_poles': 3, 'infinite_zeros': 0}
    /// ```
    ///
    /// Roots come back as Python complex numbers in rad/s; divide by `2*pi` for
    /// Hz.  `in_neg` and `out_neg` default to ground, `drive` to `"vol"` and
    /// `want` to `"pz"`.  Refuses circuits past a dense-eigensolver size limit
    /// rather than running for an unbounded time.
    #[pyo3(signature = (**kwargs))]
    pub fn pz<'py>(
        &self,
        py: Python<'py>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let (nl, registry, opts) = self.prepare(kwargs)?;
        let (ip, ineg, op, oneg, drive, want) = pz_args(&nl, kwargs)?;
        let r = py
            .allow_threads(|| {
                fairchild_core::pole_zero(
                    &nl, &registry, &opts, &ip, &ineg, &op, &oneg, drive, want,
                )
            })
            .map_err(sim_err)?;
        pz_dict(py, &r)
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
    ///   analysis: `"op"` for DC, `"tran"` for transient, `"ac"` for AC sweep,
    ///   `"noise"`, or `"dc_sweep"`.  (`"dc"` is an alias for `"op"` unless a
    ///   `src` kwarg makes it a sweep; a deck's `.dc` card is reached through
    ///   `"dc_sweep"`.)
    ///
    ///   `"tf"`, `"sens"` and `"pz"` also work here and are the same call as
    ///   `ckt.tf()` / `ckt.sens()` / `ckt.pz()` — but they return a dict (or a
    ///   list of dicts), **not** a `SimResult`, because a report has no time or
    ///   frequency axis to index.
    ///
    ///   Pass no analysis parameters and the deck's matching card is adopted
    ///   **whole**: `run("tran")` takes `step`, `stop`, `tstart`, `tmax` and
    ///   `UIC` off the `.tran` line.  Pass any one of them and the card is not
    ///   used at all — a card is never half-applied, so the numbers in a run
    ///   always come from one place.  See `analyses` for what a deck declares.
    ///
    ///   For `"tran"`: `stop` (s) and `step` (s), or a `.tran` card.
    ///
    ///   For `"ac"`: `fstart` (Hz), `fstop` (Hz), `points` (int, default 20),
    ///   `variation` (`"dec"`, `"oct"`, `"lin"`, default `"dec"`), or a `.ac`
    ///   card; `src` (excitation source name, default `None` = first V source)
    ///   is not on the card and stays yours either way.
    ///
    ///   Solver options (apply to all analyses): `reltol`, `abstol`, `vntol`,
    ///   `vmax`, `gmin`, `itl1`, `itl4`, `maxstep`,
    ///   `method` (`"be"` | `"tr"` | `"gear"`), `uic`, `temp` (°C),
    ///   `solver` (`"dense"` | `"sparse"` | `"auto"` | `"klu"`),
    ///   `variable_step` (bool, enables LTE-controlled variable-step transient),
    ///   `waveguide_delay` (bool), `cond_estimate` (bool, prints κ(A)),
    ///   `equilibrate` (bool, two-sided matrix scaling before LU),
    ///   `lambda_center_nm`, `enable_bidirectional`, `sanity_check`, and more.
    ///   These overlay any `.options` directives from the netlist.
    ///   Unrecognised kwargs raise `RuntimeError` immediately.
    ///
    ///   `params={"element.param": value}` overrides element parameters for
    ///   this call only, leaving `set_param` and the deck untouched — the inner
    ///   loop of a sweep or a fit. A name that matches no element raises, so a
    ///   typo cannot quietly return the un-overridden answer. A deck `.param`
    ///   is *not* addressable this way: it is substituted at parse time, so
    ///   change the deck text and re-load.
    #[pyo3(signature = (analysis, **kwargs))]
    #[allow(clippy::useless_conversion)]
    pub fn run(
        &self,
        py: Python<'_>,
        analysis: &str,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<PyObject> {
        // The small-signal reports are reachable by either spelling, and are
        // *the same call* either way — `run` delegates rather than reimplements,
        // so the two can never come to disagree. What they do not do is arrive
        // wrapped in a `SimResult`: that class is arrays indexed by signal name,
        // and a report has neither an axis nor signals, so every accessor would
        // have to answer with an empty array — which is also what a `SimResult`
        // returns when something went wrong.
        match analysis.to_lowercase().as_str() {
            "tf" => return Ok(self.tf(py, kwargs)?.into_py(py)),
            "sens" => return Ok(self.sens(py, kwargs)?.into_py(py)),
            "pz" => return Ok(self.pz(py, kwargs)?.into_py(py)),
            _ => {}
        }

        // Prologue extracted so `tf`/`sens`/`pz` get the same one.  Four copies
        // of "apply the overrides, then the params, then the alter" is four
        // chances for one entry point to honour something the others do not,
        // which is how `params=` came to be silently dropped once already.
        let (nl, registry, mut opts) = self.prepare(kwargs)?;

        let analysis_lc = analysis.to_lowercase();
        // "dc" with a src kwarg is a sweep; without one it's an op-point alias.
        let is_dc_sweep = analysis_lc == "dc_sweep"
            || (analysis_lc == "dc"
                && kwargs
                    .and_then(|kw| kw.get_item("src").ok().flatten())
                    .is_some());

        if is_dc_sweep {
            let p = parse_dc_kwargs(&nl, kwargs)?;
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
            }
            .into_py(py));
        }

        match analysis_lc.as_str() {
            "op" | "dc" => {
                let result = py
                    .allow_threads(|| dc_op_nr_with_registry_opts(&nl, &registry, &opts))
                    .map_err(sim_err)?;
                Ok(SimResult {
                    inner: SimResultInner::Dc(result),
                    measurements: Vec::new(),
                }
                .into_py(py))
            }
            "tran" | "transient" => {
                let (stop, step) = parse_tran_kwargs(&nl, kwargs, &mut opts)?;
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
                }
                .into_py(py))
            }
            "ac" => {
                let (freqs, src) = parse_ac_kwargs(&nl, kwargs)?;
                let result = py
                    .allow_threads(|| {
                        ac_analysis_opts(&nl, &freqs, src.as_deref(), &registry, &opts)
                    })
                    .map_err(sim_err)?;
                Ok(SimResult {
                    inner: SimResultInner::Ac(result),
                    measurements: Vec::new(),
                }
                .into_py(py))
            }
            "noise" => {
                let (freqs, out_pos, out_neg, input_src) = parse_noise_kwargs(&nl, kwargs)?;
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
                }
                .into_py(py))
            }
            other => Err(PyRuntimeError::new_err(format!(
                "unknown analysis '{}'; use 'op', 'tran', 'ac', 'noise', 'dc_sweep', \
                 or one of the small-signal reports 'tf', 'sens', 'pz' (which are also \
                 ckt.tf(), ckt.sens() and ckt.pz())",
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
    /// integrator, so `variable_step=True` (from a kwarg or from the deck's
    /// `.options`) is an error rather than a silent downgrade — the co-state
    /// recursion replays the forward step sequence, and an adaptive controller
    /// would re-decide that sequence under a perturbed parameter.
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

        let mut nl = netlist.clone();
        apply_overrides(&mut nl, &self.overrides)?;
        apply_source_overrides(&mut nl, &self.source_overrides);
        apply_param_kwarg(&mut nl, kwargs)?;

        let registry = build_registry(&nl, self.netlist_dir.as_ref())?;
        let mut opts = build_sim_options(&nl, kwargs)?;
        let (stop, step) = parse_tran_kwargs(&nl, kwargs, &mut opts)?;

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

    /// Operating-point sensitivities: every probe's value and its gradient
    /// with respect to every parameter, from one solve.
    ///
    /// The DC counterpart of `tran_adjoint` / `ac_adjoint`. It is a single call
    /// rather than a forward/backward pair because an operating point has one
    /// solve to differentiate — there is no trajectory to seed.
    ///
    /// ```python
    /// r = ckt.dc_adjoint(probes={"p": ("power", "bar", 0)},
    ///                    wrt=["Vh.dc"])
    /// r.values["p"]        # W
    /// r.grad["p"]          # dP/dVh, one entry per parameter
    /// ```
    ///
    /// Raises if a parameter reaches nothing, for the same reason the other two
    /// do: a placeholder zero is indistinguishable from a real insensitivity.
    /// `wrt` is the list of `"element.param"` names to differentiate against.
    /// It is deliberately not called `params`: that kwarg is already taken, on
    /// all three of these entry points, by the per-run value overrides an
    /// optimiser's inner loop uses.
    #[pyo3(signature = (probes, wrt, **kwargs))]
    pub fn dc_adjoint(
        &self,
        py: Python<'_>,
        probes: &Bound<'_, PyDict>,
        wrt: Vec<String>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<DcAdjointResult> {
        let netlist = self
            .netlist
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("no netlist loaded; call load() first"))?;
        let mut nl = netlist.clone();
        apply_overrides(&mut nl, &self.overrides)?;
        apply_source_overrides(&mut nl, &self.source_overrides);
        apply_param_kwarg(&mut nl, kwargs)?;
        let registry = build_registry(&nl, self.netlist_dir.as_ref())?;
        let opts = build_sim_options(&nl, kwargs)?;

        let declared: Vec<(String, Output)> = probes
            .iter()
            .map(|(k, v)| Ok((k.extract::<String>()?, parse_probe(&v)?)))
            .collect::<PyResult<_>>()?;
        if declared.is_empty() {
            return Err(PyRuntimeError::new_err(
                "dc_adjoint needs at least one probe; nothing else can be differentiated",
            ));
        }
        let outs: Vec<Output> = declared.iter().map(|(_, o)| o.clone()).collect();
        let refs: Vec<ParamRef> = wrt
            .iter()
            .map(|p| parse_param(p))
            .collect::<PyResult<_>>()?;

        let s = dc_sensitivity(&nl, &registry, &opts, &outs, &refs).map_err(sim_err)?;

        let unreached: Vec<&String> = wrt
            .iter()
            .zip(s.reached.iter())
            .filter(|(_, ok)| !**ok)
            .map(|(p, _)| p)
            .collect();
        if !unreached.is_empty() {
            return Err(PyRuntimeError::new_err(format!(
                "these parameters reach nothing in the equations, so their gradient is a \
                 placeholder zero rather than a computed one: {unreached:?}"
            )));
        }
        let shaky: Vec<(&String, f64)> = wrt
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
                    "fairchild: the gradient for {shaky:?} could not be resolved to better \
                     than the relative error shown; pass a step explicitly if you know one."
                ),),
            )?;
        }

        Ok(DcAdjointResult {
            names: declared.into_iter().map(|(n, _)| n).collect(),
            values: s.values,
            grad: s.grad,
        })
    }

    /// Run an `.ac` sweep that can be differentiated, and return the run.
    ///
    /// The frequency-domain counterpart of `tran_adjoint`, and the one a filter
    /// design wants: "put the resonance here", "flatten this passband" are
    /// least-squares fits of a response against a target.
    ///
    /// ```python
    /// run = ckt.ac_adjoint(node="out0_re_0", fstart=1e8, fstop=1e11,
    ///                      points=40, src="Vd")
    /// y = run.response                              # |V|^2 per frequency
    /// g = run.backward(2 * (y - target),            # dL/dy per frequency
    ///                  ["Xps.l_um", "Vd.dc"])       # -> (2,) numpy
    /// ```
    ///
    /// `node` is the net to read. `quantity` picks what is read off it:
    /// `"mag2"` (default, `|V|^2` — smooth even at a null, which `|V|` is not),
    /// `"re"`, or `"im"`.
    ///
    /// Everything else is the same kwargs as `run("ac", ...)`.
    #[pyo3(signature = (node, quantity=None, **kwargs))]
    pub fn ac_adjoint(
        &self,
        node: &str,
        quantity: Option<&str>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<AcAdjointRun> {
        let netlist = self
            .netlist
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("no netlist loaded; call load() first"))?;

        let mut nl = netlist.clone();
        apply_overrides(&mut nl, &self.overrides)?;
        apply_source_overrides(&mut nl, &self.source_overrides);
        apply_param_kwarg(&mut nl, kwargs)?;

        let registry = build_registry(&nl, self.netlist_dir.as_ref())?;
        let opts = build_sim_options(&nl, kwargs)?;
        let (freqs, src) = parse_ac_kwargs(&nl, kwargs)?;

        let out = match quantity.unwrap_or("mag2") {
            "mag2" | "mag_squared" | "power" => AcOutput::MagSquared { node: node.into() },
            "re" | "real" => AcOutput::Real { node: node.into() },
            "im" | "imag" => AcOutput::Imag { node: node.into() },
            other => {
                return Err(PyRuntimeError::new_err(format!(
                    "unknown quantity '{other}'; use 'mag2', 're' or 'im'"
                )))
            }
        };

        let inner =
            AcAdjoint::run(&nl, &registry, &opts, &freqs, src.as_deref()).map_err(sim_err)?;
        Ok(AcAdjointRun {
            inner,
            registry,
            out,
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
            apply_overrides(&mut nl, &self.overrides)?;
            apply_source_overrides(&mut nl, &self.source_overrides);
            let sweep_override: HashMap<String, f64> =
                [(param.to_lowercase(), val)].into_iter().collect();
            apply_overrides(&mut nl, &sweep_override)?;

            let registry = build_registry(&nl, self.netlist_dir.as_ref())?;
            let mut opts = build_sim_options(&nl, kwargs)?;

            let result = match analysis.to_lowercase().as_str() {
                "op" | "dc" => {
                    let r = dc_op_nr_with_registry_opts(&nl, &registry, &opts).map_err(sim_err)?;
                    SimResult {
                        inner: SimResultInner::Dc(r),
                        measurements: Vec::new(),
                    }
                }
                "tran" | "transient" => {
                    let (stop, step) = parse_tran_kwargs(&nl, kwargs, &mut opts)?;
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
                    let (freqs, src) = parse_ac_kwargs(&nl, kwargs)?;
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

/// Load built-in models, then every model file the netlist named — `.va`
/// sources compiled on the way in, `.osdi` artefacts as-is.
///
/// The resolving and loading is `fairchild_osdi::load_libraries`, shared with
/// the CLI. There is no command line here to carry `--openvaf`, so the compiler
/// and cache come from `FAIRCHILD_OPENVAF` / `FAIRCHILD_VA_CACHE`.
/// Parse twice, so the registry decides WDM dispatch (#52).
///
/// Pass one is permissive and only its `.model` cards and model-file paths are
/// used — neither depends on how bundles expand. The registry built from those
/// then places every instance by what its name really resolves to, which for a
/// card-named device the parser can never work out on its own.
fn two_pass(
    netlist_dir: Option<&PathBuf>,
    parse: impl Fn(&dyn ArityOracle) -> PyResult<Netlist>,
) -> PyResult<Netlist> {
    // Pass one's warnings are pass two's; emitting both would double every one.
    let was_quiet = fairchild_parser::warn::quiet();
    fairchild_parser::warn::set_quiet(true);
    let probe = parse(&PermissiveArity).and_then(|n| {
        let reg = build_registry(&n, netlist_dir)?;
        Ok(reg)
    });
    fairchild_parser::warn::set_quiet(was_quiet);
    match probe {
        Ok(reg) => parse(&reg),
        // Pass one failed: let the honest oracle produce the error the user sees.
        Err(_) => parse(&fairchild_parser::StaticArity),
    }
}

fn build_registry(netlist: &Netlist, netlist_dir: Option<&PathBuf>) -> PyResult<DeviceRegistry> {
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&netlist.models);

    #[cfg(feature = "osdi")]
    {
        fairchild_osdi::load_libraries_with_widths(
            &netlist.osdi_paths,
            &netlist.va_sources,
            netlist_dir.map(|p| p.as_path()),
            &fairchild_osdi::VaOptions::from_env(),
            &mut registry,
            &fairchild_parser::instantiated_widths(netlist),
            fairchild_parser::wires_per_channel(netlist),
        )
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        // `.model <card> <module> (...)` cards naming a descriptor we just
        // loaded. After the libraries, so the descriptors exist to alias.
        registry.register_loaded_model_cards(&netlist.models);
    }

    #[cfg(not(feature = "osdi"))]
    if !netlist.osdi_paths.is_empty() || !netlist.va_sources.is_empty() {
        return Err(PyRuntimeError::new_err(
            "netlist references .osdi/.va model files but this build was compiled without OSDI \
             support; rebuild with --features osdi or remove those references",
        ));
    }

    Ok(registry)
}

// ---------------------------------------------------------------------------
// Helper: parse tran kwargs
// ---------------------------------------------------------------------------

/// True when the caller passed none of `keys`.
///
/// This is the all-or-nothing test behind every card adoption below: a deck
/// card is taken whole or not at all, so one kwarg from the card's own set is
/// enough to hand the whole decision to the caller.
fn none_of(kwargs: Option<&Bound<'_, PyDict>>, keys: &[&str]) -> PyResult<bool> {
    let Some(kw) = kwargs else { return Ok(true) };
    for k in keys {
        if kw.get_item(k)?.is_some() {
            return Ok(false);
        }
    }
    Ok(true)
}

/// The deck's one card of a given kind, or `None` if it declares none.
///
/// Two cards of the same kind is an error rather than a silent pick: which one
/// the caller meant is unknowable, and guessing produces a run nobody asked for.
fn sole_card<'a>(
    netlist: &'a Netlist,
    kind: &str,
    is_kind: impl Fn(&Analysis) -> bool,
) -> PyResult<Option<&'a Analysis>> {
    let mut matching = netlist.analyses.iter().filter(|a| is_kind(a));
    let first = matching.next();
    if first.is_some() && matching.next().is_some() {
        return Err(PyRuntimeError::new_err(format!(
            "deck declares more than one .{kind} card, so run(\"{kind}\") cannot tell \
             which you mean: pass the parameters as kwargs, or read ckt.analyses and \
             pass the one you want"
        )));
    }
    Ok(first)
}

/// Transient timing: the caller's kwargs, or else the deck's `.tran` card taken
/// whole — `step`, `stop`, `tstart`, `tmax` and `UIC` together.
///
/// Never a mix of the two. `tstart` and `tmax` used to reach `opts` from
/// `SimOptions::from_netlist` even when the caller supplied its own `step` and
/// `stop`, which applied half a card to a run that had not asked for any of it.
fn parse_tran_kwargs(
    netlist: &Netlist,
    kwargs: Option<&Bound<'_, PyDict>>,
    opts: &mut SimOptions,
) -> PyResult<(f64, f64)> {
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

    if stop.is_none() && step.is_none() {
        if let Some(Analysis::Tran {
            step,
            stop,
            tstart,
            tmax,
            uic,
        }) = sole_card(netlist, "tran", |a| matches!(a, Analysis::Tran { .. }))?
        {
            opts.apply_tran_card(*tstart, *tmax, *uic);
            return Ok((*stop, *step));
        }
        return Err(PyRuntimeError::new_err(
            "tran needs timing: pass step= and stop= (seconds), or give the deck a \
             .tran card — run(\"tran\") then takes the whole card",
        ));
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
/// The three report shapes, in one place each.
///
/// `ckt.tf()` and a `run_all()` row must hand back the same dict — otherwise
/// code that reads one breaks on the other, and the deck-vs-caller consistency
/// #33 bought would only hold for the waveform analyses.
fn tf_dict<'py>(py: Python<'py>, r: &fairchild_core::TfResult) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new_bound(py);
    d.set_item("gain", r.gain)?;
    d.set_item("r_in", r.r_in)?;
    d.set_item("r_out", r.r_out)?;
    d.set_item("out_value", r.out_value)?;
    Ok(d)
}

fn sens_rows<'py>(
    py: Python<'py>,
    r: &fairchild_core::SensResult,
) -> PyResult<Vec<Bound<'py, PyDict>>> {
    r.rows
        .iter()
        .map(|row| {
            let d = PyDict::new_bound(py);
            d.set_item("param", &row.name)?;
            d.set_item("nominal", row.nominal)?;
            d.set_item("sensitivity", row.sensitivity)?;
            d.set_item("normalised", row.normalised)?;
            d.set_item("reached", row.reached)?;
            d.set_item("fd_error", row.fd_error)?;
            Ok(d)
        })
        .collect()
}

fn pz_dict<'py>(py: Python<'py>, r: &fairchild_core::PzResult) -> PyResult<Bound<'py, PyDict>> {
    // Python complex, so a caller can hand the list straight to numpy and take
    // `abs()` / `angle()` without unpacking pairs first.
    let to_py = |roots: &[fairchild_core::Root]| -> Vec<Bound<'py, PyComplex>> {
        roots
            .iter()
            .map(|k| PyComplex::from_doubles_bound(py, k.re, k.im))
            .collect()
    };
    let d = PyDict::new_bound(py);
    d.set_item("poles", to_py(&r.poles))?;
    d.set_item("zeros", to_py(&r.zeros))?;
    d.set_item("infinite_poles", r.infinite_poles)?;
    d.set_item("infinite_zeros", r.infinite_zeros)?;
    Ok(d)
}

/// Any analysis's result, so `run_all` can carry a heterogeneous deck's output
/// out of the parallel fan-out before any of it becomes a Python object.
///
/// Deliberately *not* `SimResult`: building one needs the GIL, and the whole
/// point of the fan-out is that it holds no GIL. The conversion happens once,
/// afterwards, in [`AnyResult::into_py_object`].
enum AnyResult {
    Dc(fairchild_core::NrResult),
    Tran(fairchild_core::TranResult, Vec<(String, f64)>),
    Ac(fairchild_core::AcResult),
    Noise(fairchild_core::NoiseResult),
    DcSweep(fairchild_core::DcSweepResult),
    Tf(fairchild_core::TfResult),
    Sens(fairchild_core::SensResult),
    Pz(fairchild_core::PzResult),
}

impl AnyResult {
    /// The waveform analyses become a `SimResult`; the reports become the same
    /// dict / list of dicts their own methods return, so a `run_all` row and a
    /// direct `ckt.tf()` are the same shape.
    fn into_py_object(self, py: Python<'_>) -> PyResult<PyObject> {
        let sim = |inner, measurements| {
            SimResult {
                inner,
                measurements,
            }
            .into_py(py)
        };
        Ok(match self {
            AnyResult::Dc(r) => sim(SimResultInner::Dc(r), Vec::new()),
            AnyResult::Tran(r, m) => sim(SimResultInner::Tran(r), m),
            AnyResult::Ac(r) => sim(SimResultInner::Ac(r), Vec::new()),
            AnyResult::Noise(r) => sim(SimResultInner::Noise(r), Vec::new()),
            AnyResult::DcSweep(r) => sim(SimResultInner::DcSweep(r), Vec::new()),
            AnyResult::Tf(r) => tf_dict(py, &r)?.into_py(py),
            AnyResult::Sens(r) => sens_rows(py, &r)?.into_py(py),
            AnyResult::Pz(r) => pz_dict(py, &r)?.into_py(py),
        })
    }
}

/// Run every analysis one corner declares, in deck order.
///
/// Returns `(kind, result)` pairs. Errors come back as strings rather than
/// `PyErr` so the caller can attribute them to a corner — a failure at 125 °C
/// and a failure in the `slow` block are different bugs, and "DC op failed" on
/// its own does not say which.
fn run_one_corner(
    corner: &fairchild_core::Corner,
    netlist_dir: Option<&PathBuf>,
) -> Result<Vec<(String, AnyResult)>, String> {
    let nl = &corner.netlist;
    let opts = &corner.opts;
    let label = || {
        format!(
            "corner alter={} temp={:.1}C",
            corner.alter_label,
            corner.temp_c()
        )
    };
    let registry = build_registry(nl, netlist_dir).map_err(|e| format!("{}: {e}", label()))?;
    let fail = |what: &str, e: SimError| format!("{} at {}: {e}", what, label());

    let mut out = Vec::new();
    for analysis in &nl.analyses {
        let entry = match analysis {
            Analysis::Op => (
                "op".to_string(),
                AnyResult::Dc(
                    dc_op_nr_with_registry_opts(nl, &registry, opts).map_err(|e| fail(".op", e))?,
                ),
            ),
            Analysis::Tran {
                step,
                stop,
                tstart,
                tmax,
                uic,
            } => {
                // Each card's tstart/tmax/UIC belong to that card's run only —
                // the same rule the CLI follows, and the leak #33 closed.
                let mut local = opts.clone();
                local.apply_tran_card(*tstart, *tmax, *uic);
                let r = if local.variable_step {
                    tran_nr_with_registry_var_opts(nl, *step, *stop, &registry, &local)
                } else {
                    tran_nr_with_registry_opts(nl, *step, *stop, &registry, &local)
                }
                .map_err(|e| fail(".tran", e))?;
                let meas = evaluate_measurements(&nl.measurements, &r)
                    .into_iter()
                    .map(|m| (m.name, m.value))
                    .collect();
                ("tran".to_string(), AnyResult::Tran(r, meas))
            }
            Analysis::Ac {
                variation,
                points,
                fstart,
                fstop,
            } => {
                let freqs = freq_points(*variation, *fstart, *fstop, *points);
                (
                    "ac".to_string(),
                    AnyResult::Ac(
                        ac_analysis_opts(nl, &freqs, None, &registry, opts)
                            .map_err(|e| fail(".ac", e))?,
                    ),
                )
            }
            Analysis::Dc {
                src,
                start,
                stop,
                step,
                nested,
            } => {
                let nested = nested
                    .as_ref()
                    .map(|n| (n.src.clone(), n.start, n.stop, n.step));
                (
                    "dc_sweep".to_string(),
                    AnyResult::DcSweep(
                        dc_sweep_with_registry_opts(
                            nl,
                            src,
                            *start,
                            *stop,
                            *step,
                            nested.as_ref().map(|(s, a, b, c)| (s.as_str(), *a, *b, *c)),
                            &registry,
                            opts,
                        )
                        .map_err(|e| fail(".dc", e))?,
                    ),
                )
            }
            Analysis::Noise {
                out_pos,
                out_neg,
                input_src,
                variation,
                points,
                fstart,
                fstop,
            } => {
                let freqs = freq_points(*variation, *fstart, *fstop, *points);
                (
                    "noise".to_string(),
                    AnyResult::Noise(
                        fairchild_core::noise_analysis(
                            nl, &freqs, out_pos, out_neg, input_src, &registry, opts,
                        )
                        .map_err(|e| fail(".noise", e))?,
                    ),
                )
            }
            Analysis::Tf { out, input_src } => (
                "tf".to_string(),
                AnyResult::Tf(
                    fairchild_core::transfer_function(nl, &registry, opts, out, input_src)
                        .map_err(|e| fail(".tf", e))?,
                ),
            ),
            Analysis::Sens { out, params } => (
                "sens".to_string(),
                AnyResult::Sens(
                    fairchild_core::sensitivity(nl, &registry, opts, out, params)
                        .map_err(|e| fail(".sens", e))?,
                ),
            ),
            Analysis::Pz {
                in_pos,
                in_neg,
                out_pos,
                out_neg,
                drive,
                want,
            } => (
                "pz".to_string(),
                AnyResult::Pz(
                    fairchild_core::pole_zero(
                        nl, &registry, opts, in_pos, in_neg, out_pos, out_neg, *drive, *want,
                    )
                    .map_err(|e| fail(".pz", e))?,
                ),
            ),
        };
        out.push(entry);
    }
    Ok(out)
}

/// `DEC`/`OCT`/`LIN` to a frequency list, in one place so `.ac` and `.noise`
/// cannot disagree about what a card's variation means.
fn freq_points(variation: AcVariation, fstart: f64, fstop: f64, points: usize) -> Vec<f64> {
    match variation {
        AcVariation::Dec => fairchild_core::freq_decade(fstart, fstop, points),
        AcVariation::Oct => fairchild_core::freq_oct(fstart, fstop, points),
        AcVariation::Lin => fairchild_core::freq_linear(fstart, fstop, points),
    }
}

/// A `v(a,b)` / `i(v1)` string from a caller, through the parser's own reader.
fn outvar_from_str(s: &str) -> PyResult<OutVar> {
    fairchild_parser::parse_outvar(s, 0)
        .map_err(|e| PyRuntimeError::new_err(format!("bad output '{s}': {e}")))
}

/// An `OutVar` back as the string a caller would have written.
fn outvar_str(out: &OutVar) -> String {
    match out {
        OutVar::NodeVoltage { pos, neg } if neg == "0" => format!("v({pos})"),
        OutVar::NodeVoltage { pos, neg } => format!("v({pos},{neg})"),
        OutVar::BranchCurrent(name) => format!("i({name})"),
    }
}

fn param_str(p: &ParamName) -> String {
    match &p.param {
        Some(x) => format!("{}.{x}", p.element),
        None => p.element.clone(),
    }
}

/// Was any of `keys` passed?  This is what decides between "adopt the deck's
/// card whole" and "the caller supplied its own", and it deliberately ignores
/// the solver options and `alter`/`params`, which are orthogonal to both.
fn any_given(kwargs: Option<&Bound<'_, PyDict>>, keys: &[&str]) -> PyResult<bool> {
    let Some(kw) = kwargs else { return Ok(false) };
    for k in keys {
        if kw.get_item(*k)?.is_some() {
            return Ok(true);
        }
    }
    Ok(false)
}

fn str_kwarg(kwargs: Option<&Bound<'_, PyDict>>, key: &str) -> PyResult<Option<String>> {
    let Some(kw) = kwargs else { return Ok(None) };
    match kw.get_item(key)? {
        Some(v) => Ok(Some(v.extract::<String>()?)),
        None => Ok(None),
    }
}

/// `.tf`'s two fields: the caller's, or else the deck's card taken whole.
fn tf_args(nl: &Netlist, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<(OutVar, String)> {
    if !any_given(kwargs, &["out", "src"])? {
        return match sole_card(nl, "tf", |a| matches!(a, Analysis::Tf { .. }))? {
            Some(Analysis::Tf { out, input_src }) => Ok((out.clone(), input_src.clone())),
            _ => Err(PyRuntimeError::new_err(
                "tf() needs either a .tf card in the deck or out=… and src=… \
                 (e.g. ckt.tf(out=\"v(out)\", src=\"Vin\"))",
            )),
        };
    }
    let out = str_kwarg(kwargs, "out")?
        .ok_or_else(|| PyRuntimeError::new_err("tf(src=…) also needs out=\"v(node)\""))?;
    let src = str_kwarg(kwargs, "src")?
        .ok_or_else(|| PyRuntimeError::new_err("tf(out=…) also needs src=\"Vin\""))?;
    Ok((outvar_from_str(&out)?, src.to_lowercase()))
}

/// `.sens`'s output and parameter list.
fn sens_args(
    nl: &Netlist,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<(OutVar, Vec<ParamName>)> {
    if !any_given(kwargs, &["out", "wrt"])? {
        return match sole_card(nl, "sens", |a| matches!(a, Analysis::Sens { .. }))? {
            Some(Analysis::Sens { out, params }) => Ok((out.clone(), params.clone())),
            _ => Err(PyRuntimeError::new_err(
                "sens() needs either a .sens card in the deck or out=… \
                 (e.g. ckt.sens(out=\"v(out)\"), optionally wrt=[\"r1\", \"m1.w\"])",
            )),
        };
    }
    let out = str_kwarg(kwargs, "out")?
        .ok_or_else(|| PyRuntimeError::new_err("sens(wrt=…) also needs out=\"v(node)\""))?;
    let mut params = Vec::new();
    if let Some(kw) = kwargs {
        if let Some(v) = kw.get_item("wrt")? {
            for name in v.extract::<Vec<String>>()? {
                let lc = name.to_lowercase();
                let (element, param) = match lc.split_once('.') {
                    Some((e, p)) => (e.to_string(), Some(p.to_string())),
                    None => (lc, None),
                };
                params.push(ParamName { element, param });
            }
        }
    }
    Ok((outvar_from_str(&out)?, params))
}

/// `.pz`'s ports and keywords.  Unlike `.tf`, most fields have a defensible
/// default (ground for the far side of each port, a voltage drive, both root
/// sets), so only the two live nodes are required.
#[allow(clippy::type_complexity)]
fn pz_args(
    nl: &Netlist,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<(String, String, String, String, PzDrive, PzWant)> {
    const KEYS: &[&str] = &["in_pos", "in_neg", "out_pos", "out_neg", "drive", "want"];
    if !any_given(kwargs, KEYS)? {
        return match sole_card(nl, "pz", |a| matches!(a, Analysis::Pz { .. }))? {
            Some(Analysis::Pz {
                in_pos,
                in_neg,
                out_pos,
                out_neg,
                drive,
                want,
            }) => Ok((
                in_pos.clone(),
                in_neg.clone(),
                out_pos.clone(),
                out_neg.clone(),
                *drive,
                *want,
            )),
            _ => Err(PyRuntimeError::new_err(
                "pz() needs either a .pz card in the deck or in_pos=… and out_pos=… \
                 (e.g. ckt.pz(in_pos=\"in\", out_pos=\"out\"))",
            )),
        };
    }
    let need = |k: &str| -> PyResult<String> {
        str_kwarg(kwargs, k)?
            .ok_or_else(|| PyRuntimeError::new_err(format!("pz() needs {k}=\"<node>\"")))
    };
    let drive = match str_kwarg(kwargs, "drive")?.as_deref() {
        None | Some("vol") => PzDrive::Vol,
        Some("cur") => PzDrive::Cur,
        Some(o) => {
            return Err(PyRuntimeError::new_err(format!(
                "pz(drive='{o}') — expected 'vol' or 'cur'"
            )))
        }
    };
    let want = match str_kwarg(kwargs, "want")?.as_deref() {
        None | Some("pz") => PzWant::Both,
        Some("pol") => PzWant::Poles,
        Some("zer") => PzWant::Zeros,
        Some(o) => {
            return Err(PyRuntimeError::new_err(format!(
                "pz(want='{o}') — expected 'pol', 'zer' or 'pz'"
            )))
        }
    };
    Ok((
        need("in_pos")?,
        str_kwarg(kwargs, "in_neg")?.unwrap_or_else(|| "0".into()),
        need("out_pos")?,
        str_kwarg(kwargs, "out_neg")?.unwrap_or_else(|| "0".into()),
        drive,
        want,
    ))
}

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
            "out",    // noise / tf / sens
            "params", // tran_adjoint per-run parameter overrides
            "in_pos", // pz
            "in_neg",
            "drive",
            "want",
            // `.sens`'s parameter list.  Not spelled `params`: that name is
            // already taken by the per-run `element.param` override dict, and
            // one kwarg meaning two things depending on which method reads it
            // is a bug waiting for someone to write a fitting loop.
            "wrt",
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

/// DC-sweep parameters: the caller's kwargs, or else the deck's `.dc` card
/// taken whole (both sweeps of a nested card included).
fn parse_dc_kwargs(netlist: &Netlist, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<DcKwargs> {
    const CARD: &[&str] = &[
        "src", "start", "stop", "step", "src2", "start2", "stop2", "step2",
    ];
    if none_of(kwargs, CARD)? {
        if let Some(Analysis::Dc {
            src,
            start,
            stop,
            step,
            nested,
        }) = sole_card(netlist, "dc", |a| matches!(a, Analysis::Dc { .. }))?
        {
            return Ok(DcKwargs {
                src: src.clone(),
                start: *start,
                stop: *stop,
                step: *step,
                nested: nested
                    .as_ref()
                    .map(|n| (n.src.clone(), n.start, n.stop, n.step)),
            });
        }
    }
    let kw = kwargs.ok_or_else(|| {
        PyRuntimeError::new_err(
            "dc sweep needs src, start, stop, step kwargs, or a .dc card in the deck",
        )
    })?;

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

/// AC sweep: the caller's frequency kwargs, or else the deck's `.ac` card taken
/// whole.  `src` is not on the card — it names the excitation source and stays
/// the caller's either way.
fn parse_ac_kwargs(
    netlist: &Netlist,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<(Vec<f64>, Option<String>)> {
    let mut fstart: Option<f64> = None;
    let mut fstop: Option<f64> = None;
    let mut points: usize = 20;
    let mut variation = AcVariation::Dec;
    let mut src: Option<String> = None;

    if none_of(kwargs, &["fstart", "fstop", "points", "variation"])? {
        if let Some(Analysis::Ac {
            variation,
            points,
            fstart,
            fstop,
        }) = sole_card(netlist, "ac", |a| matches!(a, Analysis::Ac { .. }))?
        {
            let src = match kwargs {
                Some(kw) => kw.get_item("src")?.and_then(|v| v.extract().ok()),
                None => None,
            };
            return Ok((ac_freqs(*variation, *fstart, *fstop, *points), src));
        }
    }

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

    let fstart = fstart.ok_or_else(|| {
        PyRuntimeError::new_err("ac requires 'fstart' kwarg (Hz), or a .ac card in the deck")
    })?;
    let fstop = fstop.ok_or_else(|| {
        PyRuntimeError::new_err("ac requires 'fstop' kwarg (Hz), or a .ac card in the deck")
    })?;

    Ok((ac_freqs(variation, fstart, fstop, points), src))
}

fn variation_name(v: AcVariation) -> &'static str {
    match v {
        AcVariation::Dec => "dec",
        AcVariation::Oct => "oct",
        AcVariation::Lin => "lin",
    }
}

/// The frequency vector one `DEC|OCT|LIN points fstart fstop` sweep describes.
/// Shared so a card and a kwarg set can never disagree about what a sweep means.
fn ac_freqs(variation: AcVariation, fstart: f64, fstop: f64, points: usize) -> Vec<f64> {
    match variation {
        AcVariation::Dec => freq_decade(fstart, fstop, points),
        AcVariation::Oct => freq_oct(fstart, fstop, points),
        AcVariation::Lin => freq_linear(fstart, fstop, points),
    }
}

/// Noise sweep: the caller's kwargs, or else the deck's `.noise` card taken
/// whole — observation nodes, input source and frequency sweep together.
fn parse_noise_kwargs(
    netlist: &Netlist,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<(Vec<f64>, String, String, String)> {
    const CARD: &[&str] = &[
        "out",
        "out_pos",
        "out_neg",
        "src",
        "fstart",
        "fstop",
        "points",
        "variation",
    ];
    if none_of(kwargs, CARD)? {
        if let Some(Analysis::Noise {
            out_pos,
            out_neg,
            input_src,
            variation,
            points,
            fstart,
            fstop,
        }) = sole_card(netlist, "noise", |a| matches!(a, Analysis::Noise { .. }))?
        {
            return Ok((
                ac_freqs(*variation, *fstart, *fstop, *points),
                out_pos.to_lowercase(),
                out_neg.to_lowercase(),
                input_src.to_lowercase(),
            ));
        }
    }
    let kw = kwargs.ok_or_else(|| {
        PyRuntimeError::new_err(
            "noise requires kwargs: out (or out_pos+out_neg), src, fstart, fstop \
             — or a .noise card in the deck",
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

    Ok((
        ac_freqs(variation, fstart, fstop, points),
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

/// A parameter spec from Python: `"El.param"`, or `("El.param", step)` to pin
/// the finite-difference step `∂G/∂p` is taken with.
///
/// The step exists because the automatic choice cannot always work. `∂G/∂p` is
/// differenced numerically, and the default `∛ε·|p|` assumes the residual
/// varies on the parameter's own scale — an optical length does not, since it
/// moves the propagation phase by ~17 rad/µm and a power objective is then a
/// near-total cancellation between two much larger terms. The gradient warns
/// when it could not resolve one; without this there was no way to act on that.
fn parse_param_spec(spec: &Bound<'_, PyAny>) -> PyResult<ParamRef> {
    if let Ok(name) = spec.extract::<String>() {
        return parse_param(&name);
    }
    let (name, step): (String, f64) = spec.extract().map_err(|_| {
        PyRuntimeError::new_err(
            "a parameter is 'element.param', or a ('element.param', step) tuple",
        )
    })?;
    Ok(parse_param(&name)?.with_step(step))
}

/// The display name for a spec, for error and warning messages.
fn param_spec_name(spec: &Bound<'_, PyAny>) -> String {
    spec.extract::<String>()
        .or_else(|_| spec.extract::<(String, f64)>().map(|(n, _)| n))
        .unwrap_or_else(|_| "<unparsable>".to_string())
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

#[pyclass]
pub struct DcAdjointResult {
    names: Vec<String>,
    values: Vec<f64>,
    grad: Vec<Vec<f64>>,
}

#[pymethods]
impl DcAdjointResult {
    /// `{probe: value}` at the operating point.
    #[getter]
    fn values<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new_bound(py);
        for (n, v) in self.names.iter().zip(self.values.iter()) {
            d.set_item(n, *v)?;
        }
        Ok(d)
    }

    /// `{probe: array}` — `dprobe/dp`, one entry per parameter, in the order
    /// the parameters were given.
    #[getter]
    fn grad<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new_bound(py);
        for (n, g) in self.names.iter().zip(self.grad.iter()) {
            d.set_item(n, PyArray1::from_slice_bound(py, g))?;
        }
        Ok(d)
    }
}

#[pyclass(unsendable)]
pub struct AcAdjointRun {
    inner: AcAdjoint,
    registry: DeviceRegistry,
    out: AcOutput,
}

#[pymethods]
impl AcAdjointRun {
    /// The swept frequencies, in Hz.
    #[getter]
    fn freqs<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice_bound(py, self.inner.freqs())
    }

    /// The chosen quantity at each frequency.
    #[getter]
    fn response<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let r = self.inner.response(&self.out).map_err(sim_err)?;
        Ok(PyArray1::from_vec_bound(py, r))
    }

    /// `dL/dp` for each named parameter, given `dL/d(response)` per frequency.
    ///
    /// `cotangent` is an array of `len(freqs)` — the derivative of your loss
    /// with respect to the response at each frequency. Build the loss in numpy
    /// however you like; this only needs its derivative. For a least-squares
    /// fit against a target that is `2*(response - target)`.
    ///
    /// Raises if a parameter reaches nothing in the equations: a silent zero is
    /// indistinguishable from a real insensitivity and would stall an optimiser
    /// somewhere that merely looks stationary. Warns rather than silently
    /// returning a number when the finite difference behind `∂A/∂p` could not
    /// be resolved well.
    fn backward<'py>(
        &self,
        py: Python<'py>,
        cotangent: PyReadonlyArray1<f64>,
        params: Vec<Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let weights: Vec<f64> = cotangent.as_array().to_vec();
        let n_f = self.inner.freqs().len();
        if weights.len() != n_f {
            return Err(PyRuntimeError::new_err(format!(
                "cotangent has length {} but the sweep has {n_f} frequencies",
                weights.len()
            )));
        }
        let (_, seeds) = self.inner.weighted(&self.out, &weights).map_err(sim_err)?;
        let refs: Vec<ParamRef> = params
            .iter()
            .map(parse_param_spec)
            .collect::<PyResult<_>>()?;
        let names: Vec<String> = params.iter().map(param_spec_name).collect();
        let s = self
            .inner
            .gradient(&self.registry, &seeds, &refs)
            .map_err(sim_err)?;

        let unreached: Vec<&String> = names
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

        let shaky: Vec<(&String, f64)> = names
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
    m.add_class::<AcAdjointRun>()?;
    m.add_class::<DcAdjointResult>()?;
    Ok(())
}
