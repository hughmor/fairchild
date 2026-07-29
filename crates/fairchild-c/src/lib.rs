//! C ABI for fairchild — `libfairchild`, in the spirit of `libngspice`.
//!
//! Two layers, both over the same core:
//!
//! * **Batch** — load a netlist, run `.op` / `.tran`, read result vectors.
//!   `fc_sim_new` … `fc_run_tran` … `fc_signal`.
//! * **Stepping** — hold a transient open between timesteps so the host program
//!   drives the clock: read node voltages, write source values, step.
//!   `fc_stepper_new` … `fc_get_node` / `fc_set_source` / `fc_step`.
//!
//! The stepping layer is what mixed-signal co-simulation needs, and unlike
//! libngspice there is no process-global state: every handle is independent, so
//! a host can run as many simulations concurrently as it has threads.
//!
//! Conventions:
//!
//! * Every function that can fail returns `int` — `FC_OK` (0) or an `FC_ERR_*`
//!   code.  `fc_error` / `fc_stepper_error` return the matching message,
//!   owned by the handle and valid until the next call on that handle.
//! * Strings in are NUL-terminated UTF-8.  A NULL handle or NULL string
//!   argument is `FC_ERR_ARG`, never a crash.
//! * Rust panics are caught at the boundary and reported as `FC_ERR_PANIC`
//!   rather than unwinding into C (which would be undefined behaviour).
//! * A handle is not thread-safe.  Use one per thread; no locking is done.
//!
//! # Safety
//!
//! Every `extern "C"` function here shares one contract, so it is stated once
//! rather than repeated on each of them:
//!
//! * Handle arguments are either NULL (reported as `FC_ERR_ARG`) or a pointer
//!   returned by `fc_sim_new` / `fc_stepper_new` and not yet freed.
//! * `const char *` arguments are either NULL or NUL-terminated UTF-8.
//! * Out-pointers are either NULL or writable for their type.
//! * Array arguments are valid for the stated element count.
//! * Pointers handed back by `fc_error`, `fc_signal`, and `fc_signal_name`
//!   borrow the handle's storage and die with the next call on that handle.
//!
//! Violating any of those is undefined behaviour in the usual C way; nothing
//! else about calling these functions is.

// The contract above covers every function in this module; per-function
// `# Safety` sections would restate it 20 times.
#![allow(clippy::missing_safety_doc)]

use std::ffi::{c_char, c_double, c_int, CStr, CString};
use std::panic::AssertUnwindSafe;
use std::path::PathBuf;

use fairchild_core::{
    dc_op_nr_with_registry_opts, tran_nr_with_registry_opts, tran_nr_with_registry_var_opts,
    DeviceRegistry, NrResult, SimError, SimOptions, TranResult, TranStepper,
};
use fairchild_parser::{parse_spice, parse_spice_file, Netlist};

// ---------------------------------------------------------------------------
// Status codes — keep in sync with include/fairchild.h
// ---------------------------------------------------------------------------

pub const FC_OK: c_int = 0;
/// Bad argument: NULL pointer, non-UTF-8 string, no netlist loaded yet.
pub const FC_ERR_ARG: c_int = 1;
/// The netlist could not be read or parsed.
pub const FC_ERR_PARSE: c_int = 2;
/// The simulation failed (no convergence, singular matrix, floating node …).
pub const FC_ERR_SIM: c_int = 3;
/// The requested node, source, or signal does not exist.
pub const FC_ERR_NOT_FOUND: c_int = 4;
/// A panic was caught at the boundary.  The handle must be freed, not reused.
pub const FC_ERR_PANIC: c_int = 5;

struct ApiError {
    code: c_int,
    msg: String,
}

impl ApiError {
    fn arg(msg: impl Into<String>) -> Self {
        ApiError {
            code: FC_ERR_ARG,
            msg: msg.into(),
        }
    }
    fn not_found(msg: impl Into<String>) -> Self {
        ApiError {
            code: FC_ERR_NOT_FOUND,
            msg: msg.into(),
        }
    }
}

impl From<SimError> for ApiError {
    fn from(e: SimError) -> Self {
        // UnknownNode is a lookup miss, not a solver failure — the host asked
        // for something that isn't in the circuit, and NOT_FOUND says so.
        let code = match e {
            SimError::UnknownNode(_) | SimError::UnknownModel(_) => FC_ERR_NOT_FOUND,
            SimError::Parse(_) => FC_ERR_PARSE,
            _ => FC_ERR_SIM,
        };
        ApiError {
            code,
            msg: e.to_string(),
        }
    }
}

type ApiResult = Result<(), ApiError>;

// ---------------------------------------------------------------------------
// Boundary plumbing
// ---------------------------------------------------------------------------

/// Anything with an error slot the boundary can write into.
trait ErrSlot {
    fn store(&mut self, msg: Option<CString>);
}

/// Run `f` against the handle, catching panics and recording any error message.
fn entry<H: ErrSlot, F: FnOnce(&mut H) -> ApiResult>(handle: *mut H, f: F) -> c_int {
    if handle.is_null() {
        return FC_ERR_ARG;
    }
    let res = std::panic::catch_unwind(AssertUnwindSafe(|| f(unsafe { &mut *handle })));
    // Re-borrow after the closure so the error branch can write to the slot.
    let h = unsafe { &mut *handle };
    match res {
        Ok(Ok(())) => {
            h.store(None);
            FC_OK
        }
        Ok(Err(e)) => {
            h.store(CString::new(e.msg).ok());
            e.code
        }
        Err(_) => {
            h.store(
                CString::new("panic inside fairchild — free this handle, do not reuse it").ok(),
            );
            FC_ERR_PANIC
        }
    }
}

/// Borrow a `*const c_char` as `&str`.
unsafe fn cstr<'a>(p: *const c_char) -> Result<&'a str, ApiError> {
    if p.is_null() {
        return Err(ApiError::arg("NULL string argument"));
    }
    CStr::from_ptr(p)
        .to_str()
        .map_err(|_| ApiError::arg("string argument is not valid UTF-8"))
}

fn out<'a, T>(p: *mut T) -> Result<&'a mut T, ApiError> {
    unsafe { p.as_mut() }.ok_or_else(|| ApiError::arg("NULL output pointer"))
}

fn err_ptr(slot: &Option<CString>) -> *const c_char {
    match slot {
        Some(s) => s.as_ptr(),
        None => std::ptr::null(),
    }
}

// ---------------------------------------------------------------------------
// Registry construction (mirrors the Python binding, including OSDI)
// ---------------------------------------------------------------------------

fn build_registry(netlist: &Netlist, dir: Option<&PathBuf>) -> Result<DeviceRegistry, ApiError> {
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&netlist.models);

    #[cfg(feature = "osdi")]
    for osdi_path in &netlist.osdi_paths {
        let path = if std::path::Path::new(osdi_path).is_absolute() {
            PathBuf::from(osdi_path)
        } else if let Some(d) = dir {
            d.join(osdi_path)
        } else {
            PathBuf::from(osdi_path)
        };
        let lib = unsafe { fairchild_osdi::OsdiLibrary::open(&path) }.map_err(|e| ApiError {
            code: FC_ERR_PARSE,
            msg: format!("failed to load OSDI library '{}': {e}", path.display()),
        })?;
        std::sync::Arc::new(lib).register_into(&mut registry);
    }

    #[cfg(not(feature = "osdi"))]
    {
        let _ = dir;
        if !netlist.osdi_paths.is_empty() {
            return Err(ApiError {
                code: FC_ERR_PARSE,
                msg: "netlist references .osdi files but this build has no OSDI support; \
                      rebuild with --features osdi"
                    .into(),
            });
        }
    }

    Ok(registry)
}

// ---------------------------------------------------------------------------
// fc_sim — netlist + options + batch results
// ---------------------------------------------------------------------------

/// Opaque simulation handle.  Created by `fc_sim_new`, freed by `fc_sim_free`.
pub struct FcSim {
    netlist: Option<Netlist>,
    /// Directory of a file-loaded netlist, for resolving relative `.osdi` paths.
    dir: Option<PathBuf>,
    /// `.options`-style overrides applied on top of the netlist's own.
    opt_overrides: Vec<(String, String)>,
    tran: Option<TranResult>,
    op: Option<NrResult>,
    /// Signal names for `fc_signal_name`, rebuilt after each run.
    signals: Vec<CString>,
    err: Option<CString>,
}

impl ErrSlot for FcSim {
    fn store(&mut self, msg: Option<CString>) {
        self.err = msg;
    }
}

impl FcSim {
    fn netlist(&self) -> Result<&Netlist, ApiError> {
        self.netlist
            .as_ref()
            .ok_or_else(|| ApiError::arg("no netlist loaded — call fc_load_file or fc_load_string"))
    }

    fn netlist_mut(&mut self) -> Result<&mut Netlist, ApiError> {
        self.netlist
            .as_mut()
            .ok_or_else(|| ApiError::arg("no netlist loaded — call fc_load_file or fc_load_string"))
    }

    /// Netlist `.options`, then the overrides set through `fc_set_option`.
    fn options(&self) -> Result<SimOptions, ApiError> {
        let netlist = self.netlist()?;
        let mut opts = SimOptions::from_netlist(netlist);
        for (k, v) in &self.opt_overrides {
            if !opts.set(k, v) {
                return Err(ApiError::arg(format!("unrecognised option '{k}'")));
            }
        }
        Ok(opts)
    }

    fn loaded(&mut self, netlist: Netlist, dir: Option<PathBuf>) {
        self.netlist = Some(netlist);
        self.dir = dir;
        self.tran = None;
        self.op = None;
        self.signals.clear();
    }
}

/// Library version string (static, never NULL).
#[no_mangle]
pub extern "C" fn fc_version() -> *const c_char {
    concat!(env!("CARGO_PKG_VERSION"), "\0").as_ptr() as *const c_char
}

/// Create an empty simulation handle.  Returns NULL only if allocation fails.
#[no_mangle]
pub extern "C" fn fc_sim_new() -> *mut FcSim {
    Box::into_raw(Box::new(FcSim {
        netlist: None,
        dir: None,
        opt_overrides: Vec::new(),
        tran: None,
        op: None,
        signals: Vec::new(),
        err: None,
    }))
}

/// Free a handle.  NULL is a no-op.  Invalidates every pointer previously
/// returned by `fc_signal` / `fc_signal_name` / `fc_error` for this handle.
#[no_mangle]
pub unsafe extern "C" fn fc_sim_free(sim: *mut FcSim) {
    if !sim.is_null() {
        drop(Box::from_raw(sim));
    }
}

/// Last error message for this handle, or NULL if the last call succeeded.
/// Valid until the next call on the same handle.
#[no_mangle]
pub unsafe extern "C" fn fc_error(sim: *const FcSim) -> *const c_char {
    match sim.as_ref() {
        Some(s) => err_ptr(&s.err),
        None => std::ptr::null(),
    }
}

/// Parse a netlist file.  Clears any previous results.
#[no_mangle]
pub unsafe extern "C" fn fc_load_file(sim: *mut FcSim, path: *const c_char) -> c_int {
    entry(sim, |s| {
        let path = std::path::Path::new(cstr(path)?);
        let netlist = parse_spice_file(path).map_err(|e| ApiError {
            code: FC_ERR_PARSE,
            msg: e.to_string(),
        })?;
        let dir = path.parent().map(|p| p.to_path_buf());
        s.loaded(netlist, dir);
        Ok(())
    })
}

/// Parse a netlist from a string.  Clears any previous results.
#[no_mangle]
pub unsafe extern "C" fn fc_load_string(sim: *mut FcSim, text: *const c_char) -> c_int {
    entry(sim, |s| {
        let netlist = parse_spice(cstr(text)?).map_err(|e| ApiError {
            code: FC_ERR_PARSE,
            msg: e.to_string(),
        })?;
        s.loaded(netlist, None);
        Ok(())
    })
}

/// Set a solver option by name, as `.options KEY=VALUE` would
/// (`reltol`, `method`, `solver`, `variable_step`, `maxstep`, …).
#[no_mangle]
pub unsafe extern "C" fn fc_set_option(
    sim: *mut FcSim,
    key: *const c_char,
    value: *const c_char,
) -> c_int {
    entry(sim, |s| {
        let (k, v) = (cstr(key)?.to_string(), cstr(value)?.to_string());
        // Validate immediately against a throwaway SimOptions so a typo is
        // reported here rather than at the next run.
        if !SimOptions::default().set(&k, &v) {
            return Err(ApiError::arg(format!("unrecognised option '{k}'")));
        }
        s.opt_overrides.retain(|(existing, _)| *existing != k);
        s.opt_overrides.push((k, v));
        Ok(())
    })
}

/// Retarget an element parameter: `fc_set_param(sim, "R1", "value", 2e3)`.
/// Passives take `value` or their physical name; sources take `dc`/`value`;
/// MOSFET and OSDI instances take any instance parameter.
#[no_mangle]
pub unsafe extern "C" fn fc_set_param(
    sim: *mut FcSim,
    element: *const c_char,
    param: *const c_char,
    value: c_double,
) -> c_int {
    entry(sim, |s| {
        let (element, param) = (cstr(element)?, cstr(param)?);
        let netlist = s.netlist_mut()?;
        if !fairchild_core::set_element_param(netlist, element, param, value) {
            return Err(ApiError::not_found(format!(
                "no element '{element}' with a settable parameter '{param}'"
            )));
        }
        Ok(())
    })
}

/// Replace a source's waveform with a piecewise-linear table of `n` points.
/// `t` must be ascending.  The arrays are copied; the caller keeps ownership.
///
/// This is the offline path: use it when the stimulus is known up front.  For a
/// value decided during the run, use `fc_set_source` on a stepper.
#[no_mangle]
pub unsafe extern "C" fn fc_set_source_pwl(
    sim: *mut FcSim,
    name: *const c_char,
    t: *const c_double,
    v: *const c_double,
    n: usize,
) -> c_int {
    entry(sim, |s| {
        let name = cstr(name)?;
        if n == 0 {
            return Err(ApiError::arg("PWL table needs at least one point"));
        }
        if t.is_null() || v.is_null() {
            return Err(ApiError::arg("NULL PWL array"));
        }
        let ts = std::slice::from_raw_parts(t, n);
        let vs = std::slice::from_raw_parts(v, n);
        if ts.windows(2).any(|w| w[1] < w[0]) {
            return Err(ApiError::arg("PWL times must be ascending"));
        }
        let points: Vec<(f64, f64)> = ts.iter().copied().zip(vs.iter().copied()).collect();
        let netlist = s.netlist_mut()?;
        if !fairchild_core::set_source_pwl(netlist, name, points) {
            return Err(ApiError::not_found(format!("no source named '{name}'")));
        }
        Ok(())
    })
}

// ---------------------------------------------------------------------------
// Batch analyses
// ---------------------------------------------------------------------------

/// Run the DC operating point.  Read values with `fc_op_node` / `fc_op_current`.
#[no_mangle]
pub unsafe extern "C" fn fc_run_op(sim: *mut FcSim) -> c_int {
    entry(sim, |s| {
        let opts = s.options()?;
        let netlist = s.netlist()?;
        let registry = build_registry(netlist, s.dir.as_ref())?;
        let r = dc_op_nr_with_registry_opts(netlist, &registry, &opts)?;
        s.op = Some(r);
        s.tran = None;
        s.signals.clear();
        Ok(())
    })
}

/// Run a transient from 0 to `stop` with timestep `step`.  Read waveforms with
/// `fc_signal`.  Honours `variable_step` if set via `fc_set_option`.
#[no_mangle]
pub unsafe extern "C" fn fc_run_tran(sim: *mut FcSim, step: c_double, stop: c_double) -> c_int {
    entry(sim, |s| {
        // Rejects NaN and infinity too, which `<= 0.0` would let through.
        if !step.is_finite() || step <= 0.0 || !stop.is_finite() || stop <= 0.0 {
            return Err(ApiError::arg("step and stop must be finite and positive"));
        }
        let opts = s.options()?;
        let netlist = s.netlist()?;
        let registry = build_registry(netlist, s.dir.as_ref())?;
        let r = if opts.variable_step {
            tran_nr_with_registry_var_opts(netlist, step, stop, &registry, &opts)?
        } else {
            tran_nr_with_registry_opts(netlist, step, stop, &registry, &opts)?
        };
        s.signals = std::iter::once("time".to_string())
            .chain(r.node_voltages.keys().map(|n| format!("V({n})")))
            .chain(r.vsrc_currents.keys().map(|n| format!("I({n})")))
            .filter_map(|n| CString::new(n).ok())
            .collect();
        s.tran = Some(r);
        s.op = None;
        Ok(())
    })
}

/// Number of signals available from the last transient run (0 otherwise).
#[no_mangle]
pub unsafe extern "C" fn fc_signal_count(sim: *const FcSim) -> usize {
    sim.as_ref().map_or(0, |s| s.signals.len())
}

/// Name of signal `i`, or NULL if out of range.  Valid until the next run.
#[no_mangle]
pub unsafe extern "C" fn fc_signal_name(sim: *const FcSim, i: usize) -> *const c_char {
    sim.as_ref()
        .and_then(|s| s.signals.get(i))
        .map_or(std::ptr::null(), |c| c.as_ptr())
}

/// Borrow a transient result vector by name: `"time"`, `"V(out)"`, `"I(V1)"`.
///
/// `*data` points into the handle's own storage — no copy, and no free.  It is
/// invalidated by the next run on this handle or by `fc_sim_free`.
#[no_mangle]
pub unsafe extern "C" fn fc_signal(
    sim: *mut FcSim,
    name: *const c_char,
    data: *mut *const c_double,
    len: *mut usize,
) -> c_int {
    entry(sim, |s| {
        let name = cstr(name)?;
        let data = out(data)?;
        let len = out(len)?;
        let r = s
            .tran
            .as_ref()
            .ok_or_else(|| ApiError::arg("no transient result — call fc_run_tran first"))?;

        let series: &[f64] = if name.eq_ignore_ascii_case("time") {
            &r.time
        } else if let Some(node) = strip_call(name, "v") {
            r.node_voltages
                .get(node)
                .ok_or_else(|| ApiError::not_found(format!("unknown node '{node}'")))?
        } else if let Some(vsrc) = strip_call(name, "i") {
            r.vsrc_currents
                .get(vsrc)
                .ok_or_else(|| ApiError::not_found(format!("unknown voltage source '{vsrc}'")))?
        } else {
            return Err(ApiError::arg(format!(
                "unrecognised signal '{name}'; use \"time\", \"V(node)\", or \"I(vsrc)\""
            )));
        };
        *data = series.as_ptr();
        *len = series.len();
        Ok(())
    })
}

/// Operating-point node voltage, after `fc_run_op`.
#[no_mangle]
pub unsafe extern "C" fn fc_op_node(
    sim: *mut FcSim,
    node: *const c_char,
    value: *mut c_double,
) -> c_int {
    entry(sim, |s| {
        let node = cstr(node)?;
        let value = out(value)?;
        let r =
            s.op.as_ref()
                .ok_or_else(|| ApiError::arg("no operating point — call fc_run_op first"))?;
        *value = r.node_voltage(node)?;
        Ok(())
    })
}

/// Operating-point current through a voltage source, after `fc_run_op`.
#[no_mangle]
pub unsafe extern "C" fn fc_op_current(
    sim: *mut FcSim,
    vsrc: *const c_char,
    value: *mut c_double,
) -> c_int {
    entry(sim, |s| {
        let vsrc = cstr(vsrc)?;
        let value = out(value)?;
        let r =
            s.op.as_ref()
                .ok_or_else(|| ApiError::arg("no operating point — call fc_run_op first"))?;
        *value = r.vsrc_current(vsrc)?;
        Ok(())
    })
}

/// `V(x)` / `I(x)` → `x`, case-insensitive on the prefix.
fn strip_call<'a>(s: &'a str, f: &str) -> Option<&'a str> {
    let rest = s.strip_suffix(')')?;
    let (head, inner) = rest.split_once('(')?;
    head.trim().eq_ignore_ascii_case(f).then_some(inner.trim())
}

// ---------------------------------------------------------------------------
// fc_stepper — host-driven transient
// ---------------------------------------------------------------------------

/// Opaque handle to a transient paused between timesteps.
pub struct FcStepper {
    inner: TranStepper,
    err: Option<CString>,
}

impl ErrSlot for FcStepper {
    fn store(&mut self, msg: Option<CString>) {
        self.err = msg;
    }
}

/// Open a transient on `sim` with fixed timestep `step`, solving the operating
/// point (or applying `.ic` under UIC) so the handle starts at t = 0.
///
/// The stepper snapshots the netlist: later `fc_set_param` calls on `sim` do not
/// affect it, and it stays valid after `fc_sim_free`.  Returns NULL on failure —
/// call `fc_error(sim)` for the reason.
#[no_mangle]
pub unsafe extern "C" fn fc_stepper_new(sim: *mut FcSim, step: c_double) -> *mut FcStepper {
    let mut created: *mut FcStepper = std::ptr::null_mut();
    entry(sim, |s| {
        if !step.is_finite() || step <= 0.0 {
            return Err(ApiError::arg("step must be finite and positive"));
        }
        let opts = s.options()?;
        let netlist = s.netlist()?;
        let registry = build_registry(netlist, s.dir.as_ref())?;
        let inner = TranStepper::new(netlist.clone(), &registry, &opts, step)?;
        created = Box::into_raw(Box::new(FcStepper { inner, err: None }));
        Ok(())
    });
    created
}

/// Free a stepper.  NULL is a no-op.
#[no_mangle]
pub unsafe extern "C" fn fc_stepper_free(st: *mut FcStepper) {
    if !st.is_null() {
        drop(Box::from_raw(st));
    }
}

/// Last error message for this stepper, or NULL if the last call succeeded.
#[no_mangle]
pub unsafe extern "C" fn fc_stepper_error(st: *const FcStepper) -> *const c_char {
    match st.as_ref() {
        Some(s) => err_ptr(&s.err),
        None => std::ptr::null(),
    }
}

/// Advance exactly one timestep.  `t_out` may be NULL.
///
/// On `FC_ERR_SIM` (no convergence) the stepper still holds the last accepted
/// timepoint, so the host can back off a drive level and retry.
#[no_mangle]
pub unsafe extern "C" fn fc_step(st: *mut FcStepper, t_out: *mut c_double) -> c_int {
    entry(st, |s| {
        let t = s.inner.step()?;
        if let Some(o) = t_out.as_mut() {
            *o = t;
        }
        Ok(())
    })
}

/// Step until the simulation time reaches `t_target`.  Because the step size is
/// fixed, this lands on the first grid point at or past `t_target`; `t_out`
/// (may be NULL) reports where.  Already past `t_target` ⇒ no steps taken.
#[no_mangle]
pub unsafe extern "C" fn fc_advance_to(
    st: *mut FcStepper,
    t_target: c_double,
    t_out: *mut c_double,
) -> c_int {
    entry(st, |s| {
        let t = s.inner.advance_to(t_target)?;
        if let Some(o) = t_out.as_mut() {
            *o = t;
        }
        Ok(())
    })
}

/// Current simulation time in seconds, or a negative value if `st` is NULL.
#[no_mangle]
pub unsafe extern "C" fn fc_time(st: *const FcStepper) -> c_double {
    st.as_ref().map_or(-1.0, |s| s.inner.time())
}

/// The fixed timestep, after clamping to `maxstep`.  Negative if `st` is NULL.
#[no_mangle]
pub unsafe extern "C" fn fc_step_size(st: *const FcStepper) -> c_double {
    st.as_ref().map_or(-1.0, |s| s.inner.step_size())
}

/// Node voltage at the current timepoint — the analog → digital direction.
#[no_mangle]
pub unsafe extern "C" fn fc_get_node(
    st: *mut FcStepper,
    node: *const c_char,
    value: *mut c_double,
) -> c_int {
    entry(st, |s| {
        let node = cstr(node)?;
        *out(value)? = s.inner.node(node)?;
        Ok(())
    })
}

/// Current through a voltage source at the current timepoint.
#[no_mangle]
pub unsafe extern "C" fn fc_get_current(
    st: *mut FcStepper,
    vsrc: *const c_char,
    value: *mut c_double,
) -> c_int {
    entry(st, |s| {
        let vsrc = cstr(vsrc)?;
        *out(value)? = s.inner.vsrc_current(vsrc)?;
        Ok(())
    })
}

/// Hold a voltage or current source at `value` from the next step on — the
/// digital → analog direction.  Zero-order hold, like ngspice's `GetVSRCData`.
#[no_mangle]
pub unsafe extern "C" fn fc_set_source(
    st: *mut FcStepper,
    name: *const c_char,
    value: c_double,
) -> c_int {
    entry(st, |s| {
        let name = cstr(name)?;
        s.inner.set_source(name, value)?;
        Ok(())
    })
}

/// Number of solvable nodes, for enumerating with `fc_node_name`.
#[no_mangle]
pub unsafe extern "C" fn fc_node_count(st: *const FcStepper) -> usize {
    st.as_ref().map_or(0, |s| s.inner.node_names().count())
}

/// Copy the name of node `i` into `buf` (NUL-terminated, truncated to `cap`).
/// Returns the number of bytes needed excluding the NUL, so a caller can size a
/// buffer; `FC_ERR_NOT_FOUND` maps to a negative return.
#[no_mangle]
pub unsafe extern "C" fn fc_node_name(
    st: *const FcStepper,
    i: usize,
    buf: *mut c_char,
    cap: usize,
) -> isize {
    let Some(s) = st.as_ref() else { return -1 };
    let Some(name) = s.inner.node_names().nth(i) else {
        return -1;
    };
    let bytes = name.as_bytes();
    if !buf.is_null() && cap > 0 {
        let n = bytes.len().min(cap - 1);
        std::ptr::copy_nonoverlapping(bytes.as_ptr() as *const c_char, buf, n);
        *buf.add(n) = 0;
    }
    bytes.len() as isize
}

#[cfg(test)]
mod tests {
    use super::*;

    const RC: &str = "* rc\nV1 in 0 DC 0\nR1 in out 1k\nC1 out 0 1p\n.tran 1n 100n\n.end\n";

    /// Drive every call through the real C entry points, as a C caller would.
    unsafe fn load(text: &str) -> *mut FcSim {
        let sim = fc_sim_new();
        let c = CString::new(text).unwrap();
        assert_eq!(fc_load_string(sim, c.as_ptr()), FC_OK);
        sim
    }

    #[test]
    fn batch_roundtrip() {
        unsafe {
            let sim = load(RC);
            let vin = CString::new("V1").unwrap();
            assert_eq!(
                fc_set_param(sim, vin.as_ptr(), CString::new("dc").unwrap().as_ptr(), 1.0),
                FC_OK
            );
            assert_eq!(fc_run_tran(sim, 1e-9, 100e-9), FC_OK);

            let mut data: *const c_double = std::ptr::null();
            let mut len: usize = 0;
            let key = CString::new("V(out)").unwrap();
            assert_eq!(fc_signal(sim, key.as_ptr(), &mut data, &mut len), FC_OK);
            assert!(len > 50, "got {len} timepoints");
            let series = std::slice::from_raw_parts(data, len);
            // 1 V through 1k into 1p: τ = 1 ns, so by 100 ns it is settled.
            assert!(
                (series[len - 1] - 1.0).abs() < 1e-3,
                "final {}",
                series[len - 1]
            );

            assert!(fc_signal_count(sim) >= 3);
            assert!(!fc_signal_name(sim, 0).is_null());
            assert!(fc_signal_name(sim, 9999).is_null());
            fc_sim_free(sim);
        }
    }

    #[test]
    fn stepping_reads_and_drives() {
        unsafe {
            let sim = load(RC);
            let st = fc_stepper_new(sim, 1e-10);
            assert!(!st.is_null(), "stepper_new failed");

            let out_node = CString::new("out").unwrap();
            let v1 = CString::new("V1").unwrap();
            let mut v = f64::NAN;
            assert_eq!(fc_get_node(st, out_node.as_ptr(), &mut v), FC_OK);
            assert!(v.abs() < 1e-12, "starts discharged, got {v}");

            // Bang-bang around 0.5 V — the mixed-signal loop in miniature.
            assert_eq!(fc_set_source(st, v1.as_ptr(), 1.0), FC_OK);
            for _ in 0..2000 {
                assert_eq!(fc_step(st, std::ptr::null_mut()), FC_OK);
                assert_eq!(fc_get_node(st, out_node.as_ptr(), &mut v), FC_OK);
                assert_eq!(
                    fc_set_source(st, v1.as_ptr(), if v < 0.5 { 1.0 } else { 0.0 }),
                    FC_OK
                );
            }
            assert!(
                (v - 0.5).abs() < 0.05,
                "should hold near threshold, got {v}"
            );

            let mut t = f64::NAN;
            assert_eq!(fc_advance_to(st, fc_time(st) + 5e-10, &mut t), FC_OK);
            assert!(t > fc_time(st) - 1e-18);

            fc_stepper_free(st);
            fc_sim_free(sim);
        }
    }

    /// The stepper snapshots the netlist, so it must outlive the sim handle.
    #[test]
    fn stepper_outlives_its_sim() {
        unsafe {
            let sim = load(RC);
            let st = fc_stepper_new(sim, 1e-10);
            fc_sim_free(sim);
            assert_eq!(fc_step(st, std::ptr::null_mut()), FC_OK);
            fc_stepper_free(st);
        }
    }

    #[test]
    fn bad_input_is_an_error_not_a_crash() {
        unsafe {
            // NULL handles.
            assert_eq!(fc_run_op(std::ptr::null_mut()), FC_ERR_ARG);
            assert_eq!(
                fc_step(std::ptr::null_mut(), std::ptr::null_mut()),
                FC_ERR_ARG
            );
            assert!(fc_error(std::ptr::null()).is_null());
            assert_eq!(fc_signal_count(std::ptr::null()), 0);
            fc_sim_free(std::ptr::null_mut());
            fc_stepper_free(std::ptr::null_mut());

            let sim = fc_sim_new();
            // No netlist yet.
            assert_eq!(fc_run_op(sim), FC_ERR_ARG);
            assert!(!fc_error(sim).is_null());
            // NULL string.
            assert_eq!(fc_load_string(sim, std::ptr::null()), FC_ERR_ARG);
            // Unparseable.
            let junk = CString::new("R1 dangling\n").unwrap();
            assert_ne!(fc_load_string(sim, junk.as_ptr()), FC_OK);

            assert_eq!(
                fc_load_string(sim, CString::new(RC).unwrap().as_ptr()),
                FC_OK
            );
            // Unknown option, element, signal.
            let (k, v) = (
                CString::new("nosuchopt").unwrap(),
                CString::new("1").unwrap(),
            );
            assert_eq!(fc_set_option(sim, k.as_ptr(), v.as_ptr()), FC_ERR_ARG);
            let (e, p) = (CString::new("R99").unwrap(), CString::new("value").unwrap());
            assert_eq!(
                fc_set_param(sim, e.as_ptr(), p.as_ptr(), 1.0),
                FC_ERR_NOT_FOUND
            );
            // Signal before any run.
            let mut d: *const c_double = std::ptr::null();
            let mut n: usize = 0;
            let key = CString::new("V(out)").unwrap();
            assert_eq!(fc_signal(sim, key.as_ptr(), &mut d, &mut n), FC_ERR_ARG);
            assert_eq!(fc_run_tran(sim, 1e-9, 100e-9), FC_OK);
            let bogus = CString::new("V(nope)").unwrap();
            assert_eq!(
                fc_signal(sim, bogus.as_ptr(), &mut d, &mut n),
                FC_ERR_NOT_FOUND
            );
            // Non-positive step.
            assert_eq!(fc_run_tran(sim, 0.0, 1e-9), FC_ERR_ARG);
            assert!(fc_stepper_new(sim, -1.0).is_null());

            let st = fc_stepper_new(sim, 1e-10);
            let nope = CString::new("vnope").unwrap();
            assert_eq!(fc_set_source(st, nope.as_ptr(), 1.0), FC_ERR_SIM);
            let nonode = CString::new("nonode").unwrap();
            let mut x = 0.0;
            assert_eq!(fc_get_node(st, nonode.as_ptr(), &mut x), FC_ERR_NOT_FOUND);
            assert_eq!(
                fc_get_node(
                    st,
                    CString::new("out").unwrap().as_ptr(),
                    std::ptr::null_mut()
                ),
                FC_ERR_ARG
            );
            fc_stepper_free(st);
            fc_sim_free(sim);
        }
    }

    #[test]
    fn node_name_reports_required_size_and_truncates() {
        unsafe {
            let sim = load(RC);
            let st = fc_stepper_new(sim, 1e-10);
            assert!(fc_node_count(st) >= 2);
            // NULL buffer just measures.
            let need = fc_node_name(st, 0, std::ptr::null_mut(), 0);
            assert!(need > 0);
            let mut buf = vec![0 as c_char; need as usize + 1];
            assert_eq!(fc_node_name(st, 0, buf.as_mut_ptr(), buf.len()), need);
            assert_eq!(CStr::from_ptr(buf.as_ptr()).to_bytes().len(), need as usize);
            // Truncation still NUL-terminates.
            let mut small = vec![0 as c_char; 2];
            assert_eq!(fc_node_name(st, 0, small.as_mut_ptr(), 2), need);
            assert_eq!(CStr::from_ptr(small.as_ptr()).to_bytes().len(), 1);
            assert_eq!(fc_node_name(st, 9999, std::ptr::null_mut(), 0), -1);
            fc_stepper_free(st);
            fc_sim_free(sim);
        }
    }

    #[test]
    fn pwl_injection_validates_and_drives() {
        unsafe {
            let sim = load(RC);
            let name = CString::new("V1").unwrap();
            let t = [0.0, 50e-9, 100e-9];
            let v = [0.0, 1.0, 1.0];
            assert_eq!(
                fc_set_source_pwl(sim, name.as_ptr(), t.as_ptr(), v.as_ptr(), 3),
                FC_OK
            );
            assert_eq!(fc_run_tran(sim, 1e-9, 100e-9), FC_OK);
            let mut d: *const c_double = std::ptr::null();
            let mut n: usize = 0;
            assert_eq!(
                fc_signal(
                    sim,
                    CString::new("V(out)").unwrap().as_ptr(),
                    &mut d,
                    &mut n
                ),
                FC_OK
            );
            let series = std::slice::from_raw_parts(d, n);
            assert!(
                (series[n - 1] - 1.0).abs() < 1e-2,
                "final {}",
                series[n - 1]
            );

            // Descending times and empty tables are rejected.
            let bad_t = [1e-9, 0.0];
            assert_eq!(
                fc_set_source_pwl(sim, name.as_ptr(), bad_t.as_ptr(), v.as_ptr(), 2),
                FC_ERR_ARG
            );
            assert_eq!(
                fc_set_source_pwl(sim, name.as_ptr(), t.as_ptr(), v.as_ptr(), 0),
                FC_ERR_ARG
            );
            assert_eq!(
                fc_set_source_pwl(sim, name.as_ptr(), std::ptr::null(), v.as_ptr(), 2),
                FC_ERR_ARG
            );
            fc_sim_free(sim);
        }
    }
}
