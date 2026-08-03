use std::ffi::{CStr, CString};
use std::path::Path;
use std::sync::Arc;

use fairchild_core::device_registry::{DeviceRegistry, ParamSet};
use fairchild_core::Device;

use crate::device::OsdiDevice;
use crate::error::OsdiError;
use crate::ffi::{OsdiDescriptor, OsdiLimFunction};

/// A loaded `.osdi` shared library (RAII wrapper around dlopen).
///
/// The library is kept open for the lifetime of this struct; all pointers
/// into the descriptor array are valid as long as `OsdiLibrary` is alive.
pub struct OsdiLibrary {
    handle: *mut libc::c_void,
    pub version: (u32, u32),
    pub num_descriptors: usize,
    pub descriptor_size: usize,
    /// Byte pointer to the first element of `OSDI_DESCRIPTORS`.
    descriptors_base: *const u8,
}

// SAFETY: The handle is only used through shared references that read stable,
// library-lifetime data. Callers are responsible for not racing on model state.
unsafe impl Send for OsdiLibrary {}
unsafe impl Sync for OsdiLibrary {}

impl OsdiLibrary {
    /// Load an `.osdi` shared library from `path`.
    ///
    /// # Safety
    /// Loading a shared library executes its static initialisers (arbitrary
    /// foreign code). The caller must ensure the library is a well-formed
    /// OSDI v0.4 shared object.
    pub unsafe fn open(path: &Path) -> Result<Self, OsdiError> {
        let path_c = path_to_cstring(path)?;

        let handle = libc::dlopen(path_c.as_ptr(), libc::RTLD_LAZY | libc::RTLD_LOCAL);
        if handle.is_null() {
            let msg = CStr::from_ptr(libc::dlerror())
                .to_string_lossy()
                .into_owned();
            return Err(OsdiError::DlOpen(msg));
        }

        match Self::init(handle) {
            Ok(lib) => Ok(lib),
            Err(e) => {
                libc::dlclose(handle);
                Err(e)
            }
        }
    }

    /// Finish initialisation after the handle is opened.
    unsafe fn init(handle: *mut libc::c_void) -> Result<Self, OsdiError> {
        macro_rules! sym_u32 {
            ($name:literal) => {{
                let ptr = dlsym_or_err(handle, $name)?;
                *(ptr as *const u32)
            }};
        }

        let major = sym_u32!(b"OSDI_VERSION_MAJOR\0");
        let minor = sym_u32!(b"OSDI_VERSION_MINOR\0");

        if major != 0 || minor < 4 {
            return Err(OsdiError::Version { major, minor });
        }

        let num_descriptors = sym_u32!(b"OSDI_NUM_DESCRIPTORS\0") as usize;
        let descriptor_size = sym_u32!(b"OSDI_DESCRIPTOR_SIZE\0") as usize;

        if descriptor_size < std::mem::size_of::<OsdiDescriptor>() {
            return Err(OsdiError::DescriptorSizeMismatch {
                expected: std::mem::size_of::<OsdiDescriptor>(),
                got: descriptor_size,
            });
        }

        // dlsym returns the address of the first element of OSDI_DESCRIPTORS.
        let descriptors_base = dlsym_or_err(handle, b"OSDI_DESCRIPTORS\0")? as *const u8;

        install_lim_table(handle);

        Ok(Self {
            handle,
            version: (major, minor),
            num_descriptors,
            descriptor_size,
            descriptors_base,
        })
    }

    /// Register every descriptor in this library into a `DeviceRegistry`.
    ///
    /// Each descriptor is keyed by its `name` field. The factory closure
    /// co-owns the library via `Arc`, so the library stays loaded as long as
    /// any of its devices are alive.
    pub fn register_into(self: &Arc<Self>, registry: &mut DeviceRegistry) {
        for (i, desc) in self.descriptors().enumerate() {
            let name = unsafe { CStr::from_ptr(desc.name) }
                .to_str()
                .unwrap_or("")
                .to_string();
            if name.is_empty() {
                continue;
            }
            let lib = Arc::clone(self);
            registry.register(name, move |terminals, params: &ParamSet, ctx| {
                let mut dev = OsdiDevice::from_library(Arc::clone(&lib), i)
                    .expect("descriptor index must be valid");
                dev.setup_model(ctx);
                dev.setup_instance(terminals, ctx);
                params.apply(&mut dev);
                Box::new(dev)
            });
        }
    }

    /// Iterate over the descriptors, striding by `descriptor_size` bytes.
    ///
    /// Striding by the exported size (rather than `sizeof`) keeps the iterator
    /// forward-compatible if a future OSDI version appends fields.
    pub fn descriptors(&self) -> impl Iterator<Item = &OsdiDescriptor> + '_ {
        (0..self.num_descriptors).map(move |i| unsafe {
            let ptr = self.descriptors_base.add(i * self.descriptor_size) as *const OsdiDescriptor;
            &*ptr
        })
    }
}

impl Drop for OsdiLibrary {
    fn drop(&mut self) {
        unsafe { libc::dlclose(self.handle) };
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn path_to_cstring(path: &Path) -> Result<CString, OsdiError> {
    use std::os::unix::ffi::OsStrExt;
    CString::new(path.as_os_str().as_bytes())
        .map_err(|_| OsdiError::DlOpen(format!("path contains null byte: {path:?}")))
}

unsafe fn dlsym_or_err(
    handle: *mut libc::c_void,
    name: &'static [u8],
) -> Result<*mut libc::c_void, OsdiError> {
    // Clear any pending error before calling dlsym.
    libc::dlerror();
    let ptr = libc::dlsym(handle, name.as_ptr() as *const libc::c_char);
    if ptr.is_null() {
        let err_ptr = libc::dlerror();
        let detail = if err_ptr.is_null() {
            "symbol resolved to null".to_owned()
        } else {
            CStr::from_ptr(err_ptr).to_string_lossy().into_owned()
        };
        // name is always a 'static byte literal from our call sites, so the
        // resulting &str is also 'static.
        let sym_name: &'static str = std::str::from_utf8(&name[..name.len() - 1]).unwrap_or("?");
        Err(OsdiError::Symbol {
            symbol: sym_name,
            detail,
        })
    } else {
        Ok(ptr)
    }
}

// ── $limit() support ────────────────────────────────────────────────────────

/// Limiting functions fairchild implements, keyed the way OSDI asks for them.
///
/// `num_args` is the count of arguments *beyond* the fixed `(init, check,
/// vnew, vold)` prefix — i.e. what the model passes after the name, and what
/// OpenVAF records in the table (`num_args - 2` of the lowered call's arity).
/// It is part of the ABI contract: install a 6-parameter function against an
/// entry the model calls with 5 and the callee reads a register that was never
/// set. So the arity is matched, not assumed.
///
/// Adding a limiter is one row here plus its `extern "C"` function. There is no
/// fixed set to complete: OpenVAF does not validate the name at all — it
/// forwards whatever string literal the model wrote (see
/// `hir_lower/src/expr.rs`, `BuiltIn::limit`), so this list can never be
/// exhaustive and [`osdi_no_limit`] is what makes that safe.
static LIMITERS: &[(&str, u32, LimiterFn)] = &[("pnjlim", 2, osdi_pnjlim)];

/// Signature shared by every entry in [`LIMITERS`].
///
/// Declared with the widest argument list we implement; entries with fewer
/// extra arguments simply ignore the tail. `osdi_no_limit` relies on the same
/// property — see its docs.
type LimiterFn = unsafe extern "C" fn(bool, *mut bool, f64, f64, f64, f64) -> f64;

/// Install fairchild's limiting functions into a freshly-loaded library's
/// `OSDI_LIM_TABLE`.
///
/// A library whose models call `$limit()` exports one table entry per distinct
/// limiting function, each with `func_ptr` **null**, because the simulator is
/// expected to supply the implementations. OpenVAF guards the call with the
/// `ENABLE_LIM` eval flag, so a null entry is harmless right up until the
/// simulator opts into limiting — which fairchild now does whenever a model
/// declares limit state. That makes "no entry is left null" a hard invariant
/// rather than a nicety: violate it and the first `eval` jumps to address 0.
///
/// Hence [`osdi_no_limit`] for anything unrecognised: an unimplemented limiter
/// degrades to *no limiting*, which costs convergence robustness but keeps the
/// model numerically valid, instead of killing the process.
unsafe fn install_lim_table(handle: *mut libc::c_void) {
    let table = libc::dlsym(handle, c"OSDI_LIM_TABLE".as_ptr()) as *mut OsdiLimFunction;
    let len_ptr = libc::dlsym(handle, c"OSDI_LIM_TABLE_LEN".as_ptr()) as *const u32;
    if table.is_null() || len_ptr.is_null() {
        return; // no model in this library calls $limit()
    }

    for entry in std::slice::from_raw_parts_mut(table, *len_ptr as usize) {
        let name = CStr::from_ptr(entry.name)
            .to_str()
            .unwrap_or("<invalid utf-8>");
        let chosen = LIMITERS
            .iter()
            .find(|(n, args, _)| *n == name && *args == entry.num_args);

        match chosen {
            Some((_, _, f)) => entry.func_ptr = *f as *mut libc::c_void,
            None => {
                // Distinguish "we have never heard of this" from "we have it
                // but the model calls it differently" — the second is a real
                // ABI mismatch and the more surprising of the two.
                let known_name = LIMITERS.iter().find(|(n, _, _)| *n == name);
                match known_name {
                    Some((_, expected, _)) => eprintln!(
                        "warning: OSDI model calls $limit(…, \"{name}\", …) with \
                         {} extra argument(s); fairchild implements it with {expected}. \
                         Running without limiting for this call — convergence may suffer.",
                        entry.num_args
                    ),
                    None => eprintln!(
                        "warning: OSDI model calls $limit(…, \"{name}\", …), which fairchild \
                         does not implement. Running without limiting for this call — \
                         convergence may suffer, results are unaffected. Add it to \
                         `LIMITERS` in fairchild-osdi/src/loader.rs to support it."
                    ),
                }
                entry.func_ptr = osdi_no_limit as LimiterFn as *mut libc::c_void;
            }
        }
    }
}

/// The identity limiter: hand back the proposed value untouched.
///
/// Stands in for any limiting function fairchild does not implement, so no
/// table entry is ever left null. Safe at *any* arity, which is what lets it
/// substitute for a function whose signature we do not know: it reads only the
/// four fixed leading parameters, and on every ABI fairchild targets
/// (AArch64 AAPCS, x86-64 SysV, Windows x64) surplus arguments the callee
/// never touches are harmless.
///
/// Deliberately does not set `*check`: it limited nothing, so it must not tell
/// the solver to withhold convergence.
unsafe extern "C" fn osdi_no_limit(
    _init: bool,
    _check: *mut bool,
    vnew: f64,
    _vold: f64,
    _a: f64,
    _b: f64,
) -> f64 {
    vnew
}

/// SPICE `pnjlim` — logarithmic compression of a forward junction step.
///
/// Matches `models::diode::ShockleyDiode::pnjlim` rather than SPICE3's exact
/// expression, deliberately: a Verilog-A junction and a native fairchild
/// junction then take the same path to the same answer. Not shared as one
/// function because the native one is a method carrying the device's own
/// `vcrit`, while OSDI passes `vcrit` per call.
///
/// Limiting changes the Newton path, not the solution — verified: the
/// `$limit`-using diode in `examples/verilog_a/models/va_diode.va` converges to
/// 0.6333213 V against 0.6333214 V for the same model without it.
unsafe extern "C" fn osdi_pnjlim(
    init: bool,
    check: *mut bool,
    vnew: f64,
    vold: f64,
    vt: f64,
    vcrit: f64,
) -> f64 {
    // `check` tells the solver this iterate was modified, so it must not
    // declare convergence on this step.
    let limited = |v: f64| {
        if !check.is_null() {
            *check = true;
        }
        v
    };
    if init {
        // First evaluation: no meaningful previous iterate to limit against,
        // so start at the critical voltage.
        return limited(vcrit);
    }
    if vnew > vcrit && (vnew - vold).abs() > 2.0 * vt {
        if vnew > vold {
            limited(vold + vt * ((vnew - vold) / vt + 1.0).ln())
        } else {
            limited(vold - vt * ((vold - vnew) / vt + 1.0).ln())
        }
    } else {
        vnew
    }
}
