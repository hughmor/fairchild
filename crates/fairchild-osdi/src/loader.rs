use std::ffi::{CString, CStr};
use std::path::Path;
use std::sync::Arc;

use fairchild_core::device_registry::DeviceRegistry;
use fairchild_core::Device;

use crate::device::OsdiDevice;
use crate::error::OsdiError;
use crate::ffi::OsdiDescriptor;

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
            let msg = CStr::from_ptr(libc::dlerror()).to_string_lossy().into_owned();
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
        let descriptors_base =
            dlsym_or_err(handle, b"OSDI_DESCRIPTORS\0")? as *const u8;

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
            registry.register(name, move |terminals, ctx| {
                let mut dev = OsdiDevice::from_library(Arc::clone(&lib), i)
                    .expect("descriptor index must be valid");
                dev.setup_model(ctx);
                dev.setup_instance(terminals, ctx);
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
        let sym_name: &'static str =
            std::str::from_utf8(&name[..name.len() - 1]).unwrap_or("?");
        Err(OsdiError::Symbol { symbol: sym_name, detail })
    } else {
        Ok(ptr)
    }
}
