//! OsdiDevice: wraps a loaded OsdiDescriptor and implements the Device trait.
//!
//! Uses the copy-based Jacobian path (write_jacobian_array_resist) rather than
//! the pointer-aliasing optimisation (jacobian_ptr_resist_offset). The aliasing
//! path requires MnaMatrix to have stable storage across NR iterations; the
//! current Vec<Vec<f64>> layout reallocates every iteration. Upgrade to the
//! aliasing path once MnaMatrix switches to a single contiguous buffer.
//!
//! Safety invariant for load_residual / load_jacobian (&self):
//!   The OSDI functions called (load_spice_rhs_dc, write_jacobian_array_resist)
//!   MUST only READ from `inst` memory. Correctly implemented OSDI v0.4 models
//!   satisfy this. If a model violates it, bump the trait methods to `&mut self`.

use std::os::raw::c_void;
use std::sync::Arc;

use fairchild_core::device::{Device, EvalFlags, NodeId, SimContext};
use fairchild_core::mna::MnaMatrix;

use crate::ffi::{
    OsdiDescriptor, OsdiInitInfo, OsdiSimInfo, OsdiSimParas,
    ANALYSIS_DC, CALC_RESIST_JACOBIAN, CALC_RESIST_RESIDUAL,
};
use crate::loader::OsdiLibrary;

/// Device backed by an OSDI v0.4 descriptor loaded from a `.osdi` shared library.
pub struct OsdiDevice {
    /// Keeps the library alive so the descriptor pointer stays valid.
    _lib: Option<Arc<OsdiLibrary>>,
    /// Stable pointer into the OsdiLibrary descriptor array.
    descriptor: *const OsdiDescriptor,
    /// Per-model state: descriptor.model_size bytes, 8-byte aligned.
    model: Vec<u64>,
    /// Per-instance state: descriptor.instance_size bytes, 8-byte aligned.
    instance: Vec<u64>,
    /// MNA solution-vector index for each OSDI terminal. None = ground.
    mna_nodes: Vec<NodeId>,
    /// Solution vector from the last eval() call; forwarded as prev_solve.
    x_cache: Vec<f64>,
}

// SAFETY: descriptor is read-only after construction; Vec storage is thread-safe.
unsafe impl Send for OsdiDevice {}
unsafe impl Sync for OsdiDevice {}

impl OsdiDevice {
    /// Construct an OsdiDevice from a descriptor pointer.
    ///
    /// Construct from a raw descriptor pointer (for tests / mock usage).
    ///
    /// # Safety
    /// `descriptor` must point to a valid `OsdiDescriptor` that remains valid
    /// for the entire lifetime of the returned `OsdiDevice`.
    pub unsafe fn new(descriptor: *const OsdiDescriptor) -> Self {
        let desc = &*descriptor;
        let model_u64s = ((desc.model_size as usize) + 7) / 8;
        let inst_u64s = ((desc.instance_size as usize) + 7) / 8;
        OsdiDevice {
            _lib: None,
            descriptor,
            model: vec![0u64; model_u64s.max(1)],
            instance: vec![0u64; inst_u64s.max(1)],
            mna_nodes: Vec::new(),
            x_cache: Vec::new(),
        }
    }

    /// Construct from a loaded library by descriptor index.
    ///
    /// The returned `OsdiDevice` co-owns an `Arc<OsdiLibrary>`, ensuring the
    /// library stays loaded for the device's lifetime.
    ///
    /// Returns `None` if `model_index` is out of range.
    pub fn from_library(lib: Arc<OsdiLibrary>, model_index: usize) -> Option<Self> {
        let descriptor = lib.descriptors().nth(model_index)? as *const OsdiDescriptor;
        let desc = unsafe { &*descriptor };
        let model_u64s = ((desc.model_size as usize) + 7) / 8;
        let inst_u64s = ((desc.instance_size as usize) + 7) / 8;
        Some(OsdiDevice {
            _lib: Some(lib),
            descriptor,
            model: vec![0u64; model_u64s.max(1)],
            instance: vec![0u64; inst_u64s.max(1)],
            mna_nodes: Vec::new(),
            x_cache: Vec::new(),
        })
    }

    #[inline]
    fn desc(&self) -> &OsdiDescriptor {
        unsafe { &*self.descriptor }
    }

    /// Raw *mut to model buffer (C-compatible; model functions must not alias).
    #[inline]
    fn model_ptr(&self) -> *mut c_void {
        self.model.as_ptr() as *mut c_void
    }

    /// Raw *mut to instance buffer (C-compatible; load_* must only read it).
    #[inline]
    fn inst_ptr(&self) -> *mut c_void {
        self.instance.as_ptr() as *mut c_void
    }

    /// Test-only: expose raw pointers for diagnostic assertions.
    #[cfg(test)]
    pub fn inst_ptr_test(&self) -> *mut c_void { self.inst_ptr() }
    #[cfg(test)]
    pub fn model_ptr_test(&self) -> *mut c_void { self.model_ptr() }
}

impl Device for OsdiDevice {
    fn num_terminals(&self) -> usize {
        self.desc().num_terminals as usize
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        // Cache function pointer before any mutable borrow.
        let setup_fn = self.desc().setup_model;
        if let Some(f) = setup_fn {
            let mut paras = null_sim_paras();
            let mut res =
                OsdiInitInfo { flags: 0, num_errors: 0, errors: std::ptr::null_mut() };
            unsafe {
                f(
                    std::ptr::null_mut(), // handle (unused by most models)
                    self.model_ptr(),
                    &mut paras,
                    &mut res,
                );
            }
        }
        let _ = ctx; // temperature injection deferred — will go into OsdiSimParas
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        self.mna_nodes = terminals.to_vec();

        // Cache all descriptor reads before taking the mutable borrow on instance.
        // (desc() borrows self; as_mut_ptr() is a conflicting &mut borrow.)
        let node_mapping_offset = self.desc().node_mapping_offset as usize;
        let setup_fn = self.desc().setup_instance;

        // Write the MNA↔OSDI node mapping into instance memory.
        // The model reads node_mapping[i] from (inst + node_mapping_offset) to
        // find which solution-vector index corresponds to its i-th node.
        // UINT32_MAX is the sentinel for ground (NodeId = None).
        let map_ptr = unsafe {
            (self.instance.as_mut_ptr() as *mut u8).add(node_mapping_offset) as *mut u32
        };
        for (i, &node) in terminals.iter().enumerate() {
            unsafe {
                *map_ptr.add(i) = node.map(|n| n as u32).unwrap_or(u32::MAX);
            }
        }

        if let Some(f) = setup_fn {
            let mut paras = null_sim_paras();
            let mut res =
                OsdiInitInfo { flags: 0, num_errors: 0, errors: std::ptr::null_mut() };
            unsafe {
                f(
                    std::ptr::null_mut(),
                    self.inst_ptr(),
                    self.model_ptr(),
                    ctx.temperature,
                    terminals.len() as u32,
                    &mut paras,
                    &mut res,
                );
            }
        }
    }

    fn eval(&mut self, x: &[f64], flags: EvalFlags, _ctx: &SimContext) {
        // Prepend a 0.0 guard element so OSDI's ldpsw sign-extension of the
        // ground sentinel (UINT32_MAX → -1 as i64) reads x_cache[0] = 0.0
        // (ground voltage) instead of whatever is 8 bytes before the buffer.
        self.x_cache.resize(x.len() + 1, 0.0);
        self.x_cache[0] = 0.0;
        self.x_cache[1..].copy_from_slice(x);

        let eval_fn = self.desc().eval;
        if let Some(f) = eval_fn {
            let mut osdi_flags = ANALYSIS_DC;
            if flags.resistive {
                osdi_flags |= CALC_RESIST_RESIDUAL | CALC_RESIST_JACOBIAN;
            }
            let mut info = OsdiSimInfo {
                paras: null_sim_paras(),
                abstime: 0.0,
                // Pass ptr+1: x_cache[1..] mirrors x[0..], and x_cache[0]=0.0
                // acts as a guard for OSDI's out-of-bounds index -1 (ground).
                prev_solve: unsafe { self.x_cache.as_ptr().add(1) as *mut f64 },
                prev_state: std::ptr::null_mut(),
                next_state: std::ptr::null_mut(),
                flags: osdi_flags,
            };
            unsafe {
                f(
                    std::ptr::null_mut(),
                    self.inst_ptr(),
                    self.model_ptr(),
                    &mut info,
                );
            }
        }
    }

    fn load_residual(&self, b: &mut [f64]) {
        if let Some(f) = self.desc().load_spice_rhs_dc {
            // Padded dst buffer: pass ptr+1 so OSDI's dst[-1] write (ground
            // sentinel index -1 from ldpsw) lands in temp[0] not heap metadata.
            let mut temp = vec![0.0f64; b.len() + 1];
            let prev = unsafe { self.x_cache.as_ptr().add(1) as *mut f64 };
            unsafe {
                f(self.inst_ptr(), self.model_ptr(), temp.as_mut_ptr().add(1), prev);
            }
            for i in 0..b.len() {
                b[i] += temp[i + 1];
            }
        }
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let desc = self.desc();
        let n_resist = desc.num_resistive_jacobian_entries as usize;
        if n_resist == 0 {
            return;
        }

        let f = desc.write_jacobian_array_resist.expect(
            "OsdiDevice: write_jacobian_array_resist is None — \
             aliasing Jacobian path not implemented; model must provide this function",
        );

        // Copy Jacobian values out of instance memory into a temp buffer.
        let mut jac_buf = vec![0.0f64; n_resist];
        unsafe {
            f(self.inst_ptr(), self.model_ptr(), jac_buf.as_mut_ptr());
        }

        // jacobian_entries[0..num_resistive_jacobian_entries] are the resistive entries,
        // in the same order that write_jacobian_array_resist writes to jac_buf.
        let entries = unsafe {
            std::slice::from_raw_parts(
                desc.jacobian_entries,
                desc.num_jacobian_entries as usize,
            )
        };

        for (i, entry) in entries.iter().take(n_resist).enumerate() {
            let osdi_r = entry.nodes.node_1 as usize;
            let osdi_c = entry.nodes.node_2 as usize;
            // Map OSDI terminal index → MNA matrix row/col (skip ground).
            if let (Some(mr), Some(mc)) = (
                self.mna_nodes.get(osdi_r).copied().flatten(),
                self.mna_nodes.get(osdi_c).copied().flatten(),
            ) {
                mat.a[mr][mc] += jac_buf[i];
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

fn null_sim_paras() -> OsdiSimParas {
    OsdiSimParas {
        names: std::ptr::null_mut(),
        vals: std::ptr::null_mut(),
        names_str: std::ptr::null_mut(),
        vals_str: std::ptr::null_mut(),
    }
}
