//! OsdiDevice: wraps a loaded OsdiDescriptor and implements the Device trait.
//!
//! Jacobian path: copy-based (write_jacobian_array_resist / write_jacobian_array_react).
//! Reactive entries are identified by react_ptr_off != u32::MAX in jacobian_entries;
//! write_jacobian_array_react writes values in that same traversal order.
//!
//! The aliasing path (load_jacobian_resist via jacobian_ptr_resist_offset) crashes
//! with the OpenVAF-compiled nmos_l1/pmos_l1 models — it likely requires additional
//! simulator-side initialization that we haven't determined. The copy path is sufficient.
//!
//! Safety invariant for load_residual / load_jacobian (&self):
//!   The OSDI functions called (load_spice_rhs_*, write_jacobian_array_*)
//!   MUST only READ from `inst` memory. Correctly implemented OSDI v0.4 models
//!   satisfy this.

use std::ffi::CStr;
use std::os::raw::c_void;
use std::sync::Arc;

use fairchild_core::device::{Device, EvalFlags, NodeId, SimContext};
use fairchild_core::mna::MnaMatrix;

use crate::ffi::{
    OsdiDescriptor, OsdiInitInfo, OsdiParamOpvar, OsdiSimInfo, OsdiSimParas, ANALYSIS_DC,
    ANALYSIS_TRAN, CALC_REACT_JACOBIAN, CALC_REACT_RESIDUAL, CALC_RESIST_JACOBIAN,
    CALC_RESIST_RESIDUAL, PARA_KIND_INST, PARA_KIND_MODEL,
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
    /// Solution vector from the last eval() call; used as prev_solve in eval.
    x_cache: Vec<f64>,
    /// Solution at the last committed (accepted) timestep; used as prev_solve in
    /// load_spice_rhs_tran so reactive history is pinned to the accepted state.
    /// Empty until commit_timestep() is first called.
    x_tprev: Vec<f64>,
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
        let model_u64s = (desc.model_size as usize).div_ceil(8);
        let inst_u64s = (desc.instance_size as usize).div_ceil(8);
        OsdiDevice {
            _lib: None,
            descriptor,
            model: vec![0u64; model_u64s.max(1)],
            instance: vec![0u64; inst_u64s.max(1)],
            mna_nodes: Vec::new(),
            x_cache: Vec::new(),
            x_tprev: Vec::new(),
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
        let model_u64s = (desc.model_size as usize).div_ceil(8);
        let inst_u64s = (desc.instance_size as usize).div_ceil(8);
        Some(OsdiDevice {
            _lib: Some(lib),
            descriptor,
            model: vec![0u64; model_u64s.max(1)],
            instance: vec![0u64; inst_u64s.max(1)],
            mna_nodes: Vec::new(),
            x_cache: Vec::new(),
            x_tprev: Vec::new(),
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

    /// Expose raw instance/model pointers for integration-test diagnostics.
    pub fn inst_ptr_raw(&self) -> *mut c_void {
        self.inst_ptr()
    }
    pub fn model_ptr_raw(&self) -> *mut c_void {
        self.model_ptr()
    }

    /// Read a model param's current value via access(READ), and report the pointer
    /// offset from model_ptr.  Returns (value, byte_offset_from_model_base) or None.
    pub fn probe_model_param(&self, name: &str) -> Option<(f64, isize)> {
        use crate::ffi::{ACCESS_FLAG_READ, PARA_KIND_MODEL};
        let desc = self.desc();
        let access_fn = desc.access?;
        let n_total = desc.num_params as usize;
        let n_inst = desc.num_instance_params as usize;
        if n_total == 0 || desc.param_opvar.is_null() {
            return None;
        }
        let params = unsafe { std::slice::from_raw_parts(desc.param_opvar, n_total) };
        // enumerate() preserves the absolute index — access() expects the absolute param_opvar index.
        for (i, param) in params.iter().enumerate().skip(n_inst) {
            if osdi_param_name_matches(param, name) {
                let id = PARA_KIND_MODEL | i as u32;
                let ptr = unsafe {
                    access_fn(std::ptr::null_mut(), self.model_ptr(), id, ACCESS_FLAG_READ)
                };
                if ptr.is_null() {
                    return None;
                }
                let value = unsafe { *(ptr as *const f64) };
                let offset =
                    unsafe { (ptr as *const u8).offset_from(self.model.as_ptr() as *const u8) };
                return Some((value, offset));
            }
        }
        None
    }
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
            let mut res = OsdiInitInfo {
                flags: 0,
                num_errors: 0,
                errors: std::ptr::null_mut(),
            };
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
        // Pre-size mna_nodes to num_nodes so internal slots are at least
        // ground (UINT32_MAX).  `bind_extra_nodes` later overwrites them
        // with real allocated row indices.
        let num_nodes = self.desc().num_nodes as usize;
        let num_terminals = self.desc().num_terminals as usize;
        self.mna_nodes = terminals
            .iter()
            .copied()
            .chain(std::iter::repeat_n(None, num_nodes.saturating_sub(terminals.len())))
            .take(num_nodes)
            .collect();

        // Cache all descriptor reads before taking the mutable borrow on instance.
        // (desc() borrows self; as_mut_ptr() is a conflicting &mut borrow.)
        let node_mapping_offset = self.desc().node_mapping_offset as usize;
        let setup_fn = self.desc().setup_instance;

        // Write the MNA↔OSDI node mapping into instance memory.
        // The model reads node_mapping[i] from (inst + node_mapping_offset) to
        // find which solution-vector index corresponds to its i-th node.
        // UINT32_MAX is the sentinel for ground (NodeId = None).
        // We write ALL num_nodes slots (terminals + internals).
        let map_ptr =
            unsafe { (self.instance.as_mut_ptr() as *mut u8).add(node_mapping_offset) as *mut u32 };
        for i in 0..num_nodes {
            let node = self.mna_nodes.get(i).copied().flatten();
            unsafe {
                *map_ptr.add(i) = node.map(|n| n as u32).unwrap_or(u32::MAX);
            }
        }

        if let Some(f) = setup_fn {
            let mut paras = null_sim_paras();
            let mut res = OsdiInitInfo {
                flags: 0,
                num_errors: 0,
                errors: std::ptr::null_mut(),
            };
            unsafe {
                f(
                    std::ptr::null_mut(),
                    self.inst_ptr(),
                    self.model_ptr(),
                    ctx.temperature,
                    num_terminals as u32,
                    &mut paras,
                    &mut res,
                );
            }
        }
    }

    fn num_extra_nodes(&self) -> usize {
        let desc = self.desc();
        (desc.num_nodes as usize).saturating_sub(desc.num_terminals as usize)
    }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        let num_terminals = self.desc().num_terminals as usize;
        let num_nodes = self.desc().num_nodes as usize;
        // mna_nodes is sized to num_nodes already; overwrite the trailing
        // internal slots with the allocated MNA row indices.
        for i in num_terminals..num_nodes {
            let offset = i - num_terminals;
            self.mna_nodes[i] = Some(first_idx + offset);
        }
        self.refresh_instance();
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
            let mut osdi_flags = if flags.transient {
                ANALYSIS_TRAN
            } else {
                ANALYSIS_DC
            };
            if flags.resistive {
                osdi_flags |= CALC_RESIST_RESIDUAL | CALC_RESIST_JACOBIAN;
            }
            if flags.transient {
                osdi_flags |= CALC_REACT_RESIDUAL | CALC_REACT_JACOBIAN;
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
                f(
                    self.inst_ptr(),
                    self.model_ptr(),
                    temp.as_mut_ptr().add(1),
                    prev,
                );
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
            std::slice::from_raw_parts(desc.jacobian_entries, desc.num_jacobian_entries as usize)
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

    fn load_residual_tran(&self, b: &mut [f64], alpha: f64) {
        if let Some(f) = self.desc().load_spice_rhs_tran {
            let mut temp = vec![0.0f64; b.len() + 1];
            // Use x_tprev (previous accepted timestep) as the reactive history term.
            // Fall back to x_cache (current iterate) when x_tprev is not yet populated
            // (i.e., before the first commit_timestep call in the fixed-step path).
            let prev_src = if self.x_tprev.is_empty() {
                &self.x_cache
            } else {
                &self.x_tprev
            };
            let prev = unsafe { prev_src.as_ptr().add(1) as *mut f64 };
            unsafe {
                f(
                    self.inst_ptr(),
                    self.model_ptr(),
                    temp.as_mut_ptr().add(1),
                    prev,
                    alpha,
                );
            }
            for i in 0..b.len() {
                b[i] += temp[i + 1];
            }
        } else {
            self.load_residual(b);
        }
    }

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, alpha: f64) {
        // Resistive part (same as DC).
        self.load_jacobian(mat);

        // Reactive part: stamp alpha * C entries.
        // write_jacobian_array_react writes n_react values in the traversal order of
        // jacobian_entries where react_ptr_off != u32::MAX. We match that order.
        let desc = self.desc();
        let n_react = desc.num_reactive_jacobian_entries as usize;
        if n_react == 0 {
            return;
        }
        let f = match desc.write_jacobian_array_react {
            Some(f) => f,
            None => return,
        };

        let mut jac_buf = vec![0.0f64; n_react];
        unsafe {
            f(self.inst_ptr(), self.model_ptr(), jac_buf.as_mut_ptr());
        }

        let n_total = desc.num_jacobian_entries as usize;
        let entries = unsafe { std::slice::from_raw_parts(desc.jacobian_entries, n_total) };

        // Walk all entries; for each one with a reactive pointer (react_ptr_off != MAX),
        // consume the next value from jac_buf in order.
        let mut react_idx = 0;
        for entry in entries.iter() {
            if entry.react_ptr_off == u32::MAX {
                continue;
            }
            if react_idx >= n_react {
                break;
            }
            let osdi_r = entry.nodes.node_1 as usize;
            let osdi_c = entry.nodes.node_2 as usize;
            if let (Some(mr), Some(mc)) = (
                self.mna_nodes.get(osdi_r).copied().flatten(),
                self.mna_nodes.get(osdi_c).copied().flatten(),
            ) {
                mat.a[mr][mc] += alpha * jac_buf[react_idx];
            }
            react_idx += 1;
        }
    }

    fn commit_timestep(&mut self, x: &[f64]) {
        // Snapshot x (with guard element prepended) into x_tprev.
        // The guard at x_tprev[0] = 0.0 handles the OSDI ground-sentinel (-1 index).
        self.x_tprev.resize(x.len() + 1, 0.0);
        self.x_tprev[0] = 0.0;
        self.x_tprev[1..].copy_from_slice(x);
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        let desc = self.desc();
        let access_fn = match desc.access {
            Some(f) => f,
            None => return false,
        };
        // param_opvar layout: [inst_params(0..n_inst) | model_params(n_inst..n_total) | opvars]
        // The access() id is the ABSOLUTE param_opvar index (not relative within kind).
        let n_total = desc.num_params as usize;
        let n_inst = desc.num_instance_params as usize;
        if n_total == 0 || desc.param_opvar.is_null() {
            return false;
        }
        let params = unsafe { std::slice::from_raw_parts(desc.param_opvar, n_total) };

        // Instance params: absolute indices 0..n_inst, kind = PARA_KIND_INST.
        for (i, param) in params.iter().enumerate().take(n_inst) {
            if osdi_param_name_matches(param, name) {
                let id = PARA_KIND_INST | i as u32;
                let ptr = unsafe {
                    access_fn(
                        self.inst_ptr(),
                        self.model_ptr(),
                        id,
                        crate::ffi::ACCESS_FLAG_SET,
                    )
                };
                if !ptr.is_null() {
                    unsafe {
                        *(ptr as *mut f64) = value;
                    }
                    return true;
                }
            }
        }

        // Model params: absolute indices n_inst..n_total, kind = PARA_KIND_MODEL.
        for (i, param) in params.iter().enumerate().skip(n_inst) {
            if osdi_param_name_matches(param, name) {
                let id = PARA_KIND_MODEL | i as u32;
                let ptr = unsafe {
                    access_fn(
                        std::ptr::null_mut(),
                        self.model_ptr(),
                        id,
                        crate::ffi::ACCESS_FLAG_SET,
                    )
                };
                if !ptr.is_null() {
                    unsafe {
                        *(ptr as *mut f64) = value;
                    }
                    // Re-run setup_instance so the instance struct picks up the new model value.
                    // OpenVAF caches model-param-derived quantities in the instance during
                    // setup_instance; eval() reads from instance, not directly from model.
                    self.refresh_instance();
                    return true;
                }
            }
        }

        false
    }
}

impl OsdiDevice {
    /// Re-run the OSDI setup_instance call with the current mna_nodes and model state.
    /// Required after writing a model param via access(SET) so that setup_instance can
    /// propagate the new value into the instance struct where eval() reads it.
    fn refresh_instance(&mut self) {
        let node_mapping_offset = self.desc().node_mapping_offset as usize;
        let setup_fn = self.desc().setup_instance;
        // OSDI's `num_terminals` argument is the count of *external* nodes;
        // internal flow-branch nodes are implicit (the rest of num_nodes).
        let num_terminals = self.desc().num_terminals;
        let num_nodes = self.desc().num_nodes as usize;
        let temperature = SimContext::default().temperature;

        let map_ptr =
            unsafe { (self.instance.as_mut_ptr() as *mut u8).add(node_mapping_offset) as *mut u32 };
        for i in 0..num_nodes {
            let node = self.mna_nodes.get(i).copied().flatten();
            unsafe {
                *map_ptr.add(i) = node.map(|n| n as u32).unwrap_or(u32::MAX);
            }
        }

        if let Some(f) = setup_fn {
            let mut paras = null_sim_paras();
            let mut res = OsdiInitInfo {
                flags: 0,
                num_errors: 0,
                errors: std::ptr::null_mut(),
            };
            unsafe {
                f(
                    std::ptr::null_mut(),
                    self.inst_ptr(),
                    self.model_ptr(),
                    temperature,
                    num_terminals,
                    &mut paras,
                    &mut res,
                );
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Returns true if the primary name or any alias in `param` matches `target`.
///
/// OSDI convention: `name` is an array of `num_alias + 1` C-string pointers.
/// `name[0]` is the primary name; `name[1..=num_alias]` are aliases.
fn osdi_param_name_matches(param: &OsdiParamOpvar, target: &str) -> bool {
    if param.name.is_null() {
        return false;
    }
    let n = param.num_alias as usize + 1; // primary + aliases
    let names = unsafe { std::slice::from_raw_parts(param.name, n) };
    for &name_ptr in names {
        if name_ptr.is_null() {
            continue;
        }
        let s = unsafe { CStr::from_ptr(name_ptr) };
        if s.to_str().unwrap_or("").eq_ignore_ascii_case(target) {
            return true;
        }
    }
    false
}

fn null_sim_paras() -> OsdiSimParas {
    OsdiSimParas {
        names: std::ptr::null_mut(),
        vals: std::ptr::null_mut(),
        names_str: std::ptr::null_mut(),
        vals_str: std::ptr::null_mut(),
    }
}
