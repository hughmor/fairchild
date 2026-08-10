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

use fairchild_core::device::Discretisation;
use fairchild_core::reactive::charge_current;

use crate::ffi::{
    OsdiDescriptor, OsdiInitInfo, OsdiParamOpvar, OsdiSimInfo, OsdiSimParas, ANALYSIS_DC,
    ANALYSIS_TRAN, CALC_REACT_JACOBIAN, CALC_REACT_RESIDUAL, CALC_RESIST_JACOBIAN,
    CALC_RESIST_RESIDUAL, ENABLE_LIM, INIT_LIM, PARA_KIND_INST, PARA_KIND_MODEL,
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
    /// Solution at the last committed (accepted) timestep.  Supplies the
    /// reactive history term in `load_residual_tran` — never the point the
    /// resistive residual is linearised about, which is always `x_cache`.
    /// Empty until commit_timestep() is first called.
    x_tprev: Vec<f64>,
    /// `$limit()` state, `descriptor.num_states` long each.  `eval` writes
    /// `next`; the following `eval` limits against it as `prev`, so the two
    /// swap after every call.
    ///
    /// These MUST be non-null whenever the model declares states: OpenVAF
    /// emits the state write unconditionally, so passing null segfaults on any
    /// model that calls `$limit` — which is every foundry compact model.
    lim_state_prev: Vec<f64>,
    lim_state_next: Vec<f64>,
    /// Cleared until the first `eval`, which runs with `INIT_LIM` to seed the
    /// state rather than limit against uninitialised memory.
    lim_initialised: bool,
    /// Reactive charge `q` at the last accepted timepoint, one slot per OSDI
    /// node. Read back through `load_residual_react`, which is the only
    /// sanctioned way to see a Verilog-A `ddt` contribution as a charge rather
    /// than as a Jacobian.
    ///
    /// Kept per *node* rather than per MNA row: the device has a handful of
    /// nodes and the circuit may have thousands of rows.
    q_prev: Vec<f64>,
    /// `q` one timepoint further back. `None` until two steps have been
    /// accepted, which is what gates BDF-2 for this device.
    q_prev2: Option<Vec<f64>>,
    /// Reactive current at the last accepted timepoint. Only Trapezoidal reads
    /// it; zero from a DC operating point, where the charge is static.
    i_react_prev: Vec<f64>,
    /// The integrator's discretisation, captured during `eval` because
    /// `load_*_tran` receives only `alpha`.
    disc: Option<Discretisation>,
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
            lim_state_prev: vec![0.0; desc.num_states as usize],
            lim_state_next: vec![0.0; desc.num_states as usize],
            lim_initialised: false,
            q_prev: vec![0.0; desc.num_nodes as usize],
            q_prev2: None,
            i_react_prev: vec![0.0; desc.num_nodes as usize],
            disc: None,
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
            lim_state_prev: vec![0.0; desc.num_states as usize],
            lim_state_next: vec![0.0; desc.num_states as usize],
            lim_initialised: false,
            q_prev: vec![0.0; desc.num_nodes as usize],
            q_prev2: None,
            i_react_prev: vec![0.0; desc.num_nodes as usize],
            disc: None,
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

    /// Write this instance's `$limit()` state indices into instance memory.
    ///
    /// The model reads `state_idx[k]` from `inst + state_idx_off` and then
    /// indexes `prev_state[state_idx[k]]` / `next_state[…]`. OSDI is written
    /// for a simulator holding one global state array and giving each instance
    /// a slice of it; we give each device its own two buffers instead, so the
    /// indices are just 0..num_states. Same contract, no shared allocation.
    ///
    /// Must be re-run whenever the node mapping is (see `refresh_instance`) —
    /// both live in the same instance struct that `setup_instance` reads.
    fn write_state_indices(&mut self) {
        let num_states = self.desc().num_states as usize;
        if num_states == 0 {
            return;
        }
        let off = self.desc().state_idx_off as usize;
        // i32: OpenVAF types state_idx with its generic `int` (LLVMInt32).
        let ptr = unsafe { (self.instance.as_mut_ptr() as *mut u8).add(off) as *mut i32 };
        for k in 0..num_states {
            unsafe { *ptr.add(k) = k as i32 };
        }
    }

    /// Call `f(mna_row, mna_col, dq/dx)` for every reactive Jacobian entry.
    ///
    /// `write_jacobian_array_react` writes `num_reactive_jacobian_entries`
    /// values in the traversal order of `jacobian_entries` filtered to those
    /// with `react_ptr_off != u32::MAX`; this walks that same order. Entries
    /// touching ground are dropped, as in the resistive path.
    ///
    /// Shared by the transient and frequency-domain stamps so the two cannot
    /// drift apart — they differ only in the factor applied (α vs ω).
    fn for_each_react_entry(&self, mut f: impl FnMut(usize, usize, f64)) {
        let desc = self.desc();
        let n_react = desc.num_reactive_jacobian_entries as usize;
        if n_react == 0 {
            return;
        }
        let Some(write) = desc.write_jacobian_array_react else {
            return;
        };

        let mut jac_buf = vec![0.0f64; n_react];
        unsafe {
            write(self.inst_ptr(), self.model_ptr(), jac_buf.as_mut_ptr());
        }

        let n_total = desc.num_jacobian_entries as usize;
        let entries = unsafe { std::slice::from_raw_parts(desc.jacobian_entries, n_total) };

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
                f(mr, mc, jac_buf[react_idx]);
            }
            react_idx += 1;
        }
    }

    /// Whether this model lets us build a charge-based companion — it declares
    /// reactive Jacobian entries *and* exposes the charge behind them.
    ///
    /// Both stamps consult this so they cannot disagree: stamping `scale·∂q/∂x`
    /// without the matching history is a companion with no history term, which
    /// is simply a wrong circuit. OpenVAF always emits both, so the `false` arm
    /// is for hand-written or partial libraries.
    fn has_charge_history(&self) -> bool {
        self.desc().num_reactive_jacobian_entries > 0 && self.desc().load_residual_react.is_some()
    }

    /// Read the device's reactive charge per OSDI node, via
    /// `load_residual_react`.
    ///
    /// `mna_len` is the MNA row count: OSDI writes into a solution-shaped
    /// buffer (same convention and same ground-sentinel guard as
    /// `load_spice_rhs_*`), so we gather from it into per-node slots.
    ///
    /// `None` when the model declares no reactive contribution, which is the
    /// common case and worth not allocating for.
    fn charge_per_node(&self, mna_len: usize) -> Option<Vec<f64>> {
        let f = self.desc().load_residual_react?;
        if self.desc().num_reactive_jacobian_entries == 0 {
            return None;
        }
        let mut dst = vec![0.0f64; mna_len + 1];
        unsafe {
            f(self.inst_ptr(), self.model_ptr(), dst.as_mut_ptr().add(1));
        }
        Some(
            self.mna_nodes
                .iter()
                .map(|n| n.map_or(0.0, |row| dst[row + 1]))
                .collect(),
        )
    }

    /// `Σ_j ∂q_row/∂x_j · x_j` per OSDI node, at the vector `x` currently
    /// cached by `eval`.
    ///
    /// This is what turns the charge companion into a SPICE-form RHS: the
    /// stamped Jacobian is `scale · ∂q/∂x`, so the row's linearisation
    /// contributes `scale · (∂q/∂x · x_k)` and the residual subtracts the
    /// actual current.
    fn react_jacobian_times(&self, x_padded: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0f64; self.mna_nodes.len()];
        // MNA row -> OSDI node, for scattering the row-indexed walk back.
        self.for_each_react_entry(|r, c, v| {
            for (node, slot) in self.mna_nodes.iter().enumerate() {
                if *slot == Some(r) {
                    out[node] += v * x_padded[c + 1];
                }
            }
        });
        out
    }

    /// Accumulate `scale * load_spice_rhs_{dc,tran}(prev_solve)` into `b`.
    ///
    /// OpenVAF's `load_spice_rhs_dc` writes `J_resist · prev_solve − f_resist`;
    /// `load_spice_rhs_tran` writes that *plus* `alpha · J_react · prev_solve`
    /// (it calls the DC body first — see `osdi/src/load.rs::load_spice_rhs`).
    /// `prev_solve` names the vector the linearisation is taken about, which is
    /// the current Newton iterate — not the previous timestep.
    fn accumulate_spice_rhs(&self, b: &mut [f64], x: &[f64], alpha: Option<f64>, scale: f64) {
        // Padded dst buffer: pass ptr+1 so OSDI's dst[-1] write (ground
        // sentinel index -1 from ldpsw) lands in temp[0] not heap metadata.
        let mut temp = vec![0.0f64; b.len() + 1];
        let prev = unsafe { x.as_ptr().add(1) as *mut f64 };
        unsafe {
            match alpha {
                Some(a) => match self.desc().load_spice_rhs_tran {
                    Some(f) => f(
                        self.inst_ptr(),
                        self.model_ptr(),
                        temp.as_mut_ptr().add(1),
                        prev,
                        a,
                    ),
                    None => return,
                },
                None => match self.desc().load_spice_rhs_dc {
                    Some(f) => f(
                        self.inst_ptr(),
                        self.model_ptr(),
                        temp.as_mut_ptr().add(1),
                        prev,
                    ),
                    None => return,
                },
            }
        }
        for i in 0..b.len() {
            b[i] += scale * temp[i + 1];
        }
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
            .chain(std::iter::repeat_n(
                None,
                num_nodes.saturating_sub(terminals.len()),
            ))
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
        self.write_state_indices();

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

    fn eval(&mut self, x: &[f64], flags: EvalFlags, ctx: &SimContext) {
        // `load_*_tran` gets only `alpha`, which cannot express Trapezoidal or
        // BDF-2; the method comes through the context instead.
        self.disc = ctx.discretisation;
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
            // `$limit()` limiting.  The state buffers must be real pointers
            // whenever the model declares states — OpenVAF emits the state
            // write unconditionally, so null segfaults.  INIT_LIM on the first
            // call seeds the state; after that we limit against the previous
            // iterate.
            if self.desc().num_states > 0 {
                osdi_flags |= ENABLE_LIM;
                if !self.lim_initialised {
                    osdi_flags |= INIT_LIM;
                }
            }
            let mut info = OsdiSimInfo {
                paras: null_sim_paras(),
                // `$abstime`.  SimContext::time_s is the transient clock the
                // solver advances; it stays 0 for DC/AC, which is what
                // Verilog-A expects there.
                abstime: ctx.time_s,
                // Pass ptr+1: x_cache[1..] mirrors x[0..], and x_cache[0]=0.0
                // acts as a guard for OSDI's out-of-bounds index -1 (ground).
                prev_solve: unsafe { self.x_cache.as_ptr().add(1) as *mut f64 },
                prev_state: self.lim_state_prev.as_mut_ptr(),
                next_state: self.lim_state_next.as_mut_ptr(),
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
            // What this eval limited to becomes what the next one limits
            // against.
            std::mem::swap(&mut self.lim_state_prev, &mut self.lim_state_next);
            self.lim_initialised = true;
        }
    }

    fn load_residual(&self, b: &mut [f64]) {
        self.accumulate_spice_rhs(b, &self.x_cache, None, 1.0);
    }

    /// The model's `white_noise()` / `flicker_noise()` contributions.
    ///
    /// OSDI splits this in two: the descriptor's `noise_sources` array names
    /// each generator and the node pair it drives, and `load_noise` fills an
    /// array of `num_noise_src` densities at a given frequency. So the node
    /// pairing is static and only the magnitude is evaluated per point, which
    /// is exactly the shape `Device::noise_sources` wants.
    ///
    /// Densities come back in A²/Hz — Verilog-A noise contributions are
    /// current contributions on a branch — which is what the caller expects,
    /// so nothing is converted.
    ///
    /// Before this existed a model could call `white_noise()`, compile, run,
    /// and contribute nothing to either noise analysis. That is worse than
    /// having no noise support at all: the model says it has noise and the
    /// simulator quietly disagrees.
    fn noise_sources(&self, _ctx: &SimContext, freq: f64) -> Vec<(NodeId, NodeId, f64)> {
        let desc = self.desc();
        let n = desc.num_noise_src as usize;
        let (Some(load), true) = (desc.load_noise, n > 0) else {
            return Vec::new();
        };
        let mut dens = vec![0.0f64; n];
        unsafe {
            load(self.inst_ptr(), self.model_ptr(), freq, dens.as_mut_ptr());
        }
        (0..n)
            .filter_map(|i| {
                let src = unsafe { &*desc.noise_sources.add(i) };
                let p = *self.mna_nodes.get(src.nodes.node_1 as usize)?;
                let q = *self.mna_nodes.get(src.nodes.node_2 as usize)?;
                let s = dens[i];
                (s.is_finite() && s > 0.0).then_some((p, q, s))
            })
            .collect()
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
        if self.desc().load_spice_rhs_dc.is_none() {
            self.load_residual(b);
            return;
        }
        // Resistive part: J_r·x_k − f_r(x_k), about the current iterate.
        //
        // Evaluating it about x_{n−1} instead (which this once did, by handing
        // OSDI's `tran` entry point x_tprev) leaves the Newton linearisation
        // inconsistent for any nonlinear model: it converged only while the
        // operating point sat still, and diverged on the first moving source.
        self.accumulate_spice_rhs(b, &self.x_cache, None, 1.0);

        // Reactive part, per node:  scale·(∂q/∂x · x_k) − i_n
        //
        // `charge_current` interprets the integration method — the same
        // function the native branch stamper uses, so a Verilog-A `ddt` and a
        // discrete C are integrated identically. Falling back to `alpha`
        // (Backward Euler) only when no discretisation reached us, which means
        // a caller outside the transient loop.
        if !self.has_charge_history() {
            // No charge to build a companion from. Fall back to the legacy
            // form, which derives the history from `load_spice_rhs_tran` and is
            // Backward Euler only — but is at least self-consistent with the
            // `alpha` the Jacobian stamp uses in the same situation.
            if self.desc().num_reactive_jacobian_entries > 0 {
                let hist = if self.x_tprev.is_empty() {
                    &self.x_cache
                } else {
                    &self.x_tprev
                };
                self.accumulate_spice_rhs(b, hist, Some(alpha), 1.0);
                self.accumulate_spice_rhs(b, hist, None, -1.0);
            }
            return;
        }
        let Some(q_new) = self.charge_per_node(b.len()) else {
            return;
        };
        let jq_x = self.react_jacobian_times(&self.x_cache);
        for (node, row) in self.mna_nodes.iter().enumerate() {
            let Some(row) = *row else { continue };
            let (i_n, scale) = match self.disc {
                Some(d) => charge_current(
                    d.mode,
                    d.h,
                    d.gear2_h_prev,
                    q_new[node],
                    self.q_prev[node],
                    self.q_prev2.as_ref().map(|q| q[node]),
                    self.i_react_prev[node],
                ),
                None => (alpha * (q_new[node] - self.q_prev[node]), alpha),
            };
            b[row] += scale * jq_x[node] - i_n;
        }
    }

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, alpha: f64) {
        // Resistive part (same as DC).
        self.load_jacobian(mat);
        // Reactive part: the method's coefficient times ∂q/∂x. `conductance`
        // is linear in the branch value, so evaluating it at 1.0 gives exactly
        // the scalar a matrix-valued charge needs — and keeps the Jacobian in
        // step with the residual above by construction.
        let scale = match self.disc.filter(|_| self.has_charge_history()) {
            Some(d) => fairchild_core::reactive::conductance(
                fairchild_core::device::ReactiveKind::Capacitor,
                1.0,
                d.mode,
                d.h,
                d.gear2_h_prev,
            ),
            None => alpha,
        };
        self.for_each_react_entry(|r, c, v| mat.a[r][c] += scale * v);
    }

    /// `.ac` / `.noise` counterpart of the reactive half of
    /// `load_jacobian_tran`: the same dq/dx entries in the same positions, for
    /// the caller to scale by ω instead of α.
    ///
    /// This is why `OsdiDevice` deliberately does not implement
    /// `small_signal_reactances` — a Verilog-A charge is a general matrix, and
    /// squeezing it into reciprocal two-terminal branches would silently drop
    /// transcapacitance (∂q_d/∂v_g ≠ ∂q_g/∂v_d), which is exactly what a
    /// BSIM-class model is made of.
    fn load_reactive_jacobian(&self, c_mat: &mut [Vec<f64>]) {
        self.for_each_react_entry(|r, c, v| c_mat[r][c] += v);
    }

    fn commit_timestep(&mut self, x: &[f64]) {
        // Snapshot x (with guard element prepended) into x_tprev.
        // The guard at x_tprev[0] = 0.0 handles the OSDI ground-sentinel (-1 index).
        self.x_tprev.resize(x.len() + 1, 0.0);
        self.x_tprev[0] = 0.0;
        self.x_tprev[1..].copy_from_slice(x);

        // Roll the charge history.
        //
        // The cached charge is from the last `eval`, which is one NR iterate
        // behind the converged solution — within `reltol`, but `reltol` is 1e-3
        // by default and this history feeds every later step. So correct it to
        // the converged point with the reactive Jacobian we already have:
        //
        //   q(x) ≈ q(x_eval) + ∂q/∂x · (x − x_eval)
        //
        // Exact for a linear charge (which is what makes a Verilog-A
        // `ddt(C*V)` match a discrete C bit-for-bit), first-order otherwise.
        // The native diode gets the same effect for free by recomputing its
        // charge analytically from the converged solution.
        let Some(mut q_now) = self.charge_per_node(x.len()) else {
            return;
        };
        let mut delta = vec![0.0f64; x.len() + 1];
        for (i, xi) in x.iter().enumerate() {
            delta[i + 1] = xi - self.x_cache.get(i + 1).copied().unwrap_or(0.0);
        }
        for (q, dq) in q_now.iter_mut().zip(self.react_jacobian_times(&delta)) {
            *q += dq;
        }
        // The current *entering* this accepted step becomes Trapezoidal's
        // history for the next one. Computed before q_prev is overwritten.
        if let Some(d) = self.disc {
            for node in 0..q_now.len() {
                let (i_n, _) = charge_current(
                    d.mode,
                    d.h,
                    d.gear2_h_prev,
                    q_now[node],
                    self.q_prev[node],
                    self.q_prev2.as_ref().map(|q| q[node]),
                    self.i_react_prev[node],
                );
                self.i_react_prev[node] = i_n;
            }
        }
        self.q_prev2 = Some(std::mem::replace(&mut self.q_prev, q_now));
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
        self.write_state_indices();

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
