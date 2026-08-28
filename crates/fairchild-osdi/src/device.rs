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

use std::ffi::{CStr, CString};
use std::os::raw::c_void;
use std::sync::Arc;

use fairchild_core::device::{Device, EvalFlags, NodeId, SimContext};
use fairchild_core::mna::MnaMatrix;
use fairchild_core::warn_user;

use fairchild_core::device::Discretisation;
use fairchild_core::reactive::charge_current;

use crate::ffi::{
    OsdiDescriptor, OsdiInitInfo, OsdiParamOpvar, OsdiSimInfo, OsdiSimParas, ANALYSIS_DC,
    ANALYSIS_TRAN, CALC_REACT_JACOBIAN, CALC_REACT_RESIDUAL, CALC_RESIST_JACOBIAN,
    CALC_RESIST_RESIDUAL, ENABLE_LIM, INIT_LIM,
};
use crate::loader::OsdiLibrary;

/// Device backed by an OSDI v0.4 descriptor loaded from a `.osdi` shared library.
pub struct OsdiDevice {
    /// Keeps the library alive so the descriptor pointer stays valid.
    _lib: Option<Arc<OsdiLibrary>>,
    /// Stable pointer into the OsdiLibrary descriptor array.
    descriptor: *const OsdiDescriptor,
    /// Per-model state: descriptor.model_size bytes, 8-byte aligned.
    ///
    /// 8 is enough, and the reference host's `max_align_t` (16) is not needed: an
    /// OpenVAF-generated struct holds doubles, 32-bit integers and pointers,
    /// none of which want more. This was worth checking once — a crash that
    /// looked like misalignment turned out to be [#42]'s log handle — so it is
    /// written down rather than re-tried: widening the element type without
    /// changing `div_ceil(8)` allocates half the bytes and corrupts state
    /// silently, which is what an experiment here looks like when it "works".
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
    /// The `$simparam` table, built once in `setup_model`. See [`SimParaTable`].
    sim_paras: Option<SimParaTable>,
    /// Whether the model's last `eval` returned `EVAL_RET_FLAG_FATAL`.
    ///
    /// The flags used to be discarded. `FATAL` is a model saying "this evaluation
    /// is not valid", and stamping the result anyway invents an answer out of
    /// numbers the model disowned. Reported through `Device::eval_status`.
    eval_fatal: bool,
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
    /// The OSDI `handle`: a NUL-terminated instance name the model dereferences
    /// while building a diagnostic. Passing null here segfaults any model that
    /// emits one — which every foundry compact model does about its defaults.
    handle: CString,
    /// Set only for a device elaborated from the bundle-port dialect: which of
    /// its terminals carry a wavelength label, and how a label moves between
    /// them. `None` for ordinary Verilog-A, which has no bundle structure to
    /// derive either from.
    bundle_lambda: Option<BundleLambda>,
    /// `wl_<k>` values the deck wrote on this instance, in the order they
    /// arrived. Recorded so `set_resolved_lambda` can say when a hand-set
    /// wavelength disagrees with the one the deck's own sources imply, instead
    /// of overwriting it in silence.
    wl_given: Vec<(usize, f64)>,
    /// The terminal rows the netlist gave this instance, kept because the node
    /// mapping is re-derived from scratch every time the model's collapse
    /// decisions could have changed.
    terminals: Vec<NodeId>,
    /// Whether `setup_instance` has run yet.
    ///
    /// Gates the re-setup inside `set_real_param`: OSDI's own order is *write
    /// the parameters, then set up*, so a parameter applied before setup needs
    /// no re-derivation and must not trigger one — calling `setup_instance`
    /// before the node mapping exists is not defined.
    setup_done: bool,
    /// First MNA row of this device's block of internal rows, once
    /// `bind_extra_nodes` has allocated them. `None` before that.
    extra_first: Option<usize>,
    /// The temperature `setup_instance` was last given.
    ///
    /// Re-running setup after a parameter write has to hand the model the same
    /// temperature the deck asked for. `SimContext::default()` was standing in
    /// here, so any deck with `.temp` had its devices re-derived at 300.15 K the
    /// moment a parameter was applied — every parameter a deck sets is applied
    /// after setup, so that was every device.
    temperature: f64,
}

/// Where a bundle-dialect model's wavelength labels live, computed once at
/// elaboration from the module header and the width the deck asked for.
///
/// This exists because λ is resolved before the solve, and resolution can only
/// speak for nets some device *declares*. Without it a bundle model in the
/// middle of an optical path leaves everything downstream of it unreached, and
/// an unreached port falls back to the band centre — a wrong wavelength with no
/// diagnostic, which is the exact failure mode this codebase refuses.
#[derive(Debug, Clone)]
pub struct BundleLambda {
    /// Terminals carrying a label — see `BundleModule::lambda_terminals`.
    pub terminals: Vec<usize>,
    /// `(from, to)` label routing — see `BundleModule::lambda_routing`.
    pub routing: Vec<(usize, usize)>,
    /// Channel count this instance was generated for.
    pub channels: usize,
    /// Wires per channel (3 or 5).
    pub wpc: usize,
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
            sim_paras: None,
            eval_fatal: false,
            q_prev: vec![0.0; desc.num_nodes as usize],
            q_prev2: None,
            i_react_prev: vec![0.0; desc.num_nodes as usize],
            disc: None,
            // Named after the module until the builder knows the instance name:
            // the model dereferences this while formatting a diagnostic, so it
            // must be a real string from the start.
            handle: descriptor_name(descriptor),
            bundle_lambda: None,
            wl_given: Vec::new(),
            terminals: Vec::new(),
            setup_done: false,
            extra_first: None,
            temperature: SimContext::default().temperature,
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
            sim_paras: None,
            eval_fatal: false,
            q_prev: vec![0.0; desc.num_nodes as usize],
            q_prev2: None,
            i_react_prev: vec![0.0; desc.num_nodes as usize],
            disc: None,
            // Named after the module until the builder knows the instance name:
            // the model dereferences this while formatting a diagnostic, so it
            // must be a real string from the start.
            handle: descriptor_name(descriptor),
            bundle_lambda: None,
            wl_given: Vec::new(),
            terminals: Vec::new(),
            setup_done: false,
            extra_first: None,
            temperature: SimContext::default().temperature,
        })
    }

    /// Attach the λ geometry of a bundle-dialect elaboration. Called by the
    /// registrar, before setup, because it changes nothing the setup reads.
    pub fn set_bundle_lambda(&mut self, lambda: BundleLambda) {
        self.bundle_lambda = Some(lambda);
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

    /// The MNA row each of this device's OSDI nodes ended up on, after node
    /// collapsing. `None` is the reference node.
    pub fn mna_nodes(&self) -> &[NodeId] {
        &self.mna_nodes
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
                // Same width rule as `set_real_param`: an `integer` parameter is
                // an i32, and reading it as f64 answers with garbage.
                let value = match param.flags & crate::ffi::PARA_TY_MASK {
                    crate::ffi::PARA_TY_INT => f64::from(unsafe { *(ptr as *const i32) }),
                    crate::ffi::PARA_TY_STR => return None,
                    _ => unsafe { *(ptr as *const f64) },
                };
                let offset =
                    unsafe { (ptr as *const u8).offset_from(self.model.as_ptr() as *const u8) };
                return Some((value, offset));
            }
        }
        None
    }

    /// Which node each of this device's nodes was collapsed into, as the model
    /// decided during `setup_instance`.
    ///
    /// `repr[i] == None` means node `i` is shorted to the global reference and
    /// needs no unknown at all; otherwise `repr[i]` is the node whose MNA row
    /// node `i` shares.
    ///
    /// # Why this exists
    ///
    /// A compact model declares more internal nodes than it always uses: BSIM4
    /// has `di`/`si` for the source/drain series resistance, `gi`/`gm` for the
    /// gate network, `dbulk`/`sbulk` for the body diodes. When the parameter that
    /// would give one of them a finite resistance is zero — the default for every
    /// one of them — the model asks the simulator to *collapse* that node into
    /// its neighbour, and then contributes nothing across the short.
    ///
    /// A simulator that ignores the request leaves those nodes with rows nothing
    /// stamps. `d` is then joined to the channel only through `di`, which is
    /// floating, so the drain current is whatever `gmin` allows and the
    /// transistor conducts nothing whatever its bias (#66). Every foundry model
    /// does this; the fixtures in this tree happen not to, which is why the
    /// descriptor's `collapsible`/`collapsed_offset` pair sat unread.
    ///
    /// The union is biased toward the lowest node index, so a terminal always
    /// wins over an internal node — terminals come first in OSDI's numbering.
    pub fn collapse_repr(&self) -> Vec<Option<usize>> {
        let desc = self.desc();
        let n = desc.num_nodes as usize;
        let num_pairs = desc.num_collapsible as usize;
        let mut parent: Vec<usize> = (0..n).collect();
        let mut grounded = vec![false; n];

        fn find(parent: &mut [usize], mut i: usize) -> usize {
            while parent[i] != i {
                parent[i] = parent[parent[i]];
                i = parent[i];
            }
            i
        }

        if num_pairs > 0 && !desc.collapsible.is_null() {
            // C `bool`: one byte per pair, written by `setup_instance`.
            let flags = unsafe {
                std::slice::from_raw_parts(
                    (self.instance.as_ptr() as *const u8).add(desc.collapsed_offset as usize),
                    num_pairs,
                )
            };
            let pairs = unsafe { std::slice::from_raw_parts(desc.collapsible, num_pairs) };
            for (pair, &flag) in pairs.iter().zip(flags) {
                if flag == 0 {
                    continue;
                }
                let a = pair.node_1 as usize;
                if a >= n {
                    continue;
                }
                // `node_2 == UINT32_MAX` is OpenVAF's spelling of "collapse this
                // one to ground" — an implicit equation the model retired.
                if pair.node_2 == u32::MAX || pair.node_2 as usize >= n {
                    let ra = find(&mut parent, a);
                    grounded[ra] = true;
                    continue;
                }
                let b = pair.node_2 as usize;
                let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                if ra != rb {
                    let (keep, drop) = if ra < rb { (ra, rb) } else { (rb, ra) };
                    parent[drop] = keep;
                    grounded[keep] |= grounded[drop];
                }
            }
        }

        (0..n)
            .map(|i| {
                let r = find(&mut parent, i);
                (!grounded[r]).then_some(r)
            })
            .collect()
    }

    /// Rebuild `mna_nodes` from the terminals, the allocated internal rows, and
    /// the model's collapse decisions.
    ///
    /// The dense layout is what `push_device` assumes: terminal `i` on the row
    /// the netlist gave it, internal node `i` on `extra_first + i - num_terminals`.
    /// Collapsing only *redirects* nodes onto a row that layout already holds, so
    /// a collapsed node's own row survives in the matrix with nothing stamped
    /// into it and is pinned by `stamp_gmin`'s empty-row floor. That wastes a row
    /// per collapsed node rather than renumbering, which would put this function
    /// and `push_device`'s thermal-row arithmetic in disagreement — two places
    /// interpreting one layout. Worth revisiting only if a netlist large enough
    /// for the row count to matter turns up.
    fn resolve_mna_nodes(&mut self) {
        let num_nodes = self.desc().num_nodes as usize;
        let num_terminals = self.desc().num_terminals as usize;
        let extra_first = self.extra_first;
        let dense: Vec<NodeId> = (0..num_nodes)
            .map(|i| match i.checked_sub(num_terminals) {
                None => self.terminals.get(i).copied().flatten(),
                Some(off) => extra_first.map(|f| f + off),
            })
            .collect();
        let repr = self.collapse_repr();
        self.mna_nodes = repr.iter().map(|&r| r.and_then(|r| dense[r])).collect();
    }

    /// Write `mna_nodes` into the instance's `node_mapping` array.
    ///
    /// The model reads `node_mapping[i]` from `inst + node_mapping_offset` to
    /// find which solution-vector index its `i`-th node lives on. `UINT32_MAX`
    /// is the sentinel for ground.
    fn write_node_mapping(&mut self) {
        let off = self.desc().node_mapping_offset as usize;
        let num_nodes = self.desc().num_nodes as usize;
        let ptr = unsafe { (self.instance.as_mut_ptr() as *mut u8).add(off) as *mut u32 };
        for i in 0..num_nodes {
            let node = self.mna_nodes.get(i).copied().flatten();
            unsafe {
                *ptr.add(i) = node.map(|n| n as u32).unwrap_or(u32::MAX);
            }
        }
    }

    /// Every operating-point variable the model exposes, in descriptor order.
    ///
    /// These are the model's own account of what it computed — `vth`, `ids`,
    /// `gm`, `weff` and the rest — and on a platform where the formatted-output
    /// tasks have to be stripped (see `portability`) they are the *only* channel
    /// through which a compact model can say what it thinks it is doing.
    pub fn opvar_names(&self) -> Vec<String> {
        let desc = self.desc();
        let n_params = desc.num_params as usize;
        let n_opvars = desc.num_opvars as usize;
        if n_opvars == 0 || desc.param_opvar.is_null() {
            return Vec::new();
        }
        let all = unsafe { std::slice::from_raw_parts(desc.param_opvar, n_params + n_opvars) };
        all[n_params..]
            .iter()
            .map(|p| {
                if p.name.is_null() {
                    return String::new();
                }
                let first = unsafe { *p.name };
                if first.is_null() {
                    return String::new();
                }
                unsafe { CStr::from_ptr(first) }
                    .to_string_lossy()
                    .into_owned()
            })
            .collect()
    }

    /// Read one operating-point variable by name, as of the last `eval`.
    ///
    /// Opvar ids are the plain absolute `param_opvar` index, continuing past the
    /// parameters — the same numbering `set_real_param` uses. The storage is in
    /// the *instance*, so both pointers have to go in.
    pub fn read_opvar(&self, name: &str) -> Option<f64> {
        let desc = self.desc();
        let access_fn = desc.access?;
        let n_params = desc.num_params as usize;
        let n_opvars = desc.num_opvars as usize;
        if n_opvars == 0 || desc.param_opvar.is_null() {
            return None;
        }
        let all = unsafe { std::slice::from_raw_parts(desc.param_opvar, n_params + n_opvars) };
        for (i, p) in all.iter().enumerate().skip(n_params) {
            if !osdi_param_name_matches(p, name) {
                continue;
            }
            let ptr = unsafe {
                access_fn(
                    self.inst_ptr(),
                    self.model_ptr(),
                    i as u32,
                    crate::ffi::ACCESS_FLAG_READ,
                )
            };
            if ptr.is_null() {
                return None;
            }
            return Some(match p.flags & crate::ffi::PARA_TY_MASK {
                crate::ffi::PARA_TY_INT => f64::from(unsafe { *(ptr as *const i32) }),
                crate::ffi::PARA_TY_STR => return None,
                _ => unsafe { *(ptr as *const f64) },
            });
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

/// The module's own name, as a C string, for use as the OSDI handle.
fn descriptor_name(descriptor: *const OsdiDescriptor) -> CString {
    let raw = unsafe { (*descriptor).name };
    let name = if raw.is_null() {
        "osdi_model".to_string()
    } else {
        unsafe { CStr::from_ptr(raw).to_string_lossy().into_owned() }
    };
    CString::new(name).unwrap_or_else(|_| CString::new("osdi_model").unwrap())
}

/// `INIT_ERR_OUT_OF_BOUNDS` from `osdi_0_4.h`.
const INIT_ERR_OUT_OF_BOUNDS: u32 = 1;

impl OsdiDevice {
    /// Name of the parameter an init error points at, when the descriptor has one.
    fn param_name(&self, id: u32) -> String {
        let desc = self.desc();
        if id >= desc.num_params {
            return format!("#{id}");
        }
        unsafe {
            let entry = desc.param_opvar.add(id as usize);
            if entry.is_null() || (*entry).name.is_null() {
                return format!("#{id}");
            }
            let first = *(*entry).name;
            if first.is_null() {
                return format!("#{id}");
            }
            CStr::from_ptr(first).to_string_lossy().into_owned()
        }
    }

    /// Surface whatever a model reported from `setup_model` / `setup_instance`.
    ///
    /// These used to be requested with a null array and never read, so a model
    /// saying "this parameter is out of range" was discarded — the model knew,
    /// and the run continued on a value it had rejected.
    fn report_init(&self, res: &OsdiInitInfo, phase: &str) {
        if res.num_errors == 0 || res.errors.is_null() {
            return;
        }
        let errors = unsafe { std::slice::from_raw_parts(res.errors, res.num_errors as usize) };
        for e in errors {
            match e.code {
                INIT_ERR_OUT_OF_BOUNDS => {
                    let id = unsafe { e.payload.parameter_id };
                    warn_user!(
                        "model '{}' rejected parameter '{}' during {phase}: value out of \
                         the range the model declares, and the model's own default is \
                         being used instead",
                        self.handle.to_string_lossy(),
                        self.param_name(id)
                    );
                }
                other => warn_user!(
                    "model '{}' reported init error code {other} during {phase}",
                    self.handle.to_string_lossy()
                ),
            }
        }
    }
}

impl Device for OsdiDevice {
    fn num_terminals(&self) -> usize {
        self.desc().num_terminals as usize
    }

    fn setup_model(&mut self, ctx: &SimContext) {
        // The table the model's `$simparam` calls read. Held on the device so
        // `call_setup_instance` and `refresh_instance` pass the same one, and so
        // `eval` does not rebuild it per iteration.
        self.sim_paras = Some(SimParaTable::new(ctx));
        self.call_setup_model();
    }

    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        self.terminals = terminals.to_vec();
        self.temperature = ctx.temperature;
        self.setup_done = true;
        self.write_state_indices();
        self.call_setup_instance();
        // Only now does the instance hold the model's collapse decisions, so the
        // mapping is derived after the call rather than before it.
        self.resolve_mna_nodes();
        self.write_node_mapping();
    }

    /// Nodes the model declared `thermal`, read straight off the descriptor.
    ///
    /// OSDI carries a node's discipline through as its `units` string, so a
    /// model that writes `thermal h;` arrives here as `units == "K"` with no
    /// help from the deck and nothing to keep in sync. Both halves of the node
    /// list are covered by the same loop: a thermal *port* other devices share,
    /// and the self-heating internal node most electro-thermal models keep to
    /// themselves.
    fn thermal_nodes(&self) -> Vec<usize> {
        let desc = self.desc();
        (0..desc.num_nodes as usize)
            .filter(|&i| {
                let node = unsafe { &*desc.nodes.add(i) };
                !node.is_flow
                    && !node.units.is_null()
                    && unsafe { CStr::from_ptr(node.units) }.to_bytes() == b"K"
            })
            .collect()
    }

    fn num_extra_nodes(&self) -> usize {
        let desc = self.desc();
        (desc.num_nodes as usize).saturating_sub(desc.num_terminals as usize)
    }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        self.extra_first = Some(first_idx);
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
            // The real table whenever `setup_model` has run, which is every
            // path that reaches an eval. `null_sim_paras` stays the fallback so a
            // device built without setup cannot pass a dangling pointer.
            let mut table = self.sim_paras.take();
            let paras = table.as_mut().map_or_else(null_sim_paras, |t| t.paras());
            let mut info = OsdiSimInfo {
                paras,
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
            let ret = unsafe {
                f(
                    std::ptr::null_mut(),
                    self.inst_ptr(),
                    self.model_ptr(),
                    &mut info,
                )
            };
            self.sim_paras = table;
            // `EVAL_RET_FLAG_LIM` is routine — it says `$limit` moved a voltage,
            // which is what limiting is for. `FATAL` is not.
            self.eval_fatal = ret & crate::ffi::EVAL_RET_FLAG_FATAL != 0;
            // What this eval limited to becomes what the next one limits
            // against.
            std::mem::swap(&mut self.lim_state_prev, &mut self.lim_state_next);
            self.lim_initialised = true;
        }
    }

    /// The timestep the model asked for through `$bound_step`, if any.
    ///
    /// OpenVAF writes the bound into the instance at `bound_step_offset`, as an
    /// `f64`, and initialises it to infinity — so "no request" and "a request"
    /// are the same read with different values, and a finite positive value is
    /// the only thing worth reporting.
    ///
    /// A model asks for this when it knows the solution is about to move faster
    /// than the current step can follow. Ignoring it left the LTE controller as
    /// the only thing watching, and LTE only sees the error after the step that
    /// caused it.
    fn requested_max_timestep(&self) -> Option<f64> {
        let off = self.desc().bound_step_offset;
        // `u32::MAX` is how a descriptor spells "this model never calls
        // `$bound_step`" — measured across every fixture in this tree, which all
        // report 0xffffffff, against 104 for one that does call it. Reading it as
        // a byte offset regardless is a misaligned dereference, and that is what
        // the first version of this did.
        //
        // The alignment and range checks are belt and braces: an `f64` read from
        // a raw offset supplied by a shared library is worth two comparisons.
        let size = self.desc().instance_size as usize;
        let off = off as usize;
        if off == u32::MAX as usize
            || !self.setup_done
            || !off.is_multiple_of(std::mem::align_of::<f64>())
            || off + std::mem::size_of::<f64>() > size
        {
            return None;
        }
        let dt = unsafe { *((self.instance.as_ptr() as *const u8).add(off) as *const f64) };
        // OpenVAF seeds the slot with infinity, so "no request this eval" and "a
        // request" are the same read with different values.
        (dt.is_finite() && dt > 0.0).then_some(dt)
    }

    /// `EVAL_RET_FLAG_FATAL` from the last eval — see [`Device::eval_status`].
    fn eval_status(&self) -> Result<(), String> {
        if self.eval_fatal {
            return Err(format!(
                "model '{}' returned EVAL_RET_FLAG_FATAL: it is telling this \
                 simulator that the bias point it was handed is not one it can \
                 evaluate. Common causes are a parameter outside the range the \
                 model checks, or a node voltage far enough outside normal \
                 operation that an internal expression went non-finite.",
                self.handle.to_string_lossy()
            ));
        }
        Ok(())
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

        // `write_jacobian_array_resist` writes a *packed* array: it walks every
        // entry in `jacobian_entries` order and stores a value only for the ones
        // that have a resistive part. So `jac_buf[k]` belongs to the k-th
        // resistive entry, which is `jacobian_entries[k]` only when no
        // reactive-only entry comes before it.
        //
        // Reading it as a prefix — `jacobian_entries[0..n_resist]` — silently
        // pairs conductances with the wrong node pair on any model that has a
        // capacitance the resistive Jacobian does not also touch. The linearised
        // system then disagrees with the residual it was built from, and Newton
        // converges to a point where the residual is *not* zero: BSIM4 answered
        // 4.05 mA where the model itself computes 3.43 mA, with no diagnostic
        // (#66). The reactive walk in `for_each_react_entry` already keys off the
        // entry's own flag; this now does the same.
        let entries = unsafe {
            std::slice::from_raw_parts(desc.jacobian_entries, desc.num_jacobian_entries as usize)
        };

        let resistive = entries
            .iter()
            .filter(|e| e.flags & crate::ffi::JACOBIAN_ENTRY_RESIST != 0);
        for (i, entry) in resistive.enumerate() {
            let Some(&value) = jac_buf.get(i) else {
                // The descriptor's own two counts disagree. Stamping past the
                // end would be reading the model's memory at random.
                warn_user!(
                    "{}: the model reports {n_resist} resistive Jacobian entries but at                      least {} carry the resistive flag; the extra ones are not stamped",
                    self.handle.to_string_lossy(),
                    i + 1
                );
                break;
            };
            let osdi_r = entry.nodes.node_1 as usize;
            let osdi_c = entry.nodes.node_2 as usize;
            // Map OSDI node index → MNA matrix row/col (skip ground). Two nodes
            // the model collapsed share a row, so their entries accumulate.
            if let (Some(mr), Some(mc)) = (
                self.mna_nodes.get(osdi_r).copied().flatten(),
                self.mna_nodes.get(osdi_c).copied().flatten(),
            ) {
                mat.a[mr][mc] += value;
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
    fn load_reactive_jacobian(&self, c_mat: &mut [fairchild_core::mna::SparseRow]) {
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

    fn lambda_terminals(&self) -> Vec<usize> {
        self.bundle_lambda
            .as_ref()
            .map(|b| b.terminals.clone())
            .unwrap_or_default()
    }

    fn lambda_routing(&self) -> Vec<(usize, usize)> {
        self.bundle_lambda
            .as_ref()
            .map(|b| b.routing.clone())
            .unwrap_or_default()
    }

    /// Fill the generated `wl_<k>` parameters from the resolved wavelengths.
    ///
    /// `LAMBDA(p, k)` expands to `wl_<k>` whatever `p` is — a bundle model has
    /// one channel grid, so slot `k` *is* a wavelength — which is why one value
    /// per slot is the whole answer. It is read off the first bundle port's slot
    /// `k`; after resolution every bundle port of this instance carries the same
    /// label on that slot, because `lambda_routing` says so.
    ///
    /// A deck may still write `wl_0=…` on the instance, and it still takes
    /// effect for a channel no source reaches. Where a source *does* reach it
    /// and disagrees, resolution wins and says so: two answers for one
    /// wavelength is a deck bug, and picking one silently is how a model ends up
    /// evaluating a passband at a colour that is nowhere in the circuit.
    fn set_resolved_lambda(&mut self, per_terminal: &[f64]) {
        let Some(b) = self.bundle_lambda.clone() else {
            return;
        };
        let given = std::mem::take(&mut self.wl_given);
        let lam = b.wpc - 1;
        for k in 0..b.channels {
            // The first bundle's offset is `terminals[0] - lam`; slot k sits
            // `wpc·k` further on.
            let Some(&first) = b.terminals.first() else {
                return;
            };
            let Some(&resolved) = per_terminal.get(first - lam + b.wpc * k + lam) else {
                continue;
            };
            if let Some(&(_, prev)) = given.iter().find(|&&(i, _)| i == k) {
                if (prev - resolved).abs() > 1e-18 {
                    fairchild_core::warn_user!(
                        "{}: instance sets wl_{k}={prev:e} m but the deck's sources put \
                         {resolved:e} m on that channel; the resolved wavelength wins",
                        self.handle.to_string_lossy()
                    );
                }
            }
            self.set_real_param(&format!("wl_{k}"), resolved);
        }
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        let desc = self.desc();
        let access_fn = match desc.access {
            Some(f) => f,
            None => return false,
        };
        // `param_opvar` layout: [inst_params(0..n_inst) | model_params(n_inst..) | opvars],
        // and the `access()` id is that **plain absolute index** — no
        // `PARA_KIND_*` bits. The generated function switches on the id and
        // decides *where* to look from `ACCESS_FLAG_INSTANCE`: set means the
        // instance struct, clear means the model's copy.
        //
        // Getting that wrong made every instance parameter unsettable: an id of
        // `PARA_KIND_INST | i` (0x4000_0000 | i) matched no case, so `access`
        // returned null, the value was dropped, and the caller warned "unknown
        // parameter" about one the model does declare — a MOSFET's `W` and `L`, a
        // `$mfactor`. Model parameters escaped because `PARA_KIND_MODEL` is 0.
        // The mock had no parameters at all, which is why this stood.
        let n_total = desc.num_params as usize;
        let n_inst = desc.num_instance_params as usize;
        if n_total == 0 || desc.param_opvar.is_null() {
            return false;
        }
        let params = unsafe { std::slice::from_raw_parts(desc.param_opvar, n_total) };

        for (i, param) in params.iter().enumerate() {
            if !osdi_param_name_matches(param, name) {
                continue;
            }
            // An instance parameter is written into instance memory; a model
            // parameter into the model's. Either way the write sets the model's
            // "given" flag, so re-running `setup_instance` keeps the new value
            // instead of overwriting it with the default — and re-running is
            // required, because OpenVAF caches parameter-derived quantities in
            // the instance during setup and `eval` reads only those.
            let (inst, flags) = if i < n_inst {
                (
                    self.inst_ptr(),
                    crate::ffi::ACCESS_FLAG_SET | crate::ffi::ACCESS_FLAG_INSTANCE,
                )
            } else {
                (std::ptr::null_mut(), crate::ffi::ACCESS_FLAG_SET)
            };
            let ptr = unsafe { access_fn(inst, self.model_ptr(), i as u32, flags) };
            if !ptr.is_null() {
                // Write the width the model declared, not the width we happen to
                // hold. `integer type = 1` stored as f64 reads back as 0 — a
                // BSIM4 that is neither NMOS nor PMOS and conducts nothing.
                match param.flags & crate::ffi::PARA_TY_MASK {
                    crate::ffi::PARA_TY_INT => {
                        let rounded = value.round();
                        if (value - rounded).abs() > 1e-9 {
                            fairchild_core::warn_user!(
                                "{}: parameter '{name}' is an integer in the model; \
                                 {value} was rounded to {rounded}",
                                self.handle.to_string_lossy()
                            );
                        }
                        unsafe { *(ptr as *mut i32) = rounded as i32 }
                    }
                    crate::ffi::PARA_TY_STR => {
                        fairchild_core::warn_user!(
                            "{}: parameter '{name}' is a string in the model and cannot \
                             be set from a numeric deck value; it keeps its default",
                            self.handle.to_string_lossy()
                        );
                        return false;
                    }
                    _ => unsafe { *(ptr as *mut f64) = value },
                }
                if let Some(k) = name
                    .strip_prefix("wl_")
                    .and_then(|d| d.parse::<usize>().ok())
                {
                    self.wl_given.push((k, value));
                }
                if self.setup_done {
                    self.refresh_instance();
                }
                return true;
            }
        }

        false
    }
}

impl OsdiDevice {
    /// Re-derive everything a parameter write invalidated: `setup_model`, then
    /// `setup_instance`, then the node mapping.
    ///
    /// Reached from `bind_extra_nodes`, and from a `set_real_param` on a device
    /// that is already set up — `set_resolved_lambda` writing a bundle's `wl_k`,
    /// or any caller holding a built device. A deck's own parameters no longer
    /// come this way: the registry applies them before setup, which is OSDI's
    /// own order.
    ///
    /// **Both halves, not just the instance.** This used to redo only
    /// `setup_instance`, which happens to be enough for every model in this tree
    /// — OpenVAF derives into instance storage even for a parameter that lives
    /// on the model. It is not enough in general, because `model_data` exists and
    /// a model that caches a card-level derivation there would keep its default
    /// in silence. Redoing both is one extra call and removes the exception.
    ///
    /// Not the #66 fix, though it was written while looking for it: reverting
    /// this alone leaves BSIM4 answering correctly.
    fn refresh_instance(&mut self) {
        self.call_setup_model();
        self.write_state_indices();
        self.call_setup_instance();
        self.resolve_mna_nodes();
        self.write_node_mapping();
    }

    /// The `setup_model` call itself, with nothing around it.
    fn call_setup_model(&mut self) {
        let setup_fn = self.desc().setup_model;
        if let Some(f) = setup_fn {
            // The device's own `$simparam` table. Taken and put back so the
            // storage it points into outlives the call without a second
            // mutable borrow of `self`.
            let mut table = self.sim_paras.take();
            let mut paras = table.as_mut().map_or_else(null_sim_paras, |t| t.paras());
            let mut res = OsdiInitInfo {
                flags: 0,
                num_errors: 0,
                errors: std::ptr::null_mut(),
            };
            unsafe {
                f(
                    self.handle.as_ptr() as *mut c_void,
                    self.model_ptr(),
                    &mut paras,
                    &mut res,
                );
            }
            self.sim_paras = table;
            self.report_init(&res, "model setup");
        }
    }

    /// The `setup_instance` call itself, with nothing around it.
    ///
    /// Split out so `setup_instance` and `refresh_instance` cannot drift on what
    /// they pass the model — the temperature argument in particular, which used
    /// to be `SimContext::default()` on the refresh path.
    fn call_setup_instance(&mut self) {
        let setup_fn = self.desc().setup_instance;
        // OSDI's `num_terminals` argument is the count of *external* nodes;
        // internal flow-branch nodes are implicit (the rest of num_nodes).
        let num_terminals = self.desc().num_terminals;
        let temperature = self.temperature;
        if let Some(f) = setup_fn {
            // The device's own `$simparam` table. Taken and put back so the
            // storage it points into outlives the call without a second
            // mutable borrow of `self`.
            let mut table = self.sim_paras.take();
            let mut paras = table.as_mut().map_or_else(null_sim_paras, |t| t.paras());
            let mut res = OsdiInitInfo {
                flags: 0,
                num_errors: 0,
                errors: std::ptr::null_mut(),
            };
            unsafe {
                f(
                    self.handle.as_ptr() as *mut c_void,
                    self.inst_ptr(),
                    self.model_ptr(),
                    temperature,
                    num_terminals,
                    &mut paras,
                    &mut res,
                );
            }
            self.sim_paras = table;
            self.report_init(&res, "instance setup");
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
    // A SPICE deck spells the instance multiplier `m=`; Verilog-A calls it
    // `$mfactor`, and OpenVAF gives every module one implicitly. They are the
    // same quantity — the model scales its own contributions by it — so the deck's
    // spelling has to reach it. Without this, `m=` on a compiled device warned
    // "unknown parameter" and the factor was lost, which is a wrong answer
    // exactly the size of the factor.
    let target = if target.eq_ignore_ascii_case("m") {
        "$mfactor"
    } else {
        target
    };
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

/// An empty simulator-parameter table.
///
/// `names` and `names_str` are **NUL-terminated lists**, so "no parameters" is a
/// pointer to a single null entry — not a null pointer. A null list is what a
/// model walks off the end of while looking up `$simparam`, and it crashes inside
/// the message it was trying to format about the lookup failing. OpenVAF's own
/// runtime spells this `names: &mut ptr::null_mut()`.
fn empty_name_list() -> *mut *mut std::os::raw::c_char {
    // One permanently-leaked null entry, shared by every call. Stored as a usize
    // because a raw pointer is not `Sync` and this genuinely is: it is read-only
    // and lives for the process.
    static LIST: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    let addr = *LIST.get_or_init(|| {
        let slot: &'static mut *mut std::os::raw::c_char =
            Box::leak(Box::new(std::ptr::null_mut()));
        slot as *mut *mut std::os::raw::c_char as usize
    });
    addr as *mut *mut std::os::raw::c_char
}

fn null_sim_paras() -> OsdiSimParas {
    OsdiSimParas {
        names: empty_name_list(),
        vals: std::ptr::null_mut(),
        names_str: empty_name_list(),
        vals_str: std::ptr::null_mut(),
    }
}

/// The `$simparam` table a model sees, with the storage it points into.
///
/// `OsdiSimParas` is four raw pointers, so the `CString`s and the pointer arrays
/// have to outlive the call. Keeping them in one struct beside the paras makes
/// that a borrow-checker fact rather than a comment.
///
/// # What is in it, and what is deliberately not
///
/// `gmin` is the one that matters: a foundry model asks for it, and until this
/// existed it got whatever default the model's own call site wrote — so
/// `.options gmin=` reached the native models and not the compiled ones. The same
/// fault the `const GMIN` in each native model had.
///
/// A name this table does not carry falls back to the model's own default, which
/// is what every name did before, so adding names is monotone and leaving one out
/// cannot break a model. Left out on purpose:
///
/// * `iteration` / `sourceScaleFactor` — these change per Newton iteration, and
///   building a table per eval to carry them would cost an allocation on the
///   hottest path for something no model in the OpenVAF-Reloaded suite reads.
/// * `imax` / `imelt` — current limits fairchild does not enforce, so answering
///   would be a claim it does not keep.
struct SimParaTable {
    _names: Vec<std::ffi::CString>,
    name_ptrs: Vec<*mut std::os::raw::c_char>,
    vals: Vec<f64>,
}

impl SimParaTable {
    fn new(ctx: &SimContext) -> Self {
        // `scale` and `shrink` are geometry multipliers a PDK may read; fairchild
        // applies neither, and 1.0 is what "neither" means rather than a claim.
        let entries: &[(&str, f64)] = &[
            ("gmin", ctx.gmin),
            ("scale", 1.0),
            ("shrink", 0.0),
            ("simulatorSubversion", 0.0),
        ];
        let names: Vec<std::ffi::CString> = entries
            .iter()
            .map(|(n, _)| std::ffi::CString::new(*n).expect("no interior NUL"))
            .collect();
        let mut name_ptrs: Vec<*mut std::os::raw::c_char> = names
            .iter()
            .map(|c| c.as_ptr() as *mut std::os::raw::c_char)
            .collect();
        // NUL-terminated: a model walks this list until it hits a null entry.
        name_ptrs.push(std::ptr::null_mut());
        let vals = entries.iter().map(|(_, v)| *v).collect();
        Self {
            _names: names,
            name_ptrs,
            vals,
        }
    }

    /// Borrows this table's storage. The result must not outlive `self`.
    fn paras(&mut self) -> OsdiSimParas {
        OsdiSimParas {
            names: self.name_ptrs.as_mut_ptr(),
            vals: self.vals.as_mut_ptr(),
            // No string-valued parameters. Terminated, not null — the same
            // reason `null_sim_paras` uses `empty_name_list`.
            names_str: empty_name_list(),
            vals_str: std::ptr::null_mut(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// "No simulator parameters" is a pointer to a null entry, not a null pointer.
    ///
    /// `names` and `names_str` are NUL-terminated lists, and a model walks them
    /// looking up `$simparam`. Handing it a null list is what it walks off, and it
    /// crashes inside the diagnostic it was trying to format about the failure.
    /// OpenVAF's own runtime spells the empty case `&mut ptr::null_mut()`.
    #[test]
    fn an_empty_sim_para_list_is_terminated_not_null() {
        let paras = null_sim_paras();
        assert!(!paras.names.is_null(), "the list itself must exist");
        assert!(!paras.names_str.is_null());
        unsafe {
            assert!(
                (*paras.names).is_null(),
                "its first entry must terminate it"
            );
            assert!((*paras.names_str).is_null());
        }
        // Values may be null: a model reads a value only after finding its name.
        assert!(paras.vals.is_null());
    }

    /// The list is shared and stable, so a device built later sees the same one.
    #[test]
    fn the_empty_list_is_stable_across_calls() {
        assert_eq!(empty_name_list(), empty_name_list());
    }
}
