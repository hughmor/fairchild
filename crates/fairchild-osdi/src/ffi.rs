//! `#[repr(C)]` types matching the OSDI v0.4 ABI.
//!
//! Source: OpenVAF-Reloaded `openvaf/osdi/header/osdi_0_4.h` (mob branch, GitHub).
//! All struct layouts verified against the C header field-by-field.
//! Expected size on 64-bit: see compile-time asserts at the bottom of this file.

use std::ffi::c_char;
use std::os::raw::c_void;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

pub const OSDI_VERSION_MAJOR_CURR: u32 = 0;
pub const OSDI_VERSION_MINOR_CURR: u32 = 4;

// CALC_* flags passed to eval()
pub const CALC_RESIST_RESIDUAL: u32 = 1;
pub const CALC_REACT_RESIDUAL: u32 = 2;
pub const CALC_RESIST_JACOBIAN: u32 = 4;
pub const CALC_REACT_JACOBIAN: u32 = 8;
pub const CALC_NOISE: u32 = 16;
pub const CALC_OP: u32 = 32;
pub const CALC_RESIST_LIM_RHS: u32 = 64;
pub const CALC_REACT_LIM_RHS: u32 = 128;
/// Apply `$limit()` limiting this evaluation, against `prev_state`.
pub const ENABLE_LIM: u32 = 256;
/// First evaluation of a solve: seed the limiting state instead of limiting
/// against a `prev_state` that means nothing yet.
pub const INIT_LIM: u32 = 512;
pub const ANALYSIS_DC: u32 = 2048;
pub const ANALYSIS_AC: u32 = 4096;
pub const ANALYSIS_TRAN: u32 = 8192;

pub const EVAL_RET_FLAG_LIM: u32 = 1;
pub const EVAL_RET_FLAG_FATAL: u32 = 2;

// `OsdiJacobianEntry::flags` bits from `osdi_0_4.h`.
/// The entry has a resistive (∂f/∂x) part. `write_jacobian_array_resist` writes
/// one value per entry carrying this bit and **skips the rest**, so the array it
/// fills is packed: the k-th value belongs to the k-th *resistive* entry, not to
/// `jacobian_entries[k]`.
pub const JACOBIAN_ENTRY_RESIST: u32 = 4;
/// The entry has a reactive (∂q/∂x) part; `react_ptr_off != u32::MAX` says the
/// same thing, and the reactive walk keys off that.
pub const JACOBIAN_ENTRY_REACT: u32 = 8;

// Parameter kind/type flags
pub const PARA_TY_REAL: u32 = 0;
pub const PARA_TY_INT: u32 = 1;
pub const PARA_TY_STR: u32 = 2;
/// Low bits of `OsdiParamOpvar::flags` holding the `PARA_TY_*` storage type.
///
/// The width of the slot `access()` hands back follows this, not the caller's
/// wishes: an `integer` parameter is an `i32` in the model struct, so writing a
/// `f64` over it stores the low half of an IEEE bit pattern (1.0 becomes 0) and
/// tramples the next four bytes as well.
pub const PARA_TY_MASK: u32 = 0x3;
pub const PARA_KIND_MODEL: u32 = 0 << 30;
pub const PARA_KIND_INST: u32 = 1 << 30;
pub const PARA_KIND_OPVAR: u32 = 2 << 30;

// Access flags for the `access` function pointer
pub const ACCESS_FLAG_READ: u32 = 0;
pub const ACCESS_FLAG_SET: u32 = 1;
/// Take the pointer out of the *instance* struct rather than the model's.
///
/// An instance parameter lives in instance memory, and `access()` returns null
/// for one unless this flag says where to look. Without it every instance
/// parameter — a MOSFET's `W` and `L`, a `$mfactor` — silently fails to apply.
pub const ACCESS_FLAG_INSTANCE: u32 = 4;

// ---------------------------------------------------------------------------
// Simple structs
// ---------------------------------------------------------------------------

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct OsdiNodePair {
    pub node_1: u32,
    pub node_2: u32,
}

/// A terminal or internal node of the device.
#[repr(C)]
pub struct OsdiNode {
    pub name: *mut c_char,
    pub units: *mut c_char,
    pub residual_units: *mut c_char,
    pub resist_residual_off: u32,
    pub react_residual_off: u32,
    pub resist_limit_rhs_off: u32,
    pub react_limit_rhs_off: u32,
    /// true → flow quantity (current), false → potential (voltage)
    pub is_flow: bool,
    // 7 bytes implicit C padding to align struct to 8 bytes
}

#[repr(C)]
pub struct OsdiJacobianEntry {
    pub nodes: OsdiNodePair,
    pub react_ptr_off: u32,
    pub flags: u32,
}

#[repr(C)]
pub struct OsdiParamOpvar {
    /// Null-terminated name aliases array (length = num_alias).
    pub name: *mut *mut c_char,
    pub num_alias: u32,
    // 4 bytes C padding before next pointer
    pub description: *mut c_char,
    pub units: *mut c_char,
    pub flags: u32,
    pub len: u32,
}

#[repr(C)]
pub struct OsdiNoiseSource {
    pub name: *mut c_char,
    pub nodes: OsdiNodePair,
}

/// One entry of a library's `OSDI_LIM_TABLE`.
///
/// The library exports the table with every `func_ptr` **null**; the simulator
/// is expected to walk it and install its own implementation of each named
/// limiting function before any `eval`. OpenVAF emits the call
/// unconditionally, so leaving a `func_ptr` null means a jump to address 0 the
/// first time the model evaluates.
#[repr(C)]
pub struct OsdiLimFunction {
    pub name: *mut c_char,
    /// Extra arguments beyond `(init, check, vnew, vold)`. `pnjlim` reports 2:
    /// `vt` and `vcrit`.
    pub num_args: u32,
    // 4 bytes C padding before next pointer
    pub func_ptr: *mut c_void,
}

/// `pnjlim` as OSDI expects it — see `OsdiLimFunction`. Returns the limited
/// voltage and sets `*check` when it intervened, which tells the simulator the
/// iteration is not yet converged.
pub type FnPnjlim = unsafe extern "C" fn(
    init: bool,
    check: *mut bool,
    vnew: f64,
    vold: f64,
    vt: f64,
    vcrit: f64,
) -> f64;

#[repr(C)]
pub struct OsdiSimParas {
    pub names: *mut *mut c_char,
    pub vals: *mut f64,
    pub names_str: *mut *mut c_char,
    pub vals_str: *mut *mut c_char,
}

#[repr(C)]
pub struct OsdiSimInfo {
    pub paras: OsdiSimParas,
    pub abstime: f64,
    pub prev_solve: *mut f64,
    pub prev_state: *mut f64,
    pub next_state: *mut f64,
    pub flags: u32,
    // 4 bytes C padding to align struct to 8 bytes
}

#[repr(C)]
pub union OsdiInitErrorPayload {
    pub parameter_id: u32,
}

#[repr(C)]
pub struct OsdiInitError {
    pub code: u32,
    pub payload: OsdiInitErrorPayload,
}

#[repr(C)]
pub struct OsdiInitInfo {
    pub flags: u32,
    pub num_errors: u32,
    pub errors: *mut OsdiInitError,
}

// ---------------------------------------------------------------------------
// Function pointer type aliases
// ---------------------------------------------------------------------------

pub type FnAccess =
    unsafe extern "C" fn(inst: *mut c_void, model: *mut c_void, id: u32, flags: u32) -> *mut c_void;

pub type FnSetupModel = unsafe extern "C" fn(
    handle: *mut c_void,
    model: *mut c_void,
    sim_params: *mut OsdiSimParas,
    res: *mut OsdiInitInfo,
);

pub type FnSetupInstance = unsafe extern "C" fn(
    handle: *mut c_void,
    inst: *mut c_void,
    model: *mut c_void,
    temperature: f64,
    num_terminals: u32,
    sim_params: *mut OsdiSimParas,
    res: *mut OsdiInitInfo,
);

pub type FnEval = unsafe extern "C" fn(
    handle: *mut c_void,
    inst: *mut c_void,
    model: *mut c_void,
    info: *mut OsdiSimInfo,
) -> u32;

pub type FnLoadNoise =
    unsafe extern "C" fn(inst: *mut c_void, model: *mut c_void, freq: f64, noise_dens: *mut f64);

/// Shared signature for load_residual_resist / load_residual_react /
/// load_limit_rhs_resist / load_limit_rhs_react.
pub type FnLoadResidual =
    unsafe extern "C" fn(inst: *mut c_void, model: *mut c_void, dst: *mut f64);

pub type FnLoadSpiceRhsDc = unsafe extern "C" fn(
    inst: *mut c_void,
    model: *mut c_void,
    dst: *mut f64,
    prev_solve: *mut f64,
);

pub type FnLoadSpiceRhsTran = unsafe extern "C" fn(
    inst: *mut c_void,
    model: *mut c_void,
    dst: *mut f64,
    prev_solve: *mut f64,
    alpha: f64,
);

/// Shared signature for load_jacobian_resist.
pub type FnLoadJacobian = unsafe extern "C" fn(inst: *mut c_void, model: *mut c_void);

/// Shared signature for load_jacobian_react / load_jacobian_tran.
pub type FnLoadJacobianAlpha =
    unsafe extern "C" fn(inst: *mut c_void, model: *mut c_void, alpha: f64);

/// Shared signature for given_flag_model / given_flag_instance.
pub type FnGivenFlag = unsafe extern "C" fn(ptr: *mut c_void, id: u32) -> u32;

/// Shared signature for write_jacobian_array_resist / write_jacobian_array_react.
pub type FnWriteJacobianArray =
    unsafe extern "C" fn(inst: *mut c_void, model: *mut c_void, destination: *mut f64);

/// Shared signature for load_jacobian_with_offset_resist / load_jacobian_with_offset_react.
pub type FnLoadJacobianWithOffset =
    unsafe extern "C" fn(inst: *mut c_void, model: *mut c_void, offset: usize);

// ---------------------------------------------------------------------------
// OsdiDescriptor — the central export structure
// ---------------------------------------------------------------------------

/// Top-level descriptor for one compiled Verilog-A module.
///
/// An `.osdi` shared library exports `OSDI_DESCRIPTORS[0..OSDI_NUM_DESCRIPTORS]`,
/// each one of these structs, strided by `OSDI_DESCRIPTOR_SIZE` bytes.
#[repr(C)]
pub struct OsdiDescriptor {
    pub name: *mut c_char,

    pub num_nodes: u32,
    pub num_terminals: u32,
    pub nodes: *mut OsdiNode,

    pub num_jacobian_entries: u32,
    // 4 bytes C padding before next pointer
    pub jacobian_entries: *mut OsdiJacobianEntry,

    pub num_collapsible: u32,
    // 4 bytes C padding before next pointer
    pub collapsible: *mut OsdiNodePair,

    pub collapsed_offset: u32,
    // 4 bytes C padding before next pointer
    pub noise_sources: *mut OsdiNoiseSource,

    pub num_noise_src: u32,
    pub num_params: u32,
    pub num_instance_params: u32,
    pub num_opvars: u32,
    pub param_opvar: *mut OsdiParamOpvar,

    pub node_mapping_offset: u32,
    pub jacobian_ptr_resist_offset: u32,
    pub num_states: u32,
    pub state_idx_off: u32,
    pub bound_step_offset: u32,
    pub instance_size: u32,
    pub model_size: u32,
    // 4 bytes C padding before first function pointer (align 8)
    pub access: Option<FnAccess>,
    pub setup_model: Option<FnSetupModel>,
    pub setup_instance: Option<FnSetupInstance>,
    pub eval: Option<FnEval>,
    pub load_noise: Option<FnLoadNoise>,
    pub load_residual_resist: Option<FnLoadResidual>,
    pub load_residual_react: Option<FnLoadResidual>,
    pub load_limit_rhs_resist: Option<FnLoadResidual>,
    pub load_limit_rhs_react: Option<FnLoadResidual>,
    pub load_spice_rhs_dc: Option<FnLoadSpiceRhsDc>,
    pub load_spice_rhs_tran: Option<FnLoadSpiceRhsTran>,
    pub load_jacobian_resist: Option<FnLoadJacobian>,
    pub load_jacobian_react: Option<FnLoadJacobianAlpha>,
    pub load_jacobian_tran: Option<FnLoadJacobianAlpha>,

    // v0.4 additions (not present in v0.3)
    pub given_flag_model: Option<FnGivenFlag>,
    pub given_flag_instance: Option<FnGivenFlag>,
    pub num_resistive_jacobian_entries: u32,
    pub num_reactive_jacobian_entries: u32,
    pub write_jacobian_array_resist: Option<FnWriteJacobianArray>,
    pub write_jacobian_array_react: Option<FnWriteJacobianArray>,
    pub num_inputs: u32,
    // 4 bytes C padding before next pointer
    pub inputs: *mut OsdiNodePair,
    pub load_jacobian_with_offset_resist: Option<FnLoadJacobianWithOffset>,
    pub load_jacobian_with_offset_react: Option<FnLoadJacobianWithOffset>,
}

// OSDI_DESCRIPTORS is written once by the model's static initialiser and is
// read-only from the simulator's perspective.
unsafe impl Sync for OsdiDescriptor {}
unsafe impl Send for OsdiDescriptor {}

// ---------------------------------------------------------------------------
// Compile-time layout assertions (catches any padding miscounts)
// ---------------------------------------------------------------------------

const _: () = assert!(std::mem::size_of::<OsdiDescriptor>() == 312);
const _: () = assert!(std::mem::size_of::<OsdiNode>() == 48);
const _: () = assert!(std::mem::size_of::<OsdiJacobianEntry>() == 16);
const _: () = assert!(std::mem::size_of::<OsdiParamOpvar>() == 40);
const _: () = assert!(std::mem::size_of::<OsdiNoiseSource>() == 16);
const _: () = assert!(std::mem::size_of::<OsdiLimFunction>() == 24);
const _: () = assert!(std::mem::size_of::<OsdiSimParas>() == 32);
const _: () = assert!(std::mem::size_of::<OsdiSimInfo>() == 72);
