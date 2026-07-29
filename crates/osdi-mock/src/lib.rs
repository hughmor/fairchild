//! OSDI v0.4 mock library for integration testing.
//!
//! Exports one model named "test_conductance" with two terminals: a linear
//! 1 mS conductance in parallel with a linear 1 nF capacitance.  The
//! capacitance is what gives the reactive-Jacobian path (transient companion
//! and the `.ac` / `.noise` susceptance block) any CI coverage at all — every
//! other reactive OSDI model in the tree needs OpenVAF to build.
//!   - setup_model: writes gd = 1e-3 S and c = 1e-9 F into model memory
//!   - setup_instance: no-op (node mapping is written by the simulator)
//!   - eval: records OsdiSimInfo.abstime (what `$abstime` reads) and returns 0
//!   - load_spice_rhs_dc: no-op (Jeq = 0 for a linear element)
//!   - write_jacobian_array_resist: writes [+gd, -gd, -gd, +gd] for the 4 entries
//!   - write_jacobian_array_react:   writes [+c,  -c,  -c,  +c ] for the same 4
//!
//! Instance memory layout (instance_size = 8):
//!   bytes 0-7: node_mapping[2] (two u32 MNA indices) at node_mapping_offset = 0
//!
//! Model memory layout (model_size = 24):
//!   bytes 0-7:   gd (f64, conductance in siemens)
//!   bytes 8-15:  abstime seen by the last eval (f64) — test observability only
//!   bytes 16-23: c (f64, capacitance in farads)
//!
//! Jacobian entries (4, between OSDI nodes 0 and 1; a parallel G and C share
//! the same sparsity, so the same 4 entries carry both):
//!   [0] (0,0) → +gd/+c,  [1] (0,1) → -gd/-c,
//!   [2] (1,0) → -gd/-c,  [3] (1,1) → +gd/+c

use std::ffi::c_char;
use std::os::raw::c_void;

use fairchild_osdi::ffi::{
    OsdiDescriptor, OsdiInitInfo, OsdiJacobianEntry, OsdiLimFunction, OsdiNodePair, OsdiSimInfo,
    OsdiSimParas,
};

// ---------------------------------------------------------------------------
// Model function implementations
// ---------------------------------------------------------------------------

unsafe extern "C" fn setup_model_impl(
    _handle: *mut c_void,
    model: *mut c_void,
    _sim_params: *mut OsdiSimParas,
    _res: *mut OsdiInitInfo,
) {
    *(model as *mut f64) = 1e-3; // gd = 1 mS
    *((model as *mut u8).add(CAP_OFFSET) as *mut f64) = 1e-9; // c = 1 nF
}

unsafe extern "C" fn setup_instance_impl(
    _handle: *mut c_void,
    _inst: *mut c_void,
    _model: *mut c_void,
    _temperature: f64,
    _num_terminals: u32,
    _sim_params: *mut OsdiSimParas,
    _res: *mut OsdiInitInfo,
) {
    // Node mapping is written by the simulator before this call.
    // No instance-local initialisation needed for a linear conductance.
}

/// Byte offset of the abstime slot in model memory (after `gd`).
pub const ABSTIME_OFFSET: usize = 8;
/// Byte offset of the capacitance in model memory.
pub const CAP_OFFSET: usize = 16;
/// The conductance this mock stamps, in siemens.
pub const MOCK_GD: f64 = 1e-3;
/// The capacitance this mock stamps, in farads.
pub const MOCK_C: f64 = 1e-9;

unsafe extern "C" fn eval_impl(
    _handle: *mut c_void,
    _inst: *mut c_void,
    model: *mut c_void,
    info: *mut OsdiSimInfo,
) -> u32 {
    // Record what a Verilog-A `$abstime` would have read, so a test can assert
    // the simulator passes real transient time rather than a hardcoded 0.0.
    // Model memory, not a Rust static: the mock is dlopen'd, so a static in
    // the loaded copy is not the one a test linking the rlib would observe.
    if !info.is_null() && !model.is_null() {
        *((model as *mut u8).add(ABSTIME_OFFSET) as *mut f64) = (*info).abstime;
    }
    0 // success; no voltage-limiting flag set
}

unsafe extern "C" fn load_spice_rhs_dc_impl(
    _inst: *mut c_void,
    _model: *mut c_void,
    _dst: *mut f64,
    _prev_solve: *mut f64,
) {
    // Linear conductance: Jeq = Id - gd*Vd = gd*Vd - gd*Vd = 0.
    // No contribution to the RHS.
}

unsafe extern "C" fn write_jacobian_array_resist_impl(
    _inst: *mut c_void,
    model: *mut c_void,
    destination: *mut f64,
) {
    let gd = *(model as *const f64);
    // 4 entries matching JACOBIAN_ENTRIES order: (0,0), (0,1), (1,0), (1,1).
    *destination.add(0) = gd; // ∂F[0]/∂V[0] = +gd
    *destination.add(1) = -gd; // ∂F[0]/∂V[1] = -gd
    *destination.add(2) = -gd; // ∂F[1]/∂V[0] = -gd
    *destination.add(3) = gd; // ∂F[1]/∂V[1] = +gd
}

unsafe extern "C" fn write_jacobian_array_react_impl(
    _inst: *mut c_void,
    model: *mut c_void,
    destination: *mut f64,
) {
    let c = *((model as *const u8).add(CAP_OFFSET) as *const f64);
    // Same 4 entries, same order as the resistive write.
    *destination.add(0) = c;
    *destination.add(1) = -c;
    *destination.add(2) = -c;
    *destination.add(3) = c;
}

// ---------------------------------------------------------------------------
// Jacobian entry table
// ---------------------------------------------------------------------------

static JACOBIAN_ENTRIES: [OsdiJacobianEntry; 4] = [
    OsdiJacobianEntry {
        nodes: OsdiNodePair {
            node_1: 0,
            node_2: 0,
        },
        react_ptr_off: 0,
        flags: 0,
    },
    OsdiJacobianEntry {
        nodes: OsdiNodePair {
            node_1: 0,
            node_2: 1,
        },
        react_ptr_off: 0,
        flags: 0,
    },
    OsdiJacobianEntry {
        nodes: OsdiNodePair {
            node_1: 1,
            node_2: 0,
        },
        react_ptr_off: 0,
        flags: 0,
    },
    OsdiJacobianEntry {
        nodes: OsdiNodePair {
            node_1: 1,
            node_2: 1,
        },
        react_ptr_off: 0,
        flags: 0,
    },
];

// ---------------------------------------------------------------------------
// Limiting table
// ---------------------------------------------------------------------------

/// A `$limit(..., "pnjlim", ...)` request, exported exactly as OpenVAF does:
/// `func_ptr` starts null and the simulator is expected to install its own
/// implementation before any `eval`. This mock never calls it — the point is
/// that a test can check the loader filled it in, because a null here is a
/// jump to address 0 in a real compiled model.
#[no_mangle]
pub static mut OSDI_LIM_TABLE: [OsdiLimFunction; 1] = [OsdiLimFunction {
    name: c"pnjlim".as_ptr().cast_mut(),
    num_args: 2,
    func_ptr: std::ptr::null_mut(),
}];

#[no_mangle]
pub static OSDI_LIM_TABLE_LEN: u32 = 1;

// ---------------------------------------------------------------------------
// OSDI exports
// ---------------------------------------------------------------------------

#[no_mangle]
pub static OSDI_VERSION_MAJOR: u32 = 0;
#[no_mangle]
pub static OSDI_VERSION_MINOR: u32 = 4;
#[no_mangle]
pub static OSDI_NUM_DESCRIPTORS: u32 = 1;
#[no_mangle]
pub static OSDI_DESCRIPTOR_SIZE: u32 = std::mem::size_of::<OsdiDescriptor>() as u32;

#[no_mangle]
pub static mut osdi_log: Option<unsafe extern "C" fn(*mut c_void, u32, *const c_char)> = None;

#[no_mangle]
pub static OSDI_DESCRIPTORS: [OsdiDescriptor; 1] = [OsdiDescriptor {
    name: c"test_conductance".as_ptr().cast_mut(),

    num_nodes: 2,
    num_terminals: 2,
    nodes: std::ptr::null_mut(), // OsdiNode array not needed for copy-based path

    num_jacobian_entries: 4,
    jacobian_entries: JACOBIAN_ENTRIES.as_ptr().cast_mut(),

    num_collapsible: 0,
    collapsible: std::ptr::null_mut(),

    collapsed_offset: 0,
    noise_sources: std::ptr::null_mut(),

    num_noise_src: 0,
    num_params: 0, // gd is hardcoded; no user-settable parameters
    num_instance_params: 0,
    num_opvars: 0,
    param_opvar: std::ptr::null_mut(),

    node_mapping_offset: 0,        // [u32; 2] at byte 0 of instance buffer
    jacobian_ptr_resist_offset: 8, // unused (copy-based path); must not alias node_mapping
    num_states: 0,
    state_idx_off: 0,
    bound_step_offset: 0,
    instance_size: 8, // 2 × u32 = 8 bytes
    model_size: 24,   // gd + the recorded abstime + c = 3 × f64

    access: None,
    setup_model: Some(setup_model_impl),
    setup_instance: Some(setup_instance_impl),
    eval: Some(eval_impl),
    load_noise: None,
    load_residual_resist: None,
    load_residual_react: None,
    load_limit_rhs_resist: None,
    load_limit_rhs_react: None,
    load_spice_rhs_dc: Some(load_spice_rhs_dc_impl),
    load_spice_rhs_tran: None,
    load_jacobian_resist: None,
    load_jacobian_react: None,
    load_jacobian_tran: None,

    given_flag_model: None,
    given_flag_instance: None,
    num_resistive_jacobian_entries: 4,
    num_reactive_jacobian_entries: 4,
    write_jacobian_array_resist: Some(write_jacobian_array_resist_impl),
    write_jacobian_array_react: Some(write_jacobian_array_react_impl),
    num_inputs: 0,
    inputs: std::ptr::null_mut(),
    load_jacobian_with_offset_resist: None,
    load_jacobian_with_offset_react: None,
}];
