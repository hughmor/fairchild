//! OSDI v0.4 mock library for integration testing.
//!
//! Exports one model named "test_conductance" with two terminals.
//! All function pointers are implemented with a simple linear 1 mS conductance:
//!   - setup_model: writes gd = 1e-3 S into model memory
//!   - setup_instance: no-op (node mapping is written by the simulator)
//!   - eval: no-op (linear element, nothing to linearise)
//!   - load_spice_rhs_dc: no-op (Jeq = 0 for a linear element)
//!   - write_jacobian_array_resist: writes [+gd, -gd, -gd, +gd] for the 4 entries
//!
//! Instance memory layout (instance_size = 8):
//!   bytes 0-7: node_mapping[2] (two u32 MNA indices) at node_mapping_offset = 0
//!
//! Model memory layout (model_size = 8):
//!   bytes 0-7: gd (f64, conductance in siemens)
//!
//! Jacobian entries (4 resistive, representing G between OSDI nodes 0 and 1):
//!   [0] (0,0) → +gd,  [1] (0,1) → -gd,  [2] (1,0) → -gd,  [3] (1,1) → +gd

use std::ffi::c_char;
use std::os::raw::c_void;

use fairchild_osdi::ffi::{
    OsdiDescriptor, OsdiInitInfo, OsdiJacobianEntry, OsdiNodePair, OsdiSimInfo, OsdiSimParas,
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

unsafe extern "C" fn eval_impl(
    _handle: *mut c_void,
    _inst: *mut c_void,
    _model: *mut c_void,
    _info: *mut OsdiSimInfo,
) -> u32 {
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
    *destination.add(0) = gd;   // ∂F[0]/∂V[0] = +gd
    *destination.add(1) = -gd;  // ∂F[0]/∂V[1] = -gd
    *destination.add(2) = -gd;  // ∂F[1]/∂V[0] = -gd
    *destination.add(3) = gd;   // ∂F[1]/∂V[1] = +gd
}

// ---------------------------------------------------------------------------
// Jacobian entry table
// ---------------------------------------------------------------------------

static JACOBIAN_ENTRIES: [OsdiJacobianEntry; 4] = [
    OsdiJacobianEntry { nodes: OsdiNodePair { node_1: 0, node_2: 0 }, react_ptr_off: 0, flags: 0 },
    OsdiJacobianEntry { nodes: OsdiNodePair { node_1: 0, node_2: 1 }, react_ptr_off: 0, flags: 0 },
    OsdiJacobianEntry { nodes: OsdiNodePair { node_1: 1, node_2: 0 }, react_ptr_off: 0, flags: 0 },
    OsdiJacobianEntry { nodes: OsdiNodePair { node_1: 1, node_2: 1 }, react_ptr_off: 0, flags: 0 },
];

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
    num_params: 0,        // gd is hardcoded; no user-settable parameters
    num_instance_params: 0,
    num_opvars: 0,
    param_opvar: std::ptr::null_mut(),

    node_mapping_offset: 0,         // [u32; 2] at byte 0 of instance buffer
    jacobian_ptr_resist_offset: 8,  // unused (copy-based path); must not alias node_mapping
    num_states: 0,
    state_idx_off: 0,
    bound_step_offset: 0,
    instance_size: 8,   // 2 × u32 = 8 bytes
    model_size: 8,      // 1 × f64 = 8 bytes

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
    num_reactive_jacobian_entries: 0,
    write_jacobian_array_resist: Some(write_jacobian_array_resist_impl),
    write_jacobian_array_react: None,
    num_inputs: 0,
    inputs: std::ptr::null_mut(),
    load_jacobian_with_offset_resist: None,
    load_jacobian_with_offset_react: None,
}];
