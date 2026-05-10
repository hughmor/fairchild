//! Minimal OSDI v0.4 mock library for integration testing.
//!
//! Exports one model named "test_diode" with two terminals.
//! All function pointers are null (None) — this is only used for registry-walk tests.

use std::ffi::c_char;
use std::os::raw::c_void;

use fairchild_osdi::ffi::{OsdiDescriptor, OsdiSimParas};

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
    name: c"test_diode".as_ptr().cast_mut(),
    num_nodes: 2,
    num_terminals: 2,
    nodes: std::ptr::null_mut(),
    num_jacobian_entries: 0,
    jacobian_entries: std::ptr::null_mut(),
    num_collapsible: 0,
    collapsible: std::ptr::null_mut(),
    collapsed_offset: 0,
    noise_sources: std::ptr::null_mut(),
    num_noise_src: 0,
    num_params: 1,
    num_instance_params: 0,
    num_opvars: 0,
    param_opvar: std::ptr::null_mut(),
    node_mapping_offset: 0,
    jacobian_ptr_resist_offset: 0,
    num_states: 0,
    state_idx_off: 0,
    bound_step_offset: 0,
    instance_size: 8,
    model_size: 8,
    access: None,
    setup_model: None,
    setup_instance: None,
    eval: None,
    load_noise: None,
    load_residual_resist: None,
    load_residual_react: None,
    load_limit_rhs_resist: None,
    load_limit_rhs_react: None,
    load_spice_rhs_dc: None,
    load_spice_rhs_tran: None,
    load_jacobian_resist: None,
    load_jacobian_react: None,
    load_jacobian_tran: None,
    given_flag_model: None,
    given_flag_instance: None,
    num_resistive_jacobian_entries: 0,
    num_reactive_jacobian_entries: 0,
    write_jacobian_array_resist: None,
    write_jacobian_array_react: None,
    num_inputs: 0,
    inputs: std::ptr::null_mut(),
    load_jacobian_with_offset_resist: None,
    load_jacobian_with_offset_react: None,
}];

// Suppress unused import warning for OsdiSimParas (used in function pointer types
// that are None here; included so the import doesn't get flagged).
const _: () = {
    let _ = std::mem::size_of::<OsdiSimParas>();
};
