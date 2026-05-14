pub mod ac;
pub mod connectivity;
pub mod dc_sweep;
pub mod device;
pub mod device_registry;
pub mod error;
pub mod mna;
pub mod models;
pub mod newton;
pub mod options;
pub mod solver;
pub mod tran;

pub use ac::{ac_analysis, ac_analysis_opts, freq_decade, freq_linear, freq_oct, AcResult};
pub use connectivity::check_connectivity;
pub use dc_sweep::{dc_sweep_with_registry, dc_sweep_with_registry_opts, DcSweepResult, SweepAxis};
pub use device::{Device, EvalFlags, NodeId, SimContext};
pub use device_registry::DeviceRegistry;
pub use error::SimError;
pub use mna::CircuitTopology;
pub use models::{Mosfet1, ShockleyDiode};
pub use newton::{
    build_devices, dc_op_nr, dc_op_nr_opts,
    dc_op_nr_with_devices, dc_op_nr_with_devices_opts,
    dc_op_nr_with_registry, dc_op_nr_with_registry_opts,
    NrResult,
};
pub use options::SimOptions;
pub use solver::lu_solve;
pub use tran::{
    run_tran, run_tran_tr,
    tran_nr, tran_nr_tr, tran_nr_var,
    tran_nr_with_registry, tran_nr_with_registry_opts,
    tran_nr_with_registry_tr, tran_nr_with_registry_var,
    tran_nr_with_registry_var_opts,
    IntegratorMode, TranResult,
};
