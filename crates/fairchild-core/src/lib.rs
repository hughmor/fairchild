pub mod ac;
pub mod behavioral;
pub mod connectivity;
pub mod dc_sweep;
pub mod device;
pub mod measure;
pub mod device_registry;
pub mod error;
pub mod mna;
pub mod models;
pub mod newton;
pub mod noise;
pub mod options;
pub mod sanity;
pub mod solver;
pub mod tran;

pub use ac::{ac_analysis, ac_analysis_opts, freq_decade, freq_linear, freq_oct, AcResult};
pub use connectivity::check_connectivity;
pub use dc_sweep::{dc_sweep_with_registry, dc_sweep_with_registry_opts, DcSweepResult, SweepAxis};
pub use measure::{evaluate_measurements, MeasureResult};
pub use noise::{noise_analysis, NoiseResult};
pub use device::{Device, EvalFlags, NodeId, SimContext};
pub use device_registry::DeviceRegistry;
pub use error::SimError;
pub use mna::CircuitTopology;
pub use models::{GummelPoonBjt, Mosfet1, ShockleyDiode};
pub use newton::{
    build_devices, dc_op_nr, dc_op_nr_opts,
    dc_op_nr_with_devices, dc_op_nr_with_devices_opts,
    dc_op_nr_with_registry, dc_op_nr_with_registry_opts,
    NrResult,
};
pub use options::SimOptions;
pub use sanity::check_netlist_sanity;
pub use solver::{lu_solve, SolverKind};
pub use tran::{
    run_tran, run_tran_tr,
    tran_nr, tran_nr_tr, tran_nr_var,
    tran_nr_with_registry, tran_nr_with_registry_opts,
    tran_nr_with_registry_tr, tran_nr_with_registry_var,
    tran_nr_with_registry_var_opts,
    IntegratorMode, TranResult,
};
