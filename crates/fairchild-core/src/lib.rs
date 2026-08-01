pub mod ac;
pub mod behavioral;
pub mod connectivity;
pub mod dc_sweep;
pub mod delay;
pub mod device;
pub mod device_registry;
pub mod error;
pub mod measure;
pub mod mna;
pub mod models;
pub mod netlist_edit;
pub mod newton;
pub mod noise;
pub mod options;
pub mod reactive;
pub mod sanity;
pub mod solver;
pub mod tolerance;
pub mod tran;
pub mod tran_step;

pub use ac::{ac_analysis, ac_analysis_opts, freq_decade, freq_linear, freq_oct, AcResult};
pub use connectivity::check_connectivity;
pub use dc_sweep::{dc_sweep_with_registry, dc_sweep_with_registry_opts, DcSweepResult, SweepAxis};
pub use delay::DelayLine;
pub use device::{Device, EvalFlags, NodeId, SimContext};
pub use device_registry::{DeviceRegistry, ModelFactory, ParamSet};
pub use error::SimError;
pub use measure::{evaluate_measurements, MeasureResult};
pub use mna::CircuitTopology;
pub use models::{GummelPoonBjt, Mosfet1, ShockleyDiode};
pub use netlist_edit::{set_element_param, set_source_pwl};
pub use newton::{
    build_devices, dc_op_nr, dc_op_nr_opts, dc_op_nr_with_devices, dc_op_nr_with_devices_opts,
    dc_op_nr_with_registry, dc_op_nr_with_registry_opts, NrResult,
};
pub use noise::{noise_analysis, NoiseResult};
pub use options::SimOptions;
pub use sanity::check_netlist_sanity;
pub use solver::{lu_solve, SolverKind};
pub use tran::{
    run_tran, run_tran_tr, tran_nr, tran_nr_tr, tran_nr_var, tran_nr_with_registry,
    tran_nr_with_registry_opts, tran_nr_with_registry_tr, tran_nr_with_registry_var,
    tran_nr_with_registry_var_opts, IntegratorMode, TranResult,
};
pub use tran_step::TranStepper;
