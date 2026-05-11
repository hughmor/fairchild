pub mod ac;
pub mod device;
pub mod device_registry;
pub mod error;
pub mod mna;
pub mod models;
pub mod newton;
pub mod solver;
pub mod tran;

pub use ac::{ac_analysis, freq_decade, freq_linear, AcResult};
pub use device::{Device, EvalFlags, NodeId, SimContext};
pub use device_registry::DeviceRegistry;
pub use error::SimError;
pub use models::ShockleyDiode;
pub use newton::{dc_op_nr, dc_op_nr_with_registry, NrResult};
pub use solver::{dc_op, OpResult};
pub use tran::{
    run_tran, run_tran_tr,
    tran_nr, tran_nr_tr, tran_nr_var, tran_nr_with_registry, tran_nr_with_registry_tr,
    tran_nr_with_registry_var,
    IntegratorMode, TranResult,
};
