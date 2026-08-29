// A doc link that does not resolve is a comment describing code that is not
// there, and this crate has shipped two of them: a link to a `LambdaSelect` that
// was deleted, and a `build` method on a type alias for a function. `cargo doc`
// found both and nobody was reading it, so the class is denied here instead.
//
// Private targets are allowed on purpose. The item really is there, and a link to
// it helps whoever is reading the source — which for a simulator is most readers.
// Only the rendered page loses, and only for an item nobody outside can call.
#![deny(rustdoc::broken_intra_doc_links)]
#![allow(rustdoc::private_intra_doc_links)]

pub mod ac;
pub mod adjoint;
pub mod adjoint_ac;
pub mod adjoint_tran;
pub mod behavioral;
pub mod binning;
pub mod connectivity;
pub mod corner;
pub mod dc_sweep;
pub mod delay;
pub mod device;
pub mod device_registry;
pub mod error;
pub mod lambda;
pub mod measure;
pub mod mna;
pub mod models;
pub mod netlist_edit;
pub mod newton;
pub mod noise;
pub mod options;
pub mod pz;
pub mod reactive;
pub mod sanity;
pub mod sens;
pub mod solver;
pub mod temperature;
pub mod tf;
pub mod tolerance;
pub mod tran;
pub mod tran_step;
pub mod unmodelled;

// The warning switch lives in the parser (the crate everything else depends on)
// and is re-exported here so a frontend needs one import, not two, to make
// `--quiet` mean what it says.
pub use fairchild_parser::warn;
pub use fairchild_parser::warn::set_quiet;
pub use fairchild_parser::warn_user;

pub use ac::{ac_analysis, ac_analysis_opts, freq_decade, freq_linear, freq_oct, AcResult};
pub use adjoint::{dc_sensitivity, Output, ParamRef, Sensitivities};
pub use adjoint_tran::{TranAdjoint, TranSensitivities};
pub use connectivity::check_connectivity;
pub use corner::{expand_corners, Corner, CornerGrid};
pub use dc_sweep::{dc_sweep_with_registry, dc_sweep_with_registry_opts, DcSweepResult, SweepAxis};
pub use delay::DelayLine;
pub use device::{Device, EvalFlags, NodeId, SimContext};
pub use device_registry::{ArityDecl, DeviceRegistry, ModelFactory, ParamSet};
pub use error::SimError;
pub use measure::{evaluate_measurements, MeasureResult};
pub use mna::CircuitTopology;
pub use models::{GummelPoonBjt, Mosfet1, ShockleyDiode};
pub use netlist_edit::{get_element_param, set_element_param, set_source_pwl};
pub use newton::{
    build_devices, dc_op_nr, dc_op_nr_opts, dc_op_nr_with_devices, dc_op_nr_with_devices_opts,
    dc_op_nr_with_registry, dc_op_nr_with_registry_opts, NrResult,
};
pub use noise::{noise_analysis, NoiseResult};
pub use options::SimOptions;
pub use pz::{pole_zero, PzResult, Root, MAX_PZ_SIZE};
pub use sanity::check_netlist_sanity;
pub use sens::{sensitivity, SensResult, SensRow};
pub use solver::SolverKind;
pub use tf::{transfer_function, TfResult};
pub use tran::{
    tran_nr, tran_nr_configured, tran_nr_tr, tran_nr_var, tran_nr_with_registry,
    tran_nr_with_registry_opts, tran_nr_with_registry_tr, tran_nr_with_registry_var,
    tran_nr_with_registry_var_opts, IntegratorMode, TranResult,
};
pub use tran_step::TranStepper;
