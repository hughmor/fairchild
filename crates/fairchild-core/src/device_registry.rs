use std::collections::HashMap;

use fairchild_parser::ModelCard;

use crate::device::{Device, NodeId, SimContext};
use crate::models::ShockleyDiode;

type Factory =
    Box<dyn Fn(&[NodeId], &SimContext) -> Box<dyn Device> + Send + Sync + 'static>;

/// Maps model names to device factory closures.
///
/// Factories receive the MNA node mapping and SimContext, and return a fully
/// initialised (setup_model + setup_instance already called) boxed Device.
///
/// Built-in models are registered via `register_builtin_diodes`.
/// External models (e.g. OSDI) register themselves by capturing an Arc to
/// their library, keeping it alive for the device's lifetime.
pub struct DeviceRegistry {
    factories: HashMap<String, Factory>,
}

impl DeviceRegistry {
    pub fn new() -> Self {
        Self { factories: HashMap::new() }
    }

    /// Register a factory for `name`. Overwrites any previous entry.
    pub fn register(
        &mut self,
        name: impl Into<String>,
        factory: impl Fn(&[NodeId], &SimContext) -> Box<dyn Device> + Send + Sync + 'static,
    ) {
        self.factories.insert(name.into(), Box::new(factory));
    }

    /// Populate the registry from `.model D` cards using the built-in Shockley diode.
    pub fn register_builtin_diodes(&mut self, cards: &[ModelCard]) {
        for card in cards {
            if !card.kind.starts_with('d') {
                continue;
            }
            let params: Vec<(String, f64)> = card.params.clone();
            self.register(card.name.clone(), move |terminals, ctx| {
                let mut dev = Box::new(ShockleyDiode::from_params(&params));
                dev.setup_model(ctx);
                dev.setup_instance(terminals, ctx);
                dev
            });
        }
    }

    /// Look up a factory by model name.
    pub fn get(&self, name: &str) -> Option<&Factory> {
        self.factories.get(name)
    }
}

impl Default for DeviceRegistry {
    fn default() -> Self {
        Self::new()
    }
}
