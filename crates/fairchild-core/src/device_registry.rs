use std::collections::HashMap;

use fairchild_parser::ModelCard;

use crate::device::{Device, NodeId, SimContext};
use crate::models::{
    Mosfet1, NativeDirectionalCoupler, NativePhotodetector, NativePnPhaseShifter,
    NativeSplitter, NativeThermalPhaseShifter, NativeWaveguide, ShockleyDiode,
};

type Factory =
    Box<dyn Fn(&[NodeId], &SimContext) -> Box<dyn Device> + Send + Sync + 'static>;

/// Maps model names to device factory closures.
///
/// Factories receive the MNA node mapping and SimContext, and return a fully
/// initialised (setup_model + setup_instance already called) boxed Device.
///
/// Built-in models are registered via `register_builtin_diodes` /
/// `register_builtin_mosfets`. External models (e.g. OSDI) register themselves
/// by capturing an Arc to their library, keeping it alive for the device's lifetime.
pub struct DeviceRegistry {
    factories: HashMap<String, Factory>,
    /// MOSFET model cards stored for W/L instance-param injection in build_devices.
    pub(crate) mosfet_cards: HashMap<String, (bool, Vec<(String, f64)>)>,
}

impl DeviceRegistry {
    pub fn new() -> Self {
        let mut reg = Self {
            factories: HashMap::new(),
            mosfet_cards: HashMap::new(),
        };
        // Native photonic passives are always available — no .model card or
        // .osdi import required to instantiate `fc_waveguide`, `fc_dcoupler`,
        // `fc_splitter`.
        reg.register_native_photonics();
        reg
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
            if !card.kind.to_lowercase().starts_with('d') {
                continue;
            }
            let params: Vec<(String, f64)> = card.params.clone();
            // Warn once per model card about params that aren't yet implemented.
            let (_, unknown) = ShockleyDiode::from_params(&params);
            if !unknown.is_empty() {
                eprintln!(
                    "warning: diode model '{}' params not yet implemented (using defaults): {}",
                    card.name,
                    unknown.join(", ")
                );
            }
            self.register(card.name.clone(), move |terminals, ctx| {
                let (mut dev, _) = ShockleyDiode::from_params(&params);
                dev.setup_model(ctx);
                dev.setup_instance(terminals, ctx);
                Box::new(dev)
            });
        }
    }

    /// Populate the registry from `.model NMOS` / `.model PMOS` cards.
    ///
    /// MOSFET factories do not accept instance W/L here; those are injected by
    /// `build_devices` at instantiation time using the stored `mosfet_cards` map.
    pub fn register_builtin_mosfets(&mut self, cards: &[ModelCard]) {
        for card in cards {
            let kind = card.kind.to_lowercase();
            let is_pmos = match kind.as_str() {
                "nmos" => false,
                "pmos" => true,
                _ => continue,
            };
            // Warn once per model card about unrecognised params.
            let (_, unknown) = Mosfet1::from_model_params(is_pmos, &card.params);
            if !unknown.is_empty() {
                eprintln!(
                    "warning: MOSFET model '{}' params not yet implemented (using defaults): {}",
                    card.name,
                    unknown.join(", ")
                );
            }
            self.mosfet_cards.insert(
                card.name.clone(),
                (is_pmos, card.params.clone()),
            );
        }
    }

    /// Register the always-available native photonic devices.
    /// All addressable by model name in any X-element instance line:
    ///
    /// Passives (B3):
    /// - `fc_waveguide`     — straight waveguide.
    /// - `fc_dcoupler`      — 2×2 directional coupler.
    /// - `fc_splitter`      — 1×2 Y-junction (3 dB lossless).
    ///
    /// Actives (B4):
    /// - `fc_photodetector` — PIN photodetector with linear responsivity.
    /// - `fc_thermal_ps`    — thermal phase shifter (Joule heating → φ = π·P/P_pi).
    /// - `fc_pn_ps`         — PN-junction phase shifter (Δn_eff = dn_dv·V).
    pub fn register_native_photonics(&mut self) {
        self.register("fc_waveguide", |terminals, ctx| {
            let mut d = NativeWaveguide::new();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_dcoupler", |terminals, ctx| {
            let mut d = NativeDirectionalCoupler::new();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_splitter", |terminals, ctx| {
            let mut d = NativeSplitter::new();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_photodetector", |terminals, ctx| {
            let mut d = NativePhotodetector::new();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_thermal_ps", |terminals, ctx| {
            let mut d = NativeThermalPhaseShifter::new();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_pn_ps", |terminals, ctx| {
            let mut d = NativePnPhaseShifter::new();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
    }

    /// Build a `Mosfet1` for `model_name` with specific instance params (W, L).
    pub(crate) fn build_mosfet(
        &self,
        model_name: &str,
        instance_params: &[(String, f64)],
        terminals: &[NodeId],
        ctx: &SimContext,
    ) -> Option<Box<dyn Device>> {
        let (is_pmos, model_params) = self.mosfet_cards.get(model_name)?;
        let (mut dev, _) = Mosfet1::from_model_params(*is_pmos, model_params);
        dev.set_instance_params(instance_params);
        dev.setup_model(ctx);
        dev.setup_instance(terminals, ctx);
        Some(Box::new(dev))
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
