use std::collections::HashMap;
use std::sync::Arc;

use fairchild_parser::ModelCard;

use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;
use crate::models::{
    pn_phase_shifter, pn_phase_shifter_cap, pn_phase_shifter_full, pn_phase_shifter_inj,
    pn_thermal_phase_shifter, pn_thermal_phase_shifter_cap, pn_thermal_phase_shifter_full,
    pn_thermal_phase_shifter_inj, thermal_phase_shifter, thermal_rc_phase_shifter, GummelPoonBjt,
    Mosfet1, NativeCirculator, NativeCwLaser, NativeDemux, NativeDirectionalCoupler,
    NativeGratingCoupler, NativeMux, NativeMzm, NativePhotodetector, NativeSplitter,
    NativeWaveguide, ShockleyDiode,
};

// Factory closures are `Arc` so the alias mechanism (B6) can clone a target
// factory into an outer wrapper that performs parameter-name translation.
type Factory = Arc<dyn Fn(&[NodeId], &SimContext) -> Box<dyn Device> + Send + Sync + 'static>;

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
    /// BJT model cards: model_name → (is_pnp, params).
    pub(crate) bjt_cards: HashMap<String, (bool, Vec<(String, f64)>)>,
}

impl DeviceRegistry {
    pub fn new() -> Self {
        let mut reg = Self {
            factories: HashMap::new(),
            mosfet_cards: HashMap::new(),
            bjt_cards: HashMap::new(),
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
        self.factories.insert(name.into(), Arc::new(factory));
    }

    /// Register `new_name` as an alias of an existing factory, optionally
    /// translating parameter names on `set_real_param` calls.
    ///
    /// This is the **PDK-adapter hook** (Phase B6): a user-supplied PDK
    /// mapping table can register foundry device names against native
    /// devices without leaking PDK-specific code into master.  Example:
    ///
    /// ```rust,ignore
    /// use std::collections::HashMap;
    /// let mut reg = DeviceRegistry::new();
    /// let mut remap = HashMap::new();
    /// remap.insert("waveguide_length_um".into(), "l_um".into());
    /// remap.insert("group_index".into(),        "n_g".into());
    /// reg.register_alias("pdk_foo_waveguide", "fc_waveguide", remap).unwrap();
    /// ```
    ///
    /// After registration, the netlist token `pdk_foo_waveguide` builds an
    /// `fc_waveguide` instance, with PDK-named parameters routed to the
    /// native device's underlying parameter names through the map.  Unknown
    /// parameter names pass through unchanged.
    pub fn register_alias(
        &mut self,
        new_name: impl Into<String>,
        target_name: &str,
        param_remap: HashMap<String, String>,
    ) -> Result<(), String> {
        let target_factory = self.factories.get(target_name).cloned().ok_or_else(|| {
            format!(
                "register_alias: unknown target factory '{target_name}' \
                 (register it before creating aliases)"
            )
        })?;
        let remap = Arc::new(param_remap);
        let aliased_factory = Arc::new(move |terminals: &[NodeId], ctx: &SimContext| {
            let inner = target_factory(terminals, ctx);
            Box::new(AliasedDevice {
                inner,
                remap: Arc::clone(&remap),
            }) as Box<dyn Device>
        });
        self.factories.insert(new_name.into(), aliased_factory);
        Ok(())
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
            self.mosfet_cards
                .insert(card.name.clone(), (is_pmos, card.params.clone()));
        }
    }

    /// Populate the registry from `.model NPN` / `.model PNP` cards.
    pub fn register_builtin_bjts(&mut self, cards: &[ModelCard]) {
        for card in cards {
            let kind = card.kind.to_lowercase();
            let is_pnp = match kind.as_str() {
                "npn" => false,
                "pnp" => true,
                _ => continue,
            };
            let (_, unknown) = GummelPoonBjt::from_model_params(is_pnp, &card.params);
            if !unknown.is_empty() {
                eprintln!(
                    "warning: BJT model '{}' params not yet implemented (using defaults): {}",
                    card.name,
                    unknown.join(", ")
                );
            }
            self.bjt_cards
                .insert(card.name.clone(), (is_pnp, card.params.clone()));
        }
    }

    /// Build a `GummelPoonBjt` instance for a `Q` element, injecting the
    /// stored model-card parameters. Returns `None` if the model name is unknown.
    pub(crate) fn build_bjt(
        &self,
        model_name: &str,
        terminals: &[NodeId],
        ctx: &SimContext,
    ) -> Option<Box<dyn Device>> {
        let (is_pnp, model_params) = self.bjt_cards.get(model_name)?;
        let (mut dev, _) = GummelPoonBjt::from_model_params(*is_pnp, model_params);
        dev.setup_model(ctx);
        dev.setup_instance(terminals, ctx);
        Some(Box::new(dev))
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
        self.register("fc_grating_coupler", |terminals, ctx| {
            let mut d = NativeGratingCoupler::new();
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
            let mut d = thermal_phase_shifter();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_thermal_ps_rc", |terminals, ctx| {
            let mut d = thermal_rc_phase_shifter();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_pn_ps", |terminals, ctx| {
            // Collapsed onto ActiveOpticalDevice (OpticalSegment + PnDrive).
            let mut d = pn_phase_shifter();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_pn_ps_cap", |terminals, ctx| {
            let mut d = pn_phase_shifter_cap();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_pn_th_ps", |terminals, ctx| {
            let mut d = pn_thermal_phase_shifter();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_pn_th_ps_cap", |terminals, ctx| {
            let mut d = pn_thermal_phase_shifter_cap();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_pn_ps_inj", |terminals, ctx| {
            let mut d = pn_phase_shifter_inj();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_pn_th_ps_inj", |terminals, ctx| {
            let mut d = pn_thermal_phase_shifter_inj();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_pn_ps_full", |terminals, ctx| {
            let mut d = pn_phase_shifter_full();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_pn_th_ps_full", |terminals, ctx| {
            let mut d = pn_thermal_phase_shifter_full();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_mzm", |terminals, ctx| {
            let mut d = NativeMzm::new();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_circulator", |terminals, ctx| {
            let mut d = NativeCirculator::new();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_cw_laser", |terminals, ctx| {
            let mut d = NativeCwLaser::new();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_mux", |terminals, ctx| {
            let mut d = NativeMux::new();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            Box::new(d) as Box<dyn Device>
        });
        self.register("fc_demux", |terminals, ctx| {
            let mut d = NativeDemux::new();
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

// ─── PDK adapter (B6) ────────────────────────────────────────────────────

/// Wraps another `Device`, translating parameter names on `set_real_param`
/// through a fixed remap table.  Used by `DeviceRegistry::register_alias`
/// to surface foundry-specific device names with native fairchild devices
/// underneath.  All other Device-trait methods forward verbatim.
struct AliasedDevice {
    inner: Box<dyn Device>,
    remap: Arc<HashMap<String, String>>,
}

impl Device for AliasedDevice {
    fn num_terminals(&self) -> usize {
        self.inner.num_terminals()
    }
    fn setup_model(&mut self, ctx: &SimContext) {
        self.inner.setup_model(ctx)
    }
    fn setup_instance(&mut self, terminals: &[NodeId], ctx: &SimContext) {
        self.inner.setup_instance(terminals, ctx)
    }
    fn eval(&mut self, x: &[f64], flags: EvalFlags, ctx: &SimContext) {
        self.inner.eval(x, flags, ctx)
    }
    fn load_residual(&self, b: &mut [f64]) {
        self.inner.load_residual(b)
    }
    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        self.inner.load_jacobian(mat)
    }
    fn load_residual_tran(&self, b: &mut [f64], alpha: f64) {
        self.inner.load_residual_tran(b, alpha)
    }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, alpha: f64) {
        self.inner.load_jacobian_tran(mat, alpha)
    }
    fn commit_timestep(&mut self, x: &[f64]) {
        self.inner.commit_timestep(x)
    }
    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        // Look up the remap (case-insensitive); fall through unchanged if
        // not found so users can mix PDK-named and native-named params.
        let key = name.to_lowercase();
        let translated = self.remap.get(&key).map(String::as_str).unwrap_or(name);
        self.inner.set_real_param(translated, value)
    }
    fn num_extra_nodes(&self) -> usize {
        self.inner.num_extra_nodes()
    }
    fn bind_extra_nodes(&mut self, first_idx: usize) {
        self.inner.bind_extra_nodes(first_idx)
    }
    fn noise_sources(&self, ctx: &SimContext) -> Vec<(NodeId, NodeId, f64)> {
        self.inner.noise_sources(ctx)
    }
    fn small_signal_reactances(&self) -> Vec<crate::device::ReactiveBranchSpec> {
        self.inner.small_signal_reactances()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

    /// Register a PDK-style alias mapping `pdk_widget_wg` →
    /// `fc_waveguide` with parameter-name remapping, then verify a netlist
    /// using the PDK name + PDK param names builds the right device and
    /// produces the same numbers as the native form.
    #[test]
    fn register_alias_pdk_name_with_param_remap() {
        let netlist = parse_spice(
            "* aliased waveguide\n\
             V_re in_re 0 DC 1.0\n\
             V_im in_im 0 DC 0.0\n\
             V_wl in_wl 0 DC 1.55e-6\n\
             X1 in_re in_im in_wl out_re out_im out_wl pdk_widget_wg \
                wg_len_um=100 mode_index=4.2 prop_loss_dB_cm=2.0\n\
             .op\n.end\n",
        )
        .unwrap();
        let mut registry = DeviceRegistry::new();
        let mut remap = HashMap::new();
        remap.insert("wg_len_um".to_string(), "l_um".to_string());
        remap.insert("mode_index".to_string(), "n_g".to_string());
        remap.insert("prop_loss_db_cm".to_string(), "alpha_db_cm".to_string());
        registry
            .register_alias("pdk_widget_wg", "fc_waveguide", remap)
            .expect("alias should register");
        let r = crate::newton::dc_op_nr_with_registry(&netlist, &registry)
            .expect("DC OP should converge");
        // Same numbers as the un-aliased fc_waveguide test elsewhere.
        let v_re = r.node_voltage("out_re").unwrap();
        let v_im = r.node_voltage("out_im").unwrap();
        let amp = (v_re * v_re + v_im * v_im).sqrt();
        let expected = (-23.0258509_f64 * 100e-6 / 2.0).exp();
        assert!(
            (amp - expected).abs() < 1e-5,
            "|A_out|={amp:.6} expected={expected:.6}"
        );
    }

    #[test]
    fn register_alias_unknown_target_errors() {
        let mut reg = DeviceRegistry::new();
        let res = reg.register_alias("pdk_widget_wg", "nonexistent_target", HashMap::new());
        assert!(res.is_err());
    }
}
