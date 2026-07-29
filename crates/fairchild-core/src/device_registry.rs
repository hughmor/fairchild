use std::collections::HashMap;
use std::sync::Arc;

use fairchild_parser::{Expr, ModelCard};

use crate::device::{Device, NodeId, SimContext};
use crate::models::{
    expr_phase_shifter, pn_phase_shifter, pn_phase_shifter_cap, pn_phase_shifter_full,
    pn_phase_shifter_inj, pn_thermal_phase_shifter, pn_thermal_phase_shifter_cap,
    pn_thermal_phase_shifter_full, pn_thermal_phase_shifter_inj, thermal_phase_shifter,
    thermal_rc_phase_shifter, ActiveOpticalDevice, GummelPoonBjt, Mosfet1, NativeAwgr,
    NativeCirculator, NativeCwLaser, NativeDemux, NativeDirectionalCoupler, NativeGratingCoupler,
    NativeMux, NativeMzm, NativeOptical2x2, NativePhotodetector, NativeSplitter, NativeWaveguide,
    ShockleyDiode, SpectrumTable,
};

/// The set of instance parameters from an `X…` element line, threaded into a
/// [`ModelFactory::build`] call. Keys are case-insensitive; the set tracks which
/// keys a device actually consumed so unrecognised (typo'd) params can be
/// reported, and can rename keys for the PDK-alias path.
///
/// Created per device instantiation and passed by reference; it is *not* stored
/// in the registry, so its interior mutability (consumed-tracking) needs no
/// `Sync`.
pub struct ParamSet {
    /// Lower-cased (key, value) pairs in netlist order.
    params: Vec<(String, f64)>,
    consumed: std::cell::RefCell<Vec<bool>>,
}

impl ParamSet {
    /// Build from raw instance params (keys are lower-cased).
    pub fn new(params: &[(String, f64)]) -> Self {
        let params: Vec<(String, f64)> =
            params.iter().map(|(k, v)| (k.to_lowercase(), *v)).collect();
        let n = params.len();
        ParamSet {
            params,
            consumed: std::cell::RefCell::new(vec![false; n]),
        }
    }

    /// An empty set (devices with no instance params, e.g. diodes/passives).
    pub fn empty() -> Self {
        ParamSet::new(&[])
    }

    /// Look up a single param by name (case-insensitive), marking it consumed.
    pub fn get(&self, name: &str) -> Option<f64> {
        let nl = name.to_lowercase();
        for (i, (k, v)) in self.params.iter().enumerate() {
            if *k == nl {
                self.consumed.borrow_mut()[i] = true;
                return Some(*v);
            }
        }
        None
    }

    /// The `LEVEL` selector, if present (for construction-time dispatch).
    pub fn level(&self) -> Option<i64> {
        self.get("level").map(|v| v.round() as i64)
    }

    /// Apply every param to `dev` via `set_real_param`, marking each consumed
    /// iff the device accepted it. The single, uniform instance-param
    /// application path — preserves netlist order, so a later param overrides an
    /// earlier model-card default exactly as before.
    pub fn apply(&self, dev: &mut dyn Device) {
        let mut consumed = self.consumed.borrow_mut();
        for (i, (k, v)) in self.params.iter().enumerate() {
            if dev.set_real_param(k, *v) {
                consumed[i] = true;
            }
        }
    }

    /// Keys not consumed by `get`/`apply` — i.e. params the device did not
    /// recognise (likely typos). Used to warn the user.
    pub fn unconsumed(&self) -> Vec<String> {
        let consumed = self.consumed.borrow();
        self.params
            .iter()
            .zip(consumed.iter())
            .filter(|(_, c)| !**c)
            .map(|((k, _), _)| k.clone())
            .collect()
    }

    /// A copy with keys renamed through `remap` (the PDK-alias path): the
    /// translation happens on the `ParamSet` *before* the device sees the
    /// params, so it covers construction-time params too, not just
    /// `set_real_param`.
    pub fn renamed(&self, remap: &HashMap<String, String>) -> ParamSet {
        let renamed: Vec<(String, f64)> = self
            .params
            .iter()
            .map(|(k, v)| (remap.get(k).cloned().unwrap_or_else(|| k.clone()), *v))
            .collect();
        ParamSet::new(&renamed)
    }
}

/// A device factory: constructs a device from its terminal nodes, instance
/// [`ParamSet`], and the [`SimContext`]. The returned device is fully set up
/// (`setup_model` + `setup_instance` done) with the instance params applied; the
/// caller (`build_devices`) handles extra-node allocation and unconsumed-param
/// warnings. Stored behind `Arc` so the alias mechanism (B6) can clone a target
/// factory into a wrapper that performs parameter-name translation. (A named
/// `ModelFactory` trait could wrap this later for a dlopen plugin ABI; the
/// closure form covers every in-tree need — OSDI captures `Arc<library>`, the
/// LEVEL/expr factories capture their tables/expressions.)
pub type ModelFactory =
    dyn Fn(&[NodeId], &ParamSet, &SimContext) -> Box<dyn Device> + Send + Sync + 'static;
type Factory = Arc<ModelFactory>;

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
        factory: impl Fn(&[NodeId], &ParamSet, &SimContext) -> Box<dyn Device> + Send + Sync + 'static,
    ) {
        self.factories.insert(name.into(), Arc::new(factory));
    }

    /// Register a default-constructible device: `T::default()` + the standard
    /// setup + instance-param application. Collapses the boilerplate for every
    /// passive/simple native device.
    pub fn register_default<T: Device + Default + 'static>(&mut self, name: impl Into<String>) {
        self.register(name, |terminals, params: &ParamSet, ctx| {
            let mut d = T::default();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            params.apply(&mut d);
            Box::new(d) as Box<dyn Device>
        });
    }

    /// Register a device built by a constructor function (the active-photonic
    /// `ActiveOpticalDevice` builders, which are not `Default`). Same setup +
    /// instance-param flow as [`register_default`](Self::register_default).
    fn register_ctor(&mut self, name: impl Into<String>, ctor: fn() -> ActiveOpticalDevice) {
        self.register(name, move |terminals, params: &ParamSet, ctx| {
            let mut d = ctor();
            d.setup_model(ctx);
            d.setup_instance(terminals, ctx);
            params.apply(&mut d);
            Box::new(d) as Box<dyn Device>
        });
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
        let target_factory: Factory =
            self.factories.get(target_name).cloned().ok_or_else(|| {
                format!(
                    "register_alias: unknown target factory '{target_name}' \
                 (register it before creating aliases)"
                )
            })?;
        let remap = Arc::new(param_remap);
        // Translate the ParamSet's keys BEFORE the target factory consumes them,
        // so the remap covers construction-time params too — not just the
        // post-construction set_real_param path (the old AliasedDevice limitation).
        self.register(new_name, move |terminals, params: &ParamSet, ctx| {
            let translated = params.renamed(&remap);
            target_factory(terminals, &translated, ctx)
        });
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
            self.register(card.name.clone(), move |terminals, ps: &ParamSet, ctx| {
                let (mut dev, _) = ShockleyDiode::from_params(&params);
                dev.setup_model(ctx);
                dev.setup_instance(terminals, ctx);
                ps.apply(&mut dev);
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

    /// Register active-photonic `.model` cards with a `LEVEL` selector, à la
    /// MOSFET LEVEL. A card `.model myps fc_pn_ps LEVEL=2` registers `myps` as a
    /// depletion-cap PN phase shifter; the model-card params (L_um, dn_dv, …)
    /// are baked in, and per-instance params on the `X…` line still apply on
    /// top. The base type (`fc_pn_ps` / `fc_thermal_ps` / `fc_pn_th_ps`) is the
    /// family; LEVEL selects the electrical sophistication. LEVEL 1 (or absent)
    /// is the plain device, so a bare `.model` card is always valid.
    ///
    /// | base type | LEVEL | device |
    /// |---|---|---|
    /// | `fc_pn_ps` | 1 / 2 / 3 / 4 | linear / +Cj / injection / full |
    /// | `fc_thermal_ps` | 1 / 2 | instantaneous / thermal-RC |
    /// | `fc_pn_th_ps` | 1 / 2 / 3 / 4 | +heater on each PN level |
    ///
    /// Cards whose kind is not a photonic family are ignored (handled by the
    /// diode/MOSFET/BJT registrars).
    pub fn register_photonic_models(&mut self, cards: &[ModelCard]) {
        for card in cards {
            let kind = card.kind.to_lowercase();

            // An AWG router backed by measured spectra. The file path is a
            // string, and an X-line's instance params are numeric only, so a
            // `.model` card is the only route in:
            //
            //   .model awg8 fc_awgr sfile="awgr8.csv"
            //   Xr in0 … in7 out0 … out7 awg8
            //
            // The table is read once here, not per instance.
            if kind == "fc_awgr" {
                let path = card
                    .expr_params
                    .iter()
                    .find(|(k, _)| k.eq_ignore_ascii_case("sfile"))
                    .map(|(_, v)| v.clone());
                let numeric: Vec<(String, f64)> = card.params.clone();
                let name = card.name.clone();
                self.register(card.name.clone(), move |terminals, ps: &ParamSet, ctx| {
                    let mut d = NativeAwgr::new();
                    d.setup_model(ctx);
                    d.setup_instance(terminals, ctx);
                    if let Some(p) = &path {
                        match std::fs::read_to_string(p)
                            .map_err(|e| e.to_string())
                            .and_then(|t| SpectrumTable::from_csv(&t, d.n_ports()))
                        {
                            Ok(table) => d.set_table(table),
                            Err(e) => eprintln!(
                                "warning: fc_awgr model '{name}' could not load sfile=\"{p}\" \
                                 ({e}); falling back to the analytic response"
                            ),
                        }
                    }
                    for (k, v) in &numeric {
                        d.set_real_param(k, *v);
                    }
                    ps.apply(&mut d);
                    Box::new(d) as Box<dyn Device>
                });
                continue;
            }

            // Tier-1: a declarative expression-driven phase shifter. The
            // constitutive maps live in the card's expr_params; numeric params
            // (geometry, g_pn) apply on top.
            if kind == "fc_phase_shifter_expr" {
                let parse_map = |name: &str| -> Option<Expr> {
                    card.expr_params
                        .iter()
                        .find(|(k, _)| k.eq_ignore_ascii_case(name))
                        .and_then(|(_, src)| match Expr::parse(src) {
                            Ok(e) => Some(e),
                            Err(err) => {
                                eprintln!(
                                    "warning: photonic model '{}' {name} expression \
                                     failed to parse ({err:?}); treating as 0",
                                    card.name
                                );
                                None
                            }
                        })
                };
                let dneff = parse_map("dneff");
                let dalpha = parse_map("dalpha");
                let g_pn = card
                    .params
                    .iter()
                    .find(|(k, _)| k == "g_pn")
                    .map(|(_, v)| *v)
                    .unwrap_or(1e-3);
                let numeric: Vec<(String, f64)> = card
                    .params
                    .iter()
                    .filter(|(k, _)| !k.eq_ignore_ascii_case("level"))
                    .cloned()
                    .collect();
                self.register(card.name.clone(), move |terminals, ps: &ParamSet, ctx| {
                    let mut d = expr_phase_shifter(dneff.clone(), dalpha.clone(), g_pn);
                    d.setup_model(ctx);
                    d.setup_instance(terminals, ctx);
                    // Model-card numeric params first, then instance params on top.
                    for (k, v) in &numeric {
                        d.set_real_param(k, *v);
                    }
                    ps.apply(&mut d);
                    Box::new(d) as Box<dyn Device>
                });
                continue;
            }

            let level = card
                .params
                .iter()
                .find(|(k, _)| k.eq_ignore_ascii_case("level"))
                .map(|(_, v)| v.round() as i64)
                .unwrap_or(1);
            let ctor: fn() -> ActiveOpticalDevice = match (kind.as_str(), level) {
                ("fc_pn_ps", 1) => pn_phase_shifter,
                ("fc_pn_ps", 2) => pn_phase_shifter_cap,
                ("fc_pn_ps", 3) => pn_phase_shifter_inj,
                ("fc_pn_ps", 4) => pn_phase_shifter_full,
                ("fc_thermal_ps", 1) => thermal_phase_shifter,
                ("fc_thermal_ps", 2) => thermal_rc_phase_shifter,
                ("fc_pn_th_ps", 1) => pn_thermal_phase_shifter,
                ("fc_pn_th_ps", 2) => pn_thermal_phase_shifter_cap,
                ("fc_pn_th_ps", 3) => pn_thermal_phase_shifter_inj,
                ("fc_pn_th_ps", 4) => pn_thermal_phase_shifter_full,
                // Not a photonic family (or unknown LEVEL) — leave for others;
                // warn on a recognised family with a bad LEVEL.
                (k @ ("fc_pn_ps" | "fc_thermal_ps" | "fc_pn_th_ps"), bad) => {
                    eprintln!(
                        "warning: photonic model '{}' has unsupported {k} LEVEL={bad}; \
                         using LEVEL=1",
                        card.name
                    );
                    match k {
                        "fc_thermal_ps" => thermal_phase_shifter,
                        "fc_pn_th_ps" => pn_thermal_phase_shifter,
                        _ => pn_phase_shifter,
                    }
                }
                _ => continue,
            };
            // Model-card params (minus LEVEL) bake in at construction; instance
            // params on the X-line apply on top (via the ParamSet).
            let params: Vec<(String, f64)> = card
                .params
                .iter()
                .filter(|(k, _)| !k.eq_ignore_ascii_case("level"))
                .cloned()
                .collect();
            self.register(card.name.clone(), move |terminals, ps: &ParamSet, ctx| {
                let mut d = ctor();
                d.setup_model(ctx);
                d.setup_instance(terminals, ctx);
                for (k, v) in &params {
                    d.set_real_param(k, *v);
                }
                ps.apply(&mut d);
                Box::new(d) as Box<dyn Device>
            });
        }
    }

    /// Register every `.model`-card-derived built-in family from a netlist's
    /// cards in one call: diodes, MOSFETs, BJTs, and active-photonic LEVEL
    /// models. The single entry point every analysis uses, so a model card is
    /// honoured uniformly regardless of which analysis builds the registry.
    pub fn register_builtin_models(&mut self, cards: &[ModelCard]) {
        self.register_builtin_diodes(cards);
        self.register_builtin_mosfets(cards);
        self.register_builtin_bjts(cards);
        self.register_photonic_models(cards);
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
        // Passive / simple natives (Default-constructible).
        self.register_default::<NativeWaveguide>("fc_waveguide");
        self.register_default::<NativeDirectionalCoupler>("fc_dcoupler");
        self.register_default::<NativeSplitter>("fc_splitter");
        self.register_default::<NativeGratingCoupler>("fc_grating_coupler");
        self.register_default::<NativePhotodetector>("fc_photodetector");
        self.register_default::<NativeMzm>("fc_mzm");
        self.register_default::<NativeCirculator>("fc_circulator");
        self.register_default::<NativeCwLaser>("fc_cw_laser");
        self.register_default::<NativeMux>("fc_mux");
        self.register_default::<NativeDemux>("fc_demux");
        self.register_default::<NativeOptical2x2>("fc_optical_2x2");
        self.register_default::<NativeAwgr>("fc_awgr");

        // Active phase shifters (built by constructors → ActiveOpticalDevice).
        self.register_ctor("fc_thermal_ps", thermal_phase_shifter);
        self.register_ctor("fc_thermal_ps_rc", thermal_rc_phase_shifter);
        self.register_ctor("fc_pn_ps", pn_phase_shifter);
        self.register_ctor("fc_pn_ps_cap", pn_phase_shifter_cap);
        self.register_ctor("fc_pn_th_ps", pn_thermal_phase_shifter);
        self.register_ctor("fc_pn_th_ps_cap", pn_thermal_phase_shifter_cap);
        self.register_ctor("fc_pn_ps_inj", pn_phase_shifter_inj);
        self.register_ctor("fc_pn_th_ps_inj", pn_thermal_phase_shifter_inj);
        self.register_ctor("fc_pn_ps_full", pn_phase_shifter_full);
        self.register_ctor("fc_pn_th_ps_full", pn_thermal_phase_shifter_full);
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
//
// PDK device-name aliasing with parameter-name translation is implemented in
// `register_alias`: it wraps the target factory in one that renames the
// `ParamSet` keys before `build`, so the translation covers construction-time
// params as well as `set_real_param`. (The old `AliasedDevice` Device-wrapper —
// which could only translate the post-construction `set_real_param` path — is
// gone.)

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
        let expected = (-46.0517019_f64 * 100e-6 / 2.0).exp();
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
