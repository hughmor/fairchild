use crate::warn_user;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use fairchild_parser::{ArityOracle, ArityQuery, BundleArity, Expr, ModelCard};

use crate::device::{Device, NodeId, SimContext};
use crate::models::{
    expr_phase_shifter, pn_phase_shifter, pn_phase_shifter_cap, pn_phase_shifter_full,
    pn_phase_shifter_inj, pn_thermal_phase_shifter, pn_thermal_phase_shifter_cap,
    pn_thermal_phase_shifter_full, pn_thermal_phase_shifter_inj, thermal_phase_shifter,
    thermal_rc_phase_shifter, ActiveOpticalDevice, GummelPoonBjt, Mosfet1, NativeAwgr,
    NativeCirculator, NativeCwLaser, NativeDemux, NativeDirectionalCoupler, NativeDrivenLaser,
    NativeFacet, NativeGratingCoupler, NativeMux, NativeMzm, NativeOptical2x2, NativePhotodetector,
    NativeSplitter, NativeWaveguide, ShockleyDiode, SpectrumTable,
};

/// The set of instance parameters from an `X…` element line, threaded into a
/// [`ModelFactory`] call. Keys are case-insensitive; the set tracks which
/// keys a device actually consumed so unrecognised (typo'd) params can be
/// reported, and can rename keys for the PDK-alias path.
///
/// Created per device instantiation and passed by reference; it is *not* stored
/// in the registry, so its interior mutability (consumed-tracking) needs no
/// `Sync`.
pub struct ParamSet {
    /// Lower-cased (key, value) pairs in netlist order.
    params: Vec<(String, f64)>,
    /// Where each key came from: `None` for the element line, `Some(card)` for a
    /// `.model` card default merged in by [`ParamSet::with_defaults`].
    ///
    /// Kept so an unrecognised parameter is reported against the thing that
    /// actually wrote it. Without it, merging card defaults into the set would
    /// turn `.model 'nch': unknown parameter 'tox'` into a complaint about the
    /// element line, which is not where the reader should look.
    origin: Vec<Option<Arc<str>>>,
    consumed: std::cell::RefCell<Vec<bool>>,
}

/// A parameter no device recognised, and who wrote it.
pub struct UnknownParam {
    pub key: String,
    /// The `.model` card that supplied it, or `None` for the element line.
    pub card: Option<String>,
}

impl ParamSet {
    /// Build from raw instance params (keys are lower-cased).
    pub fn new(params: &[(String, f64)]) -> Self {
        let params: Vec<(String, f64)> =
            params.iter().map(|(k, v)| (k.to_lowercase(), *v)).collect();
        let n = params.len();
        ParamSet {
            params,
            origin: vec![None; n],
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

    /// Whether `name` is present, **without** marking it consumed.
    ///
    /// For deciding whether a model-card default is overridden by an instance
    /// param: the instance param is applied by the target factory, so counting
    /// it consumed here as well would be a lie either way round.
    pub fn contains(&self, name: &str) -> bool {
        let nl = name.to_lowercase();
        self.params.iter().any(|(k, _)| *k == nl)
    }

    /// Keys not consumed by `get`/`apply` — i.e. params the device did not
    /// recognise (likely typos), each with the card it came from if it was not
    /// on the element line. Used to warn the user.
    pub fn unconsumed(&self) -> Vec<UnknownParam> {
        let consumed = self.consumed.borrow();
        self.params
            .iter()
            .zip(self.origin.iter())
            .zip(consumed.iter())
            .filter(|(_, c)| !**c)
            .map(|(((k, _), o), _)| UnknownParam {
                key: k.clone(),
                card: o.as_ref().map(|c| c.to_string()),
            })
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

    /// This set plus a `.model` card's parameters, for the keys the element line
    /// did not give. The element line always wins.
    ///
    /// The card's parameters are *merged into the set* rather than written onto
    /// the device afterwards, and that ordering is the point: a device validates
    /// itself once, at the end of construction, and it can only do that if every
    /// parameter it will ever be given is already in hand. Patching defaults in
    /// after the factory returned meant a card could complete a configuration
    /// that had already been judged incomplete — or break one that had been
    /// judged fine — with nothing to catch either.
    pub fn with_defaults(&self, defaults: &[(String, f64)], card: &str) -> ParamSet {
        let card: Arc<str> = Arc::from(card);
        let mut params = self.params.clone();
        let mut origin = self.origin.clone();
        for (k, v) in defaults {
            let k = k.to_lowercase();
            if params.iter().any(|(p, _)| *p == k) {
                continue; // element line wins
            }
            params.push((k, *v));
            origin.push(Some(Arc::clone(&card)));
        }
        let n = params.len();
        ParamSet {
            params,
            origin,
            consumed: std::cell::RefCell::new(vec![false; n]),
        }
    }
}

/// Which unknown-parameter warnings an element should emit, given what earlier
/// elements in the same build already reported.
///
/// Two different things are being reported, and they deserve different
/// frequencies:
///
/// - **An element-line parameter** is reported every time. Each element line is
///   a separate thing someone wrote, so two elements with the same typo are two
///   typos and both want naming.
/// - **A `.model` card parameter** is reported once, however many elements name
///   the card. A card is one line: a netlist with 500 diodes on one card has one
///   thing wrong with it, not 500. `register_builtin_diodes` has always followed
///   that rule for the parameters it validates at registration; the card path
///   did not, and printed the same line once per instance.
///
/// `seen` is the dedup memory, keyed on (card, parameter). Scoped to one
/// `build_devices` call rather than to the process: a warning is cosmetic, but a
/// library that silently stopped mentioning something because an earlier
/// simulation had mentioned it would be worse than a repeated line.
///
/// Returns the messages rather than printing them, so the rule above is
/// testable — `warn_user!` writes to stderr and a test cannot see it.
pub fn unknown_param_reports(
    ps: &ParamSet,
    elem: &str,
    model_name: &str,
    seen: &mut HashSet<(String, String)>,
) -> Vec<String> {
    let mut out = Vec::new();
    for p in ps.unconsumed() {
        match &p.card {
            Some(card) => {
                // Deliberately does not name an element: the card is what is
                // wrong, and naming one arbitrary instance out of 500 would
                // point the reader at a line that is fine.
                if seen.insert((card.clone(), p.key.clone())) {
                    // Not "…by '{model_name}'": for a card instance the two
                    // are the same name, and the message read
                    // ".model 'mirror' … not honoured by 'mirror'".
                    out.push(format!(
                        ".model '{card}': parameter '{}' is not honoured by this \
                         model and was ignored",
                        p.key
                    ));
                }
            }
            None => out.push(format!(
                "{elem} ('{model_name}'): instance parameter '{}' is not honoured \
                 by this model and was dropped",
                p.key
            )),
        }
    }
    out
}

/// A card's parameters without `LEVEL`, which selects the device rather than
/// configuring it and is consumed by the registrar.
fn without_level(params: &[(String, f64)]) -> Vec<(String, f64)> {
    params
        .iter()
        .filter(|(k, _)| !k.eq_ignore_ascii_case("level"))
        .cloned()
        .collect()
}

/// Set up, apply every parameter, then let the device judge itself.
///
/// The one construction order in this file. It exists because the order is
/// load-bearing and used to be spelled out at each of nine registration sites:
/// `setup_model`, `setup_instance`, parameters, validation. Get it wrong in one
/// place and that device alone sees defaults where the deck gave values, which
/// is the failure mode #31 describes — plausible numbers, no diagnostic.
fn finish<D: Device + 'static>(
    mut d: D,
    terminals: &[NodeId],
    params: &ParamSet,
    ctx: &SimContext,
) -> Result<Box<dyn Device>, String> {
    d.setup_model(ctx);
    d.setup_instance(terminals, ctx);
    apply_and_validate(d, params)
}

/// [`finish`] from `setup_instance` onwards, for the one device that has work to
/// do in between: `fc_awgr` reads its spectrum file only once `setup_instance`
/// has told it how many ports it has.
fn apply_and_validate<D: Device + 'static>(
    mut d: D,
    params: &ParamSet,
) -> Result<Box<dyn Device>, String> {
    params.apply(&mut d);
    d.validate()?;
    Ok(Box::new(d))
}

/// A device factory: constructs a device from its terminal nodes, instance
/// [`ParamSet`], and the [`SimContext`]. The returned device is fully set up
/// (`setup_model` + `setup_instance` done), has every parameter applied, and has
/// passed [`Device::validate`]; the caller (`build_devices`) handles extra-node
/// allocation and unconsumed-param warnings.
///
/// **Fallible**, because that is the whole point of doing the work here: a
/// device that cannot be built from the parameters it was given says so, once,
/// with the element named, instead of asserting from inside `eval` halfway
/// through a solve. The `Err` describes the device's own problem; the caller
/// adds the element and model name. Stored behind `Arc` so the alias mechanism (B6) can clone a target
/// factory into a wrapper that performs parameter-name translation. (A named
/// `ModelFactory` trait could wrap this later for a dlopen plugin ABI; the
/// closure form covers every in-tree need — OSDI captures `Arc<library>`, the
/// LEVEL/expr factories capture their tables/expressions.)
pub type ModelFactory = dyn Fn(&[NodeId], &ParamSet, &SimContext) -> Result<Box<dyn Device>, String>
    + Send
    + Sync
    + 'static;
type Factory = Arc<ModelFactory>;

/// Name every model-card parameter the model accepts and does not model.
///
/// Once per card, and it says what the parameter *would* have done: `IKF
/// ignored: high-injection roll-off is not modelled` tells the user whether it
/// matters to them, which `unknown parameter IKF` does not. The classification
/// itself is [`crate::unmodelled`]'s — this only prints it.
fn warn_unmodelled(
    kind: &str,
    name: &str,
    table: crate::unmodelled::Unmodelled,
    params: &[(String, f64)],
) {
    for line in crate::unmodelled::report(table, params) {
        warn_user!("{kind} model '{name}': {line}");
    }
}

/// Name every instance parameter a device could not honour.
///
/// One place, for the same reason `crate::unmodelled` is one place: a parameter
/// that parses and changes nothing is the failure these warnings exist for, and
/// fixing `AREA` while leaving three others silent only makes it rarer and
/// harder to find. Per instance rather than per card — an instance parameter is
/// written per instance, so there is nothing to collapse.
fn warn_dropped_instance_params(element: &str, model_name: &str, unknown: Vec<String>) {
    for key in unknown {
        warn_user!(
            "{element} ('{model_name}'): instance parameter '{key}' is not \
             honoured by this model and was dropped"
        );
    }
}

/// Maps model names to device factory closures.
///
/// Factories receive the MNA node mapping and SimContext, and return a fully
/// initialised (setup_model + setup_instance already called) boxed Device.
///
/// Built-in models are registered via `register_builtin_diodes` /
/// `register_builtin_mosfets`. External models (e.g. OSDI) register themselves
/// by capturing an Arc to their library, keeping it alive for the device's lifetime.
/// How a registered model's WDM dispatch is decided.
///
/// The parser used to hold this as a hardcoded list of `fc_*` names, which
/// could only ever match a name written literally on an X-line — so no
/// `.model`-card-named instance was ever found there, and the tier names drifted
/// out of it twice.  Recording it next to the registration keeps one place
/// interpreting the concept (#52).
#[derive(Clone, Copy, Debug)]
pub enum ArityDecl {
    /// Declared outright, for a device whose channel count is derived from
    /// `terminals.len()` at setup — every native photonic model.
    Fixed(BundleArity),
    /// A fixed terminal count, from an OSDI descriptor.  The instance is placed
    /// by shape: if flattening the referenced bundles hits this count the model
    /// takes the whole bus, and if one channel each hits it, it does not.  Same
    /// rule a `.subckt` instance already follows.
    Terminals(usize),
    /// A model written against the bundle-port dialect, which has no single
    /// terminal count: it is generated for whatever width the deck asks for.
    /// Any shape `scalars + per_channel·N` fits, for `N >= 1`.
    Bundle {
        /// Ports that are not part of a bundle.
        scalars: usize,
        /// Terminals every extra channel adds: `bundles × wires_per_channel`.
        per_channel: usize,
    },
}

pub struct DeviceRegistry {
    factories: HashMap<String, Factory>,
    /// WDM dispatch per model name, consulted through the `ArityOracle` impl.
    arities: HashMap<String, ArityDecl>,
    /// Numeric `.model` card parameters, by card name.
    ///
    /// Held here rather than captured in the card's factory closure so the
    /// *caller* can merge them into the one `ParamSet` it owns. A closure that
    /// merges them itself has to build a second set, and then consumption
    /// tracking happens on a set nobody else can see: the card's unrecognised
    /// parameters get reported by nobody, and an element-line parameter the
    /// device *did* consume is reported as dropped. Both of those were live for
    /// exactly as long as it took to write a test for the warning.
    card_defaults: HashMap<String, Arc<Vec<(String, f64)>>>,
    /// MOSFET model cards stored for W/L instance-param injection in build_devices.
    /// One entry per model *name*, which for a binned PDK card is the name
    /// without its `.N` suffix. See [`crate::binning`] for the selection rule.
    pub(crate) mosfet_cards: HashMap<String, (bool, crate::binning::BinGroup)>,
    /// BJT model cards: model_name → (is_pnp, params).
    pub(crate) bjt_cards: HashMap<String, (bool, Vec<(String, f64)>)>,
    /// Switch model cards: model_name → (is_current_controlled, params).
    /// Like the MOSFET/BJT maps, these are consumed by `build_devices` rather
    /// than by a factory closure: a switch needs its instance's `ON`/`OFF`
    /// keyword, and a `W` needs a controlling branch row the factory has no
    /// way to resolve.
    pub(crate) switch_cards: HashMap<String, (bool, Vec<(String, f64)>)>,
}

/// Kinds that read a card's expression params as device constitutive maps.
///
/// Every other kind reads numbers only, so an expression that reached it would be
/// dropped without a word: `.model nm NMOS (VTO={vt})` used to arrive here as an
/// expression param and the MOSFET path never looks at those, so the threshold
/// silently defaulted. The parser now evaluates `{…}` and `'…'` before the card
/// gets here, which leaves this as the backstop for whatever it could not.
const EXPR_PARAM_KINDS: &[&str] = &["fc_phase_shifter_expr", "fc_awgr"];

/// Warn about card values that nothing downstream can read.
fn warn_unusable_expr_params(cards: &[ModelCard]) {
    for card in cards {
        if card.expr_params.is_empty() {
            continue;
        }
        let kind = card.kind.to_lowercase();
        if EXPR_PARAM_KINDS.contains(&kind.as_str()) {
            continue;
        }
        let keys: Vec<&str> = card.expr_params.iter().map(|(k, _)| k.as_str()).collect();
        warn_user!(
            ".model '{}' ({kind}): parameter(s) {} are expressions, and a {kind} card \
             takes numbers — they are being ignored, not evaluated. A value that \
             depends on a .param must be written {{…}} or '…' so it is evaluated at \
             parse time; \"…\" is reserved for a device constitutive map over the \
             device's own bias, which only {} accept",
            card.name,
            keys.iter()
                .map(|k| format!("'{k}'"))
                .collect::<Vec<_>>()
                .join(", "),
            EXPR_PARAM_KINDS.join(" / ")
        );
    }
}

impl DeviceRegistry {
    pub fn new() -> Self {
        let mut reg = Self {
            factories: HashMap::new(),
            arities: HashMap::new(),
            card_defaults: HashMap::new(),
            mosfet_cards: HashMap::new(),
            bjt_cards: HashMap::new(),
            switch_cards: HashMap::new(),
        };
        // Native photonic passives are always available — no .model card or
        // .osdi import required to instantiate `fc_waveguide`, `fc_dcoupler`,
        // `fc_splitter`.
        reg.register_native_photonics();
        reg
    }

    /// Every model name this registry can build.
    ///
    /// Exists so a test can check the registry against the parser's
    /// `bundle_arity_for`. The parser cannot ask the registry what a model is
    /// (`fairchild-core` depends on `fairchild-parser`, so the dependency runs
    /// the wrong way) and the arity list is therefore a hand-maintained second
    /// list of the same facts. A test can see both, so the disagreement is
    /// catchable even though the dispatch cannot be unified.
    pub fn registered_names(&self) -> impl Iterator<Item = &str> {
        self.factories.keys().map(String::as_str)
    }

    /// Record how `name` handles WDM bundles.  Call it beside the registration
    /// of any model that can carry an optical or electrical bundle port.
    pub fn declare_arity(&mut self, name: impl Into<String>, decl: ArityDecl) {
        self.arities.insert(name.into().to_lowercase(), decl);
    }

    /// Record the numeric parameters of a `.model` card, for the caller to
    /// merge into the instance `ParamSet`. See [`Self::card_defaults`].
    fn declare_card_defaults(&mut self, card: &str, params: Vec<(String, f64)>) {
        if params.is_empty() {
            return;
        }
        self.card_defaults
            .insert(card.to_lowercase(), Arc::new(params));
    }

    /// The `.model` card parameters recorded for `name`, if it is a card.
    pub fn card_defaults(&self, name: &str) -> Option<&[(String, f64)]> {
        self.card_defaults
            .get(&name.to_lowercase())
            .map(|v| v.as_slice())
    }

    /// The instance parameters for an element, with its `.model` card's
    /// parameters merged in underneath. The element line wins.
    ///
    /// The single place a device's full parameter set is assembled, which is what
    /// makes one `unconsumed()` report cover both halves.
    pub fn params_for(&self, model_name: &str, instance: &[(String, f64)]) -> ParamSet {
        let ps = ParamSet::new(instance);
        match self.card_defaults(model_name) {
            Some(defaults) => ps.with_defaults(defaults, model_name),
            None => ps,
        }
    }

    /// Declared arity for `name`, if any.
    pub fn arity_decl(&self, name: &str) -> Option<ArityDecl> {
        self.arities.get(&name.to_lowercase()).copied()
    }

    /// Register a factory for `name`. Overwrites any previous entry.
    ///
    /// The name is folded to lower case, because that is what every lookup here
    /// asks with: a deck's model name arrives from the parser already
    /// lowercased, SPICE being case-insensitive. A registrar that hands over a
    /// name with any upper case in it — a Verilog-A module, which preserves the
    /// case the author wrote — would otherwise put an entry in the map that no
    /// deck can name. `PSP102VA`, `DIODE_CMC`, `JUNCAP200` and `hicumL2va` all
    /// loaded successfully and then reported `unknown model` at the element
    /// that used them. `declare_arity` had always folded; only this half did
    /// not, which is the whole bug.
    pub fn register(
        &mut self,
        name: impl Into<String>,
        factory: impl Fn(&[NodeId], &ParamSet, &SimContext) -> Result<Box<dyn Device>, String>
            + Send
            + Sync
            + 'static,
    ) {
        self.factories
            .insert(name.into().to_lowercase(), Arc::new(factory));
    }

    /// Register a default-constructible device: `T::default()` + the standard
    /// setup + instance-param application. Collapses the boilerplate for every
    /// passive/simple native device.
    pub fn register_default<T: Device + Default + 'static>(&mut self, name: impl Into<String>) {
        self.register(name, |terminals, params: &ParamSet, ctx| {
            finish(T::default(), terminals, params, ctx)
        });
    }

    /// Register a device built by a constructor function (the active-photonic
    /// `ActiveOpticalDevice` builders, which are not `Default`). Same setup +
    /// instance-param flow as [`register_default`](Self::register_default).
    fn register_ctor(&mut self, name: impl Into<String>, ctor: fn() -> ActiveOpticalDevice) {
        self.register(name, move |terminals, params: &ParamSet, ctx| {
            finish(ctor(), terminals, params, ctx)
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
        let target_factory: Factory = self
            .factories
            .get(&target_name.to_lowercase())
            .cloned()
            .ok_or_else(|| {
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

    /// Register every `.model <card> <kind> (...)` whose `kind` names a factory
    /// that is already in the registry, as an alias of that factory with the
    /// card's parameters as defaults.
    ///
    /// This is the `.model`-card indirection for **OSDI models**, and it is the
    /// idiom every foundry PDK ships:
    ///
    /// ```spice
    /// .osdi  bsim4.osdi
    /// .model nch  bsim4 (tox=3n vth0=0.4 …)
    /// M1 d g s b nch W=1u L=100n
    /// ```
    ///
    /// Call it *after* `register_builtin_models` and after every `.osdi`
    /// library has been registered. Cards whose name is already registered are
    /// skipped, so the native card handlers (which do construction-time work
    /// this cannot — `LEVEL` dispatch, expression parsing) keep ownership of
    /// their cards; this only fills the gap they leave.
    ///
    /// Card params are applied *after* construction and only for keys the
    /// instance line did not give, so an instance param always wins. The
    /// consequence, and the reason this is not the native card path: a device
    /// that reads params at construction time (via [`ParamSet::get`]) will not
    /// see these defaults. OSDI models have no construction-time params —
    /// everything routes through `set_real_param` — so the OSDI case is exact.
    pub fn register_loaded_model_cards(&mut self, cards: &[ModelCard]) {
        for card in cards {
            let card_name = card.name.to_lowercase();
            let kind = card.kind.to_lowercase();
            if card_name == kind || self.factories.contains_key(&card_name) {
                continue;
            }
            let Some(target) = self.factories.get(&kind).cloned() else {
                continue;
            };
            self.declare_card_defaults(&card.name, card.params.clone());
            // A plain alias. The card's parameters are merged by `params_for`
            // before this is ever called, so there is nothing to do here — which
            // is the point: one set, one consumption record, one report.
            self.register(card_name, move |terminals, params: &ParamSet, ctx| {
                target(terminals, params, ctx)
            });
        }
    }

    /// Populate the registry from `.model D` cards using the built-in Shockley diode.
    pub fn register_builtin_diodes(&mut self, cards: &[ModelCard]) {
        for card in cards {
            if !card.kind.to_lowercase().starts_with('d') {
                continue;
            }
            let params: Vec<(String, f64)> = card.params.clone();
            // Warn once per model card, however many instances it has: a netlist
            // with 500 diodes on one card has one thing wrong with it, not 500.
            let (_, unknown) = ShockleyDiode::from_params(&params);
            if !unknown.is_empty() {
                warn_user!(
                    "diode model '{}': unknown parameter(s) {}",
                    card.name,
                    unknown.join(", ")
                );
            }
            warn_unmodelled("diode", &card.name, crate::unmodelled::DIODE, &params);
            self.register(card.name.clone(), move |terminals, ps: &ParamSet, ctx| {
                let (dev, _) = ShockleyDiode::from_params(&params);
                finish(dev, terminals, ps, ctx)
            });
            // `build_devices` names whatever `apply` left unconsumed.
        }
    }

    /// Record `.model … SW` / `.model … CSW` cards for `build_devices`.
    ///
    /// Unlike diodes there is no factory closure: the element line carries an
    /// `ON`/`OFF` keyword and, for `W`, a controlling voltage-source name that
    /// only the builder can turn into an MNA row.
    pub fn register_builtin_switches(&mut self, cards: &[ModelCard]) {
        for card in cards {
            let is_current = match card.kind.to_lowercase().as_str() {
                "sw" | "vswitch" => false,
                "csw" | "iswitch" => true,
                _ => continue,
            };
            // Warn once per card, matching the diode/MOSFET/BJT convention.
            match crate::models::Switch::from_model_params(is_current, &card.params, false) {
                Ok((_, unknown)) if !unknown.is_empty() => warn_user!(
                    "switch model '{}': unknown parameter(s) {}",
                    card.name,
                    unknown.join(", ")
                ),
                // A bad RON/ROFF is reported when the instance is built, where
                // there is an error path to return it on.
                Ok(_) | Err(_) => {}
            }
            warn_unmodelled(
                "switch",
                &card.name,
                crate::unmodelled::SWITCH,
                &card.params,
            );
            self.switch_cards
                .insert(card.name.clone(), (is_current, card.params.clone()));
        }
    }

    /// Populate the registry from `.model NMOS` / `.model PMOS` cards.
    ///
    /// MOSFET factories do not accept instance W/L here; those are injected by
    /// `build_devices` at instantiation time using the stored `mosfet_cards` map.
    pub fn register_builtin_mosfets(&mut self, cards: &[ModelCard]) {
        for card in cards {
            let kind = card.kind.to_lowercase();

            // A card-named instance is looked up by the CARD's name, never by
            // its kind, so the card must inherit its kind's WDM dispatch.
            // Without this no card-based photonic device could carry a bundle
            // at all: `.model awg2 fc_awgr` was refused on the very buses an
            // AWG router exists to route, and a `.model … fc_pn_ps LEVEL=4`
            // could not take a WDM bus that bare `fc_pn_ps_full` accepted (#52).
            if let Some(decl @ ArityDecl::Fixed(_)) = self.arity_decl(&kind) {
                self.declare_arity(card.name.clone(), decl);
            }
            let is_pmos = match kind.as_str() {
                "nmos" => false,
                "pmos" => true,
                _ => continue,
            };
            // `LEVEL` gets its own warning. Folding it into the generic
            // unimplemented-params list reads as one defaulted coefficient, when
            // what actually happened is that a different device model was asked
            // for and Level 1 was simulated instead.
            if let Some((_, level)) = card
                .params
                .iter()
                .find(|(k, _)| k.eq_ignore_ascii_case("level"))
            {
                if (level - 1.0).abs() > 1e-9 {
                    warn_user!(
                        "MOSFET model '{}' asks for LEVEL={} — fairchild \
                         implements Level 1 (Shichman-Hodges) only, and is simulating \
                         this card as Level 1. Currents and capacitances will differ \
                         from the intended model, not merely in the unset parameters.",
                        card.name,
                        level
                    );
                }
            }
            // Warn once per model card about unrecognised params. The four
            // geometry selectors are not model parameters — they choose a bin —
            // so a binned card would otherwise warn about all four.
            let (_, unknown) = Mosfet1::from_model_params(is_pmos, &card.params);
            let unknown: Vec<_> = unknown
                .into_iter()
                .filter(|k| !k.eq_ignore_ascii_case("level") && !crate::binning::IS_SELECTOR(k))
                .collect();
            if !unknown.is_empty() {
                warn_user!(
                    "MOSFET model '{}': unknown parameter(s) {}",
                    card.name,
                    unknown.join(", ")
                );
            }
            warn_unmodelled(
                "MOSFET",
                &card.name,
                crate::unmodelled::MOSFET,
                &card.params,
            );
            // `nch.1` is bin 1 of `nch`, and the element line asks for `nch`.
            // An unbinned card is a group of one with an unbounded window, so
            // both take the same path out of `build_mosfet`.
            let (base, bin) = crate::binning::classify(&card.name, &card.params);
            self.mosfet_cards
                .entry(base)
                .or_insert_with(|| (is_pmos, crate::binning::BinGroup::default()))
                .1
                .push(bin);
        }
        // Once per name, after every card is in: a set whose windows overlap in
        // their interiors is malformed, and there the choice of bin really does
        // change the answer.
        for (name, (_, group)) in &self.mosfet_cards {
            for (a, b) in group.interior_overlaps() {
                warn_user!(
                    "MOSFET model '{name}': bins {a} and {b} overlap in their \
                     interiors, not merely at an edge. One geometry matches both, \
                     and the tighter window wins — check LMIN/LMAX/WMIN/WMAX on \
                     those two cards."
                );
            }
        }
    }

    /// Populate the registry from `.model NPN` / `.model PNP` cards.
    pub fn register_builtin_bjts(&mut self, cards: &[ModelCard]) {
        for card in cards {
            let kind = card.kind.to_lowercase();

            // A card-named instance is looked up by the CARD's name, never by
            // its kind, so the card must inherit its kind's WDM dispatch.
            // Without this no card-based photonic device could carry a bundle
            // at all: `.model awg2 fc_awgr` was refused on the very buses an
            // AWG router exists to route, and a `.model … fc_pn_ps LEVEL=4`
            // could not take a WDM bus that bare `fc_pn_ps_full` accepted (#52).
            if let Some(decl @ ArityDecl::Fixed(_)) = self.arity_decl(&kind) {
                self.declare_arity(card.name.clone(), decl);
            }
            let is_pnp = match kind.as_str() {
                "npn" => false,
                "pnp" => true,
                _ => continue,
            };
            let (_, unknown) = GummelPoonBjt::from_model_params(is_pnp, &card.params);
            if !unknown.is_empty() {
                warn_user!(
                    "BJT model '{}': unknown parameter(s) {}",
                    card.name,
                    unknown.join(", ")
                );
            }
            warn_unmodelled("BJT", &card.name, crate::unmodelled::BJT, &card.params);
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

            // A card-named instance is looked up by the CARD's name, never by
            // its kind, so the card must inherit its kind's WDM dispatch.
            // Without this no card-based photonic device could carry a bundle
            // at all: `.model awg2 fc_awgr` was refused on the very buses an
            // AWG router exists to route, and a `.model … fc_pn_ps LEVEL=4`
            // could not take a WDM bus that bare `fc_pn_ps_full` accepted (#52).
            if let Some(decl @ ArityDecl::Fixed(_)) = self.arity_decl(&kind) {
                self.declare_arity(card.name.clone(), decl);
            }

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
                self.declare_card_defaults(&card.name, card.params.clone());
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
                            Err(e) => warn_user!(
                                "fc_awgr model '{name}' could not load sfile=\"{p}\" \
                                 ({e}); falling back to the analytic response"
                            ),
                        }
                    }
                    apply_and_validate(d, ps)
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
                                warn_user!(
                                    "photonic model '{}' {name} expression \
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
                self.declare_card_defaults(&card.name, without_level(&card.params));
                self.register(card.name.clone(), move |terminals, ps: &ParamSet, ctx| {
                    let d = expr_phase_shifter(dneff.clone(), dalpha.clone(), g_pn);
                    finish(d, terminals, ps, ctx)
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
                    warn_user!(
                        "photonic model '{}' has unsupported {k} LEVEL={bad}; \
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
            // `LEVEL` selected the constructor above and is not a device
            // parameter, so it does not go in the table — otherwise every
            // instance would report it as unrecognised.
            self.declare_card_defaults(&card.name, without_level(&card.params));
            self.register(card.name.clone(), move |terminals, ps: &ParamSet, ctx| {
                finish(ctor(), terminals, ps, ctx)
            });
        }
    }

    /// Register every `.model`-card-derived built-in family from a netlist's
    /// cards in one call: diodes, MOSFETs, BJTs, and active-photonic LEVEL
    /// models. The single entry point every analysis uses, so a model card is
    /// honoured uniformly regardless of which analysis builds the registry.
    pub fn register_builtin_models(&mut self, cards: &[ModelCard]) {
        warn_unusable_expr_params(cards);
        self.register_builtin_diodes(cards);
        self.register_builtin_mosfets(cards);
        self.register_builtin_bjts(cards);
        self.register_builtin_switches(cards);
        self.register_photonic_models(cards);
    }

    /// Build a `GummelPoonBjt` instance for a `Q` element, injecting the
    /// stored model-card parameters and the element's own. Returns `None` if the
    /// model name is unknown.
    ///
    /// `instance_params` used not to exist here at all: `build_devices`' BJT arm
    /// did not even destructure the element's parameter list, so `AREA` on a `Q`
    /// line reached the netlist and stopped. Anything still unhonoured is named
    /// on stderr rather than dropped — `element` is there to name it.
    pub(crate) fn build_bjt(
        &self,
        model_name: &str,
        element: &str,
        instance_params: &[(String, f64)],
        terminals: &[NodeId],
        ctx: &SimContext,
    ) -> Result<Option<Box<dyn Device>>, crate::SimError> {
        let Some((is_pnp, model_params)) = self.bjt_cards.get(model_name) else {
            return Ok(None);
        };
        let (mut dev, _) = GummelPoonBjt::from_model_params(*is_pnp, model_params);
        warn_dropped_instance_params(
            element,
            model_name,
            dev.set_instance_params(instance_params),
        );
        // After the instance params: AREA scales IS, so `vcrit` — which
        // `setup_model` derives from it — has to see the scaled value.
        dev.setup_model(ctx);
        dev.setup_instance(terminals, ctx);
        dev.validate().map_err(|why| {
            crate::SimError::ParameterError(format!("{element} ('{model_name}'): {why}"))
        })?;
        Ok(Some(Box::new(dev)))
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
    /// - `fc_driven_laser`  — laser whose power follows an electrical input.
    /// - `fc_facet`         — one-port terminator / partial reflector / mirror.
    /// - `fc_thermal_ps`    — thermal phase shifter (Joule heating → φ = π·P/P_pi).
    /// - `fc_pn_ps`         — PN-junction phase shifter (Δn_eff = dn_dv·V).
    pub fn register_native_photonics(&mut self) {
        // WDM dispatch, declared here beside the factories rather than in a
        // separate list in the parser (#52).  Every photonic device is
        // bundle-aware: pure-optical ones run independent per-channel
        // propagation, and ones with electrical state share a single physical
        // interface across all N channels.  WDM is the rule, not the exception.
        //
        // The two lasers are deliberately Scalar — one laser emits one
        // wavelength; combine them with `fc_mux` for a WDM bus.
        for name in [
            "fc_waveguide",
            "fc_dcoupler",
            "fc_splitter",
            "fc_grating_coupler",
            "fc_photodetector",
            "fc_mzm",
            "fc_circulator",
            "fc_facet",
            "fc_optical_2x2",
            "fc_awgr",
            "fc_thermal_ps",
            "fc_thermal_ps_rc",
            "fc_pn_ps",
            "fc_pn_ps_cap",
            "fc_pn_ps_inj",
            "fc_pn_ps_full",
            "fc_pn_th_ps",
            "fc_pn_th_ps_cap",
            "fc_pn_th_ps_inj",
            "fc_pn_th_ps_full",
            // A card kind rather than a factory name — its instances are named
            // after the card. Declared so the inheritance below can find it.
            "fc_phase_shifter_expr",
        ] {
            self.declare_arity(name, ArityDecl::Fixed(BundleArity::Aware));
        }
        for name in ["fc_mux", "fc_demux"] {
            self.declare_arity(name, ArityDecl::Fixed(BundleArity::Bridge));
        }
        for name in ["fc_cw_laser", "fc_driven_laser"] {
            self.declare_arity(name, ArityDecl::Fixed(BundleArity::Scalar));
        }

        // Passive / simple natives (Default-constructible).
        self.register_default::<NativeWaveguide>("fc_waveguide");
        self.register_default::<NativeDirectionalCoupler>("fc_dcoupler");
        self.register_default::<NativeSplitter>("fc_splitter");
        self.register_default::<NativeGratingCoupler>("fc_grating_coupler");
        self.register_default::<NativePhotodetector>("fc_photodetector");
        self.register_default::<NativeMzm>("fc_mzm");
        self.register_default::<NativeCirculator>("fc_circulator");
        self.register_default::<NativeCwLaser>("fc_cw_laser");
        self.register_default::<NativeDrivenLaser>("fc_driven_laser");
        self.register_default::<NativeFacet>("fc_facet");
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
    ///
    /// `Ok(None)` means this name is not a MOSFET card, so the caller falls
    /// through to the generic factory lookup. `Err` means it *is* one and the
    /// instance cannot be built — today only because its geometry falls outside
    /// every bin, which is a hard error rather than a nearest-bin guess.
    pub(crate) fn build_mosfet(
        &self,
        model_name: &str,
        element: &str,
        instance_params: &[(String, f64)],
        terminals: &[NodeId],
        ctx: &SimContext,
    ) -> Result<Option<Box<dyn Device>>, crate::SimError> {
        let Some((is_pmos, group)) = self.mosfet_cards.get(model_name) else {
            return Ok(None);
        };
        // The same defaults the device itself applies, from the same constants,
        // so bin selection and model evaluation cannot disagree about geometry.
        let geom = |key: &str, default: f64| {
            instance_params
                .iter()
                .find(|(k, _)| k.eq_ignore_ascii_case(key))
                .map_or(default, |(_, v)| *v)
        };
        let l = geom("l", crate::models::mosfet1::DEFAULT_L_M);
        let w = geom("w", crate::models::mosfet1::DEFAULT_W_M);
        let model_params = group.select(model_name, l, w)?;
        let (mut dev, _) = Mosfet1::from_model_params(*is_pmos, model_params);
        // The return used to be discarded, so `M1 … banana=3` was accepted in
        // silence on the one device family whose instance params do work.
        warn_dropped_instance_params(
            element,
            model_name,
            dev.set_instance_params(instance_params),
        );
        dev.setup_model(ctx);
        dev.setup_instance(terminals, ctx);
        // `validate` ran only on the factory path (`finish`), so the MOSFET and
        // the BJT — which have builders of their own — never reached it. A hook
        // two device families bypass is a hook that does not exist for them.
        dev.validate().map_err(|why| {
            crate::SimError::ParameterError(format!("{element} ('{model_name}'): {why}"))
        })?;
        Ok(Some(Box::new(dev)))
    }

    /// Look up a factory by model name, case-insensitively — see `register`.
    pub fn get(&self, name: &str) -> Option<&Factory> {
        self.factories.get(&name.to_lowercase())
    }
}

/// The registry answers the parser's WDM dispatch question.
///
/// This is what replaces the hardcoded name list as the authority: the registry
/// is the thing that actually knows what a model name resolves to, including a
/// `.model` card's name, which the parser can never know.  Returning `None`
/// leaves the parser on its static fallback, so a name we have not registered
/// behaves exactly as it did before.
impl ArityOracle for DeviceRegistry {
    fn arity(&self, q: &ArityQuery) -> Option<BundleArity> {
        match self.arity_decl(q.model_name)? {
            ArityDecl::Fixed(a) => Some(a),
            // Placed by shape.  Flattening wins ties: a one-channel bundle makes
            // both counts equal, and there the two dispatches are the same
            // expansion anyway.
            ArityDecl::Terminals(n) => {
                if q.flattened == n {
                    Some(BundleArity::Aware)
                } else if q.single == n {
                    Some(BundleArity::Scalar)
                } else {
                    // Neither shape fits. Say nothing and let the parser report
                    // the mismatch it already reports well.
                    None
                }
            }
            // A generated model serves any width, so the only question is
            // whether the flattened shape is one an integer channel count
            // produces at all.
            ArityDecl::Bundle {
                scalars,
                per_channel,
            } => {
                let rest = q.flattened.checked_sub(scalars)?;
                (per_channel > 0 && rest % per_channel == 0 && rest / per_channel >= 1)
                    .then_some(BundleArity::Aware)
            }
        }
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

    /// A `.model` card's bad parameter is reported once, however many elements
    /// name the card. An element line's is reported every time.
    ///
    /// The frequencies differ because the mistakes differ: a card is one line in
    /// the deck, and 500 diodes on it are 500 symptoms of one typo. Two element
    /// lines with the same typo are two typos.
    ///
    /// Tested on the rule rather than on stderr because `warn_user!` prints and
    /// a test cannot see that — the printing is a one-line wrapper over this.
    #[test]
    fn a_card_parameter_is_reported_once_and_an_instance_one_every_time() {
        let mut seen = HashSet::new();
        let card = ParamSet::empty().with_defaults(&[("bogus".into(), 1.0)], "nch");

        let first = unknown_param_reports(&card, "d1", "dm", &mut seen);
        assert_eq!(first.len(), 1, "{first:?}");
        assert!(
            first[0].contains("nch") && first[0].contains("bogus"),
            "{first:?}"
        );
        // The card name appears once, not twice: for a card instance the model
        // name *is* the card name.
        assert_eq!(first[0].matches("nch").count(), 1, "{first:?}");
        // The card is what is wrong, so the message does not name an element.
        assert!(
            !first[0].contains("d1"),
            "a card report should not name one arbitrary instance: {first:?}"
        );

        // Same card, a different element: silent.
        let again = card_set();
        assert!(
            unknown_param_reports(&again, "d2", "dm", &mut seen).is_empty(),
            "the second element must not repeat the card warning"
        );
        // …and a third, in case the dedup only held for one repeat.
        let third = card_set();
        assert!(unknown_param_reports(&third, "d3", "dm", &mut seen).is_empty());

        // A *different* bad key on the same card is a different mistake.
        let other = ParamSet::empty().with_defaults(&[("alsobogus".into(), 1.0)], "nch");
        assert_eq!(
            unknown_param_reports(&other, "d4", "dm", &mut seen).len(),
            1
        );

        // The same key on a *different* card is also a different mistake.
        let other_card = ParamSet::empty().with_defaults(&[("bogus".into(), 1.0)], "pch");
        assert_eq!(
            unknown_param_reports(&other_card, "d5", "dm", &mut seen).len(),
            1
        );

        // An element-line parameter repeats: two elements, two reports.
        let inst = || ParamSet::new(&[("bogus".into(), 1.0)]);
        assert_eq!(
            unknown_param_reports(&inst(), "x1", "m", &mut seen).len(),
            1
        );
        let second = unknown_param_reports(&inst(), "x2", "m", &mut seen);
        assert_eq!(second.len(), 1, "an element-line typo is named every time");
        assert!(second[0].contains("x2"), "{second:?}");
    }

    fn card_set() -> ParamSet {
        ParamSet::empty().with_defaults(&[("bogus".into(), 1.0)], "nch")
    }

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
             .op\n",
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
