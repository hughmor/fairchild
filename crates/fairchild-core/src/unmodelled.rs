//! The one list of parameters a model **accepts and does not model**.
//!
//! A parameter can fail three ways and from the outside they look identical:
//! the netlist runs and you get a number.  It can be unknown (a typo), it can be
//! modelled, or it can be matched by the model, stored, and never read.  The
//! third is the dangerous one — the deck *names* the effect it wants, the run
//! completes, and the answer was computed without it.  `IKF` is the case that
//! motivated this: it is what somebody adds when beta looked too flat, and they
//! got the same flat beta back with nothing on stderr.
//!
//! Each model used to carry its own list inside its parameter `match`, silently.
//! Two lists in two places are two chances to disagree, which is exactly how the
//! diode came to warn about `BV` while the BJT said nothing about nine
//! parameters.  So the lists live here, next to the diagnostic that reads them,
//! and each entry carries **what the deck loses** — `IKF ignored: high-injection
//! roll-off is not modelled, so the forward current keeps its exponential slope
//! past the knee` tells a user whether it matters to them, which `unknown
//! parameter IKF` does not.
//!
//! `docs/model_status.md` is the human-facing version of these tables, and
//! `unmodelled_matches_the_audit` in `tests/model_parameter_diagnostics.rs`
//! fails if they disagree.

/// `(parameter, what is missing because it is not modelled)`.
pub type Unmodelled = &'static [(&'static str, &'static str)];

/// Diode (`.model … D`).
pub const DIODE: Unmodelled = &[
    ("eg", "IS does not vary with temperature"),
    ("xti", "IS does not vary with temperature"),
    (
        "kf",
        "flicker (1/f) noise is not modelled, in .noise or in .tran",
    ),
    (
        "af",
        "flicker (1/f) noise is not modelled, in .noise or in .tran",
    ),
    (
        "isr",
        "the recombination current is not modelled, so the low-current \
         ideality stays at N",
    ),
    (
        "nr",
        "the recombination current is not modelled (needs ISR)",
    ),
    (
        "ikf",
        "high-injection roll-off is not modelled, so the forward current keeps \
         its exponential slope past the knee",
    ),
    (
        "tnom",
        "no parameter is re-referenced to the extraction temperature",
    ),
    ("trs1", "RS does not vary with temperature"),
    ("trs2", "RS does not vary with temperature"),
    (
        "cta",
        "the junction capacitance does not vary with temperature",
    ),
    ("vpt", "punch-through is not modelled"),
];

/// BJT (`.model … NPN|PNP`).
pub const BJT: Unmodelled = &[
    (
        "cjs",
        "the collector-substrate junction capacitance is not stamped, so the \
         substrate terminal carries no charge",
    ),
    ("vjs", "the substrate junction is not stamped (needs CJS)"),
    ("mjs", "the substrate junction is not stamped (needs CJS)"),
    ("fcs", "the substrate junction is not stamped (needs CJS)"),
    (
        "xcjc",
        "all of CJC sits outside the base resistance, so RB does not see its \
         share of the collector charge",
    ),
    (
        "rbm",
        "the base resistance is constant: it does not fall towards RBM at high \
         current",
    ),
    ("irb", "the base resistance is constant (needs RBM)"),
    (
        "xtf",
        "the forward transit time is constant: TF does not rise with bias, so \
         fT is flat in current",
    ),
    ("vtf", "the forward transit time is constant (needs XTF)"),
    ("itf", "the forward transit time is constant (needs XTF)"),
    ("ptf", "excess phase is not modelled"),
    ("xtb", "the betas do not vary with temperature"),
    ("eg", "IS does not vary with temperature"),
    ("xti", "IS does not vary with temperature"),
    (
        "kf",
        "flicker (1/f) noise is not modelled, in .noise or in .tran",
    ),
    (
        "af",
        "flicker (1/f) noise is not modelled, in .noise or in .tran",
    ),
    (
        "tnom",
        "no parameter is re-referenced to the extraction temperature",
    ),
];

/// MOSFET (`.model … NMOS|PMOS`), Level 1.
pub const MOSFET: Unmodelled = &[
    (
        "is",
        "the bulk-source and bulk-drain diodes are not stamped, so a forward-biased \
         bulk conducts nothing",
    ),
    ("js", "the bulk junction diodes are not stamped"),
    ("rd", "the drain ohmic series resistance is not stamped"),
    ("rs", "the source ohmic series resistance is not stamped"),
    ("rsh", "the ohmic series resistances are not stamped"),
    ("nsub", "nothing is derived from the substrate doping"),
    (
        "nss",
        "the surface state density does not shift the threshold",
    ),
    (
        "nfs",
        "there is no subthreshold conduction: below VTO the channel current is \
         exactly zero",
    ),
    (
        "tpg",
        "the gate material does not shift the flat-band voltage",
    ),
    (
        "uo",
        "the mobility is not used to derive KP — give KP directly, or the \
         default 2e-5 applies",
    ),
    ("ucrit", "mobility degradation with field is not modelled"),
    ("uexp", "mobility degradation with field is not modelled"),
    (
        "utra",
        "transverse-field mobility degradation is not modelled",
    ),
    (
        "vmax",
        "carrier velocity saturation is not modelled, so a short channel keeps \
         the long-channel saturation current",
    ),
    (
        "xj",
        "the short-channel threshold correction is not modelled",
    ),
    (
        "ld",
        "lateral diffusion is not subtracted, so L is the drawn length rather \
         than the effective one",
    ),
    (
        "delta",
        "the narrow-width threshold correction is not modelled",
    ),
    (
        "theta",
        "mobility degradation with gate field is not modelled",
    ),
    ("eta", "static feedback on the threshold is not modelled"),
    ("kappa", "the saturation-field factor is not modelled"),
    (
        "kf",
        "flicker (1/f) noise is not modelled, in .noise or in .tran",
    ),
    (
        "af",
        "flicker (1/f) noise is not modelled, in .noise or in .tran",
    ),
    (
        "tnom",
        "no parameter is re-referenced to the extraction temperature",
    ),
    ("php", "the sidewall junction uses PB as its potential"),
];

/// Voltage/current switch (`.model … SW|CSW`).
pub const SWITCH: Unmodelled = &[];

/// Whether `key` is on `table` — the test a model's parameter `match` uses so
/// it accepts exactly what this file says it accepts, and nothing else.
pub fn is_listed(table: Unmodelled, key: &str) -> bool {
    table.iter().any(|(k, _)| k.eq_ignore_ascii_case(key))
}

/// The diagnostics for one `.model` card: one line per parameter that is
/// accepted and does nothing, naming what would have happened if it were
/// modelled.
///
/// Returns the lines rather than printing them, so a test can assert on the
/// classification without capturing stderr.  The caller emits them — once per
/// card, never once per instance: a netlist with 500 transistors on one card has
/// one thing wrong with it, not 500.
pub fn report(table: Unmodelled, params: &[(String, f64)]) -> Vec<String> {
    params
        .iter()
        .filter_map(|(k, _)| {
            table
                .iter()
                .find(|(name, _)| name.eq_ignore_ascii_case(k))
                .map(|(name, effect)| format!("{} ignored: {effect}", name.to_uppercase()))
        })
        .collect()
}
