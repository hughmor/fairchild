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
        "fcs",
        "ngspice ignores it too: its substrate junction linearises about zero \
         bias, not about FCS·VJS, and the forward capacitance is bit-identical \
         for FCS of 0.1, 0.5, 0.9 and absent",
    ),
    (
        "rbm",
        "the base resistance is constant: it does not fall towards RBM at high \
         current",
    ),
    ("irb", "the base resistance is constant (needs RBM)"),
    ("ptf", "excess phase is not modelled"),
];

/// MOSFET (`.model … NMOS|PMOS`), Level 1.
pub const MOSFET: Unmodelled = &[
    (
        "rsh",
        "the sheet resistance needs NRD/NRS squares to become a resistance, and \
         those instance parameters are not taken — give RD/RS directly",
    ),
    ("nsub", "nothing is derived from the substrate doping"),
    (
        "nss",
        "the surface state density does not shift the threshold",
    ),
    (
        "nfs",
        "a LEVEL 2/3 parameter: there is no subthreshold conduction at LEVEL 1, \
         so below VTO the channel current is exactly zero. ngspice's LEVEL 1 \
         ignores it too",
    ),
    (
        "tpg",
        "the gate material does not shift the flat-band voltage",
    ),
    (
        "ucrit",
        "a LEVEL 2 parameter: field-dependent mobility is not part of LEVEL 1, \
         and ngspice's LEVEL 1 ignores it too",
    ),
    (
        "uexp",
        "a LEVEL 2 parameter: field-dependent mobility is not part of LEVEL 1",
    ),
    (
        "utra",
        "a LEVEL 2 parameter: transverse-field mobility degradation is not part \
         of LEVEL 1",
    ),
    (
        "vmax",
        "a LEVEL 2/3 parameter: velocity saturation is not part of LEVEL 1, so a \
         short channel keeps the long-channel saturation current. ngspice's \
         LEVEL 1 ignores it too",
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
        "a LEVEL 3 parameter: gate-field mobility degradation is not part of \
         LEVEL 1, and ngspice's LEVEL 1 ignores it too",
    ),
    (
        "eta",
        "a LEVEL 3 parameter: static feedback on the threshold is not part of \
         LEVEL 1",
    ),
    (
        "kappa",
        "a LEVEL 3 parameter: the saturation-field factor is not part of LEVEL 1",
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
