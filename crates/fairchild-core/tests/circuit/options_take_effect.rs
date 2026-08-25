//! Every solver option arrives, and every solver option does something.
//!
//! `.options method=gear` parsed correctly, landed in `SimOptions.method`
//! correctly, and then ran Backward Euler, because the fixed-step integrator
//! never asked for the extra history BDF-2 needs (#93). The output was
//! byte-identical to `method=be`. Nothing in a large test suite noticed, because
//! every test that exercised an option asserted something *about that option's
//! own feature* — and there was no test whose subject was "the option took
//! effect at all".
//!
//! That is the gap this file is for, and it is deliberately two separate
//! questions:
//!
//! * **Arrival.** `set(key, value)` changes the field it names and no other. One
//!   case per key *and per alias*, generated from a table.
//! * **Effect.** Two values of the option produce different output. This is the
//!   half that catches #93, and the half a per-feature test never covers,
//!   because a feature test that only ever runs with the feature on cannot tell
//!   you the switch is wired up.
//!
//! # What keeps this from rotting
//!
//! [`every_option_is_accounted_for`] extracts the field names from
//! `SimOptions`'s own `Debug` output and fails if any is missing from the
//! tables. Rust has no reflection, so this is the closest thing to a compiler
//! check: add a field to `SimOptions` and this test tells you to say how it is
//! observable — or to write down, in `NOT_OBSERVABLE`, why it is not. The
//! exclusions are a reviewed list rather than an absence.

use std::collections::BTreeMap;

use fairchild_core::options::SimOptions;
use fairchild_core::{tran_nr_configured, DeviceRegistry};
use fairchild_parser::parse_spice;

// ─────────────────────────────────────────────────────────────────────────────
// Arrival: `set` writes the field it names, and nothing else
// ─────────────────────────────────────────────────────────────────────────────

/// `(key as a user writes it, value, the `SimOptions` field it must change)`.
///
/// Aliases get their own row: `maxstep` and `max_step` are two spellings a deck
/// can use and each is a separate chance for the match arm to point at the wrong
/// field.
const ARRIVES: &[(&str, &str, &str)] = &[
    ("reltol", "1e-5", "reltol"),
    ("abstol", "1e-15", "abstol"),
    ("vntol", "1e-9", "vntol"),
    ("temptol", "1e-4", "temptol"),
    ("vmax", "0.25", "vmax"),
    ("gmin", "1e-13", "gmin"),
    ("itl1", "77", "itl1"),
    ("itl4", "88", "itl4"),
    ("maxstep", "1n", "max_step"),
    ("max_step", "2n", "max_step"),
    ("tstart", "1n", "tstart"),
    ("gmin_max", "1e-3", "gmin_max"),
    ("gminmax", "1e-2", "gmin_max"),
    ("srcsteps", "42", "srcsteps"),
    ("srcmax", "43", "srcsteps"),
    ("method", "be", "method"),
    ("method", "gear", "method"),
    ("solver", "sparse", "solver"),
    ("uic", "1", "uic"),
    ("nopnjlim", "1", "pnjlim"),
    ("variable_step", "1", "variable_step"),
    ("variablestep", "1", "variable_step"),
    ("trannoise", "1", "trannoise"),
    ("noiseseed", "7", "noiseseed"),
    ("noisescale", "2", "noisescale"),
    ("waveguide_delay", "1", "waveguide_delay"),
    ("cond_estimate", "1", "cond_estimate"),
    ("equilibrate", "1", "equilibrate"),
    ("sanity_check", "0", "sanity_check"),
    ("verbose", "1", "verbose"),
    ("lambda_center_nm", "1310", "lambda_center_m"),
    ("enable_bidirectional", "1", "bidirectional_propagation"),
    ("temp", "85", "temp_k"),
    ("tnom", "40", "temp_k"),
    ("max_rejections", "9", "max_rejections"),
];

/// `SimOptions`'s fields, by name, read off its own `Debug` output.
///
/// The only way to enumerate a struct's fields without reflection or a derive
/// macro, and enough for a completeness gate: `{:#?}` prints one `name: value,`
/// per line at a known indent.
fn fields(opts: &SimOptions) -> BTreeMap<String, String> {
    let text = format!("{opts:#?}");
    text.lines()
        .filter_map(|l| {
            let t = l.trim();
            let (name, value) = t.split_once(": ")?;
            // Field names only — nested struct contents are indented further and
            // their lines do not survive this, which is what we want.
            if !name
                .chars()
                .all(|c| c.is_ascii_lowercase() || c == '_' || c.is_ascii_digit())
            {
                return None;
            }
            Some((name.to_string(), value.trim_end_matches(',').to_string()))
        })
        .collect()
}

#[test]
fn set_changes_the_field_it_names_and_no_other() {
    let base = SimOptions::default();
    let base_fields = fields(&base);

    for (key, value, field) in ARRIVES {
        let mut opts = SimOptions::default();
        assert!(
            opts.set(key, value),
            "'{key}' was not recognised by SimOptions::set — an unknown key is \
             already warned about, so a row here that returns false is a typo in \
             this table or a removed option"
        );
        let after = fields(&opts);

        let changed: Vec<&String> = base_fields
            .keys()
            .filter(|k| base_fields[*k] != after[*k])
            .collect();

        assert_eq!(
            changed.len(),
            1,
            "'{key}={value}' should change exactly one field; it changed {changed:?}"
        );
        assert_eq!(
            changed[0], field,
            "'{key}={value}' changed {} instead of {field}",
            changed[0]
        );
    }
}

/// A value the option cannot parse must be reported, not silently discarded.
///
/// Every numeric arm used to read `parse_num(value).unwrap_or(self.field)`, so
/// `--opt reltol=banana` kept the default and said nothing. An unknown *key* was
/// already reported, which is what made the hole in the value half easy to miss:
/// the half-working diagnostic reads as a working one.
#[test]
fn an_unparseable_value_does_not_silently_keep_the_default() {
    // `set` still returns true — the key is real — so the check is that the
    // field did not move, and (by inspection of the warning path) that the user
    // was told. The warning itself goes to stderr, which a unit test cannot see;
    // `ci.yml` asserts the text.
    for (key, field) in [("reltol", "reltol"), ("itl1", "itl1"), ("vmax", "vmax")] {
        let base = fields(&SimOptions::default());
        let mut opts = SimOptions::default();
        opts.set(key, "banana");
        let after = fields(&opts);
        assert_eq!(
            base[field], after[field],
            "'{key}=banana' must not move {field}"
        );
    }
}

/// The `.options` route and the `--opt` / kwarg route are the same code, and
/// this is the assertion that keeps them that way.
#[test]
fn the_deck_route_and_the_key_value_route_agree() {
    for (key, value, _) in ARRIVES {
        // `.options` in a deck.
        let src = format!("* opt\nV1 a 0 DC 1\nR1 a 0 1k\n.options {key}={value}\n.op\n");
        let net = parse_spice(&src).expect("parse");
        let from_deck = SimOptions::from_netlist(&net);

        // The same key through `set`, which is what the CLI's `--opt` and every
        // Python kwarg funnel into.
        let mut from_kv = SimOptions::default();
        from_kv.set(key, value);

        assert_eq!(
            fields(&from_deck),
            fields(&from_kv),
            "'{key}={value}' differs between the .options route and the key/value route"
        );
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Effect: the option changes the answer
// ─────────────────────────────────────────────────────────────────────────────

/// Decks the effect table draws on.
///
/// Two, because the options divide that way. A linear circuit converges in one
/// Newton step whatever the tolerance, so `reltol` genuinely cannot show itself
/// on `RC` — the first version of this table claimed it could and the test said
/// otherwise, which is the mechanism working.
mod decks {
    /// One time constant, ten points across it. Sensitive to *how* it is
    /// integrated and to the step schedule; insensitive to convergence
    /// tolerances, being linear.
    pub const RC: &str =
        "* rc\nV1 in 0 PULSE(0 1 0 1p 1p 1 2)\nR1 in out 1k\nC1 out 0 1n\n.tran 100n 3u\n";

    /// A diode charging a capacitor: nonlinear, so Newton takes several
    /// iterations per step and a convergence tolerance decides where it stops.
    pub const DIODE_RC: &str = "* diode rc\n.model dm D (IS=1e-14 N=1)\n\
         V1 in 0 PULSE(0 2 0 1p 1p 1 2)\nD1 in out dm\nR1 out 0 10k\n\
         C1 out 0 1n\n.tran 100n 3u\n";
}

/// Run one deck under one option setting and return time *and* value as text.
///
/// Time is included deliberately: an option that changes only the schedule
/// (`maxstep`, `variable_step`) is a real effect and must not read as no change.
fn waveform(deck: &str, key: &str, value: &str) -> String {
    let net = parse_spice(deck).expect("parse");
    let mut opts = SimOptions::from_netlist(&net);
    assert!(opts.set(key, value), "{key} not recognised");
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&net.models);
    let r = tran_nr_configured(&net, 100e-9, 3e-6, &registry, &opts)
        .unwrap_or_else(|e| panic!("{key}={value}: {e:?}"));
    r.time
        .iter()
        .zip(r.node_voltages["out"].iter())
        .map(|(t, v)| format!("{t:.17e},{v:.17e}"))
        .collect::<Vec<_>>()
        .join(";")
}

/// `(key, value A, value B, deck)` — two settings that must not produce the same
/// waveform.
///
/// This is the table #93 would have failed. `method` is the first row for that
/// reason.
const CHANGES_THE_ANSWER: &[(&str, &str, &str, &str)] = &[
    // The one that was broken: BDF-2 was Backward Euler on this path, so
    // `be` and `gear` agreed to the last bit.
    ("method", "be", "gear", decks::RC),
    ("method", "be", "tr", decks::RC),
    ("method", "tr", "gear", decks::RC),
    // Convergence tolerances need a circuit that iterates.
    ("reltol", "1e-1", "1e-10", decks::DIODE_RC),
    ("abstol", "1e-6", "1e-16", decks::DIODE_RC),
    // gmin shunts every node.
    ("gmin", "1e-12", "1e-5", decks::RC),
    // Clamping the step changes how many points the run takes.
    ("maxstep", "100n", "10n", decks::RC),
    // A different integration schedule entirely.
    ("variable_step", "0", "1", decks::RC),
    // Junction limiting only exists where there is a junction.
    ("nopnjlim", "0", "1", decks::DIODE_RC),
    // Temperature reaches the diode's thermal voltage.
    ("temp", "27", "125", decks::DIODE_RC),
];

#[test]
fn two_settings_of_an_option_do_not_give_the_same_answer() {
    for (key, a, b, deck) in CHANGES_THE_ANSWER {
        let wa = waveform(deck, key, a);
        let wb = waveform(deck, key, b);
        assert_ne!(
            wa, wb,
            "{key}={a} and {key}={b} produced identical waveforms — the option \
             parses and lands in SimOptions, and then nothing reads it. This is \
             exactly #93: `method=gear` was Backward Euler because the fixed-step \
             integrator never asked for the history BDF-2 needs."
        );
    }
}

/// Fields whose effect this deck cannot show, and why. Reviewed, not forgotten.
///
/// Every entry is a claim that the field *is* covered elsewhere or *cannot* be
/// covered by an output comparison. If you are adding to this list, the bar is:
/// could a reader be convinced the option is wired up, from something?
const NOT_OBSERVABLE: &[(&str, &str)] = &[
    (
        "vntol",
        "the absolute floor of a voltage row's Newton bound, \
               `vntol + |x|·reltol`. Newton converges quadratically, so a step \
               sequence usually jumps clean over the gap between two vntol \
               values and lands on the same iterate — finding a pair that \
               differs is tuning a test to today's numbers. `tolerance.rs` \
               asserts `bound(row, 0) == vntol` directly, which is sharper than \
               any waveform comparison",
    ),
    (
        "vmax",
        "the Newton trust region. On a circuit that converges it changes the \
              iteration count, not the answer; where it *does* change the answer \
              that is #90 — a clamped step being reported as converged — which is \
              a bug to fix rather than behaviour to pin",
    ),
    (
        "temptol",
        "bounds thermal rows, and there are none here — \
                 `thermal_discipline.rs` covers it",
    ),
    (
        "itl1",
        "an iteration ceiling: raising it cannot change a converged answer, \
              and lowering it past convergence is a failure, not a different \
              number. `singular_is_not_nonconvergence.rs` covers the failure.",
    ),
    ("itl4", "as itl1, for the transient inner loop"),
    (
        "max_rejections",
        "only reachable on the variable-step path when LTE keeps \
                        rejecting; a linear RC never rejects",
    ),
    (
        "gmin_max",
        "the gmin-stepping homotopy start, reached only when direct \
                  Newton fails",
    ),
    (
        "srcsteps",
        "source-stepping granularity, reached only when direct Newton \
                  and gmin stepping both fail",
    ),
    (
        "tstart",
        "trims the output window rather than the solve; asserted by the \
                `.tran` card tests",
    ),
    (
        "uic",
        "needs a `.ic` card to differ from the operating point",
    ),
    (
        "solver",
        "a backend choice that must give the *same* answer — \
                `solver_klu.rs` asserts exactly that agreement",
    ),
    (
        "equilibrate",
        "row/column scaling that must not change the answer beyond \
                     round-off, so an inequality assertion would be wrong",
    ),
    (
        "cond_estimate",
        "prints a condition-number estimate; stderr, not output",
    ),
    ("verbose", "prints diagnostics; stderr, not output"),
    (
        "sanity_check",
        "gates a pre-solve netlist check that either warns or does \
                      not — `spice_support` covers the checks themselves",
    ),
    (
        "trannoise",
        "needs `.options trannoise=1` plus a noise source; \
                   `transient_noise.rs` covers it",
    ),
    (
        "noiseseed",
        "only meaningful under trannoise — see transient_noise.rs",
    ),
    (
        "noisescale",
        "only meaningful under trannoise — see transient_noise.rs",
    ),
    (
        "lambda_center_m",
        "photonic; the band centre only matters to a bundle port",
    ),
    (
        "bidirectional_propagation",
        "photonic; `bidirectional_option.rs` covers it",
    ),
    (
        "waveguide_delay",
        "photonic; opt-in group delay, `native_*` tests cover it",
    ),
];

/// The gate that keeps the two tables above honest.
///
/// Adding a field to `SimOptions` fails this until the field appears in
/// `ARRIVES` (it can be set), and in either `CHANGES_THE_ANSWER` (it does
/// something) or `NOT_OBSERVABLE` (with a reason). That is the whole mechanism:
/// a new option cannot be added silently, which is how #93 got in.
#[test]
fn every_option_is_accounted_for() {
    let all = fields(&SimOptions::default());

    let arrives: Vec<&str> = ARRIVES.iter().map(|(_, _, f)| *f).collect();
    let effect: Vec<&str> = CHANGES_THE_ANSWER.iter().map(|(k, _, _, _)| *k).collect();
    let excused: Vec<&str> = NOT_OBSERVABLE.iter().map(|(f, _)| *f).collect();

    let mut missing_arrival = Vec::new();
    let mut missing_effect = Vec::new();
    for name in all.keys() {
        if !arrives.contains(&name.as_str()) {
            missing_arrival.push(name.clone());
        }
        // `CHANGES_THE_ANSWER` is keyed by the option's user-facing name, which
        // is usually but not always the field name; accept either.
        let covered = effect.contains(&name.as_str())
            || excused.contains(&name.as_str())
            || ARRIVES
                .iter()
                .any(|(k, _, f)| f == name && effect.contains(k));
        if !covered {
            missing_effect.push(name.clone());
        }
    }

    assert!(
        missing_arrival.is_empty(),
        "SimOptions fields with no row in ARRIVES (add one, so setting the option \
         is known to reach the field): {missing_arrival:?}"
    );
    assert!(
        missing_effect.is_empty(),
        "SimOptions fields that are never shown to do anything: {missing_effect:?}. \
         Add a row to CHANGES_THE_ANSWER with two values that differ, or to \
         NOT_OBSERVABLE with the reason it cannot be shown that way. An option \
         that parses and is then never read is #93."
    );

    // And the excuses have to be about real fields, or the list is a place for
    // dead entries to accumulate.
    for (field, _) in NOT_OBSERVABLE {
        assert!(
            all.contains_key(*field),
            "NOT_OBSERVABLE names '{field}', which is not a SimOptions field"
        );
    }
}
