//! λ appears in the signal listings, under `.op` and `.tran` alike (#71).
//!
//! #59 made λ a label instead of an unknown, and every enumeration surface
//! kept iterating `node_index` — so λ was probeable by name at `.op`, listed
//! nowhere, and unreachable entirely at `.tran`. The property this file pins
//! is the one that was false whichever way #71 had been decided: λ appears in
//! **both** analyses' listings **or in neither**, and by-name probing agrees
//! with the enumeration. Option A (λ in) is the decision, matching the docs'
//! "still declared, still counted, still probeable".
//!
//! One deliberate wrinkle: a λ net hand-pinned by a voltage source
//! (`Vwl b_wl_0 0 …`, the documented way to label light entering from outside
//! the deck) stays a real MNA row. It must be listed exactly once — as the
//! row it is, not again as a label.

use fairchild_core::{
    dc_op_nr_with_registry, tran_nr_with_registry_opts, DeviceRegistry, SimOptions,
};
use fairchild_parser::parse_spice;

const LAMBDA_NM: f64 = 1310.0;

fn deck() -> String {
    format!(
        "\
* one laser, one waveguide — two λ nets
.optical_port a
.optical_port b
Xl a fc_cw_laser power_mW=1.0 wavelength_nm={LAMBDA_NM}
Xw a b fc_waveguide L_um=0.1 alpha_dB_cm=0
.op
"
    )
}

const WL: f64 = LAMBDA_NM * 1e-9;

#[test]
fn lambda_is_listed_and_probeable_at_dc() {
    let net = parse_spice(&deck()).expect("deck parses");
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP converges");

    for wl_net in ["a_wl_0", "b_wl_0"] {
        // Enumeration carries it...
        let listed = r
            .all_voltages()
            .find(|(name, _)| *name == wl_net)
            .unwrap_or_else(|| panic!("{wl_net} missing from all_voltages()"));
        // ...and agrees with the by-name probe, which is the pair that
        // disagreed before this test existed.
        let probed = r.node_voltage(wl_net).expect("probeable by name");
        assert_eq!(listed.1, probed, "{wl_net}: listing and probe disagree");
        assert!((probed - WL).abs() < 1e-18, "V({wl_net}) = {probed:e}");
    }

    // CSV and Nutmeg are the surfaces with consumers outside this repo.
    let mut csv = Vec::new();
    r.write_csv(&mut csv).unwrap();
    let csv = String::from_utf8(csv).unwrap();
    let header = csv.lines().next().unwrap();
    assert!(
        header.contains("V(a_wl_0)") && header.contains("V(b_wl_0)"),
        "{header}"
    );

    let mut raw = Vec::new();
    r.write_nutmeg(&mut raw, "t").unwrap();
    let raw = String::from_utf8(raw).unwrap();
    assert!(raw.contains("v(a_wl_0)"), "λ missing from Nutmeg:\n{raw}");
    assert_nutmeg_count_matches_block(&raw);
}

#[test]
fn lambda_is_listed_and_probeable_in_transient() {
    let net = parse_spice(&deck()).expect("deck parses");
    let opts = SimOptions::from_netlist(&net);
    let r = tran_nr_with_registry_opts(&net, 1e-12, 5e-12, &DeviceRegistry::new(), &opts)
        .expect("transient runs");

    for wl_net in ["a_wl_0", "b_wl_0"] {
        assert!(
            r.lambda.contains_key(wl_net),
            "{wl_net} missing from TranResult.lambda"
        );
        // The half that was unreachable, not merely unlisted: `.tran` answers
        // the same by-name question `.op` does.
        let probed = r
            .voltage_at(wl_net, 2e-12)
            .unwrap_or_else(|| panic!("V({wl_net}) unanswerable at .tran"));
        assert!((probed - WL).abs() < 1e-18, "V({wl_net}) = {probed:e}");
    }

    let mut csv = Vec::new();
    r.write_csv(&mut csv).unwrap();
    let csv = String::from_utf8(csv).unwrap();
    let header = csv.lines().next().unwrap();
    assert!(
        header.contains("V(a_wl_0)") && header.contains("V(b_wl_0)"),
        "{header}"
    );
    // Every λ column carries the wavelength on every row, since a label does
    // not move: spot-check the last data row.
    let cols: Vec<&str> = header.split(',').collect();
    let a_col = cols.iter().position(|c| *c == "V(a_wl_0)").unwrap();
    let last = csv.lines().last().unwrap().split(',').nth(a_col).unwrap();
    let v: f64 = last.parse().unwrap();
    assert!(
        (v - WL).abs() < 1e-18,
        "λ column reads {v:e} at the last row"
    );

    let mut raw = Vec::new();
    r.write_nutmeg(&mut raw, "t").unwrap();
    let raw = String::from_utf8(raw).unwrap();
    assert!(raw.contains("v(a_wl_0)"), "λ missing from Nutmeg:\n{raw}");
    assert_nutmeg_count_matches_block(&raw);
}

/// A λ net pinned as a real row by a hand-wired source is listed once, as the
/// row — `lambda_signals` must not repeat it as a label.
#[test]
fn a_hand_pinned_lambda_net_is_not_listed_twice() {
    let src = format!(
        "\
* external light: λ pinned by a source
.optical_port a
.optical_port b
Vwl a_wl_0 0 DC {WL:e}
Xw a b fc_waveguide L_um=0.1 alpha_dB_cm=0
.op
"
    );
    let net = parse_spice(&src).expect("deck parses");
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP converges");
    let n = r
        .all_voltages()
        .filter(|(name, _)| *name == "a_wl_0")
        .count();
    assert_eq!(n, 1, "a_wl_0 listed {n} times");
}

/// `No. Variables:` and the `Variables:` block are two statements of one
/// number; #71's fix touched both, so hold them together.
fn assert_nutmeg_count_matches_block(raw: &str) {
    let declared: usize = raw
        .lines()
        .find_map(|l| l.strip_prefix("No. Variables: "))
        .expect("count line")
        .trim()
        .parse()
        .expect("count parses");
    let listed = raw
        .lines()
        .skip_while(|l| !l.starts_with("Variables:"))
        .skip(1)
        .take_while(|l| !l.starts_with("Values:"))
        .count();
    assert_eq!(
        declared, listed,
        "variable count disagrees with the block:\n{raw}"
    );
}
