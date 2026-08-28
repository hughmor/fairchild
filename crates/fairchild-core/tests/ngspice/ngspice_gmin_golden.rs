//! Where `gmin` goes, checked against ngspice on the same deck.
//!
//! `gmin` exists because a pn junction's conductance `dI/dV` collapses to
//! nothing in reverse bias. A node reaching the rest of the circuit only through
//! reverse-biased junctions then has a row of almost zeros, the Jacobian is
//! near-singular, and Newton either cannot pivot or takes an enormous step
//! because it is dividing by ~0. A tiny fixed conductance in parallel keeps every
//! junction conducting *something*.
//!
//! **Where** it is placed is the part that is easy to get wrong, and it is not a
//! free choice — it changes answers:
//!
//! * ngspice puts it **across each junction**, between that junction's own two
//!   terminals. Nodes get none directly; a node only sees `gmin` by having a
//!   junction attached to it.
//! * fairchild used to add it to **every node's diagonal, to ground,
//!   unconditionally**, and no device model mentioned `gmin` at all.
//!
//! The two coincide when a junction has one end grounded, which is why the
//! obvious diode probe agreed and the divergence went unnoticed for so long.
//! They differed in three ways, all of which this file pins:
//!
//! | | ngspice | fairchild, before |
//! |---|---|---|
//! | node with no junction (1 mA into 1 GΩ) | 1e6 | 999001 |
//! | …with `.options gmin=1e-6` | 1e6 | 999 |
//! | reverse junction between two non-ground nodes | 1.01e-12 | 1e-14 |
//! | reverse-biased BJT junction | 2.0002e-12 | 2e-16 |
//!
//! The models' `gmin` was also *Jacobian-only*: `diode.rs` added it to
//! `gd_junction` and the Norton form `jeq = Id − gd·Vd` then cancelled it out of
//! the current exactly, so it conditioned the matrix and carried nothing. That
//! is a defensible technique and it is not what SPICE's `GMIN` is.
//!
//! One divergence used to remain: ngspice put a second `gmin` across the
//! collector-substrate junction and this simulator had no such junction, so a
//! reverse-biased BJT read half of ngspice's leakage. That is closed (#97 §3), and
//! `a_reverse_biased_bjt_junction_leaks_one_gmin_per_junction` now compares
//! against ngspice instead of documenting why it could not.
//!
//! Every case here is a *runtime* comparison rather than a stored number, so it
//! cannot go stale against a new ngspice, and so the reference is visible in the
//! deck instead of being folded into a constant.

use std::collections::HashMap;
use std::io::Write;
use std::process::Command;

use fairchild_core::options::SimOptions;
use fairchild_core::{dc_op_nr_with_registry_opts, DeviceRegistry};
use fairchild_parser::parse_spice;

use super::ngspice_golden::find_ngspice;

/// What to read out of both simulators: a node voltage or a source current.
enum Probe {
    V(&'static str),
    I(&'static str),
}

impl Probe {
    fn ngspice(&self) -> String {
        match self {
            Probe::V(n) => format!("v({n})"),
            Probe::I(n) => format!("i({n})"),
        }
    }
}

/// fairchild's answer for `probe` on `deck`.
fn fairchild(deck: &str, probe: &Probe) -> f64 {
    let net = parse_spice(deck).unwrap_or_else(|e| panic!("parse: {e:?}"));
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&net.models);
    // `from_netlist`, not `SimOptions::default()`: `dc_op_nr` uses the latter and
    // therefore ignores a deck's `.options`, which quietly compared fairchild at
    // gmin=1e-12 against ngspice at gmin=1e-6 the first time this ran.
    let opts = SimOptions::from_netlist(&net);
    let r = dc_op_nr_with_registry_opts(&net, &registry, &opts)
        .unwrap_or_else(|e| panic!("fairchild failed: {e:?}"));
    match probe {
        Probe::V(n) => r.node_voltage(n).expect("node"),
        Probe::I(n) => r.vsrc_current(n).expect("source"),
    }
}

/// ngspice's answer for `probe` on the same deck, or `None` if ngspice is absent.
///
/// The deck is shared verbatim; only the `.control` block that asks for the
/// number is appended, so there is no chance of the two simulators being given
/// different circuits.
fn ngspice(deck: &str, probe: &Probe) -> Option<f64> {
    let ng = find_ngspice()?;
    let name = probe.ngspice();
    let full = format!(
        "{}\n.control\nset numdgt=12\nop\nprint {name}\n.endc\n.end\n",
        deck.trim_end().trim_end_matches(".op").trim_end()
    );
    let mut tmp = tempfile::NamedTempFile::new().ok()?;
    tmp.write_all(full.as_bytes()).ok()?;
    let out = Command::new(ng).arg("-b").arg(tmp.path()).output().ok()?;
    let text = String::from_utf8_lossy(&out.stdout);
    let mut map = HashMap::new();
    for line in text.lines() {
        if let Some((lhs, rhs)) = line.trim().split_once('=') {
            if let Ok(v) = rhs.split_whitespace().next()?.parse::<f64>() {
                map.insert(lhs.trim().to_lowercase(), v);
            }
        }
    }
    // A missing value means ngspice ran and did not answer, which is a broken
    // test rather than an absent dependency — `find_ngspice` already settled that.
    Some(*map.get(&name).unwrap_or_else(|| {
        panic!("ngspice did not report {name}.\nDeck:\n{full}\nOutput:\n{text}")
    }))
}

/// Assert the two agree on `probe` to `rel`.
fn agree(what: &str, deck: &str, probe: Probe, rel: f64) {
    let Some(want) = ngspice(deck, &probe) else {
        eprintln!("skipping {what}: ngspice not found");
        return;
    };
    let got = fairchild(deck, &probe);
    let err = (got - want).abs() / want.abs().max(f64::MIN_POSITIVE);
    assert!(
        err < rel,
        "{what}: fairchild {got:.6e}, ngspice {want:.6e} (rel {err:.2e} > {rel:.0e})\n\
         Deck:\n{deck}"
    );
}

const DM: &str = ".model dm D (IS=1e-14 N=1)";

/// A node with no junction on it gets no `gmin`, so it reads Ohm's law.
///
/// 1 mA into 1 TΩ is 1e9 V. A nodal `gmin` of 1e-12 S is a second 1 TΩ in
/// parallel and halves it — which is what fairchild did.
#[test]
fn a_node_with_no_junction_gets_no_gmin() {
    // 1 GΩ rather than 1 TΩ: the latter's operating point is 1e9 V, which needs
    // more than the default iteration limit to walk to from a zero seed. 1 GΩ
    // still puts a nodal gmin 0.1% into the answer, which is six orders above the
    // tolerance here.
    agree(
        "1 mA into 1 GΩ",
        "* no junction\nI1 0 p DC 1m\nRs p 0 1G\n.op\n",
        Probe::V("p"),
        1e-9,
    );
}

/// …and raising `gmin` six decades must not move it either.
///
/// The sharper form of the same statement: if `gmin` reached resistor nodes, this
/// would collapse from 1e6 V to about 1e3.
#[test]
fn raising_gmin_does_not_reach_a_resistor_node() {
    agree(
        "1 mA into 1 GΩ at gmin=1e-6",
        "* no junction, big gmin\n.options gmin=1e-6\nI1 0 p DC 1m\nRs p 0 1G\n.op\n",
        Probe::V("p"),
        1e-9,
    );
}

/// A reverse-biased junction *does* get `gmin`, and it dominates `IS`.
///
/// `IS = 1e-14` against `gmin·1 V = 1e-12`: two orders apart, so this reads
/// `gmin` almost directly.
#[test]
fn a_reverse_biased_junction_leaks_gmin() {
    agree(
        "reverse diode to ground",
        &format!("* reverse diode\n{DM}\nV1 a 0 DC -1\nD1 a 0 dm\n.op\n"),
        Probe::I("v1"),
        1e-6,
    );
}

/// And the leakage tracks `gmin` when it moves.
#[test]
fn junction_leakage_tracks_gmin() {
    for g in ["1e-9", "1e-6"] {
        agree(
            &format!("reverse diode at gmin={g}"),
            &format!("* reverse diode\n.options gmin={g}\n{DM}\nV1 a 0 DC -1\nD1 a 0 dm\n.op\n"),
            Probe::I("v1"),
            1e-6,
        );
    }
}

/// The case the two policies disagree on: a junction whose terminals are both
/// off ground.
///
/// `gmin` belongs *across the diode*, so it shows up in `V1`'s current. Put it
/// from each node to ground instead and the leakage at `b` returns through `V2`,
/// leaving `V1` carrying only `IS` — a hundredfold difference.
#[test]
fn a_reverse_junction_between_two_non_ground_nodes_leaks_gmin() {
    agree(
        "reverse diode between two pinned nodes",
        &format!("* floating junction\n{DM}\nV1 a 0 DC 0\nV2 b 0 DC 1\nD1 a b dm\n.op\n"),
        Probe::I("v1"),
        1e-6,
    );
}

/// A BJT's junctions get `gmin` too — one per junction, and it agrees with
/// ngspice.
///
/// Base and emitter grounded, collector at 1 V. The collector current is the two
/// reverse-biased junctions on that node, base-collector and collector-substrate,
/// so it is `gmin`-dominated (`gmin·1 V = 1e-12` against `IS = 1e-16`).
///
/// # This test used to assert half of ngspice's answer
///
/// It read `1·gmin·V` and said so, with a comment explaining why agreement was not
/// available: the collector-substrate junction was not modelled, so there was
/// nothing for the second `gmin` to cross. Asserting ngspice's total would have
/// meant asserting a junction this simulator did not have.
///
/// The junction exists now (#97 §3), so the comparison is available and is made.
/// Which junction the second `gmin` belongs to is identified in
/// `ngspice_bjt_golden::a_reverse_biased_bjt_leaks_two_gmin_one_per_junction`, by
/// pinning the substrate at the collector potential and watching exactly one
/// `gmin·V` disappear.
#[test]
fn a_reverse_biased_bjt_junction_leaks_one_gmin_per_junction() {
    for g in [1e-12, 1e-9, 1e-6] {
        let deck = format!(
            "* bjt leakage
.options gmin={g:e}
.model qm NPN (IS=1e-16 BF=100)
V1 c 0 DC 1
Q1 c 0 0 qm
.op
"
        );
        // Two junctions on the collector node, each `gmin·1V`.
        let want = 2.0 * g;
        let got = fairchild(&deck, &Probe::I("v1")).abs();
        let rel = (got - want).abs() / want;
        assert!(
            rel < 1e-3,
            "gmin={g:e}: collector leakage {got:.6e}, expected two gmin·1V = \
             {want:.6e} (rel {rel:.2e}). Half this means the substrate junction \
             went away again; `2·IS` would mean the junction gmin is back to \
             being Jacobian-only."
        );
        let Some(ng) = ngspice(&deck, &Probe::I("v1")) else {
            continue;
        };
        let rel = (got - ng.abs()).abs() / ng.abs();
        assert!(
            rel < 1e-3,
            "gmin={g:e}: fairchild {got:.6e}, ngspice {:.6e} (rel {rel:.2e}). \
             This comparison was unavailable before the substrate junction \
             existed, and closing that gap is what made it assertable.",
            ng.abs()
        );
    }
}

/// The common case must not move: a forward-biased junction is many orders above
/// `gmin`, so where `gmin` sits cannot matter there.
///
/// Included because the fix moves `gmin` around, and a change that only ever
/// makes extreme circuits right while quietly moving ordinary ones would be a bad
/// trade. 1 mA through a diode is 1e11 times `gmin`'s contribution.
#[test]
fn a_forward_biased_junction_is_unaffected() {
    // 1e-5, not tighter: the two simulators disagree by 1.7e-6 on a forward
    // diode for reasons unrelated to gmin — a different thermal-voltage constant
    // — and gmin's own contribution here is 6.6e-13 A against 1e-3 A, eleven
    // orders down. Tightening this would be pinning the Vt constant, which is
    // `ngspice_diode_golden`'s subject, not this file's.
    agree(
        "forward diode at 1 mA",
        &format!("* forward diode\n{DM}\nIb 0 b DC 1m\nD1 b 0 dm\n.op\n"),
        Probe::V("b"),
        1e-5,
    );
    agree(
        "divider, no junctions at all",
        "* divider\nV1 a 0 DC 1\nR1 a m 1k\nR2 m 0 3k\n.op\n",
        Probe::V("m"),
        1e-9,
    );
}
