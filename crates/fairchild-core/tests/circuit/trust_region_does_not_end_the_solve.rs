//! A clamped Newton step is never the last one.
//!
//! `vmax` bounds how far an iterate may move in one step. The convergence test
//! is `|Δx| < abstol + reltol·|x|`, and that threshold *grows with `|x|`* — so on
//! a node heading for hundreds of volts the bound catches up with the clamp, and
//! a walk that is still being cut short every iteration looks finished. The
//! solver then stops wherever it happens to be and reports success (#90).
//!
//! Three things had to change together, and each is asserted below:
//!
//! 1. Convergence is judged on the **Newton step** — what the linear solve asked
//!    for — not on the clamped step actually taken. Once clamped, the step taken
//!    is `vmax` on the binding row by construction and says nothing about how far
//!    the iterate still has to go.
//! 2. A clamped iteration cannot be the last one, because a clamped step is
//!    shorter than the one that would reach the solution.
//! 3. The trust region is `vmax + 1e-3·|v|` per row rather than a flat `vmax`,
//!    so (2) is affordable: a node whose answer is far from the seed climbs
//!    geometrically instead of walking `vmax` at a time. Without this, refusing
//!    to stop early turns a wrong answer into a four-minute one.
//!
//! Every expectation here is a closed form. That matters more than usual: the
//! defect produced a *plausible* number every time, and a test that compared
//! against another run of the same solver would have agreed with it.

use fairchild_core::{dc_op_nr, dc_op_nr_opts, options::SimOptions};
use fairchild_parser::parse_spice;

/// An ideal VCCS into a load: linear, and one solve from its answer.
///
/// `G1 o 0 …` pulls `gm·V(a)` out of node `o` — SPICE's convention is that the
/// current flows from `n+` to `n-` *through* the source — so `RL` must supply it
/// and `V(o) = −gm·V(a)·RL` exactly. With
/// `gm = 10`, `V(a) = 0.1` and `RL = 1 kΩ` that is **−1000 V**, and `V(a)` is
/// pinned by a voltage source so it must read exactly 0.1.
///
/// This is #90's reproduction with the Verilog-A taken out — the original used
/// OpenVAF's `vccs.va`, which cannot live in this repository. A `G` element is
/// the same circuit.
fn vccs_deck(gm: f64, rl: f64, vin: f64) -> String {
    format!(
        "* ideal VCCS into a load\n\
         Va a 0 DC {vin}\n\
         G1 o 0 a 0 {gm}\n\
         RL o 0 {rl}\n\
         .op\n"
    )
}

fn solve(deck: &str) -> (f64, f64) {
    let net = parse_spice(deck).expect("parse");
    let r = dc_op_nr(&net).expect("this circuit is linear; it must solve");
    (
        r.node_voltage("a").expect("V(a)"),
        r.node_voltage("o").expect("V(o)"),
    )
}

/// The reported case: a node the deck pins at 0.1 V read 0.0502 V, as a
/// converged operating point, in two iterations.
#[test]
fn a_pinned_node_reads_exactly_what_the_source_says() {
    let (va, vo) = solve(&vccs_deck(10.0, 1e3, 0.1));
    // `V(a)` is a voltage source's own node. There is no tolerance argument for
    // it being anything but 0.1: the equation is `V(a) = 0.1`.
    assert!(
        (va - 0.1).abs() < 1e-9,
        "V(a) = {va}, and the deck pins it at 0.1 — the trust region ended the \
         solve early (#90)"
    );
    // …and the output followed it. −gm·V(a)·RL = −1000, less the 1 mV that gmin
    // on `o` steals at this impedance.
    assert!(
        (vo + 1000.0).abs() < 1e-2,
        "V(o) = {vo}, expected ≈ −1000 (−gm·V(a)·RL)"
    );
}

/// The pathology scaled: the further the answer is from the seed, the more
/// iterations a flat `vmax` walk needs, and the more room `reltol·|x|` has to
/// overtake it. A fixed trust region fails progressively; a relative one does
/// not.
///
/// Three decades of gain, all with the same closed form.
#[test]
fn the_answer_does_not_degrade_as_it_moves_away_from_the_seed() {
    for gm in [1.0, 10.0, 100.0, 1000.0] {
        let want = -gm * 0.1 * 1e3;
        let (va, vo) = solve(&vccs_deck(gm, 1e3, 0.1));
        assert!(
            (va - 0.1).abs() < 1e-9,
            "gm={gm}: V(a) = {va}, pinned at 0.1"
        );
        // A relative bound, because gmin's share of the current grows with the
        // node voltage. At |V(o)| = 1e5 it is worth 0.1 V, which is 1e-6
        // relative — still four orders inside this.
        let rel = (vo - want).abs() / want.abs();
        assert!(
            rel < 1e-4,
            "gm={gm}: V(o) = {vo}, expected {want} (rel {rel:.2e}). The error \
             growing with gm is the signature: a flat vmax walk takes |V(o)|/vmax \
             iterations, and reltol·|x| overtakes it sooner the larger |x| gets."
        );
    }
}

/// Tightening `reltol` must not be what fixes it.
///
/// This is the sharpest statement of the defect. The old behaviour depended on
/// `reltol` — the walk stopped when `reltol·|x|` overtook the step — so the
/// answer moved when `reltol` moved, which is not how a converged answer
/// behaves. The same operating point at three tolerances, three decades apart.
#[test]
fn the_answer_does_not_depend_on_reltol() {
    let deck = vccs_deck(10.0, 1e3, 0.1);
    let net = parse_spice(&deck).expect("parse");
    let mut answers = Vec::new();
    for reltol in [1e-3, 1e-6, 1e-9] {
        let mut opts = SimOptions::from_netlist(&net);
        opts.reltol = reltol;
        let r = dc_op_nr_opts(&net, &opts).expect("must solve");
        answers.push((reltol, r.node_voltage("o").expect("V(o)")));
    }
    let (_, first) = answers[0];
    for (reltol, v) in &answers {
        assert!(
            (v - first).abs() < 1e-6 * first.abs(),
            "V(o) at reltol={reltol} is {v}, but {first} at reltol=1e-3 — a \
             converged operating point does not move with the tolerance that \
             stopped the iteration"
        );
    }
}

/// The guard must not have been bought by refusing circuits that were fine.
///
/// A plain resistive divider, a diode, and a circuit whose node genuinely sits
/// at a large voltage: all must still solve, and to their closed forms.
#[test]
fn ordinary_circuits_still_solve() {
    // Divider: 1 V across 1k + 3k.
    let (_, _) = {
        let net = parse_spice("* divider\nVa a 0 DC 1\nR1 a m 1k\nR2 m 0 3k\n.op\n").unwrap();
        let r = dc_op_nr(&net).unwrap();
        let vm = r.node_voltage("m").unwrap();
        assert!((vm - 0.75).abs() < 1e-9, "divider gave {vm}, expected 0.75");
        (vm, 0.0)
    };
    // A diode at a known bias: V = N·Vt·ln(I/Is + 1), the same anchor
    // `model_parameter_diagnostics` uses.
    let net =
        parse_spice("* diode\n.model dm D (IS=1e-14 N=1)\nIb 0 b DC 1m\nD1 b 0 dm\n.op\n").unwrap();
    let vb = dc_op_nr(&net).unwrap().node_voltage("b").unwrap();
    let vt = 1.380649e-23 * 300.15 / 1.602176634e-19;
    let want = vt * (1e-3_f64 / 1e-14 + 1.0).ln();
    assert!(
        (vb - want).abs() < 1e-4 * want,
        "diode gave {vb}, expected {want}"
    );
}

/// A high-impedance node reaches Ohm's law, not a fraction of it.
///
/// 1 mA into a shunt: `V = I·R`. The interesting part is what `R` does to the
/// walk — at `R = 1 GΩ` the operating point is 10^6 V, which a flat `vmax` of
/// 0.5 V would need two million iterations to reach, and the convergence bound
/// `vntol + reltol·|V|` overtakes it long before then. Before this fix the
/// solver reported **501988 V** at 1 GΩ and **900538 V** at 10 GΩ — half and a
/// tenth of the answer — as converged operating points.
///
/// # About the tolerance
///
/// `I·R` to within the solver's own convergence bound (`reltol`), because that
/// is the honest floor: below it, the number is where Newton stopped rather than
/// what the circuit is. Nothing else in this circuit can move the node —
/// `gmin` lives across pn junctions now, and a shunt resistor is not one
/// (`ngspice_gmin_golden.rs`), so no term stands between the answer and Ohm's
/// law in either direction.
///
/// The failure this exists to catch — a walk cut short by the flat `vmax` trust
/// region — is a factor of two or ten, five hundred times larger.
///
/// A 1 GΩ shunt is not contrived — it is a photodiode's dark resistance, a gate
/// leakage path, an ESD structure's off state.
#[test]
fn a_high_impedance_node_reaches_ohms_law() {
    // The bound the solver itself declared convergence on.
    let tol = SimOptions::default().reltol;
    for (r_str, r) in [("1G", 1e9), ("10G", 1e10)] {
        let deck = format!(
            "* high-impedance shunt
             .optical_port a
             XPD a p n fc_photodetector r_shunt={r_str}
             I1 0 p DC 1m
             Rn n 0 1
.op
"
        );
        let net = parse_spice(&deck).expect("parse");
        let v = dc_op_nr(&net)
            .unwrap_or_else(|e| panic!("r_shunt={r_str} must solve: {e:?}"))
            .node_voltage("p")
            .expect("V(p)");
        let want = 1e-3 * r;
        let rel = (v - want).abs() / want;
        assert!(
            rel < tol,
            "r_shunt={r_str}: V(p) = {v:.6e}, expected I·R = {want:.6e}              (rel {rel:.2e}). Half or a tenth of this means the walk was cut              short and the tolerance caught up with it — #90."
        );
    }
}
