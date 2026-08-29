//! A subcircuit whose defaults are expressions, which is what a PDK wrapper is.
//!
//! Every foundry cell is an `inline subckt` with a long list of extracted
//! defaults, and many of those are expressions over the other parameters. The
//! suite had no case for one, and the gap cost a real defect: `.subckt` split its
//! header on whitespace while `.param` two lines below used a brace-aware
//! splitter, so a default containing a space became several tokens, none holding
//! an `=`, and each was counted as another port.
//!
//! A three-port cell declared eleven, and every call to it was refused with an
//! error naming the *call site*. Found on a GF 45SPCLO ESD diode, where the header
//! is forty continuation lines and the phantom ports were fragments of one
//! extracted expression (#105).

use fairchild_core::{dc_op_nr_with_registry_opts, options::SimOptions, DeviceRegistry};
use fairchild_parser::parse_spice;

/// Solve a deck and return `I(v1)`.
fn current(deck: &str) -> f64 {
    let net = parse_spice(deck).unwrap_or_else(|e| panic!("parse failed:\n{deck}\n{e}"));
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&net.models);
    let opts = SimOptions::from_netlist(&net);
    dc_op_nr_with_registry_opts(&net, &registry, &opts)
        .unwrap_or_else(|e| panic!("solve failed:\n{deck}\n{e:?}"))
        .vsrc_current("v1")
        .expect("I(v1)")
}

/// A default whose value is an expression containing spaces is one parameter, not
/// several ports.
///
/// The spaces are the whole point: `{(a-b)*(n+1)}` always worked, and
/// `{(a - b)*(n + 1)}` did not. Both spellings are here so the test states which
/// difference matters.
#[test]
fn a_subckt_default_may_be_an_expression_containing_spaces() {
    const TIGHT: &str = "* tight spelling\n\
        .model dm D (is=1e-14 n=1)\n\
        .subckt cell p n w=2 nf=3 scale={(w*1e0-0.5)*(nf+1)}\n\
        d1 p n dm area={scale}\n\
        .ends cell\n\
        v1 a 0 DC 0.6\n\
        x1 a 0 cell\n\
        .op\n";
    const SPACED: &str = "* the same header, written with spaces\n\
        .model dm D (is=1e-14 n=1)\n\
        .subckt cell p n w=2 nf=3 scale={(w * 1e0 - 0.5) * (nf + 1)}\n\
        d1 p n dm area={scale}\n\
        .ends cell\n\
        v1 a 0 DC 0.6\n\
        x1 a 0 cell\n\
        .op\n";

    let tight = current(TIGHT);
    let spaced = current(SPACED);
    assert!(
        (tight - spaced).abs() / tight.abs() < 1e-12,
        "whitespace inside a braced default must not change the circuit: \
         {tight:.9e} against {spaced:.9e}"
    );

    // And the default has to have been *used*, or this passes on two circuits
    // that both ignored it. `scale = (2 - 0.5)*4 = 6`, so `AREA=6`.
    let plain = current(
        "* no scaling\n\
         .model dm D (is=1e-14 n=1)\n\
         v1 a 0 DC 0.6\n\
         d1 a 0 dm\n\
         .op\n",
    );
    assert!(
        (tight / plain - 6.0).abs() < 1e-6,
        "the default must reach AREA: {tight:.6e} is {:.4}x the unscaled \
         {plain:.6e}, expected 6",
        tight / plain
    );
}

/// `name = value` written with spaces around the `=` is one assignment.
///
/// Foundry headers align their defaults in columns — `pwbp    =1e35` — so the `=`
/// is routinely separated from its name, its value, or both. Splitting on
/// whitespace turns that into two or three tokens, and the ones without an `=`
/// become ports.
#[test]
fn a_spaced_equals_on_a_header_is_still_one_assignment() {
    let tight = current(
        "* no spaces around =\n\
         .model dm D (is=1e-14 n=1)\n\
         .subckt cell p n area1=4\n\
         d1 p n dm area={area1}\n\
         .ends cell\n\
         v1 a 0 DC 0.6\n\
         x1 a 0 cell\n\
         .op\n",
    );
    for spelling in ["area1 = 4", "area1= 4", "area1 =4", "area1    =4"] {
        let deck = format!(
            "* spaced equals\n\
             .model dm D (is=1e-14 n=1)\n\
             .subckt cell p n {spelling}\n\
             d1 p n dm area={{area1}}\n\
             .ends cell\n\
             v1 a 0 DC 0.6\n\
             x1 a 0 cell\n\
             .op\n"
        );
        let got = current(&deck);
        assert!(
            (got - tight).abs() / tight.abs() < 1e-12,
            "'{spelling}' must read as one assignment, not as ports: \
             {got:.9e} against {tight:.9e}"
        );
    }
}

/// An expression default may name a port, which is where a whitespace split does
/// the most damage: the fragments look like plausible node names.
#[test]
fn an_expression_default_over_several_parameters_survives() {
    let deck = "* several parameters, several spaces\n\
        .model dm D (is=1e-14 n=1)\n\
        .subckt cell p n a=1 b=2 c=3 k={a + b + c} m={k * 2 - 4}\n\
        d1 p n dm area={m}\n\
        .ends cell\n\
        v1 a 0 DC 0.6\n\
        x1 a 0 cell\n\
        .op\n";
    let plain = current("* unscaled\n.model dm D (is=1e-14 n=1)\nv1 a 0 DC 0.6\nd1 a 0 dm\n.op\n");
    // k = 6, m = 8.
    let got = current(deck);
    assert!(
        (got / plain - 8.0).abs() < 1e-6,
        "chained expression defaults must resolve: {:.4}x, expected 8",
        got / plain
    );
}

/// The header and a `.param` line agree about what a token is.
///
/// They did not, and that was the defect: two splitters on one file. This states
/// the invariant directly, so a future change that reintroduces a second
/// tokeniser fails here rather than on a foundry deck.
#[test]
fn a_header_default_and_a_param_line_read_the_same_value() {
    let via_param = current(
        "* the value on a .param line\n\
         .model dm D (is=1e-14 n=1)\n\
         .param scale={(2 * 1e0 - 0.5) * (3 + 1)}\n\
         v1 a 0 DC 0.6\n\
         d1 a 0 dm area={scale}\n\
         .op\n",
    );
    let via_header = current(
        "* the same value as a subckt default\n\
         .model dm D (is=1e-14 n=1)\n\
         .subckt cell p n scale={(2 * 1e0 - 0.5) * (3 + 1)}\n\
         d1 p n dm area={scale}\n\
         .ends cell\n\
         v1 a 0 DC 0.6\n\
         x1 a 0 cell\n\
         .op\n",
    );
    assert!(
        (via_param - via_header).abs() / via_param.abs() < 1e-12,
        "a `.subckt` header and a `.param` line must tokenise identically: \
         {via_param:.9e} against {via_header:.9e}"
    );
}

/// When the arity really is wrong, the error names the ports it read.
///
/// The count alone points at the call site, and the fault is usually in the
/// header. Seeing the port list is what turns this from a puzzle into a glance.
#[test]
fn a_port_count_error_lists_the_ports_it_read() {
    let err = parse_spice(
        "* genuinely mismatched\n\
         .subckt cell p n g\n\
         r1 p n 1k\n\
         .ends cell\n\
         v1 a 0 DC 1\n\
         x1 a 0 cell\n\
         .op\n",
    )
    .expect_err("two nets against three ports must be refused");
    let msg = format!("{err}");
    for needle in ["cell", "3", "p n g"] {
        assert!(
            msg.contains(needle),
            "the refusal must name {needle}, or the reader cannot tell whether the \
             call or the header is wrong: {msg}"
        );
    }
}
