//! Regression tests for fc_circulator: the three routes, one at a time.
//!
//! Bidir-only — instantiating without `.options enable_bidirectional=1` is an
//! error naming the element.
//!
//! Wire convention is the along-chain one every other device uses: `fw` runs
//! from port 1 toward port 3, the same direction at all three ports. Port 1
//! therefore behaves like a waveguide's `in` port (it reads `fw`, drives `bw`)
//! and ports 2 and 3 like `out` ports (they drive `fw`, read `bw`). Each case
//! below drives only wires the circulator does not, and reads only wires it
//! does.
//!
//! Whether the routes compose into a circuit — and conserve power once they do
//! — is `tests/bidirectional_composition.rs`; these three pin the routing table
//! itself, one entry per test.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

fn run(drives: &str) -> fairchild_core::newton::NrResult {
    let src = format!(
        ".options enable_bidirectional=1\n\
         .optical_port p1\n.optical_port p2\n.optical_port p3\n\
         Xcirc p1 p2 p3 fc_circulator\n\
         Vwl p1_wl_0 0 DC 1.55e-6\n\
         {drives}.op\n.end\n"
    );
    let net = parse_spice(&src).expect("netlist should parse");
    dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP should converge")
}

fn v(r: &fairchild_core::newton::NrResult, wire: &str) -> f64 {
    r.node_voltage(wire)
        .unwrap_or_else(|e| panic!("{wire}: {e}"))
}

/// Light entering port 1 leaves at port 2, and nowhere else.
#[test]
fn light_in_at_port_one_leaves_at_port_two() {
    let r = run("Vre p1_re_fw_0 0 DC 1.0\nVim p1_im_fw_0 0 DC 0.0\n");
    assert!(
        (v(&r, "p2_re_fw_0") - 1.0).abs() < 1e-9,
        "port 2 should carry it forward; got {}",
        v(&r, "p2_re_fw_0")
    );
    assert!(v(&r, "p2_im_fw_0").abs() < 1e-9);
    assert!(
        v(&r, "p3_re_fw_0").abs() < 1e-9,
        "port 3 gets nothing until something comes back into port 2; got {}",
        v(&r, "p3_re_fw_0")
    );
    // λ is a label on the channel: both `out`-role ports take port 1's.
    assert!((v(&r, "p2_wl_0") - 1.55e-6).abs() < 1e-18);
    assert!((v(&r, "p3_wl_0") - 1.55e-6).abs() < 1e-18);
}

/// Light coming back into port 2 leaves at port 3 — the measurement path.
/// Port 2 is an `out`-role port, so a returning wave arrives on its `bw` wires.
#[test]
fn light_back_in_at_port_two_leaves_at_port_three() {
    let r = run("Vre p2_re_bw_0 0 DC 1.0\nVim p2_im_bw_0 0 DC 0.0\n");
    assert!(
        (v(&r, "p3_re_fw_0") - 1.0).abs() < 1e-9,
        "port 3 should carry the return; got {}",
        v(&r, "p3_re_fw_0")
    );
    assert!(v(&r, "p3_im_fw_0").abs() < 1e-9);
    assert!(
        v(&r, "p2_re_fw_0").abs() < 1e-9,
        "nothing entered port 1, so port 2's forward wire stays dark; got {}",
        v(&r, "p2_re_fw_0")
    );
}

/// Light coming back into port 3 leaves at port 1 — the third leg, and the one
/// that decides whether the cycle closes or stops at port 3.
#[test]
fn light_back_in_at_port_three_leaves_at_port_one() {
    let r = run("Vre p3_re_bw_0 0 DC 1.0\nVim p3_im_bw_0 0 DC 0.0\n");
    assert!(
        (v(&r, "p1_re_bw_0") - 1.0).abs() < 1e-9,
        "port 1 should carry it back out; got {}",
        v(&r, "p1_re_bw_0")
    );
    assert!(v(&r, "p1_im_bw_0").abs() < 1e-9);
}

/// Circulator must refuse to instantiate without bidir mode — as an error the
/// caller can handle, not a panic.
///
/// It used to `panic!` inside `setup_instance`, which runs inside the registry
/// factory: that aborts the CLI with a backtrace and crosses pyo3 as a
/// `PanicException`. Now it declines, and `build_devices_with_footprints`
/// raises a `ParameterError` naming the element. The count in that error is
/// real — 9 wires supplied against the 15 a bidirectional circulator needs —
/// but it is a symptom, so the device also warns with the actual fix.
#[test]
fn circulator_requires_bidir() {
    let netlist = "\
.optical_port p0
.optical_port p1
.optical_port p2
Xcirc p0 p1 p2 fc_circulator
.op
.end
";
    let net = parse_spice(netlist).unwrap();
    let Err(err) = dc_op_nr_with_registry(&net, &DeviceRegistry::new()) else {
        panic!("a circulator without enable_bidirectional=1 must not build");
    };
    let msg = format!("{err:?}");
    assert!(
        msg.contains("xcirc") && msg.contains("fc_circulator"),
        "the error must name the element and the model: {msg}"
    );
}
