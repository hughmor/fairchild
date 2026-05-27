//! Regression tests for fc_circulator.
//!
//! Three-port circulator: 1→2, 2→3, 3→1.  Bidir-only — instantiating
//! without `.options enable_bidirectional=1` panics at setup_instance.
//!
//! Wire convention: at every port, `re_fw`/`im_fw` flow INWARD (toward
//! the device); `re_bw`/`im_bw` flow OUTWARD.  Hence external drivers
//! place their incoming signal on a port's `_re_fw_` net and read the
//! circulator's outgoing signal from the next port's `_re_bw_` net.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

/// Inject light at port 0 (port_0.re_fw = 1.0).  Verify it appears at
/// port 1's outward wires (port_1.re_bw) and NOT at port 2's outward.
#[test]
fn circulator_routes_port0_to_port1() {
    let netlist = "\
.options enable_bidirectional=1
.optical_port p0
.optical_port p1
.optical_port p2
* Drive port 0's fw (light entering port 0)
Vp0_re p0_re_fw_0 0 DC 1.0
Vp0_im p0_im_fw_0 0 DC 0.0
Vp0_wl p0_wl_0    0 DC 1.55e-6
* The bw side of port 0 is the circulator's output back toward whatever's
* upstream — pin to 0 with a probe so the matrix has a driver on it.
* (Note: in a real schematic this is the wire that returning light from
*  port 2 → port 0 routes through.)
* Other ports: nothing driving fw on port 1 / port 2 — circulator routes
* port_0.fw into port_1.bw via internal stamps.  No external loads on
* port_1.bw / port_2.bw means those wires are pure outputs we can probe.
Xcirc p0 p1 p2 fc_circulator
.op
.end
";
    let net = parse_spice(netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    let p1_re_bw = r.node_voltage("p1_re_bw_0").unwrap();
    let p1_im_bw = r.node_voltage("p1_im_bw_0").unwrap();
    let p2_re_bw = r.node_voltage("p2_re_bw_0").unwrap();
    let p2_im_bw = r.node_voltage("p2_im_bw_0").unwrap();
    assert!(
        (p1_re_bw - 1.0).abs() < 1e-9,
        "port_1.re_bw should be 1.0 (light from port 0); got {p1_re_bw}"
    );
    assert!(
        p1_im_bw.abs() < 1e-9,
        "port_1.im_bw should be 0; got {p1_im_bw}"
    );
    assert!(
        p2_re_bw.abs() < 1e-9,
        "port_2.re_bw should be 0 (no light entering port 1); got {p2_re_bw}"
    );
    assert!(p2_im_bw.abs() < 1e-9);
}

/// Reflection round-trip: inject at port 0, "reflect" at port 1 via an
/// external wire that ties port_1.re_fw back to port_1.re_bw.  The
/// reflected wave should appear at port_2.re_bw.
#[test]
fn circulator_round_trip_to_port2() {
    let netlist = "\
.options enable_bidirectional=1
.optical_port p0
.optical_port p1
.optical_port p2
Vp0_re p0_re_fw_0 0 DC 1.0
Vp0_im p0_im_fw_0 0 DC 0.0
Vp0_wl p0_wl_0    0 DC 1.55e-6
* External feedback: tie port_1.fw to port_1.bw (perfect reflection).
* Use a 0 V source as a wire-equality stamp.
Vrefl_re p1_re_fw_0 p1_re_bw_0 DC 0.0
Vrefl_im p1_im_fw_0 p1_im_bw_0 DC 0.0
Xcirc p0 p1 p2 fc_circulator
.op
.end
";
    let net = parse_spice(netlist).unwrap();
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP");
    // Light path: port_0.fw → port_1.bw (= 1.0) → external short → port_1.fw
    // → port_2.bw (= port_1.fw = 1.0).
    let p2_re_bw = r.node_voltage("p2_re_bw_0").unwrap();
    assert!(
        (p2_re_bw - 1.0).abs() < 1e-6,
        "after perfect reflection at port 1, port_2.re_bw should be 1.0; got {p2_re_bw}"
    );
}

/// Circulator must refuse to instantiate without bidir mode.
#[test]
#[should_panic(expected = "fc_circulator requires bidirectional propagation")]
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
    let _ = dc_op_nr_with_registry(&net, &DeviceRegistry::new());
}
