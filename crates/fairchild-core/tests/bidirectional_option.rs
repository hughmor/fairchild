//! Regression tests for the `enable_bidirectional` option (C.1 plumbing only).
//!
//! At this commit the flag flows through `.options` → parser → `BundlePort`
//! wire emission AND `SimOptions::bidirectional_propagation` → `SimContext`.
//! Photonic devices haven't been refactored to USE the 5-wire bundles yet
//! (C.2), so any netlist that *enables* bidir and references a port from a
//! device will fail to instantiate.  These tests just verify the plumbing.

use fairchild_core::SimOptions;
use fairchild_parser::parse_spice;

/// Default state: bidir off, ports emit 3 wires.
#[test]
fn default_optical_port_emits_3_wires() {
    let net = parse_spice(
        "\
.optical_port ch0
.optical_port bus 4
.op
.end
",
    )
    .unwrap();
    let ch0 = net.bundle_ports.iter().find(|p| p.name == "ch0").unwrap();
    let bus = net.bundle_ports.iter().find(|p| p.name == "bus").unwrap();
    assert_eq!(
        ch0.kind,
        fairchild_parser::BundleKind::Optical {
            bidirectional: false
        }
    );
    assert_eq!(ch0.wires_per_channel(), 3);
    let names = ch0.wires_for_channel(0);
    assert_eq!(names, vec!["ch0_re_0", "ch0_im_0", "ch0_wl_0"]);
    // 4-channel bus → 12 wires total.
    assert_eq!(bus.all_wires().len(), 12);

    let opts = SimOptions::from_netlist(&net);
    assert!(!opts.bidirectional_propagation);
}

/// `.options enable_bidirectional=1` flips the parser into 5-wire mode and
/// renames the wires with `_fw_` / `_bw_` suffixes.
#[test]
fn enable_bidirectional_emits_5_wires() {
    let net = parse_spice(
        "\
.options enable_bidirectional=1
.optical_port ch0
.optical_port bus 2
.op
.end
",
    )
    .unwrap();
    let ch0 = net.bundle_ports.iter().find(|p| p.name == "ch0").unwrap();
    let bus = net.bundle_ports.iter().find(|p| p.name == "bus").unwrap();
    assert_eq!(
        ch0.kind,
        fairchild_parser::BundleKind::Optical {
            bidirectional: true
        }
    );
    assert_eq!(ch0.wires_per_channel(), 5);
    let names = ch0.wires_for_channel(0);
    assert_eq!(
        names,
        vec![
            "ch0_re_fw_0",
            "ch0_im_fw_0",
            "ch0_re_bw_0",
            "ch0_im_bw_0",
            "ch0_wl_0",
        ]
    );
    // 2-channel bus → 10 wires total.
    assert_eq!(bus.all_wires().len(), 10);

    let opts = SimOptions::from_netlist(&net);
    assert!(opts.bidirectional_propagation);
    assert_eq!(opts.sim_context().wires_per_channel(), 5);
}

/// Aliases `bidirectional=1` and `bidirectional_propagation=1` also work.
#[test]
fn enable_bidirectional_aliases() {
    for keyword in [
        "enable_bidirectional",
        "bidirectional",
        "bidirectional_propagation",
    ] {
        let src = format!(".options {keyword}=1\n.optical_port ch0\n.op\n.end\n");
        let net = parse_spice(&src).unwrap();
        let ch0 = net.bundle_ports.iter().find(|p| p.name == "ch0").unwrap();
        assert_eq!(
            ch0.kind,
            fairchild_parser::BundleKind::Optical {
                bidirectional: true
            },
            "{keyword} should enable bidir"
        );
    }
}

/// Bidir off explicitly via `enable_bidirectional=0` keeps 3-wire mode.
#[test]
fn enable_bidirectional_zero_stays_off() {
    let net = parse_spice(
        "\
.options enable_bidirectional=0
.optical_port ch0
.op
.end
",
    )
    .unwrap();
    let ch0 = net.bundle_ports.iter().find(|p| p.name == "ch0").unwrap();
    assert_eq!(
        ch0.kind,
        fairchild_parser::BundleKind::Optical {
            bidirectional: false
        }
    );
}
