//! The shipped PCell library (`examples/photonic/pcells/`) must keep working.
//!
//! These decks exercise the whole hierarchical path: `.include` of a cell file,
//! `{…}` parameter arithmetic, per-instance `.model` cards, and a subcircuit
//! that takes a whole WDM bus versus one that takes a single channel.

use fairchild_core::device_registry::DeviceRegistry;
use fairchild_core::newton::{dc_op_nr_with_registry, NrResult};
use fairchild_parser::parse_spice;
use std::path::PathBuf;

fn pcell(name: &str) -> String {
    let p: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "..",
        "..",
        "examples",
        "photonic",
        "pcells",
        name,
    ]
    .iter()
    .collect();
    p.canonicalize().unwrap_or(p).to_string_lossy().into_owned()
}

fn solve(netlist: &str) -> NrResult {
    let parsed = parse_spice(netlist).expect("deck should parse");
    // The per-instance cards a `.subckt` emits have to be registered like any
    // other model card — that is exactly what is being tested here.
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&parsed.models);
    dc_op_nr_with_registry(&parsed, &registry).expect("DC OP should converge")
}

fn power(r: &NrResult, port: &str, k: usize) -> f64 {
    let re = r.node_voltage(&format!("{port}_re_{k}")).unwrap();
    let im = r.node_voltage(&format!("{port}_im_{k}")).unwrap();
    re * re + im * im
}

/// One MRM PCell instance: parses, expands, converges, and conserves power
/// across its four ports.
fn mrm_deck(radius_m: f64, n_eff: f64, lambda_nm: f64, v_pn: f64, i_heat: f64) -> String {
    format!(
        "* mrm pcell test\n\
         .include {}\n\
         .optical_port pin\n.optical_port pth\n\
         .optical_port pad\n.optical_port pdr\n\
         Xl pin fc_cw_laser power_mW=1.0 wavelength_nm={lambda_nm}\n\
         Var pad_re 0 DC 0\nVai pad_im 0 DC 0\nVaw pad_wl 0 DC {:e}\n\
         Xr pin pth pad pdr vpn 0 hc 0 mrm radius={radius_m:e} n_eff={n_eff}\n\
         Vpn vpn 0 DC {v_pn}\nIhc 0 hc DC {i_heat:e}\n.op\n.end\n",
        pcell("mrm.sp"),
        lambda_nm * 1e-9,
    )
}

/// `{pi*radius}` must actually drive the geometry: two instances differing only
/// in radius have different free spectral ranges, so their resonance combs
/// diverge. Sweeping n_eff over one FSR must find a resonance for each.
#[test]
fn mrm_pcell_radius_sets_the_fsr() {
    let lambda = 1550.0;
    // One FSR of n_eff at fixed L is λ/L, so the scan window shrinks with radius.
    let find_resonances = |radius: f64| -> usize {
        let l = 2.0 * std::f64::consts::PI * radius;
        let span = lambda * 1e-9 / l;
        let n = 60;
        let drops: Vec<f64> = (0..n)
            .map(|i| {
                let n_eff = 2.2810 + span * (i as f64) / (n as f64);
                let r = solve(&mrm_deck(radius, n_eff, lambda, 0.0, 0.0));
                power(&r, "pdr", 0)
            })
            .collect();
        // Count interior local maxima that actually couple light out.
        let peak = drops.iter().cloned().fold(0.0_f64, f64::max);
        (1..n - 1)
            .filter(|&i| {
                drops[i] > drops[i - 1] && drops[i] > drops[i + 1] && drops[i] > 0.5 * peak
            })
            .count()
    };
    // Exactly one resonance per FSR window, for both radii.
    assert_eq!(find_resonances(8e-6), 1, "8 µm ring");
    assert_eq!(find_resonances(12e-6), 1, "12 µm ring");
}

/// Per-instance `.model` cards: two instances of the same cell with different
/// EO parameters must behave differently. Here one has the fitted injection
/// coefficient and the other has it switched off, at a forward bias where
/// injection dominates.
#[test]
fn mrm_pcell_instances_carry_independent_models() {
    let deck = format!(
        "* two rings, different dn_di\n\
         .include {}\n\
         .optical_port ain\n.optical_port ath\n.optical_port aad\n.optical_port adr\n\
         .optical_port bin\n.optical_port bth\n.optical_port bad\n.optical_port bdr\n\
         Xla ain fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
         Xlb bin fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
         Vaar aad_re 0 DC 0\nVaai aad_im 0 DC 0\nVaaw aad_wl 0 DC 1.55e-6\n\
         Vbar bad_re 0 DC 0\nVbai bad_im 0 DC 0\nVbaw bad_wl 0 DC 1.55e-6\n\
         Xra ain ath aad adr va 0 ha 0 mrm n_eff=2.2872 dn_di=3.99\n\
         Xrb bin bth bad bdr vb 0 hb 0 mrm n_eff=2.2872 dn_di=0\n\
         Vva va 0 DC 0.9\nVvb vb 0 DC 0.9\n\
         Iha 0 ha DC 0\nIhb 0 hb DC 0\n.op\n.end\n",
        pcell("mrm.sp"),
    );
    let r = solve(&deck);
    let (a, b) = (power(&r, "adr", 0), power(&r, "bdr", 0));
    assert!(
        (a - b).abs() > 0.05 * b.max(1e-12),
        "instances with dn_di=3.99 vs 0 gave nearly the same drop power \
         ({a:.6e} vs {b:.6e}) — per-instance cards are not taking effect"
    );
}

/// The source bank: one instance carries the whole 8-channel bus (its declared
/// port count matches the flattened width), each channel gets its own
/// wavelength, its laser power gates it, and its MZM drive gates it too.
#[test]
fn source_bank_pcell_drives_eight_channels() {
    let mut deck = format!(
        "* source bank test\n.include {}\n.optical_port src 8\n\
         Xsrc src d1 d2 d3 d4 d5 d6 d7 d8 0 source_bank p3=0\n",
        pcell("source_bank.sp"),
    );
    // ch0 on, ch1 off via drive = v_pi, ch2 off via zero laser power,
    // ch3 half-on at v_pi/2, rest on.
    for (k, v) in [0.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0].iter().enumerate() {
        deck.push_str(&format!("Vd{} d{} 0 DC {v}\n", k + 1, k + 1));
    }
    deck.push_str(".op\n.end\n");
    let r = solve(&deck);

    let p: Vec<f64> = (0..8).map(|k| power(&r, "src", k) / 1e-3).collect();
    assert!((p[0] - 1.0).abs() < 1e-9, "ch0 should be full on: {}", p[0]);
    assert!(p[1] < 1e-6, "ch1 driven at v_pi should be off: {}", p[1]);
    assert!(p[2] < 1e-12, "ch2 has zero laser power: {}", p[2]);
    assert!((p[3] - 0.5).abs() < 1e-9, "ch3 at v_pi/2 = half: {}", p[3]);
    for (k, pk) in p.iter().enumerate().skip(4) {
        assert!((pk - 1.0).abs() < 1e-9, "ch{k} should be full on: {pk}");
    }
    // Wavelengths must survive the mux, on the 100 GHz grid the cell defaults to.
    let want = [
        1546.12, 1546.92, 1547.72, 1548.51, 1549.32, 1550.12, 1550.92, 1551.72,
    ];
    for (k, w) in want.iter().enumerate() {
        let got = r.node_voltage(&format!("src_wl_{k}")).unwrap() * 1e9;
        assert!((got - w).abs() < 1e-6, "ch{k} λ = {got}, want {w}");
    }
}

/// A subckt whose port count matches neither the flattened bus width nor a
/// single channel is a port-count error naming the width it would need, not a
/// confusing failure after expansion.
#[test]
fn subckt_bundle_width_mismatch_names_the_required_width() {
    let err = parse_spice(
        "* wrong width\n\
         .subckt two_wires a b\n\
         R1 a b 1k\n\
         .ends\n\
         .optical_port bus 4\n\
         X1 bus two_wires\n\
         .op\n.end\n",
    )
    .unwrap_err();
    let msg = format!("{err}");
    assert!(msg.contains("declares 2 ports"), "{msg}");
    assert!(msg.contains("expand to 12 wires"), "{msg}");
    assert!(msg.contains("whole 4-channel bus"), "{msg}");
}

/// A single-channel subckt on a multi-channel bundle is refused.
///
/// This is the regression that motivated dropping replication: the parser used
/// to emit four copies of the cell, so the four copies' electrical ports all
/// landed on the same two nodes and drew four times the current, silently.
#[test]
fn single_channel_subckt_on_a_wide_bundle_is_refused() {
    let err = parse_spice(
        "* would have been replicated\n\
         .subckt cell oi_re oi_im oi_l anode cathode\n\
         R1 anode cathode 1k\n\
         .ends\n\
         .optical_port bus 4\n\
         X1 bus a k cell\n\
         .op\n.end\n",
    )
    .unwrap_err();
    let msg = format!("{err}");
    assert!(msg.contains("no WDM semantics"), "{msg}");
    assert!(msg.contains("bus(4 ch)"), "{msg}");
    assert!(msg.contains("fc_demux"), "{msg}");
}

/// The same cell on a 1-channel bundle still works — that is the case where
/// replication and flattening always agreed.
#[test]
fn single_channel_subckt_on_a_one_channel_bundle_still_expands() {
    let net = parse_spice(
        "* one channel\n\
         .subckt cell oi_re oi_im oi_l anode cathode\n\
         R1 anode cathode 1k\n\
         .ends\n\
         .optical_port bus 1\n\
         X1 bus a k cell\n\
         V1 a 0 DC 1\n\
         R2 k 0 1k\n\
         .op\n.end\n",
    )
    .expect("1-channel bundle should expand");
    // The cell's 1k in series with the external 1k across 1 V.
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&net.models);
    let r = dc_op_nr_with_registry(&net, &registry).expect("solve");
    let vk = r.node_voltage("k").unwrap();
    assert!((vk - 0.5).abs() < 1e-6, "V(k) = {vk}, want 0.5");
}
