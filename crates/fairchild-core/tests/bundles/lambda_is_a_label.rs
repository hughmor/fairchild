//! λ is propagated, never computed — the premise the whole design rests on.
//!
//! This was written while λ was still an MNA unknown, to establish that every λ
//! row either was driven by a source or copied another λ row — that the solver
//! was doing label propagation and nothing else. That is what made it safe to
//! stop solving for λ; the rows are gone now (see `crate::lambda`,
//! `lambda_is_not_an_unknown.rs`) and `V(p_wl_0)` reads the resolved label.
//!
//! It is kept because the premise still has to hold, and nothing else checks
//! it: every wavelength in a converged deck must be one some source emits, or
//! the band centre an unreached port takes. If a future device genuinely
//! *computes* a wavelength — four-wave mixing, a Raman shift — resolution
//! cannot express it, and this is the test that will say so rather than
//! quietly reporting a pump's colour on an idler wire.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

/// Every λ wire in `deck` must read one of `sources` (nm), or the band-centre
/// bootstrap an undriven wire falls back to.
fn assert_every_lambda_is_a_source_label(label: &str, deck: &str, sources: &[f64]) {
    let net = parse_spice(deck).expect("deck parses");
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP converges");

    let mut checked = 0;
    for name in net.optical_nets.iter() {
        if !fairchild_parser::is_lambda_wire(name) {
            continue;
        }
        let Ok(v) = r.node_voltage(name) else {
            continue;
        };
        let nm = v * 1e9;
        checked += 1;
        let matches_source = sources.iter().any(|s| (nm - s).abs() < 1e-6);
        // 1550 nm is `lambda_center_m`, what an undriven λ wire bootstraps to.
        let bootstrap = (nm - 1550.0).abs() < 1e-6 || nm == 0.0;
        assert!(
            matches_source || bootstrap,
            "{label}: {name} solved to {nm} nm, which is neither a source \
             wavelength {sources:?} nor the band-centre bootstrap. If a device \
             now genuinely computes a wavelength, λ is no longer a label and \
             resolving it before the solve is unsound."
        );
    }
    assert!(
        checked > 0,
        "{label}: no λ wires found — the test proved nothing"
    );
}

#[test]
fn a_wdm_bus_propagates_only_its_sources_wavelengths() {
    let deck = "\
.optical_port c0
.optical_port c1
.optical_port c2
.optical_port bus 3
.optical_port out 3
Xl0 c0 fc_cw_laser power_mW=5 wavelength_nm=1546.12
Xl1 c1 fc_cw_laser power_mW=5 wavelength_nm=1550.00
Xl2 c2 fc_cw_laser power_mW=5 wavelength_nm=1553.88
Xmux bus c0 c1 c2 fc_mux
Xpn bus out a 0 fc_pn_ps_full L_um=400
Vb a 0 DC -1.0
.op
";
    assert_every_lambda_is_a_source_label("wdm bus", deck, &[1546.12, 1550.0, 1553.88]);
}

/// A ring fed through its ADD port: the case `LambdaSelect`'s latch exists for
/// — the tag arrives from the add side and has to work its way around the loop,
/// which is also the cycle a setup-time graph walk would have to handle. The
/// value that lands is still a source's.
#[test]
fn an_add_fed_ring_still_only_carries_source_wavelengths() {
    let deck = "\
.optical_port bus_in
.optical_port bus_thru
.optical_port ring_fwd
.optical_port arc1_out
.optical_port add_in
.optical_port ring_c
.optical_port drop_out
.optical_port ring_ret
Xladd add_in fc_cw_laser power_mW=1.0 wavelength_nm=1600
Xdc1 bus_in ring_ret bus_thru ring_fwd fc_dcoupler kappa_L=0.336
Xarc1 ring_fwd arc1_out fc_waveguide L_um=25 n_eff=2.4 n_g=2.4 alpha_dB_cm=1.0
Xdc2 arc1_out add_in ring_c drop_out fc_dcoupler kappa_L=0.336
Xarc2 ring_c ring_ret fc_waveguide L_um=25 n_eff=2.4 n_g=2.4 alpha_dB_cm=1.0
.op
";
    assert_every_lambda_is_a_source_label("add-fed ring", deck, &[1600.0]);
}

/// Two different wavelengths meeting at a 2×2: the routing decision that a
/// setup-time resolution would have to reproduce.
#[test]
fn a_two_input_coupler_carries_a_source_label_on_every_output() {
    let deck = "\
.optical_port p1
.optical_port p2
.optical_port t1
.optical_port t2
Xl1 p1 fc_cw_laser power_mW=3 wavelength_nm=1310
Xl2 p2 fc_cw_laser power_mW=0 wavelength_nm=1310
Xc p1 p2 t1 t2 fc_dcoupler kappa_l=0.3
.op
";
    assert_every_lambda_is_a_source_label("2x2 coupler", deck, &[1310.0]);
}
