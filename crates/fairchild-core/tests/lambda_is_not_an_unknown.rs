//! λ is a label, so it is not a row.
//!
//! `lambda_is_a_label.rs` pins that a solved λ only ever reads a source's
//! wavelength; `lambda_resolution.rs` pins that resolving it from declared
//! routing reproduces what the matrix used to solve. This file pins the
//! consequence: the matrix no longer carries λ at all, and the physics that
//! depends on λ still comes out right.
//!
//! Two halves, because either alone would pass against a bug:
//!
//! - **Structural.** No λ net is an MNA row, and every one of them is still a
//!   *net* — `V(a_wl_0)` answers with the resolved wavelength, because the
//!   X-line ABI is positional and decks and Verilog-A models address those
//!   wires by name.
//! - **Absolute.** The propagation phase through a short waveguide is checked
//!   against `2π·n_eff(λ)·L/λ` computed here, in the test, from the laser's
//!   parameter. Nothing in the solver is consulted for the expected value. A λ
//!   that failed to reach the segment would fall back to the band centre, which
//!   this deck deliberately is not at.

use fairchild_core::{dc_op_nr_with_registry, CircuitTopology, DeviceRegistry, SimOptions};
use fairchild_parser::parse_spice;

/// 500 × 220 nm SOI strip, the `fc_waveguide` defaults.
const N_EFF: f64 = 2.445;
const N_G: f64 = 4.19;
/// Short enough that φ < 2π, so a wrong λ cannot alias onto the right answer.
const L_UM: f64 = 0.1;
/// Well away from the 1550 nm band centre an unreached λ would bootstrap to.
const LAMBDA_NM: f64 = 1310.0;

fn deck() -> String {
    format!(
        "\
* one laser, one short waveguide
.optical_port a
.optical_port b
Xl a fc_cw_laser power_mW=1.0 wavelength_nm={LAMBDA_NM}
Xw a b fc_waveguide L_um={L_UM} n_eff={N_EFF} n_g={N_G} alpha_dB_cm=0
.op
"
    )
}

#[test]
fn no_lambda_net_is_a_matrix_row_and_all_of_them_are_still_nets() {
    let net = parse_spice(&deck()).expect("deck parses");
    let reg = DeviceRegistry::new();
    let ctx = SimOptions::default().sim_context();
    let topo = CircuitTopology::build_resolved(&net, &ctx, &reg);

    let lambda_nets: Vec<&String> = net
        .optical_nets
        .iter()
        .filter(|n| fairchild_parser::is_lambda_wire(n))
        .collect();
    assert_eq!(lambda_nets.len(), 2, "expected a λ wire per port");

    let x = vec![0.0; topo.size];
    for name in &lambda_nets {
        assert!(
            !topo.node_index.contains_key(name.as_str()),
            "{name} is still an MNA row"
        );
        let probed = topo
            .node_voltage(name, &x)
            .unwrap_or_else(|e| panic!("{name} is no longer probeable: {e}"));
        assert!(
            (probed - LAMBDA_NM * 1e-9).abs() < 1e-18,
            "V({name}) = {probed:e}, want {:e}",
            LAMBDA_NM * 1e-9
        );
    }

    // The field wires are still rows — this is about λ, not about optics.
    for port in ["a", "b"] {
        for wire in ["re", "im"] {
            let n = format!("{port}_{wire}_0");
            assert!(
                topo.node_index.contains_key(&n),
                "{n} should still be an unknown"
            );
        }
    }
}

/// The propagation phase must be the one the laser's wavelength implies.
///
/// This is the half that cannot be satisfied by bookkeeping: `A_out/A_in =
/// exp(−j·2π·n_eff(λ)·L/λ)` with `n_eff(λ)` the segment's linearised
/// dispersion, all of it evaluated here from the deck's numbers. At 1310 nm the
/// phase is 1.0956 rad; at the 1550 nm band centre a lost label would land on,
/// it is 0.9188 rad — 0.18 rad apart, far outside the 1e-9 this asserts.
#[test]
fn the_phase_is_the_one_the_sources_wavelength_implies() {
    let net = parse_spice(&deck()).expect("deck parses");
    let r = dc_op_nr_with_registry(&net, &DeviceRegistry::new()).expect("DC OP converges");

    let (a_re, a_im) = (
        r.node_voltage("a_re_0").unwrap(),
        r.node_voltage("a_im_0").unwrap(),
    );
    let (b_re, b_im) = (
        r.node_voltage("b_re_0").unwrap(),
        r.node_voltage("b_im_0").unwrap(),
    );
    let amp_in = (a_re * a_re + a_im * a_im).sqrt();
    assert!(amp_in > 0.0, "the laser is dark; nothing is being measured");

    let lambda = LAMBDA_NM * 1e-9;
    // `OpticalSegment` linearises dispersion about `wl_ref_m`, which defaults to
    // the band centre: n_eff(λ) = n_eff + (λ − λ_ref)·(n_eff − n_g)/λ_ref.
    let wl_ref = 1.55e-6;
    let n_eff_lam = N_EFF + (lambda - wl_ref) * (N_EFF - N_G) / wl_ref;
    let want = 2.0 * std::f64::consts::PI * n_eff_lam * (L_UM * 1e-6) / lambda;

    // A_out = A_in · exp(−jφ) ⇒ φ = −arg(A_out / A_in).
    let got = -((b_im * a_re - b_re * a_im).atan2(b_re * a_re + b_im * a_im));
    assert!(
        (got - want).abs() < 1e-9,
        "phase through the waveguide = {got:.6} rad, but the laser's \
         {LAMBDA_NM} nm implies {want:.6} rad — the wavelength that reached the \
         segment is not the one the source emits"
    );
    // Lossless, so the amplitude must survive intact: a phase check alone would
    // pass on a segment that had quietly gone dark.
    let amp_out = (b_re * b_re + b_im * b_im).sqrt();
    assert!(
        (amp_out / amp_in - 1.0).abs() < 1e-9,
        "α=0 must be lossless: |A_out|/|A_in| = {}",
        amp_out / amp_in
    );
}
