//! The add-drop worked example, checked against the cell it replaces:
//! `examples/verilog_a/models/mrm_addrop.va` vs
//! `examples/photonic/pcells/mrm.sp`.
//!
//! An agreement invariant between two subsystems cannot catch a fault common to
//! both — so this is deliberately not that. The Verilog-A model is one closed
//! form for an add-drop ring; the cell is four native devices (two couplers, two
//! phase-shifter arcs) that the solver wires together and knows nothing about
//! rings. Nothing is shared but the parameter values, which is the point: the
//! cell is the absolute anchor, and every default in the header was copied from
//! its card. Change one and this file fails.
//!
//! What is NOT checked here is anything the cell cannot do. Its `r_th` is 0, so
//! it has no self-heating and no thermal dynamics; the tests below run at a
//! microwatt precisely so that difference stays below the tolerance, and the
//! last test asserts the difference exists rather than pretending it does not.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_osdi::{load_libraries_with_widths, VaOptions};
use fairchild_parser::{instantiated_widths, parse_spice_with_arity, PermissiveArity};
use std::collections::{BTreeMap, BTreeSet};

use crate::common;

/// The discrete cell, verbatim. Inlined rather than `.include`d so the test
/// depends on the file the giona chip actually instantiates.
const CELL: &str = include_str!("../../../../examples/photonic/pcells/mrm.sp");

fn model_dir() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/verilog_a/models")
}

fn solve(deck: &str) -> fairchild_core::NrResult {
    let probe = parse_spice_with_arity(deck, &PermissiveArity).expect("probe pass parses");
    let widths: BTreeMap<String, BTreeSet<usize>> = instantiated_widths(&probe);
    let mut reg = DeviceRegistry::new();
    // The cell's arcs come from a `.model … fc_pn_th_ps LEVEL=4` card declared
    // inside the subckt, so the card has to reach the registry before anything
    // can place them. `DeviceRegistry::new()` knows the `fc_*` families and
    // nothing a deck named.
    reg.register_builtin_models(&probe.models);
    load_libraries_with_widths(
        &[],
        &[model_dir()
            .join("mrm_addrop.va")
            .to_string_lossy()
            .to_string()],
        None,
        &VaOptions {
            cache_dir: Some(common::test_cache_dir()),
            include_dirs: vec![model_dir()],
            ..VaOptions::from_env()
        },
        &mut reg,
        &widths,
        3,
    )
    .expect("the example model generates and compiles");
    let net = parse_spice_with_arity(deck, &reg).expect("the deck parses against the model");
    dc_op_nr_with_registry(&net, &reg).expect("DC OP converges")
}

/// The instance line that puts the model in card-compatible mode: the four
/// linearisations on, the Soref-Bennett path off, and the cell's terminal
/// conventions restored.
///
/// The model's defaults are no longer the cell — it counts carriers once and
/// runs the depletion side through Soref-Bennett, which is sqrt-shaped where the
/// card is linear. That is a deliberate improvement and it is checked against
/// the May capture in `experiments/giona/va_mrm_model/va_vs_may_data.py`, not here. What is
/// checked here is that the improvement did not quietly break the port: with
/// these eleven overrides the model must still be the cell, exactly, and the
/// sabotage list in this file's header still applies through them.
const LEGACY: &str = "dn_dv=-3.62e-5 da_dv=3.29e-4 dn_dnc=-8.8e-28 \
                      dalpha_dnc=1.0212e-21 vol_dep=0 sb_dn_e=0 sb_dn_h=0 \
                      sb_da_e=0 sb_da_h=0 i_sat=1.0198e-7 vol_inj=2.7531e-17 \
                      c_j0=2.75e-13 tau_sweep=1";

/// A microwatt. Low enough that the Verilog-A model's self-heating — which the
/// cell has no way to produce — moves the resonance by under 0.1 pm, which is
/// under 0.3 % of transmission on the steepest flank of the notch.
const P_UW: f64 = 0.001;

/// The unbiased resonance of this ring, so the sampled wavelengths straddle the
/// notch instead of sitting on a flat unity baseline where every parameter in
/// the card would be free to be wrong.
const LAMBDA_RES_NM: f64 = 1549.674;

/// Both rings, one laser, one solve. `port` chooses which optical port the light
/// goes into — `in` for the through path, `ad` for the add path.
fn deck(port: &str, lambda_nm: f64, v_pn: f64, v_htr: f64, p_mw: f64, extra: &str) -> String {
    let add = port == "ad";
    let (d_bus, d_add) = if add {
        ("d_dk", "src")
    } else {
        ("src", "d_dk")
    };
    let (v_bus, v_add) = if add {
        ("v_dk", "src")
    } else {
        ("src", "v_dk")
    };
    format!(
        "{CELL}\n\
         .optical_port src\n\
         .optical_port d_dk\n.optical_port d_th\n.optical_port d_dr\n\
         .optical_port v_dk\n.optical_port v_th\n.optical_port v_dr\n\
         XL src fc_cw_laser power_mW={p_mw} wavelength_nm={lambda_nm:.6}\n\
         Xd {d_bus} d_th {d_add} d_dr pn 0 htr 0 mrm\n\
         Xv {v_bus} v_th {v_add} v_dr pn 0 htr 0 tv mrm_addrop {extra}\n\
         VPN pn 0 DC {v_pn}\n\
         VHT htr 0 DC {v_htr}\n\
         .op\n"
    )
}

fn power_w(r: &fairchild_core::NrResult, net: &str) -> f64 {
    let re = r.node_voltage(&format!("{net}_re_0")).unwrap();
    let im = r.node_voltage(&format!("{net}_im_0")).unwrap();
    re * re + im * im
}

/// Assert two powers agree to `tol` relative, on a floor of a millionth of the
/// input so a dark port cannot pass by being small.
fn agree(cell: f64, va: f64, tol: f64, what: &str) {
    let floor = P_UW * 1e-3 * 1e-6;
    let denom = cell.abs().max(floor);
    let rel = (va - cell).abs() / denom;
    assert!(
        rel < tol,
        "{what}: cell {cell:.6e} W, Verilog-A {va:.6e} W — {:.3} % apart, over {:.3} %",
        rel * 100.0,
        tol * 100.0
    );
}

/// Every mechanism the cell has, over the notch, on both output ports.
///
/// The wavelengths straddle the resonance and the biases exercise each drive in
/// turn: depletion (reverse), carrier injection (forward), and the heater. The
/// last is the one whose transcription is not a copy — the cell tunes on a
/// calibrated `p_pi_th`, the Verilog-A model on a solved temperature through
/// `r_th` — so a slip in `r_th = lambda / (2 dn_dt L p_pi_th)` shows up here and
/// nowhere else.
#[test]
fn the_verilog_a_ring_reproduces_the_discrete_cell() {
    if !common::have_compiler() {
        return;
    }
    let offsets = [-0.10, -0.04, 0.0, 0.04, 0.10];
    // `r_th` is the one parameter that is a match on one row and an addition on
    // every other. On the heater row it carries `p_pi_th` and must be its
    // shipped value or nothing tunes. On the bias rows it also lets the
    // junction's own dissipation warm the ring, which the cell — at `r_th=0` —
    // cannot do: at 0.85 V that is 60 uW, a 0.2 K rise and a 13 pm red shift,
    // worth about 3 dB on the flank. Muting it there is not making the test
    // pass, it is asking the question the cell is able to answer.
    //
    // The heater row also gets a looser tolerance, for a reason worth stating
    // rather than absorbing: `p_pi_th` is a phase per watt and is the same at
    // every wavelength, while `dn_dt * r_th` is an INDEX per watt and so tunes
    // fractionally harder the further the ring sits from `wl_ref_nm`. The two
    // agree at 1550 nm by construction and drift about a picometre per
    // nanometre of detuning. That is the more physical of the two — a
    // thermo-optic coefficient is an index — so it stays. At 1.5 V the ring has
    // walked 1315 pm, so the deficit is 0.8 pm, and 0.8 pm on the steep flank of
    // a 167 pm notch is a few per cent of transmission. Hence 5 %, which is not
    // a slack test: the same flank turns a 1 % error in `r_th` into 13 pm and
    // tens of per cent, so this still resolves `r_th` to about a tenth of one.
    let biases = [
        ("unbiased", 0.0, 0.0, "r_th=1e-6", 2e-3),
        ("depletion", -3.0, 0.0, "r_th=1e-6", 2e-3),
        ("injection", 0.85, 0.0, "r_th=1e-6", 2e-3),
        ("heater", 0.0, 1.5, "", 5e-2),
    ];
    for (what, v_pn, v_htr, extra, tol) in biases {
        // The heater walks the ring by P/p_pi_th half-FSRs; follow it, or every
        // sample lands on the baseline and the test passes on nothing.
        let walk = (v_htr * v_htr / 368.8) / 26.4e-3 * 11.38 / 2.0;
        for off in offsets {
            let lam = LAMBDA_RES_NM + walk + off;
            let r = solve(&deck(
                "in",
                lam,
                v_pn,
                v_htr,
                P_UW,
                &format!("{LEGACY} {extra}"),
            ));
            let at = format!("{what} at {lam:.3} nm");
            agree(
                power_w(&r, "d_th"),
                power_w(&r, "v_th"),
                tol,
                &format!("thru, {at}"),
            );
            agree(
                power_w(&r, "d_dr"),
                power_w(&r, "v_dr"),
                tol,
                &format!("drop, {at}"),
            );
        }
    }
}

/// The add port, which the through path never touches.
///
/// A ring with two identical couplers is symmetric: light in on `add` leaves by
/// `drop` with the through path's transfer, and by `thru` with the drop path's.
/// Getting the four-way sign convention wrong is invisible from the bus side —
/// the model would pass every test above and route the add port into the wrong
/// output, which on the giona chip is a weight going to the wrong neuron.
#[test]
fn light_into_the_add_port_lands_where_the_cell_puts_it() {
    if !common::have_compiler() {
        return;
    }
    for off in [-0.06, 0.0, 0.06] {
        let lam = LAMBDA_RES_NM + off;
        let r = solve(&deck("ad", lam, 0.0, 0.0, P_UW, LEGACY));
        agree(
            power_w(&r, "d_dr"),
            power_w(&r, "v_dr"),
            0.01,
            &format!("add->drop at {lam:.3} nm"),
        );
        agree(
            power_w(&r, "d_th"),
            power_w(&r, "v_th"),
            0.01,
            &format!("add->thru at {lam:.3} nm"),
        );
    }
}

/// Absorbed light heats the ring — the thing the cell cannot do.
///
/// This is the other half of the contract. The tests above are only meaningful
/// if the agreement is a property of the matched mechanisms and not of a
/// thermal path that was quietly disabled to make the numbers line up, so:
/// raise the input power and the two models must come APART, in the direction a
/// warmer ring goes, and by a margin the microwatt tests could not have hidden.
#[test]
fn the_solved_temperature_walks_the_resonance_and_the_cell_cannot_follow() {
    if !common::have_compiler() {
        return;
    }
    // On the blue flank of the notch, a red shift means more transmission.
    let lam = LAMBDA_RES_NM - 0.06;
    let cold = solve(&deck("in", lam, 0.0, 0.0, P_UW, ""));
    let hot = solve(&deck("in", lam, 0.0, 0.0, 2.0, ""));

    let cell_t = |r: &fairchild_core::NrResult| power_w(r, "d_th") / power_w(r, "src");
    let va_t = |r: &fairchild_core::NrResult| power_w(r, "v_th") / power_w(r, "src");

    // The cell is power-blind: its normalised transfer is the same at 1 uW and
    // 2 mW, to a part in a thousand.
    let cell_drift = (cell_t(&hot) - cell_t(&cold)).abs() / cell_t(&cold);
    assert!(
        cell_drift < 1e-3,
        "the discrete cell moved with power ({cell_drift:.2e}); it has r_th=0 and cannot"
    );

    // The Verilog-A ring is not, and moves the right way.
    let va_shift = (va_t(&hot) - va_t(&cold)) / va_t(&cold);
    assert!(
        va_shift > 0.05,
        "self-heating should push the ring off the blue flank by well over 5 % of \
         transmission at 2 mW; got {:.2} %",
        va_shift * 100.0
    );
}
