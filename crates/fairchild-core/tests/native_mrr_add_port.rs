//! An add-drop ring driven through its **add** port must work.
//!
//! The bus is dark; light enters on the add port only. Every two-input photonic
//! device used to mirror port a's λ wire onto both outputs, so with port a unlit
//! the ring's arcs saw λ = 0 — and `OpticalSegment::lambda_of` quietly
//! substitutes the band centre for an undriven λ wire, so the ring resonated
//! beautifully at the wrong wavelength. No error, no obviously wrong number.
//!
//! Two assertions, one structural and one physical:
//!   - the λ tag reaches every port downstream of the add input (it was 0),
//!   - the resonance lands where the add port's own wavelength puts it, not
//!     where the band centre does.

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

const LAMBDA_NM: f64 = 1550.0;
/// Round-trip length: two 25 µm arcs.
const ARC_UM: f64 = 25.0;

/// Add-drop ring. With `bus_laser = false` there is no source on the bus at
/// all — that is the case under test. A zero-power laser would not do: a laser
/// drives its λ wire regardless of power, which is exactly what was masking
/// this bug in decks that terminate every port with a source.
fn ring(n_eff: f64, bus_laser: bool) -> String {
    ring_at(n_eff, LAMBDA_NM, bus_laser)
}

fn ring_at(n_eff: f64, lambda_nm: f64, bus_laser: bool) -> String {
    let bus = if bus_laser {
        format!("Xlbus bus_in fc_cw_laser power_mW=1.0 wavelength_nm={lambda_nm}")
    } else {
        String::new()
    };
    format!(
        "\
* add-drop ring fed through the add port
.optical_port bus_in
.optical_port bus_thru
.optical_port ring_fwd
.optical_port arc1_out
.optical_port add_in
.optical_port ring_c
.optical_port drop_out
.optical_port ring_ret

Xladd add_in fc_cw_laser power_mW=1.0 wavelength_nm={lambda_nm}
{bus}

* a=bus_in  b=ring_ret  c=bus_thru  d=ring_fwd
Xdc1 bus_in ring_ret bus_thru ring_fwd fc_dcoupler kappa_L=0.336
Xarc1 ring_fwd arc1_out fc_waveguide L_um={ARC_UM} n_eff={n_eff} n_g={n_eff} alpha_dB_cm=1.0
* a=arc1_out  b=add_in  c=ring_c  d=drop_out
Xdc2 arc1_out add_in ring_c drop_out fc_dcoupler kappa_L=0.336
Xarc2 ring_c ring_ret fc_waveguide L_um={ARC_UM} n_eff={n_eff} n_g={n_eff} alpha_dB_cm=1.0

.op
"
    )
}

fn solve(deck: &str) -> fairchild_core::newton::NrResult {
    let parsed = parse_spice(deck).expect("deck should parse");
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&parsed.models);
    dc_op_nr_with_registry(&parsed, &registry).expect("DC OP should converge")
}

fn power(r: &fairchild_core::newton::NrResult, port: &str) -> f64 {
    let re = r.node_voltage(&format!("{port}_re_0")).unwrap();
    let im = r.node_voltage(&format!("{port}_im_0")).unwrap();
    re * re + im * im
}

/// The λ tag must reach every output reachable from the add port.
#[test]
fn add_port_lambda_reaches_the_ring_and_the_outputs() {
    let r = solve(&ring(2.4, false));
    let want = LAMBDA_NM * 1e-9;
    for port in ["drop_out", "ring_c", "arc1_out", "ring_fwd", "bus_thru"] {
        let got = r.node_voltage(&format!("{port}_wl_0")).unwrap();
        assert!(
            (got - want).abs() < 1e-15,
            "{port} λ = {got:e}, want {want:e} — the add port's wavelength tag \
             did not propagate, so the ring will resonate at the band centre"
        );
    }
}

/// The resonance must land where the add port's *own* wavelength puts it.
///
/// This is the assertion that bites. Losing the λ tag does not stop the ring
/// resonating — `OpticalSegment::lambda_of` bootstraps an undriven λ wire to
/// the band centre — so the cavity still looks like a cavity. It just resonates
/// at 1550 nm when the light is at 1600 nm, which for a WDM chain (ring N's
/// drop feeding ring N−1's add) collapses every channel onto the band centre
/// while producing entirely plausible-looking output.
#[test]
fn add_port_resonance_lands_at_the_add_wavelength() {
    // Resonance when n_eff·L = m·λ, so in n_eff the comb spacing is λ/L.
    let l_m = 2.0 * ARC_UM * 1e-6;
    let at_1600 = {
        let m = (2.40 * l_m / (1600.0 * 1e-9)).round();
        m * 1600.0 * 1e-9 / l_m
    };
    let at_1550 = {
        let m = (2.40 * l_m / (1550.0 * 1e-9)).round();
        m * 1550.0 * 1e-9 / l_m
    };
    assert!(
        (at_1600 - at_1550).abs() > 5e-3,
        "test is not discriminating: predictions {at_1600} vs {at_1550}"
    );

    // Scan n_eff around both predictions and find where add→thru peaks.
    let n = 160;
    let (lo, hi) = (2.375, 2.415);
    let mut best = (f64::MIN, 0.0);
    for i in 0..n {
        let n_eff = lo + (hi - lo) * (i as f64) / ((n - 1) as f64);
        let p = power(&solve(&ring_at(n_eff, 1600.0, false)), "bus_thru");
        if p > best.0 {
            best = (p, n_eff);
        }
    }
    let found = best.1;
    let step = (hi - lo) / ((n - 1) as f64);
    assert!(
        (found - at_1600).abs() < 2.0 * step,
        "resonance at n_eff = {found:.5}; 1600 nm predicts {at_1600:.5} and the \
         1550 nm band centre predicts {at_1550:.5} — the add port's wavelength \
         is not reaching the ring"
    );
}

/// Lighting the bus as well must not change which λ the ring uses — port a
/// wins when both inputs carry light, and both carry the same λ anyway.
#[test]
fn bus_and_add_together_keep_the_same_lambda() {
    let r = solve(&ring(2.4, true));
    let want = LAMBDA_NM * 1e-9;
    let got = r.node_voltage("drop_out_wl_0").unwrap();
    assert!((got - want).abs() < 1e-15, "drop λ = {got:e}");
}
