//! A WDM channel extinguished to **exactly** zero, mid-transient.
//!
//! The other half of the "a wire that carries no light" family: `#32` is a wire
//! two devices drive, this is a wire nobody does. An optical input port left
//! unconnected has a KCL row with nothing in it, so `gmin` alone holds it — and
//! `gmin` is twelve orders below the ±1 the reading device writes into that
//! node's *column*. Partial pivoting rejects the diagonal, eliminates the row
//! against a coupling row, and the wire comes back off a `gmin`-sized pivot
//! carrying roundoff amplified by `1/gmin`.
//!
//! While the channel is lit that noise is far inside `tolerance`'s bound and
//! nothing shows. Extinguish it and the bound collapses onto the noise, and
//! Newton stops converging at any timestep — which is why eight orders of
//! "nearly dark" all solve and only the exact zero fails (issue #47).
//!
//! **What each test is worth.** The exact-zero assertion below is the
//! regression: it fails against the bug at this deck size (the wires read
//! ~1e-17 instead of 0), and it is the contract the fix establishes rather than
//! a downstream symptom of it. The non-convergence itself needs a much larger
//! deck to appear — it was reproduced and fixed against the chip netlist the
//! issue names, which is not in this tree, and no synthetic deck this suite
//! could carry reproduced it. So the power budget here is an anchor on the
//! answer, not a second regression: do not read it as covering the symptom.

use fairchild_core::{
    options::SimOptions, tran::IntegratorMode, tran_nr_with_registry_var_opts, DeviceRegistry,
};
use fairchild_parser::parse_spice;

const N: usize = 4;
const BLOCKS: usize = 4;
const P_MW: f64 = 30.0;
/// The extinguished channel's index.
const DARK_CH: usize = 2;
/// MZM extinction ratio, dB. The residual power in an "off" channel is
/// `P·10^(−ER/10)`, which is the anchor the off-state budget uses.
const ER_DB: f64 = 200.0;

/// N channels muxed onto one bus, down a 2 mm waveguide, split to `BLOCKS`
/// weight blocks — each with its **second input port left unconnected**, which
/// is the shape the failing deck has and every reduction that converged did
/// not. Channel `DARK_CH` is on from 5 ns and extinguished again at 25 ns.
///
/// `off` is the MZM drive at the edge: `v_pi` exactly darkens the channel.
fn deck(off: f64, w: f64) -> String {
    let mut s = String::from("* WDM bus -> splitter chain -> weight blocks with dark inputs\n");
    s += ".options method=gear\n";
    for k in 0..N {
        s += &format!(".optical_port ch{k}\n.optical_port mz{k}\n");
    }
    s += &format!(".optical_port bus {N}\n.optical_port wg {N}\n.optical_port stub {N}\n");
    for b in 0..BLOCKS {
        s += &format!(
            ".optical_port win{b} {N}\n.optical_port dark{b} {N}\n\
             .optical_port thru{b} {N}\n.optical_port drop{b} {N}\n"
        );
    }
    for b in 0..BLOCKS - 1 {
        s += &format!(".optical_port sp{b} {N}\n");
    }
    for k in 0..N {
        let wl = 1546.0 + 0.8 * k as f64;
        s += &format!("Xl{k} ch{k} fc_cw_laser power_mW={P_MW} wavelength_nm={wl}\n");
        s += &format!("Xm{k} ch{k} mz{k} d{k} 0 fc_mzm v_pi=1.0 e_r_db={ER_DB}\n");
        if k == DARK_CH {
            s += &format!("Vd{k} d{k} 0 PULSE({off} 0 5n 200p 200p 20n 2)\n");
        } else {
            s += &format!("Vd{k} d{k} 0 DC 0\n");
        }
    }
    s += "Xmux bus";
    for k in 0..N {
        s += &format!(" mz{k}");
    }
    s += " fc_mux\n";
    s += "Xwg bus wg fc_waveguide l_m=2e-3 n_g=4.2 alpha_dB_cm=2.0\n";
    let mut cur = String::from("wg");
    for b in 0..BLOCKS {
        let rest = if b + 1 == BLOCKS {
            "stub".to_string()
        } else {
            format!("sp{b}")
        };
        s += &format!("Xsp{b} {cur} win{b} {rest} fc_splitter\n");
        cur = rest;
    }
    for b in 0..BLOCKS {
        // `dark{b}` is declared and wired but driven by nothing at all.
        s += &format!("Xw{b} win{b} dark{b} thru{b} drop{b}");
        for k in 0..N {
            s += &format!(" W{b}_{k}");
        }
        s += " 0 fc_optical_2x2 w=0 dw_dv=1\n";
        for k in 0..N {
            let v = if k == DARK_CH { w } else { 0.0 };
            s += &format!("VW{b}_{k} W{b}_{k} 0 DC {v}\n");
        }
        s += &format!(
            "Xpt{b} thru{b} pt{b} 0 fc_photodetector responsivity=1.0 i_dark_a=0 r_shunt=1e12\n\
             Xpd{b} drop{b} pd{b} 0 fc_photodetector responsivity=1.0 i_dark_a=0 r_shunt=1e12\n\
             Rt{b} pt{b} 0 1k\nRd{b} pd{b} 0 1k\n"
        );
    }
    s += ".tran 200p 40n\n.end\n";
    s
}

fn run(off: f64, w: f64) -> fairchild_core::tran::TranResult {
    let net = parse_spice(&deck(off, w)).expect("deck should parse");
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&net.models);
    let opts = SimOptions {
        method: IntegratorMode::Gear,
        variable_step: true,
        ..SimOptions::default()
    };
    tran_nr_with_registry_var_opts(&net, 200e-12, 40e-9, &registry, &opts)
        .expect("transient should converge")
}

fn series<'a>(r: &'a fairchild_core::tran::TranResult, node: &str) -> &'a [f64] {
    r.node_voltages
        .get(node)
        .unwrap_or_else(|| panic!("no node {node}"))
}

/// **The regression.** A wire nothing drives reads exactly zero, at every
/// timepoint, on every channel — including the one that goes dark mid-run.
///
/// Exactly zero, not "small": while `gmin` alone held these rows they came back
/// at ~1e-17 here, and at ~1e-4 on the deck in the issue. The difference is
/// only how much fill-in the elimination had to work with, so a tolerance would
/// be a number with no meaning. The equation is `V = 0` and the answer is 0.
#[test]
fn a_wire_nothing_drives_reads_exactly_zero() {
    let r = run(1.0, 0.8);
    for b in 0..BLOCKS {
        for k in 0..N {
            for wire in ["re", "im"] {
                let name = format!("dark{b}_{wire}_{k}");
                let worst = series(&r, &name)
                    .iter()
                    .fold(0.0f64, |acc, v| acc.max(v.abs()));
                assert_eq!(
                    worst, 0.0,
                    "{name} is driven by nothing and must read exactly 0; \
                     peak |v| was {worst:e}"
                );
            }
        }
    }
}

/// The power budget across the extinction edge, against a hand-computed anchor.
///
/// `fc_optical_2x2` in weight mode is unitary, so a block passes on everything
/// that reaches it: summed over the bus, `P_thru + P_drop` is the launched
/// power. Only the dark channel carries a weight, so `P_drop − P_thru` is `w`
/// times *that channel's* power alone — the other three split 50/50 and cancel.
/// Extinguishing the channel therefore has to remove exactly one channel's
/// worth from the sum and take the difference to zero, leaving only the
/// modulator's finite extinction ratio, `P·10^(−ER/10)`.
///
/// This is an anchor on the answer, not a second regression — see the module
/// comment for why the non-convergence itself is not reproducible at this size.
#[test]
fn the_extinguished_channel_closes_its_power_budget() {
    const W: f64 = 0.8;
    let r = run(1.0, W);
    // Per channel into block 0: the launch, through the waveguide, halved by
    // the first 3 dB split (block b takes one arm of the b-th split).
    let wg = 10f64.powf(-2.0 * 0.2 / 10.0); // 2 mm at 2 dB/cm, power
    let p_ch = P_MW * 1e-3 * wg * 0.5;
    let floor = 10f64.powf(-ER_DB / 10.0);
    let (mut n_lit, mut n_dark) = (0, 0);
    let (pt, pd) = (series(&r, "pt0"), series(&r, "pd0"));
    for (i, &t) in r.time.iter().enumerate() {
        // 1 V across 1 kΩ at unit responsivity is 1 mW of optical power.
        let (p_thru, p_drop) = (pt[i].abs() * 1e-3, pd[i].abs() * 1e-3);
        let (sum, diff) = (p_thru + p_drop, p_drop - p_thru);
        if (10e-9..24e-9).contains(&t) {
            let want_sum = N as f64 * p_ch;
            assert!(
                (sum - want_sum).abs() / want_sum < 1e-6,
                "t={t:.3e}: P_thru+P_drop={sum:.6e} W, launched {want_sum:.6e} W"
            );
            assert!(
                (diff - W * p_ch).abs() / p_ch < 1e-6,
                "t={t:.3e}: P_drop-P_thru={diff:.6e} W, expected {:.6e} W",
                W * p_ch
            );
            n_lit += 1;
        }
        if t > 26e-9 {
            // One channel's worth gone, bar the modulator's own floor.
            let want_sum = (N as f64 - 1.0 + floor) * p_ch;
            assert!(
                (sum - want_sum).abs() / want_sum < 1e-6,
                "t={t:.3e}: P_thru+P_drop={sum:.6e} W after extinction, \
                 expected {want_sum:.6e} W"
            );
            assert!(
                diff.abs() <= W * p_ch * floor + p_ch * 1e-9,
                "t={t:.3e}: the weighted channel is dark, so P_drop-P_thru must \
                 collapse to the extinction floor; got {diff:.6e} W"
            );
            n_dark += 1;
        }
    }
    assert!(n_lit > 0 && n_dark > 0, "the run missed an edge");
}

/// The weights the issue swept, because the non-convergence was non-monotonic
/// in `w` — 0.7 sat between two failures, so a single passing value proved
/// nothing. Here they pin that the unit row pin did not disturb the answer at
/// any of them: exact extinction and one part in ten thousand short of it must
/// agree while the channel is on, which is where both always had an answer.
#[test]
fn every_swept_weight_agrees_between_exact_and_near_extinction() {
    for w in [0.05, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8] {
        let (exact, nearly) = (run(1.0, w), run(0.9999, w));
        let at = |r: &fairchild_core::tran::TranResult, node: &str| -> f64 {
            let i = r
                .time
                .iter()
                .enumerate()
                .min_by(|a, b| {
                    (a.1 - 15e-9)
                        .abs()
                        .partial_cmp(&(b.1 - 15e-9).abs())
                        .unwrap()
                })
                .map(|(i, _)| i)
                .unwrap();
            series(r, node)[i]
        };
        for node in ["pt0", "pd0"] {
            let (e, n) = (at(&exact, node), at(&nearly, node));
            assert!(
                (e - n).abs() < 1e-6,
                "w={w} {node}: exact={e:.9e} nearly={n:.9e}"
            );
        }
    }
}
