//! `fc_awgr` — N×N cyclic arrayed-waveguide grating router.
//!
//! The routing convention under test: input port `i` channel `k` leaves on
//! output port `(i + k) mod N`, still in channel slot `k`.

use fairchild_core::device_registry::DeviceRegistry;
use fairchild_core::newton::dc_op_nr_with_registry;
use fairchild_parser::parse_spice;

/// Drive every channel of every input port of an N×N router with a distinct,
/// recognisable amplitude, so each output wire can be traced back to exactly
/// one source. Amplitude for (port i, channel k) is `1 + i + 10·k`.
fn sources(n: usize, lambda_nm: &[f64]) -> String {
    let mut s = String::new();
    for i in 0..n {
        for (k, lam) in lambda_nm.iter().enumerate().take(n) {
            let amp = 1.0 + i as f64 + 10.0 * k as f64;
            s.push_str(&format!("V{i}_{k}r in{i}_{k}_re 0 DC {amp}\n"));
            s.push_str(&format!("V{i}_{k}i in{i}_{k}_im 0 DC 0.0\n"));
            s.push_str(&format!("V{i}_{k}w in{i}_{k}_wl 0 DC {:e}\n", lam * 1e-9));
        }
    }
    s
}

/// Flat net list for one side of the router: `pfx0_0_re pfx0_0_im pfx0_0_wl …`.
fn wires(pfx: &str, n: usize) -> String {
    let mut s = String::new();
    for p in 0..n {
        for k in 0..n {
            s.push_str(&format!(" {pfx}{p}_{k}_re {pfx}{p}_{k}_im {pfx}{p}_{k}_wl"));
        }
    }
    s
}

/// A 100 GHz grid anchored at 1550 nm, expressed in nm.
fn grid_nm(n: usize) -> Vec<f64> {
    const C: f64 = 299_792_458.0;
    let f0 = C / 1550e-9;
    (0..n).map(|k| C / (f0 + k as f64 * 100e9) * 1e9).collect()
}

fn solve(netlist: &str) -> fairchild_core::newton::NrResult {
    let net = parse_spice(netlist).expect("parse");
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);
    dc_op_nr_with_registry(&net, &reg).expect("DC OP should converge")
}

/// **The convention test.** An ideal `fc_awgr` must be indistinguishable from
/// the same permutation built out of `fc_demux` + `fc_mux`, which are already
/// golden-tested primitives. If the routing convention is ever redefined, this
/// is what catches it.
#[test]
fn ideal_router_equals_the_equivalent_demux_mux_permutation() {
    let n = 4;
    let lam = grid_nm(n);

    // (a) the device.
    let mut deck = format!("* ideal AWGR\n{}", sources(n, &lam));
    deck.push_str(&format!(
        "Xr{}{} fc_awgr\n.op\n.end\n",
        wires("in", n),
        wires("out", n)
    ));
    let awgr = solve(&deck);

    // (b) the same thing out of primitives: split each input port into its N
    // single-channel bundles, then recombine so that output port j takes
    // channel k from input port (j − k) mod N.
    let mut deck = format!("* demux/mux permutation\n{}", sources(n, &lam));
    for i in 0..n {
        deck.push_str(&format!("Xd{i}"));
        for k in 0..n {
            deck.push_str(&format!(" in{i}_{k}_re in{i}_{k}_im in{i}_{k}_wl"));
        }
        for k in 0..n {
            deck.push_str(&format!(" c{i}_{k}_re c{i}_{k}_im c{i}_{k}_wl"));
        }
        deck.push_str(" fc_demux\n");
    }
    for j in 0..n {
        deck.push_str(&format!("Xm{j}"));
        for k in 0..n {
            deck.push_str(&format!(" ref{j}_{k}_re ref{j}_{k}_im ref{j}_{k}_wl"));
        }
        for k in 0..n {
            let src = (j + n - k) % n;
            deck.push_str(&format!(" c{src}_{k}_re c{src}_{k}_im c{src}_{k}_wl"));
        }
        deck.push_str(" fc_mux\n");
    }
    deck.push_str(".op\n.end\n");
    let refr = solve(&deck);

    for j in 0..n {
        for k in 0..n {
            for w in ["re", "im", "wl"] {
                let got = awgr.node_voltage(&format!("out{j}_{k}_{w}")).unwrap();
                let want = refr.node_voltage(&format!("ref{j}_{k}_{w}")).unwrap();
                assert!(
                    (got - want).abs() < 1e-12,
                    "out{j}_{k}_{w}: awgr={got} vs demux/mux={want}"
                );
            }
        }
    }
}

/// The cyclic shift itself, read straight off the amplitudes.
#[test]
fn cyclic_routing_sends_input_i_channel_k_to_output_i_plus_k() {
    let n = 4;
    let lam = grid_nm(n);
    let mut deck = format!("* cyclic routing\n{}", sources(n, &lam));
    deck.push_str(&format!(
        "Xr{}{} fc_awgr\n.op\n.end\n",
        wires("in", n),
        wires("out", n)
    ));
    let r = solve(&deck);
    for i in 0..n {
        for k in 0..n {
            let j = (i + k) % n;
            let got = r.node_voltage(&format!("out{j}_{k}_re")).unwrap();
            let want = 1.0 + i as f64 + 10.0 * k as f64;
            assert!(
                (got - want).abs() < 1e-12,
                "in{i} ch{k} → out{j} ch{k}: got {got}, want {want}"
            );
        }
    }
    // Every output port carries one wavelength from every input port, so all
    // N² output wires are accounted for and none is doubly driven: the total
    // power out equals the total power in.
    let p_in: f64 = (0..n)
        .flat_map(|i| (0..n).map(move |k| (1.0 + i as f64 + 10.0 * k as f64).powi(2)))
        .sum();
    let p_out: f64 = (0..n)
        .flat_map(|j| (0..n).map(move |k| (j, k)))
        .map(|(j, k)| {
            let re = r.node_voltage(&format!("out{j}_{k}_re")).unwrap();
            let im = r.node_voltage(&format!("out{j}_{k}_im")).unwrap();
            re * re + im * im
        })
        .sum();
    assert!(
        (p_out / p_in - 1.0).abs() < 1e-12,
        "ideal mode must be lossless: {p_out} vs {p_in}"
    );
}

/// Gauss mode at the grid centres: the routed path takes exactly the specified
/// insertion loss, and the unrouted ones sit on the crosstalk floors.
#[test]
fn gauss_mode_applies_insertion_loss_and_crosstalk_floors() {
    let n = 4;
    let lam = grid_nm(n);
    let mut deck = format!("* AWGR with a real passband\n{}", sources(n, &lam));
    deck.push_str(&format!(
        "Xr{}{} fc_awgr df_ghz=100 fwhm_ghz=40 il_db=3 xt_adj_db=-30 xt_bg_db=-40\n.op\n.end\n",
        wires("in", n),
        wires("out", n)
    ));
    let r = solve(&deck);
    let il = 10f64.powf(-3.0 / 20.0);
    // Routed path: input 2 channel 1 → output 3. Its amplitude is the source
    // amplitude times the field insertion loss, plus the crosstalk that the
    // other three inputs' channel-1 light contributes to that same slot.
    let got = r.node_voltage("out3_1_re").unwrap();
    let mut want = 0.0;
    for i in 0..n {
        let src = 1.0 + i as f64 + 10.0; // channel 1 amplitudes
        let m = (3 + n - i) % n; // this pair's passband index
        let off = ((1i64 - m as i64) + n as i64) % n as i64;
        let off = if off > n as i64 / 2 {
            off - n as i64
        } else {
            off
        };
        want += src
            * il
            * match off.abs() {
                0 => 1.0,
                1 => 1e-3f64.sqrt(),
                _ => 1e-4f64.sqrt(),
            };
    }
    assert!(
        (got / want - 1.0).abs() < 1e-9,
        "out3_1_re: got {got}, want {want}"
    );
    // Sanity on the scale: the routed term dominates. Not by the full 30 dB —
    // three other inputs leak in with amplitudes of their own, which is exactly
    // the coherent accumulation this device exists to show.
    let routed = (1.0 + 2.0 + 10.0) * il;
    assert!((got / routed - 1.0).abs() < 0.1, "{got} vs {routed}");
}

/// A laser detuned off the grid must lose power on the passband skirt. This is
/// the effect that a static-coefficient model *does* capture exactly, and the
/// reason gauss mode is worth having over ideal.
#[test]
fn detuning_a_laser_costs_it_the_passband_skirt() {
    const C: f64 = 299_792_458.0;
    let n = 2;
    let f0 = C / 1550e-9;
    // Only input port 0 is lit, so out0 ch0 carries one term and nothing else:
    // with a second port driven, its Gaussian tail also lands in that slot
    // (coherently — which is the point of the model, but not of this test).
    let detuned = C / (f0 + 20e9) * 1e9; // half a FWHM off channel 0
    let mut deck = String::from("* detuned laser\n");
    for k in 0..n {
        let l = if k == 0 {
            detuned
        } else {
            C / (f0 + 100e9) * 1e9
        };
        deck.push_str(&format!("V0_{k}r in0_{k}_re 0 DC 1.0\n"));
        deck.push_str(&format!("V0_{k}i in0_{k}_im 0 DC 0.0\n"));
        deck.push_str(&format!("V0_{k}w in0_{k}_wl 0 DC {:e}\n", l * 1e-9));
    }
    deck.push_str(&format!(
        "Xr{}{} fc_awgr df_ghz=100 fwhm_ghz=40 il_db=0 xt_adj_db=-300 xt_bg_db=-300\n.op\n.end\n",
        wires("in", n),
        wires("out", n)
    ));
    let r = solve(&deck);
    // Half a FWHM off centre is the −3 dB point: √0.5 in field.
    let got = r.node_voltage("out0_0_re").unwrap();
    assert!(
        (got - 0.5f64.sqrt()).abs() < 1e-9,
        "half-FWHM detuning should cost 3 dB: got {got}"
    );
    // On-grid channel 1 is untouched, so the penalty is the detuning and not
    // some global scaling.
    let on_grid = r.node_voltage("out1_1_re").unwrap();
    assert!((on_grid - 1.0).abs() < 1e-6, "on-grid channel: {on_grid}");
}

/// Output λ tags mirror the input comb by default, so a detuned laser stays
/// detuned for whatever resonant device sits downstream.
#[test]
fn output_wavelength_tags_follow_the_input_comb() {
    let n = 2;
    let lam = vec![1549.0, 1551.0];
    let mut deck = format!("* λ tags\n{}", sources(n, &lam));
    deck.push_str(&format!(
        "Xr{}{} fc_awgr\n.op\n.end\n",
        wires("in", n),
        wires("out", n)
    ));
    let r = solve(&deck);
    for j in 0..n {
        for (k, want) in lam.iter().enumerate() {
            let got = r.node_voltage(&format!("out{j}_{k}_wl")).unwrap();
            assert!(
                (got - want * 1e-9).abs() < 1e-15,
                "out{j}_{k}_wl: got {} nm, want {want} nm",
                got * 1e9,
            );
        }
    }
}

/// A demux is this device with N−1 inputs dark; the crosstalk that a
/// single-channel-output `fc_demux` structurally cannot represent shows up
/// here, in the slots the intended channel does not occupy.
#[test]
fn one_input_lit_behaves_as_a_demux_with_crosstalk() {
    let n = 4;
    let lam = grid_nm(n);
    // Only input port 0 is driven; the other three ports are left floating.
    let mut deck = String::from("* AWGR as a demux\n");
    for (k, l) in lam.iter().enumerate() {
        deck.push_str(&format!("V0_{k}r in0_{k}_re 0 DC 1.0\n"));
        deck.push_str(&format!("V0_{k}i in0_{k}_im 0 DC 0.0\n"));
        deck.push_str(&format!("V0_{k}w in0_{k}_wl 0 DC {:e}\n", l * 1e-9));
    }
    deck.push_str(&format!(
        "Xr{}{} fc_awgr df_ghz=100 fwhm_ghz=40 xt_adj_db=-25 xt_bg_db=-40\n.op\n.end\n",
        wires("in", n),
        wires("out", n)
    ));
    let r = solve(&deck);
    let il = 10f64.powf(-3.0 / 20.0); // gauss mode's default insertion loss
                                      // Output port j gets channel j at full strength …
    for j in 0..n {
        let p = r.node_voltage(&format!("out{j}_{j}_re")).unwrap();
        assert!((p / il - 1.0).abs() < 1e-6, "out{j} ch{j}: {p}");
        // … and the neighbouring channels at the adjacent floor, in their own
        // slots — never summed into the wanted one.
        let adj = (j + 1) % n;
        let leak = r.node_voltage(&format!("out{j}_{adj}_re")).unwrap();
        let want = il * 10f64.powf(-25.0 / 20.0);
        assert!(
            (leak / want - 1.0).abs() < 1e-6,
            "out{j} ch{adj} leakage: got {leak}, want {want}"
        );
    }
}

/// An 8×8 router — the size that actually gets used — must give the right
/// answer **on default options**.
///
/// This test used to carry `.options vntol=1e-14 reltol=1e-12`, because λ wires
/// carry ~1.55e-6 and SPICE's default absolute *voltage* tolerance of 1e-6 is
/// the same order as the entire quantity: Newton's step test was satisfied
/// while λ was still ~10 pm out, which is a real detuning for a 40 GHz
/// passband. At N ≤ 5 the first step lands accurately enough that it never
/// showed; at N = 8 it did, and the router silently reported the transmission
/// for the wrong wavelength.
///
/// λ rows now carry their own absolute tolerance (`crate::tolerance`), so the
/// deck no longer needs to know about the solver. **Keep the defaults here** —
/// running this without the override is the regression.
#[test]
fn an_eight_by_eight_router_is_exact_on_default_tolerances() {
    let n = 8;
    let lam = grid_nm(n);
    let mut deck = format!("* 8×8\n{}", sources(n, &lam));
    // A 10 GHz passband on a 100 GHz grid puts the Gaussian tails below 1e-100,
    // so each output slot carries exactly one term and the expected value is
    // closed-form. (Crosstalk accumulation is covered by the n=4 test above.)
    deck.push_str(&format!(
        "Xr{}{} fc_awgr df_ghz=100 fwhm_ghz=10 xt_adj_db=-300 xt_bg_db=-300\n",
        wires("in", n),
        wires("out", n)
    ));
    deck.push_str(".op\n.end\n");
    let r = solve(&deck);
    let il = 10f64.powf(-3.0 / 20.0);
    for i in 0..n {
        for k in 0..n {
            let j = (i + k) % n;
            let got = r.node_voltage(&format!("out{j}_{k}_re")).unwrap();
            let want = (1.0 + i as f64 + 10.0 * k as f64) * il;
            assert!(
                (got / want - 1.0).abs() < 1e-9,
                "in{i} ch{k} → out{j}: got {got}, want {want}"
            );
        }
    }
}

/// A wavelength far outside the device's band must be dark, not folded back
/// onto a passband by the FSR periodicity.
///
/// The periodicity is real but bounded — the star coupler's far field rolls off
/// after a few FSRs, and an unbounded wrap makes the model claim a full
/// passband at, say, 20 nm. That is not just unphysical: Newton's early
/// iterates put λ near zero on its way up from the initial guess, and a model
/// that reports in-band transmission there thrashes its coefficients from
/// iteration to iteration until the line search collapses.
#[test]
fn light_far_outside_the_band_is_dark_rather_than_aliased() {
    let n = 2;
    // Channel 0 at a wavelength 100 nm off the grid; channel 1 on-grid.
    let lam = vec![1450.0, grid_nm(n)[1]];
    let mut deck = format!("* out of band\n{}", sources(n, &lam));
    deck.push_str(&format!(
        "Xr{}{} fc_awgr df_ghz=100 fwhm_ghz=40\n.op\n.end\n",
        wires("in", n),
        wires("out", n)
    ));
    let r = solve(&deck);
    for j in 0..n {
        let p = r.node_voltage(&format!("out{j}_0_re")).unwrap();
        assert!(p.abs() < 1e-12, "out{j} ch0 should be dark, got {p}");
    }
    // The on-grid channel is unaffected.
    let ok = r.node_voltage("out1_1_re").unwrap();
    assert!(ok.abs() > 0.1, "on-grid channel still routes: {ok}");
}

/// A bad port/channel-count combination must fail loudly at build time rather
/// than mis-parse into a differently shaped router.
#[test]
#[should_panic(expected = "2·wpc·N²")]
fn a_terminal_count_that_is_not_two_wpc_n_squared_is_rejected() {
    // 30 wires: 30/(2·3) = 5 channels per side, and 5 is not a perfect square.
    // (24 would *not* be an error — it is exactly a 2×2 router.)
    let mut deck = String::from("* wrong shape\nV1 a_re 0 DC 1.0\n");
    deck.push_str("Xr");
    for t in 0..30 {
        deck.push_str(&format!(" n{t}"));
    }
    deck.push_str(" fc_awgr\n.op\n.end\n");
    let _ = solve(&deck);
}

/// Measured-table mode: a `.model` card pointing at an `N×N` CSV of measured
/// spectra. This is the only route for a string parameter — an X-line's
/// instance params are numeric.
#[test]
fn a_measured_table_overrides_the_analytic_response() {
    let n = 2;
    let dir = std::env::temp_dir().join("fc_awgr_table_test");
    std::fs::create_dir_all(&dir).unwrap();
    let csv = dir.join("awgr2.csv");
    // Deliberately *not* what the analytic model would produce: input 0 goes
    // to output 1 at −6 dB on channel 0, contradicting the ideal permutation.
    std::fs::write(
        &csv,
        "wavelength_nm,t_0_0_db,t_0_1_db,t_1_0_db,t_1_1_db\n\
         1549.0,-40,-6,-40,-40\n\
         1551.0,-40,-6,-40,-40\n",
    )
    .unwrap();

    let lam = vec![1550.0, 1550.0];
    let mut deck = format!("* measured AWGR\n{}", sources(n, &lam));
    deck.push_str(&format!(
        ".model awg2 fc_awgr sfile=\"{}\"\n",
        csv.to_str().unwrap()
    ));
    deck.push_str(&format!(
        "Xr{}{} awg2\n.op\n.end\n",
        wires("in", n),
        wires("out", n)
    ));
    let r = solve(&deck);
    // Input 0 (amplitude 1) reaches output 1 at −6 dB …
    let got = r.node_voltage("out1_0_re").unwrap();
    let want = 10f64.powf(-6.0 / 20.0) + 2.0 * 10f64.powf(-40.0 / 20.0); // + input 1's leak
    assert!((got / want - 1.0).abs() < 1e-6, "got {got}, want {want}");
    // … and output 0 sees only the two −40 dB terms, not the ideal route.
    let leak = r.node_voltage("out0_0_re").unwrap();
    let want_leak = 3.0 * 10f64.powf(-40.0 / 20.0);
    assert!(
        (leak / want_leak - 1.0).abs() < 1e-6,
        "table must override the permutation: got {leak}, want {want_leak}"
    );
    std::fs::remove_file(&csv).ok();
}
