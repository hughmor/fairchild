//! `fc_optical_2x2` — behavioural per-channel 2×2 optical transfer block.
//!
//! Exercises the whole path: `.optical_port` / `.electrical_port` declarations →
//! parser flattening into one bundle-aware instance → per-channel weights →
//! solved DC operating point read back off the output bundles.

use fairchild_core::device_registry::DeviceRegistry;
use fairchild_core::newton::dc_op_nr_with_registry;
use fairchild_core::{options::SimOptions, tran_nr_with_registry_var_opts};
use fairchild_parser::parse_spice;

/// Optical power |E|² on channel `k` of a bundle port.
fn power(r: &fairchild_core::newton::NrResult, port: &str, k: usize) -> f64 {
    let re = r.node_voltage(&format!("{port}_re_{k}")).unwrap();
    let im = r.node_voltage(&format!("{port}_im_{k}")).unwrap();
    re * re + im * im
}

/// One laser per channel on `bus`, plus a dark second input port. Weights are
/// whatever the caller puts in `params`; control wires are held at `v_ctl`.
fn weight_bank_netlist(n: usize, params: &str, v_ctl: &[f64]) -> String {
    let mut s = String::from("* fc_optical_2x2 weight bank\n");
    for p in ["bus", "dark", "thru", "drop"] {
        s.push_str(&format!(".optical_port {p} {n}\n"));
    }
    s.push_str(&format!(".electrical_port wctl {n}\n"));
    for k in 0..n {
        // 1 mW per channel, 1 nm apart.
        s.push_str(&format!(
            "Xl{k} bus_re_{k} bus_im_{k} bus_wl_{k} fc_cw_laser power_mW=1.0 \
             wavelength_nm={}\n",
            1550.0 + k as f64
        ));
        // The unused second input needs drivers (it is an input, not an output).
        s.push_str(&format!("Vdr{k} dark_re_{k} 0 DC 0\n"));
        s.push_str(&format!("Vdi{k} dark_im_{k} 0 DC 0\n"));
        s.push_str(&format!("Vdw{k} dark_wl_{k} 0 DC 1.55e-6\n"));
        s.push_str(&format!(
            "Vc{k} wctl_{k} 0 DC {}\n",
            v_ctl.get(k).copied().unwrap_or(0.0)
        ));
    }
    s.push_str(&format!(
        "Xwb bus dark thru drop wctl 0 fc_optical_2x2 {params}\n.op\n"
    ));
    s
}

fn solve(netlist: &str) -> fairchild_core::newton::NrResult {
    let parsed = parse_spice(netlist).expect("netlist should parse");
    dc_op_nr_with_registry(&parsed, &DeviceRegistry::new()).expect("DC OP should converge")
}

/// The headline: one instance, four wavelengths, an independent bipolar weight
/// on each, and `P_drop − P_thru = w · P_in` exactly, with power conserved.
#[test]
fn per_channel_bipolar_weights_land_exactly() {
    let want = [1.0, -1.0, 0.0, 0.5];
    let params = "w_0=1.0 w_1=-1.0 w_2=0.0 w_3=0.5";
    let r = solve(&weight_bank_netlist(4, params, &[]));
    for (k, &w) in want.iter().enumerate() {
        let (p_thru, p_drop) = (power(&r, "thru", k), power(&r, "drop", k));
        let p_in = 1e-3; // 1 mW
        assert!(
            ((p_drop - p_thru) / p_in - w).abs() < 1e-9,
            "ch{k}: weight {} expected {w}",
            (p_drop - p_thru) / p_in
        );
        assert!(
            ((p_thru + p_drop) / p_in - 1.0).abs() < 1e-9,
            "ch{k}: power {} not conserved",
            (p_thru + p_drop) / p_in
        );
    }
}

/// An unindexed parameter broadcasts to every channel.
#[test]
fn unindexed_weight_broadcasts() {
    let r = solve(&weight_bank_netlist(3, "w=0.25", &[]));
    for k in 0..3 {
        let w = (power(&r, "drop", k) - power(&r, "thru", k)) / 1e-3;
        assert!((w - 0.25).abs() < 1e-9, "ch{k}: {w}");
    }
}

/// `dw_dv` makes the weight follow its own control wire — the reason the device
/// exists, since `set_param` cannot move a weight during a transient. Each
/// channel gets a different control voltage and must respond independently.
#[test]
fn control_voltage_moves_each_weight_independently() {
    // w_k(V) = 0 + 0.5·V(wctl_k)  ⇒  −0.5, 0, +0.5, and the last clamps at +1.
    let v = [-1.0, 0.0, 1.0, 4.0];
    let r = solve(&weight_bank_netlist(4, "w=0 dw_dv=0.5", &v));
    let want = [-0.5, 0.0, 0.5, 1.0];
    for (k, &w) in want.iter().enumerate() {
        let got = (power(&r, "drop", k) - power(&r, "thru", k)) / 1e-3;
        assert!((got - w).abs() < 1e-9, "ch{k}: got {got} expected {w}");
        let total = (power(&r, "thru", k) + power(&r, "drop", k)) / 1e-3;
        assert!((total - 1.0).abs() < 1e-9, "ch{k}: clamp broke passivity");
    }
}

/// Insertion loss is power dB across both outputs.
#[test]
fn insertion_loss_applies_in_power_db() {
    let r = solve(&weight_bank_netlist(1, "w=0 il_db=3.0", &[]));
    let total = (power(&r, "thru", 0) + power(&r, "drop", 0)) / 1e-3;
    let expected = 10f64.powf(-3.0 / 10.0); // 3 dB → half the power
    assert!(
        (total - expected).abs() < 1e-6,
        "total {total} expected {expected}"
    );
}

/// Explicit-matrix mode: a plain 2×2 with a phase, bypassing weight mode.
/// s11 = 0.6, s21 = 0.8·e^{j90°} routes |0.6|² to thru and |0.8|² to drop.
#[test]
fn explicit_matrix_mode_overrides_weight() {
    let params = "w=1.0 s11_mag=0.6 s11_deg=0 s21_mag=0.8 s21_deg=90 \
                  s12_mag=0 s22_mag=0";
    let r = solve(&weight_bank_netlist(1, params, &[]));
    let (p_thru, p_drop) = (power(&r, "thru", 0) / 1e-3, power(&r, "drop", 0) / 1e-3);
    assert!((p_thru - 0.36).abs() < 1e-9, "thru {p_thru}");
    assert!((p_drop - 0.64).abs() < 1e-9, "drop {p_drop}");
}

/// The λ tag must propagate to both outputs, or a downstream wavelength-aware
/// device sees an undriven wire.
#[test]
fn wavelength_labels_pass_through_to_both_outputs() {
    let r = solve(&weight_bank_netlist(2, "w=0", &[]));
    for k in 0..2 {
        let want = (1550.0 + k as f64) * 1e-9;
        for port in ["thru", "drop"] {
            let got = r.node_voltage(&format!("{port}_wl_{k}")).unwrap();
            assert!((got - want).abs() < 1e-15, "{port} ch{k}: {got} vs {want}");
        }
    }
}

/// A control bus narrower than the optical bus is meaningless for a
/// bundle-aware device and must be caught at parse time, naming both ports.
#[test]
fn control_bus_width_mismatch_is_a_parse_error() {
    let netlist = "* mismatched widths\n\
                   .optical_port bus 4\n\
                   .optical_port dark 4\n\
                   .optical_port thru 4\n\
                   .optical_port drop 4\n\
                   .electrical_port wctl 2\n\
                   Xwb bus dark thru drop wctl 0 fc_optical_2x2 w=0\n\
                   .op\n";
    let msg = format!("{}", parse_spice(netlist).unwrap_err());
    assert!(msg.contains("same channel count"), "{msg}");
    assert!(msg.contains("bus(optical, 4 ch)"), "{msg}");
    assert!(msg.contains("wctl(electrical, 2 ch)"), "{msg}");
}

/// An explicit matrix with gain is rejected unless the user opts in — a gain
/// block in a feedback path diverges silently otherwise.
#[test]
#[should_panic(expected = "has gain")]
fn explicit_matrix_with_gain_is_rejected() {
    solve(&weight_bank_netlist(
        1,
        "s11_mag=2.0 s12_mag=0 s21_mag=0 s22_mag=0",
        &[],
    ));
}

#[test]
fn explicit_gain_allowed_when_requested() {
    let r = solve(&weight_bank_netlist(
        1,
        "allow_gain=1 s11_mag=2.0 s12_mag=0 s21_mag=0 s22_mag=0",
        &[],
    ));
    let p_thru = power(&r, "thru", 0) / 1e-3;
    assert!((p_thru - 4.0).abs() < 1e-9, "thru {p_thru}");
}

/// Bidirectional mode isn't modelled yet; fail loudly rather than leave the
/// backward wires undriven. (No laser here — in 5-wire mode its terminal names
/// differ, and the device's assert fires during circuit build anyway.)
#[test]
#[should_panic(expected = "bidirectional propagation is not supported")]
fn bidirectional_mode_is_rejected() {
    let netlist = "* bidirectional 2x2\n\
                   .options enable_bidirectional=1\n\
                   .optical_port bus 1\n\
                   .optical_port dark 1\n\
                   .optical_port thru 1\n\
                   .optical_port drop 1\n\
                   .electrical_port wctl 1\n\
                   Vc0 wctl_0 0 DC 0\n\
                   Xwb bus dark thru drop wctl 0 fc_optical_2x2 w=0\n\
                   .op\n";
    solve(netlist);
}

/// `tau_s` must actually delay the output, not just be accepted as a parameter.
///
/// What gets delayed is the *optical field*: the output is `S(t) · in(t − τ)`,
/// matching `OpticalSegment` — the matrix is evaluated at the current time while
/// the fields carry the latency. So the stimulus here is a step on the input
/// field (driven straight onto the bundle wires, since a CW laser's power is a
/// static parameter), with the weight held fixed.
///
/// With `w = 0` the block is a 50/50 split and `s21 = −j·sin(π/4)`, so a real
/// input step lands entirely on `drop_im` at −0.7071·in1_re.
#[test]
fn latency_delays_the_output() {
    let cross_time = |tau_s: f64| -> (f64, f64) {
        let netlist = format!(
            "* 2x2 latency\n\
             .optical_port bus 1\n\
             .optical_port dark 1\n\
             .optical_port thru 1\n\
             .optical_port drop 1\n\
             .electrical_port wctl 1\n\
             Vin bus_re_0 0 PULSE(0 1 1n 1p 1p 20n 40n)\n\
             Vini bus_im_0 0 DC 0\n\
             Vinw bus_wl_0 0 DC 1.55e-6\n\
             Vdr0 dark_re_0 0 DC 0\n\
             Vdi0 dark_im_0 0 DC 0\n\
             Vdw0 dark_wl_0 0 DC 1.55e-6\n\
             Vc0 wctl_0 0 DC 0\n\
             Xwb bus dark thru drop wctl 0 fc_optical_2x2 w=0 tau_s={tau_s}\n\
             .tran 50p 8n\n"
        );
        let net = parse_spice(&netlist).unwrap();
        let opts = SimOptions::from_netlist(&net);
        let r = tran_nr_with_registry_var_opts(&net, 50e-12, 8e-9, &DeviceRegistry::new(), &opts)
            .expect("tran must complete");
        let (t, d) = (&r.time, r.node_voltages.get("drop_im_0").unwrap());
        let final_v = *d.last().unwrap();
        let half = 0.5 * final_v;
        let cross = (1..d.len())
            .find(|&i| (d[i] - half).signum() != (d[i - 1] - half).signum())
            .map(|i| t[i])
            .unwrap_or(f64::NAN);
        (cross, final_v)
    };

    let (t_fast, v_fast) = cross_time(0.0);
    let (t_slow, v_slow) = cross_time(2e-9);
    assert!(
        t_fast.is_finite() && t_slow.is_finite(),
        "no transition found (fast {t_fast}, slow {t_slow})"
    );
    // −1/√2 · 1 V either way: latency shifts in time, it does not attenuate.
    let expected = -1.0 / 2f64.sqrt();
    for (label, v) in [("fast", v_fast), ("slow", v_slow)] {
        assert!(
            (v - expected).abs() < 1e-9,
            "{label}: drop_im settled at {v}, expected {expected}"
        );
    }
    // The delayed run crosses ≈2 ns later (one timestep of slack).
    let shift = t_slow - t_fast;
    assert!(
        (shift - 2e-9).abs() < 1e-10,
        "expected ≈2 ns delay, got {:.3} ns (fast {:.3} ns, slow {:.3} ns)",
        shift * 1e9,
        t_fast * 1e9,
        t_slow * 1e9
    );
}
