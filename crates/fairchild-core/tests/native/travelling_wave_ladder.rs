//! A travelling-wave modulator as an interleaved ladder, against closed forms.
//!
//! The construction (#116): chop the device into `N` slices, give each one an
//! optical segment and an electrical line section, and drive slice `i`'s phase
//! from node `i` of that line. Velocity mismatch, termination ripple and
//! forward-versus-reverse drive are then *emergent* — nothing here computes a
//! walk-off integral, and the only inputs are two refractive indices and a
//! length.
//!
//! Every assertion is against an analytic form, not against another fairchild
//! path. The whole point of the exercise is that the ladder reproduces physics
//! nobody stamped into it.
//!
//! Nothing here is a new device. The ladder is forty lines of deck built from
//! `fc_pn_ps` and `T`, which is what makes it a cheap way to validate the
//! physics *before* `fc_tw_ps` exists.
//!
//! # What is not covered
//!
//! The RF-loss case — `H(ω) = (1 − e^(−(α+jΔβ)L)) / ((α+jΔβ)L)` — needs a lossy
//! line, and `T` is lossless (`O`/LTRA is not implemented; see
//! `docs/spice_support.md`). Adding series resistors between sections gives a
//! loss that is only first-order right and brings its own dispersion, which
//! would be checking an approximation against a different approximation.

use fairchild_core::{ac_analysis, DeviceRegistry, SimOptions};
use fairchild_parser::parse_spice;

const C0: f64 = 299_792_458.0;
const N_SLICES: usize = 24;
const L_M: f64 = 3e-3; // 3 mm
const N_G: f64 = 4.2; // optical group index
const Z0: f64 = 35.0; // electrode characteristic impedance
const V_PI_L: f64 = 0.012;

/// How the RF is launched relative to the light.
#[derive(Clone, Copy, PartialEq)]
enum Drive {
    /// RF enters at the same end the light does — the useful case.
    Co,
    /// RF enters at the far end, so the two sweep past each other.
    Counter,
}

/// The interleaved ladder, as a deck.
///
/// Slice `i` holds an optical segment of length `L/N` whose phase follows node
/// `e{i}` of the electrode, and an electrode section of one-way delay
/// `n_m·dz/c`. The optical group delay per slice is the segment's own
/// `τ_g = n_g·dz/c`, which is why this needs `waveguide_delay=1`.
fn ladder(n_m: f64, drive: Drive, z_term: f64) -> String {
    let dz = L_M / N_SLICES as f64;
    let td = n_m * dz / C0;
    let l_um = dz * 1e6;
    let mut s = String::from("* travelling-wave phase shifter as an interleaved ladder\n");
    s.push_str(".options waveguide_delay=1\n");
    // Optical bundle nets: o0 (in) … oN (out).
    for i in 0..=N_SLICES {
        s.push_str(&format!(".optical_port o{i}\n"));
    }
    s.push_str("Xlas o0 fc_cw_laser power_mW=1.0 wavelength_nm=1550\n");
    // The electrode. Node e0 is the launch end for a co-propagating drive and
    // the termination for a counter-propagating one.
    let (src_node, term_node) = match drive {
        Drive::Co => (0, N_SLICES),
        Drive::Counter => (N_SLICES, 0),
    };
    s.push_str(&format!("Vrf rf 0 DC 0 AC 1\nRs rf e{src_node} {Z0}\n"));
    s.push_str(&format!("Rt e{term_node} 0 {z_term}\n"));
    for i in 0..N_SLICES {
        s.push_str(&format!("T{i} e{i} 0 e{} 0 Z0={Z0} TD={td:.9e}\n", i + 1));
    }
    // The optical slices. `g_pn=0` keeps the electrode unloaded: the whole
    // point is to watch the RF propagate, and 24 shunt conductances would be a
    // different transmission line. The electrode still has a DC path to the
    // source, because a lossless line conducts at DC.
    for i in 0..N_SLICES {
        s.push_str(&format!(
            "Xps{i} o{i} o{} e{i} 0 fc_pn_ps l_um={l_um:.6} v_pi_l={V_PI_L} \
             n_g={N_G} g_pn=0 alpha_dB_cm=0 pin_at_ref=1\n",
            i + 1
        ));
    }
    s
}

/// Magnitude of the optical field's small-signal response at the ladder output.
///
/// A phase modulator's response is in quadrature with the carrier, so this
/// takes the whole complex field rather than one quadrature — the answer is
/// then independent of where the carrier's phase happens to sit.
fn eo_response(deck: &str, freqs: &[f64]) -> Vec<f64> {
    let parsed = parse_spice(deck).expect("parse");
    let opts = SimOptions::from_netlist(&parsed);
    let ac = ac_analysis(&parsed, freqs, Some("vrf"), &DeviceRegistry::new())
        .or_else(|_| {
            fairchild_core::ac::ac_analysis_opts(
                &parsed,
                freqs,
                Some("vrf"),
                &DeviceRegistry::new(),
                &opts,
            )
        })
        .expect("ac sweep");
    let out = format!("o{N_SLICES}");
    (0..freqs.len())
        .map(|i| {
            let re = ac.magnitude(&format!("{out}_re_0"), i).expect("field re");
            let im = ac.magnitude(&format!("{out}_im_0"), i).expect("field im");
            (re * re + im * im).sqrt()
        })
        .collect()
}

/// `|sinc(x)| = |sin x / x|`, the walk-off form, with `x = ω·L·Δn/(2c)`.
fn walkoff(f: f64, delta_n: f64) -> f64 {
    let x = std::f64::consts::PI * f * L_M * delta_n / C0;
    if x.abs() < 1e-12 {
        1.0
    } else {
        (x.sin() / x).abs()
    }
}

/// Step 1 of the validation ladder: matched velocities, matched termination.
///
/// With `n_m = n_g` every slice's contribution arrives at the output in phase,
/// whatever the frequency, so the response is flat. Nothing else in this deck
/// has a pole — no junction capacitance, no shunt conductance — so a roll-off
/// here would be the interleaving being mis-phased rather than physics.
#[test]
fn velocity_matched_ladder_is_flat() {
    let freqs: Vec<f64> = (0..7).map(|k| 1e9 + k as f64 * 10e9).collect();
    let r = eo_response(&ladder(N_G, Drive::Co, Z0), &freqs);
    for (f, v) in freqs.iter().zip(&r) {
        let rel = v / r[0];
        assert!(
            (rel - 1.0).abs() < 0.02,
            "matched ladder must be flat: {rel:.4} of its low-frequency value at \
             {:.0} GHz. A tilt here is the interleaving, not the physics.",
            f / 1e9
        );
    }
}

/// Step 2: velocity mismatch, lossless. The response must follow
/// `sinc(ω·L·(n_m − n_g)/2c)` — the walk-off integral, which nothing in the
/// deck computes.
#[test]
fn velocity_mismatch_follows_the_walkoff_sinc() {
    const N_M: f64 = 2.1; // a fast electrode: Δn = 2.1, first null near 47 GHz
    let delta_n = N_G - N_M;
    // Up to the first null and a little past it.
    let freqs: Vec<f64> = (1..=12).map(|k| k as f64 * 5e9).collect();
    let r = eo_response(&ladder(N_M, Drive::Co, Z0), &freqs);
    let norm = r[0] / walkoff(freqs[0], delta_n);
    for (f, v) in freqs.iter().zip(&r) {
        let want = walkoff(*f, delta_n) * norm;
        // 4 % of the *peak*, not of the local value: near a null the closed
        // form goes through zero and a relative test there is meaningless.
        assert!(
            (v - want).abs() < 0.04 * norm,
            "at {:.0} GHz the ladder gives {v:.4e}, sinc says {want:.4e}",
            f / 1e9
        );
    }
    // The null is the sharp end of the claim: it must be deep, and it must be
    // where the closed form puts it (c/(L·Δn) = 47.6 GHz for these numbers).
    let f_null = C0 / (L_M * delta_n);
    let at_null = eo_response(&ladder(N_M, Drive::Co, Z0), &[f_null]);
    assert!(
        at_null[0] < 0.06 * norm,
        "the walk-off null at {:.1} GHz is only down to {:.3} of the peak",
        f_null / 1e9,
        at_null[0] / norm
    );
}

/// Step 4, and the sharpest of them: reverse the RF drive.
///
/// Co-propagating, the walk-off is set by `n_m − n_g`; counter-propagating, by
/// `n_m + n_g`, because the two now sweep past each other. The bandwidth
/// collapses by that ratio. Nothing in the deck knows which direction the light
/// travels — the asymmetry comes from the ladder's own topology, which is why
/// this cannot pass by accident.
#[test]
fn reversing_the_drive_collapses_the_bandwidth() {
    const N_M: f64 = 2.1;
    let f3db = |drive: Drive, delta_n: f64| {
        // Walk the sweep out to the first 3 dB point of the response.
        let freqs: Vec<f64> = (1..=200).map(|k| k as f64 * 0.5e9).collect();
        let r = eo_response(&ladder(N_M, drive, Z0), &freqs);
        let half = r[0] / std::f64::consts::SQRT_2;
        let k = r.iter().position(|&v| v < half).unwrap_or(freqs.len() - 1);
        let f = freqs[k];
        // The closed form's own 3 dB point, for the message.
        let want = 1.39 * C0 / (std::f64::consts::PI * L_M * delta_n);
        (f, want)
    };
    let (f_co, want_co) = f3db(Drive::Co, N_G - N_M);
    let (f_ct, want_ct) = f3db(Drive::Counter, N_G + N_M);

    assert!(
        (f_co / want_co - 1.0).abs() < 0.15,
        "co-propagating 3 dB at {:.1} GHz, sinc says {:.1} GHz",
        f_co / 1e9,
        want_co / 1e9
    );
    assert!(
        (f_ct / want_ct - 1.0).abs() < 0.15,
        "counter-propagating 3 dB at {:.1} GHz, sinc says {:.1} GHz",
        f_ct / 1e9,
        want_ct / 1e9
    );
    // The ratio is the physics, and it is (n_g + n_m)/(n_g − n_m) = 3.0 here.
    let ratio = f_co / f_ct;
    let want_ratio = (N_G + N_M) / (N_G - N_M);
    assert!(
        (ratio / want_ratio - 1.0).abs() < 0.2,
        "reversing the drive should cost a factor {want_ratio:.2} of bandwidth, \
         measured {ratio:.2} ({:.1} GHz vs {:.1} GHz)",
        f_co / 1e9,
        f_ct / 1e9
    );
}

/// Step 5: detune the termination and the electrode rings.
///
/// A reflection at the far end returns after `2·n_m·L/c`, so the response
/// acquires a ripple of period `c/(2·n_m·L)`. Measuring the period is what
/// distinguishes a real standing wave from any other wiggle.
#[test]
fn a_detuned_termination_ripples_at_the_round_trip_period() {
    const N_M: f64 = 4.2; // velocity matched, so the *only* structure is the ripple
    let period = C0 / (2.0 * N_M * L_M);
    // Sample four periods finely enough to locate maxima.
    let n = 160;
    let freqs: Vec<f64> = (1..=n)
        .map(|k| k as f64 * 4.0 * period / n as f64)
        .collect();
    let matched = eo_response(&ladder(N_M, Drive::Co, Z0), &freqs);
    let detuned = eo_response(&ladder(N_M, Drive::Co, Z0 / 2.0), &freqs);

    // The matched line is the control: same deck, same sweep, no ripple.
    let flatness = |r: &[f64]| {
        let (lo, hi) = r
            .iter()
            .fold((f64::MAX, 0.0f64), |(a, b), &v| (a.min(v), b.max(v)));
        (hi - lo) / hi
    };
    assert!(
        flatness(&matched) < 0.02,
        "the matched ladder must not ripple: {:.1}% peak-to-peak",
        100.0 * flatness(&matched)
    );
    assert!(
        flatness(&detuned) > 0.05,
        "a 2:1 termination mismatch must ripple: only {:.1}% peak-to-peak",
        100.0 * flatness(&detuned)
    );

    // Period: the spacing of interior maxima.
    let peaks: Vec<f64> = (1..detuned.len() - 1)
        .filter(|&i| detuned[i] > detuned[i - 1] && detuned[i] > detuned[i + 1])
        .map(|i| freqs[i])
        .collect();
    assert!(
        peaks.len() >= 2,
        "expected several ripple maxima across four periods, found {}",
        peaks.len()
    );
    let spacing = (peaks[peaks.len() - 1] - peaks[0]) / (peaks.len() - 1) as f64;
    assert!(
        (spacing / period - 1.0).abs() < 0.15,
        "ripple period {:.2} GHz, round trip says {:.2} GHz",
        spacing / 1e9,
        period / 1e9
    );
}

// ── fc_tw_ps: the same construction, built by the device ─────────────────────

/// The device deck: one element where the hand-written ladder needs `2N + 3`
/// lines. `n_slices` is set explicitly so the two are the same circuit — the
/// derived count is checked separately.
fn tw_device(n_m: f64, n_slices: usize, z_term: f64) -> String {
    format!(
        "* fc_tw_ps\n\
         .options waveguide_delay=1\n\
         .optical_port oi\n\
         .optical_port oo\n\
         Xlas oi fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
         Vrf rf 0 DC 0 AC 1\n\
         Rs rf e0 {Z0}\n\
         Rt eN 0 {z_term}\n\
         Xtw oi oo e0 eN fc_tw_ps l_um={l_um} v_pi_l={V_PI_L} n_g={N_G} \
         n_m={n_m} z0={Z0} n_slices={n_slices} alpha_dB_cm=0 pin_at_ref=1\n",
        l_um = L_M * 1e6,
    )
}

fn device_response(deck: &str, freqs: &[f64]) -> Vec<f64> {
    let parsed = parse_spice(deck).expect("parse");
    let ac = ac_analysis(&parsed, freqs, Some("vrf"), &DeviceRegistry::new()).expect("ac sweep");
    (0..freqs.len())
        .map(|i| {
            let re = ac.magnitude("oo_re_0", i).expect("field re");
            let im = ac.magnitude("oo_im_0", i).expect("field im");
            (re * re + im * im).sqrt()
        })
        .collect()
}

/// `fc_tw_ps` is the hand-written ladder, to the digit.
///
/// Not an anchor — both sides are fairchild — but it is the right test for
/// *this* claim, which is that the device builds the construction the closed
/// forms above already validated. The physics is pinned upstream; what is at
/// risk here is the wiring, and a wiring error shows up as a mismatch with the
/// deck that is known to be right.
#[test]
fn the_device_reproduces_the_hand_written_ladder() {
    const N_M: f64 = 2.1;
    let freqs: Vec<f64> = (1..=10).map(|k| k as f64 * 5e9).collect();
    let hand = eo_response(&ladder(N_M, Drive::Co, Z0), &freqs);
    let dev = device_response(&tw_device(N_M, N_SLICES, Z0), &freqs);
    let scale = dev[0] / hand[0];
    assert!(
        (scale - 1.0).abs() < 1e-6,
        "low-frequency response differs by {scale:.6}x before any walk-off"
    );
    for (i, f) in freqs.iter().enumerate() {
        let rel = (dev[i] - hand[i]).abs() / hand[0];
        assert!(
            rel < 1e-6,
            "at {:.0} GHz the device gives {:.6e} and the hand ladder {:.6e}",
            f / 1e9,
            dev[i],
            hand[i]
        );
    }
}

/// The slice count follows `f_max`, and the answer converges as it rises.
///
/// The convergence sweep is the point: `N` is a discretisation, and a
/// discretisation you cannot check is a number you are trusting. Doubling
/// `slices_per_wave` must move the answer by less than the last doubling did.
#[test]
fn the_slice_count_follows_f_max_and_the_answer_converges() {
    const N_M: f64 = 2.1;
    // Well inside the band, where the walk-off is strong enough that an
    // under-resolved ladder would show it.
    let f = 30e9;
    let mut last: Option<f64> = None;
    let mut deltas = Vec::new();
    for n in [6, 12, 24, 48] {
        let r = device_response(&tw_device(N_M, n, Z0), &[f])[0];
        if let Some(prev) = last {
            deltas.push((r - prev).abs() / r);
        }
        last = Some(r);
    }
    assert!(
        deltas.windows(2).all(|w| w[1] < w[0]),
        "refining the ladder must converge; successive changes were {deltas:?}"
    );
    assert!(
        *deltas.last().unwrap() < 0.01,
        "24 -> 48 slices still moves the answer by {:.2}%",
        100.0 * deltas.last().unwrap()
    );

    // And the derived count tracks f_max: a wider band asks for more slices,
    // which is visible as a smaller requested timestep.
    let build = |f_max: f64| {
        let deck = format!(
            "* derived slice count\n\
             .optical_port oi\n\
             .optical_port oo\n\
             Xlas oi fc_cw_laser power_mW=1.0 wavelength_nm=1550\n\
             Vrf rf 0 DC 0\n\
             Rs rf e0 {Z0}\n\
             Rt eN 0 {Z0}\n\
             Xtw oi oo e0 eN fc_tw_ps l_um={l_um} v_pi_l={V_PI_L} n_g={N_G} n_m={N_M} \
             z0={Z0} f_max={f_max:e} alpha_dB_cm=0 pin_at_ref=1\n",
            l_um = L_M * 1e6,
        );
        let parsed = parse_spice(&deck).expect("parse");
        // The requested step is tau per slice over two, so it is a direct read
        // of the slice count without exposing one.
        let opts = SimOptions::from_netlist(&parsed);
        let mut topo = fairchild_core::mna::CircuitTopology::build_resolved(
            &parsed,
            &opts.sim_context(),
            &DeviceRegistry::new(),
        );
        let devices = fairchild_core::newton::build_devices(
            &parsed,
            &mut topo,
            &opts.sim_context(),
            &DeviceRegistry::new(),
        )
        .expect("build");
        devices
            .iter()
            .filter_map(|d| d.requested_max_timestep())
            .fold(f64::INFINITY, f64::min)
    };
    let coarse = build(10e9);
    let fine = build(80e9);
    assert!(
        fine < coarse / 4.0,
        "8x the bandwidth should ask for a much finer ladder: {coarse:.3e} s vs {fine:.3e} s"
    );
}
