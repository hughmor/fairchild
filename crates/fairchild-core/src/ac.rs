/// Small-signal AC analysis.
///
/// Algorithm:
///   1. Run DC operating-point Newton-Raphson to find x0.
///   2. Evaluate device small-signal Jacobians (linearized conductances G_dev) at x0.
///   3. Build frequency-dependent admittance matrix Y(jω) = G + jωC − j/ω·L.
///   4. For each frequency, solve Y·V = I_ac for complex node voltages.
///
/// AC excitation: sources flagged with `AC <mag> [<phase_deg>]` in the netlist.
/// The current parser stores only DC waveforms; for now every voltage/current source
/// whose DC value is non-zero is treated as the AC excitation with unit amplitude and
/// zero phase.  A proper `.ac` source tag (not yet parsed) can replace this later.
///
/// Complex solver: split into the equivalent 2N×2N real block system
///   [ G  -B ] [V_re]   [I_re]
///   [ B   G ] [V_im] = [I_im]
/// where B = ωC − L/ω.
use indexmap::IndexMap;
use rayon::prelude::*;

use fairchild_parser::{Element, Netlist};

use crate::device::{Device, EvalFlags, ReactiveKind, SimContext};
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::mna::{stamp_2port_by_id, stamp_netlist_scaled, stamp_passive_2port, CircuitTopology};
use crate::newton::build_devices;
use crate::options::SimOptions;
use crate::solver::LinearSolver;

/// Result of an AC sweep: complex node voltages at each frequency.
pub struct AcResult {
    /// Swept frequencies in Hz.
    pub freq: Vec<f64>,
    /// Complex node voltages: node name → [(re, im)] at each frequency.
    pub voltages: IndexMap<String, Vec<(f64, f64)>>,
}

impl AcResult {
    /// Magnitude |V(node)| at frequency index `fi`.
    pub fn magnitude(&self, node: &str, fi: usize) -> Option<f64> {
        let v = self.voltages.get(node)?;
        let (re, im) = v[fi];
        Some((re * re + im * im).sqrt())
    }

    /// Phase ∠V(node) in degrees at frequency index `fi`.
    pub fn phase_deg(&self, node: &str, fi: usize) -> Option<f64> {
        let v = self.voltages.get(node)?;
        let (re, im) = v[fi];
        Some(im.atan2(re).to_degrees())
    }

    /// Write the AC sweep as an ngspice-compatible Nutmeg ASCII rawfile.
    ///
    /// Complex values are written as `<re>,<im>` per ngspice convention.
    pub fn write_nutmeg<W: std::io::Write>(&self, mut w: W, title: &str) -> std::io::Result<()> {
        let n_vars = 1 + self.voltages.len();
        let n_pts = self.freq.len();
        writeln!(w, "Title: {title}")?;
        writeln!(w, "Plotname: AC Analysis")?;
        writeln!(w, "Flags: complex")?;
        writeln!(w, "No. Variables: {n_vars}")?;
        writeln!(w, "No. Points: {n_pts}")?;
        writeln!(w, "Variables:")?;
        writeln!(w, "\t0\tfrequency\tfrequency")?;
        for (i, name) in self.voltages.keys().enumerate() {
            writeln!(w, "\t{}\tv({name})\tvoltage", i + 1)?;
        }
        writeln!(w, "Values:")?;
        for (fi, &f) in self.freq.iter().enumerate() {
            // Frequency is always real; write as real,0 per ngspice complex convention.
            writeln!(w, " {fi}\t{f:.6e},0")?;
            for v in self.voltages.values() {
                let (re, im) = v[fi];
                writeln!(w, "\t{re:.6e},{im:.6e}")?;
            }
        }
        Ok(())
    }

    /// Write a CSV with columns: freq_hz, mag_V(<node>), phase_deg_V(<node>), ...
    pub fn write_csv<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        write!(w, "freq_hz")?;
        for name in self.voltages.keys() {
            write!(w, ",mag_V({name}),phase_deg_V({name})")?;
        }
        writeln!(w)?;
        for (fi, &f) in self.freq.iter().enumerate() {
            write!(w, "{f:.6e}")?;
            for v in self.voltages.values() {
                let (re, im) = v[fi];
                let mag = (re * re + im * im).sqrt();
                let phase = im.atan2(re).to_degrees();
                write!(w, ",{mag:.6e},{phase:.4}")?;
            }
            writeln!(w)?;
        }
        Ok(())
    }
}

/// Generate logarithmically spaced frequency points.
///
/// `points_per_decade` — number of points per decade (10× interval).
pub fn freq_decade(start_hz: f64, stop_hz: f64, points_per_decade: usize) -> Vec<f64> {
    let log_start = start_hz.log10();
    let log_stop = stop_hz.log10();
    let n_decades = log_stop - log_start;
    let n = ((n_decades * points_per_decade as f64).ceil() as usize).max(2);
    (0..=n)
        .map(|i| 10f64.powf(log_start + (log_stop - log_start) * i as f64 / n as f64))
        .filter(|&f| f <= stop_hz * 1.0001)
        .collect()
}

/// Generate linearly spaced frequency points.
pub fn freq_linear(start_hz: f64, stop_hz: f64, n_points: usize) -> Vec<f64> {
    if n_points < 2 {
        return vec![start_hz];
    }
    (0..n_points)
        .map(|i| start_hz + (stop_hz - start_hz) * i as f64 / (n_points - 1) as f64)
        .collect()
}

/// Generate logarithmically spaced frequency points.
///
/// `points_per_octave` — number of points per octave (2× interval).
pub fn freq_oct(start_hz: f64, stop_hz: f64, points_per_octave: usize) -> Vec<f64> {
    let log2_start = start_hz.log2();
    let log2_stop = stop_hz.log2();
    let n_octaves = log2_stop - log2_start;
    let n = ((n_octaves * points_per_octave as f64).ceil() as usize).max(2);
    (0..=n)
        .map(|i| 2f64.powf(log2_start + (log2_stop - log2_start) * i as f64 / n as f64))
        .filter(|&f| f <= stop_hz * 1.0001)
        .collect()
}

/// Run a small-signal AC sweep.
///
/// `ac_source` — name of the voltage source to use as the AC stimulus (amplitude 1 V,
///               phase 0°).  Pass `None` to drive all voltage sources simultaneously
///               with 1 V (useful for transfer-function sweeps with a single source).
pub fn ac_analysis(
    netlist: &Netlist,
    freqs: &[f64],
    ac_source: Option<&str>,
    registry: &DeviceRegistry,
) -> Result<AcResult, SimError> {
    ac_analysis_opts(
        netlist,
        freqs,
        ac_source,
        registry,
        &SimOptions::from_netlist(netlist),
    )
}

/// AC analysis with explicit `SimOptions`.
pub fn ac_analysis_opts(
    netlist: &Netlist,
    freqs: &[f64],
    ac_source: Option<&str>,
    registry: &DeviceRegistry,
    opts: &SimOptions,
) -> Result<AcResult, SimError> {
    crate::connectivity::check_connectivity(netlist)?;
    let ctx = opts.sim_context();
    let mut topo = CircuitTopology::build(netlist);
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();

    // --- DC operating point ---
    let mut devices = build_devices(netlist, &mut topo, &ctx, registry)?;
    let size = topo.size;
    let dc_solver = opts.linear_solver(size);
    let x0 = dc_op(&topo, netlist, &mut devices, &ctx, opts, &*dc_solver)?;
    // The AC system is 2N×2N (real-block of complex); build a sized solver.
    let ac_solver = opts.linear_solver(2 * size);

    // --- Small-signal G matrix (real, from DC Jacobian at x0) ---
    // Re-stamp the linear passive network.
    let mat0 = stamp_netlist_scaled(&topo, netlist, 1.0, &empty, &empty);
    // Add device linearization (Jacobian at x0).
    let mut g_mat = mat0.a;
    for dev in devices.iter_mut() {
        // Evaluate with transient flags so devices populate their small-signal
        // capacitance caches at the operating point. The resistive Jacobian
        // (`load_jacobian`) is identical under dc()/tran() flags; only the
        // cached reactances differ, and those are read below via
        // `small_signal_reactances()`.
        dev.eval(&x0, EvalFlags::tran(), &ctx);
        // Use a temporary MnaMatrix to collect resistive Jacobian entries.
        let mut tmp = crate::mna::MnaMatrix::zeros(size);
        dev.load_jacobian(&mut tmp);
        // `tmp.a` is sparse now, so this walks only the cells the device
        // actually stamped instead of the whole row.
        for (g_row, t_row) in g_mat.iter_mut().zip(tmp.a.iter()) {
            for (j, t) in t_row.iter() {
                g_row[j] += t;
            }
        }
    }
    // GMIN — node rows and device-internal rows (skips vsource aux rows).
    topo.stamp_gmin(&mut g_mat, opts.gmin);

    // ponytail: `.ac` assembles G/C/L and the 2n×2n system densely — O(n²)
    // memory and O(n²) work per frequency point. The DC/transient matrix went
    // sparse in a929041 and this did not, because it is not on that hot path.
    // Upgrade path: emit `SparseRow`s directly (the stamp primitives and
    // `CircuitTopology::to_csc` already take them); do it when AC on a large
    // photonic circuit starts hurting. Tracked as task #12.
    // --- Capacitance matrix C (purely imaginary part of Y) ---
    let mut c_mat = vec![vec![0.0f64; size]; size];
    for el in &netlist.elements {
        if let Element::Capacitor {
            pos,
            neg,
            capacitance,
            ..
        } = el
        {
            stamp_passive_2port(&mut c_mat, &topo.node_index, pos, neg, *capacitance);
        }
    }

    // --- Inductance matrix L_inv (contribution = -1/ωL to imaginary part) ---
    // For inductors in DC OP they appear as short circuits; their AC stamp is 1/(jωL).
    // We track "L" values and handle 1/ω at solve time.
    let mut l_mat = vec![vec![0.0f64; size]; size];
    for el in &netlist.elements {
        if let Element::Inductor {
            pos,
            neg,
            inductance,
            ..
        } = el
        {
            stamp_passive_2port(&mut l_mat, &topo.node_index, pos, neg, 1.0 / inductance);
        }
    }

    // --- Device-internal small-signal reactances (diode Cj, MOSFET
    // Meyer/junction caps, photonic parasitics) ---
    // These are stamped by transient (load_jacobian_tran / reactive_branches)
    // but were historically absent from AC, so device caps were ignored.
    // `small_signal_reactances()` reports the same physical reactances; the
    // transient eval above populated their cached values.
    for dev in devices.iter() {
        for r in dev.small_signal_reactances() {
            match r.kind {
                ReactiveKind::Capacitor => {
                    stamp_2port_by_id(&mut c_mat, r.pos, r.neg, r.value);
                }
                ReactiveKind::Inductor if r.value != 0.0 => {
                    stamp_2port_by_id(&mut l_mat, r.pos, r.neg, 1.0 / r.value);
                }
                ReactiveKind::Inductor => {}
            }
        }
        // Devices whose reactance is a general ∂q/∂x matrix rather than a set
        // of two-terminal branches (OSDI/Verilog-A) stamp it themselves.
        dev.load_reactive_jacobian(&mut c_mat);
    }

    // --- AC excitation vector ---
    // Build the RHS for AC: voltage sources contribute to the stub row; current sources to node rows.
    // Voltage source in MNA: stamps A[vi][p]=+1, A[vi][n]=-1, A[p][vi]=+1, A[n][vi]=-1, b[vi]=V_ac.
    // For the AC analysis, we set V_ac = 1 for the chosen source(s).
    // Unit amplitude is only the fallback for a deck that declares no `AC` spec
    // anywhere; a declared spec wins. See `build_ac_rhs`.
    let (b_ac_re, b_ac_im) = build_ac_rhs(&topo, netlist, ac_source, 1.0);

    // --- Sweep (parallel across frequencies) ---
    //
    // Each frequency builds its own 2N×2N block system and solves it.  The
    // assembled G / C / L matrices, the RHS vectors, and the linear solver
    // are all read-only inside the loop, so rayon can fan out one frequency
    // per worker.  Results are collected in input order via `par_iter().map()
    // .collect::<Vec<_>>()` and then written into the result IndexMap by
    // index — same observable output as the sequential path.
    let solver_ref: &dyn LinearSolver = &*ac_solver;
    let per_freq: Vec<Result<Vec<(f64, f64)>, SimError>> = freqs
        .par_iter()
        .map(|&f| {
            let omega = 2.0 * std::f64::consts::PI * f;

            // B matrix = ωC − L/ω  (susceptance)
            let mut b_mat = vec![vec![0.0f64; size]; size];
            for i in 0..size {
                for j in 0..size {
                    b_mat[i][j] = omega * c_mat[i][j] - l_mat[i][j] / omega;
                }
            }

            // Build 2N×2N block system:
            //   [ G  -B ] [V_re]   [I_re]
            //   [ B   G ] [V_im] = [I_im]
            let n2 = 2 * size;
            let mut a2 = vec![vec![0.0f64; n2]; n2];
            let mut rhs = vec![0.0f64; n2];

            for i in 0..size {
                for j in 0..size {
                    a2[i][j] = g_mat[i][j];
                    a2[i][size + j] = -b_mat[i][j];
                    a2[size + i][j] = b_mat[i][j];
                    a2[size + i][size + j] = g_mat[i][j];
                }
                rhs[i] = b_ac_re[i];
                rhs[size + i] = b_ac_im[i];
            }

            let x = solver_ref.solve(&CircuitTopology::sparse_from_dense(&a2), &rhs)?;
            let row: Vec<(f64, f64)> = topo
                .node_index
                .values()
                .map(|&idx| (x[idx], x[size + idx]))
                .collect();
            Ok(row)
        })
        .collect();

    // Bail on the first error to mirror the historic short-circuit behaviour.
    let per_freq: Vec<Vec<(f64, f64)>> = per_freq.into_iter().collect::<Result<Vec<_>, _>>()?;

    // Transpose: per_freq is [freq_idx][node_idx], result wants [node][freq_idx].
    let mut voltages: IndexMap<String, Vec<(f64, f64)>> = topo
        .node_index
        .keys()
        .map(|k| (k.clone(), Vec::with_capacity(freqs.len())))
        .collect();
    for row in &per_freq {
        for (series, &pair) in voltages.values_mut().zip(row.iter()) {
            series.push(pair);
        }
    }

    Ok(AcResult {
        freq: freqs.to_vec(),
        voltages,
    })
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

fn dc_op(
    topo: &CircuitTopology,
    netlist: &Netlist,
    devices: &mut [Box<dyn Device>],
    ctx: &SimContext,
    opts: &SimOptions,
    solver: &dyn LinearSolver,
) -> Result<Vec<f64>, SimError> {
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();
    let n_nodes = topo.n_nodes();
    // Not every unknown is a volt — see `crate::tolerance`.
    let tol = crate::tolerance::Tolerances::build(netlist, topo, opts);
    let mut x = vec![0.0f64; topo.size];

    for _ in 0..opts.itl1 {
        let mut mat = stamp_netlist_scaled(topo, netlist, 1.0, &empty, &empty);
        for dev in devices.iter_mut() {
            dev.eval(&x, EvalFlags::dc(), ctx);
            dev.load_residual(&mut mat.b);
            dev.load_jacobian(&mut mat);
        }
        topo.stamp_gmin(&mut mat.a, opts.gmin);
        let x_new = solver.solve(&mat.a, &mat.b)?;
        let max_dv = x_new
            .iter()
            .zip(x.iter())
            .take(n_nodes)
            .map(|(n, o)| (n - o).abs())
            .fold(0.0f64, f64::max);
        let x_next: Vec<f64> = if max_dv > opts.vmax {
            let scale = opts.vmax / max_dv;
            x.iter()
                .zip(x_new.iter())
                .map(|(o, n)| o + scale * (n - o))
                .collect()
        } else {
            x_new
        };
        let converged = tol.converged(&x_next, &x);
        x = x_next;
        if converged {
            return Ok(x);
        }
    }
    Err(SimError::NoConvergence { iters: opts.itl1 })
}

/// Build the AC excitation as `(real, imaginary)` RHS vectors.
///
/// A source's `AC <mag> [phase]` is the excitation: `mag·cos φ` into the real
/// vector, `mag·sin φ` into the imaginary one. Voltage sources drive their
/// auxiliary row; current sources drive their node rows.
///
/// **Two regimes, and the deck chooses.** When any source declares an `AC` spec,
/// the deck is read strictly: declared sources drive at their own magnitude and
/// phase, and undeclared ones contribute nothing — ngspice's rule. When *no*
/// source declares one, every source (or the one named by `ac_source`) drives at
/// `fallback_mag` with zero phase, which is what fairchild has always done and
/// what a transfer-function sweep wants; there is no declared intent to
/// contradict, so nothing is lost.
///
/// The regime split matters because the old code ignored `AC` entirely and drove
/// *everything* at unit amplitude. A deck pairing an `AC 1` drive with a plain
/// DC bias source therefore excited the bias source too — wrong, and silent.
fn build_ac_rhs(
    topo: &CircuitTopology,
    netlist: &Netlist,
    ac_source: Option<&str>,
    fallback_mag: f64,
) -> (Vec<f64>, Vec<f64>) {
    let n_nodes = topo.n_nodes();
    let mut re = vec![0.0f64; topo.size];
    let mut im = vec![0.0f64; topo.size];

    let strict = netlist.elements.iter().any(|el| {
        matches!(
            el,
            Element::VoltageSource { ac: Some(_), .. } | Element::CurrentSource { ac: Some(_), .. }
        )
    });

    for el in &netlist.elements {
        let (name, ac, nodes) = match el {
            Element::VoltageSource { name, ac, .. } => (name, ac, None),
            Element::CurrentSource {
                name, ac, pos, neg, ..
            } => (name, ac, Some((pos, neg))),
            _ => continue,
        };
        if !ac_source.is_none_or(|s| s.eq_ignore_ascii_case(name)) {
            continue;
        }
        let (mag, phase_deg) = match ac {
            Some(spec) => (spec.mag, spec.phase_deg),
            // Declared nothing in a deck that declares elsewhere: not a source.
            None if strict => continue,
            None => (fallback_mag, 0.0),
        };
        let phase = phase_deg.to_radians();
        let (er, ei) = (mag * phase.cos(), mag * phase.sin());

        match nodes {
            None => {
                if let Some(&vi) = topo.vsrc_index.get(name) {
                    re[n_nodes + vi] += er;
                    im[n_nodes + vi] += ei;
                }
            }
            // SPICE: current leaves n+ and enters n-.
            Some((pos, neg)) => {
                if let Some(&p) = topo.node_index.get(pos) {
                    re[p] -= er;
                    im[p] -= ei;
                }
                if let Some(&n) = topo.node_index.get(neg) {
                    re[n] += er;
                    im[n] += ei;
                }
            }
        }
    }
    (re, im)
}

#[cfg(test)]
mod tests {
    use super::*;
    /// `AC <mag> [phase]` on a source line must reach the solver.
    ///
    /// It used to be parsed away, so `.ac` always drove at unit amplitude and
    /// zero phase: an ngspice deck written with `AC 2` came out 2x too small
    /// with no diagnostic. Values below are ngspice 46 on the same decks.
    #[test]
    fn ac_spec_on_a_source_sets_magnitude_and_phase() {
        let deck = |spec: &str| {
            format!(
                "* ac spec\nV1 in 0 {spec}\nR1 in out 1k\nC1 out 0 159.15n\n                 .ac lin 1 1k 1k\n.end\n"
            )
        };
        let run = |spec: &str| {
            let nl = parse_spice(&deck(spec)).unwrap();
            let reg = DeviceRegistry::new();
            let r = ac_analysis(&nl, &[1e3], None, &reg).expect("ac failed");
            let (re, im) = r.voltages.get("out").unwrap()[0];
            ((re * re + im * im).sqrt(), im.atan2(re).to_degrees())
        };

        // Magnitude scales; the -3 dB corner keeps its -45 degrees.
        for (spec, want_mag) in [("DC 0 AC 1", 0.7071178), ("DC 0 AC 2", 1.4142356)] {
            let (mag, ph) = run(spec);
            assert!(
                (mag - want_mag).abs() < 1e-6,
                "{spec}: |V(out)| = {mag:.7}, expected {want_mag:.7}"
            );
            assert!((ph + 45.0).abs() < 0.01, "{spec}: phase {ph:.4} != -45");
        }

        // Phase rotates the excitation: +90 deg in, -45+90 = +45 deg out.
        let (mag, ph) = run("AC 1 90");
        assert!((mag - 0.7071178).abs() < 1e-6, "|V(out)| = {mag:.7}");
        assert!((ph - 45.0).abs() < 0.01, "phase {ph:.4} != +45");

        // A deck that declares nothing keeps the legacy unit drive, which is
        // what the one in-tree `.ac` example relies on.
        let (mag, _) = run("DC 0");
        assert!(
            (mag - 0.7071178).abs() < 1e-6,
            "no-spec fallback broke: {mag:.7}"
        );
    }

    /// Once any source declares an `AC` spec the deck is read strictly, so a
    /// plain DC bias source is no longer treated as a second AC drive. That
    /// silent double-excitation is what the old unconditional unit drive did.
    #[test]
    fn a_declared_ac_spec_makes_undeclared_sources_quiet() {
        let nl = parse_spice(
            "* two sources, one AC\n             V1 in 0 DC 0 AC 1\n             Vbias b 0 DC 2.5\n             R1 in out 1k\n             Rb b out 1k\n             C1 out 0 1n\n             .ac lin 1 1k 1k\n.end\n",
        )
        .unwrap();
        let reg = DeviceRegistry::new();
        let r = ac_analysis(&nl, &[1e3], None, &reg).expect("ac failed");
        let (re, im) = r.voltages.get("b").unwrap()[0];
        assert!(
            (re * re + im * im).sqrt() < 1e-12,
            "the DC bias source is being driven as an AC source: |V(b)| = {:.3e}",
            (re * re + im * im).sqrt()
        );
    }

    use crate::device_registry::DeviceRegistry;
    use fairchild_parser::parse_spice;

    /// RC low-pass filter: R=1kΩ, C=1nF. Cutoff at 1/(2πRC) ≈ 159 kHz.
    /// At f_c: |V(out)| = 1/√2 ≈ −3 dB.
    #[test]
    fn rc_lowpass_cutoff() {
        let net = parse_spice(
            "* RC low-pass\nVin in 0 DC 1\nR1 in out 1k\nC1 out 0 1n\n.ac DEC 10 1k 10Meg\n.end\n",
        )
        .unwrap();
        let registry = DeviceRegistry::new();
        let f_c = 1.0 / (2.0 * std::f64::consts::PI * 1e3 * 1e-9); // ≈ 159 kHz
        let freqs = freq_decade(1e3, 10e6, 20);
        let result = ac_analysis(&net, &freqs, Some("Vin"), &registry).unwrap();

        // Find the frequency closest to f_c.
        let fi = freqs
            .iter()
            .enumerate()
            .min_by_key(|(_, &f)| ((f - f_c).abs() * 1e6) as u64)
            .map(|(i, _)| i)
            .unwrap();

        let mag = result.magnitude("out", fi).unwrap();
        let mag_at_dc = result.magnitude("out", 0).unwrap();
        // At f_c: |H(jω)| = 1/√2 of DC gain (DC gain ≈ 1).
        assert!(
            (mag / mag_at_dc - std::f64::consts::FRAC_1_SQRT_2).abs() < 0.05,
            "At f={:.1}kHz: |V(out)|/|V(out)_dc|={:.4} (expected 1/√2 ≈ {:.4})",
            freqs[fi] / 1e3,
            mag / mag_at_dc,
            std::f64::consts::FRAC_1_SQRT_2
        );
    }

    /// A reverse-biased diode's junction capacitance must appear in `.ac` — it
    /// was historically dropped (AC built its C matrix from netlist capacitors
    /// only). `Vac` DC-biases the cathode to +1 V (reverse) and is the AC
    /// source; `Rs` + `Cj` form a lowpass. With the cap included there is a
    /// −3 dB rolloff at f_c = 1/(2π·Rs·Cj); without it the response is flat, so
    /// this test fails on the pre-fix code.
    #[test]
    fn diode_junction_cap_appears_in_ac() {
        let net = parse_spice(
            "* reverse-biased diode Cj lowpass\n\
             Vac c 0 DC 1\n\
             Rs c out 10k\n\
             D1 0 out DMOD\n\
             .model DMOD D IS=1e-14 CJO=2p VJ=0.8 M=0.5\n\
             .end\n",
        )
        .unwrap();
        let mut registry = DeviceRegistry::new();
        registry.register_builtin_models(&net.models);
        let freqs = freq_decade(1e3, 1e9, 30);
        let result = ac_analysis(&net, &freqs, Some("Vac"), &registry).unwrap();

        let mag_dc = result.magnitude("out", 0).unwrap();
        let mag_hi = result.magnitude("out", freqs.len() - 1).unwrap();
        // Clear rolloff: at 1 GHz the cap shorts `out` toward ground.
        assert!(
            mag_hi < 0.2 * mag_dc,
            "expected high-frequency rolloff: mag_hi={mag_hi:.4} mag_dc={mag_dc:.4}"
        );

        // -3 dB near f_c. At Vd_j ≈ −1 V: Cj = CJO/(1+1/0.8)^0.5 = 2pF/1.5.
        let cj = 2e-12_f64 / (1.0 + 1.0 / 0.8_f64).powf(0.5);
        let f_c = 1.0 / (2.0 * std::f64::consts::PI * 1e4 * cj); // ≈ 11.9 MHz
        let fi = freqs
            .iter()
            .enumerate()
            .min_by_key(|(_, &f)| ((f - f_c).abs()) as u64)
            .map(|(i, _)| i)
            .unwrap();
        let ratio = result.magnitude("out", fi).unwrap() / mag_dc;
        assert!(
            (ratio - std::f64::consts::FRAC_1_SQRT_2).abs() < 0.08,
            "at f={:.2} MHz |V(out)|/|V_dc|={:.4} (expected 1/√2 ≈ 0.707); \
             a flat response (cap dropped) would read ~1.0",
            freqs[fi] / 1e6,
            ratio
        );
    }

    #[test]
    fn write_csv_ac() {
        let net = parse_spice("* RC\nVin in 0 DC 1\nR1 in out 1k\nC1 out 0 1n\n.end\n").unwrap();
        let registry = DeviceRegistry::new();
        let freqs = freq_decade(1e3, 1e6, 5);
        let result = ac_analysis(&net, &freqs, Some("Vin"), &registry).unwrap();
        let mut buf = Vec::new();
        result.write_csv(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("freq_hz"), "header: {s}");
        assert!(s.contains("mag_V(out)"), "should have V(out): {s}");
    }

    #[test]
    fn write_nutmeg_ac() {
        let net = parse_spice("* RC\nVin in 0 DC 1\nR1 in out 1k\nC1 out 0 1n\n.end\n").unwrap();
        let registry = DeviceRegistry::new();
        let freqs = freq_decade(1e3, 1e6, 5);
        let result = ac_analysis(&net, &freqs, Some("Vin"), &registry).unwrap();
        let mut buf = Vec::new();
        result.write_nutmeg(&mut buf, "RC test").unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("Plotname: AC Analysis"), "plotname: {s}");
        assert!(s.contains("Flags: complex"), "flags: {s}");
        assert!(s.contains("frequency\tfrequency"), "freq var: {s}");
        assert!(s.contains("v(out)\tvoltage"), "v(out): {s}");
        assert!(s.contains("Values:"), "values: {s}");
        // Complex values should contain a comma separator.
        let values_section = s.split("Values:").nth(1).unwrap();
        assert!(
            values_section.contains(','),
            "complex pairs should use comma: {s}"
        );
    }
}
