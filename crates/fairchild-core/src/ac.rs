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

use crate::device::{Device, EvalFlags, SimContext};
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::mna::{stamp_netlist_scaled, CircuitTopology};
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
    let n_nodes = topo.n_nodes();
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
        dev.eval(&x0, EvalFlags::dc(), &ctx);
        // Use a temporary MnaMatrix to collect Jacobian entries.
        let mut tmp = crate::mna::MnaMatrix {
            a: vec![vec![0.0; size]; size],
            b: vec![0.0; size],
        };
        dev.load_jacobian(&mut tmp);
        for i in 0..size {
            for j in 0..size {
                g_mat[i][j] += tmp.a[i][j];
            }
        }
    }
    // GMIN — node rows and OSDI internal-node rows.
    for i in 0..n_nodes {
        g_mat[i][i] += opts.gmin;
    }
    let vsrc_end = n_nodes + topo.vsrc_index.len();
    for i in vsrc_end..size {
        g_mat[i][i] += opts.gmin;
    }

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

    // --- AC excitation vector ---
    // Build the RHS for AC: voltage sources contribute to the stub row; current sources to node rows.
    // Voltage source in MNA: stamps A[vi][p]=+1, A[vi][n]=-1, A[p][vi]=+1, A[n][vi]=-1, b[vi]=V_ac.
    // For the AC analysis, we set V_ac = 1 for the chosen source(s).
    let b_ac_re = build_ac_rhs(&topo, netlist, ac_source, 1.0, 0.0); // unit amplitude, 0° phase
    let b_ac_im = vec![0.0f64; size];

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

            let x = solver_ref.solve(&a2, &rhs)?;
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
    let mut x = vec![0.0f64; topo.size];

    for _ in 0..opts.itl1 {
        let mut mat = stamp_netlist_scaled(topo, netlist, 1.0, &empty, &empty);
        for dev in devices.iter_mut() {
            dev.eval(&x, EvalFlags::dc(), ctx);
            dev.load_residual(&mut mat.b);
            dev.load_jacobian(&mut mat);
        }
        for i in 0..n_nodes {
            mat.a[i][i] += opts.gmin;
        }
        let vsrc_end = n_nodes + topo.vsrc_index.len();
        for i in vsrc_end..topo.size {
            mat.a[i][i] += opts.gmin;
        }
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
        let converged = x_next
            .iter()
            .zip(x.iter())
            .all(|(n, o)| (n - o).abs() < opts.vntol + opts.reltol * n.abs());
        x = x_next;
        if converged {
            return Ok(x);
        }
    }
    Err(SimError::NoConvergence { iters: opts.itl1 })
}

/// Stamp a 2-terminal passive element value into a matrix (G or C).
/// Same pattern as stamp_conductance but directly into a raw matrix.
fn stamp_passive_2port(
    mat: &mut Vec<Vec<f64>>,
    idx: &IndexMap<String, usize>,
    pos: &str,
    neg: &str,
    val: f64,
) {
    if let Some(&p) = idx.get(pos) {
        mat[p][p] += val;
        if let Some(&n) = idx.get(neg) {
            mat[p][n] -= val;
            mat[n][p] -= val;
        }
    }
    if let Some(&n) = idx.get(neg) {
        mat[n][n] += val;
    }
}

/// Build the AC RHS vector.
///
/// For the selected source (or all sources if ac_source=None), stamp the AC voltage/current.
/// Voltage sources stamp into the auxiliary row; current sources stamp directly into node rows.
fn build_ac_rhs(
    topo: &CircuitTopology,
    netlist: &Netlist,
    ac_source: Option<&str>,
    mag: f64,
    _phase_rad: f64,
) -> Vec<f64> {
    let n_nodes = topo.n_nodes();
    let mut b = vec![0.0f64; topo.size];
    for el in &netlist.elements {
        match el {
            Element::VoltageSource { name, .. } => {
                let drives = ac_source.map_or(true, |s| s.eq_ignore_ascii_case(name));
                if drives {
                    if let Some(&vi_idx) = topo.vsrc_index.get(name) {
                        b[n_nodes + vi_idx] += mag;
                    }
                }
            }
            Element::CurrentSource { name, pos, neg, .. } => {
                let drives = ac_source.map_or(true, |s| s.eq_ignore_ascii_case(name));
                if drives {
                    if let Some(&p) = topo.node_index.get(pos) {
                        b[p] -= mag;
                    }
                    if let Some(&n) = topo.node_index.get(neg) {
                        b[n] += mag;
                    }
                }
            }
            _ => {}
        }
    }
    b
}

#[cfg(test)]
mod tests {
    use super::*;
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
